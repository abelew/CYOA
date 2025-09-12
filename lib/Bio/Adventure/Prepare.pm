package Bio::Adventure::Prepare;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Cwd qw"abs_path getcwd cwd";
use Archive::Zip qw":ERROR_CODES :CONSTANTS";
use Bio::DB::EUtilities;
use File::Basename;
use File::Copy qw"move";
use File::Path qw"make_path rmtree";
use File::Which qw"which";
use HTML::TreeBuilder::XPath;
use List::Util qw"uniq";
use Text::CSV;
use Text::CSV_XS::TSV;
use WWW::Mechanize;
use XML::LibXML;

=head1 NAME

 Bio::Adventure::Prepare - Get raw data ready for processing.

=head1 SYNOPSIS

 use Bio::Adventure;
 my $hpgl = new Bio::Adventure;
 $hpgl->Prepare(csv => 'all_samples.csv');

=head1 Methods

=cut

## Taking the wall of code out of Phage::Filter_Kraken_Worker
## I am hoping to make this more robust and useful elsewhere.

sub Download_Ensembl_Files {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        overwrite => 0,
        species => 'mus_musculus',
    );
    my $sleeper = 3;
    my $release = 0;
    my $overwrite = $options->{overwrite};
    if (defined($options->{release})) {
        $release = $options->{release};
    } elsif (defined($options->{version})) {
        $release = $options->{version};
    }
    my $got = 0;
    my $mech = WWW::Mechanize->new(autocheck => 1);
    my $ens_root = qq"https://ftp.ensembl.org/pub/";
    unless ($release) {
        my $current_release = 0;
        my @releases = ();
        print "Checking ensembl for the latest release.\n";
        my $root_index = $mech->get($ens_root);
        $got = $mech->success();
        if ($got) {
            my @ens_index = $mech->links();
          ENS_LINKS: for my $i (@ens_index) {
                my $root_url = $i->url();
                next ENS_LINKS unless ($root_url =~ /^release\-/);
                my $release_num = $root_url;
                $release_num =~ s/^.*release\-(\d+).*$/$1/;
                push(@releases, $release_num);
                if ($release_num > $current_release) {
                    $current_release = $release_num;
                }
            }
        }
        $release = $current_release;
        print "Automatically chose release: ${release}.\n";
        sleep($sleeper);
    }
    print "Downloading files from ftp.ensembl.org for species: $options->{species}, release: ${release}.\n";

    ## Gather and concatenate the genbank annotations.
    my $gb_root = qq"${ens_root}/release-${release}/genbank/$options->{species}/";
    print "Checking for ensembl genbank flat files.\n";
    my $gb_filename = qq"$options->{species}_${release}_ens.gb.gz";
    my $gb_downloaded = 0;
    if (-f $gb_filename && !$overwrite) {
        print "The genbank file already exists and overwrite is off, skipping.\n";
    } else {
        if (-f $gb_filename) {
            rm($gb_filename);
        }
        my $gb_index = $mech->get($gb_root);
        $got = $mech->success();
        if ($got) {
            my @gb_index = $mech->links();
            my $gb_concatenated = FileHandle->new(">>${gb_filename}");
          GB_LINKS: for my $i (@gb_index) {
                my $link_url = $i->url();
                if ($link_url =~ /\.dat\.gz$/) {
                    sleep $sleeper;
                    my $full_url = qq"${gb_root}${link_url}";
                    my $dl = $mech->get($full_url, ':content_file' => $gb_concatenated);
                    $got = $mech->success();
                    if ($got) {
                        print "Downloaded chromosome ${link_url}.\n";
                        $gb_downloaded++;
                    } ## Checking success of single file download.
                } ## Checking for a link to a chromosomal gb file.
            } ## Iterating over files in the directory.
            $gb_concatenated->close();
        } ## Checking if a single file downloaded.
        if ($gb_downloaded < 1) {
            print "No genbank flat files successfully downloaded, deleting the concatenated file.\n";
            rm($gb_filename);
        }
    } ## Checking overwrite status.

    ## Download the master gff3 file.
    my $gff3_root = qq"https://ftp.ensembl.org/pub/release-${release}/gff3/$options->{species}/";
    my $gff_file = qq"$options->{species}_${release}_ens.gff3.gz";
    if (-f $gff_file && !$overwrite) {
        print "The gff3 file already exists and overwrite is off, skipping.\n";
    } else {
        print "Checking for an ensembl gff3 file.\n";
        my $gff_index = $mech->get($gff3_root);
        $got = $mech->success();
        if ($got) {
            my @gff_index = $mech->links();
          GFF_LINKS: for my $i (@gff_index) {
                my $link_url = $i->url();
                if ($link_url =~ /\.${release}\.gff3\.gz$/) {
                    sleep $sleeper;
                    my $gff_out = FileHandle->new(">${gff_file}");
                    my $full_url = qq"${gff3_root}${link_url}";
                    my $dl = $mech->get($full_url, ':content_file' => $gff_out);
                    $gff_out->close();
                    $got = $mech->success();
                    if ($got) {
                        print "Downloaded ${gff_file}.\n";
                    } else {
                        rm($gff_file);
                        print "Failed to download ${gff_file}, deleting it.\n";
                    }
                    last GFF_LINKS;
                } ## Checking links for gff file.
            } ## Iterating over links
        } else {  ## Found the gff3 directory.
            print "Unable to download the gff3 file.\n";
        }
    } ## Checking overwrite status.

    ## Download the cDNA fasta file
    my $cdna_root = qq"https://ftp.ensembl.org/pub/release-${release}/fasta/$options->{species}/cdna";
    my $cdna_file = qq"$options->{species}_${release}_ens_cdna.ffn.gz";
    if (-f $cdna_file && !$overwrite) {
        print "The cdna file exists and overwrite is off, skipping.\n";
    } else {
        print "Checking for an ensembl cdna file.\n";
        my $cdna_index = $mech->get($cdna_root);
        $got = $mech->success();
        if ($got) {
            my @cdna_index = $mech->links();
          CDNA_LINKS: for my $i (@cdna_index) {
                my $link_url = $i->url();
                if ($link_url =~ /\.cdna\.all\.fa\.gz$/) {
                    sleep $sleeper;
                    my $cdna_out = FileHandle->new(">${cdna_file}");
                    my $full_url = qq"${cdna_root}/${link_url}";
                    my $dl = $mech->get($full_url, ':content_file' => $cdna_out);
                    $cdna_out->close();
                    $got = $mech->success();
                    if ($got) {
                        print "Downloaded ${cdna_file}.\n";
                    } else {
                        rm($cdna_file);
                        print "Failed to download ${cdna_file}, deleting it.\n";
                    }
                    last CDNA_LINKS;
                }
            }
        } else {
            print "Unable to download the cdna file.\n";
        }
    }

    ## Download the cds fasta file
    my $cds_root = qq"https://ftp.ensembl.org/pub/release-${release}/fasta/$options->{species}/cds";
    my $cds_file = qq"$options->{species}_${release}_ens.ffn.gz";
    if (-f $cds_file && !$overwrite) {
        print "The CDS file already exists and overwrite is off, skipping.\n";
    } else {
        print "Checking for the ensembl CDS file.\n";
        my $cds_index = $mech->get($cds_root);
        $got = $mech->success();
        if ($got) {
            my @cds_index = $mech->links();
            my $cds_out = FileHandle->new(">${cds_file}");
          CDS_LINKS: for my $i (@cds_index) {
                my $link_url = $i->url();
                if ($link_url =~ /\.cds\.all\.fa\.gz$/) {
                    sleep $sleeper;
                    my $full_url = qq"${cds_root}/${link_url}";
                    my $dl = $mech->get($full_url, ':content_file' => $cds_out);
                    $cds_out->close();
                    $got = $mech->success();
                    if ($got) {
                        print "Downloaded ${cds_file}.\n";
                    } else {
                        rm($cds_out);
                        print "Failed to download ${cds_file}, deleting it.\n";
                    }
                    last CDS_LINKS;
                }
            }
        } else {
            print "Unable to download the CDS file.\n";
        }
    }

    ## Download the ncrna fasta file
    my $ncrna_root = qq"https://ftp.ensembl.org/pub/release-${release}/fasta/$options->{species}/ncrna";
    my $ncrna_file = qq"$options->{species}_${release}_ncrna_ens.fasta.gz";
    if (-f $ncrna_file && !$overwrite) {
        print "The ncRNA file exists and overwrite is off, skipping.\n";
    } else {
        print "Checking for the ensembl ncRNA file.\n";
        my $ncrna_index = $mech->get($ncrna_root);
        $got = $mech->success();
        if ($got) {
            my @ncrna_index = $mech->links();
          NC_LINKS: for my $i (@ncrna_index) {
                my $link_url = $i->url();
                if ($link_url =~ /\.ncrna\.fa\.gz$/) {
                    sleep $sleeper;
                    my $full_url = qq"${ncrna_root}/${link_url}";
                    my $ncrna_out = FileHandle->new(">${ncrna_file}");
                    my $dl = $mech->get($full_url, ':content_file' => $ncrna_out);
                    $ncrna_out->close();
                    $got = $mech->success();
                    if ($got) {
                        print "Downloaded ${ncrna_file}.\n";
                    } else {
                        rm($ncrna_file);
                        print "Failed to download ${ncrna_file}, deleting it.\n";
                    }
                    last NC_LINKS;
                }
            }
        } else {
            print "Unable to download the ncRNA file.\n";
        }
    }

    ## Download the peptide fasta
    my $pep_root = qq"https://ftp.ensembl.org/pub/release-${release}/fasta/$options->{species}/pep";
    my $pep_file = qq"$options->{species}_${release}_ens.faa.gz";
    if (-f $pep_file && !$overwrite) {
        print "The peptide file exists and overwrite is off, skipping.\n";
    } else {
        print "Checking for the ensembl peptide file.\n";
        my $pep_index = $mech->get($pep_root);
        $got = $mech->success();
        if ($got) {
            my @pep_index = $mech->links();
          PEP_LINKS: for my $i (@pep_index) {
                my $link_url = $i->url();
                if ($link_url =~ /\.pep\.all\.fa\.gz$/) {
                    sleep $sleeper;
                    my $full_url = qq"${pep_root}/${link_url}";
                    my $pep_out = FileHandle->new(">${pep_file}");
                    my $dl = $mech->get($full_url, ':content_file' => $pep_out);
                    $pep_out->close();
                    $got = $mech->success();
                    if ($got) {
                        print "Downloaded ${pep_file}.\n";
                    } else {
                        rm($pep_out);
                        print "Failed to download ${pep_file}, deleting it.\n";
                    }
                    last PEP_LINKS;
                }
            }
        } else {
            print "Unable to download the peptide file.\n";
        }
    }

    ## Download the genome fasta file(s)
    ## Here is the set in which I am likely interested and a little information from the README
    ## 1.  The set of alts: .dna.alt.fa.gz
    ## 2.  Primary assembly: .dna.primary_assembly.fa.gz
    ##     all toplevel sequences excluding haplotypes and patches; use this when doing analyses
    ##     which would be confused by them.  These may not always exist.
    ## 4.  Toplevel assembly: .dna.toplevel.fa.gz
    ##     toplevel assemblies includ anything flagged as toplevel in ensembl: chromosomes,
    ##     unassembled regions, N-padded haplotypes/patches.
    ## 5.  masked assembly alts: .dna_rm.alt.fa.gz
    ## 6.  masked assembly toplevel: .dna_rm.toplevel.fa.gz
    ##     Any rm files are hardmasked with repeatmasker and have interspersed and low-complexity
    ##     regions filled in with N
    ## 7.  masked assembly primary: .dna_rm.primary.fa.gz
    ## 8.  softmasked assembly alts: .dna_sm.alt.fa.gz
    ## 9.  softmasked assembly toplevel: .dna_sm.toplevel.fa.gz
    ##     Ibid except masked by lowercasing instead of N
    ## 10. softmasked assembly primary: .dna_sm.primary.fa.gz
    ## 11-n: There may also be explicit haplotype chromosomes and all the chromosomes separately.
    my %dna_files = (
        dna_alt => '\.dna\.alt\.fa\.gz',
        primary => '\.dna\.primary_assembly\.fa\.gz',
        toplevel => '\.dna\.toplevel\.fa\.gz',
        masked_alt => '\.dna_rm\.alt\.fa\.gz',
        masked_toplevel => '\.dna_rm\.toplevel\.fa\.gz',
        masked_primary => '\.dna_rm\.primary\.fa\.gz',
        soft_alt => '\.dna_sm\.alt\.fa\.gz',
        soft_toplevel => '\.dna_sm\.toplevel\.fa\.gz',
        soft_primary => '\.dna_sm\.primary\.fa\.gz',
    );
    my $dna_root = qq"https://ftp.ensembl.org/pub/release-${release}/fasta/$options->{species}/dna";
    my $dna_index = $mech->get($dna_root);
    print "Checking for a directory of genome fasta files.\n";
    $got = $mech->success();
    if ($got) {
        my @dna_index = $mech->links();
      TYPES: for my $type (keys %dna_files) {
          DNA_LINKS: for my $i (@dna_index) {
                my $link_url = $i->url();
                if ($link_url =~ /$dna_files{$type}$/) {
                    sleep $sleeper;
                    my $full_url = qq"${dna_root}/${link_url}";
                    my $type_file = qq"$options->{species}_${release}_${type}_ens.fasta.gz";
                    my $type_out = FileHandle->new(">${type_file}");
                    my $dl = $mech->get($full_url, ':content_file' => $type_out);
                    $got = $mech->success();
                    if ($got) {
                        print "Downloaded ${type_file}.\n";
                    } else {
                        rm($type_out);
                        print "Failed to download ${type_file}, deleting it.\n";
                    }
                    $type_out->close();
                    last DNA_LINKS;
                }
            }
        }
    } else {
        print "Unable to download the genome fasta files.\n";
    }
    ## Create a nice return here.
}

sub Download_NCBI_Assembly {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        required => ['input'],
        term => 'ORGN',
    );
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $output_dir = $paths->{output_dir};
    my $gbk_dir = $paths->{gbk_dir};
    ## The final output filenames, which may be invoked early
    ## if it is not possible to perform the filtering.
    my $log = FileHandle->new(">$paths->{log}");
    my $search_string = $options->{input};
    print $log "Downloading new assembly by searching for ${search_string}.\n";
    my $escaped_search = $search_string;
    $escaped_search =~ s/\s/\+/g;
    my $search_data = undef;
    my $downloaded_file;
    my $full_search = qq"${escaped_search}";
    if ($options->{term} ne 'none') {
        $full_search .= qq"[$options->{term}]";
    }
    my $fact = Bio::DB::EUtilities->new(-eutil => 'esearch',
                                        -email => $options->{email},
                                        -db => 'assembly',
                                        -usehistory => 'y',
                                        -term => $full_search,);
    my $hist = $fact->next_History;
    my $count = $fact->get_count;
    my @ids = $fact->get_ids();
    my $accession;
    if ($count) {
        print $log "The search string: ${search_string} returned ${count} ids.\n";
        print "The search string: ${search_string} returned ${count} ids.\n";
        $fact->reset_parameters(-eutil => 'esummary',
                                -history => $hist,
                                -email => $options->{email},
                                -db => 'assembly',
                                -retmax => 1,
                                -id => \@ids,);
        my $xml_string = $fact->get_Response->content;
        my $xml = XML::LibXML->new();
        my $dom = $xml->load_xml(string => $xml_string);
        my $uid;
        for my $node ($dom->findnodes('//*[@uid]')) {
            $uid = $node->getAttribute('uid');
        }
        $accession = $dom->getElementsByTagName('AssemblyAccession');
        my @test = split(/\./, $accession);
        my $url = $dom->getElementsByTagName('FtpPath_GenBank');
        print $log "Found ${accession} with uid ${uid} at ${url}.\n";
        print "Found ${accession} with uid ${uid} at ${url}.\n";
        ## This url is a ftp link with a nice series of downloadable files.
        my $mech = WWW::Mechanize->new(autocheck => 1);
        my @suffixes = ('_ani_contam_ranges.tsv', '_ani_report.txt',
                        '_assembly_report.txt', '_assembly_stats.txt',
                        '_cds_from_genomic.fna.gz', '_fcs_report.txt',
                        '_feature_count.txt', '_feature_table.txt.gz',
                        '_genomic.fna.gz', '_genomic.gbff.gz',
                        '_genomic.gtf.gz', '_protein.faa.gz',
                        '_protein.gpff.gz', '_rna_from_genomic.fna.gz',
                        '_translated_cds.faa.gz');
        my %downloads;
        for my $suffix (@suffixes) {
            my $name = $suffix;
            $name =~ s/^_//g;
            $name = basename($name, ('.gz'));
            my $url_basename = basename($url);
            my $url_suffix = qq"${url_basename}${suffix}";
            $downloads{$name} = qq"${url}/${url_suffix}";
        }
        print $log "The full download url for the genbank file is: $downloads{'genomic.gbff'}.\n";
        print "The full download url for the genbank file is: $downloads{'genomic.gbff'}.\n";
        $downloaded_file = qq"$paths->{gbk_dir}/${accession}.gbff.gz";
        print $log "Downloading to ${downloaded_file}.\n";
        $mech->get($downloads{'genomic.gbff'});
        $mech->save_content($downloaded_file);
    } else {
        print STDOUT "The search for an assembly for ${search_string} failed.\n";
    }
    if (!-r $downloaded_file) {
        print "Unable to download assembly information for: ${accession}.\n";
    }
}

=head2 C<Download_NCBI_Accession>

 Given an accession, download the fasta/genbank/etc file.

 This function expects an input argument which is the NCBI accession.
 It currently expects to find that accession in the nucleotide database
 (which should be made a parameter).

 A somewhat different implementation of this exists in Phage.pm where it
 is used to download an assembly.

=cut
sub Download_NCBI_Accession {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        ## library => 'nucleotide',
        jprefix => '11',);
    my $job_name = $class->Get_Job_Name();
    ## Make an array of the accession(s)
    my @unique = ();
    if (-f $options->{input}) {
        my $ids = FileHandle->new("<$options->{input}");
        my @id_array;
        while (my $line = <$ids>) {
            chomp $line;
            push(@id_array, $line);
        }
        $ids->close();
    } else {
        push(@unique, $options->{input});
    }
    @unique = uniq(@unique);

    ## Ordered list of databases I think might have my accession.
    my @dbs = ('nuccore', 'genome', 'bioproject', 'biosample', 'gene',
               'books', 'protein', 'snp', 'structure',
               'cdd', 'gap', 'dbvar', 'gds', 'geoprofiles',
               'homologene', 'mesh', 'toolkit', 'nlmcatalog', 'popset',
               'probe', 'proteinclusters', 'pcassay', 'pccompound',
               'pcsubstance', 'pmc', 'taxonomy',);
    my $found = 0;
    my $eutil;
  DBLOOP: for my $db (@dbs) {
        print "Submitting search for @unique to $db.\n";
        $eutil = Bio::DB::EUtilities->new(
            -eutil => 'esearch', -email => $options->{email},
            -retmax => 10000,
            -db => $db, -term => $options->{input},);
        my @found_ids = $eutil->get_ids;
        use Data::Dumper;
        if (scalar(@found_ids) >= 1) {
            print "Found @unique in ${db} with @found_ids\n";
            $found++;
            $eutil = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                              -email => $options->{email},
                                              -db => $db,
                                              -id => \@found_ids,);
            last DBLOOP;
        } else {
            print "Did not find @unique in ${db}.\n";
        }
        sleep(5);
    }
    if ($found == 0) {
        print "Unable to find the accession: @unique in any database.\n";
        return(undef);
    }

    while (my $docsum = $eutil->next_DocSum) {
        my $acc_version = '';
        my $accession = '';
      ITEMS: while (my $item = $docsum->next_Item) {
            my $item_name = $item->get_name;
            if ($item_name eq 'AccessionVersion') {
                $acc_version = $item->get_content();
            } elsif ($item_name eq 'Caption') {
                $accession = $item->get_content();
            } else {
                next ITEMS;
            }
        } ## End checking the document summary

        ## Now check if we already have this file
        if (-r qq"${acc_version}.gb") {
            print "Already have: ${acc_version}.gb\n";
        } else {
            print "Downloading ${accession}\n";
            my $download = Bio::DB::EUtilities->new(-eutil => 'efetch',
                                                    -db => $options->{library},
                                                    -rettype => 'gb',
                                                    -email => $options->{email},
                                                    -id => $accession,);
            my $output_file = qq"${acc_version}.gb";
            $download->get_Response(-file => $output_file);
            sleep(1);

            my @current_files = glob(qq"${acc_version}*");
            my ($first_acc, $first_ver, $first_ext);
            my ($second_acc, $second_ver, $second_ext);
            if (scalar(@current_files) > 1) {
                ($first_acc, $first_ver, $first_ext) = split(/\./, $current_files[0]);
                ($second_acc, $second_ver, $second_ext) = split(/\./, $current_files[1]);
                if ($first_ver > $second_ver) {
                    unlink($current_files[1]);
                } else {
                    unlink($current_files[0]);
                }
            }
        } ## Finished checking if we already have this accession
    } ## Finished iterating over every phage ID
}

=head2 C<Download_SRA_PRJNA>

  Extract information about and download every sample of a Bioproject.

  Given a PRJNA ID, this will query the entrezutilities for it and do the following:

  1. Download and extract as much metadata as possible.
  2. Extract every biosample accession in the bioproject.
  3. Download the _last_ biosample's SRA run.

The reason I explicitly say this downloads the last biosample
accession is due to the section below with the tag 'ITEMS', it
overwrites the accession variable with each iteration's sra
accession. If we want every SRA run associated with every biosample,
then this should check if $sra_info{$accession} is already defined,
and if so make another entry.

A good test bioproject: PRJNA799123 and the sample SRR17686731 vs. SRR17686732

=cut
sub Download_SRA_PRJNA {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jname => 'prjnadownload',
        overwrite => 1,
        jprefix => '00',);
    my %sra_info = ();
    my $accession = $options->{input};
    my $output_dir = qq"outputs/$options->{jprefix}$options->{jname}";
    my $output_log = qq"${output_dir}/prjna_download.stdout";
    my $output_csv = qq"${output_dir}/$options->{input}_sra.csv";
    make_path($output_dir);
    my $log = FileHandle->new(">${output_log}");
    my $csv_fh = FileHandle->new(">$output_csv");
    print $log "Beginning download of SRA samples associated with: $options->{input}.\n";
    my $eutil = Bio::DB::EUtilities->new(
        -eutil => 'esearch', -email => $options->{email},
        -retmax => 10000,
        -db => 'sra', -term => $accession,);
    my $id_count = $eutil->get_count;
    my $id_max = $eutil->get_retmax;
    if (!defined($id_count)) {
        print $log "Did not find any SRA accession for this accession.\n";
    } else {
        print $log "Initial search found ${id_count} accessions.\n";
    }
    my @ids = $eutil->get_ids;
    my $summary = Bio::DB::EUtilities->new(
        -eutil => 'esummary', -email => $options->{email},
        -retmax => 10000,
        -db => 'sra', -id => \@ids);
    my $count = 0;
    my $item_count = 0;
  DOCS: while (my $docsum = $summary->next_DocSum) {
        my $id_info = {
            uid => $ids[$count],
            num_samples => 0,
            sample_id_list => '',
        };
        my $accession;
      ITEMS: while (my $item = $docsum->next_Item) {
            ## This prints stuff like 'Runs ExtLinks CreateDate etc' followed by the data associated therein.
            $item_count++;
            $id_info->{num_samples}++;
            my $name = $item->get_name;
            if ($name eq 'Runs') {
                my $stuff = $item->get_content;
                $accession = $stuff;
                $accession =~ s/^.*acc="(\w+?)".*$/$1/g;
                my $spots = $stuff;
                $spots =~ s/.*total_spots="(\d+?)".*/$1/g;
                $id_info->{total_spots} = $spots;
                my $bases = $stuff;
                $bases =~ s/.*total_bases="(\d+?)".*/$1/g;
                $id_info->{total_bases} = $bases;
            } elsif ($name eq 'ExpXml') {
                my $stuff = $item->get_content;
                my $title = $stuff;
                $title =~ s/.*\<Title\>(.*?)\<\/Title\>.*/$1/g;
                $id_info->{title} = $title;
                my $instrument = $stuff;
                $instrument =~ s/.*Platform instrument_model="(.*?)"\>.*/$1/g;
                $id_info->{instrument} = $instrument;
                my $experiment = $stuff;
                $experiment =~ s/.*Experiment acc="(.*?)".*$/$1/g;
                $id_info->{experiment_acc} = $experiment;
                my $study = $stuff;
                $study =~ s/.*Study acc="(.*?)".*/$1/g;
                $id_info->{study_acc} = $study;
                my $name = $stuff;
                $name =~ s/.* name="(.*?)".*/$1/g;
                $id_info->{study_name} = $name;
                my $taxid = $stuff;
                $taxid =~ s/.*Organism taxid="(\d+?)".*/$1/g;
                $id_info->{taxid} = $taxid;
                my $species = $stuff;
                $species =~ s/.*ScientificName="(.*?)".*/$1/g;
                $id_info->{species} = $species;
                my $sample_acc = $stuff;
                $sample_acc =~ s/.*Sample acc="(.*?)".*/$1/g;
                $id_info->{sample_acc} = $sample_acc;
                if (defined($sra_info{$accession})) {
                    print "$accession has already been defined.\n";
                    $id_info->{sample_id_list} .= qq", ${sample_acc}";
                } else {
                    $id_info->{sample_id_list} = qq"${sample_acc}";
                }
                my $protocol = $stuff;
                $protocol =~ s/.*\<LIBRARY_CONSTRUCTION_PROTOCOL\>(.*?)\<\/LIBRARY_CON.*/$1/g;
                $id_info->{protocol} = $protocol;
            } elsif ($name eq 'ExtLinks') {
                my $stuff = $item->get_content;
                my $title = $stuff;
                $title =~ s/.*\<Title\>(.*?)\<\/Title\>.*/$1/g;
                $id_info->{links} = $title;
            } elsif ($name eq 'CreateDate') {
                my $stuff = $item->get_content;
                my $title = $stuff;
                $title =~ s/.*\<Title\>(.*?)\<\/Title\>.*/$1/g;
                $id_info->{create} = $title;
            } elsif ($name eq 'UpdateDate') {
                my $stuff = $item->get_content;
                my $title = $stuff;
                $title =~ s/.*\<Title\>(.*?)\<\/Title\>.*/$1/g;
                $id_info->{update} = $title;
            } else {
                print "This item has a name I haven't parsed yet: $name\n";
            }
        } ## End iterating over items
        $sra_info{$accession} = $id_info;
        $count++;
    } ## End iterating over the document summary

    ## Print the output csv file.
    my @order = ('accession', 'total_spots', 'total_bases', 'experiment_acc',
                 'study_acc', 'study_name', 'instrument', 'taxid', 'species',
                 'sample_acc', 'title', 'protocol');
    my $acc_count = 0;
    for my $acc (keys %sra_info) {
        $acc_count++;
        my @values = ($acc);
        my $inner = $sra_info{$acc};
        for my $k2 (@order) {
            push(@values, $inner->{$k2});
        }
        if ($acc_count == 1) {
            for my $h (@order) {
                print $csv_fh qq"${h}\t";
            }
            print $csv_fh "\n";
        }
        for my $v (@values) {
            print $csv_fh qq"${v}\t" if (defined($v));
        }
        print $csv_fh "\n";

        if ($options->{download}) {
            print "Beginning download to $options->{preprocess_dir}\n" if ($acc_count == 1);
            my $start = cwd();
            my $pre_dir = qq"${start}/$options->{preprocess_dir}";
            my $made = make_path($pre_dir);
            chdir($pre_dir);
            my $downloaded = $class->Bio::Adventure::Prepare::Fastq_Dump(input => $acc);
            chdir($start);
        } else {
            print "Not downloading accession: ${acc}.\n" if ($acc_count == 1);
        }
    }
    $csv_fh->close();
    $log->close();
}

sub Download_NCBIDatasets_Accession {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        delete => 0,
        required => ['input'],
        jprefix => '11',);
    my $job_name = $class->Get_Job_Name();
    ## Make an array of the accession(s)
    my $mech = WWW::Mechanize->new(autocheck => 1);
    my $url = qq"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$options->{input}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GBFF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED";
    print "Downloading data for $options->{input}.\n";
    my $downloaded_file = qq"$options->{input}.zip";
    if (-f $downloaded_file) {
        print "Already downloaded ${downloaded_file}.\n";
    } else {
        $mech->get($url);
        $mech->save_content($downloaded_file);
    }
    print "Extracting downloaded .zip file.\n";
    ## I was going to use Archive::Zip for this, but it was being a PITA.
    my $writer = FileHandle->new("unzip -o ${downloaded_file} |");
  EXTRACT: while (my $line = <$writer>) {
        chomp $line;
        next EXTRACT unless ($line =~ /inflating/);
        my $extracted = $line;
        $extracted =~ s/^\s+inflating:\s*//g;
        $extracted =~ s/\s+$//g;
        my $base = basename($extracted);
        if ($base eq 'protein.faa') {
            $base = qq"$options->{input}.faa";
        } elsif ($base eq 'genomic.gff') {
            $base = qq"$options->{input}.gff";
        } elsif ($base eq 'genomic.gbff') {
            $base = qq"$options->{input}.gbff";
        } elsif ($base eq 'cds_from_genomic.fna') {
            $base = qq"$options->{input}.ffn";
        } elsif ($base =~ /_genomic\.fna/) {
            $base = qq"$options->{input}.fsa";
        }
        my $moved = move($extracted, $base);
    }
    $writer->close();
    print "Deleting downloaded file: ${downloaded_file}.\n";
    my $gzipped = qx"gzip -9 $options->{input}.gbff";
    if ($options->{delete}) {
        unlink $downloaded_file;
    }
    print "Deleting ncbi_dataset directory.\n";
    my $removed = rmtree("ncbi_dataset");
}

=head2 C<Fastq_Dump>

 Invoke fastq_dump to download some data from sra.

 This expects an input argument which is the individual sample accession.
 In its current form it will dump fastq files as paired reads.

=cut
sub Fastq_Dump {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        output => undef);
    my $fastq_comment = qq"## This script should download an sra accession to local fastq.gz files.
";
    my $job_basename = $class->Get_Job_Name();
    my @inputs = split(/$options->{delimiter}/, $options->{input});
    my @outputs = ();
    my $first_output = undef;
    if (defined($options->{output})) {
        @outputs = split(/$options->{delimiter}/, $options->{output});
        $first_output = $outputs[0];
    }

    my %fastq_jobs = ();
    my $count = 0;
    my $fastq_job;
    for my $i (0 .. $#inputs) {
        $count++;
        my $in = $inputs[$i];
        print "Invoking fastq-dump for ${in}.\n";
        my $stderr = qq"${in}/${in}_fastqdump.stderr";
        my $stdout = qq"${in}/${in}_fastqdump.stdout";
        my $jstring = "";
        if (defined($outputs[$i]) || defined($first_output)) {
            $outputs[$i] = $first_output if (!defined($outputs[$i]));
            $jstring = qq"
mkdir -p $outputs[$i]
prefetch ${in} 2>${stderr} 1>${stdout}
fastq-dump --outdir $outputs[$i] \\
  --gzip --skip-technical --readids \\
  --read-filter pass --dumpbase \\
  --split-3 --clip ${in} \\
  2>>${stderr} 1>>${stdout}
";
        } else {
            $jstring = qq"
mkdir -p ${in}
prefetch ${in} 2>${stderr} 1>${stdout}
fastq-dump --outdir ${in} --gzip --skip-technical --readids \\
  --read-filter pass --dumpbase \\
  --split-3 --clip ${in} \\
  2>>${stderr} 1>>${stdout}
";
        }
        ## This works for single ended, but not paired.
        my $output_files = qq"$options->{input}/$options->{input}_pass.fastq.gz";
        my $output_paired = qq"$options->{input}/$options->{input}_pass_1.fastq.gz:$options->{input}/$options->{input}_pass_2.fastq.gz";
        my $current_fastq_job = $class->Submit(
            comment => $fastq_comment,
            input => $in,
            jdepends => $options->{jdepends},
            jname => qq"fqd_${in}",
            jstring => $jstring,
            jprefix => "01",
            jmem => 12,
            jwalltime => '6:00:00',
            output => $output_files,
            output_paired => $output_paired,
            stdout => $stdout,
            stderr => $stderr,);

        if (defined($fastq_job)) {
            $fastq_job->{$count} = $current_fastq_job;
        } else {
            $fastq_job = $current_fastq_job;
        }

        ## That is weird, I thought I implemented this, but no.
        ## Add a job which figures out if the fastq-dump results are se or paired.
        #my $decision_string = qq?
        #use Bio::Adventure::Prepare;
        #my \$result = \$h->Bio::Adventure::Prepare::Write_Input_Worker(
        #  input => '$fastq_job->{output}',
        #  input_paired => '$fastq_job->{output_paired}',
        #  output => '<input.txt');
        #?;
        #        my $input_worker = $class->Submit(
        #            input => $fastq_job->{output},
        #            input_paired => $fastq_job->{output_paired},
        #            jprefix => $options->{jprefix} + 1,
        #            jstring => $decision_string,
        #            jname => 'se_paired',
        #            output => '<input.txt');
    }   ## Foreach my $input
    return($fastq_job);
}

=head2 C<Read_Samples>

 This function currently has no real use-case.  It should be merged with the
 xlsx reader/writer which is used by the assembly annotation merger in order
 to provide a more flexible system for dealing with sample sheets/etc.

 With that in mind, this will dump an input csv file to a 2d hash and
 return that hash.

=cut
sub Read_Samples {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input']);
    my $fh = new FileHandle("<$options->{input}");
    my $csv = Text::CSV->new({binary => 1});
    my $row_count = 0;
    my @headers = ();
    my $data = {};
    while (my $row = $csv->getline($fh)) {
        $row_count++;
        if ($row_count == 1) {
            @headers = @{$row};
        } else {
            my $id = $row->[0];
            foreach my $c (1 .. $#headers) {
                my $key = $headers[$c];
                $data->{$id}->{$key} = $row->[$c];
            }
        }
    }
    $csv->eof or $csv->error_diag();
    $fh->close();
    return($data);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<fastq-dump>

=cut

1;
