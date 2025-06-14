package Bio::Adventure::Prepare;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Cwd qw"abs_path getcwd cwd";
use Bio::DB::EUtilities;
use File::Basename;
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
        library => 'nucleotide',
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

    my $eutil = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                         -email => $options->{email},
                                         -db => $options->{library},
                                         -id => \@unique,);
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
