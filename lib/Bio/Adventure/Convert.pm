package Bio::Adventure::Convert;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use feature 'try';
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::FeatureIO;
use Bio::Tools::GFF;
use Bio::Root::Exception;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Path qw"make_path remove_tree";
use File::Which qw"which";
use List::MoreUtils qw"uniq";
use Text::CSV_XS::TSV;

no warnings qw"experimental::try";

=head1 NAME

 Bio::Adventure::Convert - Perform conversions between various formats.


=head1 SYNOPSIS

 The functions here handle the various likely conversions one may perform.
 sam to bam, gff to fasta, genbank to fasta, etc.

=head1 METHODS

=cut
sub Any2Any {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    my $in_fh = FileHandle->new("less $options{input} |");
    my $input_seq = Bio::SeqIO->new(-fh => $in_fh,
                                    -format => $options{informat},);
    my $output_seq = Bio::SeqIO->new(-file => ">$options{main}",
                                     -format => $options{outformat},);
    while (my $input_seq = $input_seq->next_seq) {
        $output_seq->write_seq($input_seq);
    }
    $in_fh->close();
}

=head2 C<Bam2Bigwig>

Convert a set of aligned reads to a bigwig coverage map.

=cut
sub Bam2Bigwig {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    $options{extra_args} = $class->Passthrough_Args();
    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    my $bam2bw_job = $class->Submit();
    return($bam2bw_job);
}

=head2 C<Bam2Fastq>

Convert an aligned set of reads back to fastq via samtools and bedtools.

=cut
sub Bam2Fastq {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    my $output_name = basename($options{input}, ('.bam', '.sam'));
    $options{output_string} = q{};
    if ($options{paired}) {
        $options{output_string} = qq"-fq $options{output_dir}/${output_name}_1.fastq -fq2 $options{output_dir}/${output_name}_2.fastq";
    } else {
        $options{output_string} = qq"-fq $options{output_dir}/${output_name}_1.fastq";
    }
    $options{resorted} = qq"$options{output_dir}/resorted.bam";

    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    my $convert = $class->Submit();
    return($convert);
}

sub Extract_Subseq {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    my $in_fh = FileHandle->new("less $options{input} |");
    my $input_seq = Bio::SeqIO->new(-fh => $in_fh,
                                    -format => $options{informat},);
    my $output_seq = Bio::SeqIO->new(-file   => qq">$options{main}",
                                     -format => $options{outformat},);
    my $seq_count = 0;
    my $found_id = 0;
  SEQ_LOOP: while (my $input_seq = $input_seq->next_seq) {
        my $sequence_id = $input_seq->id;
        if (defined($options{seqid})) {
            if ($options{seqid} eq $sequence_id) {
                $found_id = 1;
            }
            ## If seqid is not provided, assume there is only 1 sequence and use the first.
        } else {
            $found_id = 1;
        }
        if ($found_id) {
            my $subseq_region = $input_seq->subseq($options{start}, $options{end});
            my $subseq_seq = Bio::Seq->new(
                -seq => $subseq_region,
                -display_id => qq"${sequence_id}_$options{start}_$options{end}");
            $output_seq->write_seq($subseq_seq);
            last SEQ_LOOP;
        }
    }
    $in_fh->close();
}

=head2 C<Extract_SeqIDs>

Yank specific Sequences by ID from an input file of arbitrary type (defaulting to fasta).

=cut
sub Extract_SeqIDs {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    my @ids = split(/$options{delimiter}/, $options{seqid});
    my $in_fh = FileHandle->new("less $options{input} |");
    my $input_seq = Bio::SeqIO->new(-fh => $in_fh,
                                    -format => $options{informat},);
    my $output_seq = Bio::SeqIO->new(-file   => qq">$options{main}",
                                     -format => $options{outformat},);
    my $seq_count = 0;
    my $found_id = 0;
  SEQ_LOOP: while (my $input_seq = $input_seq->next_seq) {
        my $sequence_id = $input_seq->id;
      WANTED: for my $wanted_id (@ids) {
            if ($sequence_id =~ /$wanted_id/) {
                $output_seq->write_seq($input_seq);
            }
        }
    }
    $in_fh->close();
}

=head2 C<Fasta2Gff>

  Convert a fasta file to a simple gff file.

  This is really only intended to make a fasta database of stuff like
  rRNAs into a simple gff file so that it is easier to count up rRNA
  reads observed per molecule per sample.  As such it makes no attempt
  to do anything useful other than write out the entry names and lengths.

=cut
sub Fasta2Gff {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    my $dir = dirname($options{input});
    my $base = basename($options{input}, ('.fasta', '.fsa', '.fna', '.fa'));
    my $gff_name = qq"${dir}/${base}.gff";
    my $gff_opened = open(GFF, ">${gff_name}");
    my $features;
    my $seqio_in = Bio::SeqIO->new(-format => 'fasta', -file => $options{input});
    my $gff_out = Bio::Tools::GFF->new(-fh => \*GFF, -gff_version => 3);
    my $features_written = 0;
  SEQS: while (my $seq = $seqio_in->next_seq) {
      my $len = $seq->length;
      my $id = $seq->id;
      print "Writing entry for $id\n";
      my $fasta_feature = Bio::SeqFeature::Generic->new(
          -primary => 'misc_feature',
          -seq_id => $id,
          -source => 'fasta_file',
          -start => 1,
          -end => $len,
          -strand => +1,
          -score => undef,
          -frame => 0,
          -tag => {
              'ID' => $id,
          },);
      $fasta_feature->display_name($id);
      $gff_out->write_feature($fasta_feature);
      $features_written++;
  }
    return($features_written);
}

=head2 C<GATK_Dedup>

 An invocation of GATK's mark duplicate function.

 Given an existing sorted bam alignment, this will attempt to mark and remove
 duplicates without the help of UMIs.

=cut
sub GATK_Dedup {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    if ($options{remove}) {
        $options{remove_string} = 'true';
    } else {
        $options{remove_string} = 'false'
    }
    $options{marked} = qq"$options{output_dir}/deduplication_stats.txt";
    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    my $gatk_job = $class->Submit();
    return($gatk_job);
}

=head2 C<Gb2Gff>

 A poorly named genbank flat file converter.

 Gb2Gff() takes a genbank genome file and splits it into:
 A genomic fasta file, CDS fasta, peptide fasta, gff file of all
 entries, CDS, and interCDS regions.

=cut
sub Gb2Gff {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    $options{perl_job} = qq"$options{basedir}/$options{script_dir}/$options{jprefix}$options{jname}.pl";
    $options{worker} = 'Gb2Gff_Worker';
    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    my $convert_job = $class->Submit();
    return($convert_job);
}

sub Gb2Gff_Worker {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    my $handle = IO::Handle->new;
    open($handle, "less $options{input} |");
    my $base = basename($options{input}, ('.xz', '.gz', '.bz2'));
    $base = basename($base, ('.gb', '.gba', '.gbk', '.genbank', '.gbff',));
    $options{output_dir} = $inputs->[0]->{directory};
    $options{output_base} = $base;
    ## By default, put all the outputs in the same directory as the input.
    print "TESTME: I expect output_dir and output_base to be passed from the original input:
  $options{output_dir} and $options{output_base}\n";
    make_path($options{output_dir}) unless(-d $options{output_dir});
    $options{fasta} = qq"$options{output_dir}/$options{output_base}.fsa";
    $options{cds_fasta} = qq"$options{output_dir}/$options{output_base}.ffn";
    $options{misc_fasta} = qq"$options{output_dir}/$options{output_base}.ffn";
    $options{nc_fasta} = qq"$options{output_dir}/$options{output_base}_ncRNA.fasta";
    $options{pep_fasta} = qq"$options{output_dir}/$options{output_base}.faa";
    $options{rrna_fasta} = qq"$options{output_dir}/$options{output_base}_rRNA.fasta";
    $options{all_gff} = qq"$options{output_dir}/$options{output_base}_all.gff";
    $options{gene_gff} = qq"$options{output_dir}/$options{output_base}.gff";
    $options{cds_gff} = qq"$options{output_dir}/$options{output_base}_cds.gff";
    $options{inter_gff} = qq"$options{output_dir}/$options{output_base}_inter_cds.gff";
    $options{misc_gff} = qq"$options{output_dir}/$options{output_base}_miscRNA.gff";
    $options{nc_gff} = qq"$options{output_dir}/$options{output_base}_ncRNA.gff";
    $options{rrna_gff} = qq"$options{output_dir}/$options{output_base}_rRNA.gff";

    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $handle);
    my $seq_count = 0;
    my $total_nt = 0;
    my $feature_count = 0;
    my $fasta = Bio::SeqIO->new(-file => qq">$options{fasta}",
                                -format => 'fasta', -flush => 0);
    my $gffout = Bio::Tools::GFF->new(-file => ">$options{all_gff}",
                                      -gff_version => 3);
    my $gene_gff = Bio::Tools::GFF->new(-file => ">$options{gene_gff}",
                                        -gff_version => 3);
    my $cds_gff = Bio::Tools::GFF->new(-file => ">$options{cds_gff}",
                                       -gff_version => 3);
    my $inter_gffout = Bio::Tools::GFF->new(-file => ">$options{inter_gff}",
                                            -gff_version => 3);
    my $cds_fasta = Bio::SeqIO->new(-file => qq">$options{cds_fasta}",
                                    -format => 'Fasta');
    my $pep_fasta = Bio::SeqIO->new(-file => qq">$options{pep_fasta}",
                                    -format => 'Fasta');
    my $nc_fasta = Bio::SeqIO->new(-file => qq">$options{nc_fasta}",
                                   -format => 'Fasta');
    my $misc_fasta = Bio::SeqIO->new(-file => qq">$options{misc_fasta}",
                                     -format => 'Fasta');
    my $nc_gffout = Bio::Tools::GFF->new(-file => ">$options{nc_gff}", -gff_version => 3);
    my $misc_gffout = Bio::Tools::GFF->new(-file => ">$options{misc_gff}", -gff_version => 3);
    my $rrna_gffout = Bio::Tools::GFF->new(-file => ">$options{rrna_gff}", -gff_version => 3);
    my $rrna_fasta = Bio::SeqIO->new(-file => qq">$options{rrna_fasta}", -format => 'Fasta');
    my $rrna_id_counter = {};
    my $cds_id_counter = {};
    my $ncrna_id_counter = {};
    my $miscrna_id_counter = {};
    my $id_count = 0;

    while (my $seq = $seqio->next_seq) {
        $seq_count++;
        print "Writing feature: $seq_count\n";
        $total_nt = $total_nt + $seq->length();
        $fasta->write_seq($seq);
        print "Wrote ${seq_count} features.\n";
        my @feature_list = ();
        ## defined a default name

        ## $feature is an object of type: Bio::SeqFeatureI
        ## Some of the things you can call on it are:
        ## display_name(), primary_tag(), has_tag(something), get_tag_values(), get_all_tags(), attach_seq(Bio::Seq),
        ## seq(), entire_seq(), seq_id(), gff_string(), spliced_seq(), primary_id(), phase(),
        for my $feature ($seq->top_SeqFeatures()) {
            $gffout->write_feature($feature);
            my $feat_count = 0;
        }
      FEAT: for my $feat_object ($seq->get_SeqFeatures) {
            $feature_count++;
            my $feat_count = 0;
            if ($feat_object->primary_tag eq 'rRNA') {
                $rrna_gffout->write_feature($feat_object);
                my $id_string = '';
                my $desc = '';
                my $id_hash = {};
                for my $thing (keys %{$feat_object->{_gsf_tag_hash}}) {
                    my @arr = @{$feat_object->{_gsf_tag_hash}->{$thing}};
                    $id_hash->{$thing} = $arr[0];
                }
                my $id = '';
                if (defined($id_hash->{standard_name})) {
                    $id = qq"$id_hash->{standard_name}";
                } elsif (defined($id_hash->{gene})) {
                    $id = qq"$id_hash->{gene}";
                }

                if (defined($rrna_id_counter->{id})) {
                    my $id_count = $rrna_id_counter->{id} + 1;
                    $rrna_id_counter->{id} = $id_count;
                    $id .= qq".${id_count}";
                } else {
                    $rrna_id_counter->{id} = 0;
                }

                if (defined($id_hash->{gene})) {
                    $desc .= "$id_hash->{gene} ; ";
                }

                if (defined($id_hash->{db_xref})) {
                    $desc .= "$id_hash->{db_xref} ; ";
                }

                if (defined($id_hash->{product})) {
                    $desc .= "$id_hash->{product}";
                }

                my $start = $feat_object->start;
                my $end = $feat_object->end;
                my $len = $feat_object->length;
                my $seq;
                try {
                    $seq = $feat_object->spliced_seq->seq;
                }
                catch ($e) {
                    print qq"Something went wrong getting the sequence.\n";
                    next FEAT;
                }
                my $size = {
                    id => $id_string,
                    start => $start,
                    end => $end,
                    len => $len,
                    feat => $feat_object,
                };
                $desc =~ s/\;\s*$//g;
                my $rrna_object = Bio::PrimarySeq->new(-id => $id_string, -seq => $seq,
                                                       description => $desc);
                $rrna_fasta->write_seq($rrna_object);
            } elsif ($feat_object->primary_tag eq 'CDS') {
                $cds_gff->write_feature($feat_object);
                $feat_count++;
                my $id_string = '';
                my $id_hash = {};
                for my $thing (keys %{$feat_object->{_gsf_tag_hash}}) {
                    my @arr = @{$feat_object->{_gsf_tag_hash}->{$thing}};
                    $id_hash->{$thing} = $arr[0];
                }
                my $id = '';
                if (defined($id_hash->{standard_name})) {
                    $id = qq"$id_hash->{standard_name}";
                } elsif (defined($id_hash->{protein_id})) {
                    $id = qq"$id_hash->{protein_id}";
                }

                if (defined($cds_id_counter->{id})) {
                    my $id_count = $cds_id_counter->{id} + 1;
                    $cds_id_counter->{id} = $id_count;
                    $id .= qq".${id_count}";
                } else {
                    $cds_id_counter->{id} = 0;
                }

                my $desc = '';
                if (defined($id_hash->{gene})) {
                    $desc .= "$id_hash->{gene};";
                }

                if (defined($id_hash->{db_xref})) {
                    $desc .= "$id_hash->{db_xref};";
                }

                if (defined($id_hash->{product})) {
                    $desc .= "$id_hash->{product}";
                }

                $desc .= "\n";
                my $start = $feat_object->start;
                my $end = $feat_object->end;
                my $len = $feat_object->length;
                ## This is in response to the puzzling error:
                ## "Error::throw("Bio::Root::Exception", "Location end (601574) exceeds length (0) of called sequence C"...) called at /sw/local/perl/5.28.1/perl5/Bio/Root/Root.pm line 449"
                my $seq = '';
                my $pep = '';

                my $e;
                try {
                    $seq = $feat_object->spliced_seq->seq;
                }
                catch ($e) {
                    print "Something went wrong getting the sequence.\n";
                    next FEAT;
                }
                try {
                    $pep = $feat_object->spliced_seq->translate->seq;
                }
                catch ($e) {
                    print "Something went wrong getting the peptide sequence.\n";
                    next FEAT;
                }
                my $size = {
                    id => $id,
                    start => $start,
                    end => $end,
                    len => $len,
                    feat => $feat_object,
                };
                push(@feature_list, $size);
                $desc =~ s/\;\s*$//g;
                my $seq_object = Bio::PrimarySeq->new(-id => $id, -seq => $seq, description => $desc);
                my $pep_object = Bio::PrimarySeq->new(-id => $id, -seq => $pep, description => $desc);
                $cds_fasta->write_seq($seq_object);
                my $ttseq = $seq_object->seq;
                $pep_fasta->write_seq($pep_object);
            } elsif ($feat_object->primary_tag eq 'ncRNA') {
                $nc_gffout->write_feature($feat_object);
                $feat_count++;
                my $id_string = '';
                my $id_hash = {};
                for my $thing (keys %{$feat_object->{_gsf_tag_hash}}) {
                    my @arr = @{$feat_object->{_gsf_tag_hash}->{$thing}};
                    $id_hash->{$thing} = $arr[0];
                }
                my $id = '';
                if (defined($id_hash->{standard_name})) {
                    $id = $id_hash->{standard_name};
                } elsif (defined($id_hash->{gene})) {
                    $id = $id_hash->{gene};
                }

                if (defined($ncrna_id_counter->{$id})) {
                    my $id_count = $ncrna_id_counter->{$id} + 1;
                    $ncrna_id_counter->{$id} = $id_count;
                    $id .= ".${id_count}";
                } else {
                    $ncrna_id_counter->{$id} = 0;
                }

                my $desc = '';
                if (defined($id_hash->{db_xref})) {
                    $desc .= "$id_hash->{db_xref};";
                }

                if (defined($id_hash->{product})) {
                    $desc .= "$id_hash->{product}";
                }
                $desc .= "\n";

                my $start = $feat_object->start;
                my $end = $feat_object->end;
                my $len = $feat_object->length;
                ## This is in response to the puzzling error:
                ## "Error::throw("Bio::Root::Exception", "Location end (601574) exceeds length (0) of called sequence C"...) called at /sw/local/perl/5.28.1/perl5/Bio/Root/Root.pm line 449"
                my $seq = '';
                my $e;
                try {
                    $seq = $feat_object->spliced_seq->seq;
                }
                catch ($e) {
                    print "Something went wrong getting the sequence.\n";
                    next FEAT;
                }
                my $size = {
                    id => $id,
                    start => $start,
                    end => $end,
                    len => $len,
                    feat => $feat_object,
                };
                push(@feature_list, $size);
                $desc =~ s/\;\s*$//g;
                my $seq_object = Bio::PrimarySeq->new(-id => $id, -seq => $seq, description => $desc);
                $nc_fasta->write_seq($seq_object);
                my $ttseq = $seq_object->seq;
            } elsif ($feat_object->primary_tag eq 'misc_RNA') {
                ## Note: misc_RNA appears to me to be an ensembl-specific type and includes the various
                ## feature types which fall under ncRNA at NCBI...
                ## This is noteworthy, because ensembl's genbank files keep the ID under the tag
                ## 'standard_name' rather than ID or gene or whatever NCBI does.
                $misc_gffout->write_feature($feat_object);
                $feat_count++;
                my $id_string = '';
                my $id_hash = {};
                ## I should convert this to use $feat_object->get_tags and $feat_object->get_tag_values
                for my $thing (keys %{$feat_object->{_gsf_tag_hash}}) {
                    my @arr = @{$feat_object->{_gsf_tag_hash}->{$thing}};
                    $id_hash->{$thing} = $arr[0];
                }
                my $id = '';
                if (defined($id_hash->{standard_name})) {
                    $id = $id_hash->{standard_name};
                } elsif (defined($id_hash->{gene})) {
                    $id = $id_hash->{gene};
                }

                if (defined($miscrna_id_counter->{$id})) {
                    my $id_count = $miscrna_id_counter->{$id} + 1;
                    $miscrna_id_counter->{$id} = $id_count;
                    $id .= ".${id_count}";
                } else {
                    $miscrna_id_counter->{$id} = 0;
                }
                my $desc = '';

                if (defined($id_hash->{db_xref})) {
                    $desc .= "$id_hash->{db_xref};";
                }

                if (defined($id_hash->{product})) {
                    $desc .= "$id_hash->{product}";
                }
                $desc .= "\n";

                my $start = $feat_object->start;
                my $end = $feat_object->end;
                my $len = $feat_object->length;
                ## This is in response to the puzzling error:
                ## "Error::throw("Bio::Root::Exception", "Location end (601574) exceeds length (0) of called sequence C"...) called at /sw/local/perl/5.28.1/perl5/Bio/Root/Root.pm line 449"
                my $seq = '';
                my $e;
                try {
                    $seq = $feat_object->spliced_seq->seq;
                }
                catch ($e) {
                    print "Something went wrong getting the sequence.\n";
                    next FEAT;
                }
                my $size = {
                    id => $id,
                    start => $start,
                    end => $end,
                    len => $len,
                    feat => $feat_object,
                };
                push(@feature_list, $size);
                $desc =~ s/\;\s*$//g;
                my $seq_object = Bio::PrimarySeq->new(-id => $id, -seq => $seq, description => $desc);
                $misc_fasta->write_seq($seq_object);
                ## Why do I have this in my previous block?
                ##my $ttseq = $seq_object->seq;
            }   ## End looking for miscRNA
        } ## End looking at every feature and putting them into the @feature_list
        my %inter_features = ();
        my ($p_st, $p_en, $c_st, $c_en, $inter_start, $inter_end) = 0;
      INTER: for my $c (0 .. $#feature_list) {
            my $p_size = {};
            if ($c == 0) {
                $p_size = $feature_list[$#feature_list];
            } else {
                $p_size = $feature_list[$c - 1];
            }

            my $c_size = $feature_list[$c];
            next INTER if (!defined($p_size->{start}));
            next INTER if (!defined($p_size->{end}));
            next INTER if (!defined($c_size->{start}));
            next INTER if (!defined($c_size->{end}));

            my $p_start = $p_size->{start};
            my $p_end = $p_size->{end};
            my $c_start = $c_size->{start};
            my $c_end = $c_size->{end};
            my $interstart = 0;
            my $interend = 0;

            if ($p_start > $p_end) {
                $interstart = $p_start;
            } else {
                $interstart = $p_end;
            }

            if ($c_start < $c_end) {
                $interend = $c_start;
            } else {
                $interend = $c_end;
            }

            my $tmp_start = $interstart;
            my $tmp_end = $interend;
            if ($tmp_start > $tmp_end) {
                $interstart = $tmp_end;
                $interend = $tmp_start;
            }

            ## I need a little logic to catch the first feature
            ## As it is currently written, it grabs the position of the last feature
            ## and uses that number rather than the beginning of the chromosome.
            ## Therefore we end up with a first feature which contains the entire chromosome.
            if ($c <= 2) { ## Arbitrarily set it to check the first 2 features
                if ($interend >= 100000) { ## If the first feature take >= 1/2 the chromosome
                    $interend = $interstart;
                    $interstart = 1;
                }
            }

            my $new_feature = $c_size->{feat};
            my $inter_location = Bio::Location::Atomic->new(-start => $interstart,
                                                            -end => $interend,
                                                            -strand => 1);
            $new_feature->{_location} = $inter_location;
            $inter_gffout->write_feature($new_feature);
        }
    } ## End while every sequence

    my $ret_stats = {
        num_sequences => $seq_count,
        total_nt => $total_nt,
        num_features => $feature_count,
        fasta_out => $options{fasta},
        all_gff_out => $options{all_gff},
        gene_gff => $options{gene_gff},
        inter_cds => $options{inter_gff},
        cds_fasta => $options{cds_fasta},
        cds_gff => $options{cds_gff},
        pep_fasta => $options{pep_fasta},
        rrna_gff => $options{rrna_gff},
        rrna_fasta => $options{rrna_fasta},
        ncrna_gff => $options{nc_gff},
        ncrna_fasta => $options{nc_fasta},
        misc_gff => $options{misc_gff},
        misc_fasta => $options{misc_fasta},
    };
    close($handle);
    $fasta->close();
    $gffout->close();
    $gene_gff->close();
    $cds_gff->close();
    $inter_gffout->close();
    $cds_fasta->close();
    $pep_fasta->close();
    $rrna_gffout->close();
    $rrna_fasta->close();
    $nc_gffout->close();
    $nc_fasta->close();
    return($ret_stats);
}

=head2 C<Gff2Fasta>

 It writes one file of amino acid sequence and one of nucleotides.
 Upon completion, it returns the number of entries written.

=over

=item C<Arguments>

 species: Add the .gff and .fasta suffixes to this and look in $libdir.
 gff(required): Gff file to read.
 input(required): Genome fasta file to read.
 gff_tag('gene_id'): ID tag to use when extracting features.
 gff_type(undef): gff type to extract, left undefined it will count
  them up and choose the most highly represented.

=item C<Invocation>

> cyoa --method gff2 --input lmajor.fasta \
   --gff lmajor.gff --gff_type CDS --gff_tag ID

=back

=cut
sub Gff2Fasta {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        gff_tag => 'gene_id',
        gff_type => 'cds',);
    unless ($options->{species} || ($options->{input} && $options->{gff})) {
        die("This requires either a species or an input gff/fasta.")
    }
    my $genome;
    my $gff;
    my $species;
    if ($options->{species}) {
        $genome = $paths->{fasta};
        $gff = $paths->{gff};
    } elsif ($options->{input} && $options->{gff}) {
        $genome = $options->{input};
        $gff = $options->{gff};
    }

    $species = basename($genome, @{$class->{suffix_array}});
    $species = basename($species, @{$class->{suffix_array}});
    ## Count IDs as we come across them.
    my $id_counter = {};
    my $nt_id_counter = {};
    my $contig_id_counter = {};
    my $aa_id_counter = {};

    my $wanted_tag = $options->{gff_tag};
    my $wanted_type = $options->{gff_type};
    print "Reading: ${genome} to extract entries of type: ${wanted_type} by tag ${wanted_tag}.\n";
    my $chromosomes = $class->Read_Genome_Fasta(genome => $genome);
    my $gff_handle = Bio::Adventure::Get_FH(input => $gff);
    my $nt_out = Bio::SeqIO->new(
        -format => 'Fasta',
        -file => qq">$outputs->{nt}",);
    my $aa_out = Bio::SeqIO->new(
        -format => 'Fasta',
        -file => qq">$outputs->{aa}",);
    my $gff_out = FileHandle->new(">$outputs->{gff}");
    ## Note that this and the next line might not be a good idea,
    ## HT_Types only looks at the first n (40,000) records and uses that as a heuristic
    ## to see that the wanted type is actually in the gff file.
    print "Counting ${wanted_type} entries in ${gff}.\n";
    die("Missing the gff file: ${gff}.") unless (-r $gff);
    my $feature_type = $class->Bio::Adventure::Count::HT_Types(
        annotation => $gff,
        feature_type => $wanted_type,);
    $feature_type = $feature_type->[0];
    my $annotation_in = Bio::Tools::GFF->new(-fh => $gff_handle, -gff_version => 3);
    my $features_written = 0;
    my $features_read = 0;
    print "Iterating over gff entries in search of type: ${wanted_type}.\n";
  LOOP: while (my $feature = $annotation_in->next_feature()) {
        $features_read++;
        my $primary_type = $feature->primary_tag();
        unless ($primary_type eq $wanted_type) {
            next LOOP;
        }
        my $start = $feature->start();
        my $end = $feature->end();
        my $strand = $feature->strand();
        my $gff_chr = $feature->seq_id();
        my @ids;
        my $e;
        try {
            @ids = $feature->each_tag_value($wanted_tag);
        }
        catch ($e) {
            print "Did not find the tag: ${wanted_tag}, perhaps check the gff file.\n";
            next LOOP;
        }
        my $gff_string = $annotation_in->gff_string($feature);
        print $gff_out qq"${gff_string}\n";
        my $id = $ids[0];
        my $genome_obj = $chromosomes->{$gff_chr}->{obj};
        my $cds = $genome_obj->trunc($start, $end);
        if ($strand == -1) {
            $cds = $cds->revcom;
        }
        my $id_string = $id;
        my $long_id_string = qq"chr_${gff_chr}_id_${id}_start_${start}_end_${end}";

        if (!defined($id_counter->{$id})) {
            $id_counter->{$id} = 1;
        } else {
            $id_counter->{$id}++;
            $id_string = $long_id_string;
        }

        if (!defined($id_string)) {
            $id_string = $long_id_string;
        }

        my $nt_sequence = Bio::Seq->new(
            -id => $id_string,
            -seq => $cds->seq());
        $nt_out->write_seq($nt_sequence);
        my $aa_sequence = Bio::Seq->new(
            -id => $id_string,
            -seq => $cds->translate->seq());
        $aa_out->write_seq($aa_sequence);
        $features_written++;
    }                           ## End LOOP
    $gff_handle->close();
    print "Wrote ${features_written} features to the aa/nt files.\n";
    my $retlist = {
        features_written => $features_written,
        outputs => $outputs,
        output => $outputs->{nt},
        output_nt => $outputs->{nt},
        output_aa => $outputs->{aa},
        output_gff => $outputs->{gff},
    };
    return($retlist);
}

=head2 C<Read_GFF>

 use Bio::Tools::GFF to gather information from a gff file.

 This pulls information of potential interest from a gff file.  In
 theory it could be trivially modified to be robust to different
 formats and other shenanigans.

=cut
sub Read_GFF {
    my ($class, %args) = @_;
    my $chromosomes = $args{chromosomes};
    my $gff = FileHandle->new("<$args{gff}");
    my $annotation_in = Bio::Tools::GFF->new(-fh => \$gff, -gff_version => 3);
    my $gff_out = {};
    print "Starting to read gff: $args{gff}\n";
  LOOP: while(my $feature = $annotation_in->next_feature()) {
        next LOOP unless ($feature->{_primary_tag} eq $args{gff_type});
        my $location = $feature->{_location};
        my $start = $location->start();
        my $end = $location->end();
        my $strand = $location->strand();
        my @ids = $feature->each_tag_value('ID');
        my $id = "";
        my $gff_chr = $feature->{_gsf_seq_id};
        my $gff_string = $annotation_in->gff_string($feature);
        if (!defined($chromosomes->{$gff_chr})) {
            print STDERR "Something is wrong with $gff_chr\n";
            next LOOP;
        }

        for my $i (@ids) {
            $i =~ s/^cds_//g;
            $i =~ s/\-\d+$//g;
            $id .= "${i} ";
        }
        $id =~ s/\s+$//g;
        my @gff_information = split(/\t+/, $gff_string);
        my $description_string = $gff_information[8];
        my $orf_chromosome = $gff_chr;
        my $annot = {
            id => $id,
            start => $start,  ## Genomic coordinate of the start codon
            end => $end,      ## And stop codon
            strand => $strand,
            description_string => $description_string,
            chromosome => $gff_chr,
        };
        $gff_out->{$gff_chr}->{$id} = $annot;
    }                     ## End looking at every gene in the gff file
    $gff->close();
    return($gff_out);
}

=head2 C<Sam2Bam>

 Sort, compress, index a sam file into a usable bam file.

 Used to invoke samtools to take the sam output from bowtie/bwa and
 convert it to an compressed-sorted-indexed bam alignment. This
 function just calls $class->Samtools(), but has a little logic
 to see if the invocation of this is for an extent .sam file or
 calling it on an existing .fastq(.gz), in which case one must
 assume it is being called on one or more files in the bowtie_out/
 directory and which start with $basename, include something like
 -trimmed-v0M1.sam.

=over

=item C<Arguments>

 input(required): Input sam file.
 species(required): Species name to hunt down the appropriate
  fasta/gff files.
 modules('samtools', 'bamtools'): Environment modules to load.

=item C<Invocation>

> cyoa --task convert --method sam2bam --input bob.sam --species hg38_91

=back

=cut
sub Sam2Bam {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],);
    my @input_list = ();
    if ($options->{input}) {
        push(@input_list, $options->{input});
    } elsif (-r $options->{input} and $options->{input} =~ /\.sam$/) {
        push(@input_list, $options->{input});
    } elsif (-r $options->{input} and $options->{input} =~ /\.fastq$/) {
        if (-d "bowtie_out") {
            find({ wanted => sub { push(@input_list, "bowtie_out/$_") if ($_ =~ /\.sam$/); }, follow => 1 }, 'bowtie_out/');

        } else {
            foreach my $k (%{$options->{bt_args}}) {
                my $output_string = qq"bowtie_out/$paths->{jbasename}-${k}.sam";
                push(@input_list, $output_string);
            }
            my $bt = $class->Bio::Adventure::Map::Bowtie(%args);
        }
    } else {
        die("I don't know what to do without a .fastq file nor a .sam file.\n");
    }

    my $loaded = $class->Module_Load(modules => $options->{modules});
    my $sam = $class->Bio::Adventure::Convert::Samtools(%args,
                                                        jdepends => $options->{jdepends},
                                                        sam => \@input_list);
    return($sam);
}

=head2 C<Samtools>

 Does the actual work of converting a sam to a bam file.

 $hpgl->Samtools() calls (in order): samtools view, samtools sort,
 and samtools index.  Upon completion, it invokes bamtools stats to
 see what the alignments looked like.

 It explicitly does not pipe one samtools invocation into the next,
 not for any real reason but because when I first wrote it, it
 seemed like the sorting was taking too long if I did not already
 have the alignments in a bam file.

=over

=item C<Arguments>

 input(required): Input sam file.
 species(required): Input species prefix for finding a genome.

=back

=cut
sub Samtools {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    my $ready;
    if (!$options{jdepends}) {
        $ready = $class->Check_Input(files => $inputs);
    }
    $options{samtools_args} = $class->Passthrough_Args(arbitrary => $options{arbitrary});

    ## There are a few options in my samtools template to check:
    ## sam_input which may just be input, main (output),
    ## paired and output_paired,
    ## fasta (the reference genome and should come from gfv),
    ## nosuffix (some version of this should come from gfv),
    ## nomismatch (defaults to 0) and output_sorted_nomismatch,
    ## mapped (defaults to 0) and output_mapped (fastq),
    ## unmapped (defaults to 0) and output_unmapped (fastq)
    $options{sam_input} = $inputs->[0]->{fullpath};
    my $sam_name = $inputs->[0]->{basename};
    my $sam_dir = $inputs->[0]->{directory};
    my $sam_ext = $inputs->[0]->{filebase_extension};
    $options{main} = qq"${sam_dir}/${sam_name}.bam";
    $options{output_sorted} = qq"${sam_dir}/${sam_name}_sorted.bam";
    $options{output_paired} = qq"${sam_dir}/${sam_name}_paired.bam";
    $options{output_nosuffix} = qq"${sam_dir}/${sam_name}";
    $options{output_nomismatch} = qq"${sam_dir}/${sam_name}_nomismatch.bam";
    $options{output_mismatch} = qq"${sam_dir}/${sam_name}_mismatch.bam";
    if ($options{samtools_mapped}) {
        $options{output_mapped} = qq"${sam_dir}/${sam_name}_samtools_mapped.fastq";
    }

    if ($options{samtools_unmapped}) {
        $options{output_unmapped} = qq"${sam_dir}/${sam_name}_samtools_unmapped.fastq";
    }

    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    print "TESTME: Samtools about to submit with main: $options{main}\n";
    my $samtools_job = $class->Submit();
    print "After submission main: $samtools_job->{main}\n";
    return($samtools_job);
}

=head2 C<Species2SF>

Given a species identifier, read the gff and fasta files and produce a SeqFeatures database from it.

=cut
sub Species2SF {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['species'],);
    ## Read the GFF
    print "Reading $paths->{gff}\n";
    my $gff_fh = Bio::Adventure::Get_FH(input => $paths->{gff}, suffix => '| grep -v "^#"');
    my $sf_by_contig = {};
    my $input_gff = new Bio::FeatureIO(-format => 'GFF', -version => 3, -fh => $gff_fh);
  INFEAT: while (my $feature = $input_gff->next_feature()) {
        my @sf_tmp_features = ();
        my $contig_id = $feature->seq_id;
        if (defined($sf_by_contig->{$contig_id})) {
            @sf_tmp_features = @{$sf_by_contig->{$contig_id}};
        } else {
            $sf_by_contig->{$contig_id} = ();
        }

        push(@sf_tmp_features, $feature);
        $sf_by_contig->{$contig_id} = \@sf_tmp_features;
    }
    $gff_fh->close();
    my @contig_ids = sort keys %{$sf_by_contig};
    my $num_contigs = scalar(@contig_ids);

    ## Read the genome
    print "Reading genome: $paths->{fasta}\n";
    my $input_genome = Bio::SeqIO->new(-file => $paths->{fasta}, -format => 'Fasta');
    my $sf_sequences = {};
    my $contig_features = {};
  GENOME: while (my $in_seq = $input_genome->next_seq()) {
        ## Create a feature containing the entire contig.
        my $id = $in_seq->id;
        my $contig_feature = new Bio::SeqFeature::Generic(
            -start => 1,
            -end => $in_seq->length,
            -id => $id,
            -strand => 0,
            -primary => 'contig',
            -display_name => $in_seq->display_name,);
        my $attached = $in_seq->add_SeqFeature($contig_feature);
        $contig_features->{$id} = $contig_feature;
        my @in_contig_features = @{$sf_by_contig->{$id}};
        my $num_gff_features = scalar(@in_contig_features);
        for my $in_feature (@in_contig_features) {
            my $in_gff_attached = $in_seq->add_SeqFeature($in_feature);
        }
        $sf_sequences->{$id} = $in_seq;
    }
    my $ret = {
        contigs => $contig_features,
        sf_sequences => $sf_sequences,
        sf_by_contig => $sf_by_contig, };
    return($ret);
}

=head2 C<Split_Fasta>

Split an input fasta file into individual single pieces.

=cut
sub Split_Fasta {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input'],);
    my $fh = Bio::Adventure::Get_FH(input => $options->{input});
    my $input = Bio::SeqIO->new(-fh => $fh, -format => 'Fasta');
    my $ret = {};
    while (my $sequence = $input->next_seq()) {
        my $id = $sequence->id;
        $id =~ s/[[:punct:]]/_/g;
        my $seq_string = $sequence->seq;
        my $out_filename = qq"${id}.fasta";
        my $writer = Bio::SeqIO->new(-format => 'fasta', -file => qq">${out_filename}");
        $writer->write_seq($sequence);
        $ret->{$id} = $sequence;
    }
    return($ret);
}

=head1 AUTHOR - atb

Email abelew@gmail.com

=head1 SEE ALSO

L<samtools> L<Bio::FeatureIO> L<Bio::Tools::GFF> L<Bio::Seq> L<Bio::SeqIO>

=cut

1;
