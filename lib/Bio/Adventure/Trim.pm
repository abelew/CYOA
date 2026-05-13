package Bio::Adventure::Trim;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::Adventure::Config;
use Cwd qw"abs_path getcwd cwd";
use Data::Printer;
use File::Basename;
use File::Path qw"make_path";
use File::ShareDir qw":ALL";
use File::Spec;
use File::Temp qw":POSIX";
use File::Which qw"which";
use FileHandle;
use Bio::SeqIO;

=head1 NAME

 Bio::Adventure::Trim - Use trimomatic/cutadapt/etc to trim libraries

=head1 SYNOPSIS

 This file is responsible for invoking the various sequence trimmers.

=head1 METHODS

=head2 C<Cogent>

 Invoke takara's cogent tool to remove UMIs from a sequencing run.

=cut
sub Cogent {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        modules => ['cogent'],
        required => ['input', 'species', 'input_umi'],
        type => 'Stranded_UMI',);
    my $cwd_name = basename(cwd());
    my @input_filenames = split(/$options->{delimiter}/, $options->{input});
    my $stdout = qq"$paths->{output_dir}/cogent.stdout";
    my $stderr = qq"$paths->{output_dir}/cogent.stderr";
    my $jstring = qq!cogent demux \\
  -i $input_filenames[0] \\
  -p $input_filenames[1] \\
  -b $options->{input_umi} \\
  -t $options->{type} \\
  -o $paths->{output_dir}/demux \\
  -n $options->{jcpu} \\
  2>${stderr} 1>${stdout}
test=\$?
if [[ "\${test}" -eq "0" ]]; then
  echo "cogent demux succeeded."
else
  echo "cogent demux failed."
  exit \$?
fi
cogent analyze \\
  -i $paths->{output_dir}/demux \\
  -g $options->{species} \\
  -o $paths->{output_dir}/analyze \\
  -t $options->{type} \\
  --threads $options->{jcpu} \\
  2>>${stderr} 1>>${stdout}
test=\$?
if [[ "\${test}" -eq "0" ]]; then
  echo "cogent analyze succeeded."
else
  echo "cogent analyze failed."
  exit \$?
fi
!;
    my $cogent = $class->Submit(
        jdepends => $options->{jdepends},
        jstring => $jstring,
        outputs => $outputs,
        stdout => $stdout,
        stderr => $stderr,);
    return($cogent);
}

=head2 C<Cutadapt>

 Invoke cutadapt on a pile of sequence.
 10.14806/ej.17.1.200

 Use biopieces/cutadapt to attempt to remove sequence adapters from a library.
 This is most common used by me for ribosome profiling libraries.

=cut
sub Cutadapt {
    my ($class, %args) = @_;
    ## Once I get all this cleaned up, I would like to create %template_args
    ## right at the beginning of each function and fill it in as we go;
    ## but until I am comfortable with the workflow let us add the variables
    ## right before Update_Spec() so that I can poke and play easily.
    ## Ideally, I will just make a copy of class->current_state as %template_args
    ## and modify the slots of interest.
    print "TESTME Starting Cutadapt current_state->cluster: $class->{current_state}->{cluster}\n";
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    ## Dereference this so that we may change them without messing with the underlying
    ## function information.
    my %options = %{$function_info};
    print "Starting cutadapt with cluster: $options{cluster}\n";
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];

    my $extra_args = $class->Passthrough_Args(arbitrary => $options{arbitrary});
    my $num_inputs = scalar(@{$inputs});

    $options{adapter_string} = q{};
    if ($options{any_adapter}) {
        $options{adapter_string} .= qq"-a $options{any_adapter} ";
    }
    if ($options{left_adapter}) {
        $options{adapter_string} .= qq"-b $options{left_adapter} ";
    }
    if ($options{right_adapter}) {
        $options{adapter_string} .= qq"-g $options{right_adapter} ";
    }
    if ($options{type} eq 'old_tnseq') {
        $options{adapter_string} = $options{old_tnseq_adapter};
        $options{maxlength} = 20;
    } elsif ($options{type} eq 'old_riboseq') {
        $options{adapter_string} = $options{old_riboseq_adapter};
        $options{minlength} = 16;
        $options{maxlength} = 42;
    }

    $options{min_flag} = q{};
    if ($options{minlength}) {
        $options{min_flag} = qq"-m $options{minlength} --too-short-output=$options{too_short_r1}";
    }
    $options{max_flag} = q{};
    if ($options{maxlength}) {
        $options{max_flag} = qq"-M $options{maxlength} --too-long-output=$options{too_long_r1}";
    }
    ## The maximum number of times an adapter may be removed.
    unless ($options{max_removed}) {
        $options{max_removed} = 3;
    }
    ## Separate trimmed/untrimmed reads into two piles?
    ## output_r1 is in functions.json under Cutadapt->paths
    $options{compress_string} = qq"xz -9e -f $options{output_r1}\n";
    $options{untrimmed_string} = q{};
    if ($options{separate}) {
        $options{untrimmed_string} = qq"--untrimmed-output=$options{untrimmed_r1}";
        $options{compress_string} .= qq"xz -9e -f $options{untrmmed_r1}\n";
    }

    ## Define input/output strings here, then fill in the exceptions.
    $options{input_string} = q{};
    my $output = qq"$options{output_r1}";
    $options{output_string} = q{};

    $options{paired} = 0;
    if ($num_inputs == 1) {
        $options{input_string} = Bio::Adventure::Get_FD(
            input => $inputs->[0]->{fullpath},
            ignored => '.gz');
        my $r1_base = $inputs->[0]->{filebase_extension};
        if ($options{compress}) {
            $output .= '.xz';
            $options{compress_string} = qq"
xz -9e -f $options{output_r1}
";
        }
        $options{output_string} = qq"-o $options{output_r1}";
    } elsif ($num_inputs == 2) {
        $options{paired} = 1;
        my $r1_fd = $class->Get_FD(input => $inputs->[0]->{fullpath});
        my $r2_fd = $class->Get_FD(input => $inputs->[1]->{fullpath});
        $options{input_string} = qq"${r1_fd} ${r2_fd}";
        $options{output_string} = qq"-o $options{output_r1} -p $options{output_r2}";
        if ($options{separate}) {
            $options{untrimmed_string} .= qq" --untrimmed-paired-output=$options{untrimmed_r2}";
            $options{compress_string} .= qq"xz -9e -f $options{untrimmed_r2}\n";
        }
        if ($options{minlength}) {
            $options{min_flag} .= qq" --too-short-paired-output=$options{too_short_r2}";
        }
        if ($options{maxlength}) {
            $options{max_flag} .= qq" --too-long-paired-output=$options{too_long_r2}";
        }
        $options{compress_string} .= qq"xz -9e -f  $options{output_r2}\n";
    } else {
        croak("This requires either 1 or 2 inputs, not ${num_inputs}.\n");
    }

    if ($options{maxerr}) {
        $options{max_err_string} = qq"-e $options{maxerr}";
    }
    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    my $cutadapt = $class->Submit();
    my $stats = $class->Bio::Adventure::Metadata::Cutadapt_Stats(
        input => $options{stdout},
        input_name => $options{output_base},
        jdepends => $cutadapt->{job_id},
        paired => $options{paired});
    $cutadapt->{stats} = $stats;
    return($cutadapt);
}

=head2 C<Fastp>

 Invoke fastp on a sequence dataset.
 10.1093/bioinformatics/bty560

 Fastp is an excellent trimomatic alternative.  A version of this invocation was
 taken from one of the kmcp tutorials and provides flags suitable for short read
 removal, adapter trimming, and the html report:

 https://bioinf.shenwei.me/kmcp/tutorial/profiling/

 Here are some of the options which may be passed along:
 * -i/-o/-I/-O: input/output r1/R2
 * -m merge R1/R2
 * -6 use phred64 instead of 33
 * -A/--disable-adapter-trimming
 * -a/--adapter_sequence R1 adapter provided
 * --adapter_sequence_r2 R2 adapter provided
 * --adapter_fasta fasta file containing adapters
 * --detect_adapter_for_pe Turn on paired end auto-adapter-detection
 * -f/--trim_front1/-F/--trim_front2  # bases to trim
 * -t/--trim_tail1/-T/--trim_tail2   # bases to trim
 * -b/--max_len1/-B/--max_len2     Maximum length
 * -D/--dedup Deduplication
 * --dup_calc_accuracy accuracy level to calculate duplicates
 * --dont_eval_duplication
 * -g/--trim_poly_g/--poly_g_min_len/-G/--disable_trim_poly_g  various polyG options
 * -x/--trim_poly_x/--poly_x_min_len various polyX
 * -5/--cut_front sliding window quality from front
 * -3/--cut_tail from back
 * -r/--cut_right
 * -W/--cut_window_size
 * -M/--cut_mean_quality
 * --cut_front(tail)_window_size/--cut_front(tail)_mean_quality/--cut_front_mean_quality
 * -Q/--disable_quality_filtering
 * -q/--qualified_quality_phred (15)
 * -u/--unqualified_percent_limit %bases < threshold to drop read
 * -n/--n_base_limit max N
 * -e/--average_qual
 * -L/--disable_length_filtering
 * -l/--length_required/--length_limit
 * -y/--low_complexity_filter
 * -Y/--complexity_threshold
 * --filter_by_index1/--filter_by_index2 1 barcode/line to filter out
 * --filter_by_index_threshold
 * -c/--correction base correction via overlaps
 * --overlap_diff_limit/--overlap_diff_percent_limit Thresholds for correction
 * -U/--umi Enable molecular identifier analysis
 * --umi_loc/--umi_len/--umi_prefix/--umi_skip
 * -p/--overrepresentation_analysis Enable over representation analysis
 * -P/--overrepresentation_sampling Sample 1/x reads
 * -j/--json json report
 * -h/--htlp html report
 * -R/--report_title Report title
 * -w/--thread # threads
 * -s/--split/-S/--split_by_lines/-d/--spit_prefix_digits Splitting output file stuff

One note: it appears that fastp is unable to successfully read from a bash subshell: <()

=cut
sub Fastp {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];
    my $extra_args = $class->Passthrough_Args(arbitrary => $options{arbitrary});
    $extra_args .= ' -D ' if ($options{deduplicate});
    $extra_args .= ' -c ' if ($options{correction});
    $extra_args .= ' -y ' if ($options{complexity});
    my $input_string = q{};
    my $num_inputs = scalar(@{$inputs});
    my $output = q{};
    my $compress_string = q{};
    $options{paired} = 0;
    $options{r1} = $inputs->[0]->{fullpath};
    my $r1 = $inputs->[0]->{fullpath};
    my $r1_base = $inputs->[0]->{filebase_extension};
    $options{r1} = Bio::Adventure::Get_FD(input => $inputs->[0]->{fullpath},
                                          ignored => '.gz');
    $options{input_r1_string} = qq"-i $options{r1}";
    $options{output_r1_string} = qq"-o $options{output_r1}";
    if ($options{compress}) {
        $output = qq"$options{output_r1}.xz";
        $compress_string = qq"
xz -9e -f $options{output_r1}
";
    }
    if ($num_inputs == 2) {
        $options{paired} = 1;
        my $r2 = $inputs->[1]->{fullpath};
        $options{r2} = Bio::Adventure::Get_FD(input => $inputs->[1]->{fullpath},
                                              ignored => '.gz');
        my $r2_base = $inputs->[1]->{filebase_extension};
        if ($options{compress}) {
            $options{main} = qq"$options{output_r1}.xz:$options{output_r2}.xz";
            $compress_string = qq"
xz -9e -f $options{output_r1}
xz -9e -f $options{output_r2}
";
        }
        $options{input_r2_string} = qq"-I $options{r2}";
        $options{output_r2_string} = qq"-O $options{output_r2}";
    } elsif ($num_inputs > 2) {
        die("An unusual number of inputs was provided: ${num_inputs}.");
    }
    my $umi_flags = '';
    if ($options{umi}) {
        $umi_flags = '-U';
    }
    my $report_flags = qq"-h $options{report}";
    my $json_flags = qq"-j $options{json}";

    $options{report_flags} = $report_flags;
    $options{json_flags} = $json_flags;
    $options{umi_flags} = $umi_flags;
    $options{input_string} = $input_string;
    $options{extra_args} = $extra_args;

    ## Set up the template here.
    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    my $fastp = $class->Submit();
    my $stats = $class->Bio::Adventure::Metadata::Fastp_Stats(
        input => $options{stderr},
        input_name => $options{output_base},
        jdepends => $fastp->{job_id},
        paired => $options{paired},);
    $fastp->{stats} = $stats;
    return($fastp);
}

sub PolyA_Extractor {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        min_A => 10,
        min_seq => 10,
        output => qq"outputs/10polyA_extracted/polyA_reads.fastq",
        required => ['input', ],);
    my $output_dir;
    my $final_output;
    if (defined($options->{output})) {
        $output_dir = dirname($options->{output});
        $final_output = $options->{output};
    } else {
        $output_dir = qq"outputs/$options->{jprefix}polyA_extracted";
        $final_output = qq"${output_dir}/poyA_reads.fastq";
    }
    my @inputs = ();
    if ($options->{input} =~ /$options->{delimiter}/) {
        @inputs = split(/$options->{delimiter}/, $options->{input});
    } else {
        push(@inputs, $options->{input});
    }
    my $stdout = qq"${output_dir}/polyA_reads.stdout";
    my $stderr = qq"${output_dir}/polyA_reads.stderr";
    my $comment = qq"## Extract polyA reads from a sequencing file. This currently assumes an unstranded library.";
    my $polyA;
    for my $input_file (@inputs) {
        my $jname = qq"polyAextract_$options->{job_basename}";
        my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Trim;
my \$result = \$h->Bio::Adventure::Trim::PolyA_Extractor_Worker(
  input => '${input_file}',
  jname => '${jname}',
  jprefix => '$options->{jprefix}',
  min_A => '$options->{min_A}',
  min_seq => '$options->{min_seq}',
  output => '${final_output}',
  output_dir => '${output_dir}',);
!;
        $polyA = $class->Submit(
            input => $input_file,
            jdepends => $options->{jdepends},
            jstring => $jstring,
            language => 'perl',
            min_A => $options->{min_A},
            min_seq => $options->{min_seq},
            output => $final_output,
            output_dir => $output_dir,
            outputs => $outputs,
            stdout => $stdout,
            stderr => $stderr,);
        $options->{jdepends} = $polyA->{job_id};
    }
    return($polyA);
}

sub PolyA_Extractor_Worker {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        length => 1000000,
        min_A => 10,
        min_seq => 10,
        required => ['input', ],);
    my $out_fh = FileHandle->new(">>$options->{output}");
    my $out_seq = Bio::SeqIO->new(-fh => $out_fh, -format => 'Fastq');
    my $log_dir = dirname($options->{stdout});
    my $log_name = basename($options->{stdout}, ('.stdout'));
    my $log = FileHandle->new(">>${log_dir}/${log_name}.log");
    my $reads_counted = 0;
    my $reads_written = 0;
    print $log "Starting to process: $options->{input}\n";
    my $fh = $class->Get_FH(input => $options->{input});
    my $in_seq = Bio::SeqIO->new(-fh => $fh, -format => 'Fastq');
    my %info = (
        reads_counted => 0,
        reads_written => 0,);
  FQ: while (my $seq = $in_seq->next_seq) {
        $info{reads_counted}++;
        ##print $log "$info{reads_counted}\n";
        my $id = $seq->id;
        my $sequence = $seq->seq;
        ##print $log "SEQ: $sequence\n";
        my $len = $seq->length;
        my $match_begin = undef;
        my $new_seq = '';
        if ($sequence =~ /A{$options->{min_A},}$/) {
            $match_begin = $-[0];
        } elsif ($sequence =~ /^T{$options->{min_A},}/) {
            $seq = $seq->revcom;
            $sequence = $seq->seq;
            $match_begin = $len - $+[0];
        } else {
            ##print $log "No match.\n";
            next FQ;
        }
        if ($match_begin < $options->{min_seq}) {
            ## Too much polyA
            ##print $log "Matched but Too short.\n";
            next FQ;
        }
        $new_seq = $seq->trunc(1, $match_begin);
        $info{reads_written}++;
        ##print $log "Writing a sequence.\n";
        my $written = $out_seq->write_seq($new_seq);
    } ## End iterating over every sequence.
    ##$in_fh->close();
    ##close(INFH);
    $fh->close();
    print $log "Counted: $info{reads_counted} and wrote: $info{reads_written} reads.\n";
    $log->close();
    $out_fh->close();
    return(\%info);
}

=head2 C<Racer>

 Use the RACER command from hitec to correct sequencer-based errors.
 10.1093/bioinformatics/btq653

=cut
sub Racer {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        compress => 1,
        length => 1000000,
        required => ['input', ],);
    my $input = $options->{input};
    my @input_list = split(/$options->{delimiter}/, $input);
    my @suffixes = ('.fastq', '.gz', '.xz');
    my $decompress_input = 0;
    my @base_list = ();
    for my $in (@input_list) {
        $decompress_input = 1 unless ($in =~ /\.fastq$/);
        my $shorted = basename($in, @suffixes);
        push(@base_list, basename($shorted, @suffixes));
    }
    my $jstring = qq'';
    my $output_dir = qq"outputs/$options->{jprefix}racer";
    my @created = make_path($output_dir);
    my @output_files;
    my $stdout = qq"${output_dir}/racer.stdout";
    my $stderr = qq"${output_dir}/racer.stderr";
    foreach my $c (0 .. $#input_list) {
        my $name = File::Temp::tempnam($output_dir, 'racer');
        my $output = qq"${output_dir}/$base_list[$c]-corrected.fastq";
        if ($decompress_input) {
            my $input_fc = Bio::Adventure::Get_FC(input => $input_list[$c]);
            $jstring .= qq"${input_fc} > ${name}.fastq\n";
        } else {
            $jstring .= qq!ln -sf "\$(pwd)"/$input_list[$c] ${name}.fastq\n!;
        }
        $jstring .= qq"RACER \\
  ${name}.fastq \\
  ${output} \\
  $options->{length} \\
  1>>${stdout} \\
  2>>${stderr}
rm ${name}.fastq
echo \"Finished correction of $input_list[$c].\" >> ${stdout}
";
        if ($options->{compress}) {
            $jstring .= qq"xz -9e -f ${output}
";
            push(@output_files, "${output}.xz");
        } else {
            push(@output_files, "${output}");
        }
    }

    my $output = '';
    for my $o (@output_files) {
        $output .= qq"${o}:";
    }
    $output =~ s/\:$//g;

    my $racer = $class->Submit(
        input => $input,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        outputs => $outputs,
        output => $output,
        stdout => $stdout,
        stderr => $stderr);
    return($racer);
}

=head2 C<Trimomatic>

 Invoke trimomatic to get sequencing data ready to play with.
 10.1093/bioinformatics/btu170

 Call trimomatic to remove adapters/low quality sequences.  If $args{input} has a
 ':' or ',' then this will assume the input is comprised of two pairwise files
 and will call 'Trimomatic_Pairwise()', otherwise 'Trimomatic_Single()'.

=cut
sub Trimomatic {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];
    my $extra_args = $class->Passthrough_Args(arbitrary => $options{arbitrary});
    my $num_inputs = scalar(@{$inputs});
    my $trim;
    if ($num_inputs == 1) {
        $options{paired} = 0;
        $trim = $class->Bio::Adventure::Trim::Trimomatic_Single(%options);
    } elsif ($num_inputs == 2) {
        $options{paired} = 1;
        $trim = $class->Bio::Adventure::Trim::Trimomatic_Pairwise(%options);
    } else {
        die("I do not understand this number of inputs: ${num_inputs}.");
    }
    return($trim);
}

=head2 C<Trimomatic_Pairwise>

 Invoke trimomatic with parameters suitable for pairwise sequence libraries.

=cut
sub Trimomatic_Pairwise {
    my ($class, %options) = @_;
    my $all_inputref = $class->inputs;
    my $inputs = $all_inputref->[$#$all_inputref];
    my $found_exe = 0;
    my $exe = q{};
    my %exe_list = (trimomatic => 'trimomatic PE',
                    TrimmomaticPE => 'TrimmomaticPE',
                    TrimomaticPE => 'TrimmomaticPE',
                    trimmomatic => 'trimmomatic PE');
    for my $test_exe (keys %exe_list) {
        if (which($test_exe)) {
            $exe = $exe_list{$test_exe};
        }
    }
    if (!defined($exe)) {
        die('Unable to find the trimomatic executable.');
    } else {
        $options{trimomatic_exe} = $exe;
    }
    unless ($options{adapter_file}) {
        $options{adapter_file} = dist_file('Bio-Adventure', 'genome/adapters.fa');
    }

    my @suff = ('.fastq', '.gz', '.xz', '.bz2');
    my $r1 = $inputs->[0]->{filename};
    my $r2 = $inputs->[1]->{filename};
    my $r1b = basename($r1, @suff);
    $r1b = basename($r1b, @suff);
    my $r2b = basename($r2, @suff);
    $r2b = basename($r2b, @suff);
    $options{input_string} = qq"${r1} ${r2}";
    if ($r1 =~ /\.fastq\.xz$/) {
        my $r1_fd = $class->Get_FD(input => $r1);
        my $r2_fd = $class->Get_FD(input => $r2);
        $options{input_string} = qq"${r1_fd} ${r2_fd}";
    }
    $options{thread_flags} = qq"-threads $options{jcpu}";

    my $r1o = qq"$options{output_dir}/${r1b}-trimmed.fastq";
    $options{r1o} = $r1o;
    my $r1op = qq"$options{output_dir}/${r1b}-trimmed_paired.fastq";
    $options{r1op} = $r1op;
    my $r1ou = qq"$options{output_dir}/${r1b}-trimmed_unpaired.fastq";
    $options{r1ou} = $r1ou;

    my $r2o = qq"$options{output_dir}/${r2b}-trimmed.fastq";
    $options{r2o} = $r2o;
    my $r2op = qq"$options{output_dir}/${r2b}-trimmed_paired.fastq";
    $options{r2op} = $r2op;
    my $r2ou = qq"$options{output_dir}/${r2b}-trimmed_unpaired.fastq";
    $options{r2ou} = $r2ou;

    $options{leader_flags} = '';
    if ($options{task} eq 'dnaseq') {
        $options{leader_flags} = 'HEADCROP:$options{min_quality} LEADING:3 TRAILING:3';
    }
    my $suffix_trim = '';
    if (defined($options{arbitrary}) &&
        $options{arbitrary} =~ /^CROP/) {
        $suffix_trim = $options{arbitrary};
    }

    my $output = qq"${r1op}:${r2op}";
    my $output_unpaired = qq"${r1ou}:${r2ou}";
    my $compress_string = '';
    if ($options{compress}) {
        $output = qq"${r1op}.xz:${r2op}.xz";
        $output_unpaired = qq"${r1ou}.xz:${r2ou}.xz";
        $options{compress_string} = qq"
## Recompress the unpaired reads, this should not take long.
xz -9e -f ${r1ou}
xz -9e -f ${r2ou}
## Recompress the paired reads.
xz -9e -f ${r1op}
xz -9e -f ${r2op}
";
    }

    if (-r $options{adapter_file}) {
        $options{adapter_flags} = qq"ILLUMINACLIP:$options{adapter_file}:2:$options{min_quality}:10:2";
    } else {
        carp("The adapter file: $options{adapter_file} is not readable.");
        $options{adapter_flags} = qq"ILLUMINACLIP:TruSeq3-PE.fa:2:$options{min_quality}:10:2";
    }
    $options{suffix_flags} = qq"${suffix_trim}";
    $options{window_flags} = qq"SLIDINGWINDOW:4:$options{min_quality} MINLEN:$options{min_length}";
    my $stdout = qq"$options{output_dir}/$options{basename}-trimomatic.stdout";
    my $stderr = qq"$options{output_dir}/$options{basename}-trimomatic.stderr";
    $options{template} = $options{template_pair};
    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    my $trim = $class->Submit();
    my $stats = $class->Bio::Adventure::Metadata::Trimomatic_Stats(
        input => $options{stderr},
        input_name => $options{output_base},
        jdepends => $trim->{job_id},
        paired => 1);
    $trim->{stats} = $stats;
    return($trim);
}

=head2 C<Trimomatic_Single>

 Invoke trimomatic with parameters suitable for single-read sequence libraries.

=cut
sub Trimomatic_Single {
    my ($class, %options) = @_;
    my $all_inputref = $class->inputs;
    my $inputs = $all_inputref->[$#$all_inputref];
    my $found_exe = 0;
    my $exe = q{};
    my %exe_list = (trimomatic => 'trimomatic SE',
                    TrimmomaticSE => 'TrimmomaticSE',
                    TrimomaticSE => 'TrimmomaticSE',
                    trimmomatic => 'trimmomatic SE');
    for my $test_exe (keys %exe_list) {
        if (which($test_exe)) {
            $exe = $exe_list{$test_exe};
        }
    }
    if (!defined($exe)) {
        die('Unable to find the trimomatic executable.');
    } else {
        $options{trimomatic_exe} = $exe;
    }
    unless ($options{adapter_file}) {
        $options{adapter_file} = dist_file('Bio-Adventure', 'genome/adapters.fa');
    }
    $options{leader_flags} = '';
    if ($options{task} eq 'dnaseq') {
        $options{leader_flags} = 'HEADCROP:20 LEADING:3 TRAILING:3';
    }
    my $suffix_trim = '';
    if (defined($options{arbitrary}) &&
        $options{arbitrary} =~ /^CROP/) {
        $suffix_trim = $options{arbitrary};
    }

    if (-r $options{adapter_file}) {
        $options{adapter_flags} = qq"ILLUMINACLIP:$options{adapter_file}:2:$options{min_quality}:10:2";
    } else {
        carp("The adapter file: $options{adapter_file} is not readable.");
        $options{adapter_flags} = qq"ILLUMINACLIP:TruSeq3-PE.fa:2:$options{min_quality}:10:2";
    }
    $options{suffix_flags} = qq"${suffix_trim}";
    $options{window_flags} = qq"SLIDINGWINDOW:4:$options{min_quality} MINLEN:$options{min_length}";
    my $input = $inputs->[0]->{filename};
    my $basename = $input;
    $basename = basename($basename, ('.gz', '.xz', '.bz2'));
    $basename = basename($basename, ('.fastq'));
    $options{input_string} = $input;
    if ($input =~ /\.fastq\.xz$/) {
        $options{input_string} = $class->Get_FD(input => $input);
    }
    $options{thread_flags} = qq"-threads $options{jcpu}";
    $options{compress_string} = '';
    if ($options{compress}) {
        $options{compress_string} = qq"
## Compress the trimmed reads.
xz -9e -f $options{output_r1}
";
    }
    $options{template} = $options{template_single};
    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    my $trim = $class->Submit();
    my $trim_stats = $class->Bio::Adventure::Metadata::Trimomatic_Stats(
        basename => $basename,
        input => $options{stderr},
        input_name => $options{output_base},
        jdepends => $trim->{job_id},
        paired => 0,
    );
    $trim->{stats} = $trim_stats;
    return($trim);
}

sub Umi_Tools_Dedup {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        count => 1,);
    my $jname = qq"umi_dedup";
    my $stdout = qq"$paths->{output_dir}/umi_dedup.stdout";
    my $stderr = qq"$paths->{output_dir}/umi_dedup.stderr";
    my $output = qq"$paths->{output_dir}/umi_tools_deduplicated.bam";
    my $jstring = qq!umi_tools dedup \\
  --output-stats=$paths->{output_dir}/umi_dedup_stats.txt \\
  --buffer-whole-contig \\
  -I $options->{input} \\
  -S ${output} \\
  2>${stdout} 1>${stderr}
samtools index ${output}
!;
    my $umi_job = $class->Submit(
        jdepends => $options->{jdepends},
        input => $options->{input},
        jstring => $jstring,
        outputs => $outputs,
        output => $output,
        stderr => $stderr,
        stdout => $stdout,);

  my $htmulti;
  my $new_jprefix = qq"$options->{jprefix}_1";
  if ($options->{count}) {
    $htmulti = $class->Bio::Adventure::Count::HT_Multi(
      input => $umi_job->{output},
      jdepends => $umi_job->{job_id},
      jname => qq"htmulti_umi_dedup",
      jprefix => $new_jprefix,);
   }
   $umi_job->{recount} = $htmulti;
   return($umi_job);
}

sub Umi_Tools_Extract {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        extract_method => 'regex',
        extract_string => '^(?P<umi_1>.{8})(?P<discard_1>.{6}).*',
        required => ['input'],);
    my $jname = qq"umi_tools";
    my $test_file;
    my $stranded = 'yes';
    my $output = qq"$paths->{output_dir}/r1_extracted.fastq.gz";
    my $umi_input = qq"-I $options->{input}";
    my $umi_output = qq"-S ${output}";
    my $paired = 0;
    my @pair_listing = ();
    if ($options->{input} =~ /$options->{delimiter}/) {
        @pair_listing = split(/$options->{delimiter}/, $options->{input});
        $paired = 1;
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        $umi_input = qq"-I $pair_listing[0] -2 $pair_listing[1] ";
        ## I expect umi_tools has similar problems reading from an anonymous bash
        ## handle as per other python programs.
        if ($pair_listing[0] =~ /\.[x|b]z$/) {
            my $r1_fd = $class->Get_FD(input => $pair_listing[0]);
            my $r2_fd = $class->Get_FD(input => $pair_listing[1]);
            $umi_input = qq"-I ${r1_fd} --read2-in ${r2_fd}";
            $umi_output = qq"-S $paths->{output_dir}/r1_extracted.fastq.gz --read2-out $paths->{output_dir}/r2_extracted.fastq.gz";
        } else {
            $umi_input = qq"-I $pair_listing[0] --read2-in $pair_listing[1]";
            $umi_output = qq"-S $paths->{output_dir}/r1_extracted.fastq.gz --read2-out $paths->{output_dir}/r2_extracted.fastq.gz";
        }
        $output .= qq":$paths->{output_dir}/r2_extracted.fastq.gz";
    } else {
        $stranded = 'no';
        $test_file = File::Spec->rel2abs($options->{input});
        if ($test_file =~ /\.[x|b]z$/) {
            ## It is noteworthy that I modified hisat2 on my computer so this is never necessary.
            my $r1_fd = $class->Get_FD(input => $test_file);
            $umi_input = qq" -I ${r1_fd}";
        }
    }
    my $stdout = qq"$paths->{output_dir}/umi_tools.stdout";
    my $stderr = qq"$paths->{output_dir}/umi_tools.stderr";
    my $jstring = qq!mkdir -p $paths->{output_dir}
umi_tools extract \\
 ${umi_input} \\
 ${umi_output} \\
 --extract-method $options->{extract_method} \\
 --bc-pattern2 '$options->{extract_string}' \\
 --log $paths->{output_dir}/extract.log \\
 2>${stderr} \\
 1>${stdout}
!;
    my $umi_job = $class->Submit(
        jdepends => $options->{jdepends},
        input => $umi_input,
        jstring => $jstring,
        outputs => $outputs,
        output => $output,
        stderr => $stderr,
        stdout => $stdout,);
    return($umi_job);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<trimomatic> L<cutadapt>

=cut

1;
