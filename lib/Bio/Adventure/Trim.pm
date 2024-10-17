package Bio::Adventure::Trim;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::Adventure::Config;
use Cwd qw"abs_path getcwd cwd";
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
    my $options = $class->Get_Vars(
        args => \%args,
        modules => ['cogent'],
        required => ['input', 'species', 'input_umi'],
        jcpu => 4,
        jprefix => '01',
        jmem => 24,
        jwalltime => 36,
        type => 'Stranded_UMI',);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $cwd_name = basename(cwd());
    my @input_filenames = split(/$options->{delimiter}/, $options->{input});
    my $output_dir = qq"outputs/$options->{jprefix}cogent";
    my $comment = qq!## This is a cogent demultiplexing/trimming job.
!;
    my $stdout = qq"${output_dir}/cogent.stdout";
    my $stderr = qq"${output_dir}/cogent.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
cogent demux \\
  -i $input_filenames[0] \\
  -p $input_filenames[1] \\
  -b $options->{input_umi} \\
  -t $options->{type} \\
  -o ${output_dir}/demux \\
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
  -i ${output_dir}/demux \\
  -g $options->{species} \\
  -o ${output_dir}/analyze \\
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
        comment => $comment,
        jcpu => $options->{jcpu},
        jdepends => $options->{jdepends},
        jname => "cogent_${job_name}",
        jqueue => 'large',
        jstring => $jstring,
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
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => undef,
        maxlength => 42,
        maxerr => 0.1,
        maxremoved => 3,
        minlength => 8,
        left => undef,
        right => undef,
        either => undef,
        type => 'rnaseq',
        jmem => 12,
        jwalltime => '48:00:00',
        jprefix => '12',);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $minlength = $options->{minlength};
    my $maxlength = $options->{maxlength};
    my $arbitrary = '';
    if (defined($options->{arbitrary})) {
        $arbitrary = $options->{arbitrary};
    }

    my $input = $options->{input};
    my @suffixes = split(/\,/, $options->{suffixes});
    my $basename = basename($input, @suffixes);
    $basename = basename($basename, @suffixes);
    my $adapter_flags = "";

    if ($options->{type} eq 'old_tnseq') {
        $adapter_flags = ' -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATC -a ACAGTCCCCGGTCTGACACATCTCCCTAT -a ACAGTCCNCGGTCTGACACATCTCCCTAT ';
        $maxlength = 20;
    }
    elsif ($options->{type} eq 'riboseq') {
        $adapter_flags = ' -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC ';
        $minlength = 16;
        $maxlength = 42;
    }
    else {
        $minlength = 30;
        $maxlength = undef;
        $adapter_flags = ' -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC ';
    }

    my $comment = qq!## This script makes use of biopieces and cutadapt to trim away adapters
## and separate the sequence file into a few pieces depending on size
## and adapter status.  It also performs some simple graphs of the data.!;
    my $out_dir = qq"outputs/$options->{jprefix}cutadapt";
    my $stderr = qq"${out_dir}/cutadapt.stderr";
    my $stdout = qq"${out_dir}/cutadapt.stdout";
    my $type_flag = '';
    my $out_suffix = 'fastq';
    my $input_fc = Bio::Adventure::Get_FC(input => $input);
    my $input_flags = qq"${input_fc} | cutadapt - ";
    if ($input =~ /\.csfasta/) {
        $type_flag = '-c -t --strip-f3';
        $options = $class->Get_Vars(
            args => \%args,
            required => ['qual'],
        );
        $input_flags = qq"${input_fc} | cutadapt - $options->{qual} "; ## If we are keeping quality files
        $out_suffix = 'fastq';
    }
    my $minarg = '';
    my $too_short = undef;
    if ($minlength) {
        $too_short = qq"${out_dir}/${basename}_too_short.${out_suffix}";
        $minarg = qq" -m ${minlength} --too-short-output=${too_short}";
    }
    my $maxarg = '';
    my $too_long = undef;
    if ($maxlength) {
        $too_long = qq"${out_dir}/${basename}_too_long.${out_suffix}";
        $maxarg = qq" -M ${maxlength} --too-long-output=${too_long}";
    }
    my $maxerr = '';
    if ($options->{maxerr}) {
        $maxerr = qq" -e $options->{maxerr}";
    }
    my $output = qq"${out_dir}/${basename}-trimmed_ca.${out_suffix}";
    my $compressed_out = qq"${output}.xz";
    my $jstring = qq!## I am putting the args in an odd order in case some are not defined.
mkdir -p ${out_dir}
${input_flags} \\
  ${type_flag} ${adapter_flags} ${arbitrary} \\
  -o ${output} ${minarg} \\
  -n $options->{maxremoved} ${maxerr} ${maxarg} \\
  --untrimmed-output=${out_dir}/${basename}_untrimmed.${out_suffix} \\
  2>${stderr} \\
  1>${stdout}
xz -9e -f ${output}
xz -9e -f ${out_dir}/${basename}_untrimmed.${out_suffix}
!;
    if (defined($too_short)) {
        $jstring .= qq!
xz -9e -f ${too_short}
!;
    }
    if (defined($too_long)) {
        $jstring .= qq!
xz -9e -f ${too_long}
!;
    }

    my $cutadapt = $class->Submit(
        comment => $comment,
        input => $input,
        jmem => $options->{jmem},
        jname => qq"cutadapt_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'workstation',
        jstring => $jstring,
        jwalltime => '8:00:00',
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        output => $compressed_out,);
    if ($options->{type} eq 'old_tnseq') {
        my $ta_check = $class->Bio::Adventure::TNSeq::TA_Check(
            comment => '## Check that TAs exist.',
            input => $compressed_out,
            jdepends => $cutadapt->{job_id},
            jname => qq"tacheck_${job_name}",
            jprefix => $options->{jprefix} + 1,);
        $cutadapt->{tacheck} = $ta_check;
    }
    my $stats = $class->Bio::Adventure::Metadata::Cutadapt_Stats(
        input => $cutadapt->{stdout},
        jcpu => 1,
        jmem => 1,
        jwalltime => '00:03:00',
        jdepends => $cutadapt->{job_id},);
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
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => undef,
        jmem => 12,
        jwalltime => '24:00:00',
        jprefix => '12',);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});

    my $extra_args = $class->Passthrough_Args(arbitrary => $options->{arbitrary});
    $extra_args .= ' -D ' if ($options->{deduplication});
    $extra_args .= ' -c ' if ($options->{correction});
    $extra_args .= ' -y ' if ($options->{complexity});

    my $comment = qq!## Run fastp on raw data
!;
    my $out_dir = qq"outputs/$options->{jprefix}fastp";
    my $stderr = qq"${out_dir}/fastp.stderr";
    my $stdout = qq"${out_dir}/fastp.stdout";
    my $input_flags = '';
    my $num_inputs = scalar(@{$inputs});
    if ($num_inputs == 2) {
        my $r1 = $inputs->[0]->{filename};
        my $r2 = $inputs->[1]->{filename};
        my $r1_base = $inputs->[0]->{filebase_extension};
        my $r2_base = $inputs->[1]->{filebase_extension};
        my $output_dir = $inputs->[0]->{directory};
        my $output_r1 = qq"${output_dir}/${r1_base}-fastp.fastq";
        my $output_r2 = qq"${output_dir}/${r2_base}-fastp.fastq";
        $input_flags = qq" -i ${r1} -o ${output_r1} \\
  -I ${r2} -O ${output_r2} ";
    }
    elsif ($num_inputs == 1) {
        my $r1 = $inputs->[0]->{filename};
        my $r1_base = $inputs->[0]->{filebase_extension};
        my $output_dir = $inputs->[0]->{directory};
        my $output_r1 = qq"${output_dir}/${r1_base}-fastp.fastq";
        $input_flags = qq" -i ${r1} -o ${output_r1} ";
    }
    else {
        die("An unusual number of inputs was provided: ${num_inputs}.");
    }
    my $umi_flags = '';
    if ($options->{do_umi}) {
        $umi_flags = ' -U ';
    }
    my $report_flags = qq"-h ${out_dir}/fastp_report.html -j ${out_dir}/fastp_report.json";

    my $jstring = qq!
mkdir -p ${out_dir}
fastp ${umi_flags} ${input_flags} \\
  ${report_flags} ${extra_args} \\
  2>${stderr} \\
  1>${stdout}
!;

    my $fastp = $class->Submit(
        comment => $comment,
        input => $options->{input},
        jmem => $options->{jmem},
        jname => qq"fastp_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => '8:00:00',
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($fastp);
}

sub PolyA_Extractor {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        min_A => 10,
        min_seq => 10,
        jmem => 24,
        jprefix => '10',
        jwalltime => '40:00:00',
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
    my $paths = $class->Get_Paths($final_output);
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
        my $job_name = $class->Get_Job_Name(input => $input_file);
        print "TESTME: $job_name and $final_output\n";
        my $jname = qq"polyAextract_${job_name}";
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
            comment => $comment,
            input => $input_file,
            jdepends => $options->{jdepends},
            jmem => $options->{jmem},
            jname => $jname,
            jprefix => $options->{jprefix},
            jstring => $jstring,
            language => 'perl',
            min_A => $options->{min_A},
            min_seq => $options->{min_seq},
            output => $final_output,
            output_dir => $output_dir,
            stdout => $stdout,
            stderr => $stderr,);
        $options->{jdepends} = $polyA->{job_id};
    }
    return($polyA);
}

sub PolyA_Extractor_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        length => 1000000,
        jmem => 24,
        jprefix => '10',
        jwalltime => '40:00:00',
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
        }
        elsif ($sequence =~ /^T{$options->{min_A},}/) {
            $seq = $seq->revcom;
            $sequence = $seq->seq;
            $match_begin = $len - $+[0];
        }
        else {
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
    }                           ## End iterating over every sequence.
    ##$in_fh->close();
    ##close(INFH);
    $fh->close();
    print $log "Counted: ${reads_counted} and wrote: ${reads_written} reads.\n";
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
    my $options = $class->Get_Vars(
        args => \%args,
        compress => 1,
        length => 1000000,
        jmem => 24,
        jprefix => '10',
        jwalltime => '40:00:00',
        required => ['input', ],);
    my $job_name = $class->Get_Job_Name();
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
    my $comment = qq!## This calls RACER to try to remove
## arbitrary errors in sequencing data.
!;
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
        }
        else {
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
        }
        else {
            push(@output_files, "${output}");
        }
    }

    my $output = '';
    for my $o (@output_files) {
        $output .= qq"${o}:";
    }
    $output =~ s/\:$//g;

    my $racer = $class->Submit(
        comment => $comment,
        input => $input,
        jcpu => 4,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => "racer_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'workstation',
        jstring => $jstring,
        jwalltime => '12:00:00',
        prescript => $options->{prescript},
        postscript => $options->{postscript},
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
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        compress => 1,
        input_paired => undef,
        length => 50,
        quality => '20',
        jcpu => 4,
        jmem => 24,
        jprefix => '01',
        jwalltime => '36:00:00',);
    my $trim;
    if ($options->{input} =~ /$options->{delimiter}/) {
        $trim = $class->Bio::Adventure::Trim::Trimomatic_Pairwise(%{$options});
    }
    else {
        $trim = $class->Bio::Adventure::Trim::Trimomatic_Single(%{$options});
    }
    return($trim);
}

=head2 C<Trimomatic_Pairwise>

 Invoke trimomatic with parameters suitable for pairwise sequence libraries.

=cut
sub Trimomatic_Pairwise {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        compress => 1,
        jcpu => 4,
        jmem => 24,
        jprefix => '01',
        jwalltime => '48:00:00',
        length => 50,
        quality => '20',
        required => ['input',],);
    my $output_dir = qq"outputs/$options->{jprefix}trimomatic";
    my $job_name = $class->Get_Job_Name();
    my $exe = undef;
    my $found_exe = 0;
    my %modules = Bio::Adventure::Get_Modules(caller => 1);
    my $loaded = $class->Module_Loader(%modules);
    my %exe_list = (trimomatic => 'trimomatic PE',
                    TrimmomaticSE => 'TrimmomaticPE',
                    TrimomaticSE => 'TrimmomaticPE',
                    trimmomatic => 'trimmomatic PE');
    for my $test_exe (keys %exe_list) {
        if (which($test_exe)) {
            $exe = $exe_list{$test_exe};
        }
    }
    if (!defined($exe)) {
        die('Unable to find the trimomatic executable.');
    }
    my $adapter_file = dist_file('Bio-Adventure', 'genome/adapters.fa');
    my $input = $options->{input};
    my @input_list = split(/$options->{delimiter}/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = $class->Bio::Adventure::Trim::Trimomatic_Single(input => $input,);
        return($ret);
    }
    my $r1 = $input_list[0];
    my $r2 = $input_list[1];
    my @suff = ('.fastq', '.gz', '.xz');

    my $basename = basename($r1, @suff);
    $basename = basename($basename, @suff);
    $basename =~ s/_R1$//g;
    $basename =~ s/_forward$//g;

    my $r1b = basename($r1, @suff);
    $r1b = basename($r1b, @suff);
    my $r2b = basename($r2, @suff);
    $r2b = basename($r2b, @suff);
    my $reader = qq"${r1} ${r2}";
    if ($r1 =~ /\.fastq\.xz$/) {
        my $r1_fd = $class->Get_FD(input => $r1);
        my $r2_fd = $class->Get_FD(input => $r2);
        $reader = qq"${r1_fd} ${r2_fd}";
    }
    my $r1o = qq"${r1b}-trimmed.fastq";
    my $r1op = qq"${r1b}-trimmed_paired.fastq";
    my $r1ou = qq"${r1b}-trimmed_unpaired.fastq";

    my $r2o = qq"${r2b}-trimmed.fastq";
    my $r2op = qq"${r2b}-trimmed_paired.fastq";
    my $r2ou = qq"${r2b}-trimmed_unpaired.fastq";

    my $leader_trim = '';
    #if ($options->{task} eq 'dnaseq') {
    #    $leader_trim = 'HEADCROP:20 LEADING:3 TRAILING:3';
    #}
    my $suffix_trim = '';
    if (defined($options->{arbitrary}) &&
        $options->{arbitrary} =~ /^CROP/) {
        $suffix_trim = $options->{arbitrary};
    }

    my $output = qq"${r1o}:${r2o}";
    my $output_unpaired = qq"${r1ou}:${r2ou}";
    my $compress_string = '';
    if ($options->{compress}) {
        $output = qq"${r1o}.xz:${r2o}.xz";
        $output_unpaired = qq"${r1ou}.xz:${r2ou}.xz";
        $compress_string = qq"
## Recompress the unpaired reads, this should not take long.
xz -9e -f ${r1ou}
xz -9e -f ${r2ou}
## Recompress the paired reads.
xz -9e -f ${r1o}
xz -9e -f ${r2o}
ln -sf ${r1o}.xz r1_trimmed.fastq.xz
ln -sf ${r2o}.xz r2_trimmed.fastq.xz
";
    }

    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $stdout = qq"${output_dir}/${basename}-trimomatic.stdout";
    my $stderr = qq"${output_dir}/${basename}-trimomatic.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
## Note that trimomatic prints all output and errors to STDERR, so send both to output
${exe} \\
  -threads 1 \\
  -phred33 \\
  ${reader} \\
  ${r1op} ${r1ou} \\
  ${r2op} ${r2ou} \\
  ${leader_trim} ILLUMINACLIP:${adapter_file}:2:$options->{quality}:10:2:keepBothReads ${suffix_trim} \\
  SLIDINGWINDOW:4:$options->{quality} MINLEN:$options->{length} \\
  1>${stdout} \\
  2>${stderr}
excepted=\$( { grep "Exception" "${output_dir}/${basename}-trimomatic.stdout" || test \$? = 1; } )
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "\${excepted}" \!= "" ]]; then
  ${exe} \\
    -threads 1 \\
    -phred33 \\
    ${reader} \\
    ${r1op} ${r1ou} \\
    ${r2op} ${r2ou} \\
    ${leader_trim} SLIDINGWINDOW:4:25 MINLEN:$options->{length} \\
    1>${stdout} \\
    2>${stderr}
fi
sleep 10
mv ${r1op} ${r1o}
mv ${r2op} ${r2o}
!;

    if ($options->{compress}) {
        $jstring .= qq"${compress_string}\n";
    }

    ## Example output from trimomatic:
    ## Input Read Pairs: 10000 Both Surviving: 9061 (90.61%) Forward Only Surviving: 457 (4.57%) Reverse Only Surviving: 194 (1.94%) Dropped: 288 (2.88%)
    ## Perhaps I can pass this along to Get_Stats()
    my $trim = $class->Submit(
        args => \%args,
        comment => $comment,
        input => $input,
        jcpu => $options->{jcpu},
        jmem => $options->{jmem},
        jname => qq"trim_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'workstation',
        jstring => $jstring,
        jwalltime => $options->{jwalltime},
        length => $options->{length},
        output => $output,
        output_unpaired => $output_unpaired,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stdout => $stdout,
        stderr => $stderr);
    my $new_prefix = qq"$options->{jprefix}_1";
    my $trim_stats = $class->Bio::Adventure::Metadata::Trimomatic_Stats(
        basename => $basename,
        jdepends => $trim->{job_id},
        jcpu => 1,
        jprefix => $new_prefix,
        jname => "trst_${job_name}",
        jwalltime => '00:03:00',
        pairwise => 1,
        input => $stderr,
        output_dir => $output_dir,);
    $trim->{stats} = $trim_stats;
    return($trim);
}

=head2 C<Trimomatic_Single>

 Invoke trimomatic with parameters suitable for single-read sequence libraries.

=cut
sub Trimomatic_Single {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        compress => 1,
        jmem => 24,
        jcpu => 4,
        jprefix => '01',
        jwalltime => '48:00:00',
        length => 50,
        quality => '20',
        required => ['input',],);
    my $exe = undef;
    my $found_exe = 0;
    my %modules = Bio::Adventure::Get_Modules(caller => 1);
    my $loaded = $class->Module_Loader(%modules);
    my %exe_list = (trimomatic => 'trimomatic SE',
                    TrimmomaticSE => 'TrimmomaticSE',
                    TrimomaticSE => 'TrimmomaticSE',
                    trimmomatic => 'trimmomatic SE');
    for my $test_exe (keys %exe_list) {
        if (which($test_exe)) {
            $exe = $exe_list{$test_exe};
        }
    }
    my $output_dir = qq"outputs/$options->{jprefix}trimomatic";
    if (!defined($exe)) {
        die('Unable to find the trimomatic executable.');
    }
    my $adapter_file = dist_file('Bio-Adventure', 'genome/adapters.fa');
    my $leader_trim = "";
    #if (defined($options->{task}) && $options->{task} eq 'dnaseq') {
    #    $leader_trim = 'HEADCROP:20 LEADING:3 TRAILING:3';
    #}

    my $input = $options->{input};
    my $basename = $input;
    $basename = basename($basename, ('.gz', '.xz', '.bz2'));
    $basename = basename($basename, ('.fastq'));
    my $job_name = $class->Get_Job_Name();
    my $output = qq"${basename}-trimmed.fastq";
    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $stdout = qq"${output_dir}/${basename}-trimomatic.stdout";
    my $stderr = qq"${output_dir}/${basename}-trimomatic.stderr";
    my $r1_fd = $class->Get_FD(input => $input);
    my $jstring = qq!mkdir -p ${output_dir}
## Note that trimomatic prints all output and errors to STDERR, so send both to output
${exe} \\
  -phred33 \\
  ${r1_fd} \\
  ${output} \\
  ${leader_trim} ILLUMINACLIP:${adapter_file}:2:30:10 \\
  SLIDINGWINDOW:4:25 MINLEN:$options->{length} \\
  1>${stdout} 2>${stderr}
!;
    my $compress_string = '';
    if ($options->{compress}) {
        $compress_string = qq"
## Compress the trimmed reads.
xz -9e -f ${output}
ln -sf ${output}.xz r1_trimmed.fastq.xz
";
    }
    $jstring .= $compress_string;
    my $trim = $class->Submit(
        comment => $comment,
        input => $input,
        jcpu => $options->{jcpu},
        jmem => $options->{jmem},
        jname => qq"trim_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => $options->{jwalltime},
        length => $options->{length},
        output => $output,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    my $trim_stats = $class->Bio::Adventure::Metadata::Trimomatic_Stats(
        basename => $basename,
        input => $stderr,
        jcpu => 1,
        jdepends => $trim->{job_id},
        jname => qq"trst_${job_name}",
        jprefix => '06',
        jwalltime => '00:03:00',
        output_dir => $output_dir,
    );
    $trim->{stats} = $trim_stats;
    return($trim);
}

sub Umi_Tools_Dedup {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 48,
        jwalltime => '36:00:00',
        jprefix => '04',);
    my $jname = qq"umi_dedup";
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $stdout = qq"$paths->{output_dir}/umi_dedup.stdout";
    my $stderr = qq"$paths->{output_dir}/umi_dedup.stderr";
    my $output = qq"$paths->{output_dir}/umi_tools_deduplicated.bam";
    my $jstring = qq!mkdir -p $paths->{output_dir}
umi_tools dedup \\
  --output-stats=$paths->{output_dir}/umi_dedup_stats.txt \\
  --buffer-whole-contig \\
  -I $options->{input} \\
  -S ${output} \\
  2>${stdout} 1>${stderr}
samtools index ${output}
!;
    my $comment = qq!## This is a umi_tools deduplication script
!;
    my $umi_job = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        input => $options->{input},
        jname => $jname,
        jstring => $jstring,
        output => $output,
        stderr => $stderr,
        stdout => $stdout,);
    return($umi_job);
}

sub Umi_Tools_Extract {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        extract_method => 'regex',
        extract_string => '^(?P<umi_1>.{8})(?P<discard_1>.{6}).*',
        required => ['input'],
        jmem => 20,
        jprefix => '01',
        jwalltime => '8:00:00');
    my $jname = qq"umi_tools";
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
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
        }
        else {
            $umi_input = qq"-I $pair_listing[0] --read2-in $pair_listing[1]";
            $umi_output = qq"-S $paths->{output_dir}/r1_extracted.fastq.gz --read2-out $paths->{output_dir}/r2_extracted.fastq.gz";
        }
        $output .= qq":$paths->{output_dir}/r2_extracted.fastq.gz";
    }
    else {
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
    my $comment = qq!## This is a umi_tools extraction (soon to be deduplication) script
## currently it supports only the Takara UMI scheme.
!;
    my $umi_job = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        input => $umi_input,
        jname => $jname,
        jstring => $jstring,
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
