package Bio::Adventure::QA;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use feature 'try';
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use File::Basename;
use File::Which qw"which";

=head1 NAME

 Bio::Adventure::QA - Use trimomatic/cutadapt/etc to trim libraries

=head1 SYNOPSIS

 The functions here invoke various quality assurance programs.

=head1 METHODS

=head2 C<Biopieces_Graph>

 Biopieces is an older but fun toolkit.
 https://github.com/maasha/biopieces

 Reads in a fastq file and uses biopieces to make some graphs describing the
 sequences therein.

=cut
sub Biopieces_Graph {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 8,
        jprefix => '02',);
    my $check = which('read_fastq');
    die("Could not find biopieces in your PATH.") unless($check);
    my $input = $options->{input};
    my $jname = qq"biop_$options->{job_basename}";
    my @inputs = split(/$options->{delimiter}/, $input);
    my $comment = '## This script uses biopieces to draw some simple graphs of the sequence.';
    my $output_dir = 'output/biopieces';
    my $stdout = qq"${output_dir}/biopieces.stdout";
    my $stderr = qq"${output_dir}/biopieces.stderr";
    my $bp;
    if (scalar(@inputs) > 1) { ## multiple comma/colon/semicolon inputs were provided.
        foreach my $in (@inputs) {
            my $short_in = basename($in, ('.gz','.fastq'));
            $short_in = basename($short_in, ('.gz','.fastq'));
            my $input_fc = Bio::Adventure::Get_FC(input => $in);
            my $jstring = qq!
## Do not forget that _only_ the last command in a biopieces string is allowed to have the -x.
mkdir -p outputs/biopieces
${input_fc} | read_fastq -i - -e base_$options->{phred} |\\
 plot_scores -T 'Quality Scores' -t svg -o outputs/biopieces/${short_in}_quality_scores.svg |\\
 plot_nucleotide_distribution -T 'NT. Distribution' -t svg -o outputs/biopieces/${short_in}_ntdist.svg |\\
 plot_lendist -T 'Length Distribution' -k SEQ_LEN -t svg -o outputs/biopieces/${short_in}_lendist.svg |\\
 analyze_gc |\\
     bin_vals -b 20 -k 'GC%' |\\
     plot_distribution -k 'GC%_BIN' -t svg -o outputs/biopieces/${short_in}_gcdist.svg |\\
 analyze_gc |\\
     mean_vals -k 'GC%' -o outputs/biopieces/${short_in}_gc.txt |\\
 count_records -o outputs/biopieces/${short_in}_count.txt -x
!;
            $bp = $class->Submit(
                input => $in,
                comment => $comment,
                jdepends => $options->{jdepends},
                jmem => $options->{jmem},
                jname => qq"${jname}_${in}",
                jprefix => $options->{jprefix},
                jstring => $jstring,
                stderr => $stderr,
                stdout => $stdout,
                prescript => $args{prescript},
                postscript => $args{postscript},);
        }
    } else { ## A single input was provided
        my $input_fc = Bio::Adventure::Get_FC(input => $input);
        my $jstring = qq!
## Do not forget that _only_ the last command in a biopieces string is allowed to have the -x.
mkdir -p outputs/biopieces
${input_fc} | read_fastq -i - -e base_33 |\\
 plot_scores -T 'Quality Scores' -t svg -o outputs/biopieces/${jname}_quality_scores.svg |\\
 plot_nucleotide_distribution -T 'NT. Distribution' -t svg -o outputs/biopieces/${jname}_ntdist.svg |\\
 plot_lendist -T 'Length Distribution' -k SEQ_LEN -t svg -o outputs/biopieces/${jname}_lendist.svg |\\
 analyze_gc |\\
     bin_vals -b 20 -k 'GC%' |\\
     plot_distribution -k 'GC%_BIN' -t svg -o outputs/biopieces/${jname}_gcdist.svg |\\
 analyze_gc |\\
     mean_vals -k 'GC%' -o outputs/biopieces/${jname}_gc.txt |\\
 count_records -o outputs/biopieces/${jname}_count.txt -x
!;
        $bp = $class->Submit(
            input => $input,
            comment => $comment,
            jdepends => $options->{jdepends},
            jmem => $options->{jmem},
            jname => 'biop',
            jprefix => $options->{jprefix},
            jstring => $jstring,
            stderr => $stderr,
            stdout => $stdout,
            prescript => $options->{prescript},
            postscript => $options->{postscript},);
    }
    return($bp);
}

=head2 C<Fastqc>

 Invoke fastqc on some sequence files.
 https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

=cut
sub Fastqc {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        filtered => 'unfiltered',
        required => ['input',],);
    my $outdir = $paths->{output_dir};
    my $dirname = $inputs->[0]->{dirname};

    my $fastqc_job;
    my $input_file_string = '';
    my $subshell = 0;
    my $modified_input = undef;
    if (scalar(@{$inputs}) > 1) {
        my @input_files;
        for my $element (@{$inputs}) {
            push(@input_files, $element->{fullpath});
        }
        for my $in (@input_files) {
            $modified_input = basename($in, ('.gz', '.bz2', '.xz')) unless(defined($modified_input));
            $input_file_string = qq"$input_file_string ${in} ";
            if ($in =~ /\.xz$|\.bz2$/) {
                $subshell = 1;
                my $r1_fd = $class->Get_FD(input => $in);
                $input_file_string = qq"$input_file_string ${r1_fd} ";
            }
        }
    } elsif ($options->{input} =~ /\.xz$|\.bz2$/) {
        $modified_input = basename($options->{input}, ('.gz', '.bz2', '.xz')) unless ($modified_input);
        $subshell = 1;
        my $r1_fd = $class->Get_FD(input => $options->{input});
        $input_file_string = qq"${r1_fd} ";
    } else {
        $modified_input = basename($options->{input}, ('.gz')) unless ($modified_input);
        $input_file_string = qq"$options->{input} ";
    }

    ## This is where fastqc should put its various output files
    ## with one important exception: It tries to be smart and take its own basename of the input
    ## but since I am using a bash subshell <(), it will get /dev/fd/xxx and so put the outputs into
    ## outputs/${jprefix}fastqc/xxx_fastqc...
    my $jstring = qq!which perl 2>$outputs->{stderr} 1>&2
which fastqc 2>>$outputs->{stderr} 1>&2
## Even if fastqc finishes happily using a subshell, it might exit with SIGERR
trap - ERR
finished=\$(fastqc --extract \\
  -o $paths->{output_dir} \\
  ${input_file_string} \\
  2>>$outputs->{stderr} \\
  1>>$outputs->{stdout} || /bin/true)
echo "fastqc finished with \$? andor \${finished}."
!;
    if ($subshell) {
        $jstring .= qq!
## Note that if this is using a subshell, fastqc will assume that the inputs
## are /dev/fd/xx (usually 60 something).
## We can likely cheat and get the subshell fd with this:
read_num=1
for badname in  \$(cd $paths->{output_dir} && /bin/ls -d *_fastqc); do
  echo \${badname}
  mv $paths->{output_dir}/\${badname}.html $paths->{output_dir}/r\${read_num}_trimmed_fastqc.html 2>/dev/null
  echo "html move finished with: \$?"
  mv $paths->{output_dir}/\${badname}.zip $paths->{output_dir}/r\${read_num}_trimmed_fastqc.zip 2>/dev/null
  echo "zip move finished with: $?"
  mv $paths->{output_dir}/\${badname} outputs/01fastqc/r\${read_num}_trimmed_fastqc 2>/dev/null
  echo "directory move finished with: \$?"
  let "read_num++"
done
!;
    }
    my $fqc = $class->Submit(
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        outputs => $outputs,);
    my $stats = $class->Bio::Adventure::Metadata::Fastqc_Stats(
        input => $outputs->{main},
        jdepends => $fqc->{job_id},);
    $fqc->{stats} = $stats;
    return($fqc);
}

sub MultiQC {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        input => 'preprocessing',
        jprefix => '12',);
    my $comment = '## A Multiqc run!';
    my $stderr = 'multiqc.stderr';
    my $stdout = 'multiqc.stdout';
    my $jstring = qq!
cd $options->{input}
multiqc --no-ansi . 2>${stderr} 1>${stdout}
!;
    my $multiqc = $class->Submit(
        comment => $comment,
        input => $options->{input},
        jname => qq"multiqc_$options->{job_basename}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($multiqc);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<fastqc> L<biopieces>

=cut

1;
