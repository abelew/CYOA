package Bio::Adventure::Splicing;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use feature 'try';
no warnings 'experimental::try';

use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Capture::Tiny qw":all";
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Spec;
use File::Path qw"make_path rmtree";
use File::Which qw"which";
use File::ShareDir qw":ALL";
use IPC::Open2;
use List::MoreUtils qw"uniq";
use Scalar::Util;
use Spreadsheet::Read;
use String::Approx qw"amatch";
use Symbol qw"gensym";


=head2 C<Suppa>

  Set up and invoke Suppa, a differential transcript and splicing
  event calculator.  10.1186/s13059-018-1417-1

  I am putting this function in this file because it is relevant to
  the gene structure, which is admittedly a bit of a stretch.

  This function should be able to invoke the tool which creates the
  catalog of splicing events (suppa.py generateEvents), create the
  table of TPM values by transcript.  I am thinking to use the sample
  sheet as the input variable with a parameter containing the
  salmon/etc filename and have this generate the associated table.
  Given that information, this should be able to run the psiPerIsoform
  function and/or psiPerEvent.  Note, it turns out that suppa provides
  a helper function for creating the expression table 'joinFiles'; so
  that will simplify things quite a bit.

=cut
sub Suppa {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        type => 'SE',
        condition_column => 'drug',
        file_column => 'hg38100salmonfile',
        mapper => 'salmon',
        diff_method => 'classical',
        paired => 0,
        jprefix => '90',);

    my $jname = qq"suppa_$options->{species}_$options->{condition_column}";
    my $suppa_dir = qq"outputs/$options->{jprefix}${jname}";
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $gtf_file = $paths->{gtf};
    die("Suppa requires a gtf file: ${gtf_file}, which is missing.") unless (-r $gtf_file);
    if (!-r $options->{input}) {
        die("Cannot open input spreadsheet: $options->{input}");
    }

    ## The following basenames will define the locations of the output files and logs
    ## for each step performed.
    my $event_base = qq"${suppa_dir}/events";
    my $tpm_base = qq"${suppa_dir}/tpm";
    my $psi_base = qq"${suppa_dir}/psi";
    my $diff_base = qq"${suppa_dir}/diff";
    my $cluster_base = qq"${suppa_dir}/cluster";

    ## suppa seems to assume htseq-like inputs: e.g. a file with 2
    ## columns: transcript IDs (with no version ID suffix), and TPM
    ## values.  Given that I am likely to use salmon at least
    ## sometimes, I will need some logic to reformat the salmon output
    ## into an input suitable for suppa.  I chose to use awk.
    my $awk_string = '';
    my $directory_string = qq!
mkdir -p ${event_base}
mkdir -p ${tpm_base}
mkdir -p ${psi_base}
mkdir -p ${diff_base}
mkdir -p ${cluster_base}
!;
    my $event_string = '';
    ## I think we can use the same invocation string for every event
    ## type?  I will at least endeavor to do so.
    my @event_types = ('SE', 'SS', 'MX', 'RI', 'FL');
    my %mapped_types = (
        'SE' => 'SE',
        'SS' => 'A3',
        'MX' => 'MX',
        'RI' => 'RI',
        'FL' => 'AF');
    my @alt_types = ('SE', 'A3', 'A5', 'AF', 'AL', 'MX', 'RI');
    ## Read in the gtf (NOT GTF) annotation file and produce a series
    ## of event-type specific annotations.
    for my $ty (@event_types) {
        my $check_type = $mapped_types{$ty};
        $event_string .= qq!
if [[ -r "${event_base}/local_as_events_${check_type}_strict.ioe" ]]; then
  echo "The ${ty} events already exist." >> ${event_base}.stdout
else
  echo "Generating ${ty} events." >> ${event_base}.stdout
  suppa.py generateEvents -i ${gtf_file} -f ioe -e ${ty} -o ${event_base}/local_as_events \\
    2>>${event_base}.stderr 1>>${event_base}.stdout
fi
!;
    }
    ## In addition, create one file for all local events.
    $event_string .= qq!
if [[ -r "${event_base}/transcript_events.ioi" ]]; then
  echo "The ioi events already exist." >> ${event_base}.stdout
else
  echo "Generating the ioi events." >> ${event_base}.stdout
  suppa.py generateEvents -i ${gtf_file} -f ioi -o ${event_base}/transcript_events \\
    2>>${event_base}.stderr 1>>${event_base}.stdout
fi
!;

    ## The following few lines read in the sample sheet and extract
    ## the samples associated with each experimental condition.
    ## Reminder to self: if the various spreadsheet reader modules are
    ## not installed, this just returns a somewhat unhelpful undef.
    my $reader = Spreadsheet::Read->new($options->{input});
    my $sheet = $reader->sheet(1);
    ## Strangely, Spreadsheet::Read only uses numeric values for columns, so grab the first row
    ## and get the number of my column from it...
    my $file_number = undef;
    my $condition_number = undef;
    my $numeric_condition = Scalar::Util::looks_like_number($options->{condition_column});
    my $numeric_file = Scalar::Util::looks_like_number($options->{file_column});
    if ($numeric_condition) {
        $condition_number = $options->{condition_column};
    }
    if ($numeric_file) {
        $file_number = $options->{file_column};
    }
    if (!defined($condition_number) && !defined($file_number)) {
        my @column_names = $sheet->cellrow(1);
        my $count = 0;
      COLUMNS: for my $name (@column_names) {
            $count++;
            if ($name eq $options->{file_column}) {
                print "Found $options->{file_column} as number: $count\n";
                $file_number = $count;
            }
            if ($name eq $options->{condition_column}) {
                print "Found $options->{condition_column} as number: $count\n";
                $condition_number = $count;
            }
        }
    }

    ## Keep in mind that cellcolumn is 1-indexed and it provides an
    ## array which includes the header as the first element.  Thus if
    ## I want the sample IDs, I pull cellcolumn(1) and if I do not
    ## want to include the header (which I don't), I will need to
    ## shift (or just skip it) the first element off the result.
    ## Either way, I want arrays containing the sample IDs, filenames,
    ## and conditions.
    my @samplenames = $sheet->cellcolumn(1);
    shift @samplenames;
    my @filenames = $sheet->cellcolumn($file_number);
    my $filename_heading = shift @filenames;
    my @conditions = $sheet->cellcolumn($condition_number);
    my $condition_heading = shift @conditions;
    ## Map samples to files and conditions
    my %samples_to_files = ();
    my %samples_to_conditions = ();
    ## Mape conditions to lists of samples/files.
    my %files_by_condition = ();
    my %samples_by_condition = ();
    my $sample_count = -1;
  SAMPLENAMES: for my $s (@samplenames) {
        $sample_count++;
        my $f = $filenames[$sample_count];
        my $cond = $conditions[$sample_count];
        if (!defined($s)) {
            print "The samplename at row ${sample_count} is not defined, skipping it.\n";
            next SAMPLENAMES;
        }
        if (!defined($f)) {
            print "This sample: ${s} does not have a file, skipping it.\n";
            next SAMPLENAMES;
        }
        if (!-r $f) {
            print "The file for ${s} is not readable, skipping it.\n";
            next SAMPLENAMES;
        }
        if (defined($options->{numerator}) && defined($options->{denominator})) {
            if ($cond eq $options->{numerator}) {
                print "Adding ${f} sample ${s} as numerator: $options->{numerator}.\n";
            }
            elsif ($cond eq $options->{denominator}) {
                print "Adding ${f} sample ${s} as denominator: $options->{denominator}.\n";
            }
            else {
                print "Skipping sample ${s}, it is neither numerator nor denominator.\n";
                next SAMPLENAMES;
            }
        }
        $samples_to_files{$s} = $f;
        $samples_to_conditions{$s} = $cond;
        print "Working on $cond\n";
        my @cond_samples = ();
        my @cond_filenames = ();
        if (defined($files_by_condition{$cond})) {
            @cond_filenames = @{$files_by_condition{$cond}};
            @cond_samples = @{$samples_by_condition{$cond}};
        }
        push(@cond_filenames, $f);
        push(@cond_samples, $s);
        $files_by_condition{$cond} = \@cond_filenames;
        $samples_by_condition{$cond} = \@cond_samples;
    }                           ## End iterating over the samples

    if ($filename_heading ne $options->{file_column}) {
        print "Something seems wrong with the name of the file entry.\n";
    }
    if ($condition_heading ne $options->{condition_column}) {
        print "Something seems wrong with the name of the condition entry.\n";
    }

    ## Now start creating the invocation strings for the various steps
    ## suppa will perform.  The first step is already complete: the
    ## creation of the event inputs.  Next we need to join the tpm's
    ## into one file / condition; then we create the values/event (psi).
    my $join_string = '';
    my $psi_string = '';
    my $psi_local_string = '';
    ## Thus we will iterate over every condition _and_ file/condition
    ## in order to: 1. copy/reformat the TPM data 2. join them into one file/condition
    my @condition_names = sort keys %files_by_condition;
    for my $cond (@condition_names) {
        my $file_string = '';
        my @sample_names = @{$samples_by_condition{$cond}};
        my $count = -1;
      FILES: for my $f (@{$files_by_condition{$cond}}) {
            next FILES unless(defined($f));
            $count++;
            my $sample = $sample_names[$count];
            ## The following if handles salmon/kallisto inputs; htseq
            ## will require an explicit tpm conversion.
            if ($options->{mapper} eq 'salmon') {
                $awk_string .= qq!
if [[ -r "${tpm_base}/${sample}.tpm" ]]; then
  echo "The tpm file for ${sample} already exists."
else
  awk '{printf("\%s\\t\%s\\n", \$1, \$4)}' ${f} > ${tpm_base}/${sample}.tpm
  perl -p -i.unmodified -e 's/^(\\w+)\\.\\d+(\\s+.*\$)/\$1\$2/g' ${tpm_base}/${sample}.tpm
fi
!;
            }
            elsif ($options->{mapper} eq 'kallisto') {
                ## I think I can just copy outputs from rsem too? (check that)
                $awk_string = qq!cp ${f} ${tpm_base}/${sample}.tpm
!;
            }
            else {
                die("I need a tpm conversion here and have not added it yet.");
            }
            $file_string .= qq!${tpm_base}/${sample}.tpm !;
        }
        ## Now use suppa's joinFiles function to reformat the inputs
        ## into its (bizarre to me) expected format.
        $join_string .= qq!
if [[ -r "${tpm_base}/{cond}.tpm" ]]; then
  echo "The ${cond}.tpm already exists." >>${tpm_base}.stdout
else
  echo "Joining files for condition: ${cond}." >> ${tpm_base}.stdout
  suppa.py joinFiles -f tpm -i ${file_string} -o ${tpm_base}/${cond} \\
    2>>${tpm_base}.stderr 1>>${tpm_base}.stdout
fi
!;

        ## Given the full set of annotations and tpm inputs, create a
        ## PSI for all events.
        $psi_string .= qq!
if [[ -r "${psi_base}/${cond}_isoform.psi" ]]; then
  echo "The ${cond} psi file already exists." >>${psi_base}.stdout
else
  echo "Performing psiPerIsoform for ${cond}." >>${psi_base}.stdout
  suppa.py psiPerIsoform -g ${gtf_file} -e ${tpm_base}/${cond}.tpm \\
   -o ${psi_base}/${cond} \\
    2>>${psi_base}.stderr 1>>${psi_base}.stdout
fi
!;
        ## But also create per-type specific PSI files.  The output
        ## filenames get a prefix of 'local', which may be a poor choice?
        for my $ty (@alt_types) {
            $psi_local_string .= qq!
if [[ -r "${psi_base}/${cond}_local_${ty}.psi" ]]; then
  echo "The ${cond} ${ty} psi file already exists." >>${psi_base}.stdout
else
  echo "Performing psiPerEvent: ${ty} for ${cond}." >>${psi_base}.stdout
  suppa.py psiPerEvent --ioe-file ${event_base}/local_as_events_${ty}_strict.ioe \\
    --expression-file ${tpm_base}/${cond}.tpm \\
    -o ${psi_base}/${cond}_local_${ty} \\
    2>>${psi_base}.stderr 1>>${psi_base}.stdout
fi
!;
        }
    }

    ## The differential splicing invocations require the we iterate
    ## over the 2chooseN pairs of conditions.
    my $diffsplice_ioi_string = '';
    my $diffsplice_ioe_string = '';
    my $cluster_string = '';
    ## Note, suppa has a syntax for performing all pairwise
    ## comparisons, but it doesn't make much sense to me; as a result
    ## I think I will need to move the files with the suffix dpsi.temp.0
    ## to .dpsi because I am explicitly performing each pairwise
    ## comparison myself and therefore the little code block in
    ## diff_tools.py is not getting called.
    my $end = $#condition_names;
    print "Examining @condition_names\n";
    for my $c (0 .. ($end - 1)) {
        my $next = $c + 1;
        my $numerator_sample_number = 0;
        my $denominator_sample_number = 0;
        for my $d ($next .. $end) {
            ## I do not actually yet know for 100% certain that suppa
            ## puts the denominator second.  Check this.
            my $denominator = $condition_names[$c];
            my $numerator = $condition_names[$d];
            my @denominator_samples = @{$samples_by_condition{$denominator}};
            my @numerator_samples = @{$samples_by_condition{$numerator}};
            my $first_start = 1;
            my $first_end = scalar(@numerator_samples);
            my $second_start = $first_end + 1;
            my $second_end = $first_end + scalar(@denominator_samples);
            my $group_string = qq"${first_start}-${first_end},${second_start}-${second_end}";

            ## Note, these extra flags should be handled via
            ## arbitrary_args; but for now just fill them in until I
            ## understand them.
            my $extra_flags = '--area 1000 --lower-bound 0.05 -gc ';
            if ($options->{paired}) {
                $extra_flags .= '-pa ';
            }
            my $extra_cluster_flags = '-c DBSCAN -m euclidean --sig-threshold 0.05 --eps 0.05 ';
            $diffsplice_ioi_string .= qq!
if [[ -r "${diff_base}/${numerator}_${denominator}_ioi.dpsi" ]]; then
  echo "The dpsi file exists for ${numerator} vs ${denominator}." >>${diff_base}.stdout
else
  echo "Starting ioi diffSplice between ${numerator} and ${denominator}." >>${diff_base}.stdout
  suppa.py diffSplice -s --method empirical --input ${event_base}/transcript_events.ioi \\
    --psi ${psi_base}/${numerator}_isoform.psi ${psi_base}/${denominator}_isoform.psi \\
    --tpm ${tpm_base}/${numerator}.tpm ${tpm_base}/${denominator}.tpm \\
    ${extra_flags} -o ${diff_base}/${numerator}_${denominator}_ioi_empirical \\
    2>>${diff_base}.stderr 1>>${diff_base}.stdout
  if [[ -r "${diff_base}/${numerator}_${denominator}_ioi_empirical.dpsi.temp.0" ]]; then
    mv ${diff_base}/${numerator}_${denominator}_ioi_empirical.dpsi.temp.0 \\
      ${diff_base}/${numerator}_${denominator}_ioi_empirical.dpsi
  fi
  suppa.py diffSplice --method classical --input ${event_base}/transcript_events.ioi \\
    --psi ${psi_base}/${numerator}_isoform.psi ${psi_base}/${denominator}_isoform.psi \\
    --tpm ${tpm_base}/${numerator}.tpm ${tpm_base}/${denominator}.tpm \\
    ${extra_flags} -o ${diff_base}/${numerator}_${denominator}_ioi_classical \\
    2>>${diff_base}.stderr 1>>${diff_base}.stdout
  if [[ -r "${diff_base}/${numerator}_${denominator}_ioi_classical.dpsi.temp.0" ]]; then
    mv ${diff_base}/${numerator}_${denominator}_ioi_classical.dpsi.temp.0 \\
      ${diff_base}/${numerator}_${denominator}_ioi_classical.dpsi
  fi
fi
!;
            $cluster_string .= qq!
suppa.py clusterEvents --dpsi ${diff_base}/${numerator}_${denominator}_ioi.dpsi \\
  --psivec ${diff_base}/${numerator}_${denominator}_ioi.psivec \\
  ${extra_cluster_flags} --groups ${group_string} \\
  -o ${numerator}_${denominator}_ioi.clusters \\
  2>>${cluster_base}.stderr 1>>${cluster_base}.stdout
!;
            for my $ty (@alt_types) {
                $diffsplice_ioe_string .= qq!
if [[ -f "${diff_base}/${numerator}_${denominator}_${ty}.dpsi" ]]; then
  echo "The ${numerator} vs ${denominator} ${ty} dpsi file exists." >>${diff_base}.stdout
else
  echo "Starting ioe diffSplice between ${numerator} and ${denominator}, type: ${ty}." >>${diff_base}.stdout
  suppa.py diffSplice -s --method $options->{diff_method} --input ${event_base}/local_as_events_${ty}_strict.ioe \\
    --psi ${psi_base}/${numerator}_local_${ty}.psi ${psi_base}/${denominator}_local_${ty}.psi \\
    --tpm ${tpm_base}/${numerator}.tpm ${tpm_base}/${denominator}.tpm \\
    ${extra_flags} -o ${diff_base}/${numerator}_${denominator}_${ty} \\
    2>>${diff_base}.stderr 1>>${diff_base}.stdout
  if [[ -r "${diff_base}/${numerator}_${denominator}_${ty}.dpsi.temp.0" ]]; then
    mv ${diff_base}/${numerator}_${denominator}_${ty}.dpsi.temp.0 ${diff_base}/${numerator}_${denominator}_${ty}.dpsi
  fi
fi
!;
                $cluster_string .= qq!
suppa.py clusterEvents --dpsi ${diff_base}/${numerator}_${denominator}_${ty}.dpsi \\
  --psivec ${diff_base}/${numerator}_${denominator}_${ty}.psivec \\
  ${extra_cluster_flags} --groups ${group_string} \\
  -o ${numerator}_${denominator}_${ty}.clusters \\
  2>>${cluster_base}.stderr 1>>${cluster_base}.stdout
!;
            }
        }
    }

    my $comment = qq"# Generate the event types for suppa.";
    my $jstring = qq!${directory_string}

${event_string}

${awk_string}

${join_string}

${psi_string}

${psi_local_string}

${diffsplice_ioi_string}

${diffsplice_ioe_string}

${cluster_string}
!;
    my $suppa_job = $class->Submit(
        comment => $comment,
        cpus => $options->{jcpu},
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => qq"${suppa_dir}/suppa_generate.stdout",
        stdout => qq"${suppa_dir}/suppa_generate.stdout",
        stderr => qq"${suppa_dir}/suppa_generate.stderr",);
    return($suppa_job);
}

=head2 C<SLSearch>

 Search a pile of reads for the trypanosome spliced leader.

 Use some simple pattern matching on a pile of reads to look for
 sequences of interest.  By default this looks for a portion of the
 spliced leader sequence from Leishmania major.  Having written and
 used this, I realized that it is ... dumb.  I should have used
 jellyfish and simply counted up the hits across a range of the SL.
 Doing this with jellyfish would also let me count up each window of
 the SL and plot a histogram of the occurrences, thus showing the
 optimal starting point to discriminate the SL as opposed to my ad hoc
 choice of position generated via repeated invocations of 'grep | wc'.

 With that in mind, this counts up the number of times a string of
 interest appears in the input files in the forward and RC directions,
 and prints that in an easy-to-read format.

=over

=item C<Arguments>

 input(required): Set of input files to read.
 search('AGTTTCTGTACTTTATTGG'): SL substring to search.
 jmem(24): Memory to allocate for this task.
 jprefix(50): Default jobname/output prefix.

=back

=cut
sub SLSearch {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => 'lmajor',
        jmem => 24,
        jprefix => '50',
        seed => 15,);

    my $paths = $class->Bio::Adventure::Config::Get_Paths();

    ## Ideally, I would like to run this function on the mapping
    ## directories and search the aligned/unaligned sequences
    ## separately.
    ## Unfortunately, in our cleaning, we deleted those files.
    ## Therefore it (for the moment) will instead just read over the
    ## trimmed files.
    ## I think I will write the following: the read IDs for every
    ## read which matches the forward or RC versions of the SL and
    ## a summary txt file of the result.

    ## I am going to collect the full SL sequences from various
    ## species here.  These are going to exclude the first few bases
    ## with the assumption that the various methylations observed on
    ## them will leave them more susceptible to sequencing
    ## shenanigans. (AACTA in Keith's document, I assume a G is on
    ## there?)
    my $paired = 0;
    my $observed_sl = {
        tcruzi => ['ACGCTATTATTGATACAGTTTCTGTACTATATTG'],
        lmajor => ['ACGCTATATAAGTATCAGTTTCTGTACTTTATTG'],
    };
    my $input_string = qq"";
    my @input_lst = ();
    if ($options->{input} =~ /$options->{delimiter}/) {
        $paired = 1;
        @input_lst = split(/$options->{delimiter}/, $options->{input});
        for my $i (@input_lst) {
            $input_string .= "${i} ";
        }
    } else {
        if (-r $options->{input}) {
            push(@input_lst, $options->{input});
            $input_string = $options->{input};
        }
    }
    if (scalar(@input_lst) == 0) {
        return(undef);
    }

    my @fwd_sl = ();
    my @rev_sl = ();
    my @full_observed = @{$observed_sl->{$options->{species}}};
    for my $sl (@full_observed) {
        my $fwd = substr($sl, ($options->{seed} * -1));
        push(@fwd_sl, $fwd);
        my $rev = reverse($fwd);
        $rev =~ tr/AGTCUagtcu/tcagatcaga/;
        push(@rev_sl, $rev);
    }

    my $log_file = qq"$paths->{output_dir}/slsearch_log.txt";
    print "TESTME: $log_file and $options->{input} and paired: $paired\n";
    my $log_fh = FileHandle->new(">${log_file}");
    print $log_fh qq"Searching for putative SL/polyA containing reads
in the file(s): $options->{input}.\n";
    ## Now start the main loop, open file handles for a main log
    ## and for the per-input outputs.  Create a couple of global counters.

    my $fwd_sl_string = qq"";
    for my $fw (uniq @fwd_sl) {
        $fwd_sl_string .= "${fw} ";
    }
    my $rev_sl_string = qq"";
    for my $rv (uniq @rev_sl) {
        $rev_sl_string .= "${rv} ";
    }

    my $fwd_polya_string = qq"AAAAAAAAAAAA";
    my $rev_polya_string = qq"TTTTTTTTTTTT";
    my $fwd_cmd = qq"
cutadapt -a ${fwd_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_fwdorient_threepr1_r1.fastq.gz -p $paths->{output_dir}/fwdsl_fwdorient_threepr1_r2.fastq.gz \\
   --interleaved  ${input_string}
cutadapt -g ${fwd_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_revorient_fivepr1_r1.fastq.gz -p $paths->{output_dir}/fwdsl_revorient_fivepr1_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -A ${fwd_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_fwdorient_threepr2_r1.fastq.gz -p $paths->{output_dir}/fwdsl_fwdorient_threepr2_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -G ${fwd_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_revorient_fivepr2_r1.fastq.gz -p $paths->{output_dir}/fwdsl_revorient_fivepr2_r2.fastq.gz \\
   --interleaved ${input_string}

cutadapt -a ${fwd_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_fwdorient_threepr1_r1.fastq.gz -p $paths->{output_dir}/fwdpolya_fwdorient_threepr1_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -g ${fwd_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_revorient_fivepr1_r1.fastq.gz -p $paths->{output_dir}/fwdpolya_revorient_fivepr1_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -A ${fwd_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_fwdorient_threepr2_r1.fastq.gz -p $paths->{output_dir}/fwdpolya_fwdorient_threepr2_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -G ${fwd_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_revorient_fivepr2_r1.fastq.gz -p $paths->{output_dir}/fwdpolya_revorient_fivepr2_r2.fastq.gz \\
   --interleaved ${input_string}

";
    my $rev_cmd = qq"
cutadapt -a ${rev_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revsl_fwdorient_threepr1_r1.fastq.gz -p $paths->{output_dir}/revsl_fwdorient_threepr1_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -g ${rev_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revsl_revorient_fivepr1_r1.fastq.gz -p $paths->{output_dir}/revsl_revorient_fivepr1_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -A ${rev_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revsl_fwdorient_threepr2_r1.fastq.gz -p $paths->{output_dir}/revsl_fwdorient_threepr2_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -G ${rev_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revsl_revorient_fivepr2_r1.fastq.gz -p $paths->{output_dir}/revsl_revorient_fivepr2_r2.fastq.gz \\
   --interleaved ${input_string}

cutadapt -a ${rev_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revpolya_fwdorient_threepr1_r1.fastq.gz -p $paths->{output_dir}/revpolya_fwdorient_threepr1_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -g ${rev_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revpolya_revorient_fivepr1_r1.fastq.gz -p $paths->{output_dir}/revpolya_revorient_fivepr1_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -A ${rev_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revpolya_fwdorient_threepr2_r1.fastq.gz -p $paths->{output_dir}/revpolya_fwdorient_threepr2_r2.fastq.gz \\
   --interleaved ${input_string}
cutadapt -G ${rev_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revpolya_revorient_fivepr2_r1.fastq.gz -p $paths->{output_dir}/revpolya_revorient_fivepr2_r2.fastq.gz \\
   --interleaved ${input_string}
";
    unless ($paired) {
        $fwd_cmd = qq"
cutadapt -b ${fwd_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_fwdorient_r1.fastq.gz \\
   ${input_string}
cutadapt -B ${fwd_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_revorient_r1.fastq.gz \\
   ${input_string}

cutadapt -b ${fwd_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_fwdorient_r1.fastq.gz \\
   ${input_string}
cutadapt -B ${fwd_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_revorient_r1.fastq.gz \\
   ${input_string}
";
        $rev_cmd = qq"
cutadapt -b ${rev_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revsl_fwdorient_r1.fastq.gz \\
   ${input_string}
cutadapt -B ${rev_sl_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revsl_revorient_r1.fastq.gz \\
   ${input_string}

cutadapt -b ${rev_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revpolya_fwdorient_r1.fastq.gz \\
   ${input_string}
cutadapt -B ${rev_polya_string} -b ${rev_polya_string} -e 1 --trimmed-only \\
   -o $paths->{output_dir}/revpolya_revorient_r1.fastq.gz \\
   ${input_string}
";
    }

    my $stderr = qq"$paths->{output_dir}/slsearch_cutadapt.stderr";
    my $stdout = qq"$paths->{output_dir}/slsearch_cutadapt.stdout";
    my $comment = qq!## This cutadapt invocation seeks to extract SL/polyA reads.!;

    my $jstring = qq!mkdir -p $paths->{output_dir}
${fwd_cmd} 2>>${stderr} 1>>${stdout}
${rev_cmd} 2>>${stderr} 1>>${stdout}
!;

    my $jname = qq"slsearch_$options->{species}";
    my $sl_extract = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        input => $options->{input},
        jname => $jname,
        jmem => $options->{jmem},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        stderr => $stderr,
        stdout => $stdout,
    );

}

=head2 C<SL_UTR>

 Search a pile of reads for UTR boundaries with the SL/polyA.

 I have been rereading Keith's code for gathering UTRs from SL sites.
 I like it quite a lot and would like to use it to improve my python
 style, but it will require quite a bit of retooling to work with the
 newer generation of tools/methods.  As a result I am thinking I will
 first port the logic here, then once I am reasonably confident I can
 get a 'correct' answer I will turn back to that code.

=over

=item C<Arguments>

 input(required): Set of input files to read.
 search('AGTTTCTGTACTTTATTGG'): SL substring to search.
 jmem(24): Memory to allocate for this task.
 jprefix(50): Default jobname/output prefix.

=back

=cut
sub SL_UTR {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 24,
        jprefix => '50',
        search => 'AGTTTCTGTACTTTATTGG',);
    my $output_dir =  qq"outputs/$options->{jprefix}SL_UTR";
    my $output_made = make_path($output_dir);
    my $comment = '## Search for SL sub-sequences.';
    my $stdout = qq"${output_dir}/sl_utr.stdout";
    my $stderr = qq"${output_dir}/sl_utr.stderr";
    my $jstring = qq?
my \$result = \$h->Bio::Adventure::Splicing::SL_UTR_Worker(
  input => '$options->{input}',
  output => '${output_dir}',
  jprefix => '$options->{jprefix}',
  jname => 'slsearch',);
?;
    my $sl_utr = $class->Submit(
        comment => $comment,
        input => $options->{input},
        jcpu => 1,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'sl_utr',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        output => $output_dir,
        stderr => $stderr,
        stdout => $stdout,);
    return($sl_utr);
}

=head2 C<SL_UTR_Worker>

 This function does the actual work for SL_UTR().

 A good sample for testing: HPGL0453: This is an amazonensis procyclic sample.

=cut
sub SL_UTR_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'output'],
        jmem => 24,
        jprefix => '50',
        polya_search => 'AAAAAAAAA',
        sl_search => 'AGTTTCTGTACTTTAT',);
    ## Ideally, I would like to run this function on the mapping
    ## directories and search the aligned/unaligned sequences
    ## separately.
    ## Unfortunately, in our cleaning, we deleted those files.
    ## Therefore it (for the moment) will instead just read over the
    ## trimmed files.
    ## I think I will write the following: the read IDs for every
    ## read which matches the forward or RC versions of the SL and
    ## a summary txt file of the result.
    my @input_lst = ();
    if ($options->{input} =~ /$options->{delimiter}/) {
        my @tmp_lst = split(/$options->{delimiter}/, $options->{input});
        for my $i (@tmp_lst) {
            push (@input_lst, $i) if (-r $i);
        }
    }
    else {
        if (-r $options->{input}) {
            push(@input_lst, $options->{input});
        }
    }
    if (scalar(@input_lst) == 0) {
        return(undef);
    }

    my $sl_fwd_search = $options->{sl_search};
    my $sl_rc_search = reverse($sl_fwd_search);
    $sl_rc_search =~ tr/ATGCUatgcu/TACGAtacga/;
    my $search_length = length($sl_fwd_search);
    my $polya_fwd_search = $options->{polya_search};
    ## Doing the explicit reverse complement is silly, yes, but I might extend this to
    ## stuff like polyY
    my $polya_rc_search = reverse($polya_fwd_search);
    $polya_rc_search =~ tr/ATGCUatgcu/TACGAtacga/;

    my $log_fh = FileHandle->new(">$options->{output}/slsearch_log.txt");
    print $log_fh qq"Starting search for a portion of the SL sequence: $options->{search}
in the file(s): $options->{input}.\n";
    ## Now start the main loop, open file handles for a main log
    ## and for the per-input outputs.  Create a couple of global counters.

    my %global_search_result = (
        found => 0,     ## The total number of observed SL.
        fwd_found => 0, ## The number found in the forward orientation.
        rc_found => 0,  ## The number of revcomp found.
        searched => 0); ## The total number of sequences searched.
    ## One may reasonably ask why all_found is not just the sum of the next two.
    ## I am thinking to catch the pathological case where a single input sequence
    ## has both the forward and RC sequences.

    for my $i (@input_lst) {
        my $ind_name = basename($i, ('.gz', '.bz2', '.xz'));
        my $format = 'Fastq';
        $format = 'Fasta' if ($i =~ /\.fasta/);

        $ind_name = basename($ind_name, ('.fastq', '.fasta'));
        my $ind_result = FileHandle->new(">$options->{output}/${ind_name}.tsv");
        my $output_sl_reads = FileHandle->new(">$options->{output}/${ind_name}_sl_reads.fasta");
        my $output_polya_reads = FileHandle->new(">$options->{output}/${ind_name}_polya_reads.fasta");
        ## I want to save the position information too so that I can see if I my search string
        ## is too close to the beginning of the read.
        print $ind_result "ReadID\tposition\torientation\n";
        my %ind_search_result = (
            found => 0,
            sl_fwd_found => 0,
            sl_rc_found => 0,
            polya_fwd_found => 0,
            polya_rc_found => 0,
            searched => 0);
        my $reader = FileHandle->new("less ${i} |");
        my $seqio = Bio::SeqIO->new(-format => $format, -fh => $reader);
      FSA: while (my $seq = $seqio->next_seq) {
            $ind_search_result{searched}++;
            $global_search_result{searched}++;
            my $seq_id = $seq->id;
            my $read_seq = $seq->seq;
            my $fwd_end = undef;
            if ($read_seq =~ m/$sl_fwd_search/g) {
                print "TESTME: Found forward SL: $read_seq\n";
                $ind_search_result{sl_fwd_found}++;
                $fwd_end = pos($read_seq);
                my $new_read = $read_seq;
                $new_read =~ s/^.*$sl_fwd_search//g;
                print $output_sl_reads ">${seq_id} fwd
$new_read\n";
            } elsif ($read_seq =~ m/$polya_fwd_search/g) {
                print "TESTME: Found forward polyA: $read_seq\n";
                $ind_search_result{polya_fwd_found}++;
                $fwd_end = pos($read_seq);
                my $new_read = $read_seq;
                $new_read =~ s/^.*$polya_fwd_search//g;
                print $output_polya_reads ">${seq_id} fwd
$new_read\n";
            }
            my $rc_end = undef;
            if ($read_seq =~ m/^.*$sl_rc_search/g) {
                print "TESTME: Found reverse SL: $read_seq\n";
                $ind_search_result{sl_rev_found}++;
                $rc_end = pos($read_seq);
                my $new_read = reverse($read_seq);
                $new_read =~ tr/AGCTU/tcgaa/;
                $new_read =~ s/^.*$sl_fwd_search//g;
                print $output_sl_reads ">${seq_id} rev
${new_read}\n";
            } elsif ($read_seq =~ m/^.*$polya_rc_search/g) {
                print "TESTME: Found reverse polya: $read_seq\n";
                $ind_search_result{polya_rev_found}++;
                $rc_end = pos($read_seq);
                my $new_read = reverse($read_seq);
                $new_read =~ tr/AGCTU/tcgaa/;
                $new_read =~ s/^.*$polya_fwd_search//g;
                print $output_polya_reads ">${seq_id} rev
${new_read}\n";
            }
            ## Get out if we do not find the SL portion.
            if (!defined($fwd_end) && !defined($rc_end)) {
                next FSA;
            }

            if ($fwd_end) {
                my $fwd_start = $fwd_end - ($search_length - 1);
                print $ind_result "${seq_id}\t${fwd_start}\tFWD\n";
                $ind_search_result{found}++;
                $ind_search_result{fwd_found}++;
                $global_search_result{found}++;
                $global_search_result{fwd_found}++;
            }

            if ($rc_end) {
                my $rc_start = $rc_end - ($search_length - 1);
                print $ind_result "${seq_id}\t${rc_start}\tRC\n";
                $ind_search_result{found}++;
                $ind_search_result{rc_found}++;
                $global_search_result{found}++;
                $global_search_result{rc_found}++;
            }
        }  ## End reading the input fastq/fasta file.

        print $log_fh qq"
${ind_name} results:
  sequences searched: $ind_search_result{searched}
  subsequences observed: $ind_search_result{found}
  forward observed: $ind_search_result{fwd_found}
  reverse-complement observed: $ind_search_result{rc_found}\n";
        $ind_result->close();
        $reader->close();
    }                           ## End of the input loop
    print $log_fh qq"
Total results:
  sequences searched: $global_search_result{searched}
  subsequences observed: $global_search_result{found}
  forward observed: $global_search_result{fwd_found}
  reverse-complement observed: $global_search_result{rc_found}\n";
    $log_fh->close();
}



1;
