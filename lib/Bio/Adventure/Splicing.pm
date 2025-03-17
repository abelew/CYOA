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


sub RMats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        type => 'SE',
        numerator => 't4h',
        denominator => 'no',
        condition_column => 'drug',
        file_column => 'bamfile',
        mapper => 'hisat',
        paired => 0,
        jcpu => 16,
        jprefix => '90',
        read_length => 100,);
    my $jname = qq"rmats_$options->{species}_$options->{condition_column}";
    my $rmats_dir = qq"outputs/$options->{jprefix}${jname}";
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $gtf_file = $paths->{gtf};
    die("Rmats requires a gtf file: ${gtf_file}, which is missing.") unless (-r $gtf_file);
    if (!-r $options->{input}) {
        die("Cannot open input spreadsheet: $options->{input}");
    }

    ## The following basenames will define the locations of the output files and logs
    ## for each step performed.
    my $rmats_tmp = qq"${rmats_dir}/tmp";
    my $bam_group1 = qq"${rmats_dir}/group1.txt";
    my $bam_group2 = qq"${rmats_dir}/group2.txt";
    my $made = make_path($rmats_tmp);

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
    my $group1_string = '';
    my $group2_string = '';
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
                $group1_string .= qq"${f},";
            } elsif ($cond eq $options->{denominator}) {
                print "Adding ${f} sample ${s} as denominator: $options->{denominator}.\n";
                $group2_string .= qq"${f},";
            } else {
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
    } ## End iterating over the samples
    $group1_string =~ s/\,$//g;
    $group2_string =~ s/\,$//g;
    my $g1_fh = FileHandle->new(qq">${bam_group1}");
    print $g1_fh qq"${group1_string}\n";
    $g1_fh->close();
    my $g2_fh = FileHandle->new(qq">${bam_group2}");
    print $g2_fh qq"${group2_string}\n";
    $g2_fh->close();
    my $stdout = qq"${rmats_dir}/rmats.stdout";
    my $stderr = qq"${rmats_dir}/rmats.stderr";
    my $comment = qq"## Invoke rMats.";
    my $jstring = qq!mkdir -p ${rmats_tmp}
rmats.py \\
  --b1 ${bam_group1} \\
  --b2 ${bam_group2} \\
  --gtf ${gtf_file} \\
  -t paired --nthread $options->{jcpu} --od ${rmats_dir} --tmp ${rmats_tmp} \\
  --readLength $options->{read_length} --variable-read-length \\
  2>${stderr} 1>${stdout}
!;
    my $rmats_job = $class->Submit(
        comment => $comment,
        cpus => $options->{jcpu},
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => qq"${rmats_dir}/rmats.stdout",
        stdout => qq"${rmats_dir}/rmats.stdout",
        stderr => qq"${rmats_dir}/rmats.stderr",);
    return($rmats_job);
}

sub Spladder {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        type => 'SE',
        numerator => 't4h',
        denominator => 'no',
        condition_column => 'drug',
        file_column => 'bamfile',
        mapper => 'hisat',
        paired => 0,
        jcpu => 16,
        jprefix => '90',
        read_length => 100,);
    my $jname = qq"spladder_$options->{species}_$options->{condition_column}";
    my $spladder_dir = qq"outputs/$options->{jprefix}${jname}";
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $gtf_file = $paths->{gtf};
    die("Spladder requires a gtf file: ${gtf_file}, which is missing.") unless (-r $gtf_file);
    if (!-r $options->{input}) {
        die("Cannot open input spreadsheet: $options->{input}");
    }

    my $made = make_path($spladder_dir);
    my $stdout = qq"${spladder_dir}/spladder.stdout";
    my $stderr = qq"${spladder_dir}/spladder.stderr";
    my $info = $class->Bio::Adventure::Metadata::Get_Metadata_Column(
        input => $options->{input},
        condition_column => $options->{condition_column},
        wanted_column => $options->{file_column},
        type => 'file');
    my @groups = keys %{$info->{wanted_by_cond}};
    my %group_strings = ();
    my $all_bam = '';
    for my $g (@groups) {
        my $group_string = join(',', @{$info->{wanted_by_cond}{$g}});
        $group_strings{$g} = $group_string;
        $all_bam .= qq"${group_string},";
    }
    $all_bam =~ s/\,$//g;
    my $build_string .= qq"
spladder build \\
  --annotation ${gtf_file} \\
  --bams ${all_bam} \\
  --outdir ${spladder_dir} \\
  --output-gff3 \\
  --parallel 8 --verbose --readlen $options->{read_length} \\
  2>${stderr} 1>${stdout}
";
    my $test_string = qq"
spladder test \\
  --conditionA $group_strings{$options->{numerator}} \\
  --conditionB $group_strings{$options->{denominator}} \\
  --outdir ${spladder_dir} \\
  --labelA $options->{numerator} \\
  --labelB $options->{denominator} \\
  --out-tag $options->{numerator}_vs_$options->{denominator} \\
  --diagnose-plots --plot-format png --verbose \\
  --parallel 8 --high-memory \\
  2>>${stderr} 1>>${stdout}
";
    my $comment = qq"## Invoke spladder.";
    my $jstring = qq!mkdir -p ${spladder_dir}
${build_string}
${test_string}
!;
    my $spladder_job = $class->Submit(
        comment => $comment,
        cpus => $options->{jcpu},
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => qq"${spladder_dir}/spladder.stdout",
        stdout => qq"${spladder_dir}/spladder.stdout",
        stderr => qq"${spladder_dir}/spladder.stderr",);
    return($spladder_job);
}


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
        error_rate => 0.1,
        jmem => 24,
        jprefix => '50',
        seed => 25,
        allowed_missing => 3,);

    ## Set the minimum overlap parameter with the desired seed minus
    ## the allowed missing number of nucleotides.
    my $min_overlap = $options->{seed} - $options->{allowed_missing};

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
        lpanamensis => ['AAGTATCAGTTTCTGTACTTTATTG'],
    };
    my $input_string = qq"";
    my @input_lst = ();
    if ($options->{input} =~ /$options->{delimiter}/) {
        $paired = 1;
        @input_lst = split(/$options->{delimiter}/, $options->{input});
        for my $i (@input_lst) {
            $input_string .= "${i} ";
        }
    }
    else {
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

    my $fwd_polya_string = qq"AAAAAAAAAAAAAAAAA";
    my $rev_polya_string = qq"TTTTTTTTTTTTTTTTT";
    my $polya_fivep_stderr = qq"$paths->{output_dir}/slsearch_polya_fivep_cutadapt.stderr";
    my $polya_fivep_stdout = qq"$paths->{output_dir}/slsearch_polya_fivep_cutadapt.stdout";
    my $polya_threep_stderr = qq"$paths->{output_dir}/slsearch_polya_threep_cutadapt.stderr";
    my $polya_threep_stdout = qq"$paths->{output_dir}/slsearch_polya_threep_cutadapt.stdout";
    my $sl_fivep_stderr = qq"$paths->{output_dir}/slsearch_sl_fivep_cutadapt.stderr";
    my $sl_fivep_stdout = qq"$paths->{output_dir}/slsearch_sl_fivep_cutadapt.stdout";
    my $sl_threep_stderr = qq"$paths->{output_dir}/slsearch_sl_threep_cutadapt.stderr";
    my $sl_threep_stdout = qq"$paths->{output_dir}/slsearch_sl_threep_cutadapt.stdout";

    my $fwd_cmd = qq"
cutadapt -a ${fwd_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_fwdorient_threepr1_r1.fastq.gz -p $paths->{output_dir}/fwdsl_fwdorient_threepr1_r2.fastq.gz \\
   --interleaved  ${input_string} \\
   2>${sl_threep_stderr} 1>${sl_threep_stdout}
cutadapt -g ${fwd_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_revorient_fivepr1_r1.fastq.gz -p $paths->{output_dir}/fwdsl_revorient_fivepr1_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>${sl_fivep_stderr} 1>${sl_fivep_stdout}
cutadapt -A ${fwd_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_fwdorient_threepr2_r1.fastq.gz -p $paths->{output_dir}/fwdsl_fwdorient_threepr2_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${sl_threep_stderr} 1>>${sl_threep_stdout}
cutadapt -G ${fwd_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_revorient_fivepr2_r1.fastq.gz -p $paths->{output_dir}/fwdsl_revorient_fivepr2_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${sl_fivep_stderr} 1>>${sl_fivep_stdout}

cutadapt -a ${fwd_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_fwdorient_threepr1_r1.fastq.gz -p $paths->{output_dir}/fwdpolya_fwdorient_threepr1_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>${polya_threep_stderr} 1>${polya_threep_stdout}
cutadapt -g ${fwd_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_revorient_fivepr1_r1.fastq.gz -p $paths->{output_dir}/fwdpolya_revorient_fivepr1_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>${polya_fivep_stderr} 1>${polya_fivep_stdout}
cutadapt -A ${fwd_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_fwdorient_threepr2_r1.fastq.gz -p $paths->{output_dir}/fwdpolya_fwdorient_threepr2_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${polya_threep_stderr} 1>>${polya_threep_stdout}
cutadapt -G ${fwd_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_revorient_fivepr2_r1.fastq.gz -p $paths->{output_dir}/fwdpolya_revorient_fivepr2_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${polya_fivep_stderr} 1>>${polya_fivep_stdout}

";
    my $rev_cmd = qq"
cutadapt -a ${rev_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revsl_fwdorient_threepr1_r1.fastq.gz -p $paths->{output_dir}/revsl_fwdorient_threepr1_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${sl_threep_stderr} 1>>${sl_threep_stdout}
cutadapt -g ${rev_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revsl_revorient_fivepr1_r1.fastq.gz -p $paths->{output_dir}/revsl_revorient_fivepr1_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${sl_fivep_stderr} 1>>${sl_fivep_stdout}
cutadapt -A ${rev_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revsl_fwdorient_threepr2_r1.fastq.gz -p $paths->{output_dir}/revsl_fwdorient_threepr2_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${sl_threep_stderr} 1>>${sl_threep_stdout}
cutadapt -G ${rev_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revsl_revorient_fivepr2_r1.fastq.gz -p $paths->{output_dir}/revsl_revorient_fivepr2_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${sl_fivep_stderr} 1>>${sl_fivep_stdout}

cutadapt -a ${rev_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revpolya_fwdorient_threepr1_r1.fastq.gz -p $paths->{output_dir}/revpolya_fwdorient_threepr1_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${polya_threep_stderr} 1>>${polya_threep_stdout}
cutadapt -g ${rev_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revpolya_revorient_fivepr1_r1.fastq.gz -p $paths->{output_dir}/revpolya_revorient_fivepr1_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${polya_fivep_stderr} 1>>${polya_fivep_stdout}
cutadapt -A ${rev_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revpolya_fwdorient_threepr2_r1.fastq.gz -p $paths->{output_dir}/revpolya_fwdorient_threepr2_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${polya_threep_stderr} 1>>${polya_threep_stdout}
cutadapt -G ${rev_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revpolya_revorient_fivepr2_r1.fastq.gz -p $paths->{output_dir}/revpolya_revorient_fivepr2_r2.fastq.gz \\
   --interleaved ${input_string} \\
   2>>${polya_fivep_stderr} 1>>${polya_fivep_stdout}
";
    unless ($paired) {
        $fwd_cmd = qq"
cutadapt -b ${fwd_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_fwdorient_r1.fastq.gz \\
   ${input_string} \\
   2>>${sl_fivep_stderr} 1>>${sl_fivep_stdout}
cutadapt -B ${fwd_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdsl_revorient_r1.fastq.gz \\
   ${input_string} \\
   2>>${sl_threep_stderr} 1>>${sl_threep_stdout}

cutadapt -b ${fwd_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_fwdorient_r1.fastq.gz \\
   ${input_string} \\
   2>>${polya_fivep_stderr} 1>>${polya_fivep_stdout}
cutadapt -B ${fwd_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/fwdpolya_revorient_r1.fastq.gz \\
   ${input_string} \\
   2>>${polya_threep_stderr} 1>>${polya_threep_stdout}
";
        $rev_cmd = qq"
cutadapt -b ${rev_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revsl_fwdorient_r1.fastq.gz \\
   ${input_string} \\
   2>>${sl_fivep_stderr} 1>>${sl_fivep_stdout}
cutadapt -B ${rev_sl_string} -e $options->{error_rate} -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revsl_revorient_r1.fastq.gz \\
   ${input_string} \\
   2>>${sl_threep_stderr} 1>>${sl_threep_stdout}

cutadapt -b ${rev_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revpolya_fwdorient_r1.fastq.gz \\
   ${input_string} \\
   2>>${polya_fivep_stderr} 1>>${polya_fivep_stdout}
cutadapt -B ${rev_polya_string} -b ${rev_polya_string} --no-indels -O ${min_overlap} --trimmed-only \\
   -o $paths->{output_dir}/revpolya_revorient_r1.fastq.gz \\
   ${input_string} \\
   2>>${polya_threep_stderr} 1>>${polya_threep_stdout}
";
    }
    my $comment = qq!## This cutadapt invocation seeks to extract SL/polyA reads.!;

    my $jstring = qq!mkdir -p $paths->{output_dir}
${fwd_cmd}
${rev_cmd}
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
        stderr => $sl_fivep_stderr,
        stdout => $sl_fivep_stdout,
        polya_fivep_stderr =>  $polya_fivep_stderr,
        polya_fivep_stdout =>  $polya_fivep_stdout,
        polya_threep_stderr =>  $polya_threep_stderr,
        polya_threep_stdout =>  $polya_threep_stdout,
        sl_fivep_stderr =>  $sl_fivep_stderr,
        sl_fivep_stdout =>  $sl_fivep_stdout,
        sl_threep_stderr =>  $sl_threep_stderr,
        sl_threep_stdout =>  $sl_threep_stdout,
    );

}

=head2 C<SLSeq_Recorder>

  Use the reads from a SLSeq experiment to define the 5' UTRs of parasite genes.
  This is taking ideas from 10.1038/s41598-017-03987-0 and hopefully improving
  the logic a bit.

  There is an important caveat: this code currently assumes/requires paired-end reads.

  The general idea: Read in the extant gene annotations and make a contig-keyed hash
  where each key is comprised of an array of gene features including the relevant
  information about the current gene and the 'next' gene (reading from beginning to end
  of each contig).  Then read each aligned read from a position-sorted bam file from hisat
  or whatever and keep only the reads which are: a) pointing in the same direction as the ORF
  b) positioned in the 5' UTR of an ORF (thus the information about the next gene).  Given that,
  record every SL position with respect to the AUG and count how many reads are observed at every
  relative position.  Upon completion, write up a new gff file with new features of type 'SLSeq'
  that have a new start for the + and new end for the - strand features.

=over

=item C<Arguments>

  input(required): Sorted/indexed bam from hisat/bowtie/bwa/etc
  species: Find the original genome/gff annotations with this.
  gff_type(protein_coding_gene): Use this feature type to create new features.
  output_type(SLSeq_gene): Create new features of this type.
  id_tag(ID): Identify genes using this feature tag.

=back

=cut
sub SLSeq_Recorder {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        gff_type => 'protein_coding_gene',
        output_type => 'SLSeq_gene',
        jprefix => '50',
        id_tag => 'ID',);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $output_dir = $paths->{output_dir};
    make_path($output_dir) unless (-d $output_dir);
    my $output_gff = qq"${output_dir}/$options->{species}_maximum_utr.gff";
    my $output_tsv = qq"${output_dir}/$options->{species}_recorded_sl.tsv";
    my $jname = qq"$options->{jprefix}SLSeq_$options->{species}";
    my $jstring = qq!use Bio::Adventure::Splicing;
my \$result = \$h->Bio::Adventure::Splicing::SLSeq_Recorder_Worker(
  gff_type => '$options->{gff_type}',
  id_tag => '$options->{id_tag}',
  input => '$options->{input}',
  jname => '${jname}',
  output => '${output_gff}',
  output_tsv => '${output_tsv}',
  output_dir => '${output_dir}',
  output_type => '$options->{output_type}',
  species => '$options->{species}',);
!;
    my $comment = qq!## Seek SL positions.!;
    my $slseq = $class->Submit(
        comment => $comment,
        input => $options->{input},
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $output_gff,
        output_dir => $output_dir,
        output_tsv => $output_tsv,
        language => 'perl',);
    return($slseq);
}

sub SLSeq_Recorder_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        gff_type => 'protein_coding_gene',
        output_type => 'SLSeq_5pUTR',
        gff_tag => 'ID',);
    my $input_gff = qq"$options->{libpath}/$options->{libtype}/gff/$options->{species}.gff";
    my $fasta = qq"$options->{libpath}/$options->{libtype}/fasta/$options->{species}.fasta";
    unless (-r $options->{input}) {
        die("Unable to find input bam alignment.");
    }
    my $output_dir = $options->{output_dir};
    make_path($output_dir) unless (-d $output_dir);
    ## Set up the output filehandles:
    ##  The output csv file.
    my $output_csv = FileHandle->new(qq">${output_dir}/$options->{species}_recorded_sl.csv");
    print $output_csv qq"$options->{gff_tag}\tstart\tend\tstrand\trelative_pos\tnum_reads\n";
    ## The output gff file, this will receive individual gff strings:
    my $gff_output_filename = qq"${output_dir}/$options->{species}_modified_genes.gff";
    my $gff_out = FileHandle->new(">${gff_output_filename}");
    ## Use Bio::Tools::GFF to create the gff strings.
    my $gffout_io = Bio::Tools::GFF->new(-noparse => 1, -gff_version => 3);
    ## The log file:
    my $log_file = qq"${output_dir}/slseq_recorder.log";
    my $log = FileHandle->new(">${log_file}");
    print "Starting SLSeq Recorder, reading $options->{gff_tag} tags from $options->{gff_type}
entries of ${input_gff}.
Writing results to ${gff_output_filename}.\n";
    my ($counters, $observed) = $class->Bio::Adventure::Parsers::Count_Extract_GFF(
        gff => $input_gff, log => $log,
        gff_type => $options->{gff_type},
        gff_tag => $options->{gff_tag},);
    print $log "Beginning to read alignments from: $options->{input}.\n";
    my $sam = Bio::DB::Sam->new(-bam => $options->{input},
                                -fasta => ${fasta},);
    my @targets = $sam->seq_ids;
    my $num = scalar(@targets);
    my $bam = Bio::DB::Bam->open($options->{input});
    my $header = $bam->header;
    my $target_count = $header->n_targets;
    my $target_names = $header->target_name;
    my $align_count = 0;
    my $million_aligns = 0;
    my $alignstats = qx"samtools idxstats $options->{input}";
    my @alignfun = split(/\n/, $alignstats);
    my @aligns = split(/\t/, $alignfun[0]);
    my @unaligns = split(/\t/, $alignfun[1]);
    my $number_reads = $aligns[2] + $unaligns[3];
    my $output_name = qq"$options->{input}.txt";
    my $summary_string = qq"There are ${number_reads} alignments in $options->{input}.
There are $aligns[2] aligned reads and $unaligns[3] unaligned reads.\n";
    print $summary_string;
    my $new_contig = 'start';
    my $previous_contig = 'start';
  BAMLOOP: while (my $align = $bam->read1) {
        $counters->{reads}++;
        ## last BAMLOOP if ($counters->{quantified_reads} > 10000);
        my $read_seqid = $target_names->[$align->tid];
        ## my $start = $align->pos + 1;
        ## my $end = $align->calend;
        my $read_start = $align->pos;
        my $read_end = $align->calend - 1;
        my $read_strand = $align->strand;
        if (!defined($read_strand)) {
            $counters->{unstranded_reads}++;
        } elsif ($read_strand eq '1') {
            $counters->{plus_reads}++;
        } elsif ($read_strand eq '-1') {
            $counters->{minus_reads}++;
        } else {
            $counters->{notplusminus_reads}++;
        }
        ##my $read_seq = $align->query->dna;
        ##my $read_qual = $align->qual;
        ##my $read_pairedp = $align->paired;
        my $read_unmappedp = $align->unmapped;
        if ($read_unmappedp) {
            $counters->{unmapped_reads}++;
        }
        my $read_mate_unmappedp = $align->munmapped;
        ##my $read_reversedp = $align->reversed;
        ##my $read_mate_reversedp = $align->mreversed;
        ##my $read_mate_id = $align->mtid;
        ##my $read_mate_start = $align->mpos;
        ##my $read_properp = $align->proper_pair;
        if ($read_unmappedp && $read_mate_unmappedp) {
            $counters->{unmapped_pairs}++;
            next BAMLOOP;
        }
        my @tags = $align->get_all_tags;
        ##my $score_ref = $align->qscore;
        ##my $cigar = $align->cigar_str;
        my $not_primaryp = $align->get_tag_values('NOT_PRIMARY');
        unless ($not_primaryp) {
            $counters->{primary_reads}++;
        }
        my $supplementalp = $align->get_tag_values('SUPPLEMENTARY');
        if ($supplementalp) {
            $counters->{supplemental_reads}++;
            next BAMLOOP;
        }
        if ($not_primaryp) {
            $counters->{secondary_reads}++;
            next BAMLOOP;
        }
        if ($read_seqid ne $new_contig) {
            print "Starting new contig: ${read_seqid}\n";
            $previous_contig = $new_contig;
            if ($previous_contig ne 'start') {
                print "Writing data for ${previous_contig}\n";
                Write_SLData(observed => $observed,
                             contig => $previous_contig,
                             output_csv => $output_csv,
                             gff_tag => $options->{gff_tag},
                             gffio => $gffout_io,
                             gff_out => $gff_out,);
            }
            $new_contig = $read_seqid;
        }
        ## If the read strand is +1, then I want it associated with a feature on the +1 with a start which is > this read start.
        ## If the read strand is -1, then I want it associated with a feature on the -1 with a end which is < this read end.
        my $num_features = undef;
        next BAMLOOP unless (defined($observed->{$read_seqid}}));
        my $type = ref($observed->{$read_seqid});
        my $test_features = scalar(@{$observed->{$read_seqid}}) - 1;
        if (defined($test_features)) {
            unless (defined($num_features)) {
                $num_features = $test_features;
            }
        } else {
            $counters->{nogene_contigs}++;
            print "There appear to be no features for ${read_seqid}\n";
            next BAMLOOP;
        }
        my $checker = 0;
      FEAT_SEARCH: for my $c (0 .. $num_features) {
            my $checker++;
            my $current_feat = $observed->{$read_seqid}[$c];
            my $current_tag = $current_feat->{primary_tag};
            my $current_start = $current_feat->{start};
            my $current_end = $current_feat->{end};
            my $current_strand = $current_feat->{strand};
            my $next_distance = $current_feat->{next_distance};
            my $next_start = $current_feat->{next_start};
            my $next_strand = $current_feat->{next_strand};
            my $current_name = $current_feat->{gene_name};
            last FEAT_SEARCH if ($checker > $num_features);
            if (!defined($current_strand)) {
                print $log "In feature search, somehow strand is undefined for: ${current_name}\n";
                next FEAT_SEARCH;
            }
            if ($current_strand eq $read_strand &&
                $current_start <= $read_start &&
                $current_end > $read_start) {
                $counters->{same_strand_in_gene}++;
                next BAMLOOP;
            }
            ## Handle the (rare) case where reads are after a + strand gene and
            ## before a - strand gene
            if ($current_strand eq '1' &&
                $current_end < $read_start &&
                $next_strand eq '-1' &&
                $next_start > $read_start) {
                print "Contig:${read_seqid} Str:${current_strand} St:${read_start} vs. last end:${current_end} next start:${next_start}, unannotated ORF?\n";
                $counters->{unannotated_read}++;
                next BAMLOOP;
            }
            ## Sadly, I do not think I can exclude the opposite orientation,
            ## if the current feature is -1 and the next is +1, then the current reads
            ## may be associated with both

            next FEAT_SEARCH if (!defined($current_strand));

            if ($current_strand ne $read_strand &&
                $next_strand ne $read_strand) {
                $counters->{opposite_strand}++;
                next BAMLOOP;
            }
            if ($read_strand eq '1') {
                my $relative_position = $current_start - $read_start;
                next FEAT_SEARCH if ($relative_position >= $next_distance);
                if ($current_start > $read_start) {
                    $counters->{quantified_reads}++;
                    ## print "Got a +1 hit: $counters->{quantified_reads}\n";
                    if (defined($observed->{$read_seqid}[$c]->{observed}->{$relative_position})) {
                        $observed->{$read_seqid}[$c]->{observed}->{$relative_position}++;
                    } else {
                        $observed->{$read_seqid}[$c]->{observed}->{$relative_position} = 1;
                    }
                    if ($relative_position > $observed->{$read_seqid}[$c]->{most_upstream}) {
                        $observed->{$read_seqid}[$c]->{most_upstream} = $relative_position;
                    }
                    last FEAT_SEARCH;
                } else {
                    next FEAT_SEARCH;
                }
            } elsif ($read_strand eq '-1') {
                my $relative_position = $read_end - $current_end;
                next FEAT_SEARCH if ($relative_position >= $next_distance);
                if ($current_end < $read_end) {
                    $counters->{quantified_reads}++;
                    ## print "Got a -1 hit: $counters->{quantified_reads}\n";
                    if (defined($observed->{$read_seqid}[$c]->{observed}->{$relative_position})) {
                        $observed->{$read_seqid}[$c]->{observed}->{$relative_position}++;
                    } else {
                        $observed->{$read_seqid}[$c]->{observed}->{$relative_position} = 1;
                    }
                    if ($relative_position > $observed->{$read_seqid}[$c]->{most_upstream}) {
                        $observed->{$read_seqid}[$c]->{most_upstream} = $relative_position;
                    }
                    last FEAT_SEARCH;
                } else {
                    next FEAT_SEARCH;
                }
            } else { ## End checking the - strand
                print "What what? Neither +1 nor -1: $read_strand\n";
                next BAMLOOP;
            }
        } ## End iterating over every feature
    } ## End iterating over the alignments
    print "Writing data for ${new_contig}\n";
    Write_SLData(observed => $observed,
                 contig => $new_contig,
                 output_csv => $output_csv,
                 gff_tag => $options->{gff_tag},
                 gffio => $gffout_io,
                 gff_out => $gff_out,);
    use Data::Dumper;
    print Dumper $counters;
    $gff_out->close();
    $log->close();
    return($observed);
}

sub Write_SLData {
    my %args = @_;
    my $observed = $args{observed};
    my $contig = $args{contig};
    my $output_csv = $args{output_csv};
    my $gff_tag = $args{gff_tag};
    my $gffio = $args{gffio};
    my $gff_out = $args{gff_out};
    my $observed_genes = scalar(@{$observed->{$contig}}) - 1;
    for my $c (0 .. $observed_genes) {
        my $feat_info = $observed->{$contig}[$c];
        my $observations = $feat_info->{observed};
        my @observed_positions = keys %{$observations};
        my $gene_name = $feat_info->{gene_name};
        my $gene_start = $feat_info->{start};
        my $new_start = $gene_start;
        my $gene_end = $feat_info->{end};
        my $new_end = $gene_end;
        my $gene_strand = $feat_info->{strand};
        my $farthest = 0;
        my $num_positions = scalar(@observed_positions);
        if ($num_positions > 0) {
            for my $obs (sort {$a <=> $b } @observed_positions) {
                $farthest = $obs;
                print $output_csv qq"${gene_name}\t${gene_start}\t${gene_end}\t${gene_strand}\t${obs}\t$observations->{$obs}\n";
            }
            if ($gene_strand eq '1') {
                $new_start = $gene_start - $farthest;
                $new_end = $gene_start;
            } elsif ($gene_strand eq '-1') {
                $new_start = $gene_end;
                $new_end = $gene_end + $farthest;
            } else {
                die("This should not happen, neither +1/-1.\n");
            }
            if ($new_start < 0) {
                print "This is a failed entry: $contig $gene_name with start: $new_start\n";
                $new_start = 0;
            }
            my $standardized = Bio::SeqFeature::Generic->new(
                -primary_tag => 'SLSeq_gene',
                -display_name => $gene_name,
                -seq_id => $contig,
                -start => $new_start,
                -end => $new_end,
                -strand => $gene_strand,
                -score => 0,
                -tag => {
                    $gff_tag => $gene_name,
                    inference => 'SLSeq',
                    old_start => $gene_start,
                    old_end => $gene_end,
                });
            my $gff_string = $gffio->gff_string($standardized);
            print $gff_out "${gff_string}\n";
        }
    }
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
        found => 0,         ## The total number of observed SL.
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
        my $reader = Bio::Adventure::Get_FH(input => $i);
        my $seqio = Bio::SeqIO->new(-format => $format, -fh => $reader);
      FSA: while (my $seq = $seqio->next_seq) {
            $ind_search_result{searched}++;
            $global_search_result{searched}++;
            my $seq_id = $seq->id;
            my $read_seq = $seq->seq;
            my $fwd_end = undef;
            if ($read_seq =~ m/$sl_fwd_search/g) {
                $ind_search_result{sl_fwd_found}++;
                $fwd_end = pos($read_seq);
                my $new_read = $read_seq;
                $new_read =~ s/^.*$sl_fwd_search//g;
                print $output_sl_reads ">${seq_id} fwd
$new_read\n";
            }
            elsif ($read_seq =~ m/$polya_fwd_search/g) {
                $ind_search_result{polya_fwd_found}++;
                $fwd_end = pos($read_seq);
                my $new_read = $read_seq;
                $new_read =~ s/^.*$polya_fwd_search//g;
                print $output_polya_reads ">${seq_id} fwd
$new_read\n";
            }
            my $rc_end = undef;
            if ($read_seq =~ m/^.*$sl_rc_search/g) {
                $ind_search_result{sl_rev_found}++;
                $rc_end = pos($read_seq);
                my $new_read = reverse($read_seq);
                $new_read =~ tr/AGCTU/tcgaa/;
                $new_read =~ s/^.*$sl_fwd_search//g;
                print $output_sl_reads ">${seq_id} rev
${new_read}\n";
            }
            elsif ($read_seq =~ m/^.*$polya_rc_search/g) {
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
        }                   ## End reading the input fastq/fasta file.

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
