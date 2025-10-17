package Bio::Adventure::Map;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use File::Basename;
use File::Copy qw"cp";
use File::Spec;
use File::Which qw"which";

=head1 NAME

 Bio::Adventure::Map - Perform highthroughput sequence alignments with tools like bowtie/tophat/etc

=head1 SYNOPSIS

 use Bio::Adventure;
 my $hpgl = new Bio::Adventure;
 $hpgl->Bowtie();

=head1 METHODS

=head2 C<Bowtie>

 Run the OG (except bwa) short read aligner.
 10.1186/gb-2009-10-3-r25

 Perform a bowtie alignment.  Unless instructed otherwise, it will do so with 0
 mismatches and with multi-matches randomly placed 1 time among the
 possibilities (options -v 0 -M 1).  It will then convert the resulting sam alignment
 to a sorted-compressed-indexed bam, count it with htseq-count, compress the various
 output fastq files, and collect a few alignment statistics.

=over

 This requires the arguments: 'input' and 'species'.  The input is likely a
 colon-separated pair of (compressed)fastq files.  The species will be used
 to look for bowtie indexes in ${libdir}/${libtype}/indexes/${species}.

 The argument bt_type(v0M1: no mismatches, randomly place multi-matches in 1 location)
 defines the mismatch and multimatch parameters; I
 pre-defined a few likely option sets for these rather important options.

 The count(1: e.g. yes) argument defines whether htseq-count will be performed.

 libtype(genome: as opposed to rRNA or contaminants etc) defines the type of
 index to search against.

 gff_type(gene) defines the type of feature to count with htseq-count.  This is effectively
 the third column of a gff file.

 gff_tag(ID: ID is common for eukaryotic organisms, locus_tag is common for bacteria, most
 other species follow their own arbitrary rules) defines the ID type for htseq-count. These
 are the tags in the last column of a gff file.

 jprefix(10): Used all over the place to define the prefix job number.

 modules(bowtie1)

=back

=cut
sub Bowtie {
    my ($class, %args) = @_;
    say('Recall that you can change the bowtie arguments via "bt_type".');
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        bt_type => 'v0M1',
        count => 1,
        libtype => 'genome',
        jmem => 12,
        jprefix => '10',);
    my $start_species = $options->{species};
    my $species = $start_species;
    my $paths = $class->Bio::Adventure::Config::Get_Paths(species => $species,
                                                          bt_type => $options->{bt_type});
    if ($species =~ /\:/) {
        my @species_lst = split(/:/, $species);
        my @result_lst = ();
        foreach my $sp (@species_lst) {
            print "Invoking bowtie on ${sp}\n";
            $class->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::Bowtie(%{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }
    my $bt_type = $options->{bt_type};
    my $bt_args = $options->{bt_args}->{$bt_type};
    $bt_args = ' --best -v 0 -M 1 ' if (!defined($bt_args));
    my $bt_input = $options->{input};
    my $stranded = $options->{stranded};
    my $paired = 0;
    my $test_file = "";
    if ($bt_input =~ /$options->{delimiter}/) {
        $paired = 1;
        my @pair_listing = split(/$options->{delimiter}/, $bt_input);
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        my $r1_fd = $class->Get_FD(input => $pair_listing[0]);
        my $r2_fd = $class->Get_FD(input => $pair_listing[1]);
        $bt_input = qq" ${r1_fd} ${r2_fd} ";
        $test_file = $pair_listing[0];
    } else {
        $test_file = File::Spec->rel2abs($bt_input);
        my $r1_fd = $class->Get_FD(input => $test_file);
        $bt_input = qq" ${r1_fd} ";
        $stranded = 'no';
    }

    my $jname = qq"bt${bt_type}_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $libtype = $options->{libtype};

    my $count = $options->{count};
    my $bt_dir = qq"outputs/$options->{jprefix}bowtie_${species}";
    $bt_dir = $options->{bt_dir} if ($options->{bt_dir});

    ## Check that the indexes exist
    my $bt_reftest = $paths->{index_file};
    my $current_prefix = 10;
    if (defined($options->{jprefix})) {
        $current_prefix = $options->{jprefix};
    }
    my $index_job;
    if (!-r $bt_reftest && !$options->{bt1_indexjobs}) {
        $class->{bt1_indexjobs} = 1;
        my $genome_input = qq"$options->{libpath}/$options->{libtype}/${species}.fasta";
        $index_job = $class->Bio::Adventure::Index::BT1_Index(
            input => $paths->{fasta},
            jdepends => $options->{jdepends},
            jprefix => $current_prefix - 1,
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $bowtie_input_flag = "-q"; ## fastq by default
    $bowtie_input_flag = "-f" if ($options->{input} =~ /\.fasta/);
    my $comment = qq!## This is a bowtie1 alignment of ${bt_input} against
## $paths->{index} using arguments: ${bt_args}.!;
    my $aligned_filename = qq"${bt_dir}/$options->{jbasename}-${bt_type}_aligned_${species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/$options->{jbasename}-${bt_type}_unaligned_${species}.fastq";
    my $sam_filename = qq"${bt_dir}/$options->{jbasename}-${bt_type}.sam";
    my $jstring = qq!mkdir -p ${bt_dir}
bowtie \\
  $paths->{index} \\
  ${bt_args} \\
  -p $options->{jcpu} \\
  ${bowtie_input_flag} ${bt_input} \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>$paths->{stderr} \\
  1>$paths->{stdout}
!;

    my $bt_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $sam_filename,
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        stderr => $paths->{stderr},
        stdout => $paths->{stdout},
        unaligned => $unaligned_filename,);
    if (defined($index_job)) {
        $bt_job->{index} = $index_job;
    }

    my $expected_hours = $class->Bio::Adventure::Config::Estimate_Time_Composite(process => 'bt1');
    my $time_string = qq"${expected_hours}:00:00";
    my $compress_files = qq"${bt_dir}/$options->{jbasename}-${bt_type}_unaligned_${species}.fastq:${bt_dir}/$options->{jbasename}-${bt_type}_aligned_${species}.fastq";
    my $comp = $class->Bio::Adventure::Compress::Compress(
        comment => '## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt_args}.',
        jdepends => $bt_job->{job_id},
        jname => qq"${jname}_xzun",
        jprefix => $options->{jprefix} + 1,
        jwalltime => $time_string,
        input => $compress_files,);
    $bt_job->{compression} = $comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    my $trim_output_file = qq"outputs/trimomatic_stats.csv";
    my $sam_job = $class->Bio::Adventure::Convert::Samtools(
        jcpu => 1,
        input => $sam_filename,
        jdepends => $bt_job->{job_id},
        jname => qq"s2b_${jname}",
        jprefix => $options->{jprefix} + 3,
        paired => $paired,
        species => $options->{species},);
    $bt_job->{samtools} = $sam_job;

    my $htmulti;
    if ($count) {
        if ($libtype eq 'rRNA') {
            $htmulti = $class->Bio::Adventure::Count::HTSeq(
                input => $sam_job->{output},
                jdepends => $sam_job->{job_id},
                jname => qq"ht_${jname}",
                jprefix => $options->{jprefix} + 4,
                libtype => $libtype,
                mapper => 'bowtie1',
                stranded => $stranded,
                suffix => $bt_type,);
            $bt_job->{rRNA_count} = $htmulti;
        } else {
            $htmulti = $class->Bio::Adventure::Count::HT_Multi(
                input => $sam_job->{output},
                jdepends => $sam_job->{job_id},
                jname => qq"ht_${jname}",
                jprefix => $options->{jprefix} + 4,
                libtype => $libtype,
                mapper => 'bowtie1',
                stranded => $stranded,
                suffix => $bt_type,);
            $bt_job->{htseq} = $htmulti;
        }
    }  ## End if ($count)

    my $stats = $class->Bio::Adventure::Metadata::BT1_Stats(
        input => $paths->{stderr},
        bt_type => $bt_type,
        count_table => $bt_job->{htseq}->[0]->{output},
        jdepends => $bt_job->{job_id},
        jname => qq"${jname}_stats",
        jprefix => $options->{jprefix} + 2,
        trim_input => ${trim_output_file},);
    $bt_job->{stats} = $stats;
    $bt_job->{job_id} = $stats->{job_id};
    return($bt_job);
}

=head2 C<Bowtie2>

 Perform a bowtie2 alignment.

 This is pretty much a twin to the Bowtie function above with the
 obvious caveat that it uses bowtie2. It converts the resulting sam
 alignment to a sorted-compressed-indexed bam, count it with
 htseq-count, compress the various output fastq files, and collect a
 few alignment statistics.

=over

 input(required): colon-separated pair of (compressed)fastq files.
 species(required): used to look for bowtie indexes in
  ${libdir}/${libtype}/indexes/${species}.
 bt_type(v0M1: no mismatches, randomly place multi-matches in 1 location):
  defines the mismatch and multimatch parameters; I pre-defined a few
  likely option sets for these rather important options.
 count(1: e.g. yes): argument defines whether htseq-count will be performed.
 libtype(genome: as opposed to rRNA or contaminants etc): defines the type of
  index to search against.
 gff_type(gene): defines the type of feature to count with
  htseq-count.  This is effectively the third column of a gff file.
 gff_tag(ID: ID is common for eukaryotic organisms, locus_tag is common for bacteria, most
  other species follow their own arbitrary rules): defines the ID type for htseq-count. These
  are the tags in the last column of a gff file.
 jprefix(10): Prefix for the job/output directory
 modules(bowtie2): Load this environment module.

=back

=cut
sub Bowtie2 {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        count => 1,
        jmem => 28,
        jprefix => '20',);
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking bowtie2 on ${sp}\n";
            $class->{species} = $sp;
            my $result = Bio::Adventure::Map::Bowtie2($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $ready = $class->Check_Input(files => $options->{input},);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my %bt_jobs = ();
    my $libtype = 'genome';
    my $bt2_args = $options->{bt2_args};
    my $prefix_name = qw"bt2";
    my $bt2_name = qq"${prefix_name}_$options->{species}";
    my $suffix_name = $prefix_name;
    if ($options->{jname}) {
        $bt2_name .= qq"_$options->{jname}";
        $suffix_name .= qq"_$options->{jname}";
    }

    my $bt_dir = $paths->{output_dir};
    if ($args{bt_dir}) {
        $bt_dir = $args{bt_dir};
    }
    my $stranded = $options->{stranded};
    my $bt_input = $options->{input};
    my $test_file = '';
    if ($bt_input =~ /$options->{delimiter}/) {
        my @pair_listing = split(/$options->{delimiter}/, $bt_input);
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        my $r1_fd = $class->Get_FD(input => $pair_listing[0]);
        my $r2_fd = $class->Get_FD(input => $pair_listing[1]);
        $bt_input = qq" -1 ${r1_fd} -2 ${r2_fd} ";
        $test_file = $pair_listing[0];
    } else {
        $test_file = File::Spec->rel2abs($bt_input);
        my $r1_fd = $class->Get_FD(input => $test_file);
        $bt_input = qq" ${r1_fd} ";
        $stranded = 'no';
    }

    ## Check that the indexes exist
    my $job_suffix = 0;
    my $bt_reflib = $paths->{index_shell};
    my $bt_reftest = $paths->{index_file_shell};
    my $bt_reftest_large = $paths->{index_file2_shell};
    my $index_job;
    my $cluster = $options->{cluster};
    if (!-r $bt_reftest && !-r $bt_reftest_large) {
        $job_suffix++;
        print "Hey! The Indexes do not appear to exist, check this out: ${bt_reftest}\n";
        sleep(20);
        my $genome_input = $paths->{fasta};
        $index_job = $class->Bio::Adventure::Index::BT2_Index(
            input => $genome_input,
            cluster => 0,
            jdepends => $options->{jdepends},
            jprefix => qq"$options->{jprefix}_${job_suffix}",
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }
    my $bowtie_input_flag = '-q '; ## fastq by default
    $bowtie_input_flag = '-f ' if (${bt_input} =~ /\.fasta$/);

    my $jcpu = $options->{jcpu};
    my $stderr = qq"${bt_dir}/$options->{jbasename}.stderr";
    my $stdout = qq"${bt_dir}/$options->{jbasename}.stdout";
    my $comment = qq!## This is a bowtie2 alignment of ${bt_input} against
## ${bt_reflib} using arguments: ${bt2_args}.
!;
    my $aligned_filename = qq"${bt_dir}/$options->{jbasename}_aligned_$options->{species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/$options->{jbasename}_unaligned_$options->{species}.fastq";
    my $sam_filename = qq"${bt_dir}/$options->{jbasename}.sam";
    my $jstring = qq!mkdir -p ${bt_dir}
bowtie2 -x ${bt_reflib} ${bt2_args} \\
    -p ${jcpu} \\
    ${bowtie_input_flag} ${bt_input} \\
    --un ${unaligned_filename} \\
    --al ${aligned_filename} \\
    -S ${sam_filename} \\
    2>${stderr} \\
    1>${stdout}
!;
    my $bt2_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        cluster => $cluster,
        input => $bt_input,
        jname => $bt2_name,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => $options->{jmem},
        output => $sam_filename,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stderr => $stderr,
        stdout => $stdout,
        unaligned => $unaligned_filename,);
    my $compression_files = qq"${bt_dir}/$options->{jbasename}_unaligned_$options->{species}.fastq:${bt_dir}/$options->{jbasename}_aligned_$options->{species}.fastq";
    $job_suffix++;
    my $comp = $class->Bio::Adventure::Compress::Recompress(
        comment => '## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt2_args}.',
        cluster => $cluster,
        input => $compression_files,
        jdepends => $bt2_job->{job_id},
        jname => qq"xzun_${suffix_name}",
        jprefix => qq"$options->{jprefix}_${job_suffix}",);
    $bt2_job->{compression} = $comp;
    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    ## my $trim_output_file = qq"outputs/$options->{jbasename}-trimomatic.out";
    $job_suffix++;
    my $stats = $class->Bio::Adventure::Metadata::BT2_Stats(
        input => $stderr,
        cluster => $cluster,
        count_table => qq"$options->{jbasename}.count.xz",
        jdepends => $bt2_job->{job_id},
        jname => qq"bt2st_${suffix_name}",
        jprefix => qq"$options->{jprefix}_${job_suffix}",);
    $bt2_job->{stats} = $stats;
    $job_suffix++;
    my $sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $sam_filename,
        cluster => $cluster,
        jdepends => $bt2_job->{job_id},
        jname => qq"s2b_${suffix_name}",
        jprefix => qq"$options->{jprefix}_${job_suffix}",);
    $bt2_job->{samtools} = $sam_job;
    my $htseq_input = $sam_job->{output};
    my $htmulti;
    if ($options->{count}) {
        if ($libtype eq 'rRNA') {
            $job_suffix++;
            $htmulti = $class->Bio::Adventure::Count::HTSeq(
                cluster => $cluster,
                input => $sam_job->{output},
                jdepends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => qq"$options->{jprefix}_${job_suffix}",
                libtype => $libtype,
                mapper => 'bowtie2',
                stranded => $stranded,);
        } else {
            $job_suffix++;
            $htmulti = $class->Bio::Adventure::Count::HT_Multi(
                cluster => $cluster,
                input => $sam_job->{output},
                jdepends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => qq"$options->{jprefix}_${job_suffix}",
                libtype => $libtype,
                mapper => 'bowtie2',
                stranded => $stranded,);
            $bt2_job->{htseq} = $htmulti;
        }
    }
    my $last_value = $#$htmulti;
    my $last_job = $htmulti->[$last_value];
    $bt2_job->{job_id} = $last_job->{job_id};
    return($bt2_job);
}

=head2 C<Bowtie_RRNA>

 Search for rRNA using bowtie.

 This function is probably extraneous at this point.  It simply calls Bowtie()
 with the libtype set to 'rRNA' in order to get it to look for ribosomal reads
 instead of the default, genomic reads.

=cut
sub Bowtie_RRNA {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        exclude => 0,
        required => ["species"],);
    my $job_name = qq"rRNA_$options->{jbasename}";
    my $exclude = $options->{exclude};
    my $species = $options->{species};
    my ${bt_dir} = qq"outputs/bowtie_$options->{species}";
    my $job = $class->Bio::Adventure::Map::Bowtie(
        jdepends => $options->{jdepends},
        jname => $job_name,
        libtype => 'rRNA',
        prescript => $args{prescript},
        postscript => $args{postscript},);
    return($job);
}

=head2 C<BT_Multi>

 Run multiple bowtie parameters at once.

 Attempts to run multiple bowtie1 runs for a given species.  One run is performed
 for each of a few parameter sets which are kept in the global 'bt_args' variable.
 and generally include: 0 mismatch, 1 mismatch, 2 mismatches, 1 randomly placed
 hit, 0 randomly placed hits, or default options.

 This should either be removed or modified to work more generally with other tools.

=cut
sub BT_Multi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],);
    my $bt_input = $options->{input};
    my $species = $options->{species};
    my %bt_types = %{$options->{bt_args}};
    my $count = 0;
    my $job;
    foreach my $type (keys %bt_types) {
        my $jname = qq"bt${type}_${species}";
        my $bt_job = $class->Bio::Adventure::Map::Bowtie(
            input => $bt_input,
            bt_type => $type,
            jdepends => $options->{jdepends},
            jname => $jname,
            prescript => $args{prescript},
            postscript => $args{postscript},);
        if ($count == 0) {
            $job = $bt_job;
        } else {
            $job->{$count} = $bt_job;
        }
        $count++;
    }
    return($job);
}

=head2 C<BWA>

 The other OG short read aligner!
 10.1093/bioinformatics/btp324

 Perform a bwa alignment using both the sam(s|p)e and aln algorithms.  It then
 converts the output (when appropriate) to sorted/indexed bam and passes them to
 htseq.

=over

 input(required): likely a colon-separated pair of (compressed)fastq files.
 species(required): will be used to look for bwa indexes in
  ${libdir}/${libtype}/indexes/${species}.
 count(1: e.g. yes): argument defines whether htseq-count will be performed.
 libtype(genome: as opposed to rRNA or contaminants etc): defines the type of
  index to search against.
 gff_type(gene): defines the type of feature to count with
  htseq-count.  This is effectively the third column of a gff file.
 gff_tag(ID: ID is common for eukaryotic organisms, locus_tag is common for bacteria, most
  other species follow their own arbitrary rules): defines the ID type for htseq-count. These
  are the tags in the last column of a gff file.
 jprefix(30): Prefix for jobname and output directory.
 modules(bwa): load this environment module

=back

=cut
sub BWA {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jcpu => 8,
        jmem => 40,
        jprefix => 30,
        jwalltime => '72:00:00',
        required => ['input', 'species'],
        samtools_mapped => 1,
        samtools_unmapped => 1);
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking bwa on ${sp}\n";
            $options->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::BWA(%{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $bwa_input = $options->{input};
    my $stranded = $options->{stranded};
    my $test_file = '';
    my $forward_reads = '';
    my $reverse_reads = undef;
    if ($bwa_input =~ /$options->{delimiter}/) {
        my @pair_listing = split(/$options->{delimiter}/, $bwa_input);
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        $test_file = $pair_listing[0];
        if ($test_file =~ /\.bz2|\.xz/) {
            $forward_reads = $class->Get_FD(input => $pair_listing[0]);
            $reverse_reads = $class->Get_FD(input => $pair_listing[1]);
        } elsif ($test_file =~ /\.fastq|\.gz|\.fq/) {
            $forward_reads = qq"$pair_listing[0]";
            $reverse_reads = qq"$pair_listing[1]";
        } else {
            print "These inputs seem strange and may not work.\n";
            $forward_reads = qq"<(less $pair_listing[0])";
            $reverse_reads = qq"<(less $pair_listing[1])";
        }
        $bwa_input = qq" ${forward_reads} ${reverse_reads} ";
    } else {
        $stranded = 'no';
        $test_file = File::Spec->rel2abs($bwa_input);
        if ($test_file =~ /\.bz2|\.xz/) {
            my $reads = $class->Get_FD(input => $test_file);
            $bwa_input = qq" ${reads} ";
        } elsif ($test_file =~ /\.fastq|\.gz|\.fq/) {
            $bwa_input = qq" ${bwa_input} ";
        } else {
            print "These inputs seem strange and may not work.\n";
            $bwa_input = qq" <(less ${bwa_input} ";
        }
        $forward_reads = $bwa_input;
    }

    my $jname = qq"bwa_$options->{species}";
    $jname = $options->{jname} if ($options->{jname});

    my $bwa_dir = $paths->{output_dir};
    $bwa_dir = $options->{bwa_dir} if ($options->{bwa_dir});

    my $uncompress_jobid = undef;
    my $index_jobid = undef;
    ## Check that the indexes exist
    ## NOTE: BWA Might require the file to be named ...fa instead of ...fasta
    my $bwa_reflib = qq"$paths->{index_dir}/$options->{species}.fa";
    my $bwa_reftest = qq"${bwa_reflib}.sa";
    my $index_job;
    my $cluster = $options->{cluster};
    if (!-r $bwa_reftest) {
        my $index_extras = {
            input => $paths->{fasta_shell},
            jprefix => $options->{jprefix} - 1,
        };

        my %index_args = $class->Extra_Options(options => $options, extras => $index_extras, cluster => 0);
        $index_job = $class->Bio::Adventure::Index::BWA_Index(%index_args);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $method = $options->{bwa_method};
    my @bwa_methods = split(/[:;,]/, $options->{bwa_method});
    my $job_string = '';
    my $comment = '## A series of bwa commands.';
    my $sam_outs = '';
    my @sam_files = ();
    my $stderr = qq"${bwa_dir}/bwa.stderr";
    my $stdout = qq"${bwa_dir}/bwa.stdout";
    my $htmulti;

    for my $method (@bwa_methods) {
        my $sam_out = qq"${bwa_dir}/$options->{jbasename}_${method}.sam";
        push(@sam_files, $sam_out);
        $sam_outs = qq"${sam_out}:";
        my $extra_args = qq" $options->{arbitrary} -t $options->{jcpu} ";
        if ($method eq 'mem') {
            $extra_args = qq" -M ${extra_args} ";
        }

        if ($method eq 'mem') {
            $job_string .= qq!mkdir -p ${bwa_dir}
bwa ${method} \\
  ${bwa_reflib} \\
  ${bwa_input} \\
  -o ${sam_out} \\
  ${extra_args} \\
  1>${stdout} \\
  2>>${stderr}
test=\$?
if [[ "\${test}" -eq "0" ]]; then
  echo "bwa ${method} finished happily with return 0." >> ${bwa_dir}/bwa_${method}.stderr
else
  echo "bwa ${method} finished sadly with return \${test}." >> ${bwa_dir}/bwa_${method}.stderr
  exit \${test}
fi
!;
        } elsif ($method eq 'aln') {
            $job_string .= qq!mkdir -p ${bwa_dir}
bwa ${method} ${extra_args} \\
  ${bwa_reflib} \\
  ${forward_reads} \\
  1>${bwa_dir}/$options->{jbasename}_aln-r1.sai \\
  2>>${bwa_dir}/bwa_${method}_r1.stderr
test=\$?
if [[ "\${test}" -eq "0" ]]; then
  echo "bwa ${method} finished happily with return 0." >> ${bwa_dir}/bwa_${method}.stderr
else
  echo "bwa ${method} finished sadly with return \${test}." >> ${bwa_dir}/bwa_${method}.stderr
  exit \${test}
fi
!;
            if (defined($reverse_reads)) {
                $job_string .= qq!mkdir -p ${bwa_dir}
bwa ${method} ${extra_args} \\
  ${bwa_reflib} \\
  ${reverse_reads} \\
  1>${bwa_dir}/$options->{jbasename}_aln-r2.sai \\
  2>>${bwa_dir}/bwa_${method}_r1.stderr
test=\$?
if [[ "\${test}" -eq "0" ]]; then
  echo "bwa ${method} reverse finished happily with return 0." >> ${bwa_dir}/bwa_${method}.stderr
else
  echo "bwa ${method} reverse finished sadly with return \${test}." >> ${bwa_dir}/bwa_${method}.stderr
  exit \${test}
fi

bwa sampe \\
  ${bwa_reflib} \\
  ${bwa_dir}/$options->{jbasename}_aln-r1.sai \\
  ${bwa_dir}/$options->{jbasename}_aln-r2.sai \\
  ${bwa_input} \\
  1>${sam_out} \\
  2>>${bwa_dir}/bwa_sampe.stderr
test=\$?
if [[ "\${test}" -eq "0" ]]; then
  echo "bwa sampe finished happily with return 0." >> ${bwa_dir}/bwa_${method}.stderr
else
  echo "bwa sampe finished sadly with return \${test}." >> ${bwa_dir}/bwa_${method}.stderr
  exit \${test}
fi
!;
            } else {
                $job_string .= qq!
bwa samse \\
  ${bwa_reflib} \\
  ${bwa_dir}/$options->{jbasename}_aln-r1.sai \\
  ${bwa_input} \\
  1>${sam_out} \\
  2>>${bwa_dir}/bwa_samse.stderr
test=\$?
if [[ "\${test}" -eq "0" ]]; then
  echo "bwa ${method} reverse finished happily with return 0." >> ${bwa_dir}/bwa_${method}.stderr
else
  echo "bwa ${method} reverse finished sadly with return \${test}." >> ${bwa_dir}/bwa_${method}.stderr
  exit \${test}
fi

!;
            }
        } elsif ($method eq 'samse') {
                $job_string .= qq!
bwa samse ${extra_args} \\
  ${bwa_reflib} \\
  ${bwa_dir}/$options->{jbasename}_aln-r1.sai \\
  ${bwa_input} \\
  1>${sam_out} \\
  2>>${bwa_dir}/bwa_samse.stderr
test=\$?
if [[ "\${test}" -eq "0" ]]; then
  echo "bwa ${method} reverse finished happily with return 0." >> ${bwa_dir}/bwa_${method}.stderr
else
  echo "bwa ${method} reverse finished sadly with return \${test}." >> ${bwa_dir}/bwa_${method}.stderr
  exit \${test}
fi
!;
        } elsif ($method eq 'sampe') {
            $job_string .= qq!
bwa sampe ${extra_args} \\
  ${bwa_reflib} \\
  ${bwa_dir}/$options->{jbasename}_aln-r1.sai \\
  ${bwa_dir}/$options->{jbasename}_aln-r2.sai \\
  ${bwa_input} \\
  1>${sam_out} \\
  2>>${bwa_dir}/bwa_sampe.stderr
test=\$?
if [[ "\${test}" -eq "0" ]]; then
  echo "bwa ${method} reverse finished happily with return 0." >> ${bwa_dir}/bwa_${method}.stderr
else
  echo "bwa ${method} reverse finished sadly with return \${test}." >> ${bwa_dir}/bwa_${method}.stderr
  exit \${test}
fi
!;
        }
    } ## End iterating over potential bwa methods.
    $sam_outs =~ s/:$//g;
    my %submit_extras = (
        cluster => $cluster,
        comment => $comment,
        input => $bwa_input,
        jname => qq"bwa_$options->{species}",
        jstring => $job_string,
        stdout => $stdout,
        stderr => $stderr,);
    my $bwa_job = $class->Submit(%submit_extras);
    my @samtools_jobs = ();
    my $sam_count = 0;
    for my $sam (@sam_files) {
        my $sam_method = $bwa_methods[$sam_count];
        my %sam_extras = (
            input => $sam,
            jdepends => $bwa_job->{job_id},
            jname => qq"s2b_${sam_method}",
            jmem => '30',);
        my $sam_job = $class->Bio::Adventure::Convert::Samtools(%sam_extras);
        my $sam_name = qq"samtools_${sam_method}";
        my $htseq_name = qq"htseq_${sam_method}";
        $bwa_job->{$sam_name} = $sam_job;
        if ($options->{count}) {
            my $htseq_extras = {
                input => $sam_job->{output},
                jdepends => $sam_job->{job_id},
                jname => $htseq_name,
                mapper => 'bwa',
                stranded => $stranded,
            };
            my %htseq_args = $class->Extra_Options(options => $options, extras => $htseq_extras);
            $htmulti = $class->Bio::Adventure::Count::HTSeq(%htseq_args);
            $bwa_job->{$htseq_name} = $htmulti;
        }
        $sam_count++;
    } ## Done running samtools/htseq on the various bwa output files.
    $bwa_job->{job_id} = $htmulti->{job_id};
    return($bwa_job);
}

=head2 C<Downsample_Guess_Strand>

  Downsample the data and use salmon to guestimate the strandedness of the library.

=over

 input(required): Trimmed fastq files.
 species(required): Used to find the salmon indexes.
 reads(1e6): Pull out this many reads.

=back

=cut
sub Downsample_Guess_Strand {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        jmem => 24,
        reads => 1000000, ## Needs to be big enough for salmon to smaple
        jprefix => '45',);
    my $depends = $options->{jdepends};
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking salmon on ${sp}\n";
            $options->{species} = $sp;
            my $result = Bio::Adventure::Map::Salmon($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $ready = $class->Check_Input(files => $options->{input});
    my %sa_jobs = ();
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = $options->{species};

    my $jname = qq"sal_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $outdir = qq"outputs/$options->{jprefix}salmon_${species}";
    my $stderr = qq"${outdir}/salmon_${species}.stderr";
    my $stdout = qq"${outdir}/salmon_${species}.stdout";

    my $sa_args = '';
    my $sa_input = $options->{input};
    my $input_name = $sa_input;
    my $sa_rm = qq"rm ${outdir}/r1_downsampled.fastq
";
    my $downsample_commands = '';
    if ($sa_input =~ /$options->{delimiter}/) {
        my @pair_listing = split(/$options->{delimiter}/, $sa_input);
        $sa_args .= qq" -1 ${outdir}/r1_downsampled.fastq -2 ${outdir}/r2_downsampled.fastq ";
        $input_name = $pair_listing[0];
        $downsample_commands .= qq"sampled_r1=\$( { less $pair_listing[0] ||:; } |\\
  head -n $options->{reads} > ${outdir}/r1_downsampled.fastq )
sampled_r2=\$( { less $pair_listing[1] ||:; } |\\
  head -n $options->{reads} > ${outdir}/r2_downsampled.fastq )
";
        my $sa_rm .= qq"rm ${outdir}/r2_downsampled.fastq
";
    } else {
        $sa_args .= qq" -r ${outdir}/r1_downsampled.fastq ";
        $downsample_commands .= qq"sampled_r1=\$( { less $sa_input} ||:; } |\\
  head -n $options->{reads} > ${outdir}/r1_downsampled.fastq )
";
    }

    ## Check that the indexes exist
    my $sa_reflib = qq"$options->{libpath}/${libtype}/indexes/$options->{species}_salmon_index";
    my $index_job;
    if (!-r $sa_reflib) {
        my $transcript_file = qq"$options->{libpath}/${libtype}/$options->{species}_cds.fasta";
        $index_job = $class->Bio::Adventure::Index::Salmon_Index(
            input => $transcript_file,
            depends => $options->{jdepends},
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $comment = qq!## This is a salmon pseudoalignment of ${sa_input} against
## ${sa_reflib}.
!;
    my $jstring = qq!mkdir -p ${outdir}
## For reasons passing my understanding, the less downsample is throwing a ERR
## when running in the script; but the outputs are created correctly.
${downsample_commands}
salmon quant -i ${sa_reflib} \\
  -l A --gcBias --validateMappings  \\
  ${sa_args} \\
  -o ${outdir} \\
  2>${stderr} \\
  1>${stdout}
${sa_rm}
!;
    my $salmon = $class->Submit(
        comment => $comment,
        input => $sa_input,
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        output => qq"${outdir}/quant.sf",
        stderr => $stderr,
        stdout => $stdout,
        prescript => $args{prescript},
        postscript => $args{postscript},);
    $salmon->{index_job} = $index_job;
    my $strand = $class->Bio::Adventure::Metadata::Salmon_Strand(
        input => $stderr,
        jdepends => $salmon->{job_id},
        jname => qq"sastrand_$options->{species}",
        jprefix => $options->{jprefix},);
    $salmon->{strand} = $strand;
    $salmon->{job_id} = $strand->{job_id};
    return($salmon);
}

=head2 C<Hisat2>

 Invoke hisat2!
 10.1038/s41587-019-0201-4

 Hisat2 is currently my favorite aligner.

=over

 input(required): Colon separated fastq input files.
 species(required): Set the indexes.
 gff_type('gene'): Define the htseq-count type parameter.
 gff_tag('ID'): Define the htseq-count id attribute parameter.
 count(1): Count after aligning?
 libtype('genome'): Change this for different index types (contaminants/rRNA/genome).
 jmem(24): Expected memory required.
 jprefix('40'): Set the job prefix and output directory.
 modules('hisat2','samtools','htseq','bamtools'): Load these environment modules.

=back

=cut
sub Hisat2 {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        compress => 1,
        count => 1,
        counter => 'featureCounts',
        jcpu => 8,
        jmem => 48,
        jname => 'hisat',
        jprefix => '40',
        jwalltime => '36:00:00',
        get_insertsize => 1,
        maximum => undef,
        output_dir => undef,
        output_unaligned => undef,
        unaligned_discordant => undef,
        required => ['species', 'input',],
        samtools => 1,);
    if ($options->{species} =~ /:/) {
        my $start_species = $options->{species};
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        foreach my $sp (@species_lst) {
            print "Invoking hisat2 on ${sp}\n";
            $options->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::Hisat2(%{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    unless ($options->{introns}) {
        print "It appears that introns are disabled for this hisat run, is this correct?\n";
        sleep(5);
    }

    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    print "In Map::Hisat2 the output directory is: $paths->{output_dir}\n";
    print "Hisat output directory does not exist!\n" unless (-d $paths->{output_dir});
    my $ready;
    if (!$options->{jdepends}) {
        $ready = $class->Check_Input(files => $options->{input},);
    }
    my $hisat_args = $class->Passthrough_Args(arbitrary => $options->{hisat_args});
    $hisat_args .= qq" -k $options->{maximum} " if (defined($options->{maximum}));
    if ($options->{task} eq 'dnaseq' or !$options->{introns}) {
        $hisat_args .= qq" --no-spliced-alignment ";
    }
    my $prefix_name = 'hisat';
    $prefix_name .= qq"_k$options->{maximum}" if (defined($options->{maximum}));
    my $hisat_name = qq"${prefix_name}_$options->{species}_$options->{libtype}";
    my $xz_jname = qq"xz_${hisat_name}";
    my $sam_jname = qq"s2b_${hisat_name}";
    my $counter_jname = qq"counter_${hisat_name}";
    my $stat_jname = qq"${hisat_name}_stat";

    my $hisat_input = $options->{input};
    my $test_file = '';
    my $paired = 0;
    my $stranded = $options->{stranded};
    my @pair_listing;
    if ($hisat_input =~ /$options->{delimiter}/) {
        @pair_listing = split(/$options->{delimiter}/, $hisat_input);
        $paired = 1;
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        ## After years of working without problem, suddenly my lesspipe
        ## process substitution pre-filter for arbitrarily compressed files
        ## stopped working and might well give me an aneurysm trying to figure out.
        ## It turns out that there is a race condition somewhere which is triggered when
        ## less buffers its output -- so make sure the environment variable 'LESS'
        ## contains --unbuffered
        $hisat_input = qq" -1 $pair_listing[0] -2 $pair_listing[1] ";
        if ($pair_listing[0] =~ /\.[x|g|b]z$/) {
            ## It is noteworthy that I modified hisat2 on my computer
            ## so this is never necessary.
            ## my $r1_fd = $class->Get_FD(input => $pair_listing[0]);
            ## my $r2_fd = $class->Get_FD(input => $pair_listing[1]);
            $hisat_input = qq" -1 $pair_listing[0] -2 $pair_listing[1] ";
            ## $hisat_input = qq" -1 <(less $pair_listing[0]) -2 <(less $pair_listing[1]) ";
        }
        $test_file = $pair_listing[0];
    } else {
        $stranded = 'no';
        $options->{get_insertsize} = 0;
        $test_file = File::Spec->rel2abs($hisat_input);
        $hisat_input = qq" -U ${test_file} ";
        if ($test_file =~ /\.[x|g|b]z$/) {
            ## It is noteworthy that I modified hisat2 on my computer so this is never necessary.
            $hisat_input = qq" -U ${test_file} ";
        }
    }
    ## Check that the indexes exist
    my $cluster = $options->{cluster};
    if (!-r $paths->{index_file} && !-r $paths->{index_file2}) {
        my $genome_fasta = $paths->{fasta};
        my $index_job = $class->Bio::Adventure::Index::Hisat2_Index(
            cluster => 0,
            input => $genome_fasta,
            jprefix => qq"$options->{jprefix}_0",
            output_dir => $paths->{output_dir},
            jdepends => $options->{jdepends},
            libtype => $options->{libtype},);
        ## The following line inserts the indexer into the dependency chain
        ## Do not forget this, it is rather important.
    }

    my $hisat_input_flag = '-q '; ## fastq by default
    $hisat_input_flag = '-f ' if (${hisat_input} =~ /\.fasta$/);
    my $jcpu = $options->{jcpu};

    my $stderr = qq"$paths->{output_dir}/hisat_$options->{species}_$options->{libtype}";
    my $stdout = $stderr;
    my $metrics = $stderr;
    if (defined($options->{jbasename})) {
        $stderr .= "_$options->{jbasename}.stderr";
        $stdout .= "_$options->{jbasename}.stdout";
        $metrics .= "_$options->{jbasename}.metrics";
    } else {
        $stderr .= ".stderr";
        $stdout .= ".stdout";
        $metrics .= ".metrics";
    }
    my $comment = qq!## This is a hisat2 alignment of ${hisat_input} against $paths->{index_shell}
!;
    $comment .= qq"## This alignment is using arguments: ${hisat_args}.\n" unless ($hisat_args eq '');
    my $aligned_discordant_filename = qq"$paths->{output_dir}/aldis_$options->{species}_$options->{libtype}.fastq";
    my $unaligned_discordant_filename = qq"$paths->{output_dir}/unaldis_$options->{species}_$options->{libtype}.fastq";
    if (defined($options->{unaligned_discordant})) {
        $unaligned_discordant_filename = $options->{unaligned_discordant};
    }
    my $aligned_concordant_filename = qq"$paths->{output_dir}/alcon_$options->{species}_$options->{libtype}.fastq";
    my $unaligned_concordant_filename = qq"$paths->{output_dir}/unalcon_$options->{species}_$options->{libtype}.fastq";
    if (defined($options->{unaligned_output})) {
        $unaligned_concordant_filename = $options->{unaligned_output};
    }
    my $sam_filename = qq"$paths->{output_dir}/$options->{species}_$options->{libtype}.sam";
    my $jstring = qq!mkdir -p $paths->{output_dir}
mapped=1
{
  /usr/bin/time -v -o ${stdout}.time -a hisat2 -x $paths->{index_shell} ${hisat_args} \\
    -p ${jcpu} \\
    ${hisat_input_flag} ${hisat_input} \\
    --phred$options->{phred} \\
    --un ${unaligned_discordant_filename} \\
    --al ${aligned_discordant_filename} \\
!;
    ## Record the aligned/unaligned filenames both before and after compression.
    my $aligned_filenames = $aligned_discordant_filename;
    my $aligned_xz = qq"${aligned_filenames}.xz";
    my $unaligned_filenames = $unaligned_discordant_filename;
    my $unaligned_xz = qq"${unaligned_filenames}.xz";
    if ($paired) {
        ## The concordant flag tries to send reads to xxx.1 and xxx.2 even if
        ## the data is not paired.
        $jstring .= qq!    --un-conc ${unaligned_concordant_filename} \\
    --al-conc ${aligned_concordant_filename} \\
!;
        my $aligned_base = basename($aligned_concordant_filename, ('.fastq', '.fasta'));
        my $unaligned_base = basename($unaligned_concordant_filename, ('.fastq', '.fastq'));
        my $the_dirname = dirname($aligned_concordant_filename);
        ## For the set of all aligned filenames, there should now be 3, the original discordant and two more.
        ## For the xz filenames, there should be 2, just the concordant; because these will be used downstream.
        $aligned_filenames .= qq":${the_dirname}/${aligned_base}.1.fastq:${the_dirname}/${aligned_base}.2.fastq";
        $aligned_xz = qq"${the_dirname}/${aligned_base}.1.fastq.xz:${the_dirname}/${aligned_base}.2.fastq.xz";

        $unaligned_filenames = qq"${the_dirname}/${unaligned_base}.1.fastq:${the_dirname}/${unaligned_base}.2.fastq:${unaligned_discordant_filename}";
        $unaligned_xz = qq"${the_dirname}/${unaligned_base}.1.fastq.xz:${the_dirname}/${unaligned_base}.2.fastq.xz:${unaligned_discordant_filename}.xz";
    }
    my $all_filenames = qq"${aligned_filenames}:${unaligned_filenames}";
    $jstring .= qq!    -S ${sam_filename} \\
    --met-file ${metrics} \\
    2>${stderr} \\
    1>${stdout}
} || {
  echo 'hisat2 failed.' >> ${stderr}
}
## It seems like small files on NFS can cause a race condition leading to failure to properly
## read the sam/un(aligned) outputs.
sync
!;
    ## Explicitly setting the cluster state to whatever it was before checking for indexes
    ## Because I started setting the index creation in a bash job to avoid
    ## accidental overwriting of indexes by multiple jobs.
    my $hisat_job = $class->Submit(
        cluster => $cluster,
        comment => $comment,
        jdepends => $options->{jdepends},
        input => $hisat_input,
        jname => $hisat_name,
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => $options->{jmem},
        paired => $paired,
        output => $sam_filename,
        output_metrics => $metrics,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        aligned => $aligned_filenames,
        aligned_comp => $aligned_xz,
        stderr => $stderr,
        stdout => $stdout,
        unaligned => $unaligned_filenames,
        unaligned_dis => $unaligned_discordant_filename,
        unaligned_comp => $unaligned_xz,);
    my $new_jprefix = qq"$options->{jprefix}_1";
    ## HT1_Stats also reads the trimomatic output, which perhaps it should not.
    ## my $trim_output_file = qq"outputs/$options->{jbasename}-trimomatic.out";
    my $sam_jprefix = qq"$options->{jprefix}_2";
    my $sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $sam_filename,
        jcpu => 1,
        jdepends => $hisat_job->{job_id},
        jname => $sam_jname,
        jprefix => $sam_jprefix,
        paired => $paired,
        species => $options->{species},);
    $hisat_job->{samtools} = $sam_job;
    if ($options->{compress}) {
        my $comp = $class->Bio::Adventure::Compress::Compress(
            jname => $xz_jname,
            jprefix => $new_jprefix,
            input => qq"${aligned_filenames}:${unaligned_filenames}",
            jdepends => $sam_job->{job_id});
        $hisat_job->{compression} = $comp;
        ## Sneak the compression job's ID in place as hisat's.
        ## No, I changed my mind, let the samtools jump ahead.
        ## $hisat_job->{job_id} = $comp->{job_id};
    }

    if ($options->{get_insertsize}) {
        my $sizes = $class->Bio::Adventure::Count::Insert_Size(
            input => $sam_job->{output},
            jprefix => qq"$options->{jprefix}_3",
            jdepends => $sam_job->{job_id},);
        $hisat_job->{insertsizes} = $sizes;
    }
    if ($options->{duplicate}) {
        my $duplication = $class->Bio::Adventure::Convert::GATK_Dedup(
            input => $sam_job->{output},
            remove => $options->{remove},
            jprefix => qq"$options->{jprefix}_6",
            jdepends => $sam_job->{job_id},);
    }
    $new_jprefix = qq"$options->{jprefix}_3";
    my $counter_input;
    if ($paired) {
        $counter_input = $sam_job->{paired_output};
    } else {
        $counter_input = $sam_job->{output};
    }
    my $counter_out = '';
    my $counter;
    if ($options->{count} && $options->{counter} eq 'htseq') {
        if ($options->{libtype} eq 'rRNA') {
            $counter = $class->Bio::Adventure::Count::HTSeq(
                input => $counter_input,
                jdepends => $sam_job->{job_id},
                jname => $counter_jname,
                jprefix => $new_jprefix,
                libtype => $options->{libtype},
                mapper => 'hisat2',
                paired => $paired,
                stranded => $stranded,);
        } else {
            $counter = $class->Bio::Adventure::Count::HT_Multi(
                input => $counter_input,
                jdepends => $sam_job->{job_id},
                jname => $counter_jname,
                jprefix => $new_jprefix,
                libtype => $options->{libtype},
                mapper => 'hisat2',
                paired => $paired,
                stranded => $stranded,);
        }
        $hisat_job->{htseq} = $counter;
        if (ref($hisat_job->{htseq}) eq 'ARRAY') {
            $counter_out = $hisat_job->{htseq}->[0]->{output};
        } else {
            $counter_out = $hisat_job->{htseq}->{output};
        }
        ## If this is not the final job in a chain, then make sure
        ## that any jobs queued after it do not start until
        ## samtools/htseq/etc are finished.
    } elsif ($options->{count} && $options->{counter} eq 'featureCounts') {
        $counter = $class->Bio::Adventure::Count::Feature_Counts(
            input => $counter_input,
            jdepends => $sam_job->{job_id},
            jname => $counter_jname,
            jprefix => $new_jprefix,
            libtype => $options->{libtype},
            mapper => 'hisat2',
            paired => $paired,
            stranded => $stranded,);
        ## If this is not the final job in a chain, then make sure
        ## that any jobs queued after it do not start until
        ## samtools/htseq/etc are finished.
        $hisat_job->{job_id} = $counter->{job_id};
    }
    $hisat_job->{counter} = $counter;
    my $stats = $class->Bio::Adventure::Metadata::HT2_Stats(
        input => $stderr,
        count_table => $counter_out,
        jdepends => $hisat_job->{job_id},
        jname => $stat_jname,
        jprefix => $new_jprefix,
        output_dir => $paths->{output_dir},);
    $hisat_job->{stats} = $stats;
    $hisat_job->{job_id} = $stats->{job_id};
    return($hisat_job);
}

=head2 C<Kallisto>

 Perform a kallisto transcript quantification.
 10.1038/nbt.3519

 Kallisto and salmon are my two favorite 'voting' based aligners.

=over

 input(required): Colon separated fastq input file(s).
 species(required): Define the location of the indexes.
 jmem(24): Expected memory requirements.
 jprefix('46'): Set the output directory/job prefix.
 modules('kallisto'): Load this environment module.

=back

=cut
sub Kallisto {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jcpu => 8,
        jwalltime => '8:00:00',
        jmem => 24,
        jprefix => '46',
        required => ['input', 'species',],);
    die('No species provided.') unless($options->{species});
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking kallisto on ${sp}\n";
            $options->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::Kallisto(%{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my %ka_jobs = ();
    my $ka_depends_on = '';
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = $options->{species};
    my $jname = qq"kall_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $ka_args = '';
    my $ka_input = $options->{input};
    my $input_name = $ka_input;
    if ($ka_input =~ /$options->{delimiter}/) {
        my @pair_listing = split(/$options->{delimiter}/, $ka_input);
        $ka_args .= ' --bias ';
        if ($options->{stranded} != 0) {
            $ka_args .= " --$options->{stranded} ";
        }
        my $r1_fd = $class->Get_FD(input => $pair_listing[0]);
        my $r2_fd = $class->Get_FD(input => $pair_listing[1]);

        $ka_input = qq" ${r1_fd} ${r2_fd} ";
        $input_name = $pair_listing[0];
    } else {
        my $r1_fd = $class->Get_FD(input => $ka_input);
        $ka_input = qq" ${r1_fd} ";
        $ka_args .= qq' --bias --single -l 40 -s 10 ';
    }

    ## Check that the indexes exist
    my $ka_reflib = $paths->{index_file};
    my $index_job;
    if (!-r $ka_reflib) {
        my $transcriptome_fasta = $paths->{fasta};
        $index_job = $class->Bio::Adventure::Index::Kallisto_Index(
            input => $transcriptome_fasta,
            jdepends => $options->{jdepends},
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $outdir = $paths->{output_dir};
    my $stderr = $paths->{stderr};
    my $stdout = $paths->{stdout};
    my $output_sam = qq"${outdir}/kallisto_${species}.sam";
    my $output_bam = qq"${outdir}/kallisto_${species}.bam";
    my $output_stats = qq"${outdir}/kallisto_${species}.stats";
    my $sorted_bam = qq"${outdir}/kallisto_${species}-sorted";
    my $comment = qq!## This is a kallisto pseudoalignment of ${ka_input} against
## ${ka_reflib}.
## Other candidates for making a pretty count table include:
##  perl -F'\\t' -a -n -e 'print "\$F[0] \$F[3]\\n"' ${outdir}/abundance.tsv > ${outdir}/abundance.count
##   awk '{printf("%s %s\\n", \$1, \$4)}' ${outdir}/abundance.tsv > ${outdir}/abundance.count
## The sam->bam conversion is copy/pasted from Sam2Bam() but I figured why start another job
## because kallisto is so fast
!;
    my $dropped_args = ' --pseudobam ';
    my $jstring = qq!mkdir -p ${outdir}
kallisto quant ${ka_args} \\
  --plaintext -t 4 -b 100 \\
  -o ${outdir} \\
  -i ${ka_reflib} \\
  ${ka_input} \\
  2>${stderr} \\
  1>${output_sam} && \\
  cut -d "	" -f 1,4 ${outdir}/abundance.tsv > ${outdir}/${input_name}_abundance.count && \\
  xz -9e -f ${outdir}/${input_name}_abundance.count
  echo "Finished kallisto." >${stdout}
!;
    ## Newer kallisto does not seem to do this well anymore...
##  samtools view -u -t ${ka_reflib} -S ${output_sam} \\
##    2>${output_bam}.err \\
##    1>${output_bam} && \\
##  samtools sort -l 9 ${output_bam} ${sorted_bam} \\
##    2>${sorted_bam}.err \\
##    1>${sorted_bam}.out
##mv ${sorted_bam}.bam ${output_bam} && samtools index ${output_bam} && \\
##  bamtools stats -in ${output_bam} 2>${output_stats} 1>&2
##!;
    my $kallisto = $class->Submit(
        comment => $comment,
        input => $ka_input,
        depends => $options->{jdepends},
        jname => qq"${jname}",
        jprefix => '30',
        jstring => $jstring,
        jmem => $options->{jmem},
        output => qq"${outdir}/abundance.tsv",
        count => qq"${outdir}/${input_name}_abundance.count",
        stderr => $stderr,
        stdout => $stdout,
        prescript => $args{prescript},
        postscript => $args{postscript},);
    $kallisto->{index} = $index_job;
    return($kallisto);
}

=head2 C<PolyA_Recorder>

  This is a sister function to Bio::Adventure::Splicing::SLSeq_Recorder()
  One might reasonably ask why this is here and that there, I decided to separate
  them because the SLSeq is really only interesting in the context of trans-splicing
  in that relatively small cohort of organisms weird enough to use it.  In contrast,
  mapping the locations of polyA reads feels to me like it might have use in other
  contexts.  It should also be noted that this will likely only be useful following
  Bio::Adventure::Trim::PolyA_Extractor() and Hisat/BWA/etc.  With that in mind, it is
  likely that I will put these pieces together into a trans-splicing pipeline.

  The general idea: Read in the extant gene annotations and make a contig-keyed hash
  where each key is comprised of an array of gene features including the relevant
  information about the current gene and the 'next' gene (reading from beginning to end
  of each contig).  Then read each aligned read from a position-sorted bam file from hisat
  or whatever and keep only the reads which are: a) pointing in the same direction as the ORF
  b) positioned in the 3' UTR of an ORF (thus the information about the next gene).  Given that,
  record every polyA position with respect to the stop codon and count how many reads are
  observed at every relative position.  Upon completion, write up a new gff file with new
  features of type 'polyA' that have a new end for the + and new start for the - strand genes.

  Given the context I am using this, it should be unsurprising that it makes trypanosome-specific
  assumptions.

=over

  input(required): Sorted/indexed bam from hisat/bowtie/bwa/etc
  species: Find the original genome/gff annotations with this.
  gff_type(protein_coding_gene): Use this feature type to create new features.
  output_type(polyA): Create new features of this type.
  id_tag(ID): Identify genes using this feature tag.

=back

=cut
sub PolyA_Recorder {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        input_gff => '',
        gff_type => 'protein_coding_gene',
        gff_tag => 'ID',
        output_type => 'polyA_gene',
        jprefix => '50',);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $output_dir = $paths->{output_dir};
    make_path($output_dir) unless (-d $output_dir);
    my $output_gff = qq"${output_dir}/$options->{species}_maximum_3p_utr.gff";
    my $output_tsv = qq"${output_dir}/$options->{species}_recorded_polyA.tsv";
    my $jname = qq"polyA_$options->{species}";
    my $jstring = qq!use Bio::Adventure::Map;
my \$result = \$h->Bio::Adventure::Map::PolyA_Recorder_Worker(
  gff_type => '$options->{gff_type}',
  gff_tag => '$options->{gff_tag}',
  input => '$options->{input}',
  jname => '${jname}',
  output => '${output_gff}',
  output_tsv => '${output_tsv}',
  output_dir => '${output_dir}',
  output_type => '$options->{output_type}',
  species => '$options->{species}',);
!;
    my $comment = qq!## Seek polyA positions.!;
    my $polyA = $class->Submit(
        comment => $comment,
        gff_type => $options->{gff_type},
        gff_tag => $options->{gff_tag},
        input => $options->{input},
        input_gff => $options->{input_gff},
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $output_gff,
        output_dir => $output_dir,
        output_tsv => $output_tsv,
        output_type => $options->{output_type},
        language => 'perl',);
    return($polyA);
}

=head2 C<PolyA_Recorder_Worker>

  Does the actual work for PolyA_Recorder()

  This currently assumes a fairly specific read orientation and uses a couple of regexes
  to seek out and extract any polyA reads.  It trims off the As and optionally
  reverse-complements the sequence so that (hopefully) any valid polyA-derived sequence
  is in the same orientation as its upstream CDS.

=cut
sub PolyA_Recorder_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        input_gff => '',
        gff_type => 'protein_coding_gene',
        output_type => 'SLSeq_3pUTR',
        gff_tag => 'ID',);
    ## A couple of notes:
    ## 1.  These are all extracted from larger polyA containing reads
    ##     Thus I should not need to think about R1 vs R2 nor the read strand, only the
    ##     relative position of the nearest stop codon.

    my $input_gff = $options->{input_gff};
    unless ($input_gff) {
        $input_gff = qq"$options->{libpath}/$options->{libtype}/gff/$options->{species}.gff";
        print "Input gff file not provided, using:
${input_gff}.\n";
    }
    my $fasta = qq"$options->{libpath}/$options->{libtype}/fasta/$options->{species}.fasta";
    unless (-r $options->{input}) {
        die("Unable to find input bam alignment.");
    }
    my $output_dir = $options->{output_dir};
    make_path($output_dir) unless (-d $output_dir);
    ## Set up the output filehandles:
    ##  The output csv file.
    my $output_tsv_filename = qq"${output_dir}/$options->{species}_recorded_polya.tsv";
    my $output_tsv = FileHandle->new(qq">${output_dir}/$options->{species}_recorded_sl.tsv");
    print $output_tsv qq"$options->{gff_tag}\tstart\tend\tstrand\trelative_pos\tnum_reads\n";
    ## The output gff file, this will receive individual gff strings:
    my $gff_output_filename = qq"${output_dir}/$options->{species}_3p_modified_genes.gff";
    my $gff_out = FileHandle->new(">${gff_output_filename}");
    ## Use Bio::Tools::GFF to create the gff strings.
    my $gffout_io = Bio::Tools::GFF->new(-noparse => 1, -gff_version => 3);
    ## The log file:
    my $log_file = qq"${output_dir}/polya_recorder.log";
    my $log = FileHandle->new(">${log_file}");
    print $log "Starting PolyA Recorder, reading $options->{gff_tag} tags from $options->{gff_type}
entries of ${input_gff}.
Writing results to ${gff_output_filename}.\n";
    my ($counters, $observed) = $class->Bio::Adventure::Parsers::Count_Extract_GFF(
        gff => $input_gff, log => $log, gff_type => $options->{gff_type},
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
    $counters->{inside_gene} = 0;
  BAMLOOP: while (my $align = $bam->read1) {
        $counters->{reads}++;
        ##last BAMLOOP if ($counters->{quantified_reads} > 20000);
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
        my $read_unmappedp = $align->unmapped;
        if ($read_unmappedp) {
            $counters->{unmapped_reads}++;
            next BAMLOOP;
        }
        my @tags = $align->get_all_tags;
        my $not_primaryp = $align->get_tag_values('NOT_PRIMARY');
        unless ($not_primaryp) {
            $counters->{primary_reads}++;
        }
        my $supplementalp = $align->get_tag_values('SUPPLEMENTARY');
        ## I am going to take another look at the supplemental alignments
        ## My thought is that since these are only single ended and have lost the polyA region,
        ## perhaps they will align relatively poorly.
        if ($supplementalp) {
            $counters->{supplemental_reads}++;
        }
        if ($not_primaryp) {
            $counters->{secondary_reads}++;
        }
        if ($read_seqid ne $new_contig) {
            print "Starting new contig: ${read_seqid}\n";
            $previous_contig = $new_contig;
            if ($previous_contig ne 'start') {
                print "Writing data for ${previous_contig}\n";
                Write_PolyA_Data(observed => $observed,
                                 contig => $previous_contig,
                                 output_tsv => $output_tsv,
                                 gff_tag => $options->{gff_tag},
                                 gffio => $gffout_io,
                                 gff_out => $gff_out,);
            }
            $new_contig = $read_seqid;
        }
        ## Since these reads came from all manner of library types,
        ## I no longer want to consider the relative strand of the closest gene and the
        ## strand of the read.
        my $num_features = undef;

        ## Recall, that the read_seqid may include scaffolds which do not have annotated ORFs
        next BAMLOOP unless (defined($observed->{$read_seqid}));
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
        ## Now iterate over the features on this chromosome/contig and look for the following:
        ## 1.  If it is ---------->, then I want reads which are > the end.
        ## 2.  If it is <----------, then I want reads which are < the end.
        my $feature_start = 0;
      FEAT_SEARCH: for my $c (0 .. $num_features) {
            my $checker++;
            next FEAT_SEARCH if ($feature_start > $c);
            my $current_feat = $observed->{$read_seqid}[$c];
            my $current_tag = $current_feat->{primary_tag};
            my $current_start = $current_feat->{start};
            my $current_end = $current_feat->{end};
            my $current_strand = $current_feat->{strand};
            my $next_distance = $current_feat->{next_distance};
            my $next_start = $current_feat->{next_start};
            my $next_strand = $current_feat->{next_strand};
            my $previous_distance = $current_feat->{previous_distance};
            my $previous_end = $current_feat->{previous_end};
            my $previous_strand = $current_feat->{previous_strand};
            my $current_name = $current_feat->{gene_name};
            last FEAT_SEARCH if ($checker > $num_features);
            if (!defined($current_strand)) {
                print $log "In feature search, somehow strand is undefined for: ${current_name}\n";
                next FEAT_SEARCH;
            }
            ## Handle the (rare) case where reads are after a + strand gene and
            ## before a - strand gene
            if ($current_strand eq '1' && $current_end < $read_start &&
                $next_strand eq '-1' && $next_start > $read_start) {
                print $log "Contig:${read_seqid} Str:${current_strand} St:${read_start} vs. last end:${current_end} next start:${next_start}, unannotated ORF?\n";
                $counters->{unannotated_read}++;
                next BAMLOOP;
            }
            if ($current_strand ne $read_strand && $next_strand ne $read_strand) {
                $counters->{opposite_strand}++;
            }

            ## Don't bother with reads inside the current gene.
            if ($current_start >= $read_start &&
                $current_end <= $read_start) {
                $counters->{inside_gene}++;
                next BAMLOOP;
            }

            if ($current_strand eq '1') {
                my $relative_position = $read_end - $current_end;
                ## When we get to > 10000 nt past the first gene in the list
                ## stop looking at the first gene and move forward
                if ($relative_position < -10000) {
                    $feature_start++;
                }
                if ($read_end < $next_start) {
                    $counters->{quantified_reads}++;
                    ## print "Got a +1 hit: $counters->{quantified_reads}\n";
                    if (defined($observed->{$read_seqid}[$c]->{observed}->{$relative_position})) {
                        $observed->{$read_seqid}[$c]->{observed}->{$relative_position}++;
                    } else {
                        $observed->{$read_seqid}[$c]->{observed}->{$relative_position} = 1;
                    }
                    if ($relative_position > $observed->{$read_seqid}[$c]->{most_downstream}) {
                        $observed->{$read_seqid}[$c]->{most_downstream} = $relative_position;
                    }
                    last FEAT_SEARCH;
                } else {
                    next FEAT_SEARCH;
                }
            } elsif ($current_strand eq '-1') {
                my $relative_position = $current_start - $read_start;
                if ($read_end < $current_start) {
                    $counters->{quantified_reads}++;
                    if (defined($observed->{$read_seqid}[$c]->{observed}->{$relative_position})) {
                        $observed->{$read_seqid}[$c]->{observed}->{$relative_position}++;
                    } else {
                        $observed->{$read_seqid}[$c]->{observed}->{$relative_position} = 1;
                    }
                    if ($relative_position > $observed->{$read_seqid}[$c]->{most_downstream}) {
                        $observed->{$read_seqid}[$c]->{most_downstream} = $relative_position;
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
    Write_PolyA_Data(observed => $observed,
                     contig => $new_contig,
                     output_tsv => $output_tsv,
                     gff_tag => $options->{gff_tag},
                     gffio => $gffout_io,
                     gff_out => $gff_out,);
    print $log "Reads observed: $counters->{reads}
Unstranded reads: $counters->{unstranded_reads}
Plus reads: $counters->{plus_reads}
Minus reads: $counters->{minus_reads}
Unstranded: $counters->{notplusminus_reads}
Unmapped: $counters->{unmapped_reads}
Unmapped pair: $counters->{unmapped_pairs}
Primary: $counters->{primary_reads}
Supplemental: $counters->{supplemental_reads}
Secondary: $counters->{secondary_reads}
Same strand inside gene: $counters->{same_strand_in_gene}
Unannotated: $counters->{unannotated_read}
Opposite strand: $counters->{opposite_strand}
Quantified: $counters->{quantified_reads}
";
    $gff_out->close();
    $log->close();
    return($observed);
}

=head2 C<Write_PolyA_Data>

 Sister to the Write_SL_Data function.

 Use this to record information about each observed polyA-derived read.

=cut
sub Write_PolyA_Data {
    my %args = @_;
    my $observed = $args{observed};
    my $contig = $args{contig};
    my $output_tsv = $args{output_tsv};
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
                print $output_tsv qq"${gene_name}\t${gene_start}\t${gene_end}\t${gene_strand}\t${obs}\t$observations->{$obs}\n";
            }
            if ($gene_strand eq '1') {
                $new_end = $gene_end + $farthest;
                $new_start = $gene_end;
            } elsif ($gene_strand eq '-1') {
                $new_start = $gene_start - $farthest;
                $new_end = $gene_start;
            } else {
                die("This should not happen, neither +1/-1.\n");
            }
            if ($new_end < 0) {
                print "This is a failed entry: $contig $gene_name with end: $new_end\n";
                $new_end = 0;
            }
            my $standardized = Bio::SeqFeature::Generic->new(
                -primary_tag => 'polyA',
                -display_name => $gene_name,
                -seq_id => $contig,
                -start => $new_start,
                -end => $new_end,
                -strand => $gene_strand,
                -score => 0,
                -tag => {
                    $gff_tag => $gene_name,
                    inference => 'polyA',
                    old_start => $gene_start,
                    old_end => $gene_end,
                });
            my $gff_string = $gffio->gff_string($standardized);
            print $gff_out "${gff_string}\n";
        }
    }
}

=head2 C<RSEM>

 Invoke RSEM.
 10.1186/1471-2105-12-323

 The most accurate and slow transcript quantification method.

=over

 input(required): Trimmed fastq input file(s).
 species(required): Used to find the indexes

=back

Currently everything else is left as the default.

=cut
sub RSEM {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        jmem => 24,);
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking RSEM on ${sp}\n";
            $options->{species} = $sp;
            my $result = Bio::Adventure::Map::RSEM($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $jbasename = $class->Get_Job_Name();
    my $cds = qq"$options->{libpath}/$options->{libtype}/$options->{species}_cds_nt.fasta";
    my $idx = qq"$options->{libpath}/$options->{libtype}/indexes/rsem/$options->{species}";
    my $test_idx = qq"${idx}.transcripts.fa";
    my $index_job;
    unless (-r $test_idx) {
        print "Need to create the RSEM indexes, invoking RSEM_Index().\n";
        unless (-r $cds) {
            die("RSEM_Index requires a cds fasta file at: ${cds}, create it with a cyoa2 conversion.");
        }
        $index_job = $class->Bio::Adventure::Index::RSEM_Index(
            input => $cds,
            index => $idx,
            jdepends => $options->{jdepends},
            jname => qq"rsidx_${jbasename}",
            libtype => $options->{libtype},
            modules => $options->{modules},);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $rsem_input = $options->{input};
    my $test_file = "";
    if ($rsem_input =~ /$options->{delimiter}/) {
        my @pair_listing = split(/$options->{delimiter}/, $rsem_input);
        my $r1_fd = $class->Get_FD(input => $pair_listing[0]);
        my $r2_fd = $class->Get_FD(input => $pair_listing[1]);
        $rsem_input = qq"--paired-end ${r1_fd} ${r2_fd}";
        $test_file = $pair_listing[0];
    } else {
        my $r1_fd = $class->Get_FD(input => $rsem_input);
        $rsem_input = qq" ${r1_fd} ";
    }

    my $rsem_dir = qq"outputs/rsem_$options->{species}";
    my $stderr = qq"${rsem_dir}/$options->{species}.stderr";
    my $stdout = qq"${rsem_dir}/$options->{species}.stdout";
    my $rsem_comment = '## This is a rsem invocation script.';
    my $output_file = qq"$options->{species}_something.txt";
    my $jstring = qq"mkdir -p ${rsem_dir} && rsem-calculate-expression --bowtie2 \\
  --calc-ci ${rsem_input} \\
  ${idx} \\
  ${rsem_dir}/$options->{species} \\
  2>${stderr} \\
  1>${stdout}
";
    my $rsem = $class->Submit(
        comment => $rsem_comment,
        input => $rsem_input,
        jname => qq"rsem_${jbasename}",
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jprefix => '28',
        mem => $options->{jmem},
        output => $output_file,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $rsem->{index_job} = $index_job;
    return($rsem);
}

=head2 C<Salmon>

 Perform a salmon quantification of transcript abundances.
 10.1038/nmeth.4197

 My favorite transcript aware quantification method.

=over

 input(required): Trimmed file(s) containing the reads.
 species(required): Used to locate the indexes.

=back

=cut
sub Salmon {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        jcpu => 8,
        libtype => 'CDS',
        jwalltime => '4:00:00',
        jmem => 24,
        stranded => 'A',
        jprefix => '45',);
    if ($options->{libtype} ne 'CDS') {
        print "The library type for salmon generally needs to be CDS unless this is something special.
Waiting 10 seconds to see if you change your mind.\n";
        sleep(10);
    }
    my $depends = $options->{jdepends};
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking salmon on ${sp}\n";
            $options->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::Salmon(%{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $ready = $class->Check_Input(files => $options->{input});
    my %sa_jobs = ();
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = $options->{species};

    my $jname = qq"sal_${species}_$options->{libtype}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $sa_args = '';
    my $sa_input = $options->{input};
    my $input_name = $sa_input;
    my $num_inputs = 1;
    if ($sa_input =~ /$options->{delimiter}/) {
        my @pair_listing = split(/$options->{delimiter}/, $sa_input);
        my $r1_fd = $class->Get_FD(input => $pair_listing[0]);
        my $r2_fd = $class->Get_FD(input => $pair_listing[1]);
        $sa_args .= qq" -1 ${r1_fd} -2 ${r2_fd} ";
        $input_name = $pair_listing[0];
        $num_inputs = 2;
    } else {
        my $r1_fd = $class->Get_FD(input => $sa_input);
        $sa_args .= qq" -r ${r1_fd} ";
    }

    ## A reminder from the salmon documentation regarding the various strand declarations:
    ## 'A': autodetect
    ##   Orientations:
    ## I = inward
    ## O = outward
    ## M = matching
    ##   Stranded vs. not: S/U
    ##   read1 is on the forward/reverse strand: F/R
    ## Thus, a paired end, stranded, away from each other library <-> with r1 on the reverse strand:
    ## 'OSR'
    my $strand_string = $options->{stranded};
    if (!defined($strand_string)) {
        $strand_string = 'A';
    } elsif ($strand_string eq '0' || 'no') {
        if ($num_inputs == 1) {
            $strand_string = 'U';
        } else {
            $strand_string = 'IU';
        }
    } elsif ($strand_string eq '1' || 'forward') {
        $strand_string = 'ISF';
    } elsif ($strand_string eq '-1' || 'reverse') {
        $strand_string = 'ISR';
    } elsif ($strand_string =~ /S|U|F|R|I|O|M/) {
        print "The strand string: ${strand_string} was provided to cyoa.\n";
    } else {
        $strand_string = 'A';
    }

    ## Check that the indexes exist
    my $sa_reflib = $paths->{index_dir};
    my $index_job;
    my $cluster = $options->{cluster};
    if (!-d $sa_reflib) {
        my $transcript_file = qq"$options->{libpath}/${libtype}/fasta/$options->{species}.fasta";
        my $index_extras = {
            input => $transcript_file,
        };
        my %index_args = $class->Extra_Options(options => $options, extras => $index_extras);
        $index_job = $class->Bio::Adventure::Index::Salmon_Index(%index_args, cluster => 0);
    }
    my $outdir = $paths->{output_dir};
    my $stderr = $paths->{stderr};
    my $stdout = $paths->{stdout};
    my $comment = qq!## This is a salmon pseudoalignment of ${sa_input} against
## ${sa_reflib}.
!;
    my $jstring = qq!mkdir -p ${outdir}
mapped=1
{
  /usr/bin/time -v -o ${stdout}.time -a \\
    salmon quant -i ${sa_reflib} \\
      -l ${strand_string} --gcBias --validateMappings  \\
      ${sa_args} \\
      -o ${outdir} \\
      2>${stderr} \\
      1>${stdout}
} || {
  echo 'salmon failed.' >> ${stderr}
}
!;
    ## Explicitly setting the cluster state to whatever it was before
    ## testing for the indexes, because I forced the index creation
    ## to run outside of the context of a cluster.
    my $salmon = $class->Submit(
        cluster => $cluster,
        comment => $comment,
        input => $sa_input,
        jdepends => $options->{jdepends},
        jname => $jname,
        jstring => $jstring,
        output => $paths->{output},
        stderr => $stderr,
        stdout => $stdout,
    );
    $salmon->{index_job} = $index_job;
    my $stats = $class->Bio::Adventure::Metadata::Salmon_Stats(
        input => qq"${outdir}/lib_format_counts.json",
        jdepends => $salmon->{job_id},
        jprefix => qq"$options->{jprefix}_1",
        jname => qq"sastats_$options->{species}",
    );
    $salmon->{stats} = $stats;
    $salmon->{job_id} = $stats->{job_id};
    return($salmon);
}

=head2 C<STAR>

 Invoke STAR for transcript abundances.
 10.1093/bioinformatics/bts635

=cut
sub STAR {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        jmem => 48,);
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking STAR on ${sp}\n";
            $options->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::STAR(%{$options});
            push (@result_lst, $result);
        }
        return(@result_lst);
    }

    my $ready = $class->Check_Input(files => $options->{input});
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = $options->{species};
    my $jname = qq"star_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $star_inputstring = qq"";
    my $star_input = $options->{input};
    my $input_name = $star_input;
    if ($star_input =~ /$options->{delimiter}/) {
        my @pair_listing = split(/$options->{delimiter}/, $star_input);
        $star_inputstring = qq"$pair_listing[0],$pair_listing[1]";
        $input_name = $pair_listing[0];
    } else {
        $star_inputstring = qq"${star_input}";
    }

    ## Check that the indexes exist
    my $star_refdir = qq"$options->{libpath}/${libtype}/indexes/$options->{species}_star_index";
    my $star_reflib = qq"${star_refdir}/SAindex";
    my $index_job;
    if (!-r $star_reflib) {
        my $genome_file = qq"$options->{libpath}/${libtype}/$options->{species}.fasta";
        $index_job = $class->Bio::Adventure::Index::STAR_Index(
            input => $genome_file,
            jdepends => $options->{jdepends},
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }
    my $outdir = qq"outputs/star_${species}";
    my $stderr = qq"${outdir}/star_${species}.stderr";
    my $stdout = qq"${outdir}/star_${species}.stdout";
    my $comment = qq!## This is a star pseudoalignment of ${star_input} against
## ${star_reflib}.
## This jobs depended on: $options->{jdepends}
## Currently, this only works with the module star/git_201803
!;
    my $jstring = qq!mkdir -p ${outdir}
STAR \\
  --genomeDir ${star_refdir} \\
  --outFileNamePrefix outputs/star_$options->{species}/${input_name} \\
  --outSAMtype BAM SortedByCoordinate \\
  --outBAMcompression 10 \\
  --chimOutType WithinBAM \\
  --quantMode GeneCounts \\
  --readFilesIn ${star_inputstring} \\
  --readFilesCommand /usr/bin/lesspipe.sh \\
  --runThreadN 6 \\
  2>${stderr} \\
  1>${stdout}

STAR \\
  --genomeDir ${star_refdir} \\
  --outFileNamePrefix outputs/star_$options->{species}/${input_name}_tx \\
  --outSAMtype BAM SortedByCoordinate \\
  --outBAMcompression 10 \\
  --chimOutType WithinBAM \\
  --quantMode TranscriptomeSAM \\
  --readFilesIn ${star_inputstring} \\
  --readFilesCommand /usr/bin/lesspipe.sh \\
  --runThreadN 6 \\
  2>>${stderr} \\
  1>>${stdout}
!;

    my $star_job = $class->Submit(
        comment => $comment,
        input => $star_input,
        jdepends => $options->{jdepends},
        jname => qq"${jname}",
        jprefix => '33',
        jstring => $jstring,
        jmem => $options->{jmem},
        prescript => $args{prescript},
        postscript => $args{postscript},
        stderr => $stderr,
        stdout => $stdout,);
    $star_job->{index_job} = $index_job;
    return($star_job);
}

=head2 C<Tophat>

 Invokes tophat!
 10.1186/gb-2013-14-4-r36

 It also sorts/indexes the accepted_hits file, collects some
 statistics, and passes the hits to htseq-count.

=cut
sub Tophat {
    my ($class, %args) = @_;
    my $check = which('tophat');
    die('Could not find tophat in your PATH.') unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input', 'gff_type'],);
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking tophat on ${sp}\n";
            $options->{species} = $sp;
            my $result = Bio::Adventure::Map::Tophat($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $tophat_cpus = 4;
    my $inputs = $options->{input};
    my @in = split(/:/, $inputs);
    $inputs =~ s/:/ /g;
    my $paired = 0;
    if (scalar(@in) > 1) {
        my $r1_fd = $class->Get_FD(input => $in[0]);
        my $r2_fd = $class->Get_FD(input => $in[1]);
        $inputs = qq" ${r1_fd} ${r2_fd} ";
        $paired = 1;
    } else {
        my $r1_fd = $class->Get_FD(input => $in[0]);
        $inputs = qq" ${r1_fd} ";
    }
    my $tophat_args = ' -g 1 --microexon-search --b2-very-sensitive ';
    if ($options->{tophat_args}) {
        $tophat_args = $options->{tophat_args};
    }
    ## $tophat_args .= ' --no-mixed --no-discordant ' if (scalar(@in) > 1);
    ## $tophat_args .= ' ' if (scalar(@in) > 1);

    my $tophat_walltime = '18:00:00';
    my $tophat_mem = 8;
    if ($options->{species} eq 'hsapiens' or $options->{species} eq 'mmusculus') {
        $tophat_walltime =  '144:00:00';
        $tophat_mem = 20;
    }

    my $tophat_dir = qq"outputs/tophat_$options->{species}";
    if ($options->{tophat_dir}) {
        $tophat_dir = $options->{tophat_dir};
    }
    my $libtype = $options->{libtype};
    my $bt_reflib = qq"$options->{libpath}/${libtype}/indexes/$options->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    my $index_job = undef;
    if (!-r $bt_reftest) {
        print "Did not find the index for $options->{species} at: ${bt_reflib}, indexing now.\n";
        my $genome_file = qq"$options->{libpath}/${libtype}/$options->{species}.fasta";
        $index_job = $class->Bio::Adventure::Index::BT2_Index(
            input => $genome_file,
            jdepends => $options->{jdepends},);
        $options->{jdepends} = $index_job->{job_id};
    }
    my $gtf_file = qq"$options->{libpath}/genome/$options->{species}.gtf";
    if (!-r $gtf_file) {
        print "Missing the gtf file for $options->{species}\n";
        print "Using the gff file.\n";
        $gtf_file =~ s/\.gtf/\.gff/;
    }

    my $stderr = qq"${tophat_dir}/tophat.stderr";
    my $stdout = qq"${tophat_dir}/tophat.stdout";
    my $spliced = 0;
    if ($options->{spliced}) {
        $spliced = 1;
    }
    my $jname = qq"th_$options->{species}";
    my $jstring = qq!
mkdir -p ${tophat_dir}
tophat ${tophat_args} \\
  -G ${gtf_file} \\
  -p ${tophat_cpus} -o ${tophat_dir} \\
!;
    if ($spliced) {
        $jstring .= qq!  --no-novel-juncs \\
!;
    }
    $jstring .= qq!  $options->{libpath}/genome/indexes/$options->{species} \\
  ${inputs} \\
  2>${stderr} \\
  1>${stdout}
samtools sort -l 9 -n ${tophat_dir}/accepted_hits.bam -o ${tophat_dir}/accepted_sorted.bam \\
  2>${tophat_dir}/samtools_sort.stderr 1>${tophat_dir}/samtools_sort.stdout
samtools index ${tophat_dir}/accepted_hits.bam \\
  2>${tophat_dir}/samtools_index.stderr 1>${tophat_dir}/samtools_index.stdout
samtools sort -l 9 -n ${tophat_dir}/unmapped.bam -o ${tophat_dir}/unmapped_sorted.bam \\
  2>>${tophat_dir}/samtools_sorted.stderr 1>>${tophat_dir}/samtools_sort.stdout
samtools index ${tophat_dir}/unmapped.bam \\
  2>>${tophat_dir}/samtools_index.stderr 1>>${tophat_dir}/samtools_sort.stdout
!;
    if ($paired) {
        $jstring .= qq!
if [ -r "${tophat_dir}/accepted_hits.bam" ]; then
  samtools view -b -f 2 ${tophat_dir}/accepted_hits.bam > ${tophat_dir}/accepted_paired.bam
  samtools index ${tophat_dir}/accepted_paired.bam
fi
!;
    }
    my $comment = qq!## I still have no clue what I am doing when I use tophat...
## However, I know that -g 1 will allow only 1 hit in the case of multihits, but randomly place it
## From the manual:  "If there are more alignments with the same score than this
## number, TopHat will randomly report only this many alignments"
## -N 1 will discard anything with >1 mismatch (default is 2)
## -r adjusts the allowable mean distance between the paired reads
## --mate-std-dev sets the deviation of -r
## --microexon-search will tell it to search short exons for reads >=50
!;
    my $tophat = $class->Submit(
        comment => $comment,
        jcpu => $tophat_cpus,
        jdepends => $options->{jdepends},
        jname => $options->{jname},
        jprefix => '31',
        jstring => $jstring,
        jmem => $tophat_mem,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jwalltime => $tophat_walltime,);
    ## Set the input for htseq
    my $accepted = qq"${tophat_dir}/accepted_hits.bam";
    $accepted = $options->{accepted_hits} if ($options->{accepted_hits});
    my $count_table = 'accepted_hits.count';
    $count_table = $options->{count_table} if ($options->{count_table});
    my $htmulti = $class->Bio::Adventure::Count::HT_Multi(
        htseq_input => $accepted,
        jdepends => $tophat->{job_id},
        jname => qq"hts_$options->{species}",
        jprefix => '32',
        mapper => 'tophat',);
    $tophat->{htseq} = $htmulti;
    ## Perform a separate htseq run using only the successfully paired hits
    if ($paired) {
        my $ht_paired = $class->Bio::Adventure::Count::HT_Multi(
            htseq_input => qq"${tophat_dir}/accepted_paired.bam",
            jdepends => $tophat->{job_id},
            jname => qq"htsp_$options->{species}",
            jprefix => '32',
            mapper => 'tophat',);
    }
    ## Tophat_Stats also reads the trimomatic output, which perhaps it should not.
    my $unaccepted = $accepted;
    $unaccepted =~ s/accepted_hits/unmapped/g;
    my $input_read_info = $accepted;
    $input_read_info =~ s/accepted_hits\.bam/prep_reads\.info/g;
    my $stats = $class->Bio::Adventure::Metadata::Tophat_Stats(
        accepted_input => $accepted,
        count_table => qq"${count_table}.xz",
        jdepends => $tophat->{job_id},
        jname => qq"tpstats_$options->{species}",
        jprefix => '33',
        prep_input => $input_read_info,
        unaccepted_input => $unaccepted,);
    $tophat->{stats} = $stats;
    $tophat->{job_id} = $stats->{job_id};
    return($tophat);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;
