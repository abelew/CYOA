package Bio::Adventure::RNASeq_Map;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::Which qw"which";

=head1 NAME

    Bio::Adventure::RNASeq_Map - Perform highthroughput sequence alignments with tools like bowtie/tophat/etc

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->Bowtie();

=head2 Methods

=over 4

=item C<Bowtie>

    $hpgl->Bowtie() performs a bowtie alignment.  Unless instructed
    otherwise, it will do so with 0 mismatches and with multi-matches
    randomly placed 1 time among the possibilities. (options -v 0 -M 1)


    It checks to see if a bowtie1 compatible index is in
    $libdir/$libtype/indexes/$species, if not it attempts to create
    them.

    It will continue on to convert the bowtie sam output to a
    compressed, sorted, indexed bam file, and pass that to htseq-count
    using a gff file of the same species.

=cut
sub Bowtie {
    my ($class, %args) = @_;
    my $check = which('bowtie-build');
    die("Could not find bowtie in your PATH.") unless($check);

    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => 'lmajor', 'bt_type' => 'v0M1', count => 1,
        libtype => 'genome',
    );
    my $species = $options->{species};
    my $bt_type = $options->{bt_type};
    my $bt_args = $options->{bt_args}->{$bt_type};
    $bt_args = ' --best -v 0 -M 1 ' if (!defined($bt_args));

    my $sleep_time = 3;
    my %bt_jobs = ();
    my $bt_input = $options->{input};
    my $bt_depends_on;
    $bt_depends_on = $options->{job_depends} if ($options->{job_depends});
    my $job_basename = $options->{job_basename};

    my $job_name = qq"bt${bt_type}_${species}";
    $job_name = $options->{job_name} if ($options->{job_name});
    my $libtype = $options->{libtype};
    my $count = $options->{count};

    my $bt_dir = qq"outputs/bowtie_${species}";
    $bt_dir = $options->{bt_dir} if ($options->{bt_dir});

    my $uncompress_jobid = undef;
    my $index_jobid = undef;
    if ($bt_input =~ /\.gz$|\.bz2$|\.xz$/ ) {
        print "The input needs to be uncompressed, doing that now.\n" if ($options->{debug});
        my $uncomp = Bio::Adventure::Compress::Uncompress(
            $class,
            input => $bt_input,
            job_depends => $bt_depends_on,
            job_name => 'uncomp',
            job_prefix => '09',
        );
        $bt_input = basename($bt_input, ('.gz', '.bz2', '.xz'));
        $bt_jobs{uncompress} = $uncomp;
        $options = $class->Set_Vars(input => $bt_input);
        $uncompress_jobid = $uncomp->{job_id};
    }

    ## Check that the indexes exist
    my $bt_reflib = "$options->{libdir}/${libtype}/indexes/${species}";
    my $bt_reftest = qq"${bt_reflib}.1.ebwt";
    if (!-r $bt_reftest && !$options->{bt1_indexjobs}) {
        $options = $class->Set_Vars(bt1_indexjobs => 1);
        my $index_job = Bio::Adventure::RNASeq_Map::BT1_Index(
            $class,
            job_depends => $bt_depends_on,
            libtype => $libtype,
        );
        $bt_jobs{index} = $index_job;
        $index_jobid = $index_job->{job_id};
    }

    ## Make a depends string containing the uncompress, indexer, or both
    if (defined($index_jobid) && defined($uncompress_jobid)) {
        $bt_depends_on = qq"${index_jobid}:${uncompress_jobid}";
    } elsif (defined($index_jobid)) {
        $bt_depends_on = $index_jobid;
    } elsif (defined($uncompress_jobid)) {
        $bt_depends_on = $uncompress_jobid;
    } else {
        $bt_depends_on = "";
    }

    my $bowtie_input_flag = "-q"; ## fastq by default
    $bowtie_input_flag = "-f" if ($options->{input} =~ /\.fasta/);

    my $cpus = $options->{cpus};
    my $error_file = qq"${bt_dir}/${job_basename}-${bt_type}.err";
    my $comment = qq!## This is a bowtie1 alignment of ${bt_input} against
## ${bt_reflib} using arguments: ${bt_args}.
!;
    if ($bt_depends_on) {
        $comment .= qq!## This jobs depended on: ${bt_depends_on}.
!;
    }
    my $aligned_filename = qq"${bt_dir}/${job_basename}-${bt_type}_aligned_${species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/${job_basename}-${bt_type}_unaligned_${species}.fastq";
    my $sam_filename = qq"${bt_dir}/${job_basename}-${bt_type}.sam";
    my $job_string = qq!mkdir -p ${bt_dir} && sleep ${sleep_time} && bowtie ${bt_reflib} ${bt_args} \\
  -p ${cpus} ${bowtie_input_flag} ${bt_input} \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>${bt_dir}/${job_basename}-${bt_type}.out
!;

    my $bt_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        input => $bt_input,
        job_depends => $bt_depends_on,
        job_name => $job_name,
        job_output => $sam_filename,
        job_prefix => "10",
        job_string => $job_string,
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        queue => 'workstation',
        unaligned => $unaligned_filename,
    );
    $bt_jobs{bowtie} = $bt_job;

    my $un_comp = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt_args}.\n",
        job_depends => $bt_job->{job_id},
        job_name => "xzun",
        job_prefix => "11",
        xz_input => "${bt_dir}/${job_basename}-${bt_type}_unaligned_${species}.fastq",
    );
    $bt_jobs{unaligned_compression} = $un_comp;

    my $al_comp = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the sequences which successfully aligned against ${bt_reflib} using options ${bt_args}.",
        job_depends => $bt_job->{job_id},
        job_name => "xzal",
        job_prefix => "11",
        xz_input => "${bt_dir}/${job_basename}-${bt_type}_aligned_${species}.fastq",
    );
    $bt_jobs{aligned_compression} = $al_comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    my $trim_output_file = qq"outputs/trimomatic_stats.csv";

    my $sam_job = Bio::Adventure::Convert::Samtools(
        $class,
        input => $sam_filename,
        job_depends => $bt_job->{job_id},
        job_name => "${job_name}_s2b",
        job_prefix => "13",
    );
    $bt_jobs{samtools} = $sam_job;
    $options = $class->Set_Vars(job_output => $sam_job->{job_output});
    my $htmulti;
    if ($count) {
        if ($libtype eq 'rRNA') {
            $htmulti = Bio::Adventure::RNASeq_Count::HTSeq(
                $class,
                htseq_id => $options->{htseq_id},
                htseq_input => $sam_job->{job_output},
                htseq_type => $options->{htseq_type},
                job_depends => $sam_job->{job_id},
                job_name => "ht_${job_name}",
                job_prefix => '14',
                libtype => $libtype,
                queue => 'workstation',
                suffix => $bt_type,
            );
        } else {
            $htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
                $class,
                htseq_id => $options->{htseq_id},
                htseq_input => $sam_job->{job_output},
                htseq_type => $options->{htseq_type},
                job_depends => $sam_job->{job_id},
                job_name => "ht_${job_name}",
                job_prefix => '14',
                libtype => $libtype,
                queue => 'workstation',
                suffix => $bt_type,
            );
            $bt_jobs{htseq} = $htmulti;
        }
    }  ## End if ($count)

    my $stats = Bio::Adventure::RNASeq_Map::BT1_Stats(
        $class, %args,
        bt_input => $error_file,
        bt_type => $bt_type,
        count_table => qq"${job_basename}-${bt_type}.count.xz",
        job_depends => $bt_job->{job_id},
        job_name => "${job_name}_stats",
        job_prefix => "12",
        trim_input => ${trim_output_file},
    );
    $bt_jobs{stats} = $stats;

    return(\%bt_jobs);
}

=item C<Bowtie2>

    $hpgl->Bowtie2() performs a bowtie2 alignment.  Unless instructed
    otherwise, it will do so with 0 mismatches and with multi-matches
    randomly placed 1 time among the possibilities. (options -v 0 -M 1)

    It checks to see if a bowtie2 compatible index is in
    $libdir/$libtype/indexes/$species, if not it attempts to create
    them.

    It will continue on to convert the bowtie sam output to a
    compressed, sorted, indexed bam file, and pass that to htseq-count
    using a gff file of the same species.

=cut
sub Bowtie2 {
    my ($class, %args) = @_;
    my $check = which('bowtie2-build');
    die("Could not find bowtie2 in your PATH.") unless($check);

    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input', 'htseq_type'],
        do_htseq => 1,
    );
    my $sleep_time = 3;
    my %bt_jobs = ();
    my $libtype = 'genome';
    my $bt_depends_on = "";
    $bt_depends_on = $options->{job_depends} if ($options->{job_depends});
    my $bt2_args = $options->{bt2_args};
    my $job_name = qq"bt2_$options->{species}";
    $job_name = $options->{job_name} if ($options->{job_name});

    my $job_basename = $options->{job_basename};

    my $bt_dir = qq"outputs/bowtie2_$options->{species}";
    if ($args{bt_dir}) {
        $bt_dir = $args{bt_dir};
    }
    my $bt_input = $options->{input};
    if ($bt_input =~ /\.gz$|\.bz2$|\.xz$/ ) {
        print "The input needs to be uncompressed, doing that now.\n" if ($options->{debug});
        my $uncomp = Bio::Adventure::Compress::Uncompress(
            $class,
            input => $bt_input,
            job_depends => $bt_depends_on,
        );
        $bt_input =~ s/\:|\;|\,|\s+/ /g;
        $bt_input =~ s/\.gz|\.bz|\.xz//g;
        ##$bt_input = $bt_inputbasename($bt_input, ('.gz', '.bz2', '.xz'));
        $options = $class->Set_Vars(input => $bt_input);
        $bt_depends_on = $uncomp->{job_id};
    }

    my $test_file = "";
    if ($bt_input =~ /\s+/) {
        my @pair_listing = split(/\s+/, $bt_input);
        $bt_input = qq" -1 $pair_listing[0] -2 $pair_listing[1] ";
        $test_file = $pair_listing[0];
    } else {
        $test_file = $bt_input;
    }

    ## Check that the indexes exist
    my $bt_reflib = "$options->{libdir}/$options->{libtype}/indexes/$options->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    if (!-r $bt_reftest) {
        print "Hey! The Indexes do not appear to exist, check this out: ${bt_reftest}\n";
        sleep(20);
        my $index_job = Bio::Adventure::RNASeq_Map::BT2_Index(
            $class,
            job_depends => $bt_depends_on,
            libtype => $libtype,
        );
        $bt_jobs{index} = $index_job;
        $bt_depends_on = $index_job->{job_id};
    }
    my $bowtie_input_flag = "-q "; ## fastq by default
    $bowtie_input_flag = "-f " if (${bt_input} =~ /\.fasta$/);

    my $cpus = $options->{cpus};
    my $error_file = qq"${bt_dir}/${job_basename}.err";
    my $comment = qq!## This is a bowtie2 alignment of ${bt_input} against
## ${bt_reflib} using arguments: ${bt2_args}.
## This jobs depended on: ${bt_depends_on}
!;
    my $aligned_filename = qq"${bt_dir}/${job_basename}_aligned_$options->{species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/${job_basename}_unaligned_$options->{species}.fastq";
    my $sam_filename = qq"${bt_dir}/${job_basename}.sam";
    my $job_string = qq!mkdir -p ${bt_dir} && sleep ${sleep_time} && bowtie2 -x ${bt_reflib} ${bt2_args} \\
  -p ${cpus} \\
  ${bowtie_input_flag} ${bt_input} \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>${bt_dir}/${job_basename}.out
!;

    ## Last chance before submitting, make sure the input file exists.
    if (!-r "${test_file}") {
        die("Could not find the bowtie2 input file(s): ${bt_input}\n");
    }

    my $bt2_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        input => $bt_input,
        job_name => $job_name,
        job_depends => $bt_depends_on,
        job_string => $job_string,
        job_prefix => "15",
        output => $sam_filename,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        unaligned => $unaligned_filename,
    );
    $bt_jobs{bowtie} = $bt2_job;

    my $un_comp = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt2_args}\n",
        input => "${bt_dir}/${job_basename}_unaligned_$options->{species}.fastq",
        job_depends => $bt2_job->{job_id},
        job_name => "xzun",
        job_prefix => "16",
        );
    $bt_jobs{unaligned_compression} = $un_comp;

    my $al_comp = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the sequences which successfully aligned against ${bt_reflib} using options ${bt2_args}",
        input => "${bt_dir}/${job_basename}_aligned_$options->{species}.fastq",
        job_name => "xzal",
        job_prefix => "17",
        job_depends => $bt2_job->{job_id},
    );
    $bt_jobs{aligned_compression} = $al_comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    ## my $trim_output_file = qq"outputs/${job_basename}-trimomatic.out";
    my $stats = Bio::Adventure::RNASeq_Map::BT2_Stats(
        $class,
        bt_input => $error_file,
        count_table => qq"${job_basename}.count.xz",
        job_depends => $bt2_job->{job_id},
        job_name => "${job_name}_stats",
        job_prefix => "18",
        ## trim_input => ${trim_output_file},
    );
    my $sam_job = Bio::Adventure::Convert::Samtools(
        $class,
        input => $sam_filename,
        job_depends => $bt2_job->{job_id},
        job_name => "${job_name}_s2b",
        job_prefix => "19",
    );
    $bt_jobs{samtools} = $sam_job;
    my $htseq_input = $sam_job->{job_output};
    my $htmulti;
    if ($options->{do_htseq}) {
        if ($libtype eq 'rRNA') {
            $htmulti = Bio::Adventure::RNASeq_Count::HTSeq(
                $class,
                htseq_input => $sam_job->{job_output},
                job_depends => $sam_job->{job_id},
                job_prefix => "20",
                libtype => $libtype,
            );
        } else {
            $htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
                $class,
                htseq_input => $sam_job->{job_output},
                job_depends => $sam_job->{job_id},
                job_prefix => "21",
                libtype => $libtype,
            );
            $bt_jobs{htseq} = $htmulti;
        }
    }
    return(\%bt_jobs);
}

=item C<BT_Multi>

    $hpgl->BT_Multi() attempts to run multiple bowtie1 runs for a
    given species.  One run is performed for each of a few parameters
    which are kept in the variable '$hpgl->{bt_types}' and generally
    include: 0 mismatch, 1 mismatch, 2 mismatches, 1 randomly placed
    hit, 0 randomly placed hits, or default options.

=cut
sub BT_Multi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species","input", "htseq_type"],
    );
    my $job_basename = $options->{job_basename};
    my $species = $options->{species};
    my $depends_on = $options->{job_depends};
    my %bt_types = %{$options->{bt_args}};
    my @jobs = ();
    foreach my $type (keys %bt_types) {
        my $job_name = qq"bt${type}_${species}";
        my $job = Bio::Adventure::RNASeq_Map::Bowtie(
            $class,
            bt_type => $type,
            job_depends => $depends_on,
            job_name => $job_name,
            prescript => $args{prescript},
            postscript => $args{postscript},
        );
        push(@jobs, $job);
    }
    return(\@jobs);
}

=item C<Bowtie_RRNA>

    $hpgl->Bowtie_RRNA() performs an alignment against a home-curated
    set of ribosomal RNA/tRNA sequences.  The alignment requires a
    fastq input library and fasta library found in
    'libraries/rRNA/$class->{species}.fasta'

  Example:
    my $rrna = $hpgl->Bowtie_RRNA();
    ## If you want to exclude the rRNA sequences from future alignments:
    my $rrna = $hpgl->Bowtie_RRNA(exclude => 1);

=cut
sub Bowtie_RRNA {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species"],
    );
    my $job_basename = $options->{job_basename};
    $job_basename = qq"rRNA_${job_basename}";
    my $exclude = 0;
    $exclude = $options->{exclude} if ($options->{exclude});
    my $species = $options->{species};
    my $depends_on = $options->{depends};
    my ${bt_dir} = qq"outputs/bowtie_$options->{species}";
    if ($exclude) {
        my $in = $options->{input};
        $options = $class->Set_Vars(
            postscript => qq"mv $in includerrna_${in} && mv ${bt_dir}/${job_basename}-rRNA_unaligned_${species}.fastq $in",
        );
    }
    my $job = Bio::Adventure::RNASeq_Map::Bowtie(
        $class,
        job_depends => $depends_on,
        job_name => qq"btrrna",
        libtype => 'rRNA',
        prescript => $args{prescript},
        postscript => $args{postscript},
    );
    ## Return the basename back to normal so that future tasks don't
    ## get confuseled.
    return($job);
}

=item C<BT1_Index>

    $hpgl->BT1_Index() creates a bowtie1 index using
    $hpgl->{species}.fasta and leaves it in the indexes/ directory.

=cut
sub BT1_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args, required => ["species"],
                                   job_depends => "");
    my $job_basename = $options->{job_basename};
    my $job_string = qq!bowtie-build $options->{libdir}/$options->{libtype}/$options->{species}.fasta \\
  $options->{libdir}/$options->{libtype}/indexes/$options->{species}
!;
    my $comment = qq!## Generating bowtie1 indexes for species: $options->{species} in $options->{libdir}/$options->{libtype}/indexes!;
    my $bt1_index = $class->Submit(
        comment => $comment,
        job_name => "bt1idx",
        job_depends => $options->{job_depends},
        job_string => $job_string,
        job_prefix => "10",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
    );
    return($bt1_index);
}

=item C<BT2_Index>

    $hpgl->BT2_Index() creates a bowtie2 index using
    $hpgl->{species}.fasta and leaves it in the indexes/ directory.

=cut
sub BT2_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args, required => ["species"]);
    my $job_basename = $options->{job_basename};
    my $dep = "";
    $dep = $options->{job_depends};
    my $libtype = $options->{libtype};
    my $libdir = File::Spec->rel2abs($options->{libdir});
    my $job_string = qq!
if [ \! -r "${libdir}/genome/$options->{species}.fa" ]; then
  ln -s ${libdir}/genome/$options->{species}.fasta ${libdir}/genome/$options->{species}.fa
fi
if [ \! -r "${libdir}/genome/indexes/$options->{species}.fa" ]; then
  ln -s ${libdir}/genome/$options->{species}.fasta ${libdir}/genome/indexes/$options->{species}.fa
fi

bowtie2-build $options->{libdir}/genome/$options->{species}.fasta \\
  $options->{libdir}/${libtype}/indexes/$options->{species}
!;
    my $comment = qq!## Generating bowtie2 indexes for species: $options->{species} in $options->{libdir}/${libtype}/indexes!;
    my $indexer = $class->Submit(
        comment => $comment,
        job_depends => $dep,
        job_name => "bt2idx",
        job_prefix => "15",
        job_string => $job_string,
        prescript => $args{prescript},
        postscript => $args{postscript},
    );
    return($indexer);
}

=item C<BWA>

    $hpgl->BWA() performs a bwa alignment using both the sam(s|p)e and
    aln algorithms.  It then converts the output (when appropriate) to
    sorted/indexed bam and passes them to htseq.

=cut
sub BWA {
    my ($class, %args) = @_;
    my $check = which('bwa');
    die("Could not find bwa in your PATH.") unless($check);

    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => 'lmajor',
        libtype => 'genome',
    );

    my $sleep_time = 3;
    my %bwa_jobs = ();
    my $bwa_input = $options->{input};
    my $bwa_depends_on;
    $bwa_depends_on = $options->{job_depends} if ($options->{job_depends});
    my $job_basename = $options->{job_basename};

    my $job_name = qq"bwa_$options->{species}";
    $job_name = $options->{job_name} if ($options->{job_name});
    my $libtype = $options->{libtype};

    my $bwa_dir = qq"outputs/bwa_$options->{species}";
    $bwa_dir = $options->{bwa_dir} if ($options->{bwa_dir});

    my $uncompress_jobid = undef;
    my $index_jobid = undef;
    if ($bwa_input =~ /\.gz$|\.bz2$|\.xz$/) {
        print "The input needs to be uncompressed, doing that now.\n" if ($options->{debug});
        my $uncomp = Bio::Adventure::Compress::Uncompress(
            $class,
            input => $bwa_input,
            job_depends => $bwa_depends_on,
        );
        $bwa_input = basename($bwa_input, ('.gz', '.bz2', '.xz'));
        $bwa_jobs{uncompress} = $uncomp;
        $options = $class->Set_Vars(input => $bwa_input);
        $uncompress_jobid = $uncomp->{job_id};
    }

    ## Check that the indexes exist
    my $bwa_reflib = "$options->{libdir}/${libtype}/indexes/$options->{species}.fa";
    my $bwa_reftest = qq"${bwa_reflib}.sa";
    if (!-r $bwa_reftest) {
        my $index_job = Bio::Adventure::RNASeq_Map::BWA_Index(
            $class,
            job_depends => $bwa_depends_on,
            libtype => $libtype,
        );
        $bwa_jobs{index} = $index_job;
        $bwa_depends_on = $index_job->{job_id};
    }

    ## Make a depends string containing the uncompress, indexer, or both
    if (defined($index_jobid) && defined($uncompress_jobid)) {
        $bwa_depends_on = qq"${index_jobid}:${uncompress_jobid}";
    } elsif (defined($index_jobid)) {
        $bwa_depends_on = $index_jobid;
    } elsif (defined($uncompress_jobid)) {
        $bwa_depends_on = $uncompress_jobid;
    } else {
        $bwa_depends_on = "";
    }

    ## Include logic to deal with paired end/single end reads.
    my $test_file = "";
    my @pair_listing = ();
    my ($forward_reads, $reverse_reads);
    if ($bwa_input =~ /\s+/) {
        ($forward_reads, $reverse_reads) = split(/\s+/, $bwa_input);
    } else {
        $forward_reads = $bwa_input;
    }

    my $aln_sam = qq"${bwa_dir}/${job_basename}_aln.sam";
    my $mem_sam = qq"${bwa_dir}/${job_basename}_mem.sam";
    my $job_string = qq!mkdir -p ${bwa_dir}
bwa mem -a ${bwa_reflib} ${bwa_input} 2>${bwa_dir}/bwa.err 1>${mem_sam}
!;
    my $reporter_string = qq"bwa samse ${bwa_reflib} \\
  ${bwa_dir}/${job_basename}_aln-forward.sai $bwa_input \\
  1>${aln_sam} 2>${bwa_dir}/${job_basename}.samerr";
    my $aln_string = qq"bwa aln ${bwa_reflib} \\
  ${forward_reads} 1>${bwa_dir}/${job_basename}_aln-forward.sai \\
  2>${bwa_dir}/${job_basename}_aln-forward.err";
    if (defined($reverse_reads)) {
        $aln_string = qq"${aln_string}
bwa aln ${bwa_reflib} \\
  ${reverse_reads} 1>${bwa_dir}/${job_basename}_aln-reverse.sai \\
  2>${bwa_dir}/${job_basename}_aln-reverse.err";
        $reporter_string = qq"bwa sampe ${bwa_reflib} \\
  ${bwa_dir}/${job_basename}_aln-forward.sai ${bwa_dir}/${job_basename}_aln-reverse.sai \\
  ${forward_reads} ${reverse_reads} \\
  1>${aln_sam} 2>${bwa_dir}/${job_basename}.samerr";
    }

    my $comment = qq!## This is a BWA alignment of ${bwa_input} against
## ${bwa_reflib}.
## It will perform a separate bwa mem run and bwa aln run.!;

    $job_string = qq!${job_string}
${aln_string}
${reporter_string}
!;
    my $bwa_job = $class->Submit(
        comment => $comment,
        input => $bwa_input,
        job_depends => $bwa_depends_on,
        job_name => "bwa_$options->{species}",
        job_output => qq"${bwa_dir}/${job_basename}_mem.sam",
        job_prefix => "20",
        job_string => $job_string,
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        queue => 'workstation',
    );
    $bwa_jobs{bwa} = $bwa_job;

    my $mem_sam_job = Bio::Adventure::Convert::Samtools(
        $class,
        input => $bwa_job->{job_output},
        job_depends => $bwa_job->{job_id},
        job_name => "s2b_mem",
        job_prefix => "21",
    );
    $bwa_jobs{samtools_mem} = $mem_sam_job;

    my $aln_sam_job = Bio::Adventure::Convert::Samtools(
        $class,
        input => ${aln_sam},
        job_depends => $mem_sam_job->{job_id},
        job_name => "s2b_aln",
        job_prefix => "21",
    );
    $bwa_jobs{samtools_mem} = $aln_sam_job;

    my $mem_htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
        $class,
        htseq_id => $options->{htseq_id},
        htseq_input => $mem_sam_job->{job_output},
        htseq_type => $options->{htseq_type},
        job_depends => $mem_sam_job->{job_id},
        job_name => "htmem_${job_name}",
        job_prefix => '22',
    );
    $bwa_jobs{htseq_mem} = $mem_htmulti;

    my $aln_htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
        $class,
        htseq_id => $options->{htseq_id},
        htseq_input => $aln_sam_job->{job_output},
        htseq_type => $options->{htseq_type},
        job_depends => $aln_sam_job->{job_id},
        job_name => "htaln_${job_name}",
        job_prefix => '22',
    );
    $bwa_jobs{htseq_aln} = $aln_htmulti;

    my $bwa_stats = Bio::Adventure::RNASeq_Map::BWA_Stats(
        $class,
        job_depends => $mem_sam_job->{job_id},
        job_name => 'bwastats',
        job_prefix => "25",
        aln_output => $aln_sam_job->{job_output},
        mem_output => $mem_sam_job->{job_output},
    );
    $bwa_jobs{stats} = $bwa_stats;

    return(\%bwa_jobs);
}

=item C<BWA_Stats>

    $hpgl->BWA_Stats() collects some alignment statistics from bwa.

=cut
sub BWA_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $job_basename = $options->{job_basename};
    my $aln_input = $options->{aln_output};
    $aln_input = qq"${aln_input}.stats";
    my $mem_input = $options->{mem_output};
    $mem_input = qq"${mem_input}.stats";
    my $stat_output = qq"outputs/bwa_stats.csv";

    my $job_depends = $options->{job_depends};
    my $job_name = "bwa_stats";
    $job_name = $options->{job_name} if ($options->{job_name});
    my $jobid = qq"${job_basename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $job_string = qq!
if [ \! -r "${stat_output}" ]; then
    echo "# original reads, reads used, aln-aligned reads, mem-aligned reads, rpm" > ${stat_output}
fi
original_reads=0
if [ -r "outputs/trimomatic_stats.csv" ]; then
    original_reads_tmp=\$(tail -n 1 outputs/trimomatic_stats.csv | awk -F, '{print \$2}')
    original_reads=\${original_reads_tmp:-0}
fi
reads_tmp=\$(grep "^Total reads: " ${aln_input} | awk '{print \$3}' | sed 's/ //g')
reads=\${reads_tmp:-0}
aln_aligned_tmp=\$(grep "^Mapped reads" ${aln_input} | awk '{print \$3}' | sed 's/ .*//g')
aln_aligned=\${aln_aligned_tmp:-0}
mem_aligned_tmp=\$(grep "^Mapped reads" ${mem_input} | awk '{print \$3}' | sed 's/ .*//g')
mem_aligned=\${mem_aligned_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${aligned})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${job_basename},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aln_aligned}" "\${mem_aligned}" "\$rpm")
echo "\${stat_string}" >> ${stat_output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $aln_input,
        job_depends => $job_depends,
        job_name => $job_name,
        job_prefix => $options->{job_prefix},
        job_string => $job_string,
        cpus => 1,
        mem => 1,
        queue => "throughput",
        wall => "00:10:00",
    );
    return($stats);
}

=item C<BWA_Index>

    $hpgl->BWA_Index() creates bwa indexes.

=cut
sub BWA_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $job_basename = $options->{job_basename};
    my $job_string = qq!
if [ \! -r "$options->{libdir}/genome/$options->{species}.fa" ]; then
  ln -s $options->{libdir}/genome/$options->{species}.fasta $options->{libdir}/genome/$options->{species}.fa
fi
start=\$(pwd)
cd $options->{libdir}/$options->{libtype}/indexes &&
  bwa index $options->{species}.fa
cd \$start
!;
    my $comment = qq!## Generating bwa indexes for species: $options->{species} in $options->{libdir}/$options->{libtype}/indexes!;
    my $bwa_index = $class->Submit(
        comment => $comment,
        job_depends => $options->{job_depends},
        job_name => "bwaidx",
        job_prefix => "26",
        job_string => $job_string,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
    );
    return($bwa_index);
}

=item C<BT2_Stats>

    $hpgl->BT2_Stats() collects alignment statistics from bowtie 2.

=cut
sub BT2_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $bt_input = $options->{bt_input};
    my $job_basename = $options->{job_basename};
    my $bt_type = $options->{bt_type};
    my $job_name = "bt2_stats";
    $job_name = $options->{job_name} if ($options->{job_name});
    my $jobid = qq"${job_basename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $output = "outputs/bowtie2_stats.csv";
    my $job_string = qq!
if [ \! -r "${output}" ]; then
    echo "original reads, single hits, failed reads, multi-hits, rpm" > ${output}
fi
original_reads_tmp=\$(grep " reads; of these" "${bt_input}" 2>/dev/null | awk '{print \$1}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
one_align_tmp=\$(grep " aligned exactly 1 time" "${bt_input}" | awk '{print \$1}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep " aligned 0 times" "${bt_input}" | awk '{print \$1}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep " aligned >1 times" "${bt_input}" | awk '{print \$1}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \$(( \${one_align} + \${sampled} )) ) " 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${job_basename},${bt_type},%s,%s,%s,%s,%s" "\${original_reads}" "\${one_align}" "\${failed}" "\${sampled}" "\${rpm}")
echo "\$stat_string" >> ${output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $bt_input,
        job_name => $job_name,
        job_depends => $options->{job_depends},
        job_prefix => $args{job_prefix},
        job_string => $job_string,
        cpus => 1,
        mem => 1,
        queue => "throughput",
        wall => "00:10:00",
    );
    return($stats);
}

=item C<Kallisto_Index

    $hpgl->Kallisto_Index() uses kallisto and an annotated_CDS fasta sequence library to do the expected.

=cut
sub Kallisto_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species", "genome"],
    );
    my $job_basename = $options->{job_basename};
    my $dep = "";
    $dep = $options->{job_depends};
    my $libtype = $options->{libtype};
    my $genome = File::Spec->rel2abs($options->{genome});
    unless (-r $genome) {
        die("The indexing operation for kallisto will fail because the $options->{species} genome does not exist.")
    }

    my $job_string = qq!
kallisto index -i $options->{libdir}/${libtype}/indexes/$options->{species}.idx ${genome}!;
    my $comment = qq!## Generating kallisto indexes for species: $options->{species} in $options->{libdir}/${libtype}/indexes!;
    my $jobid = $class->Submit(
        comment => $comment,
        job_depends => $dep,
        job_string => $job_string,
        job_name => "kalidx",
        job_prefix => $options->{job_prefix},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
    );
    return($jobid);
}

=item C<Kallisto>

    $hpgl->Kallisto() should perform a kallisto alignment, I have not
    yet properly tested this  TODO!

=cut
sub Kallisto {
    my ($class, %args) = @_;
    my $check = which('kallisto');
    die("Could not find kallisto in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species","input"],
    );
    my %ka_jobs = ();
    my $ka_depends_on;
    $ka_depends_on = $options->{job_depends} if ($options->{job_depends});
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $job_basename = $options->{job_basename};
    my $species = $options->{species};
    my $ka_input = $options->{input};

    my $sleep_time = 3;
    my $job_name = qq"kall_${species}";
    $job_name = $options->{job_name} if ($options->{job_name});
    my $ka_args = qq"";

    if ($ka_input =~ /\.bz2$|\.xz$/ ) {
        my $uncomp = Bio::Adventure::Compress::Uncompress(
            $class,
            input => $ka_input,
            job_depends => $ka_depends_on,
        );
        $ka_input = basename($ka_input, ('.gz','.bz2','.xz'));
        $options = $class->Set_Vars(input => $ka_input);
        $ka_depends_on = $uncomp->{job_id};
    }
    my $input_name = $ka_input;
    $input_name = basename($input_name, ('.fastq'));
    my @input_test = split(/\:/, $ka_input);
    if (scalar(@input_test) == 1) {
        $ka_args .= " --single -l 40 -s 10 ";
    } elsif (scalar(@input_test) == 2) {
        ## Do nothing
        print "There are 2 inputs\n";
    } else {
        print "There are not 1 nor 2 inputs.\n";
    }
    $ka_input =~ s/\:/ /g;

    ## Check that the indexes exist
    my $ka_reflib = "$options->{libdir}/${libtype}/indexes/$options->{species}.idx";
    if (!-r $ka_reflib) {
        my $index_job = Bio::Adventure::RNASeq_Map::Kallisto_Index(
            $class,
            job_depends => $ka_depends_on,
            libtype => $libtype,
        );
        $ka_jobs{index} = $index_job;
        $ka_depends_on = $index_job->{job_id};
    }

    my $outdir = qq"outputs/kallisto_${species}";
    my $error_file = qq"${outdir}/kallisto_${species}.stderr";
    my $output_sam = qq"${outdir}/kallisto_${species}.sam";
    my $output_bam = qq"${outdir}/kallisto_${species}.bam";
    my $output_stats = qq"${outdir}/kallisto_${species}.stats";
    my $sorted_bam = qq"${outdir}/kallisto_${species}-sorted";
    my $comment = qq!## This is a kallisto pseudoalignment of ${ka_input} against
## ${ka_reflib}.
## This jobs depended on: ${ka_depends_on}
## Other candidates for making a pretty count table include:
##  perl -F'\\t' -a -n -e 'print "\$F[0] \$F[3]\\n"' ${outdir}/abundance.tsv > ${outdir}/abundance.count
##   awk '{printf("%s %s\\n", \$1, \$4)}' ${outdir}/abundance.tsv > ${outdir}/abundance.count
## The sam->bam conversion is copy/pasted from Sam2Bam() but I figured why start another job
## because kallisto is so fast
!;
    my $job_string = qq!mkdir -p ${outdir} && sleep ${sleep_time} && \\
kallisto quant ${ka_args} --plaintext --pseudobam -t 4 -b 100 -o ${outdir} -i ${ka_reflib} \\
  ${ka_input} 2>${error_file} 1>${output_sam} && \\
  cut -d "	" -f 1,4 ${outdir}/abundance.tsv > ${outdir}/${input_name}_abundance.count && \\
  gzip ${outdir}/${input_name}_abundance.count && \\
  samtools view -u -t ${ka_reflib} -S ${output_sam} 1>${output_bam} && \\
  samtools sort -l 9 ${output_bam} ${sorted_bam} && \\
  rm ${output_bam} && mv ${sorted_bam}.bam ${output_bam} && samtools index ${output_bam} && \\
  bamtools stats -in ${output_bam} 2>${output_stats} 1>&2
!;
    my $ka_job = $class->Submit(
        comment => $comment,
        input => $ka_input,
        job_depends => $ka_depends_on,
        job_name => qq"${job_name}",
        job_prefix => "30",
        job_string => $job_string,
        prescript => $args{prescript},
        postscript => $args{postscript},
        queue => "workstation",
    );
    $ka_jobs{kallisto} = $ka_job;

    return(\%ka_jobs);
}

=item C<Tophat>

    $hpgl->Tophat() runs... guess... tophat!  It also sorts/indexes
    the accepted_hits file, collects some statistics, and passes the
    hits to htseq-count.

=cut
sub Tophat {
    my ($class, %args) = @_;
    my $check = which('tophat');
    die("Could not find tophat in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species", "input", "htseq_type"],
    );
    my $job_depends = $options->{job_depends};
    my $tophat_cpus = 4;
    my $inputs = $options->{input};
    my @in = split(/:/, $inputs);
    $inputs =~ s/:/ /g;
    my $number_inputs = scalar(@in);
    my $tophat_args = ' -g 1 --microexon-search --b2-very-sensitive ';
    if ($options->{tophat_args}) {
        $tophat_args = $options->{tophat_args};
    }
    ## $tophat_args .= ' --no-mixed --no-discordant ' if (scalar(@in) > 1);
    ## $tophat_args .= ' ' if (scalar(@in) > 1);

    my $tophat_queue = $options->{queue};
    my $tophat_walltime = '18:00:00';
    my $tophat_mem = 8;
    if ($options->{species} eq 'hsapiens' or $options->{species} eq 'mmusculus') {
        $tophat_queue = 'workstation';
        $tophat_walltime =  '144:00:00';
        $tophat_mem = 20;
    }

    my $tophat_dir = qq"outputs/tophat_$options->{species}";
    if ($options->{tophat_dir}) {
        $tophat_dir = $options->{tophat_dir};
    }
    my $job_basename = $options->{job_basename};
    my $libtype = $options->{libtype};
    my $bt_reflib = "$options->{libdir}/${libtype}/indexes/$options->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    my $index_job = undef;
    if (!-r $bt_reftest) {
        print "Did not find the index for $options->{species} at: ${bt_reflib}, indexing now.\n";
        $index_job = Bio::Adventure::RNASeq_Map::BT2_Index(
            $class,
            job_depends => $job_depends,
        );
        ## Use a colon separated append to make tophat depend on multiple jobs
        ## $job_depends .= qq":$index_job->{jobid}";
        ## Or just replace the dependency string with this job's and insert it into the stack
        $job_depends = $index_job->{job_id};
    }
    my $gtf_file = qq"$options->{libdir}/genome/$options->{species}.gtf";
    if (!-r $gtf_file) {
        print "Missing the gtf file for $options->{species}\n";
        print "Using the gff file.\n";
        $gtf_file =~ s/\.gtf/\.gff/;
        ##my $written = $class->Gff2Gtf(gff => "$class->{libdir}/genome/$class->{species}.gff");
    }

    my $spliced = 0;
    if ($options->{spliced}) {
        $spliced = 1;
    }
    my $job_name = qq"th_$options->{species}";
    my $job_string = qq!
mkdir -p ${tophat_dir} && tophat ${tophat_args} \\
  -G ${gtf_file} \\
  -p ${tophat_cpus} -o ${tophat_dir} \\
!;
    if ($spliced) {
        $job_string .= qq!  --no-novel-juncs \\
!;
    }
    $job_string .= qq!  $options->{libdir}/genome/indexes/$options->{species} \\
  ${inputs} 2>outputs/tophat.out 1>&2 && \\
 samtools sort -l 9 -n ${tophat_dir}/accepted_hits.bam -o ${tophat_dir}/accepted_sorted.bam && \\
 samtools index ${tophat_dir}/accepted_hits.bam && \\
 samtools sort -l 9 -n ${tophat_dir}/unmapped.bam -o ${tophat_dir}/unmapped_sorted.bam && \\
 samtools index ${tophat_dir}/unmapped.bam
!;
    if ($number_inputs > 1) {
        $job_string .= qq!
if [ -r "${tophat_dir}/accepted_hits.bam" ]; then
  samtools view -b -f 2 ${tophat_dir}/accepted_hits.bam > ${tophat_dir}/accepted_paired.bam && \\
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
        comment => ${comment},
        cpus => ${tophat_cpus},
        job_depends => ${job_depends},
        job_name => ${job_name},
        job_prefix => "31",
        job_string => ${job_string},
        mem => ${tophat_mem},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => ${tophat_queue},
        wall => ${tophat_walltime},
    );

    ## Set the input for htseq
    my $accepted = "${tophat_dir}/accepted_hits.bam";
    $accepted = $options->{accepted_hits} if ($options->{accepted_hits});
    my $count_table = "accepted_hits.count";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
        $class,
        htseq_input => $accepted,
        htseq_id => $options->{htseq_id},
        htseq_type => $options->{htseq_type},
        job_depends => $tophat->{job_id},
        job_name => qq"hts_$options->{species}",
        job_prefix => "32",
    );
    $tophat->{htseq} = $htmulti;
    ## Perform a separate htseq run using only the successfully paired hits
    if ($number_inputs > 1) {
        my $ht_paired = Bio::Adventure::RNASeq_Count::HT_Multi(
            $class,
            htseq_input => qq"${tophat_dir}/accepted_paired.bam",
            htseq_id => $options->{htseq_id},
            htseq_type => $options->{htseq_type},
            job_depends => $tophat->{job_id},
            job_name => qq"htsp_$options->{species}",
            job_prefix => "32",
        );
    }
    ## Tophat_Stats also reads the trimomatic output, which perhaps it should not.
    my $unaccepted = $accepted;
    $unaccepted =~ s/accepted_hits/unmapped/g;
    my $input_read_info = $accepted;
    $input_read_info =~ s/accepted_hits\.bam/prep_reads\.info/g;
    my $stats = Bio::Adventure::RNASeq_Map::Tophat_Stats(
        $class,
        accepted_input => $accepted,
        count_table => qq"${count_table}.xz",
        job_depends => $tophat->{job_id},
        job_name => "tpstats_$options->{species}",
        job_prefix => "33",
        prep_input => $input_read_info,
        unaccepted_input => $unaccepted,
    );
    $tophat->{stats} = $stats;

    return($tophat);
}

=item C<BT1_Stats>

    $hpgl->BT1_Stats() collects some alignment statistics from
    bowtie1.

=cut
sub BT1_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $bt_input = $options->{bt_input};
    my $job_basename = $options->{job_basename};
    my $bt_type = "";
    $bt_type = $options->{bt_type} if ($options->{bt_type});
    my $job_depends = $options->{job_depends};
    my $job_name = "stats";
    $job_name = $options->{job_name} if ($options->{job_name});
    my $jobid = qq"${job_basename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $stat_output = qq"outputs/bowtie_stats.csv";
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $job_string = qq!
if [ \! -r "${stat_output}" ]; then
  echo "name,type,original_reads,reads,one_hits,failed,samples,rpm,count_table" > ${stat_output}
fi
original_reads=0
if [ -r "outputs/trimomatic_stats.csv" ]; then
  original_reads_tmp=\$(tail -n 1 outputs/trimomatic_stats.csv | awk -F, '{print \$2}')
  original_reads=\${original_reads_tmp:-0}
fi
reads_tmp=\$(grep "^# reads processed" ${bt_input} | awk -F: '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
one_align_tmp=\$(grep "^# reads with at least one reported" ${bt_input} | awk -F": " '{print \$2}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep "^# reads that failed to align" ${bt_input} | awk -F": " '{print \$2}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep "^# reads with alignments sampled" ${bt_input} | awk -F": " '{print \$2}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${one_align})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${job_basename},${bt_type},%s,%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${one_align}" "\${failed}" "\${sampled}" "\$rpm")
echo "\$stat_string" >> ${stat_output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $bt_input,
        job_depends => $job_depends,
        job_name => $job_name,
        job_prefix => $args{job_prefix},
        job_string => $job_string,
        cpus => 1,
        mem => 1,
        queue => "throughput",
        wall => "00:10:00",
    );
    return($stats);
}

=item C<Tophat_Stats>

    $hpgl->Tophat_Stats() collects alignment statistics from the
    accepted_hits.bam/unaligned.bam files generated by a tophat run.

=cut
sub Tophat_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $accepted_input = $options->{accepted_input};
    my $accepted_output = qq"${accepted_input}.stats";
    my $unaccepted_input = $options->{unaccepted_input};
    my $unaccepted_output = qq"${unaccepted_input}.stats";
    my $read_info = $options->{prep_input};
    my $job_basename = $options->{job_basename};
    my $job_depends = $options->{job_depends};

    my $job_name = "stats";
    $job_name = $options->{job_name} if ($options->{job_name});
    my $jobid = qq"${job_basename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $output = "outputs/tophat_stats.csv";
    my $comment = qq!## This is a stupidly simple job to collect tophat alignment statistics.!;
    my $job_string = qq!
if [ \! -r "${output}" ]; then
  echo "basename,species,original_reads,aligned_reads,failed_reads,rpm,count_table" > ${output}
fi
bamtools stats < "${accepted_input}" \\
    2>${accepted_output} 1>&2 && \\
  bamtools stats < "${unaccepted_input}" \\
    2>${unaccepted_output} 1>&2

original_reads=0
if [ -r "outputs/trimomatic_stats.csv" ]; then
  original_reads_tmp=\$(tail -n 1 outputs/trimomatic_stats.csv | awk -F, '{print \$2}')
  original_reads=\${original_reads_tmp:-0}
fi
reads_tmp=\$(grep "^reads_in " ${read_info} | awk -F= '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
aligned_tmp=\$(grep "^Total reads" ${accepted_output} | awk '{print \$3}' | sed 's/ .*//g')
aligned=\${aligned_tmp:-0}
failed_tmp=\$(grep "^Total reads" ${unaccepted_output} | awk '{print \$3}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${aligned})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${job_basename},$options->{species},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aligned}" "\${failed}" "\$rpm")
echo "\$stat_string" >> "${output}"!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $accepted_input,
        job_name => $job_name,
        job_depends => $job_depends,
        job_prefix => $args{job_prefix},
        job_string => $job_string,
        mem => 1,
        queue => "throughput",
        wall => "00:10:00",
    );
    return($stats);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;