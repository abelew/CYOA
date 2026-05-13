package Bio::Adventure::Job_Estimators;
use strict;

=head2 C<Estimate_Memory_Single>

Estimate the amount of memory required for a job given a single parameter and a constant.
Ideally I would like to have a file comprised of a set of functions which define the required
memory in a fashion like:
  memory = { input_file_size * constant } or
  memory = { input_lines * constant } etc

=cut
sub Estimate_Memory_Single {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars();
    ## In order to collect helpful memory usage numbers on a slurm
    ## cluster, this kind of requires that the slurmdbd maintains the
    ## state of previous jobs.

    ## I thus ran cyoa --checkjob in an arbitrarily chosen recent
    ## processing directory
    ## (rnaseq/mmusculus_iprgc_2024/preprocessing/umd_sequenced/something)
    ## and looked over the outputs/logs/jobs.csv and discovered that
    ## for a sample with r1/r2 of 997M/1.2G sizes:
    ##
    ## 1. umi_extract took 5 G memory at maximum, 23 minutes (always
    ## rounding up)
    ## 2. hisat2 vs. the mm39_112 genome took 37G memory,
    ## 18 minutes
    ## 3. samtools took 32G memory, 111 minutes
    ## 4. htseq took 1G memory, 23 minutes
    ## 5. hisat2 vs the mouse rRNA database took 24G memory, 14 minutes.
    ## 6. samtools took 26G and 106 minutes
    ## 7. htseq took 100Mb and 1 minute
    ## 8. umi_dedup took more memory than the 20G I allocated and failed.

    ## Nonetheless, I can use these numbers to get started on some memory estimates.
    my $multipliers = {
        umi_extract => 1.2,  ## Sum the input fastq and multiply.
        htseq => 4, ## Take the input gff and multiply
        samtools => 8, ## Take the input fastq files and multiply (this _should_ be using the sam,
        ## but I don't know its size)
    };

    ## I am thinking to pass this function an option 'mem_estimator'
    ## which will be the path to whatever we choose to examine.
    my $total_mem = 0.1;  ## Start with 0.1G ram and move up from
    ## there.  I actually do not know if slurm requires integer memory
    ## amounts.  If so we will just do a ceil...

    my $job_type = $options->{mem_job};
    my $job_multiplier = $multipliers->{$job_type};
    if (!defined($job_multiplier)) {
        print "No memory multiplier defined for this job type, assuming 10x the input file size.\n";
        $job_multiplier = 10;
    }

    my @inputs = split(/[:;,\s+]/, $options->{mem_estimator});
    for my $input (@inputs) {
        my $st = stat($input);
        my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size,
            $atime, $mtime, $ctime, $blksize, $blocks) = stat($input);
        my @pieces = split(/\./, $input);
        my $ext = $pieces[$#pieces];
        if (defined($multipliers->{$ext})) {
            $size = $size * $job_multiplier;
        }
        $total_mem = $total_mem + $size;
    }
    return($total_mem);
}

=head2 C<Estimate_Time_Single>

Analagous to Estimate_Memory_Single, but the required time for a given job assuming it may be
estimated via a single parameter and some constant.  If I am going to actually use this, it
will require some logic from the cluster informing the submitting job of the slowest machine
in the cluster; or I will need to just assume that machines will continue getting faster and
therefore the times I have already recorded are thus the slowest.

=cut
sub Estimate_Time_Single {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        process => 'xz',
        species => q{},
    );

    ## Approximate compression multipliers, taken from:
    ## https://community.centminmod.com/threads/round-3-compression-comparison-benchmarks-zstd-vs-brotli-vs-pigz-vs-bzip2-vs-xz-etc.17259/
    my $multipliers = {
        xz => 4.4,
        bz2 => 3.9,
        gz => 3.1,
    };
    my $mb_per_sec = {
        gz => 12.0 * 211957760 * 10000,
        bz2 => 12.9 * 211957760 * 10000,
        xz => 2.3 * 211957760 * 10000,
    };
    ## I need to acquire from slurm the relative cpu frequency, but assuming similar speeds
    ## I should be able to get reasonable time estimates...
    my @inputs = split(/[:;,\s+]/, $options->{input});
    my $total_hours = 0;
    for my $input (@inputs) {
        my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size,
            $atime, $mtime, $ctime, $blksize, $blocks) = stat($input);
        my @pieces = split(/\./, $input);
        my $ext = $pieces[$#pieces];
        if (defined($multipliers->{$ext})) {
            $size = $size * $multipliers->{$ext};
        }
        my $multi_key = $options->{process};
        if (defined($mb_per_sec->{$multi_key})) {
            my $needed = $size / $mb_per_sec->{$ext};
            my $hours = ceil($needed * 60.0 * 60.0);
            $total_hours = $total_hours + $hours;
        }
    }
    return($total_hours);
}

=head2 C<Estimate_Time_Composite>

Estimate required job time for occasions when multiple parameters are relevant.  Things like
htseq are good examples: a combination of input filesize and the complexity of the reference genome
define how long the job will take.

=cut
sub Estimate_Time_Composite {
    my ($class, %args) = @_;
    my $options = $args{options};
    my $paths = $args{paths};
    my $inputs = $args{inputs};
    my $species_file = qq"$options->{libpath}/$options->{libtype}/fasta/$options->{species}.fasta";
    my $species_size = 1;
    if (!-r $species_file) {
        print "The species file does not exist!\n";
    } else {
        my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size,
        $atime, $mtime, $ctime, $blksize, $blocks) = stat($species_file);
        $species_size = $size / 1e7;
    }

    ## I need to acquire from slurm the relative cpu frequency, but assuming similar speeds
    ## I should be able to get reasonable time estimates...
    ## I am going to assume that the less preprocessor will make xz/gz/bz/uncomp
    ## input file size differences negligible with respect to time to run.
    my $total_input_size = 0;
    for my $in (@{$inputs}) {
        if (!-r $in->{fullpath}) {
            die("This file does not exist: ${in}.");
        } else {
            my ($dev, $ino, $mode, $nlink, $uid, $gid, $rdev, $size,
                $atime, $mtime, $ctime, $blksize, $blocks) = stat($in->{fullpath});
            $total_input_size = $total_input_size + $size;
        }
    }

    # A genome multiplier should be used when the time estimate for a tool can be described as:
    ## (input_size * constant) * (genome_size * constant)
    ## Really it should be the index size, but I assume the fasta is a sufficient proxy.
    my $genome_mult = {
        bt1 => 1.0,
        hisat => 1.0,
    };
    ## Approximate compression multipliers, taken from:
    ## https://community.centminmod.com/threads/round-3-compression-comparison-benchmarks-zstd-vs-brotli-vs-pigz-vs-bzip2-vs-xz-etc.17259/
    my $multipliers = {
        bt1 => 1.0,
        hisat => 4.4,
    };

    ## My current thinking is:
    ## ((inputsize * constant) * (genomesize * constant)) -> seconds to finish job.
    ## Once I get that nailed down, then I will scale the result by the slurm cpu frequency values.
    ## For an example, a 300Mb (xz compressed) paired end sample against a 6.4Mb genome took 13 minutes on one of the old cbcb nodes.
    ## Conversely, a 980Mb (xz) paired end sample against a 12Gb genome took 12 (720 seconds) minutes on one of the new nodes.
    ## So, this is a case where I think we can assume IO is not limiting but cpu frequency is...
    ## A more relevant comparison, another smaller sample on the new nodes:
    ## 250Mb input (xz) paired end vs a 1.8Mb genome took 1 minute (60 seconds)
    ## (/ 720.9 12000000.0) -> 6.0e-5
    ## (/ 60 (* 250 18.0)) -> 0.013
    my $key = $args{process};
    $total_input_size = $total_input_size / 1e7;
    my $needed_input_mult = $total_input_size * $multipliers->{$key};
    my $needed_genome_mult = $species_size * $genome_mult->{$key};
    my $needed = ($needed_input_mult * $needed_genome_mult);
    my $total_hours = ceil($needed * 60.0 * 60.0);
    return($total_hours);
}

=head2 C<Get_Memory>

In my time estimator, I tried out some heuristics and have found them
somewhat difficult to rely upon because of the huge variance in cpu speeds.
I thought I would therefore try some simpler heuristics for memory usage
especially since memory usage tends to be a bit more consistent
which of course may be an argument for a multiplier based system.

=cut
sub Get_Memory {
    my ($class, %args) = @_;
    my $function = 'Default';
    my %memory_data = (
        'Example' => {
            genes => {
                '1000' => 4,
                '10000' => 12,
                '20000' => 24,
                'default' => 30,
            },
        },
        'Example2' => {
            genome_size => {
                '1000000' => 4,
                '10000000' => 8,
                '100000000' => 12,
                'default' => 24,
            },
        },
    );

    my @function_array = split(/::/, $function);
    $function = $function_array[$#function_array];
    $function =~ s/_Worker$//g;
    my $datum = $memory_data{$function};
    my $needed_memory = 6;
    if ($datum) {
        my @expected_args = keys(%$datum);
        my $expected = shift @expected_args;
        my $actual = $args{$expected};
        my @cutoff_values = keys %{$datum->{$expected}};
        my $found_cutoff = 0;
        my $used_default = 0;
      CUTOFF: for my $cut (@cutoff_values) {
            if ($actual <= $cut) {
                $needed_memory = $datum->{$expected}->{$cut};
                $found_cutoff = 1;
                last CUTOFF;
            }
            if ($cut eq 'default') {
                $needed_memory = $datum->{$expected}->{$cut};
                $found_cutoff = 1;
                $used_default = 1;
                last CUTOFF;
            }
        }
        if ($used_default) {
            print "Used the default value: ${needed_memory}.\n";
        }
        unless ($found_cutoff) {
            print "Did not find the expected cutoff type: ${expected}, using the global default: ${needed_memory}.\n";
        }
    } else {
        print "This function does not have a memory configuration, using the global default: ${needed_memory}.\n";
    }
    return($needed_memory);
}

=head2 C<Get_Modules>

Figure out the appropriate environment module/condat environment for a job.

This is the first place I will move my code into a json/ini/etc file.  Each entry of the hash
is a function with either a list of modules or environments required for it to work.
If the exe key is provided, then it will check to see that the program exists and is executable.
Currently that check is located elsewhere, perhaps move that here?

=cut
sub Get_Modules {
    ## Once all of the information from Old_Get_Modules is moved into
    ## functions.json, delete this and it.
    say "I exist now only to fill in job specifications and be deleted.";
}

sub Old_Get_Modules {
    my %args = @_;
    my $function = $args{function};
    ## I think that with recent changes, one need no longer add 'cyoa' as a module.
    my %module_data = (
        'Example' => {
            modules => ['abricate'],
            conda => ['torsten'],
            exe => 'example' },
        'Abricate' => { ## Adding cyoa because Path::Tiny isn't in the abricate directory yet.
            modules => ['cyoa', 'any2fasta', 'abricate', 'blast', 'blastdb',],
            exe => 'abricate' },
        'Abyss' => { modules => 'abyss' },
        'Angsd_Filter' => { modules => 'angsd' },
        'Any2Any' => { modules => 'cyoa' },
        'Aragorn' => { modules => 'aragorn', exe => 'aragorn' },
        'Assembly_Coverage' => {
            modules => ['cyoa', 'hisat2', 'samtools', 'bbmap',], },
        'Bacphlip' => { modules => 'bacphlip', exe => 'bacphlip' },
        'Bam2Coverage' => { modules => ['samtools', 'bbmap',] },
        'Bam2Bigwig' => {modules => ['deeptools',] },
        'Bedtools_Coverage' =>  { modules => 'bedtools', exe => 'bedtools' },
        'Biopieces_Graph' => { modules => ['biopieces'] },
        'Bowtie' => { modules => 'bowtie1', },
        'Bowtie2' => { modules => 'bowtie2' },
        'BT1_Index' => { modules => ['bowtie1'] },
        'BT2_Index' => { modules => ['bowtie2'] },
        'BWA' => { modules => 'bwa' },
        'BWA_Index' => { modules => ['bwa'] },
        'Caical' => { modules => ['cyoa', 'caical'] },
        'Casfinder' => { modules => ['casfinder',], exe => 'casfinder.sh' },
        'Casoff' => { modules => 'casoff', exe => 'cas-offinder' },
        'CGView' => { modules => ['cyoa', 'cgview'], }, ## cyoa has the prereqs required
        'Check_Blastdb' => { modules => ['blast'], },
        'Classify_Phage' => {
            modules => ['cyoa', 'blast', 'blastdb'], exe => 'tblastx' },
        'Cogent' => { modules => ['cogent'], exe => 'cogent'},
        'Collect_Assembly' => { modules => 'cyoa' },
        'Consolidate_TAs' => { modules => 'cyoa', },
        'Cutadapt' => { modules => ['cyoa', 'cutadapt'], exe => 'cutadapt'},
        'Dantools_RNASeq' => { modules => ['dantools'] },
        'Downsample_Guess_Strand' => { modules => ['salmon'] },
        'Essentiality_TAs' => { modules => 'cyoa', },
        'Extend_Kraken_DB' => { modules => ['kraken'], exe => ['kraken2'] },
        'Fasta2Gff' => { modules => ['cyoa'] },
        'Fastqc' => { modules => ['fastqc'] },
        'Fastq_Dump' => { modules => ['cyoa', 'sra'], },
        'Fastp' => { modules => ['fastp'], exe => 'fastp' },
        'Feature_Counts' => { modules => ['subread'], exe => 'featureCounts', },
        'Filter_Host_Kraken' => {
            modules => ['cyoa', 'kraken', 'hisat2', 'htseq', 'samtools'], },
        'Fragment_Genome' => {
            modules => ['cyoa' ], },
        'Freebayes_SNP_Search' => {
            modules => ['gatk', 'freebayes', 'samtools', 'bcftools', 'vcftools'],
            exe => ['gatk', 'freebayes'], },
        'GATK_Dedup' => { modules => 'gatk' },
        'Gb2Gff' => { modules => 'cyoa' },
        'Generate_Samplesheet' => { modules => 'R' },
        'GFF_Tree' => { modules => ['cyoa'] },
        'Glimmer' => { modules => ['cyoa', 'glimmer'], exe => 'glimmer3' },
        'Glimmer_Single' => {
            modules => ['glimmer', 'cyoa'], exe => 'glimmer3' },
        'Gubbins' => { modules => 'gubbins' },
        'Guess_Strand' => { modules => ['cyoa', 'samtools'] },
        'Hisat2' => { modules => 'hisat2', exe => 'hisat2'},
        'Hisat2_Index' => { modules => ['hisat2'], exe => ['hisat2'], },
        'HT_Multi' => { modules => 'htseq' },
        'HTSeq' => { modules => 'htseq', exe => 'htseq-count' },
        'Insert_Size' => { modules => ['gatk', 'R'], },
        'Interpro_Long2Wide' => { modules => 'cyoa' },
        'Interproscan' => {
            modules => ['cyoa', 'interproscan'], exe => 'interproscan.sh' },
        'Jellyfish' => { modules => ['jellyfish'], exe => 'jellyfish' },
        'Jellyfish_Matrix' => { modules => ['cyoa'] },
        'Kallisto' => { modules => 'kallisto', exe => 'kallisto' },
        'Kallisto_Index' => { modules => ['kallisto'], exe => ['kallisto'], },
        'Kraken' => { modules => 'kraken', exe => 'kraken2', },
        'Kraken_Best_Hit' => { modules => ['cyoa', 'kraken'] },
        'Make_Blast_Job' => { modules => ['cyoa', 'blastdb', 'blast'], exe => 'tblastx' },
        'Make_Fasta_Job' => { modules => ['cyoa', 'fasta'] },
        'Make_Codon_Table' => { modules => ['cyoa'] },
        'Mash' => { modules => 'mash', exe => 'mash' },
        'MauveAligner' => { modules => 'mauve' },
        'Merge_Annotations' => { modules => ['cyoa', 'ncbi_tools/6.1'], exe => 'tbl2asn' },
        'Merge_CDS_Predictions' => { modules => ['cyoa', 'ncbi_tools/6.1'], exe => 'tbl2asn' },
        'Merge_Parse_Blast' => { modules => ['cyoa'], },
        'Mpileup' => { modules => ['samtools', 'gatk'] },
        'Mpileup_SNP_Search' => {
            modules => ['samtools/1.13', 'bcftools', 'gatk', 'vcftools'],
            exe => ['samtools', 'bcftools'], },
        'Nextflow_RNASeq' => { modules => ['nextflow'], },
        'OrthoFinder' => { modules => ['cyoa', 'orthofinder'], exe => ['orthofinder'] },
        'OrthoMCL_Pipeline' => { modules => ['orthomcl'], exe => ['orthomclPairs'] },
        'Pairwise_Similarity_Matrix' => { modules => ['blast'], },
        'Parse_Fasta_Global' => { modules => 'cyoa' },
        'Phagepromoter' => { modules => 'phagepromoter' },
        'Phageterm' => { modules => ['cyoa', 'phageterm'], exe => 'PhageTerm.py' },
        'Phanotate' => { modules => ['trnascan', 'phanotate'], exe => 'phanotate.py' },
        'Phastaf' => { modules => ['cyoa', 'phastaf'] },
        'PolyA_Extractor_Worker' => { modules => ['cyoa'] },
        'PolyA_Recorder' => { modules => ['cyoa', 'samtools'] },
        'Prodigal' => { modules => ['cyoa', 'prodigal'] },
        'Progressive_Mauve' => { modules => 'mauve' },
        'Prokka' => { ## Prokka should not need cyoa; it is getting requisite perl module from it for now
            modules => ['cyoa', 'prokka', 'blast'], exe => 'prokka'},
        'ProteinFold' => { modules => ['cyoa', 'alphafold3', 'cuda'] },
        'ProteinFold_PairIDs' => { modules => ['cyoa', 'alphafold3', 'cuda'] },
        'Racer' => { modules => ['hitec'], exe => ['RACER'], },
        'Resfinder' => { modules => 'resfinder', exe => 'run_resfinder.py' },
        'Restriction_Catalog' => { modules => 'cyoa' },
        'Rgi' => { modules => ['kma', 'jellyfish', 'bowtie', 'bwa', 'diamond', 'rgi'], },
        'Rho_Predict' => { modules => 'rhotermpredict' },
        'RMats' => { modules => ['rmats'] },
        'RNAFold_Windows' => { modules => ['cyoa', 'vienna'], exe => 'RNAfold' },
        'Rosalind_Plus' => { modules => ['cyoa', 'prodigal'] },
        'RSEM' => { modules => 'rsem', exe => 'rsem-prepare-reference' },
        'RSEM_Index' => { modules => ['rsem', 'bowtie2'], exe => ['bowtie2'] },
        'Run_Essentiality' => { modules => 'essentiality' },
        'Run_Parse_Blast' => { modules => ['cyoa', 'blast'], exe => 'blastp' },
        'Salmon' => { modules => 'salmon', exe => 'salmon', },
        'Salmon_Index' => { modules => ['salmon'], exe => ['salmon']},
        'Samtools' => { modules => ['samtools', 'bamtools'], exe => 'samtools' },
        'Shovill' => { modules => 'shovill', exe => 'shovill' },
        'SLSearch' => { modules => ['cutadapt'] },
        ## FIXME: I should only set one of these two (and for the following ones too)
        'SL_Recorder' => { modules => ['cyoa', 'samtools'] },
        'SLSeq_Recorder_Worker' => { modules => ['cyoa', 'samtools'] },
        'SL_UTR' => { modules => 'cyoa' },
        'Snippy' => { modules => ['snippy', 'gubbins', 'fasttree'],
                      exe => ['snippy', 'gubbins'], },
        'SNP_Ratio' => {
            modules => ['cyoa', 'freebayes', 'gatk',
                        'samtools', 'bcftools', 'vcftools'],
            exe => ['samtools', 'cyoa', 'freebayes'] },
        'SNP_Ratio_Intron' => {
            modules => ['cyoa', 'freebayes', 'gatk',
                        'samtools', 'bcftools', 'vcftools'], },
        'SNP_Ratio_Worker' => {
            modules => ['cyoa', 'freebayes', 'samtools',
                        'bcftools', 'vcftools'], },
        'SNP_Ratio_Intron_Worker' => {
            modules => ['cyoa', 'freebayes', 'samtools',
                        'bcftools', 'vcftools'], },
        'Sort_Indexes' => { modules => 'cyoa' },
        'Spladder' => { modules => 'spladder' },
        'Split_Align_Blast' => { modules => ['cyoa', 'blast', 'blastdb'], exe => 'blastn' },
        'Split_Align_Fasta' => { modules => ['cyoa', 'fasta',], exe => 'fasta36' },
        'Split_Fasta' => { modules => ['cyoa'] },
        'Spring' => { modules => ['spring',] },
        'STAR' => { modules => 'star' },
        'STAR_Index' => { modules => ['star'] },
        'Suppa' => { modules => ['suppa'] },
        'TA_Check' => { modules => 'cyoa' },
        'Terminase_ORF_Reorder' => { modules => ['cyoa', 'fasta', 'blast', 'blastdb'], },
        'Test' => { modules => ['bcftools'], exe => 'bcftools', },
        'Train_Prodigal' => { modules => 'prodigal' },
        'Tophat' => { modules => 'tophat' },
        'Transdecoder' => { modules => 'transdecoder' },
        'Transit_TPP' => { modules => ['bwa', 'transit', 'htseq' ], exe => 'tpp' },
        'Transposonpsi' => { modules => 'transposonpsi', exe => 'transposonPSI.pl' },
        ## I can set the exe for single vs pairwise now
        'Trimomatic_Pairwise' => { modules => ['trimomatic'], },
        'Trimomatic_Single' => { modules => ['trimomatic'], },
        'Trinity' => { modules => ['trinity', 'cyoa'], exe => 'Trinity', },
        'Trinity_Post' => { modules => ['cyoa', 'rsem'], },
        'Trinotate' => {
            modules => ['cyoa', 'divsufsort', 'transdecoder', 'blast', 'blastdb', 'signalp/4.1',
                        'hmmer', 'tmhmm', 'rnammer', 'trinotate',],
            exe => ['Trinotate'] },
        'tRNAScan' => { modules => ['infernal', 'trnascan'], exe => 'trnascan' },
        'Unicycler' => {
            modules => ['trimomatic', 'bowtie2', 'spades', 'unicycler',
                        'flash', 'shovill', 'bwa', 'pilon'],
            exe => 'unicycler', },
        'Umi_Tools_Dedup' => { modules => ['umi_tools', 'samtools'], },
        'Umi_Tools_Extract' => { modules => 'umi_tools', },
        'Unicycler_Filter_Depth' => { modules => 'cyoa', },
        'Velvet' => { modules => 'velvet', exe => 'velveth' },
        'Xref_Crispr' => { modules => 'cyoa', },
    );
    my $datum = $module_data{$function};
    if ($datum) {
        ## Sometimes I don't bother to make the module list
        ## a list, but just a scalar.
        my $module_type = ref($datum->{modules});
        my @mod_lst;
        my @conda_lst;
        if (!$module_type) {
            my $mod_name = $datum->{modules};
            if ($mod_name) {
                push(@mod_lst, $mod_name);
                $datum->{modules} = \@mod_lst;
            } else {
                $datum = {};
            }
            my $conda_name = $datum->{conda};
            if ($conda_name) {
                push(@conda_lst, $conda_name);
                $datum->{conda} = \@conda_lst;
            }
        }
    } else {
        print "The function: ${function} has no modules, adding cyoa.\n";
        $datum = {
            modules => ['cyoa'],
        };
    }
    return(%{$datum});
}

1;
