package Bio::Adventure::Config;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Exporter;
use File::Basename qw"basename dirname";
use File::Path qw"make_path";
use File::stat;
use Getopt::Long qw"GetOptionsFromArray";

our @EXPORT = qw"Get_Modules Get_Menus Get_TODOs";
our @ISA = qw"Exporter";

sub Get_Caller {
    my %args = @_;
    my $caller_number = 2;
    $caller_number = $args{caller} if ($args{caller});
    my $attempts = 4;
    my $found_sub = undef;
    my $subroutine;
    my $last_attempt = $caller_number + $attempts;
  COUNTER: for my $attempt ($caller_number .. $last_attempt) {
        my ($package, $filename, $line, $sub, $hasargs,
            $wantarray, $evaltext, $is_require, $hints, $bitmask, $hinthash) = caller($attempt);
        if ($sub =~ /^Bio::.*Submit$/) {
            next COUNTER;
        }
        if ($sub =~ /^Bio::.*Get_Caller$/) {
            next COUNTER;
        }
        ## All callers which are needing modules have the form Bio::Adventure::Something::Something
        my @pieces = split(/::/, $sub);
        if (scalar(@pieces) == 4) {
            $subroutine = $pieces[3];
            last COUNTER;
        }
    }
    return($subroutine);
}

=head2 C<Get_Menus>

This is responsible for printing out the cyoa menus and connecting the
various options to the functions in the package.

The top-level of the hash is the set of tasks, so stuff like sequence alignment,
mapping, indexing, conversion, counting, etc.

There are a few keys inside each of these; the names are the option which may be used
with --task to differentiate between contexts (e.g. if one wishes to use mapping
options specific for tnseq vs. snps, for example).  The message is just silly.
Finally, the choices are a set of strings which point to a function.
Choosing the number associated with the string will therefore lead to the
invocation of that function.

=cut
sub Get_Menus {
    my $menus = {
        Alignment => {
            name => 'alignment',
            message => 'Hari Seldon once said violence is the last refuge of the incompetent.  Go to page 6626070.',
            choices => {
                '(blastsplit): Split the input sequence into subsets and align with blast.' => \&Bio::Adventure::Align_Blast::Split_Align_Blast,
                '(fastasplit): Split the input sequence into subsets and align with fasta36.' => \&Bio::Adventure::Align_Fasta::Split_Align_Fasta,
                '(fastamismatch): Split the input sequence into subsets and align with fasta36.' => \&Bio::Adventure::Align_Fasta::Parse_Fasta_Mismatches,
                '(concat): Merge split searches into a single set of results.' => \&Bio::Adventure::Align::Concatenate_Searches,
                '(fastaparse): Parse fasta36 output into a reasonably simple table of hits.' => \&Bio::Adventure::Align_Fasta::Parse_Fasta,
                '(blastparse): Parse blast output into a reasonably simple table of hits.' => \&Bio::Adventure::Align_Blast::Parse_Blast,
                '(fastamerge): Merge and Parse fasta36 output into a reasonably simple table of hits.' => \&Bio::Adventure::Align_Fasta::Merge_Parse_Fasta,
                '(blastmerge): Merge and Parse blast output into a reasonably simple table of hits.' => \&Bio::Adventure::Align_Blast::Merge_Parse_Blast,
                '(orthomcl): Run OrthoMCL assuming a directory "input" with the amino acid faa files.' => \&Bio::Adventure::Align_Blast::OrthoMCL_Pipeline,
                '(orthofinder): Run OrthoFinder assuming a directory "input" with the amino acid faa files.' => \&Bio::Adventure::Align_Blast::OrthoFinder,
            },
        },
        Annotation => {
            name => 'annotation',
            message => 'How come Aquaman can control whales?  They are mammals!  Makes no sense.',
            choices =>  {
                '(abricate): Search for Resistance genes across databases.' => \&Bio::Adventure::Resistance::Abricate,
                '(casfinder): Search for Crispr/Cas9 cassettes/sequences in a bacterial assembly.' => \&Bio::Adventure::Annotation::Casfinder,
                '(classify_virus): Use ICTV data to classify viral sequences/contigs.' => \&Bio::Adventure::Phage::Classify_Phage,
                '(extend_kraken): Extend a kraken2 database with some new sequences.' => \&Bio::Adventure::Index::Extend_Kraken_DB,
                '(interproscan): Use interproscan to analyze ORFs.' => \&Bio::Adventure::Annotation::Interproscan,
                '(kraken2): Taxonomically classify reads.' => \&Bio::Adventure::Count::Kraken,
                '(phageterm): Invoke phageterm to hunt for likely phage ends.' => \&Bio::Adventure::Phage::Phageterm,
                '(phastaf): Invoke phastaf to attempt classifying phage sequence.' => \&Bio::Adventure::Phage::Phastaf,
                '(rhopredict): Search for rho terminators.' => \&Bio::Adventure::Feature_Prediction::Rho_Predict,
                '(terminasereorder): Reorder an assembly based on the results of a blast search.' => \&Bio::Adventure::Phage::Terminase_ORF_Reorder,
                '(prokka): Invoke prokka to annotate a genome.' => \&Bio::Adventure::Annotation::Prokka,
                '(resfinder): Search for antimicrobial resistance genes.' => \&Bio::Adventure::Resistance::Resfinder,
                '(rgi): Search for resistance genes with genecards and peptide fasta input.' => \&Bio::Adventure::Resistance::Rgi,
                '(mergeannotations): Merge annotations into a genbank file.' => \&Bio::Adventure::Metadata::Merge_Annotations_Make_Gbk,
            },
        },
        Assembly => {
            name => 'assembly',
            message => 'The wise man fears the wrath of a gentle heart. Go to page 314159.',
            choices => {
                '(abyss): Run abyss to create a new assembly.' => \&Bio::Adventure::Assembly::Abyss,
                '(assemblycoverage): Calculate Assembly coverage across contigs.' => \&Bio::Adventure::Assembly::Assembly_Coverage,
                '(extract_trinotate): Extract the most likely hits from Trinotate.' => \&Bio::Adventure::Annotation::Extract_Trinotate,
                '(filterdepth): Filter contigs based on sequencing depth.' => \&Bio::Adventure::Assembly::Unicycler_Filter_Depth,
                '(kraken2): Taxonomically classify reads.' => \&Bio::Adventure::Count::Kraken,
                '(transdecoder):  Run transdecoder on a putative transcriptome.' => \&Bio::Adventure::Assembly::Transdecoder,
                '(trinotate): Perform de novo transcriptome annotation with trinotate.' => \&Bio::Adventure::Annotation::Trinotate,
                '(trinity): Perform de novo transcriptome assembly with trinity.' => \&Bio::Adventure::Assembly::Trinity,
                '(trinitypost): Perform post assembly analyses with trinity.' => \&Bio::Adventure::Assembly::Trinity_Post,
                '(velvet): Perform de novo genome assembly with velvet.' => \&Bio::Adventure::Assembly::Velvet,
                '(unicycler): Perform de novo assembly with unicycler.' => \&Bio::Adventure::Assembly::Unicycler,
                '(shovill): Perform the shovill pre/post processing with spades.' => \&Bio::Adventure::Assembly::Shovill,
                '(terminasereorder): Reorder an existing assembly to the most likely terminase.' => \&Bio::Adventure::Phage::Terminase_ORF_Reorder,
                '(collect_assembly): Copy files generated by an assembly into one directory.' => \&Bio::Adventure::Metadata::Collect_Assembly,
                '(rosalind_plus): Make sure the Rosalind strand has the most pluses.' => \&Bio::Adventure::Annotation::Rosalind_Plus,
            },
        },
        Conversion => {
            name => 'convert',
            message => 'And it rained a fever. And it rained a silence. And it rained a sacrifice. And it rained a miracle. And it rained sorceries and saturnine eyes of the totem.  Go to page 2584981.',
            choices => {
                '(sam2bam): Convert a sam mapping to compressed/sorted/indexed bam.' => \&Bio::Adventure::Convert::Sam2Bam,
                '(bam2cov): Convert a sorted bam to tsv of coverage/base.' => \&Bio::Adventure::Convert::Bam2Coverage,
                '(gb2gff): Convert a genbank flat file to gff/fasta files.' => \&Bio::Adventure::Convert::Gb2Gff,
                '(gff2fasta): Convert a gff file to a fasta file.' => \&Bio::Adventure::Convert::Gff2Fasta,
                '(gatkdedup): Deduplicate an existing bam alignment with gatk.' => \&Bio::Adventure::Convert::GATK_Dedup,
            },
        },
        Counting => {
            name => 'count',
            message => 'Once men turned their thinking over to machines in the hope that this would set them free. But that only permitted other men with machines to enslave them.  Go to page 27812',
            choices => {
                '(guessstrand): Count reads with respect to features and report strandedness.' => \&Bio::Adventure::Count::Guess_Strand,
                '(htseq): Count mappings with htseq-count.' =>  \&Bio::Adventure::Count::HTSeq,
                '(htmulti): Use different option sets for counting with htseq.' => \&Bio::Adventure::Count::HT_Multi,
                '(insertsize): Gather insert size stats with gatk.' => \&Bio::Adventure::Count::Insert_Size,
                '(jellyfish): Perform a kmer count of some data.' => \&Bio::Adventure::Count::Jellyfish,
                '(mash): Use mash to count pairwise distances among sequences.' => \&Bio::Adventure::Count::Mash,
                '(mimap): Count mature miRNA species.' => \&Bio::Adventure::Count::Mi_Map,
                '(mpileup): Count coverage with mpileup.' => \&Bio::Adventure::Count::Mpileup,
                '(countstates): Count ribosome positions.' => \&Bio::Adventure::Riboseq::Count_States,
                '(slsearch): Count frequency of SL (or an arbitrary) sequences.' => \&Bio::Adventure::Splicing::SLSearch,
            },
        },
        FeaturePrediction => {
            name => 'features',
            message => '',
            choices => {
                '(aragorn): Search for tRNAs with aragorn.' => \&Bio::Adventure::Feature_Prediction::Aragorn,
                '(glimmer): Look for ORFs in bacterial/viral sequence.' => \&Bio::Adventure::Feature_Prediction::Glimmer,
                '(glimmersingle): Use glimmer to search for ORFs without training.' => \&Bio::Adventure::Feature_Prediction::Glimmer_Single,
                '(phanotate): Look for ORFs in bacterial/viral sequence.' => \&Bio::Adventure::Feature_Prediciton::Phanotate,
                '(prodigal): Look for ORFs in bacterial/viral sequence.' => \&Bio::Adventure::Feature_Prediction::Prodigal,
                '(rhopredict): Search for rho terminators.' => \&Bio::Adventure::Feature_Prediction::Rho_Predict,
                '(trnascan): Search for tRNA genes with trnascan.' => \&Bio::Adventure::Feature_Prediction::tRNAScan,
                '(trainprodigal): Train prodgial using sequences from a species/strain.' => \&Bio::Adventure::Feature_Prediction::Train_Prodigal,
            },
        },
        Mapping => {
            name => 'map',
            message => 'The world is dark and full of terrors, take this and go to page 6022140.',
            choices => {
                '(bowtie): Map trimmed reads with bowtie1 and count with htseq.' => \&Bio::Adventure::Map::Bowtie,
                '(bt2): Map trimmed reads with bowtie2 and count with htseq.' => \&Bio::Adventure::Map::Bowtie2,
                '(hisat): Map trimmed reads with hisat2 and count with htseq.' => \&Bio::Adventure::Map::Hisat2,
                '(btmulti): Map trimmed reads and count using multiple bowtie1 option sets.' => \&Bio::Adventure::Map::BT_Multi,
                '(bwa): Map reads with bwa and count with htseq.' => \&Bio::Adventure::Map::BWA,
                '(kallisto): Pseudo-align and count reads using kallisto.' => \&Bio::Adventure::Map::Kallisto,
                '(salmon): Pseudo-align and count reads using salmon.' => \&Bio::Adventure::Map::Salmon,
                '(star): Pseudo-align and count reads using STAR.' => \&Bio::Adventure::Map::STAR,
                '(mimap): Attempt to map reads explicitly to mature miRNA species.' => \&Bio::Adventure::Count::Mi_Map,
                '(rrnabowtie): Map rRNA reads using bowtie1.' => \&Bio::Adventure::Map::Bowtie_RRNA,
                '(rsem): Quantify reads using rsem.' => \&Bio::Adventure::Map::RSEM,
                '(tophat): Map reads using tophat2 and count with htseq.' => \&Bio::Adventure::Map::Tophat,
            },
        },
        Indexers => {
            name => 'index',
            message => '',
            choices => {
                '(extend_kraken): Extend a kraken2 database with some new sequences.' => \&Bio::Adventure::Index::Extend_Kraken_DB,
                '(indexbt1): Create bowtie1 compatible indexes.' => \&Bio::Adventure::Index::BT1_Index,
                '(indexbt2): Create bowtie2 compatible indexes.' => \&Bio::Adventure::Index::BT2_Index,
                '(indexhisat): Create hisat2 compatible indexes.' => \&Bio::Adventure::Index::Hisat2_Index,
                '(indexbwa): Create bwa compatible indexes.' => \&Bio::Adventure::Index::BWA_Index,
                '(indexkallisto): Create kallisto compatible indexes.' => \&Bio::Adventure::Index::Kallisto_Index,
                '(indexrsem): Create rsem indexes.' => \&Bio::Adventure::Index::RSEM_Index,
                '(indexsalmon): Create salmon indexes.' => \&Bio::Adventure::Index::Salmon_Index,
            },
        },
        Phage => {
            name => 'phage',
            message => '',
            choices => {
                '(filterkraken): Filter out host sequences using a kraken report.' => \&Bio::Adventure::Phage::Filter_Host_Kraken,
                '(classifyphage): Use ICTV data to classify a phage assembly.' => \&Bio::Adventure::Phage::Classify_Phage,
                '(phageterm): Invoke phageterm to search for terminal repeats.' => \&Bio::Adventure::Phage::Phageterm,
                '(phagepromoter): Search for phage promoters.' => \&Bio::Adventure::Feature_Prediction::Phagepromoter,
                '(phastaf): Search for phage regions in arbitrary assemblies.' => \&Bio::Adventure::Phage::Phastaf,
                '(restriction): Search for restriction sites.' => \&Bio::Adventure::Phage::Restriction_Catalog,
                '(terminasereorder): Reorder a phage assembly to the putative terminase.' => \&Bio::Adventure::Phage::Terminase_ORF_Reorder,
                '(caical): Calculate codon adaptation index of a phage vs. some host.' => \&Bio::Adventure::Phage::Caical,
                '(recatalog): Create a catalog of restriction enzyme hits.' => \&Bio::Adventure::Phage::Restriction_Catalog,
                '(xref_crispr): Cross reference against observed crispr sequences.' => \&Bio::Adventure::Phage::Xref_Crispr,
            },
        },
        Phylogeny => {
            name => 'phylogeny',
            message => '',
            choices => {
                '(gubbins): Run Gubbins with an input msa.' => \&Bio::Adventure::Phylogeny::Run_Gubbins,
            },
        },
        Pipeline => {
            name => 'pipeline',
            message => 'When Mr. Bilbo Baggins announced he would shortly be celebrating his eleventyfirst birthday, there was much talk and excitement in Hobbiton.  Go to page 1618033',
            choices => {
                '(pdnaseq): Preset DNAseq tasks.' => \&Bio::Adventure::Pipeline::Process_DNAseq,
                '(prnaseq): Preset RNAseq tasks.' => \&Bio::Adventure::Pipeline::Process_RNAseq,
                '(priboseq): Perform a preset pipeline of ribosome profiling tasks.' => \&Bio::Adventure::Pipeline::Riboseq,
                '(ptnseq): Perform a preset pipeline of TNSeq tasks.' => \&Bio::Adventure::Pipeline::TNSeq,
                '(pbt1): Perform a preset group of bowtie1 tasks.' => \&Bio::Adventure::Pipeline::Bowtie,
                '(pbt2): Use preset bowtie2 tasks.' => \&Bio::Adventure::Pipeline::Bowtie2,
                '(phisat): Use preset tophat tasks.' => \&Bio::Adventure::Pipeline::Hisat,
                '(pbwa): Try preset bwa tasks.' => \&Bio::Adventure::Pipeline::BWA,
                '(psalmon): Try preset salmon tasks.' => \&Bio::Adventure::Pipeline::Salmon,
                '(passemble): Preset assembly.' => \&Bio::Adventure::Pipeline::Assemble,
                '(phageassemble): Preset phage assembly.' => \&Bio::Adventure::Pipeline::Phage_Assemble,
                '(annotateassembly): Do some searches on an assembly.' => \&Bio::Adventure::Pipline::Annotation_Assembly,
            },
        },
        PopulationGenetics => {
            name => 'popgen',
            message => 'The front pattern DOES move - and no wonder! The woman behind shakes it!',
            choices => {
                '(angsd): Get the set of filtered variants given a list of bam files.' => \&Bio::Adventure::PopGen::Angsd_Filter,
            },
        },
        Prepare => {
            name => 'preparation',
            message => 'Whan that Aprille withe her shoures sote, the droughte of Marche hath perced to the rote.  go to Cantebury.',
            choices => {
                '(read_samples): Read samples using a csv file to determine the raw data locations.' => \&Bio::Adventure::Prepare::Read_Samples,
                '(copyraw): Copy data from the raw data archive to scratch.' => \&Bio::Adventure::Prepare::Copy_Raw,
                '(sradownload): Extract SRA accession from a bioproject and download.' => \&Bio::Adventure::Prepare::Download_SRA_PRJNA,
                '(fastqdump): Download data from sra.' => \&Bio::Adventure::Prepare::Fastq_Dump,
                '(ncbidownload): Download accessions from ncbi.' => \&Bio::Adventure::Prepare::Download_NCBI_Accessions,
            },
        },
        QA => {
            name => 'qa',
            message => 'There is a time when the operation of the machine becomes so odious that you cannot take part.',
            choices => {
                '(biopieces): Use biopieces to graph some metrics of the data.' => \&Bio::Adventure::QA::Biopieces_Graph,
                '(cogent): Use cogent to remove UMIs and trim data.' => \&Bio::Adventure::Trim::Cogent,
                '(cutadapt): Perform adapter trimming with cutadapt.' => \&Bio::Adventure::Trim::Cutadapt,
                '(fastp): Use fastp to trim fastq/remove UMIs.' => \&Bio::Adventure::Trim::Fastp,
                '(fastqc): Use fastqc to check the overall quality of the raw data.' => \&Bio::Adventure::QA::Fastqc,
                '(racer): Perform sequence correction with hitec/RACER.' => \&Bio::Adventure::Trim::Racer,
                '(trimomatic): Perform adapter trimming with Trimomatic.' => \&Bio::Adventure::Trim::Trimomatic,
                '(umi_tools): Perform umi extraction with umi_tools.' => \&Bio::Adventure::Trim::Umi_Tools_Extract,
                '(umi_dedup): Perform umi deduplication with umi_tools following mapping.' => \&Bio::Adventure::Trim::Umi_Tools_Dedup,
            },
        },
        RiboSeq => {
            name => 'riboseq',
            message => 'Awake Awake Fear Fire Foes!  Go to page 5291772',
            choices => {
                '(biopieces): Make some plots of the demultiplexed/trimmed reads.' => \&Bio::Adventure::QA::Biopieces_Graph,
                '(cutadapt): Use cutadapt to remove the adapters.' => \&Bio::Adventure::Trim::Cutadapt,
                '(rrnabowtie): Map rRNA reads using bowtie1.' => \&Bio::Adventure::Map::Bowtie_RRNA,
                '(btmulti): Use bowtie1 to find putative ribosomal positions.' => \&Bio::Adventure::Map::BT_Multi,
                '(calibrate): Calibrate the positions of the a/p/e sites of the ribosomes.' => \&Bio::Adventure::Riboseq::Calibrate,
                '(countstates): Count the positions of a/p/e/etc sites across the genome/transcriptome.' => \&Bio::Adventure::Riboseq::Count_States,
                '(graphreads): Plot the coverage of ribosomes across the genome by ORF.' => \&Bio::Adventure::Riboseq::Graph_Reads,
            },
        },
        SNP => {
            name => 'snp',
            message => qq"When my god comes back I'll be waiting for him with a shotgun.  And I'm keeping the last shell for myself. (inexact quote)  Go to page 667408",
            choices => {
                '(mpileup): Count coverage with mpileup.' => \&Bio::Adventure::Count::Mpileup,
                '(trim): Trim sequence with an additional rule to remove the first 10 nucleotides.' => \&Bio::Adventure::Trim::Trimomatic,
                '(bwa): Map reads with bwa and count with htseq.' => \&Bio::Adventure::Map::BWA,
                '(bowtie): Map trimmed reads with bowtie1 and count with htseq.' => \&Bio::Adventure::Map::Bowtie,
                '(bt2): Map trimmed reads with bowtie2 and count with htseq.' => \&Bio::Adventure::Map::Bowtie2,
                '(freebayes): Use freebayes to create a vcf file and filter it.' => \&Bio::Adventure::SNP::Freebayes_SNP_Search,
                '(hisat): Map trimmed reads with hisat2 and count with htseq.' => \&Bio::Adventure::Map::Hisat2,
                '(parsnp): Parse an existing bcf file and print some fun tables.' => \&Bio::Adventure::SNP::SNP_Ratio,
                '(snpsearch): Use mpileup to create a vcf file and filter it. (bam input)' => \&Bio::Adventure::SNP::Mpileup_SNP_Search,
                '(snpratio): Count the variant positions by position and create a new genome. (bcf input)' => \&Bio::Adventure::SNP::SNP_Ratio,
                '(snp): Perform alignments and search for variants. (fastq input)' => \&Bio::Adventure::SNP::Align_SNP_Search,
                '(snippy): Invoke snippy. (fastq and genbank inputs)' => \&Bio::Adventure::SNP::Snippy,
            },
        },
        TNSeq => {
            name => 'tnseq',
            message => 'You have enterred a world of jumping DNA, be ware and go to page 42.',
            choices => {
                '(sortindex): Demultiplex raw reads based on the peculiar TNSeq indexes.' => \&Bio::Adventure::TNSeq::Sort_Indexes,
                '(consolidate): Extract the reads of read-pairs with beginning TAs.' => \&Bio::Adventure::TNSeq::Consolidate_TAs,
                '(cutadapt): Use cutadapt to remove the odd tnseq adapters.' => \&Bio::Adventure::Trim::Cutadapt,
                '(tacheck): Make certain that all reads have a leading or terminal TA.' => \&Bio::Adventure::TNSeq::TA_Check,
                '(biopieces): Make some plots of the demultiplexed/trimmed reads.' => \&Bio::Adventure::QA::Biopieces_Graph,
                '(essentialityta): Count the hits/TA in preparation for essentiality.' => \&Bio::Adventure::TNSeq::Essentiality_TAs,
                '(runessentiality): Run the essentiality suite of tools.' => \&Bio::Adventure::TNSeq::Run_Essentiality,
                '(gumbel): Run the essentiality suite of tools on the ta counts.' => \&Bio::Adventure::TNSeq::Run_Essentiality,
                '(tpp): Run the transit preprocessing script.' => \&Bio::Adventure::TNSeq::Transit_TPP,
            },
        },
        Splicing => {
            name => 'splicing',
            message => 'Tools to consider splicing',
            choices => {
                '(rmats): Quantify changes in transcript splicing events across samples.' => \&Bio::Adventure::Splicing::RMats,
                '(suppa): Quantify changes in transcript splicing events across samples.' => \&Bio::Adventure::Splicing::Suppa,
                '(slseq_recorder): Record observed spliced leaders in an SLSeq experiment.' => \&Bio::Adventure::Splicing_SLSeq_Recorder,
                '(slsearch): Count frequency of SL (or an arbitrary) sequences.' => \&Bio::Adventure::Splicing::SLSearch,
                '(slutr): Search for SL junction reads and use them to define UTRs.' => \&Bio::Adventure::Splicing::SL_UTR,
            },
        },
        Structure => {
            name => 'structure',
            message => 'Intended for gene and RNA structure tools',
            choices => {
                '(vienna): Use the RNAfold method from vienna to examine potential RNA structures.' => \&Bio::Adventure::Structure::Vienna,
            },
        },
        Visualize => {
            name => 'visualize',
            message => 'When Red wins, she stands alone.  Go to page: 96485.332',
            choices =>  {
                '(cgview): Invoke cgview to visualize a genome.' => \&Bio::Adventure::Visualization::CGView,
            },
        },
        Test => {
            name => 'test',
            message => 'All happy families are happy in the same way. Go to page 5670367.',
            choices => {
                '(test): Run a test job' => \&Test_Job,
            },
        },
    };
    return($menus);
}

sub Estimate_Memory_Single {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars();
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
        my $size = $st->size;
        my @pieces = split(/\./, $input);
        my $ext = $pieces[$#pieces];
        if (defined($multipliers->{$ext})) {
            $size = $size * $job_multiplier;
        }
        $total_mem = $total_mem + $size;
    }
    return($total_mem);
}

sub Estimate_Time_Single {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        process => 'xz',
        species => '',
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
        my $st = stat($input);
        my $size = $st->size;
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

sub Estimate_Time_Composite {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        process => 'xz',
    );

    my $species_file = qq"$options->{libpath}/$options->{libtype}/indexes/$options->{species}";
    my $species_stat = stat($species_file);
    my $species_size = $species_stat->size;

    ## I need to acquire from slurm the relative cpu frequency, but assuming similar speeds
    ## I should be able to get reasonable time estimates...
    my @inputs = split(/[:;,\s+]/, $options->{input});
    ## I am going to assume that the less preprocessor will make xz/gz/bz/uncomp
    ## input file size differences negligible with respect to time to run.
    my $total_input_size = 0;
    for my $input (@inputs) {
        my $st = stat($input);
        my $size = $st->size;
        $total_input_size = $total_input_size + $size;
    }

    # A genome multiplier should be used when the time estimate for a tool can be described as:
    ## (input_size * constant) * (genome_size * constant)
    ## Really it should be the index size, but I assume the fasta is a sufficient proxy.
    my $genome_mult = {
        hisat => 1.0,
    };
    ## Approximate compression multipliers, taken from:
    ## https://community.centminmod.com/threads/round-3-compression-comparison-benchmarks-zstd-vs-brotli-vs-pigz-vs-bzip2-vs-xz-etc.17259/
    my $multipliers = {
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
    my $key = $options->{process};
    my $needed_input_mult = $total_input_size * $multipliers->{$key};
    my $needed_genome_mult = $species_size * $genome_mult->{$key};
    my $needed = ($needed_input_mult * $needed_genome_mult);
    my $total_hours = ceil($needed * 60.0 * 60.0);
    return($total_hours);
}

sub Get_Modules {
    my %args = @_;
    my $subroutine = Bio::Adventure::Config::Get_Caller();
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
        'Aragorn' => { modules => 'aragorn', exe => 'aragorn' },
        'Assembly_Coverage' => {
            modules => ['cyoa', 'hisat2', 'samtools', 'bbmap',], },
        'Bacphlip' => { modules => 'bacphlip', exe => 'bacphlip' },
        'Bam2Coverage' => { modules => ['samtools', 'bbmap',] },
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
        'Downsample_Guess_Strand' => { modules => ['salmon'] },
        'Essentiality_TAs' => { modules => 'cyoa', },
        'Extend_Kraken_DB' => { modules => ['kraken'], exe => ['kraken2'] },
        'Fasta2Gff' => { modules => ['cyoa'] },
        'Fastqc' => { modules => ['fastqc'] },
        'Fastq_Dump' => { modules => ['cyoa', 'sra'], },
        'Fastp' => { modules => ['fastp'], exe => 'fastp' },
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
        'Glimmer' => { modules => ['cyoa', 'glimmer'], exe => 'glimmer3' },
        'Glimmer_Single' => {
            modules => ['glimmer', 'cyoa'], exe => 'glimmer3' },
        'Gubbins' => { modules => 'gubbins' },
        'Guess_Strand' => { modules => ['cyoa', 'samtools'] },
        'Hisat2' => { modules => 'hisat2', exe => 'hisat2'},
        'Hisat2_Index' => { modules => ['hisat2'], exe => ['hisat2'], },
        'HT_Multi' => { modules => 'htseq' },
        'HTSeq' => { modules => 'htseq', exe => 'htseq-count' },
        'Insert_Size' => { modules => 'gatk' },
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
        'Make_Codon_Table' => { modules => ['cyoa'] },
        'Mash' => { modules => 'mash', exe => 'mash' },
        'Merge_Annotations' => { modules => ['cyoa', 'ncbi_tools/6.1'], exe => 'tbl2asn' },
        'Merge_CDS_Predictions' => { modules => ['cyoa', 'ncbi_tools/6.1'], exe => 'tbl2asn' },
        'Merge_Parse_Blast' => { modules => ['cyoa'], },
        'Mpileup' => { modules => 'samtools' },
        'Mpileup_SNP_Search' => {
            modules => ['samtools/1.13', 'bcftools', 'vcftools'],
            exe => ['samtools', 'bcftools'], },
        'OrthoFinder' => { modules => ['cyoa', 'orthofinder'], exe => ['orthofinder'] },
        'OrthoMCL_Pipeline' => { modules => ['orthomcl'], exe => ['orthomclPairs'] },
        'Parse_Fasta_Global' => { modules => 'cyoa' },
        'Phagepromoter' => { modules => 'phagepromoter' },
        'Phageterm' => { modules => ['cyoa', 'phageterm'], exe => 'PhageTerm.py' },
        'Phanotate' => { modules => ['trnascan', 'phanotate'], exe => 'phanotate.py' },
        'Phastaf' => { modules => ['cyoa', 'phastaf'] },
        'PolyA_Extractor_Worker' => { modules => ['cyoa'] },
        'Prodigal' => { modules => ['cyoa', 'prodigal'] },
        'Prokka' => { ## Prokka should not need cyoa; it is getting requisite perl module from it for now
            modules => ['cyoa', 'prokka', 'blast'], exe => 'prokka'},
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
        'SLSeq_Recorder' => { modules => ['cyoa', 'samtools'] },
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
        'STAR' => { modules => 'star' },
        'STAR_Index' => { modules => ['star'] },
        'Suppa' => { modules => ['suppa'] },
        'TA_Check' => { modules => 'cyoa' },
        'Terminase_ORF_Reorder' => { modules => ['cyoa', 'fasta', 'blast', 'blastdb'], },
        'Train_Prodigal' => { modules => 'prodigal' },
        'Tophat' => { modules => 'tophat' },
        'Transdecoder' => { modules => 'transdecoder' },
        'Transit_TPP' => { modules => ['bwa', 'transit', 'htseq' ], exe => 'tpp' },
        'Transposonpsi' => { modules => 'transposonpsi', exe => 'transposonPSI.pl' },
        ## I can set the exe for single vs pairwise now
        'Trimomatic_Pairwise' => { modules => ['trimomatic'], },
        'Trimomatic_Single' => { modules => ['trimomatic'], },
        'Trinity' => { modules => 'trinity', exe => 'Trinity', },
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
    my @function_array = split(/::/, $subroutine);
    my $function = $function_array[$#function_array];
    $function =~ s/_Worker$//g;
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
            }
            else {
                $datum = {};
            }

            my $conda_name = $datum->{conda};
            if ($conda_name) {
                push(@conda_lst, $conda_name);
                $datum->{conda} = \@conda_lst;
            }
        }
    }
    else {
        print "The function: $function does not appear to have defined modules.\n" if ($args{debug});
        $datum = {};
    }
    return(%{$datum});
}

sub Get_Paths {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $subroutine = $args{subroutine};
    $subroutine = Bio::Adventure::Config::Get_Caller() unless($subroutine);
    my $libpath_prefix = qq"$options->{libpath}/$options->{libtype}";
    my $libdir_prefix = qq"$options->{libdir}/$options->{libtype}";
    my $output_prefix = qq"$options->{output_base}/$options->{jprefix}";
    my $paths = {};
    if (defined($options->{species})) {
        $paths->{fasta} = qq"${libpath_prefix}/fasta/$options->{species}.fasta";
        $paths->{fasta_shell} = qq"${libdir_prefix}/fasta/$options->{species}.fasta";
        $paths->{gff} = qq"${libpath_prefix}/gff/$options->{species}.gff";
        $paths->{gff_shell} = qq"${libdir_prefix}/gff/$options->{species}.gff";
        $paths->{gtf} = qq"${libpath_prefix}/gtf/$options->{species}.gtf";
        $paths->{gtf_shell} = qq"${libdir_prefix}/gtf/$options->{species}.gtf";
        $paths->{gbk} = qq"${libpath_prefix}/genbank/$options->{species}.gbk";
        $paths->{gbk_shell} = qq"${libdir_prefix}/genbank/$options->{species}.gbk";
        my $fasta_check = dirname($paths->{fasta});
        make_path($fasta_check) unless (-d $fasta_check);
        my $gff_check = dirname($paths->{gff});
        make_path($gff_check) unless (-d $gff_check);
        my $gtf_check = dirname($paths->{gtf});
        make_path($gtf_check) unless (-d $gtf_check);
        my $gbk_check = dirname($paths->{gbk});
        make_path($gbk_check) unless (-d $gbk_check);
    }
    my $index_prefix = qq"${libpath_prefix}/indexes";
    my $index_prefix_shell = qq"${libdir_prefix}/indexes";
    if ($subroutine eq 'Bowtie') {
        $paths->{index} = qq"${index_prefix}/bt1/$options->{species}";
        $paths->{index_shell} = qq"${index_prefix_shell}/bt1/$options->{species}";
        $paths->{index_file} = qq"$paths->{index}.1.ebwt";
        $paths->{index_file_shell} = qq"$paths->{index_shell}.1.ebwt";
        $paths->{output_dir} = qq"${output_prefix}bowtie_$options->{species}";
    }
    elsif ($subroutine eq 'Bowtie2') {
        $paths->{index} = qq"${index_prefix}/bt2/$options->{species}";
        $paths->{index_shell} = qq"${index_prefix_shell}/bt2/$options->{species}";
        $paths->{index_file} = qq"$paths->{index}.1.bt2";
        $paths->{index_file_shell} = qq"$paths->{index_shell}.1.bt2";
        $paths->{index_file2} = qq"$paths->{index}.1.bt2l";
        $paths->{index_file2_shell} = qq"$paths->{index_shell}.1.bt2l";
        $paths->{output_dir} = qq"${output_prefix}bowtie2_$options->{species}";
    }
    elsif ($subroutine eq 'BWA') {
        $paths->{index} = qq"${index_prefix}/bwa/$options->{species}";
        $paths->{index_shell} = qq"${index_prefix_shell}/bwa/$options->{species}";
        $paths->{index_file} = qq"$paths->{index}.sa";
        $paths->{index_file_shell} = qq"$paths->{index_shell}.sa";
        $paths->{output_dir} = qq"${output_prefix}bwa_$options->{species}";
    }
    elsif ($subroutine eq 'Downsample_Guess_Strand') {
        $paths->{index} = qq"${index_prefix}/salmon/$options->{species}";
        $paths->{index_shell} = qq"${index_prefix_shell}/salmon/$options->{species}";
        $paths->{output_dir} = qq"${output_prefix}salmon_downsample_$options->{species}";
    }
    elsif ($subroutine eq 'GATK_Dedup') {
        $paths->{output_dir} = qq"${output_prefix}gatk_dedup";
    }
    elsif ($subroutine eq 'Hisat2' or $subroutine eq 'Hisat2_Index') {
        $paths->{index_dir} = qq"${index_prefix}/hisat";
        $paths->{index} = qq"$paths->{index_dir}/$options->{species}";
        $paths->{index_shell} = qq"${index_prefix_shell}/hisat/$options->{species}";
        $paths->{index_file} = qq"$paths->{index}.1.ht2";
        $paths->{index_file_shell} = qq"$paths->{index_shell}.1.ht2";
        $paths->{index_file2} = qq"$paths->{index}.1.ht2l";
        $paths->{index_file2_shell} = qq"$paths->{index_shell}.1.ht2l";
        $paths->{output_dir} = qq"${output_prefix}hisat_$options->{species}";
    }
    elsif ($subroutine eq 'Kallisto') {
        $paths->{index_file} = qq"${index_prefix}/kallisto/$options->{species}.idx";
        $paths->{index_file_shell} = qq"${index_prefix_shell}/kallisto/$options->{species}.idx";
        $paths->{output_dir} = qq"${output_prefix}kallisto_$options->{species}";
    }
    elsif ($subroutine eq 'OrthoFinder') {
        $paths->{output_dir} = qq"${output_prefix}orthofinder";
    }
    elsif ($subroutine eq 'Pairwise_Similarity_Matrix') {
        $paths->{output_dir} = qq"${output_prefix}pairwise_blast";
    }
    elsif ($subroutine eq 'ProgressiveMauve') {
        $paths->{output_dir} = qq"${output_prefix}/pmauve";
    }
    elsif ($subroutine eq 'RSEM') {
        $paths->{index} = qq"${index_prefix}/rsem/$options->{species}";
        $paths->{index_shell} = qq"${index_prefix_shell}/rsem/$options->{species}";
        $paths->{index_file} = qq"$paths->{index}.transcripts.fa";
        $paths->{index_file_shell} = qq"$paths->{index_shell}.transcripts.fa";
        $paths->{output_dir} = qq"${output_prefix}rsem_$options->{species}";
    }
    elsif ($subroutine eq 'Run_Parse_Blast') {
        my $libname = basename($options->{library}, $class->{suffixes});
        $paths->{output_dir} = qq"${output_prefix}blast_$options->{input}_vs_${libname}";
    }
    elsif ($subroutine eq 'Salmon') {
        $paths->{index} = qq"${index_prefix}/salmon/$options->{species}";
        $paths->{index_shell} = qq"${index_prefix_shell}/salmon/$options->{species}";
        $paths->{output_dir} = qq"${output_prefix}salmon_$options->{species}";
    }
    elsif ($subroutine eq 'SLSearch') {
        $paths->{output_dir} = qq"${output_prefix}slsearch_$options->{species}";
    }
    elsif ($subroutine eq 'SLSeq_Recorder') {
        $paths->{output_dir} = qq"${output_prefix}slseq_recorder_$options->{species}";
    }
    elsif ($subroutine eq 'Split_Align_Blast') {
        my $libname = basename($options->{library}, $class->{suffixes});
        $paths->{output_dir} = qq"${output_prefix}blastsplit_$options->{input}_vs_${libname}";
    }
    elsif ($subroutine eq 'STAR') {
        $paths->{index} = qq"${index_prefix}/star/$options->{species}";
        $paths->{index_shell} = qq"${index_prefix_shell}/star/$options->{species}";
        $paths->{index_file} = qq"$paths->{index}/SAindex";
        $paths->{index_file_shell} = qq"$paths->{index_shell}/SAindex";
        $paths->{output_dir} = qq"${output_prefix}star_$options->{species}";
    }
    elsif ($subroutine eq 'Tophat') {
        $paths->{index} = qq"${index_prefix}/tophat/$options->{species}";
        $paths->{index_shell} = qq"${index_prefix_shell}/tophat/$options->{species}";
        $paths->{index_file} = qq"$paths->{index}.1.bt2";
        $paths->{index_file_shell} = qq"$paths->{index_shell}.1.bt2";
        $paths->{output_dir} = qq"${output_prefix}tophat_$options->{species}";
    }
    elsif ($subroutine eq 'Umi_Tools_Extract') {
        $paths->{output_dir} = qq"${output_prefix}umi_tools";
    }
    elsif ($subroutine eq 'Umi_Tools_Dedup') {
        $paths->{output_dir} = qq"${output_prefix}umi_dedup";
    }
    if ($paths->{index_file}) {
        my $index_directory = dirname($paths->{index_file});
        make_path($index_directory) unless (-d $index_directory);
    }
    if ($paths->{output_dir}) {
        make_path($paths->{output_dir});
    }
    return($paths);
}

=head2 C<Get_TODOs>

The keys of possible_todos provide the GetOpt::Long options required
to figure out if the user is requesting the various functions.  When
using GetOpt::Long, if an options has a '+' suffix, its value is
incremented by one when it is set without a following argument.  Thus
doing --abricate will set the
todo_list{todo}{Bio::Adventure::Resistance::Abricate} value from 0 to 1.
The cyoa script will then iterate over the todo hash and look for
anything with a positive value and run its associated function.

=cut
sub Get_TODOs {
    my %args = @_;
    my $todo_list = ();
    my $possible_todos = {
        ## Run Torsten Seeman's abricate against an assembly to look for groups of genes,
        ## primarily antibiotic resistance, but any database may be performed.
        "abricate+" => \$todo_list->{todo}{'Bio::Adventure::Resistance::Abricate'},
        ## Create an assembly using Abyss.
        "abyss+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Abyss'},
        ## Filter a population genetics dataset via Angsd
        "angsdfilter+" => \$todo_list->{todo}{'Bio::Adventure::PopGen::Angsd_Filter'},
        ## Use my phage pipeline to attempt to annotate a phage assembly.
        "annotatephage+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Annotate_Phage'},
        ## Given an assembly, search for tRNA-like genes.
        "aragorn+" => \$todo_list->{todo}{'Bio::Adventure::Feature_Prediction::Aragorn'},
        ## Remap reads to a new assembly to estimate coverage on a per-baseis.
        "assemblycoverage+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Assembly_Coverage'},
        ## Use the logic from Assembly_Coverage to calculate coverage vs. an existing genome.
        "bam2cov+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Bam2Coverage'},
        ## Use bedtools to calculate coverage, partially redundant with bam2cov
        "bedcov+" => \$todo_list->{todo}{'Bio::Adventure::Count::Bedtools_Coverage'},
        ## Biopieces provides some fun utilities for QA/QC.
        "biopieces+" => \$todo_list->{todo}{'Bio::Adventure::QA::Biopieces_Graph'},
        ## Split up a large set of sequences, run blast on them, and merge the result.
        "blastmerge+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Merge_Blast_Parse'},
        ## Simplify the results from blast into a table
        "blastparse+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Blast_Parse'},
        ## Perform all of the steps noted in blastmerge
        "blastsplitalign+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Split_Align_Blast'},
        ## Run old-school bowtie
        "bowtie+" => \$todo_list->{todo}{'Bio::Adventure::Map::Bowtie'},
        "bowtie2+" => \$todo_list->{todo}{'Bio::Adventure::Map::Bowtie2'},
        "bowtierrna+" => \$todo_list->{todo}{'Bio::Adventure::Map::Bowtie_RRNA'},
        "rrnabowtie+" => \$todo_list->{todo}{'Bio::Adventure::Map::Bowtie_RRNA'},
        "bt2+" => \$todo_list->{todo}{'Bio::Adventure::Map::Bowtie2'},
        "btmulti+" => \$todo_list->{todo}{'Bio::Adventure::Map::BT_Multi'},
        "bwa+" => \$todo_list->{todo}{'Bio::Adventure::Map::BWA'},
        "caical+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Caical'},
        "calibrate+" => \$todo_list->{todo}{'Bio::Adventure::Riboseq::Calibrate'},
        "casfinder+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Casfinder'},
        "cgview+" => \$todo_list->{todo}{'Bio::Adventure::Visualization::CGView'},
        "checkslurm+" => \$todo_list->{todo}{'Bio::Adventure::Check_Job'},
        "checkjob+" => \$todo_list->{todo}{'Bio::Adventure::Check_Job'},
        "classifyphage+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Classify_Phage'},
        "cogent+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Cogent'},
        "collectassembly+" => \$todo_list->{todo}{'Bio::Adventure::Metadata::Collect_Assembly'},
        "consolidate+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Consolidate_TAs'},
        "copyraw+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Copy_Raw'},
        "countstates+" => \$todo_list->{todo}{'Bio::Adventure::Riboseq::Count_States'},
        "concat+" => \$todo_list->{todo}{'Bio::Adventure::Align::Concatenate_Searches'},
        "cutadapt+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Cutadapt'},
        "dedupgatk+" => \$todo_list->{todo}{'Bio::Adventure::Convert::GATK_Dedup'},
        "ncbidownload+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Download_NCBI_Accession'},
        "downsampleguess+" => \$todo_list->{todo}{'Bio::Adventure::Map::Downsample_Guess_Strand'},
        "essentialitytas+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Essentiality_TAs'},
        "extendkraken+" => \$todo_list->{todo}{'Bio::Adventure::Index::Extend_Kraken_DB'},
        "extracttrinotate+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Extract_Trinotate'},
        "fastp+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Fastp'},
        "splitalignfasta+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Split_Align_Fasta'},
        "fastado+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Do_Fasta'},
        "fastasplitalign+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Split_Align_Fasta'},
        "fastamerge+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Merge_Parse_Fasta'},
        "fastamismatch+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Parse_Fasta_Mismatches'},
        "fastaparse+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Parse_Fasta'},
        "fasta2gff+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Fasta2Gff'},
        "fastqct+" => \$todo_list->{todo}{'Bio::Adventure::QA::Fastqc'},
        "fastqdump+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Fastq_Dump'},
        "featureextract+" => \$todo_list->{todo}{'Bio::Adventure::Annotation_Genbank::Extract_Features'},
        "filterdepth+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Unicycler_Filter_Depth'},
        "filterkraken+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Filter_Host_Kraken'},
        "fragment+" => \$todo_list->{todo}{'Bio::Adventure::SeqMisc::Fragment_Genome'},
        "freebayes+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Freebayes_SNP_Search'},
        "gatkdedup+" => \$todo_list->{todo}{'Bio::Adventure::Convert::GATK_Dedup'},
        "gb2gff+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Gb2Gff'},
        "gff2fasta+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Gff2Fasta'},
        "glimmer+" => \$todo_list->{todo}{'Bio::Adventure::Feature_Prediction::Glimmer'},
        "glimmersingle+" => \$todo_list->{todo}{'Bio::Adventure::Feature_Prediction::Glimmer_Single'},
        "graphreads+" => \$todo_list->{todo}{'Bio::Adventure::Riboseq::Graph_Reads'},
        "gubbins+" => \$todo_list->{todo}{'Bio::Adventure::Phylogeny::Run_Gubbins'},
        "guessstrand+" => \$todo_list->{todo}{'Bio::Adventure::Count::Guess_Strand'},
        "gumbel+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Run_Essentiality'},
        "hisat+" => \$todo_list->{todo}{'Bio::Adventure::Map::Hisat2'},
        "htmulti+" => \$todo_list->{todo}{'Bio::Adventure::Count::HT_Multi'},
        "htseq+" => \$todo_list->{todo}{'Bio::Adventure::Count::HTSeq'},
        "indexhisat+" => \$todo_list->{todo}{'Bio::Adventure::Index::Hisat2_Index'},
        "indexbt1+" => \$todo_list->{todo}{'Bio::Adventure::Index::BT1_Index'},
        "indexbt2+" => \$todo_list->{todo}{'Bio::Adventure::Index::BT2_Index'},
        "indexbwa+" => \$todo_list->{todo}{'Bio::Adventure::Index::BWA_Index'},
        "indexkallisto+" => \$todo_list->{todo}{'Bio::Adventure::Index::Kallisto_Index'},
        "indexrsem+" => \$todo_list->{todo}{'Bio::Adventure::Index::RSEM_Index'},
        "indexsalmon+" => \$todo_list->{todo}{'Bio::Adventure::Index::Salmon_Index'},
        "insertsize+" => \$todo_list->{todo}{'Bio::Adventure::Count::Insert_Size'},
        "interproscan+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Interproscan'},
        "interprolong+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Interpro_Long2Wide'},
        "jellyfish+" => \$todo_list->{todo}{'Bio::Adventure::Count::Jellyfish'},
        "kallisto+" => \$todo_list->{todo}{'Bio::Adventure::Map::Kallisto'},
        "kraken+" => \$todo_list->{todo}{'Bio::Adventure::Count::Kraken'},
        "mash+" => \$todo_list->{todo}{'Bio::Adventure::Count::Mash'},
        "mergeannotations+" => \$todo_list->{todo}{'Bio::Adventure::Metadata::Merge_Annotations'},
        "mergecds+" => \$todo_list->{todo}{'Bio::Adventure::Annotation_Genbank::Merge_CDS_Predictions'},
        "mergeparse+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Merge_Parse_Blast'},
        "mergeprodigal+" => \$todo_list->{todo}{'Bio::Adventure::Metadata::Merge_Annot_Prodigal'},
        "mimap+" => \$todo_list->{todo}{'Bio::Adventure::MiRNA::Mi_Map'},
        "mpileup+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Mpileup_SNP_Search'},
        "orthomcl+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::OrthoMCL_Pipeline'},
        "orthofinder+" => \$todo_list->{todo}{'Bio::Adventure::Align::OrthoFinder'},
        "phageterm+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Phageterm'},
        "phanotate+" => \$todo_list->{todo}{'Bio::Adventure::Feature_Prediction::Phanotate'},
        "phastaf+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Phastaf'},
        "prodigal+" => \$todo_list->{todo}{'Bio::Adventure::Feature_Prediction::Prodigal'},
        "parseblast+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Parse_Blast'},
        "parsefasta+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Parse_Fasta_Global'},
        "parsebcf+" => \$todo_list->{todo}{'Bio::Adventure::SNP::SNP_Ratio'},
        "phagepromoter+" => \$todo_list->{todo}{'Bio::Adventure::Feature_Prediction_Phagepromoter'},
        "polya+" => \$todo_list->{todo}{'Bio::Adventure::Trim::PolyA_Extractor'},
        "posttrinity+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity_Post'},
        "prokka+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Prokka'},
        "queryjob+" => \$todo_list->{todo}{'Bio::Adventure::Slurm::Query_Job'},
        "racer+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Racer'},
        "readsample+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Read_Samples'},
        "recatalog+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Restriction_Catalog'},
        "recurse+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Recurse'},
        "resfinder+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Resfinder'},
        "restriction+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Restriction_Catalog'},
        "rhopredict+" => \$todo_list->{todo}{'Bio::Adventure::Feature_Prediction::Rho_Predict'},
        "rgi+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Rgi'},
        "rmats+" => \$todo_list->{todo}{'Bio::Adventure::Splicing::RMats'},
        "rsem+" => \$todo_list->{todo}{'Bio::Adventure::Map::RSEM'},
        "runessentiality+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Run_Essentiality'},
        "sam2bam+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Sam2Bam'},
        "salmon+" => \$todo_list->{todo}{'Bio::Adventure::Map::Salmon'},
        "terminasereorder+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Terminase_ORF_Reorder'},
        "shovill+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Shovill'},
        "slseq_record+" => \$todo_list->{todo}{'Bio::Adventure::Splicing::SLSeq_Recorder'},
        "slsearch+" => \$todo_list->{todo}{'Bio::Adventure::Splicing::SLSearch'},
        "slutr+" => \$todo_list->{todo}{'Bio::Adventure::Splicing::SL_UTR'},
        "snippy+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Snippy'},
        "alignsnpsearch+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Align_SNP_Search'},
        "snpsearch+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Mpileup_SNP_Search'},
        "snpratio+" => \$todo_list->{todo}{'Bio::Adventure::SNP::SNP_Ratio'},
        "snpgenome+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Make_Genome'},
        "sortindexes+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Sort_Indexes'},
        "spladder+" => \$todo_list->{todo}{'Bio::Adventure::Splicing::Spladder'},
        "splitalign+" => \$todo_list->{todo}{'Bio::Adventure::Align::Split_Align'},
        "sradownload+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Download_SRA_PRJNA'},
        "star+" => \$todo_list->{todo}{'Bio::Adventure::Map::STAR'},
        "suppa+" => \$todo_list->{todo}{'Bio::Adventure::Splicing::Suppa'},
        "tacheck+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::TA_Check'},
        "test+" => \$todo_list->{todo}{'Bio::Adventure::Test_Job'},
        "tophat+" => \$todo_list->{todo}{'Bio::Adventure::Map::Tophat'},
        "tpp+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Transit_TPP'},
        "trainprodigal+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Train_Prodigal'},
        "transdecoder+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Transdecoder'},
        "trimomatic+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Trimomatic'},
        "trinity+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity'},
        "trinitypost+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity_Post'},
        "trinotate+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Trinotate'},
        "trnascan+" => \$todo_list->{todo}{'Bio::Adventure::Feature_Prediction::tRNAScan'},
        "umidedup+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Umi_Tools_Dedup'},
        "umitools+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Umi_Tools_Extract'},
        "unicycler+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Unicycler'},
        "variantgenome+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Make_Genome'},
        "velvet+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Velvet'},
        "vienna+" => \$todo_list->{todo}{'Bio::Adventure::Structure::RNAFold_Windows'},
        "rosalindplus+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Rosalind_Plus'},
        "xrefcrispr+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Xref_Crispr'},
        ## A set of pipelines combining the above.
        "pannotate+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Annotate_Assembly'},
        "passemble+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Assemble'},
        "pbt1+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Bowtie'},
        "pbt2+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Bowtie2'},
        "pbwa+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::BWA'},
        "pdnaseq+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Process_DNASeq'},
        "phageassemble+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Phage_Assemble'},
        "phisat+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Hisat'},
        "pkallisto+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::RNAseq_Pipeline_Kallisto'},
        "priboseq+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Riboseq_Pipeline'},
        "prnaseq+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Process_RNAseq'},
        "psalmon+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Salmon'},
        "ptnseq+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::TNseq_Pipeline'},
        "ptophat+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::RNAseq_Pipeline_Tophat'},
        "help+" => \$todo_list->{todo}{'Bio::Adventure::Adventure_Help'},
    };

    if (defined($args{method})) {
        my $m = '--' . $args{method};
        my @method = ($m,);
        my $array_result = GetOptionsFromArray(\@method, %{$possible_todos});
    }

    my @methods_to_run = ();
    my $set_of_todos = $todo_list->{todo};
    foreach my $job (keys %{$set_of_todos}) {
        my $thing = $set_of_todos->{$job};
        if ($thing) {
            if (!defined($args{task})) {
                print "Going on an adventure to ${job}.\n";
            }
            else {
                print "Going on an adventure to ${job} in the $args{task} context.\n";
            }
            my $final = \&${job};
            push(@methods_to_run, $final);
        }
    }
    return(\@methods_to_run);
}

1;
