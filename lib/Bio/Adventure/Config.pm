package Bio::Adventure::Config;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use v5.14;
use Cwd qw"abs_path fastcwd";
use Data::Dumper;
use Exporter;
use File::Basename qw"basename dirname";
use FileHandle;
use File::Path qw"make_path";
use File::ShareDir qw"dist_file module_dir dist_dir";
use Getopt::Long qw"GetOptionsFromArray";
use JSON;
use POSIX qw"ceil";
use String::Interpolate qw"interpolate";
use Template;
my $template_base = dist_dir('Bio-Adventure');
my $template_dir = qq"${template_base}/jobs";

our @EXPORT = qw"Get_Modules Get_Menus Get_TODOs";
our @ISA = qw"Exporter";

sub Add_Config_File {
    my %args = @_;
    my $input = $args{input};
    my %current = %{$args{values}};
    my $config = AppConfig->new({ CREATE => '.*' , GLOBAL => {
        ARGCOUNT => AppConfig::ARGCOUNT_ONE } });
    my $info = $config->file($input);
    my $variable = q{};
    my $steps = {};
    my $defaults = {};
    my $step = '0';
    my %data = $config->varlist(".*");
    for my $config_option (keys %data) {
        my $step = 'defaults';
        $variable = $config_option;
        my $config_value = $data{$config_option};
        if ($config_option =~ /^defaults_/) {
            $variable =~ s/^defaults_//g;
        }
        if ($config_option =~ /^(\d+)_(.+)$/) {
            $step = $1;
            $variable = $2;
            if (!defined($steps->{$step})) {
                $steps->{$step} = {};
            }
            $steps->{$step}->{$variable} = $config_value;
        } else {
            $defaults->{$variable} = $config_value;
        }
    }

    my $new = Merge_Current_Ini(current => %current, steps => $steps, defaults => $defaults);


    return($new);
}

sub Get_ARGV {
    my $cmdline_args = {};
    if ($#ARGV <= 0) {
        return $cmdline_args;
    }
    my $arg_string = q{};
    my @spec_strings = ();
    foreach my $a (@ARGV) {
        $arg_string .= "$a ";
        if ($a =~ /^\-/) {
            my $spec = $a;
            $spec =~ s/^\-*//g;
            push(@spec_strings, $spec);
        }
    }

    my %argv_information = ();
    my %conf_specification = ();
    for my $spec (@spec_strings) {
        $argv_information{$spec} = undef;
        my $def_string = qq"${spec}:s";
        $conf_specification{$def_string} = \$argv_information{$spec};
    }
    my $argv_result = GetOptions(%conf_specification);
    $cmdline_args = \%argv_information;
    $cmdline_args->{argv_string} = $arg_string;
    return($cmdline_args);
}

=head2 C<Get_Caller>

Figure out the calling function so that I may fill in some configuration information.

=cut
sub Get_Caller {
    my %args = @_;
    my $caller_number = 2;
    $caller_number = $args{caller} if ($args{caller});
    my $attempts = 4;
    my $found_sub = undef;
    my $function;
    my $last_attempt = $caller_number + $attempts;
  COUNTER: for my $attempt ($caller_number .. $last_attempt) {
        my ($package, $filename, $line, $sub, $hasargs,
            $wantarray, $evaltext, $is_require, $hints, $bitmask, $hinthash) = caller($attempt);
        last COUNTER unless (defined($sub));
        if ($sub =~ /^Bio::.*Submit$/) {
            next COUNTER;
        } elsif ($sub =~ /^Bio::.*Get_Vars$/) {
            next COUNTER;
        } elsif ($sub =~ /^Bio::.*Get_Caller$/) {
            next COUNTER;
        }
        ## All callers which are needing modules have the form Bio::Adventure::Category::Function
        $function = $sub;
        $function =~ s/^(.*::)(\w+)$/$2/g;
        last COUNTER if (defined($function));
    }
    if (!defined($function)) {
        $function = 'unknown';
    }
    return($function);
}

## My thought here: Make it possible to have a set of 'global' defaults
## which may be overridden at a few levels, branching down to the indivdual sample.
## I organize my data as /type-of-project/individual-project/preprocessing/sample/
## so there are likely 4 levels at which one may most likely want to make changes.
sub Get_Ini {
    my %args = @_;
    my $filename = basename($args{input_file});
    my $ini_values = {};
    ## At the super-project level, e.g. all the RNASeq analyses may share defaults
    ## which are different from the TNSeq etc.
    my $super_super_super_file = qq"../../../${filename}";
    if (-f $super_super_super_file) {
        $ini_values = Add_Config_File(
            input => $super_super_super_file,
            values => $ini_values,
        );
    }
    ## At the project directory, which is likely redundant to the preprocessing directory
    ## but I can imagine having a separate directory for other types of analysis
    my $super_super_file = qq"../../${filename}";
    if (-f $super_super_file) {
        $ini_values = Add_Config_File(
            input => $super_super_file,
            values => $ini_values,
        );
    }
    ## At the preprocessing directory
    my $super_file = qq"../${filename}";
    if (-f $super_file) {
        $ini_values = Add_Config_File(
            input => $super_file,
            values => $ini_values,
        );
    }
    ## Change defaults at the sample level
    my $file = qq"${filename}";
    if (-f $file) {
        $ini_values = Add_Config_File(
            input => $file,
            values => $ini_values,
        );
    }
    return($ini_values);
}

sub Get_JSON {
    my %args = @_;
    if (!defined($args{input_file})) {
        return({});
    }
    my $default_base = basename($args{input_file});
    if (!defined($default_base)) {
        return({});
    }
    my $input_file;
    my $dist_dir = dist_dir('Bio-Adventure');
    my $dist_file = qq"${dist_dir}/${default_base}";
    if (-r $args{input_file}) {
        $input_file = $args{input_file};
    } elsif (-f $dist_file) {
        $input_file = $dist_file;
    } else {
        say "Unable to find global default file, returning an empty hash ref.";
        return({});
    }
    my $json_string = q{};
    my $in = FileHandle->new(qq"<${input_file}");
    while (my $line = <$in>) {
        $json_string .= $line;
    }
    $in->close();
    my $json = JSON->new->allow_nonref;
    $json = $json->relaxed([1]);
    my $information = $json->decode($json_string);
    return($information);
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
        Compression => {
            name => 'compress',
            message => 'There is a time when the operation of the machine becomes so odious,
makes you so sick at heart, that you can’t take part; you can’t even passively take part,
and you’ve got to put your bodies upon the gears and upon the wheels, upon the levers,
upon all the apparatus, and you’ve got to make it stop. Go to page 2025.',
            choices => {
                '(spring): Compress fastq data with spring.' => \&Bio::Adventure::Compress::Spring,
            },
        },
        Conversion => {
            name => 'convert',
            message => 'And it rained a fever. And it rained a silence. And it rained a sacrifice. And it rained a miracle. And it rained sorceries and saturnine eyes of the totem.  Go to page 2584981.',
            choices => {
                '(any2any): Convert from one format to another with Bio::SeqIO' => \&Bio::Adventure::Convert::Any2Any,
                '(sam2bam): Convert a sam mapping to compressed/sorted/indexed bam.' => \&Bio::Adventure::Convert::Sam2Bam,
                '(bam2cov): Convert a sorted bam to tsv of coverage/base.' => \&Bio::Adventure::Convert::Bam2Coverage,
                '(gb2gff): Convert a genbank flat file to gff/fasta files.' => \&Bio::Adventure::Convert::Gb2Gff,
                '(gff2fasta): Convert a gff file to a fasta file.' => \&Bio::Adventure::Convert::Gff2Fasta,
                '(gatkdedup): Deduplicate an existing bam alignment with gatk.' => \&Bio::Adventure::Convert::GATK_Dedup,
                '(extract_seqids): Extract arbitrary sequence IDs from a reference.' => \&Bio::Adventure::Convert::Extract_SeqIDs,
                '(extract_subseq): Extract a subsequence from a reference genome.' => \&Bio::Adventure::Convert::Extract_Subseq,
                '(split_fasta): Separate a fasta file into its pieces.' => \&Bio::Adventure::Convert::Split_Fasta,
            },
        },
        Counting => {
            name => 'count',
            message => 'Once men turned their thinking over to machines in the hope that this would set them free. But that only permitted other men with machines to enslave them.  Go to page 27812',
            choices => {
                '(fcount): Use featureCounts to count aligned reads.' => \&Bio::Adventure::Count::Feature_Counts,
                '(gfftree): Count features by type in a gff file.' => \&Bio::Adventure::Count::GFF_Tree,
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
            message => q{},
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
            message => 'It is important, when killing a nun, to ensure that you bring
 an army of sufficient size.  Go to page 299792458.',
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
            message => q{},
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
            message => q{},
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
                '(downloadassembly): Download an assembly from ncbi.' => \&Bio::Adventure::Prepare::Download_NCBI_Assembly,
            },
        },
        QA => {
            name => 'qa',
            message => 'There is a time when the operation of the machine becomes so
 odious that you cannot take part.  Go to page 91093837139',
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
            message => qq"When my god comes back I'll be waiting for him with a shotgun.
And I'm keeping the last shell for myself. (inexact quote)  Go to page 667408",
            choices => {
                '(mpileup): Count coverage with mpileup.' => \&Bio::Adventure::Count::Mpileup,
                '(trim): Trim sequence with an additional rule to remove the first 10 nucleotides.' => \&Bio::Adventure::Trim::Trimomatic,
                '(bwa): Map reads with bwa and count with htseq.' => \&Bio::Adventure::Map::BWA,
                '(bowtie): Map trimmed reads with bowtie1 and count with htseq.' => \&Bio::Adventure::Map::Bowtie,
                '(bt2): Map trimmed reads with bowtie2 and count with htseq.' => \&Bio::Adventure::Map::Bowtie2,
                '(dantools_rnseq): Invoke dantools using RNASeq data as input.' => \&Bio::Adventure::SNP::Dantools_RNASeq,
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
            message => 'I decided that Orion needed to die after the second time he saved my life.
Go to page 160669',
            choices => {
                '(rmats): Quantify changes in transcript splicing events across samples.' => \&Bio::Adventure::Splicing::RMats,
                '(suppa): Quantify changes in transcript splicing events across samples.' => \&Bio::Adventure::Splicing::Suppa,
                '(sl_recorder): Record observed spliced leaders in an SLSeq experiment.' => \&Bio::Adventure::Splicing::SL_Recorder,
                '(polya_recorder): Record observed spliced leaders in an SLSeq experiment.' => \&Bio::Adventure::Map::PolyA_Recorder,
                '(slsearch): Count frequency of SL (or an arbitrary) sequences.' => \&Bio::Adventure::Splicing::SLSearch,
                '(slutr): Search for SL junction reads and use them to define UTRs.' => \&Bio::Adventure::Splicing::SL_UTR,
            },
        },
        Structure => {
            name => 'structure',
            message => '... a form of behavior marked by obsessive preoccupation with community
decline, humiliation, or victimhood and by compensatory cults of unity, energy, and purity,
in which a mass-based party of committed nationalist militants, working in uneasy but
effective collaboration with traditional elites, abandons democratic liberties and pursues
with redemptive violence and without ethical or legal restraints goals of internal
cleansing and external expansion.  Go to page: 192210',
            choices => {
                '(alphafold): Invoke alphafold on (currently) a peptide sequence.' => \&Bio::Adventure::Structure::ProteinFold,
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

sub Get_Modulecmd {
    ## Check that the module command is available as a bash function.
    ## I was initially going to do this in BUILD(), but I think that might mess up jobs
    ## which set the jprefix parameter.
    ## This command, if run via bash, requires a controlling pty in order to find the
    ## bash function 'module'.  slurm and the various managers will not provide that
    ## so, we can use typescript to allocate one post-facto.
    my $module_test = qq!/usr/bin/script -q -c 'bash -i -c "type module"' /dev/null |!;
    my $modulecmd_handle = FileHandle->new($module_test);
    my $modulecmd_text = q{};
    MOD: while (my $line = <$modulecmd_handle>) {
          next MOD if ($line =~ /is a function/);
          $modulecmd_text .= $line;
      }
    $modulecmd_handle->close();
    if (-f 'typescript') {
        print "The stupid typescript file still got created.\n";
        rm('typescript');
    }
      return($modulecmd_text);
}

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

=head2 C<Get_Input_Info>

Given a file/directory name, provide a set of paths which will
hopefully prove useful for the various functions in CYOA.  I think
this function is a good candidate for replacing Check_Input() below.

=cut
sub Get_Input_Info {
    my ($class, @inputs) = @_;
    my %ret = ();
    return(undef) unless($inputs[0]);
    if ($inputs[0] =~ /[:;,\s]+/) {
        @inputs = split(/[:;,\s]+/, $inputs[0]);
    }
    my $num_inputs = scalar(@inputs);
    if ($num_inputs == 0) {
        die("This requires an input filename.");
    }
    my @outputs = ();
    for my $in (@inputs) {
        my $filename = basename($in);
        my $filebase_compress = basename($in, ('.gz', '.xz', '.bz2'));
        my @exts = ('.fastq', '.fasta', '.fsa', '.faa', 'fna', '.fa', '.ffn',
                    'gff', '.gff3', '.gbk', '.gbf', '.sqn', 'tbl', '.tsv', '.txt');
        my $filebase_extension = basename($filebase_compress, @exts);
        my $directory = dirname($in);
        my $dirname = basename($directory);
        ## Check that doing basename(dirname(input)) leaves us with something interesting
        if ($dirname eq '.') {
            $dirname = undef;
        }
        ## This test/path build is because abs_path only works on stuff which exists.
        if (! -e $directory) {
            make_path($directory);
        }
        my $full_path = abs_path($directory);
        $full_path .= "/${filename}";
        my $local_path = qq"${directory}/${filename}";

        ## Make a new job_basename which doesn't have some cruft
        my $jbasename = $filebase_extension;
        $jbasename =~ s/_forward.*//g;
        $jbasename =~ s/_R1.*//g;
        $jbasename =~ s/-trimmed.*//g;
        $jbasename =~ s/_trimmed.*//g;

        my %ret = (
            filename => $filename,
            filebase_compress => $filebase_compress,
            filebase_extension => $filebase_extension,
            directory => $directory,
            dirname => $dirname,
            localpath => $local_path,
            fullpath => $full_path,
            jbasename => $jbasename,
            basename => $jbasename,
            input => $local_path,);
        push(@outputs, \%ret);
    }
    return(\@outputs);
}

=head C<Get_Memory>

In my time estimator, I tried out some heuristics and have found them
somewhat difficult to rely upon because of the huge variance in cpu speeds.
I thought I would therefore try some simpler heuristics for memory usage
especially since memory usage tends to be a bit more consistent
which of course may be an argument for a multiplier based system.

=cut
sub Get_Memory {
    my %args = @_;
    my $function = Bio::Adventure::Config::Get_Caller();
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

=head2 C<Get_Outputs>

Get_Paths seeks out a set of useful paths for each function, this defines output filenames.
This largely depends on Get_Paths.  One could easily argue that it should be a part of that
function; but I separated it because I like having a conceptual separation between the various
locations for stuff and the final outputs.

=cut
sub Get_Outputs {
    my (%args) = @_;
    my $subroutine = $args{subroutine};
    my $options = $args{options};
    my $paths = $args{paths};
    if (!defined($subroutine)) {
        $subroutine = Bio::Adventure::Config::Get_Caller();
    }

    my $outputs = {};
    given ($subroutine) {
        when ('Any2Any') {
            $outputs->{main} = qq"$paths->{input_dirname}/$paths->{output_base}.$options->{outformat}";
        }
        when ('Bam2Bigwig') {
            $outputs->{main} = qq"$paths->{output_dir}/$paths->{output_base}.bw";
        }
        when ('Bedtools_Coverage') {
            ## Not currently using anything
        }
        when ($_ eq 'Bowtie' || $_ eq 'BT1_Index') {
            $outputs->{stdout} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{species}_$options->{bt_type}.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{species}_$options->{bt_type}.stderr";
            $outputs->{aligned} = qq"$paths->{output_dir}/$paths->{output_base}_aligned_$options->{species}_$options->{libtype}.fastq";
            $outputs->{unaligned} = qq"$paths->{output_dir}/$paths->{output_base}_unaligned_$options->{species}_$options->{libtype}.fastq";
            $outputs->{main} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{species}_$options->{libtype}.sam";
        }
        when ('BT1_Stats') {
            $outputs->{main} = qq"$paths->{output_dir}/bt1_stats.csv";
        }
        when($_ eq 'Bowtie2' || $_ eq 'BT2_Index') {
            $outputs->{stdout} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{species}_$options->{bt_type}.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{species}_$options->{bt_type}.stderr";
        }
        when ($_ eq 'BWA' || $_ eq 'BWA_Index') {
            $outputs->{fa} = qq"$paths->{index_dir}/$options->{species}.fa";
            $outputs->{sa} = qq"$outputs->{fa}.sa";
            $outputs->{pac} = qq"$outputs->{fa}.pac";
            $outputs->{bwt} = qq"$outputs->{fa}.bwt";
            $outputs->{ann} = qq"$outputs->{fa}.ann";
            $outputs->{amb} = qq"$outputs->{fa}.amb";
        }
        when ($_ eq 'Caical_Worker' || $_ eq 'Caical') {
            ## No output defined yet
        }
        when ('Check_Blastdb') {
            $outputs->{stderr} = 'makeblastdb.stderr';
            $outputs->{stdout} = 'makeblastdb.stdout';
        }
        when ($_ eq 'Classify_Phage' || $_ eq 'Classify_Phage_Worker') {
            ## No output defined yet
        }
        when ($_ eq 'Compress' || $_ eq 'Recompress') {
            $outputs->{stdout} = qq"$paths->{output_dir}/xz.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/xz.stderr";
        }
        when ($_ eq 'Concatenate_Searches' || $_ eq 'Concatenate_Searches_Worker') {
            ## The _Worker is wrong
            $outputs->{main} = qq"$paths->{output_dir}/split_search.txt.xz";
            $outputs->{stderr} = qq"$paths->{output_dir}/concat.stderr";
            $outputs->{stdout} = qq"$paths->{output_dir}/concat.stdout";
        }
        when ('Cutadapt') {
            $outputs->{main} = qq"$paths->{output_dir}/$paths->{output_base}$paths->{output_suffix}.$paths->{input_suffix}";
        }
        when ('Cutadapt_Stats') {
            $outputs->{main} = qq"$paths->{output_dir}/cutadapt_stats.csv";
        }
        when (/^Dantools/) {
            $outputs->{log} = qq"$paths->{output_dir}/dantools_rnaseq.log";
        }
        when ('Download_NCBI_Assembly') {
            $outputs->{log} = qq"$paths->{output_dir}/download_assembly.log";
        }
        when ('Downsample_Guess_Strand') {
            ## No output defined yet
        }
        when ('Duplicate_Remove') {
            $outputs->{main} = qq"$paths->{input_dirname}/$paths->{output_base}_removed.txt";
        }
        when ('Extend_Kraken_DB') {
            $outputs->{stdout} = qq"$paths->{output_dir}/kraken2-build.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/kraken2-build.stderr";
            $outputs->{main} = qq"$paths->{output_dir}/kraken2-build.out";
        }
        when ('Extract_Subseq') {
            $outputs->{main} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{start}_$options->{end}.$options->{outformat}";
        }
        when ('Fastp') {
            ## No output defined yet
        }
        when ('Fastp_Stats') {
            $outputs->{main} = qq"$paths->{output_dir}/fastp_stats.csv";
        }
        when ('Fastqc') {
            $outputs->{main} = qq"$paths->{output_subdir}/fastqc_data.txt";
            $outputs->{summary} = qq"$paths->{output_subdir}/summary.txt";
        }
        when ('Fastqc_Stats') {
            $outputs->{main} = qq"$paths->{output_dir}/fastqc_stats.csv";
        }
        when ('Feature_Counts') {
            ## This if statement is repeated in Count::FeatureCounts
            ## and happens because htseq uses yes/no/f/r while fc uses
            ## 0,1,2 and I do not expect the user to remember this.
            ## but I want the outputs to reflect correctly the choice.
            if (defined($options->{stranded})) {
                if ($options->{stranded} eq 'forward') {
                    $options->{stranded} = 1;
                } elsif ($options->{stranded} eq 'reverse') {
                    $options->{stranded} = 2;
                } else {
                    $options->{stranded} = 0;
                }
            } else {
                $options->{stranded} = 0;
            }
            $outputs->{main} = qq"$paths->{nosuffix}_s$options->{stranded}_$options->{gff_type}_$options->{gff_tag}_fcounts.csv";
            $outputs->{featurecount} = qq"$paths->{nosuffix}.bam.featureCounts";
            $outputs->{main_variables} = qq"$paths->{nosuffix}_s\${stranded}_\${gff_type}_\${gff_tag}_fcounts.csv";

            $outputs->{stdout} = qq"$paths->{nosuffix}_featureCounts.stdout";
            $outputs->{stderr} = qq"$paths->{nosuffix}_featureCounts.stderr";
        }
        when ($_ eq 'Filter_Host_Kraken' || $_ eq 'Filter_Kraken_Worker') {
            $outputs->{stderr} = qq"$paths->{output_dir}/filter_host.stderr";
            $outputs->{stdout} = qq"$paths->{output_dir}/filter_host.stdout";
            $outputs->{log} = qq"$paths->{output_dir}/filter_kraken.log";
        }
        when ($_ eq 'Freebayes_SNP_Search' || $_ eq 'SNP_Ratio' || $_ eq 'SNP_Ratio_Worker' ||
              $_ eq 'SNP_Ratio_Intron' || $_ eq 'SNP_Ratio_Intron_Worker') {
            ## No output defined yet
        }
        when ('GATK_Dedup') {
            $outputs->{marked} = qq"$paths->{output_dir}/deduplication_stats.txt";
            $outputs->{stdout} = qq"$paths->{output_dir}/deduplication.stdout";
            $outputs->{stder} = qq"$paths->{output_dir}/deduplication.stderr";
            $outputs->{main} = qq"$paths->{output_dir}/deduplicated.bam";
        }
        when (/^Gb2Gff/) {
            $outputs->{stdout} = qq"$paths->{output_dir}/$paths->{output_base}.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/$paths->{output_base}.stderr";
            $outputs->{fasta} = qq"$paths->{output_dir}/$paths->{output_base}.fsa";
            $outputs->{cds_fasta} = qq"$paths->{output_dir}/$paths->{output_base}.ffn";
            $outputs->{misc_fasta} = qq"$paths->{output_dir}/$paths->{output_base}_miscRNA.fasta";
            $outputs->{nc_fasta} = qq"$paths->{output_dir}/$paths->{output_base}_ncRNA.fasta";
            $outputs->{pep_fasta} = qq"$paths->{output_dir}/$paths->{output_base}.faa";
            $outputs->{rrna_fasta} = qq"$paths->{output_dir}/$paths->{output_base}_rRNA.fasta";
            $outputs->{all_gff} = qq"$paths->{output_dir}/$paths->{output_base}_all.gff";
            $outputs->{gene_gff} = qq"$paths->{output_dir}/$paths->{output_base}.gff";
            $outputs->{cds_gff} = qq"$paths->{output_dir}/$paths->{output_base}_cds.gff";
            $outputs->{inter_gff} = qq"$paths->{output_dir}/$paths->{output_base}_inter_cds.gff";
            $outputs->{misc_gff} = qq"$paths->{output_dir}/$paths->{output_base}_miscRNA.gff";
            $outputs->{nc_gff} = qq"$paths->{output_dir}/$paths->{output_base}_ncRNA.gff";
            $outputs->{rrna_gff} = qq"$paths->{output_dir}/$paths->{output_base}_rRNA.gff";
        }
        when ('Get_Substitution_Matrix') {
            $outputs->{matrix_file} = qq"$paths->{output_dir}/$options->{matrix}";
        }
        when ('GFF_Tree') {
            ## No output defined yet
        }
        when ('Gff2Fasta') {
            if (defined($options->{species})) {
                $outputs->{aa} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{species}_$options->{gff_type}_$options->{gff_tag}.faa";
                $outputs->{nt} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{species}_$options->{gff_type}_$options->{gff_tag}.ffn";
                $outputs->{gff} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{species}_$options->{gff_type}_$options->{gff_tag}.gff";
            } else {
                $outputs->{aa} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{gff_type}_$options->{gff_tag}.faa";
                $outputs->{nt} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{gff_type}_$options->{gff_tag}.ffn";
                $outputs->{gff} = qq"$paths->{output_dir}/$paths->{output_base}_$options->{gff_type}_$options->{gff_tag}.gff";
            }
        }
        when (/^Guess_Strand/) {
            $outputs->{main} = qq"$paths->{nosufix}.txt";
        }
        when ($_ eq 'Hisat2' || $_ eq 'Hisat2_Index') {
            $outputs->{index_stderr} = qq"$paths->{output_dir}/hisat2_index_$options->{species}.stderr";
            $outputs->{index_stdout} = qq"$paths->{output_dir}/hisat2_index_$options->{species}.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/hisat_$options->{species}_$options->{libtype}";
            $outputs->{stdout} = $paths->{stderr};
            $outputs->{metrics} = $paths->{stderr};
            if (defined($options->{jbasename})) {
                $outputs->{stderr} .= qq"_$options->{jbasename}.stderr";
                $outputs->{stdout} .= qq"_$options->{jbasename}.stdout";
                $outputs->{metrics} .= qq"_$options->{jbasename}.metrics";
            } else {
                $outputs->{stderr} .= '.stderr';
                $outputs->{stdout} .= '.stdout';
                $outputs->{metrics} .= '.metrics';
            }
            $outputs->{aligned_discordant} = qq"$paths->{output_dir}/aldis_$options->{species}_$options->{libtype}.fastq";
            if ($options->{unaligned_discordant}) {
                $outputs->{unaligned_discordant} = $options->{unaligned_discordant};
            } else {
                $outputs->{unaligned_discordant} = qq"$paths->{output_dir}/unaldis_$options->{species}_$options->{libtype}.fastq";
            }
            $outputs->{aligned_concordant} = qq"$paths->{output_dir}/alcon_$options->{species}_$options->{libtype}.fastq";
            if (defined($options->{unaligned_output})) {
                $outputs->{unaligned_concordant} = $options->{unaligned_output};
            } else {
                $outputs->{unaligned_concordant} = qq"$paths->{output_dir}/unalcon_$options->{species}_$options->{libtype}.fastq";
            }
            $outputs->{main} = qq"$paths->{output_dir}/$options->{species}_$options->{libtype}.sam";
        }
        when ('Hisat2_Index') {
            ## No output defined yet
        }
        when ('Hisat2_Stats') {
            $outputs->{main} = qq"$paths->{output_dir}/hisat_stats.csv";
        }
        when ($_ eq 'HT_Multi' || $_ eq 'HTSeq' || $_ eq 'HT_Types') {
            $outputs->{output_base} = qq"$paths->{output_dir}/$paths->{nosuffix}";
            $outputs->{main} = qq"$paths->{output_dir}/$paths->{output_base}";
            $outputs->{stderr} = qq"$paths->{output_dir}/$paths->{output_base}_htseq.stderr";
            $outputs->{stdout} = qq"$paths->{output_dir}/$paths->{output_base}_htseq.stdout";
        }
        when ('Insert_Size') {
            $outputs->{stderr} = qq"$paths->{output_dir}/$paths->{output_base}_insert_size.stderr";
            $outputs->{stdout} = qq"$paths->{output_dir}/$paths->{output_base}_insert_size.stdout";
            $outputs->{main} = qq"$paths->{output_dir}/$paths->{output_base}_insert_size.txt";
            $outputs->{pdf} = qq"$paths->{output_dir}/$paths->{output_base}_insert_size.pdf";
        }
        when ($_ eq 'Jellyfish' || $_ eq 'Jellyfish_Matrix') {
            $outputs->{count} = qq"$paths->{jelly_base}.count";
            $outputs->{info} = qq"$paths->{jelly_base}.info";
            $outputs->{histogram} = qq"$paths->{jelly_base}.hist";
            $outputs->{count_fasta} = qq"$paths->{jelly_base}_by_count.fasta";
            $outputs->{matrix} = qq"$paths->{jelly_base}_matrix.csv";
            $outputs->{stderr} = qq"$paths->{jelly_base}.stderr";
            $outputs->{stdout} = qq"$paths->{jelly_base}.stdout";
            $outputs->{main} = $outputs->{histogram};
        }
        when ($_ eq 'Kallisto' || $_ eq 'Kallisto_Index') {
            $outputs->{index_stderr} = qq"$paths->{index_dir}/index.stderr";
            $outputs->{index_stdout} = qq"$paths->{index_dir}/index.stdout";
        }
        when ('Kraken') {
            $outputs->{stdout} = qq"$paths->{output_dir}/kraken.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/kraken.stderr";
            $outputs->{report} = qq"$paths->{output_dir}/kraken_report.txt";
            $outputs->{main} = $outputs->{report};
            $outputs->{classified} = qq"$paths->{output_dir}/classified#.fastq.gz";
            $outputs->{unclassified} = qq"$paths->{output_dir}/unclassified#.fastq.gz";
        }
        when (/^Kraken_to_Matrix/) {
            $outputs->{main} = qq"$paths->{nosuffix}_matrix.tsv";
            $outputs->{stderr} = qq"$paths->{output_dir}/kraken_to_matrix.stderr";
            $outputs->{stdout} = qq"$paths->{output_dir}/kraken_to_matrix.stdout";
            $outputs->{log} = qq"$paths->{output_dir}/kraken_to_matrix.log";
        }
        when ('Make_Blast_Job') {
            $outputs->{stderr} = qq"$paths->{output_dir}/$paths->{output_base}_blast_job.stderr";
            $outputs->{stdout} = qq"$paths->{output_dir}/$paths->{output_base}_blast_job.stdout";
            $outputs->{main} = qq"$paths->{output_dir}/$options->{blast_tool}.out";
        }
        when ('Make_Codon_Table') {
            $outputs->{codon_table} = qq"$paths->{codon_dir}/$options->{species}.txt";
        }
        when ('Mash') {
            $outputs->{stdout} = qq"$paths->{output_dir}/mash_$paths->{cwd_name}.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/mash_$paths->{cwd_name}.stderr";
            $outputs->{dist_dir} = qq"$paths->{output_dir}/dist";
            $outputs->{sketch_dir} = qq"$paths->{output_dir}/sketch";
        }
        when ('Mauve_Aligner') {
            $outputs->{islands} = qq"$paths->{output_dir}/islands.txt";
            $outputs->{backbone} = qq"$paths->{output_dir}/backbone.txt";
            $outputs->{id_matrix} = qq"$paths->{output_dir}/id_matrix.txt";
            $outputs->{guide_tree} = qq"$paths->{output_dir}/guide_tree.txt";
            $outputs->{main} = qq"$paths->{output_dir}/alignment.txt";
        }
        when (/^Merge_Annotations/) {
            ## No output defined yet
        }
        when ($_ eq 'Mi_Map' || $_ eq 'Read_Bam_Mi' || $_ eq 'Final_Print_Mi') {
            $outputs->{pre_map} = qq"$paths->{nosuffix}_mirnadb.txt";
            $outputs->{final_map} = qq"$paths->{nosuffix}_mature.count";
        }
        when ('Mpileup') {
            $outputs->{stderr} = qq"$paths->{output_dir}/mpileup.stderr";
            $outputs->{stdout} = qq"$paths->{output_dir}/mpileup.stdout";
            $outputs->{bcferr} = qq"$paths->{output_dir}/mpileup_bcf.stdout";
            $outputs->{bam} = qq"$paths->{output_dir}/sorted.bam";
            $outputs->{vcf} = qq"$paths->{output_dir}/sorted.vcf";
            $outputs->{main} = qq"$paths->{output_dir}/sorted.bcf";
        }
        when ('Nextflow_RNASeq') {
            ## No output defined yet
        }
        when ($_ eq 'OrthoFinder' || $_ eq 'Orthofinder_Names_Worker' ||
              $_ eq 'Orthofinder_Names') {
            $outputs->{all} = qq"$paths->{output_dir}/output/Orthogroups/Orthogroups.tsv";
            $outputs->{single} = qq"$paths->{output_dir}/output/Orthogroups/Orthogroups_SingleCopyOrthologues.txt";
            $outputs->{named} = qq"$paths->{output_dir}/orthogroups_all_named.tsv";
            $outputs->{single_named} = qq"$paths->{output_dir}/orthogroups_single_named.tsv";
            $outputs->{namer_log} = qq"$paths->{output_dir}/orthonamer.log";
        }
        when ('Pairwise_Similarity_Matrix') {
            $outputs->{evalues} = qq"$paths->{output_dir}/pairwise_evalues.tsv";
            $outputs->{bit_scores} = qq"$paths->{output_dir}/pairwise_bits.tsv";
            $outputs->{pair_scores} = qq"$paths->{output_dir}/pairwise_scores.tsv";
            $outputs->{main} = qq"$paths->{output_dir}/output.txt";
        }
        when ('PolyA_Recorder') {
            ## No output defined yet
        }
        when ('Prodigal') {
            $outputs->{stdout} = qq"$paths->{output_dir}/prodigal.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/prodigal.stderr";
        }
        when ('Progressive_Mauve') {
            $outputs->{main} = qq"$paths->{output_dir}/alignment";
            $outputs->{guide_tree} = qq"$paths->{output_dir}/guide_tree";
            $outputs->{backbone} = qq"$paths->{output_dir}/backbone";
        }
        when ($_ eq 'ProteinFold' || $_ eq 'ProteinFold_Worker') {
            $outputs->{main} = qq"$paths->{output_dir}/proteinfold.txt";
        }
        when ($_ eq 'ProteinFold_PairIDs' || $_ eq 'ProteinFold_PairIDs_Worker') {
            $outputs->{main} = qq"$paths->{output_dir}/proteinfold.txt";
        }
        when (/^RNAFold_Windows/) {
            $outputs->{main} = qq"$paths->{output_dir}/$args{output_name}.tsv.xz";
            $outputs->{stdout} = qq"$paths->{output_dir}/$args{output_name}_rnafold.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/$args{output_name}_rnafold.stderr";
        }
        when ($_ eq 'RSEM' || $_ eq 'RSEM_Index') {
            ## No outputs defined yet
        }
        when ('Run_Parse_Blast') {
            $outputs->{counts} = qq"$paths->{output_dir}/counts.txt";
            $outputs->{zeros} = qq"$paths->{output_dir}/zeros.fasta";
            $outputs->{singles} = qq"$paths->{output_dir}/singles.fasta";
            $outputs->{doubles} = qq"$paths->{output_dir}/doubles.fasta";
            $outputs->{few} = qq"$paths->{output_dir}/few.fasta";
            $outputs->{many} = qq"$paths->{output_dir}/many.fasta";
            $outputs->{main} = qq"$paths->{output_dir}/blast_output.txt";
        }
        when ($_ eq 'Salmon' || $_ eq 'Salmon_Index') {
            $outputs->{main} = qq"$paths->{output_dir}/quant.sf";
            $outputs->{stdout} = qq"$paths->{output_dir}/salmon_$options->{species}.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/salmon_$options->{species}.stderr";
        }
        when ($_ eq 'Sam2Bam' || $_ eq 'Samtools') {
            $outputs->{main} = qq"$paths->{input_dirname}/$paths->{output_base}.bam";
            $outputs->{sorted} = qq"$paths->{input_dirname}/$paths->{output_base}-sorted.bam";
            $outputs->{sorted_nomismatch} = qq"$paths->{input_dirname}/$paths->{output_base}_nomismatch.bam";
            $outputs->{paired} = qq"$paths->{input_dirname}/$paths->{output_base}-paired.bam";
            $outputs->{stdout} = qq"$paths->{input_dirname}/$paths->{output_base}_samtools.stdout";
            $outputs->{stderr} = qq"$paths->{input_dirname}/$paths->{output_base}_samtools.stderr";
        }
        when ('SLSearch') {
            ## No output defined yet
        }
        when ('SL_Recorder') {
            ## No output defined yet
        }
        when ('Split_Align_Blast') {
            $outputs->{main} = qq"$paths->{output_dir}/$paths->{output_base}_vs_$paths->{libname}.txt";
        }
        when ($_ eq 'Split_Align_Fasta' || $_ eq 'Split_Align_Fasta_Worker' ||
              $_ eq 'Parse_Fasta' || $_ eq 'Parse_Fasta_Global' ||
              $_ eq 'Make_Fasta_Job' || $_ eq 'Parse_Fasta_Worker' ||
              $_ eq 'Parse_Fasta_Mismatches') {
            $outputs->{fasta} = qq"$paths->{output_dir}/$paths->{output_filename_base}.txt";
            $outputs->{fasta_comp} = qq"$outputs->{fasta}.xz";
            $outputs->{main} = qq"$paths->{output_dir}/$paths->{output_filename_base}.parsed.txt";
            $outputs->{counts} = qq"$paths->{output_dir}/$paths->{output_filename_base}.count";
            $outputs->{singles} = qq"$paths->{output_dir}/$paths->{output_filename_base}_singles.txt";
            $outputs->{doubles} = qq"$paths->{output_dir}/$paths->{output_filename_base}_doubles.txt";
            $outputs->{few} = qq"$paths->{output_dir}/$paths->{output_filename_base}_few.txt";
            $outputs->{many} = qq"$paths->{output_dir}/$paths->{output_filename_base}_many.txt";
            $outputs->{zero} = qq"$paths->{output_dir}/$paths->{output_filename_base}_zero.txt";
            $outputs->{all} = qq"$paths->{output_dir}/$paths->{output_filename_base}_all.txt";
            $outputs->{variants} = qq"$paths->{output_dir}/$paths->{output_filename_base}_mismatches.txt";
            $outputs->{variant_numbers} = qq"$paths->{output_dir}/$paths->{output_filename_base}_mismatch_nums.txt";
        }
        when ('Split_Fasta') {
            ## No output defined yet
        }
        when ('Spring') {
            $outputs->{stdout} = qq"$paths->{output_dir}/$paths->{output_base}_spring.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/$paths->{output_base}_spring.stderr";
        }
        when ($_ eq 'STAR' || $_ eq 'STAR_Index') {
            ## No output defined yet
        }
        when (/^Terminase_ORF_Reorder/) {
            $outputs->{stderr} = qq"$paths->{output_dir}/termreorder.stderr";
            $outputs->{stdout} = qq"$paths->{output_dir}/termreorder.stdout";
            $outputs->{main} = qq"$paths->{output_dir}/final_assembly.fasta";
            $outputs->{prodigal} = qq"$paths->{output_dir}/prodigal_cds.fasta";
        }
        when ('Test_Worker') {
            $outputs->{stderr} = '/tmp';
            $outputs->{stdout} = '/tmp/funkytown.stdout';
        }
        when ('Tophat') {
            ## No output defined yet
        }
        when ('Transit_TPP') {
            ## No output defined yet
        }
        when ('Trimomatic_Pairwise') {
            $outputs->{main} = join(':', @{$paths->{input_base_array}});
            $outputs->{stdout} = qq"$paths->{output_dir}/$paths->{output_base}-trimomatic.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/$paths->{output_base}-trimomatic.stderr";
        }
        when ('Trimomatic_Single') {
            $outputs->{main} = qq"$paths->{output_dir}/$paths->{output_base}-trimmed.fastq";
            $outputs->{stdout} = qq"$paths->{output_dir}/$paths->{output_base}-trimomatic.stdout";
            $outputs->{stderr} = qq"$paths->{output_dir}/$paths->{output_base}-trimomatic.stderr";
        }
        when ('Trimomatic_Stats') {
            $outputs->{main} = qq"$paths->{output_dir}/trim_stats.csv";
        }
        when ('Umi_Tools_Extract') {
            ## No output defined yet
        }
        when ('Umi_Tools_Dedup') {
            ## No output defined yet
        }
        when ('Unicycler') {
            ## No output defined yet
        }
        when ('Velvet') {
            ## No output defined yet
        }
    } ## End the given structure
    $outputs->{main} = qq"$paths->{output_dir}/default.out" unless (defined($outputs->{main}));
    $outputs->{stdout} = qq"$paths->{output_dir}/stdout" unless (defined($outputs->{stdout}));
    $outputs->{stderr} = qq"$paths->{output_dir}/stderr" unless (defined($outputs->{stderr}));
    $outputs->{log} = qq"$paths->{output_dir}/log.txt" unless (defined($outputs->{log}));
    return($outputs);
}

=head2 C<Get_Class_Paths>

This should be responsible for acquiring paths used by anything in
this instance of the Adventure class.

TODO: Once this is settled, take out the function-specific material
and move it to Get_Function_Paths() which will contain stuff
like output_dir.

=cut
sub Get_Class_Paths {
    my ($class, %args) = @_;
    my $spec = $args{spec};
    my $options = $class->current_state;
    my $libtype = 'genome';
    $libtype = $options->{libtype} if (defined($options->{libtype}));
    my $libpath_prefix = qq"$options->{libpath}/${libtype}";
    my $libdir_prefix = qq"$options->{libdir}/${libtype}";
    my $output_prefix = qq"$options->{output_base}/$spec->{jprefix}";
    ##my $cwd_name = basename(cwd());
    my $paths = {};
    my $outputs = {};
    my @input_array = ();
    my @output_array = ();
    my @output_base_array = ();
    if (defined($options->{input})) {
        @input_array = split(/$options->{delimiter}/, $options->{input});
        $paths->{input_suffix} = $options->{input};
        $paths->{input_suffix} =~ s/\.[xgb]z$//g;
        if ($paths->{input_suffix} =~ /\..*$/) {
            $paths->{input_suffix} =~ s/^.*\.(.*)$/$1/g;
        }
        ##$paths->{output_base} = basename($options->{input}, @{$options->{suffix_array}});
        ##$paths->{output_base} = basename($paths->{output_base}, @{$options->{suffix_array}});
        $paths->{input_dirname} = dirname($options->{input});
        if ($paths->{input_dirname} eq '.' || $paths->{input_dirname} eq q{}) {
            $paths->{input_dirname} = fastcwd();
        }
    } else {
        $paths->{input_dirname} = fastcwd();
        @input_array = ($options->{input});
    }
    local *basebase = sub {
        my $in = shift;
        if (!defined($in)) {
            warn("Expected a file path but received undefined.\n");
            return('undefined');
        }
        my $out = basename($in, @{$options->{suffix_array}});
        $out = basename($out, @{$options->{suffix_array}});
        return($out);
    };

    my @input_base_array = ();
    if (@input_array) {
        @input_base_array = map { basebase($_) } @input_array;
    } else {
        print "The input array was not defined in Config::Get_Paths.\n";
        @input_array = ('undefined');
    }

    $paths->{output_base} = $input_base_array[0];
    $paths->{input_array} = \@input_array;
    $paths->{input_base_array} = \@input_base_array;
    ## Ok, something in cwd/basename is causing a problem.
    $paths->{cwd_name} = basename(fastcwd());
    $paths->{nosuffix} = qq"$paths->{input_dirname}/$paths->{output_base}";
    $paths->{blast_dir} = qq"${libpath_prefix}/blast";
    $paths->{codon_dir} = qq"${libpath_prefix}/codon_tables";
    $paths->{fasta_dir} = qq"${libpath_prefix}/fasta";
    $paths->{gbk_dir} = qq"${libpath_prefix}/genbank";
    $paths->{gff_dir} = qq"${libpath_prefix}/gff";
    $paths->{gtf_dir} = qq"${libpath_prefix}/gtf";
    $paths->{hmm_dir} = qq"${libpath_prefix}/hmm";
    $paths->{index_prefix} = qq"${libpath_prefix}/indexes";
    $paths->{index_prefix_shell} = qq"${libdir_prefix}/indexes";
    $paths->{output_prefix} = $output_prefix;
    $paths->{output_dir} = qq"${output_prefix}_$spec->{jname_prefix}_$paths->{output_base}";
    $paths->{time} = qq"$paths->{output_dir}/$spec->{jname_prefix}_$paths->{output_base}.time";
    $paths->{jname} = $paths->{output_base};
    if (!defined($options->{jprefix})) {
        $options->{jprefix} = '01';
    }
    $options->{jname} = $paths->{jname};
    $paths->{time} = qq"$paths->{output_dir}/$options->{jname}.time";
    if ($spec->{jprefix}) {
        $spec->{jname} = qq"$spec->{jprefix}$spec->{jname_prefix}";
    }
    if ($options->{libtype}) {
        $options->{jname} .= "_$options->{libtype}";
    }
    make_path($paths->{output_dir}) unless (-d $paths->{output_dir});
    make_path($paths->{blast_dir}) unless (-d $paths->{blast_dir});
    make_path($paths->{codon_dir}) unless (-d $paths->{codon_dir});
    make_path($paths->{fasta_dir}) unless (-d $paths->{fasta_dir});
    make_path($paths->{gff_dir}) unless (-d $paths->{gff_dir});
    make_path($paths->{gtf_dir}) unless (-d $paths->{gtf_dir});
    make_path($paths->{gbk_dir}) unless (-d $paths->{gbk_dir});
    make_path($paths->{hmm_dir}) unless (-d $paths->{hmm_dir});
    make_path($paths->{index_prefix}) unless (-d $paths->{index_prefix});
    if (defined($ENV{BLAST_DB})) {
        if ($ENV{BLAST_DB} =~ /:/) {
            my @blast = split(/:/, $ENV{BLAST_DB});
            unless (grep(/$paths->{blast_dir}/, @blast)) {
                $ENV{BLAST_DB} = qq"$ENV{BLAST_DB}:$paths->{blast_dir}";
            }
        } else {
            unless ($paths->{blast_dir} eq $ENV{BLAST_DB}) {
                $ENV{BLAST_DB} = qq"$ENV{BLAST_DB}:$paths->{blast_dir}";
            }
        }
    } else {
        $ENV{BLAST_DB} = $paths->{blast_dir};
    }
    if (defined($options->{species})) {
        $paths->{fasta} = qq"$paths->{fasta_dir}/$options->{species}.fasta";
        $paths->{fasta_shell} = qq"${libdir_prefix}/fasta/$options->{species}.fasta";
        $paths->{gff} = qq"$paths->{gff_dir}/$options->{species}.gff";
        $paths->{gff_shell} = qq"${libdir_prefix}/gff/$options->{species}.gff";
        $paths->{gtf} = qq"$paths->{gtf_dir}/$options->{species}.gtf";
        $paths->{gtf_shell} = qq"${libdir_prefix}/gtf/$options->{species}.gtf";
        $paths->{gbk} = qq"$paths->{gbk_dir}/$options->{species}.gbk";
        $paths->{gbk_shell} = qq"${libdir_prefix}/genbank/$options->{species}.gbk";
    }
    return($paths);
}


=head2 C<Get_Paths>

Define a set of useful paths for the calling function.  My medium-term goal is to modify
this and Get_Outputs so that the function-specific paths are supplied by a config or
json file.  Before I do that, though, I want to clean out all of the miscellaneous places where
I arbitrarily set paths/filenames.  Thus, for now, these function are serving as a warehouse
for that information so that I may remove all those bits and pieces of junk before creating the
final and hopefully better way of storing this information.

=cut
sub Get_Paths {
    my ($options, %args) = @_;
    my $subroutine = $args{subroutine};
    $subroutine = Bio::Adventure::Config::Get_Caller() unless($subroutine);
    $subroutine = 'undef' if (!defined($subroutine));
    my $libtype = 'genome';
    $libtype = $options->{libtype} if (defined($options->{libtype}));
    my $libpath_prefix = qq"$options->{libpath}/${libtype}";
    my $libdir_prefix = qq"$options->{libdir}/${libtype}";
    my $output_prefix = qq"$options->{output_base}/$options->{jprefix}";
    my $cwd_name = basename(fastcwd());
    my $paths = {};
    my $outputs = {};
    my @input_array = ();
    my @output_array = ();
    my @output_base_array = ();
    if (defined($options->{input})) {
        @input_array = split(/$options->{delimiter}/, $options->{input});
        $paths->{input_suffix} = $options->{input};
        $paths->{input_suffix} =~ s/\.[xgb]z$//g;
        if ($paths->{input_suffix} =~ /\..*$/) {
            $paths->{input_suffix} =~ s/^.*\.(.*)$/$1/g;
        }
        ##$paths->{output_base} = basename($options->{input}, @{$options->{suffix_array}});
        ##$paths->{output_base} = basename($paths->{output_base}, @{$options->{suffix_array}});
        $paths->{input_dirname} = dirname($options->{input});
        if ($paths->{input_dirname} eq '.' || $paths->{input_dirname} eq q{}) {
            $paths->{input_dirname} = fastcwd();
        }
    } else {
        $paths->{input_dirname} = fastcwd();
        @input_array = ($options->{input});
    }

    local *basebase = sub {
        my $in = shift;
        if (!defined($in)) {
            warn("Expected a file path but received undefined.\n");
            return('undefined');
        }
        my $out = basename($in, @{$options->{suffix_array}});
        $out = basename($out, @{$options->{suffix_array}});
        return($out);
    };

    my @input_base_array = ();
    if (@input_array) {
        @input_base_array = map { basebase($_) } @input_array;
    } else {
        print "The input array was not defined in Config::Get_Paths.\n";
        @input_array = ('undefined');
    }

    $paths->{output_base} = $input_base_array[0];
    $paths->{input_array} = \@input_array;
    $paths->{input_base_array} = \@input_base_array;
    $paths->{cwd_name} = $cwd_name;
    $paths->{nosuffix} = qq"$paths->{input_dirname}/$paths->{output_base}";
    $paths->{time} = qq"$paths->{nosuffix}.time";
    $paths->{blast_dir} = qq"${libpath_prefix}/blast";
    $paths->{codon_dir} = qq"${libpath_prefix}/codon_tables";
    $paths->{fasta_dir} = qq"${libpath_prefix}/fasta";
    $paths->{gbk_dir} = qq"${libpath_prefix}/genbank";
    $paths->{gff_dir} = qq"${libpath_prefix}/gff";
    $paths->{gtf_dir} = qq"${libpath_prefix}/gtf";
    $paths->{hmm_dir} = qq"${libpath_prefix}/hmm";
    $paths->{index_prefix} = qq"${libpath_prefix}/indexes";
    $paths->{index_prefix_shell} = qq"${libdir_prefix}/indexes";
    $paths->{output_prefix} = $output_prefix;
    $paths->{output_dir} = qq"${output_prefix}$paths->{output_base}";  ## In case output_dir is not provided.
    $paths->{jname} = $paths->{output_base};
    if (!defined($options->{jprefix})) {
        $options->{jprefix} = '01';
    }
    $options->{jname} = $paths->{jname};
    if ($options->{jprefix}) {
        $options->{jname} = qq"$options->{jprefix}$options->{jname_prefix}";
    }
    if ($options->{libtype}) {
        $options->{jname} .= "_$options->{libtype}";
    }
    make_path($paths->{blast_dir}) unless (-d $paths->{blast_dir});
    make_path($paths->{codon_dir}) unless (-d $paths->{codon_dir});
    make_path($paths->{fasta_dir}) unless (-d $paths->{fasta_dir});
    make_path($paths->{gff_dir}) unless (-d $paths->{gff_dir});
    make_path($paths->{gtf_dir}) unless (-d $paths->{gtf_dir});
    make_path($paths->{gbk_dir}) unless (-d $paths->{gbk_dir});
    make_path($paths->{hmm_dir}) unless (-d $paths->{hmm_dir});
    make_path($paths->{index_prefix}) unless (-d $paths->{index_prefix});
    if (defined($ENV{BLAST_DB})) {
        if ($ENV{BLAST_DB} =~ /:/) {
            my @blast = split(/:/, $ENV{BLAST_DB});
            unless (grep(/$paths->{blast_dir}/, @blast)) {
                $ENV{BLAST_DB} = qq"$ENV{BLAST_DB}:$paths->{blast_dir}";
            }
        } else {
            unless ($paths->{blast_dir} eq $ENV{BLAST_DB}) {
                $ENV{BLAST_DB} = qq"$ENV{BLAST_DB}:$paths->{blast_dir}";
            }
        }
    } else {
        $ENV{BLAST_DB} = $paths->{blast_dir};
    }
    if (defined($options->{species})) {
        $paths->{fasta} = qq"$paths->{fasta_dir}/$options->{species}.fasta";
        $paths->{fasta_shell} = qq"${libdir_prefix}/fasta/$options->{species}.fasta";
        $paths->{gff} = qq"$paths->{gff_dir}/$options->{species}.gff";
        $paths->{gff_shell} = qq"${libdir_prefix}/gff/$options->{species}.gff";
        $paths->{gtf} = qq"$paths->{gtf_dir}/$options->{species}.gtf";
        $paths->{gtf_shell} = qq"${libdir_prefix}/gtf/$options->{species}.gtf";
        $paths->{gbk} = qq"$paths->{gbk_dir}/$options->{species}.gbk";
        $paths->{gbk_shell} = qq"${libdir_prefix}/genbank/$options->{species}.gbk";
    }

    print "At end of primary Get_Paths\n";
    ## Adding this so I can write down examples above for the various slots in $paths.
    print Dumper $paths;

    ## Fill this in instead of this stupid if elsif shenanigans.
    ##my $subroutine_paths = {
    ##    'Bowtie' => {
    ##        index => qq"$paths->{index_prefix}/bt1/$options->{species}",
    ##        index_shell => qq"$paths->{index_prefix_shell}/bt1/$options->{species}",
    ##        index_file => qq"$paths->{index}.1.ebwt",
    ##        index_file_shell =>  qq"$paths->{index_shell}.1.ebwt",
    ##        output_dir => qq"${output_prefix}bowtie_$options->{species}",
    ##    },
    ##    'Bowtie2' => {
    ##    },
    ##};
    given ($subroutine) {
        when ('Any2Any') {
            ## Not yet
        }
        when ('Bam2Bigwig') {
            $paths->{output_dir} = $paths->{input_dirname};
        }
        when ('Bedtools_Coverage') {
            ## Not currently using anything
        }
        when ($_ eq 'Bowtie' || $_ eq 'BT1_Index') {
            $paths->{index} = qq"$paths->{index_prefix}/bowtie/$options->{species}";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/bowtie/$options->{species}";
            $paths->{index_file} = qq"$paths->{index}.1.ebwt";
            $paths->{index_file_shell} = qq"$paths->{index_shell}.1.ebwt";
            $paths->{output_dir} = qq"${output_prefix}bowtie_$options->{species}";
        }
        when ('BT1_Stats') {
            $paths->{output_dir} = $paths->{input_dirname};
        }
        when ($_ eq 'Bowtie2' || $_ eq 'BT2_Index') {
            $paths->{index} = qq"$paths->{index_prefix}/bt2/$options->{species}";
            $paths->{index_dir} = qq"$paths->{index_prefix}/bt2";
            $paths->{index_file} = qq"$paths->{index}.1.bt2";
            $paths->{index_file2} = qq"$paths->{index}.1.bt2l";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/bt2/$options->{species}";
            $paths->{index_file_shell} = qq"$paths->{index_shell}.1.bt2";
            $paths->{index_file2_shell} = qq"$paths->{index_shell}.1.bt2l";
            $paths->{output_dir} = qq"${output_prefix}bowtie2_$options->{species}";
        }
        when ($_ eq 'BWA' || $_ eq 'BWA_Index') {
            $paths->{index} = qq"$paths->{index_prefix}/bwa/$options->{species}";
            $paths->{index_dir} = qq"$paths->{index_prefix}/bwa";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/bwa/$options->{species}";
            $paths->{index_file} = qq"$paths->{index}.sa";
            $paths->{index_file_shell} = qq"$paths->{index_shell}.sa";
            $paths->{output_dir} = qq"${output_prefix}bwa_$options->{species}";
        }
        when (/^Caical/) {
            $paths->{output_dir} = qq"${output_prefix}caical";
        }
        when ('Check_Blastdb') {
            ## Not yet
        }
        when (/Classify_Phage/) {
            $paths->{output_dir} = qq"${output_prefix}classify_phage";
            $paths->{log} = qq"$paths->{output_dir}/classify_phage.log";
        }
        when ($_ eq 'Compress' || $_ eq 'Recompress') {
            $paths->{output_dir} = qq"${output_prefix}xz";
        }
        ## The _Worker is wrong
        when (/^Concatenate_Searches/) {
            $paths->{output_dir} = $paths->{input_dirname};
        }
        when (/^Dantools/) {
            $paths->{output_dir} = qq"${output_prefix}dantools_rna";
            $paths->{log} = qq"$paths->{output_dir}/dantools_rnaseq.log";
        }
        when ('Download_NCBI_Assembly') {
            $paths->{output_dir} = qq"${output_prefix}download_assembly";
            $paths->{log} = qq"$paths->{output_dir}/download_assembly.log";
        }
        when ('Downsample_Guess_Strand') {
            $paths->{index} = qq"$paths->{index_prefix}/salmon/$options->{species}";
            $paths->{index_dir} = qq"$paths->{index_prefix}/salmon";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/salmon/$options->{species}";
            $paths->{output_dir} = qq"${output_prefix}salmon_downsample_$options->{species}";
        }
        when ('Duplicate_Remove') {
            ## Not yet
        }
        when ('Extend_Kraken_DB') {
            $paths->{output_dir} = qq"${output_prefix}extend_kraken";
        }
        when ('Extract_Subseq') {
            $paths->{output_dir} = qq"${output_prefix}subseq";
        }
        when ('Fastqc') {
            $paths->{output_dir} = qq"${output_prefix}fastqc";
            $paths->{output_subdir} = qq"$paths->{output_dir}/$paths->{output_base}_fastqc";
        }
        when ('Fastqc_Stats') {
            $paths->{output_dir} = $paths->{input_dirname};
        }
        when ('Feature_Counts') {
            $paths->{output_dir} = qq"${output_prefix}featureCounts";
        }
        when (/^Filter_Host_Kraken/) {
            $paths->{output_dir} = qq"${output_prefix}filter_kraken";
            $paths->{index_dir} = qq"$paths->{index_prefix}/hisat";
        }
        when ($_ eq 'Freebayes_SNP_Search' || $_ eq 'SNP_Ratio' || $_ eq 'SNP_Ratio_Worker' ||
              $_ eq 'SNP_Ratio_Intron' || $_ eq 'SNP_Ratio_Intron_Worker') {
            $paths->{output_dir} = qq"${output_prefix}freebayes_$options->{species}";
        }
        when ('GATK_Dedup') {
            $paths->{output_dir} = qq"${output_prefix}gatk_dedup";
        }
        when (/^Gb2Gff/) {
            $paths->{output_dir} = qq"${output_prefix}gb2gff";
        }
        when ('Get_Substitution_Matrix') {
            $paths->{output_dir} = qq"${output_prefix}matrices";
        }
        when ('GFF_Tree') {
            $paths->{output_dir} = qq"${output_prefix}gff_tree";
        }
        when ('Gff2Fasta') {
            $paths->{output_dir} = qq"${output_prefix}gff2fasta";
        }
        when (/^Guess_Strand/) {
            ## Not yet
        }
        when ($_ eq 'Hisat2' || $_ eq 'Hisat2_Stats') {
            $paths->{index_dir} = qq"$paths->{index_prefix}/hisat";
            $paths->{index} = qq"$paths->{index_dir}/$options->{species}";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/hisat/$options->{species}";
            $paths->{index_file} = qq"$paths->{index}.1.ht2";
            $paths->{index_file_shell} = qq"$paths->{index_shell}.1.ht2";
            $paths->{index_file2} = qq"$paths->{index}.1.ht2l";
            $paths->{index_file2_shell} = qq"$paths->{index_shell}.1.ht2l";
            $paths->{output_dir} = qq"${output_prefix}hisat_$options->{species}";
            $paths->{index_stderr} = qq"${output_prefix}/hisat2_index_$options->{species}.stderr";
            $paths->{index_stdout} = qq"${output_prefix}/hisat2_index_$options->{species}.stdout";
        }
        when ('Hisat2_Index') {
            if (!defined($options->{species})) {
                $options->{species} = basename($options->{input}, ('.fsa', '.fasta', '.fa'));
                $paths->{fasta} = qq"$paths->{fasta_dir}/$options->{species}.fasta";
                $paths->{fasta_shell} = qq"${libdir_prefix}/fasta/$options->{species}.fasta";
                $paths->{gff} = qq"$paths->{gff_dir}/$options->{species}.gff";
                $paths->{gff_shell} = qq"${libdir_prefix}/gff/$options->{species}.gff";
            }
            $paths->{index_dir} = qq"$paths->{index_prefix}/hisat";
            $paths->{index} = qq"$paths->{index_dir}/$options->{species}";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/hisat/$options->{species}";
            $paths->{index_file} = qq"$paths->{index}.1.ht2";
            $paths->{index_file_shell} = qq"$paths->{index_shell}.1.ht2";
            $paths->{index_file2} = qq"$paths->{index}.1.ht2l";
            $paths->{index_file2_shell} = qq"$paths->{index_shell}.1.ht2l";
            $paths->{output_dir} = qq"${output_prefix}index_hisat2_$options->{species}";
            make_path($paths->{output_dir}, {verbose => 0}) unless (-r $paths->{output_dir});
        }
        when ('Hisat2_Stats') {
            $paths->{output_dir} = $paths->{input_dirname};
        }
        when ($_ eq 'HT_Multi' || $_ eq 'HTSeq' || $_ eq 'HT_Types') {
            ## Not yet
        }
        when ('Insert_Size') {
            $paths->{output_dir} = dirname($options->{input});
        }
        when ($_ eq 'Jellyfish' || $_ eq 'Jellyfish_Matrix') {
            $paths->{output_dir} = qq"${output_prefix}jellyfish_${cwd_name}";
            $paths->{jelly_base} = qq"$paths->{output_dir}/$ {cwd_name}_$options->{length}";
        }
        when ($_ eq 'Kallisto' || $_ eq 'Kallisto_Index') {
            $paths->{index_dir} = qq"$paths->{index_prefix}/kallisto";
            $paths->{index_file} = qq"$paths->{index_prefix}/kallisto/$options->{species}.idx";
            $paths->{index_file_shell} = qq"$paths->{index_prefix_shell}/kallisto/$options->{species}.idx";
            $paths->{output_dir} = qq"${output_prefix}kallisto_$options->{species}";
            $paths->{index_dir} = qq"${output_prefix}kallisto_index_$options->{species}";
        }
        when (/^Kraken/) {
            $paths->{output_dir} = qq"${output_prefix}kraken_$options->{library}";
        }
        when ('Make_Blast_Job') {
            ## Not yet
        }
        when ('Make_Codon_Table') {
            ## Not yet
        }
        when ('Mash') {
            ## Not yet
        }
        when ('Mauve_Aligner') {
            $paths->{output_dir} = qq"${output_prefix}/mauve";
        }
        when (/^Merge_Annotations/) {
            $paths->{output_dir} = qq"${output_prefix}mergeannot";
        }
        when ($_ eq 'Mi_Map' || $_ eq 'Read_Bam_Mi' || %_ eq 'Final_Print_Mi') {
            ## Not yet
        }
        when ('Mpileup') {
            $paths->{output_dir} = qq"${output_prefix}mpileup_$options->{species}";
        }
        when ('Nextflow_RNASeq') {
            $paths->{output_dir} = qq"${output_prefix}nf_rna_$options->{species}";
        }
        when (/^OrthoFinder/) {
            $paths->{output_dir} = qq"${output_prefix}orthofinder";
            $paths->{fasta_dir} = qq"$paths->{output_dir}/input";
        }
        when ('Pairwise_Similarity_Matrix') {
            $paths->{output_dir} = qq"${output_prefix}pairwise_blast";
        }
        when ('PolyA_Recorder') {
            $paths->{output_dir} = qq"${output_prefix}polya_recorder_$options->{species}";
        }
        when ('Prodigal') {
            my $input_name = basename($options->{input}, ('.fasta', '.fsa', '.ffn'));
            $paths->{output_dir} = qq"${output_prefix}prodigal_${input_name}";
        }
        when ('Progressive_Mauve') {
            $paths->{output_dir} = qq"${output_prefix}/pmauve";
        }
        when (%_ eq 'ProteinFold' || $_ eq 'ProteinFold_Worker') {
            my $output_suffix = basename($options->{input}, ('.fsa', '.faa', '.fasta', '.fa', '.ffn'));
            $paths->{output_dir} = qq"${output_prefix}proteinfold_${output_suffix}";
        }
        when ($_ eq 'ProteinFold_PairIDs' || $_ eq 'ProteinFold_PairIDs_Worker') {
            my ($in1, $in2) = split(/$options->{delimiter}/, $options->{input});
            my $output_suffix_first = basename($in1, ('.txt', '.fsa', '.faa', '.fasta', '.fa', '.ffn'));
            my $output_suffix_second = basename($in2, ('.txt', '.fsa', '.faa', '.fasta', '.fa', '.ffn'));
            $paths->{output_dir} = qq"${output_prefix}${output_suffix_first}";
            $paths->{cif_output_dir} = qq"${output_prefix}${output_suffix_first}/${output_suffix_second}";
            $paths->{suffix_first} = $output_suffix_first;
            $paths->{suffix_second} = $output_suffix_second;
        }
        when (/^RNAFold_Windows/) {
            $paths->{output_dir} = qq"${output_prefix}rnafold";
        }
        when (/^RSEM/) {
            $paths->{index} = qq"$paths->{index_prefix}/rsem/$options->{species}";
            $paths->{index_dir} = qq"$paths->{index_prefix}/rsem";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/rsem/$options->{species}";
            $paths->{index_file} = qq"$paths->{index}.transcripts.fa";
            $paths->{index_file_shell} = qq"$paths->{index_shell}.transcripts.fa";
            $paths->{output_dir} = qq"${output_prefix}rsem_$options->{species}";
        }
        when ('Run_Parse_Blast') {
            my $libname = basename($options->{library}, $options->{suffix_array});
            $paths->{output_dir} = qq"${output_prefix}blast_$options->{input}_vs_${libname}";
        }
        when ($_ eq 'Salmon' || $_ eq 'Salmon_Index') {
            $paths->{index} = qq"$paths->{index_prefix}/salmon/$options->{species}";
            $paths->{index_dir} = qq"$paths->{index_prefix}/salmon/$options->{species}";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/salmon/$options->{species}";
            $paths->{output_dir} = qq"${output_prefix}salmon_$options->{species}_$options->{libtype}";
        }
        when ($_ eq 'Sam2Bam' || $_ eq 'Samtools') {
            $paths->{output_dir} = qq"$paths->{input_dirname}";
            $paths->{output_dir} = qq"$paths->{input_dirname}";
        }
        when ('SLSearch') {
            $paths->{output_dir} = qq"${output_prefix}slsearch_$options->{species}";
        }
        when ('SL_Recorder') {
            $paths->{output_dir} = qq"${output_prefix}sl_recorder_$options->{species}";
        }
        when ('Split_Align_Blast') {
            my $libname = basename($options->{library}, $options->{suffix_array});
            $paths->{output_dir} = qq"${output_prefix}blastsplit_$paths->{output_base}_vs_${libname}";
        }
        when ($_ eq 'Split_Align_Fasta' || $_ eq 'Split_Align_Fasta_Worker' ||
              $_ eq 'Parse_Fasta' || $_ eq 'Parse_Fasta_Global' ||
              $_ eq 'Make_Fasta_Job' || $_ eq 'Parse_Fasta_Worker' ||
              $_ eq 'Parse_Fasta_Mismatches') {
            my $libname = basename($options->{library}, @{$options->{suffix_array}});
            $libname = basename($libname, @{$options->{suffix_array}});
            my $query = basename($options->{input}, @{$options->{suffix_array}});
            $query = basename($query, @{$options->{suffix_array}});
            if ($subroutine =~ /Parse/) {
                ## The input for this is already x_vs_y
                $query =~ s/^(.*)_vs_.*$/$1/g;
            }
            $options->{output_base} = qq"outputs";
            $output_prefix = qq"$options->{output_base}/$options->{jprefix}";
            $paths->{output_dir} = qq"${output_prefix}fasta_${query}_${libname}";
            ## Note, basing the output filename on a variable is dangerous
            $paths->{output_filename_base} = qq"${query}_vs_${libname}";
        }
        when ('Split_Fasta') {
            $paths->{output_dir} = qq"${output_prefix}split";
        }
        when ('Spring') {
            my $input_name = basename($options->{input}, ('.gz', '.bz2', '.xz'));
            $input_name = basename($input_name, ('.fastq'));
            $paths->{output_dir} = qq"${output_prefix}spring";
        }
        when (/^STAR/) {
            $paths->{index} = qq"$paths->{index_prefix}/star/$options->{species}";
            $paths->{index_dir} = qq"$paths->{index_prefix}";
            $paths->{index_file} = qq"$paths->{index}/SAindex";
            $paths->{index_file_shell} = qq"$paths->{index_shell}/SAindex";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/star/$options->{species}";
            $paths->{output_dir} = qq"${output_prefix}star_$options->{species}";
        }
        when (/^Terminase_ORF_Reorder/) {
            $paths->{output_dir} = qq"${output_prefix}terminase_reorder";
        }
        when ('Test_Worker') {
            $paths->{funkytown} = qq"$paths->{index_prefix}/funkytown";
        }
        when ('Tophat') {
            $paths->{index} = qq"$paths->{index_prefix}/tophat/$options->{species}";
            $paths->{index_dir} = qq"$paths->{index_prefix}/tophat";
            $paths->{index_file} = qq"$paths->{index}.1.bt2";
            $paths->{index_file_shell} = qq"$paths->{index_shell}.1.bt2";
            $paths->{index_shell} = qq"$paths->{index_prefix_shell}/tophat/$options->{species}";
            $paths->{output_dir} = qq"${output_prefix}tophat_$options->{species}";
        }
        when ('Transit_TPP') {
            $paths->{output_dir} = qq"${output_prefix}transit_$options->{species}";
        }
        when ('Trimomatic_Pairwise') {
            $paths->{output_dir} = qq"${output_prefix}trimomatic";
        }
        when ('Trimomatic_Single') {
            $paths->{output_dir} = qq"${output_prefix}trimomatic";
        }
        when ('Trimomatic_Stats') {
            $paths->{output_dir} = $paths->{input_dirname};
        }
        when ('Umi_Tools_Extract') {
            $paths->{output_dir} = qq"${output_prefix}umi_tools";
        }
        when ('Umi_Tools_Dedup') {
            $paths->{output_dir} = qq"${output_prefix}umi_dedup";
        }
        when ('Unicycler') {
            $paths->{output_dir} = qq"${output_prefix}unicycler";
        }
        when ('Velvet') {
            $paths->{output_dir} = qq"${output_prefix}velvet";
        }
        default {
            $paths->{output_dir} = $paths->{input_dirname};
        }
    }  ## End the given
    if ($paths->{index_file}) {
        my $index_directory = dirname($paths->{index_file});
        make_path($index_directory) unless (-d $index_directory);
    }
    if (!defined($paths->{log})) {
        $paths->{log} = qq"$paths->{output_dir}/log.txt";
    }
    ## Create the output directory and set a default stdout/stderr file.
    make_path($paths->{output_dir}) unless (-d $paths->{output_dir});
    return($paths);
}

=head2 C<Get_Spec>

Define the parameters on a computing cluster required to submit a job and expect it to succeed.
This is where I would really like to start using the heuristics for memory/time consumption.

=cut
sub Get_Spec {
    my ($class, %args) = @_;
    my $function = Bio::Adventure::Config::Get_Caller();
    ## I think that with recent changes, one need no longer add 'cyoa' as a module.
    my %spec_data = (
        'Default' => {
            jcomment => 'Default placeholder comment.',
            jcpu => 1,
            jgpu => 0,
            jmem => 4,
            jname_prefix => 'default',
            jprefix => '99',
            jwalltime => '2:00:00',
        },
        'Test' => {
            jcomment => 'Comment for a test job.',
            jcpu => 1,
            jgpu => 0,
            jmem => 12,
            jname_prefix => 'test',
            jprefix => '99',
            jwalltime => '2:00:00',
        },
        'Abyss' => {
            jcomment => 'Assemble a genome via makefile shenanigans.',
            jcpu => 6,
            jmem => 12,
            jname_prefix => 'abyss',
            jprefix => '15',
            jwalltime => '4:00:00',
        },
        'Align_SNP_Search' => {
            jcomment => 'Extract variants from short reads.',
            jcpu => 4,
            jmem => 8,
            jname_prefix => 'snpsearch',
            jprefix => '80',
            jwalltime => '8:00:00',
        },
        'Aragorn' => {
            jcomment => 'Seek (m)tRNA genes with Aragorn.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'aragorn',
            jprefix => '21',
            jwalltime => '1:00:00',
        },
        'Assembly_Coverage' => {
            jcomment => 'Use hisat2 and bbmap to create a coverage map.',
            jcpu => 6,
            jmem => 18,
            jname_prefix => 'assemblecov',
            jprefix => '18',
            jwalltime => '8:00:00',
        },
        'Bacphlip' => {
            jcomment => 'ML method to seek lytic/lysogenic phage sequences.',
            jcpu => 2,
            jmem => 12,
            jname_prefix => 'bacphlip',
            jprefix => '80',
            jwalltime => '1:00:00',
        },
        'Bam2Bigwig' => {
            jcomment => 'Convert an alignment to bigwig with deeptools.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'bam2bw',
            jprefix => '60',
            jwalltime => '01:00:00',
        },
        'Bedtools_Coverage' => {
            jcomment => 'Use bedtools to create a coverage map.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'bedtools',
            jprefix => '60',
            jwalltime => '2:00:00',
        },
        'Bowtie' => {
            jcomment => 'Fast short-read aligner.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'bowtie',
            jprefix => '10',
            jwalltime => '4:00:00',
        },
        'Bowtie_RRNA' => {
            jcomment => 'Quantify rRNA reads in a short-read dataset.',
            jcpu => 2,
            jmem => 8,
            jname_prefix => 'bt_rrna',
            jprefix => '11',
            jwalltime => '3:00:00',
        },
        'Bowtie2' => {
            jcomment => 'Fast short-read aligner, version 2.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'bowtie2',
            jprefix => '20',
            jwalltime => '8:00:00',
        },
        'BT1_Index' => {
            jcomment => 'Index a genome using bowtie1.',
            jcpu => 2,
            jmem => 20,
            jname_prefix => 'bt1_index',
            jprefix => '10',
            jwalltime => '2:00:00',
        },
        'BT1_Stats' => {
            jcomment => 'Collect mapping stats for bowtie1.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'bt1_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'BT2_Index' => {
            jcomment => 'Index a genome using bowtie2.',
            jcpu => 2,
            jmem => 20,
            jname_prefix => 'bt2_index',
            jprefix => '10',
            jwalltime => '2:00:00',
        },
        'BT2_Stats' => {
            jcomment => 'Collect mapping stats for bowtie2.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'bt2_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'BWA' => {
            jcomment => 'Use the Burroughs-Wheeler aligner to quickly align short reads.',
            jcpu => 4,
            jmem => 40,
            jname_prefix => 'bwa',
            jprefix => '30',
            jwalltime => '72:00:00',
        },
        'BWA_Index' => {
            jcomment => 'Index a genome using bwa.',
            jcpu => 2,
            jmem => '12',
            jname_prefix => 'bwa_index',
            jprefix => '10',
            jwalltime => '2:00:00',
        },
        'BWA_Stats' => {
            jcomment => 'Collect mapping stats for bwa.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'bwa_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'Caical' => {
            jcomment => 'Calculate the codon adaptation index of a set of sequences.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'caical',
            jprefix => '80',
            jwalltime => '2:00:00',
        },
        'Caical_CDS' => {
            jcomment => 'Calculate the codon adaptation index of a set of sequences.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'caical',
            jprefix => '80',
            jwalltime => '2:00:00',
        },
        'Casfinder' => {
            jcomment => 'Seek Cas-like genes in bacterial genomes.',
            jcpu => 4,
            jmem => 8,
            jname_prefix => 'casfinder',
            jprefix => '22',
            jwalltime => '1:00:00',
        },
        'Casoff' => {
            jcomment => 'An incomplete, separate Cas-search.',
            jcpu => 4,
            jmem => 8,
            jname_prefix => 'casoff',
            jprefix => '22',
            jwalltime => '2:00:00',
        },
        'CGView' => {
            jcomment => 'Run cgview on an assembly.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'cgview',
            jwalltime => '2:00:00',
        },
        'Classify_Phage' => {
            jcomment => 'Viral classification via tblastx.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'phageclass',
            jprefix => '18',
            jwalltime => '2:00:00',
        },
        'Cogent' => {
            jcomment => 'Use the cogent pipeline to handle umi data.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'cogent',
            jprefix => '01',
            jwalltime => '36:00:00',
        },
        'Collect_Assembly' => {
            jcomment => 'Collect all the outputs of an assembly pipeline into 1 place.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'collect',
            jprefix => '81',
            jwalltime => '0:10:00',
        },
        'Compress' => {
            jcomment => 'Compress an input using xz.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'xz',
            jprefix => $args{jprefix},
            jwalltime => '12:00:00',
        },
        'Concatenate_Searches' => {
            jcomment => 'Concatenate BLAST output files.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'concat',
            jprefix => $args{jprefix},
            jwalltime => '1:00:00',
        },
        'Consolidate_TAs' => {
            jcomment => 'Consolidate mariner sequences from a TNSeq data set.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'consolidateta',
            jprefix => '59',
            jwalltime => '4:00:00',
        },
        'Contig_Reoder' => {
            jcomment => 'Use mauve to rearrange contigs',
            jcpu => 2,
            jmem => 24,
            jname_prefix => 'contigreorder',
            jprefix => '19',
            jwalltime => '1:00:00',
        },
        'Cutadapt' => {
            jcomment => 'Trim sequence with cutadapt.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'cutadapt',
            jprefix => '01',
            jwalltime => '30:00:00',
        },
        'Cutadapt_Stats' => {
            jcomment => 'Collect mapping stats from cutadapt.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'cutadapt_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'Dantools_RNASeq' => {
            jcomment => 'Invoke dantools on rnaseq data.',
            jcpu => 8,
            jmem => 24,
            jname_prefix => 'dantools',
            jprefix => '22',
            jwalltime => '64:00:00',
        },
        'Downsample_Guess_Strand' => {
            jcomment => 'Use salmon on the first few reads to guess the library type.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'salmonguess',
            jprefix => '45',
            jwalltime => '0:20:00',
        },
        'Essentiality_TAs' => {
            jcomment => 'Extract TA sequences from a TNSeq dataset.',
            jcpu => 2,
            jmem => 8,
            jname_prefix => 'essentiality_ta',
            jprefix => '40',
            jwalltime => '12:00:00',
        },
        'Extend_Kraken_DB' => {
            jcomment => 'Update an existing kraken2 database.',
            jcpu => 6,
            jmem => 64,
            jname_prefix => 'kraken_extend',
            jprefix => '09',
            jwalltime => '36:00:00',
        },
        'Fastqc' => {
            jcomment => 'Estimate sequence quality with fastqc.',
            jcpu => 4,
            jmem => 6,
            jname_prefix => 'fastqc',
            jprefix => '02',
            jwalltime => '03:00:00',
        },
        'Fastqc_Stats' => {
            jcomment => 'Collect results from fastqc.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'fastqc_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'Fastp' => {
        },
        'Fastp_Stats' => {
            jcomment => 'Collect mapping stats from fastp.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'fastp_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'Feature_Counts' => {
            jcomment => 'Quantify reads from an alignment with featureCounts.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'fcounts',
            jprefix => '44',
            jwalltime => '4:00:00',
        },
        'Filter_Host_Kraken' => {
            jcomment => 'Use kraken to find the likely host of a viral sequence and remove them.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'filter_host',
            jprefix => '06',
            jwalltime => '2:00:00',
        },
        'Freebayes_SNP_Search' => {
            jcomment => 'Seek variants with freebayes and write out the results.',
            jcpu => 4,
            jmem => 36,
            jname_prefix => 'freebayes',
            jprefix => '60',
            jwalltime => '24:00:00',
        },
        'GATK_Dedup' => {
            jcomment => 'Deduplicate reads via alignment positions of R1/R2.',
            jcpu => 4,
            jmem => 36,
            jname_prefix => 'gatk_dedup',
            jprefix => '55',
            jwalltime => '40:00:00',
        },
        'Gb2Gff' => {
            jcomment => 'Extract the various datatypes from a genbank file.',
            jcpu => 1,
            jmem => 6,
            jname_prefix => 'gb2gff',
            jprefix => '78',
            jwalltime => '01:00:00',
        },
        'Generate_Samplesheet' => {
            jcomment => 'Create a samplesheet from the logs produced by all these tools',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'generate_metadata',
            jprefix => '99',
            jwalltime => '2:00:00',
        },
        'Get_Split' => {
            jcomment => 'Split a fasta into pieces for searching.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'get_split',
            jprefix => $args{jprefix},
            jwalltime => '00:10:00',
        },
        'Gff2Fasta' => {  ## FIXME: Dispatch this to a worker!
            jcomment => 'Extract sequences of arbitrary types from a genome and gff file.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'gff2fasta',
            jprefix => '78',
            jwalltime => '01:00:00',
        },
        'Glimmer' => {
            jcomment => 'Seek prokaryotic genes with the OG ORF search method.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'glimmer',
            jprefix => '16',
            jwalltime => '1:00:00',
        },
        'Glimmer_Single' => {
            jcomment => 'Seek prokaryotic genes with the OG ORF search method, without training.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'glimmer',
            jprefix => '16',
            jwalltime => '1:00:00',
        },
        'Gubbins' => {
            jcomment => 'Given a fasta and outgroup, estimate their tree.',
            jcpu => 8,
            jmem => 24,
            jname_prefix => 'gubbins',
            jprefix => '65',
            jwalltime => '24:00:00',
        },
        'Guess_Strand' => {
            jcomment => 'Examine a bamfile to gather strandedness information.',
            jcpu => 1,
            jmem => 12,
            jname_prefix => 'guess_strand',
            jprefix => '46',
            jwalltime => '8:00:00',
        },
        'Hisat2' => {
            jcomment => 'Another iteration of bowtie, now with a new name!',
            jcpu => 8,
            jmem => 48,
            jname_prefix => 'hisat',
            jprefix => '40',
            jwalltime => '36:00:00',
        },
        'Hisat2_Index' => {
            jcomment => 'Index a genome using hisat2.',
            jcpu => 4,
            jmem => 20,
            jname_prefix => 'index_hisat',
            jprefix => '09',
            jwalltime => '4:00:00',
        },
        'Hisat2_Stats' => {
            jcomment => 'Collect mapping stats from hisat.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'hisat_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'HTSeq' => {
            jcomment => 'Quantify observed reads/feature with htseq.',
            jcpu => 1,
            jmem => 30,
            jname_prefix => 'htseq',
            jprefix => '42',
            jwalltime => '12:00:00',
        },
        'HT_Types' => {
            ## This does not currently dispatch to a worker.
            jcomment => 'Check the first lines of a gff file to find the most common tags/types.',
            jcpu => 1,
            jmem => 3,
            jname_prefix => 'ht_types',
            jprefix => '89',
            jwalltime => '00:02:00',
        },
        'Insert_Size' => {
            jcomment => 'Estimate the insert size of a library by looking at aligned paired reads.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'insert_size',
            jprefix => '44',
            jwalltime => '04:00:00',
        },
        'Interproscan' => {
            jcomment => 'Use the a large suite of search tools to seek similarities to known elements.',
            jcpu => 6,
            jmem => 24,
            jname_prefix => 'interpro',
            jprefix => '21',
            jwalltime => '36:00:00',
        },
        'Interpro_Long2Wide' => {
            jcomment => 'Rewrite interproscan outputs to a tabular format.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'interl2w',
            jprefix => '21',
            jwalltime => '1:00:00',
        },
        'Jellyfish' => {
            jcomment => 'Count kmers with jellyfish.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'jellyfish',
            jprefix => '17',
            jwalltime => '02:00:00',
        },
        'Jellyfish_Matrix' => {
            jcomment => 'Create a matrix with rows as kmers and elements as number of occurrences.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'jelly_mtrx',
            jprefix => '18',
            jwalltime => '01:00:00',
        },
        'Kallisto' => {
            jcomment => 'Fast quantitation method with the power of deBruijn.',
            jcpu => 8,
            jmem => 24,
            jname_pref => 'kallisto',
            jprefix => '46',
            jwalltime => '8:00:00',
        },
        'Kallisto_Index' => {
            jcomment => 'Index a genome using kallisto.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'kallisto_index',
            jprefix => '09',
            jwalltime => '4:00:00',
        },
        'Kraken' => {
            jcomment => 'Use a kmer index to estimate reads/species.',
            jcpu => 8,
            jmem => 64,
            jname_prefix => 'kraken',
            jprefix => '11',
            jwalltime => '04:00:00',
        },
        'Kraken_Best_Hit' => {
            jcomment => 'Extract the best genus/phylum/etc from a kraken result.',
            jcpu => 8,
            jmem => 64,
            jname_prefix => 'kraken_best',
            jprefix => '11',
            jwalltime => '4:00:00',
        },
        'Kraken_to_Matrix' =>  {
            jcomment => 'Convert kraken output to a tabular format.',
            jcpu => 1,
            jmem => 6,
            jname_prefix => 'kraken_mtrx',
            jprefix => '12',
            jwalltime => '00:30:00',
        },
        'Last_Stat' => {
            ## Currently there is no scenario in which this job will run on a cluster.
            cluster => 0,
            jcomment => 'Get the last line of the stats file for a job.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'last_stat',
            jprefix => '99',
            jwalltime => '0:01:00',
        },
        'Make_Blast_Job' => {
            jcomment => 'Run multiple BLAST jobs',
            jcpu => 8,
            jmem => 24,
            jname_prefix => 'blast',
            jprefix => '89',
            jwalltime => '8:00:00',
        },
        'Make_Codon_Table' => {  ## Note to self, this needs to be dispatched.
            jcomment => 'Create a table of relative codon abundances given a set of CDS.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'make_codon_table',
            jprefix => '08',
            jwalltime => '1:00:00',
        },
        'Make_Directories' => {
            jcomment => 'Create directories for splitting a fasta file.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'align_mkdir',
            jprefix => $args{jprefix},
            jwalltime => '00:10:00',
        },
        'Make_Fasta_Job' => {
            jcomment => 'Run multiple fasta36 jobs',
            jcpu => 8,
            jmem => 22224,
            jname_prefix => 'fasta36',
            jprefix => $args{jprefix},
            jwalltime => '8:00:00',
        },
        'Mash' => {
            jcomment => 'Compute distances between large sets of sequences.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'mash',
            jprefix => '21',
            jwalltime => '02:00:00',
        },
        'Mauve_Aligner' => {
            jcomment => 'Attempt to run mauveAligner using a directory of genbank/fasta files.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'mauve',
            jprefix => '20',
            jwalltime => '4:00:00',
        },
        'Merge_Annotations' => {
            jcomment => 'Gather annotation data from multiple sources.',
            jcpu => 1,
            jmem => 6,
            jname_prefix => 'merge_annot',
            jprefix => '15',
            jwalltime => '00:30:00',
        },
        'Merge_CDS_Predictions' => {
            jcomment => 'Merge CDS predictions from prodigal/glimmer/phanotate/etc.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'merge_cds',
            jprefix => '19',
            jwalltime => '0:30:00',
        },
        'Mpileup' => {
            jcomment => 'Seek counts/variants per position with samtools.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'mpileup',
            jprefix => '55',
            jwalltime => '06:00:00',
        },
        'Mpileup_SNP_Search' => {
            jcomment => 'Use samtools mpileup to seek variants.',
            jcpu => 2,
            jmem => 12,
            jname_prefix => 'mpileup',
            jprefix => '70',
            jwalltime => '36:00:00',
        },
        'Merge_Parse_Blast' => {
            jcomment => 'Merge BLAST outputs and parse the result into a simpler table.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'mergeblast',
            jprefix => 44,
            jwalltime => '1:00:00',
        },
        'OrthoFinder' => {
            jcomment => 'Attempt to run orthofinder using a faa directory.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'orthofinder',
            jprefix => '50',
            jwalltime => '8:00:00',
        },
        'OrthoMCL_Pipeline' => {
            jcomment => 'Use github.com/apetkau/orthomcl-pipeline to run orthomcl.',
            jmem => 8,
            jcpu => 8,
            jname_prefix => 'orthomcl',
            jprefix => '60',
            jwalltime => '36:00:00',
        },
        'Parse_Fasta' => {
            jcomment => 'Parse the output from fasta36.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'parsefasta',
            jprefix => $args{jprefix},
            jwalltime => '01:00:00',
        },
        'Phagepromoter' => {
            jcomment => 'Seek phage promoter sequences in sequences.',
            jcpu => 1,
            jmem => 6,
            jname_prefix => 'phagepromoter',
            jprefix => '30',
            jwalltime => '2:00:00',
        },
        'Phageterm' => {
            jcomment => 'Seek phager terminii in short read sequence.',
            jcpu => 8,
            jmem => 12,
            jname_prefix => 'phageterm',
            jprefix => '14',
            jwalltime => '1:00:00',
        },
        'Phanotate' => {
            jcomment => 'A very aggressive CDS search tuned for phage genes.',
            jcpu => 1,
            jmem => 6,
            jname_prefix => 'phagepromoter',
            jprefix => '17',
            jwalltime => '2:00:00',
        },
        'Phastaf' => {
            jcomment => 'Seek phage-derived sequences in bacterial assemblies.',
            jcpu => 4,
            jmem => 8,
            jname_prefix => 'phastaf',
            jprefix => '15',
            jwalltime => '2:00:00',
        },
        'PolyA_Extractor' => {
            jcomment => 'Extract polyA reads from a dataset.',
            jcpu => 1,
            jmem => 12,
            jname_prefix => 'polyaextract',
            jprefix => '10',
            jwalltime => '36:00:00',
        },
        'PolyA_Recorder' => {
            jcomment => 'Extract and record polyA locations from short reads.',
            jcpu => 1,
            jmem => 5,
            jname_prefix => 'polya_recorder',
            jprefix => '50',
            jwalltime => '4:00:00',
        },
        'Prodigal' => {
            jcomment => 'An easily trained CDS search for bacterial sequences.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'prodigal',
            jprefix => '17',
            jwalltime => '3:00:00',
        },
        'Progressive_Mauve' => {
            jcomment => 'Attempt to run progressive mauve using a directory of genbank/fasta files.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'pmauve',
            jprefix => '20',
            jwalltime => '4:00:00',
        },
        'Prokka' => {
            jcomment => 'Quick bacterial assembly annotation pipeline.',
            jcpu => 8,
            jmem => 12,
            jname_prefix => 'prokka',
            jprefix => '19',
            jwalltime => '3:00:00',
        },
        'ProteinFold' => {
            jcomment => 'Use alphafold3 to fold some sequences.',
            jcpu => 8,
            jgpu => 1,
            jmem => 24,
            jname_prefix => 'alphafold',
            jprefix => '79',
            jwalltime => '8:00:00',
        },
        'ProteinFold_JSON_Pairwise_OneInput' => {
            jcomment => 'Use alphafold3 to fold paired sequences.',
            jcpu => 8,
            jgpu => 1,
            jmem => 24,
            jname_prefix => 'alphafold',
            jprefix => '79',
            jwalltime => '8:00:00',
        },
        'ProteinFold_JSON_Pairwise_TwoSeq' => {
            jcomment => 'Use alphafold3 to fold two sequences.',
            jcpu => 8,
            jgpu => 1,
            jmem => 24,
            jname_prefix => 'alphafold',
            jprefix => '79',
            jwalltime => '8:00:00',
        },
        'ProteinFold_JSON_Pairwise_TwoInput' => {
            jcomment => 'Use alphafold3 to fold two sequences from two sequence databases.',
            jcpu => 8,
            jgpu => 1,
            jmem => 24,
            jname_prefix => 'alphafold',
            jprefix => '79',
            jwalltime => '8:00:00',
        },
        'ProteinFold_JSON_Separate' => {
            jcomment => 'Use alphafold3 to fold separate sequences.',
            jcpu => 8,
            jgpu => 1,
            jmem => 24,
            jname_prefix => 'alphafold',
            jprefix => '79',
            jwalltime => '8:00:00',
        },
        'ProteinFold_JSON_Together' => {
            jcomment => 'Use alphafold3 to fold paired sequences.',
            jcpu => 8,
            jgpu => 1,
            jmem => 24,
            jname_prefix => 'alphafold',
            jprefix => '79',
            jwalltime => '8:00:00',
        },
        'ProteinFold_PairIDs' => {
            jcomment => 'Use alphafold3 to fold pairs of sequences.',
            jcpu => 8,
            jgpu => 1,
            jmem => 36,
            jname_prefix => 'alphafold',
            jprefix => '79',
            jwalltime => '12:00:00',
        },
        'Racer' => {
            jcomment => 'correct reads for sequencing errors with kmers.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'racer',
            jprefix => '10',
            jwalltime => '8:00:00',
        },
        'Read_Genome_Fasta' => {
            ## Currently there is no scenario in which this job will run on a cluster.
            cluster => 0,
        },
        'Recompress' => {
            jcomment => 'Compress some files.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'recomp',
            jprefix => '80',
            jwalltime => '08:00:00',
        },
        'Restriction_Catalog' => {
            jcomment => 'Catalog the set of restriction enzymes in a sequence.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'recatalog',
            jprefix => '18',
            jwalltime => '1:00:00',
        },
        'Rho_Terminase_Predict' => {
            jcomment => 'Search for rho terminases.',
            jcpu => 2,
            jmem => 12,
            jname_prefix => 'rhoterm',
            jprefix => '51',
            jwalltime => '4:00:00',
        },
        'RMats' => {
            jcomment => 'Seek splicing changes from aligned reads with rmats.',
            jcpu => 16,
            jmem => 24,
            jname_prefix => 'rmats',
            jprefix => '90',
            jwalltime => '24:00:00',

        },
        'RNAFold_Windows' => {
            jcomment => 'Invoke rnafold on a rolling sequence window.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'rnafold',
            jprefix => '81',
            jwalltime => '8:00:00',
        },
        'Rosalind_Plus' => {
            jcomment => 'Make the strand with the most ORFs the plus strand.',
            jcpu => 2,
            jmem => 3,
            jname_prefix => 'rplus',
            jprefix => '44',
            jwalltime => '1:00:00',
        },
        'RSEM' => {
            jcomment => 'Gene/tx level quantification from counts.',
            jcpu => 8,
            jmem => 24,
            jname_prefix => 'rsem',
            jprefix => '28',
            jwalltime => '20:00:00',
        },
        'RSEM_Index' => {
            jcomment => 'Index a genome for use with RSEM.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'rsem_index',
            jprefix => '09',
            jwalltime => '10:00:00',
        },
        'Run_Essentiality' => {
            jcomment => 'Invoke the DeJesus essentiality tool',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'essentiality',
            jprefix => '44',
            jwalltime => '24:00:00'
        },
        'Salmon' => {
            jcomment => 'Exceedingly fast quantitation method.',
            jcpu => 8,
            jmem => 24,
            jname_prefix => 'salmon',
            jprefix => '45',
            jwalltime => '4:00:00',
        },
        'Salmon_Index' => {
            jcomment => 'Index a transcriptome for use with salmon.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'salmon_index',
            jprefix => '09',
            jwalltime => '1:00:00',
        },
        'Salmon_Stats' => {
            jcomment => 'Collect mapping stats from salmon.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'salmon_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'Salmon_Strand' => {
            jcomment => 'Estimate library strandedness from the first few reads from a dataset.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'salmon_strand',
            jprefix => '12',
            jwalltime => '00:05:00',
        },
        'Samtools' => {
            jcomment => 'Convert a sam alignment to bam with/out filtering.',
            jcpu => 1,
            jmem => 20,
            jname_prefix => 'samtools',
            jprefix => '21',
            jwalltime => '06:00:00',
        },
        'Shovill' => {
            jcomment => 'Iterative Spades settings by Torsten Seeman.',
            jcpu => 8,
            jmem => 12,
            jname_prefix => 'shovil',
            jprefix => '13',
            jwalltime => '12:00:00',
        },
        'SL_Recorder' => {
            jcomment => 'Extract Spliced Leader reads via pattern matching.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'sl_record',
            jprefix => '50',
            jwalltime => '12:00:00',
        },
        'SLSearch' => {
            jcomment => 'Extract Spliced Leader reads with cutadapt.',
            jcpu => 4,
            jmem => 12,
            jname => 'slsearch',
            jprefix => '50',
            jwalltime => '08:00:00',
        },
        'SL_UTR' => {
            jcomment => 'Get SL addition positions from split reads.',
            jcpu => 1,
            jmem => 24,
            jname => 'slutr',
            jprefix => '50',
            jwalltime => '08:00:00',
        },
        'Snippy' => {
            jcomment => 'Use snippy to search for variants.',
            jcpu => 4,
            jmem => 12,
            jname_prefix => 'snippy',
            jprefix => '44',
            jwalltime => '4:00:00',
        },
        'SNP_Ratio' => {
            jcomment => 'Find variants in CDSs without introns from freebayes/mpileup',
            jcpu => 1,
            jmem => 16,
            jname_prefix => 'snpratio',
            jwalltime => '8:00:00'
        },
        'SNP_Ratio_Intron' => {
            jcomment => 'Find variants in CDSs with introns from freebayes/mpileup',
            jcpu => 1,
            jmem => 24,
            jname_prefix => 'snpratiointron',
            jwalltime => '12:00:00'
        },
        'Sort_Indexes' => {
            jcomment => 'Sort input tnseq data into separate samples.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'sortindexes',
            jprefix => '10',
            jwalltime => '24:00:00',
        },
        'Spladder' => {
            jcomment => 'Quantify differential splicing via spladder.',
            jcpu => 8,
            jmem => 24,
            jname_prefix => 'spladder',
            jprefix => '90',
            jwalltime => '24:00:00',
        },
        'Split_Align_Blast' => {
            jcomment => 'Split a large set of BLAST searches into smaller pieces',
            jcpu => 4,
            jmem => 8,
            jname_prefix => 'splitblast',
            jprefix => '60',
            jwalltime => '12:00:00',
        },
        'Split_Align_Fasta' => {
            jcomment => 'Split a large set of fasta36 searches into smaller pieces',
            jcpu => 8,
            jmem => 12,
            jname_prefix => 'splitfasta',
            jprefix => '60',
            jwalltime => '12:00:00',
        },
        'STAR' => {
            jcomment => 'Transcript abundance estimator.',
            jcpu => 8,
            jmem => 48,
            jname_prefix => 'star',
            jprefix => '48',
            jwalltime => '12:00:00',
        },
        'STAR_Index' => {
            jcomment => 'Index a genome for use with STAR',
            jcpu => 4,
            jmem => 96,
            jname_prefix => 'start_index',
            jprefix => '09',
            jwalltime => '12:00:00',
        },
        'Suppa' => {
            jcomment => 'Quantify differential splicing via suppa.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'suppa',
            jprefix => '90',
            jwalltime => '12:00:00',
        },
        'TA_Check' => {
            jcomment => 'Ensure each read of a TNSeq set starts with TA.',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'tacheck',
            jprefix => '30',
            jwalltime => '12:00:00',
        },
        'Terminase_ORF_Reorder' => {
            jcomment => 'Reorder an assembly to put the putative terminase first.',
            jcpu => 4,
            jmem => 8,
            jname_prefix => 'term_reorder',
            jprefix => '30',
            jwalltime => '1:00:00',
        },
        'Tophat' => {
            jcomment => 'Early splicing aware mapper/aligner.',
            jcpu => 4,
            jmem => 20,
            jname_prefix => 'tophat',
            jprefix => '35',
            jwalltime => '36:00:00',
        },
        'Tophat_Stats' => {
            jcomment => 'Collect mapping stats from tophat.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'tophat_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'Train_Prodigal' => {
            jcomment => 'Train prodigal on a specific sequence dataset.',
            jcpu => 4,
            jmem => 8,
            jname_prefix => 'trainprodigal',
            jprefix => '13',
            jwalltime => '2:00:00',
        },
        'Transdecoder' => {
            jcomment => 'Seek viable CDS in the 6 reading frames of a set of sequences.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'trandecode',
            jprefix => '21',
            jwalltime => '2:00:00',
        },
        'Transit_TPP' => {
            jcomment => 'Invoke the Transit preprocessing pipeline.',
            jcpu => 4,
            jmem => 8,
            jname_prefix => 'tpp',
            jprefix => '09',
            jwalltime => '8:00:00',
        },
        'Transposonpsi' => {
            jcomment => 'Seek transposable elements.',
            jcpu => 2,
            jmem => 4,
            jname_prefix => 'tpsi',
            jprefix => '21',
            jwalltime => '2:00:00',
        },
        'Trimomatic_Pairwise' => {
            jcomment => 'Trim reads with the OG trimmer.',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'trimomaticp',
            jprefix => '01',
            jwalltime => '4:00:00',
        },
        'Trimomatic_Single' => {
            jcomment => 'Trim reads with the OG trimmer.',
            jcpu => 2,
            jmem => 12,
            jname_prefix => 'trimomatics',
            jprefix => '01',
            jwalltime => '4:00:00',
        },
        'Trimomatic_Stats' => {
            jcomment => 'Collect read stats from trimomatic.',
            jcpu => 1,
            jmem => 1,
            jname_prefix => 'trimomatic_stats',
            jprefix => '90',
            jwalltime => '00:01:00',
        },
        'Trinity' => {
            jcomment => 'Automated transcriptome assembly from RNASeq data.',
            jcpu => 8,
            jmem => 80,
            jname_prefix => 'trinity',
            jprefix => '60',
            jwalltime => '36:00:00',
        },
        'Trinity_Post' =>{
            jcomment => 'Post processing of trinity results with RSEM etc.',
            jcpu => 8,
            jmem => 24,
            jname_prefix => 'trinpost',
            jprefix => '61',
            jwalltime => '8:00:00',
        },
        'Trinotate' => {
            jcomment => 'Automated annotation of sequences primarily via blast.',
            jcpu => 8,
            jmem => 36,
            jname_prefix => 'trinotate',
            jprefix => '62',
            jwalltime => '36:00:00',
        },
        'tRNAScan' => {
            jcomment => 'Search (m)tRNA genes.',
            jcpu => 1,
            jmem => 6,
            jname_prefix => 'trnascan',
            jprefix => '28',
            jwalltime => '4:00:00',
        },
        'Umi_Tools_Dedup' => {
            jcomment => 'Remove umis after mapping.',
            jcpu => 4,
            jmem => 48,
            jname_prefix => 'umi_dedup',
            jprefix => '04',
            jwalltime => '24:00:00'
        },
        'Umi_Tools_Extract' => {
            jcomment => 'Use a regex to pull umi reads.',
            jcpu => 1,
            jmem => 20,
            jname_prefix => 'umi_extract',
            jprefix => '05',
            jwalltime => '08:00:00',
        },
        'Unicycler' => {
            jcomment => 'Heuristics to improve spades assemblies.',
            jcpu => 8,
            jmem => 24,
            jname_prefix => 'unicycler',
            jprefix => '13',
            jwalltime => '72:00:00',
        },
        'Unicycler_Filter_Depth' => {
            jcomment => 'Filter an assembly based on relative depth of contigs.',
            jcpu => 1,
            jmem => 4,
            jname_prefix => 'depthfilter',
            jprefix => '13',
            jwalltime => '1:00:00',
        },
        'Velvet' => {
            jcomment => 'yet another kmer assembler',
            jcpu => 4,
            jmem => 24,
            jname_prefix => 'velvet',
            jprefix => '20',
            jwalltime => '8:00:00',
        },
        'Xref_Crispr' => {
            jcomment => 'Seek crispr sequences',
            jcpu => 1,
            jmem => 8,
            jname_prefix => 'xref_crispr',
            jprefix => '80',
            jwalltime => '1:00:00',
        },
    );
    my $datum = {};
    ## This function is fed the returned_values within Get_Vars()
    ## so we can use them to see if the user asked for a specifc prefix/walltime/etc
    ## Changed this to remove the prefix 'j' so that one may override the
    ## spec by specifying the mem/gpu/etc.
    $datum->{jcomment} = $args{comment} if ($args{comment});
    $datum->{jcpu} = $args{cpu} if ($args{cpu});
    $datum->{jprefix} = $args{prefix} if ($args{prefix});
    $datum->{jgpu} = $args{gpu} if ($args{gpu});
    $datum->{jmem} = $args{mem} if ($args{mem});
    $datum->{jname_prefix} = $args{name_prefix} if ($args{name_prefix});
    $datum->{jwalltime} = $args{walltime} if ($args{walltime});
    ## If the user does not provide spec information, use the above hash of hashes
    if (defined($function)) {
        my @function_array = split(/::/, $function);
        my $function = $function_array[$#function_array];
        $function =~ s/_Worker$//g;
        my %spec_info = %{$spec_data{Default}};
        if (defined($spec_data{$function})) {
            %spec_info = %{$spec_data{$function}};
        } else {
            print "The function: $function does not have a spec, using the defaults.\n";
        }
        ## Fill in the cpu/gpu/etc if it was not specified by the job itself.
        for my $k (keys %spec_info) {
            $datum->{$k} = $spec_info{$k} unless ($datum->{$k});
        }
    }

    ## If things are still missing from the spec, fill in the defaults here.
    $datum->{jcomment} = $spec_data{Default}->{jcomment} unless ($datum->{jcomment});
    $datum->{jcpu} = $spec_data{Default}->{jcpu} unless ($datum->{jcpu});
    $datum->{jgpu} = $spec_data{Default}->{jgpu} unless ($datum->{jgpu});
    $datum->{jmem} = $spec_data{Default}->{jmem} unless ($datum->{jmem});
    $datum->{jname_prefix} = $spec_data{Default}->{jprefix} unless ($datum->{jname_prefix});
    $datum->{jprefix} = $spec_data{Default}->{jprefix} unless ($datum->{jprefix});
    $datum->{jwalltime} = $spec_data{Default}->{jwalltime} unless ($datum->{jwalltime});
    return($datum);
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
        "alphafold+" => \$todo_list->{todo}{'Bio::Adventure::Structure::ProteinFold'},
        ## Filter a population genetics dataset via Angsd
        "angsdfilter+" => \$todo_list->{todo}{'Bio::Adventure::PopGen::Angsd_Filter'},
        ## Use my phage pipeline to attempt to annotate a phage assembly.
        "annotatephage+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Annotate_Phage'},
        ## Arbitrary sequence format converter.
        "any2any+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Any2Any'},
        ## Given an assembly, search for tRNA-like genes.
        "aragorn+" => \$todo_list->{todo}{'Bio::Adventure::Feature_Prediction::Aragorn'},
        ## Remap reads to a new assembly to estimate coverage on a per-baseis.
        "assemblycoverage+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Assembly_Coverage'},
        ## Use the logic from Assembly_Coverage to calculate coverage vs. an existing genome.
        "bam2bigwig+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Bam2Bigwig'},
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
        "checkjobs+" => \$todo_list->{todo}{'Bio::Adventure::Check_Job'},
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
        "dantools_rnaseq+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Dantools_RNASeq'},
        "dedupgatk+" => \$todo_list->{todo}{'Bio::Adventure::Convert::GATK_Dedup'},
        "ncbidownload+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Download_NCBI_Accession'},
        "dlassembly+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Download_NCBI_Assembly'},
        "dldataset+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Download_NCBIDatasets_Accession'},
        "dlncbi+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Download_NCBI_Accession'},
        "dlensembl+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Download_Ensembl_Files'},
        "downsampleguess+" => \$todo_list->{todo}{'Bio::Adventure::Map::Downsample_Guess_Strand'},
        "essentialitytas+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Essentiality_TAs'},
        "extendkraken+" => \$todo_list->{todo}{'Bio::Adventure::Index::Extend_Kraken_DB'},
        "extractseqids+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Extract_SeqIDs'},
        "extractsubseq+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Extract_Subseq'},
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
        "fcount+" => \$todo_list->{todo}{'Bio::Adventure::Count::Feature_Counts'},
        "featurecount+" => \$todo_list->{todo}{'Bio::Adventure::Count::Feature_Counts'},
        "featureextract+" => \$todo_list->{todo}{'Bio::Adventure::Annotation_Genbank::Extract_Features'},
        "filterdepth+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Unicycler_Filter_Depth'},
        "filterkraken+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Filter_Host_Kraken'},
        "fragment+" => \$todo_list->{todo}{'Bio::Adventure::SeqMisc::Fragment_Genome'},
        "freebayes+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Freebayes_SNP_Search'},
        "gatkdedup+" => \$todo_list->{todo}{'Bio::Adventure::Convert::GATK_Dedup'},
        "gb2gff+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Gb2Gff'},
        "gff2fasta+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Gff2Fasta'},
        "gfftree+" => \$todo_list->{todo}{'Bio::Adventure::Count::GFF_Tree'},
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
        "juggle+" => \$todo_list->{todo}{'Bio::Adventure::Juggle'},
        "kallisto+" => \$todo_list->{todo}{'Bio::Adventure::Map::Kallisto'},
        "kraken+" => \$todo_list->{todo}{'Bio::Adventure::Count::Kraken'},
        "mash+" => \$todo_list->{todo}{'Bio::Adventure::Count::Mash'},
        "mauve+" => \$todo_list->{todo}{'Bio::Adventure::Align::MauveAligner'},
        "mergeannotations+" => \$todo_list->{todo}{'Bio::Adventure::Metadata::Merge_Annotations'},
        "mergecds+" => \$todo_list->{todo}{'Bio::Adventure::Annotation_Genbank::Merge_CDS_Predictions'},
        "mergeparse+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Merge_Parse_Blast'},
        "mergeprodigal+" => \$todo_list->{todo}{'Bio::Adventure::Metadata::Merge_Annot_Prodigal'},
        "mimap+" => \$todo_list->{todo}{'Bio::Adventure::MiRNA::Mi_Map'},
        "mpileup+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Mpileup_SNP_Search'},
        "nfrna+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Nextflow_RNASeq'},
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
        "polyarecord+" => \$todo_list->{todo}{'Bio::Adventure::Map::PolyA_Recorder'},
        "posttrinity+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity_Post'},
        "pmauve+" => \$todo_list->{todo}{'Bio::Adventure::Align::Progressive_Mauve'},
        "prokka+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Prokka'},
        "proteinfold+" => \$todo_list->{todo}{'Bio::Adventure::Structure::ProteinFold'},
        "proteinpair+" => \$todo_list->{todo}{'Bio::Adventure::Structure::ProteinFold_PairIDs'},
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
        "sl_record+" => \$todo_list->{todo}{'Bio::Adventure::Splicing::SL_Recorder'},
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
        "splitfasta+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Split_Fasta'},
        "spring+" => \$todo_list->{todo}{'Bio::Adventure::Compress::Spring'},
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
        print "TESTME: Found args{method}: $args{method}\n";
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
            } else {
                print "Going on an adventure to ${job} in the $args{task} context.\n";
            }
            my $final = \&${job};
            push(@methods_to_run, $final);
        }
    }
    return(\@methods_to_run);
}

sub Get_Function_Vars {
    my ($class, $funcall_args) = @_;
    my $current_state = $class->current_state;
    my $last_input;
    if ($funcall_args->{input}) {
        $last_input = $class->Bio::Adventure::Config::Get_Input_Info($funcall_args->{input});
        use Data::Dumper;
        print Dumper $last_input;
        my $number_inputs = $class->Append_Inputs($last_input);
    } else {
        $last_input = $class->Get_Last_Input();
    }
    print "TESTME: Get_Function_Vars after if()\n";
    print Dumper $last_input;
    my %all_function_args = %{$class->function_defaults};
    my $function = Bio::Adventure::Config::Get_Caller();
    my $function_args = $all_function_args{$function};
    my $job_spec = $function_args->{spec};
    print "TESTME: In Config::Get_Vars, about to invoke Get_Class_Paths\n";
    my $class_paths = $class->Bio::Adventure::Config::Get_Class_Paths(spec => $job_spec);
    print "TESTME: After Get_Class_Paths: output_dir: $class_paths->{output_dir}\n";
    ## Define the set of accessible interpolation variables here.
    my %available = %{$class->current_state};
    $available{output_dir} = q{};
    ## Add the results from Get_Class_Spec() to tthe set of availble
    ## interpolators, in addition add them to the function_args{paths} slot.
    print "After get_current_state: $class_paths->{output_dir}\n";
    for my $k (keys %{$class_paths}) {
        $available{$k} = $class_paths->{$k};
        $function_args->{paths}->{$k} = $class_paths->{$k};
    }
    print "TESTME: Adding Functioncall arguments to available:\n";
    for my $k (keys %{$funcall_args}) {
        print "Adding $k: $funcall_args->{$k} to available.\n";
        $available{$k} = $funcall_args->{$k};
    }

    my $fb = $last_input->[0]->{filebase_extension};
    $available{r1_base} = qq"$available{output_dir}/$last_input->[0]->{filebase_extension}";
    if ($last_input->[1]) {
        print "SETTING R2_BASE\n";
        $available{r2_base} = qq"$available{output_dir}/$last_input->[1]->{filebase_extension}";
    }
    $function_args->{caller} = $function;

    ## Some specifications may provide 2 outputs but only use 1
    ## E.g. fastp can have a r1+r2 output or just r1
    ## We can check for this case by checking r2_base.
    while (my ($key) = each(%{$function_args->{paths}})) {
        if ($key =~ /_r1$/) {
            if ($available{r1_base}) {
                print "keeping $available{r1_base}\n";
            } else {
                delete $function_args->{paths}->{$key};
            }
        } elsif ($key =~ /_r2$/) {
            if ($available{r2_base}) {
                print "Keeping $available{r2_base}\n";
            } else {
                print "DELETING $key\n";
                delete $function_args->{paths}->{$key};
            }
        }
    }

    ## Now we can create the outputs interpolation variable
    my $outputs = q{};
    for my $k (sort keys %{$function_args->{paths}}) {
        if ($k =~ /^output/) {
            $outputs .= qq"$function_args->{paths}->{$k}:";
        }
    }
    $outputs =~ s/:$//g;
    print "AFTER OUTPUTS_LOOP: ${outputs}\n";
    $available{outputs} = $outputs;

    for (values %{$function_args->{paths}}) {
        if ($_ =~ /[\$\@]/) {
            $_ = interpolate(interpolate($_, \%available), \%available);
        }
    }

    $function_args->{inputs} = $last_input;
    ## Add this set of options to the option list in the cyoa object.
    ## Thus, we can trivially pull all of those options in future jobs.
    my @option_array = @{$class->options};
    push(@option_array, $function_args);
    $class->options(\@option_array);
    return($function_args);
}

## For now, just overwrite keys from global->function->ini->json->cmdline
## in order to create the hash of the current state.
## Once I start getting into the dirty implementation details
## I think I will want to use those portions of the function hash to make this
## into something which sets the parameters of each step of a process.
## But while I am just experimenting with how I want to work with this, a simpler
## approach makes sense.
sub Merge_Args {
    my %args = @_;
    my $global = $args{global};
    my $cmdline = $args{cmdline};
    my $ini = $args{ini};
    my $json = $args{json};
    my %current;
    for my $k (keys %{$global}) {
        $current{$k} = $global->{$k};
    }
    for my $k (keys %{$ini}) {
        $current{$k} = $ini->{$k};
    }
    for my $k (keys %{$json}) {
        $current{$k} = $json->{$k};
    }
    for my $k (keys %{$cmdline}) {
        $current{$k} = $cmdline->{$k};
    }
    return(\%current);
}

sub Merge_Current_Ini {
    my %args = @_;
    my %current = $args{current};
    my $steps = $args{steps};
    my $defaults = $args{defaults};
    for my $def_key (keys %{$defaults}) {
        my $current_key = $def_key;
        while (defined($current{$current_key})) {
            $current_key .= "_";
            say "${def_key} is already defined in the current state, adding an _.";
        }
        $current{$current_key} = $defaults->{$def_key};
    }
    for my $step_key (keys %{$steps}) {
        my $current_key = $step_key;
        while (defined($current{$current_key})) {
            $current_key .= "_";
            say "${step_key} is already defined in the current state, adding an _.";
        }
        $current{$current_key} = $steps->{$step_key};
    }
    return(\%current);
}

sub Update_Spec {
    my ($class, %args) = @_;
    my $function_info = $args{function_info};
    my $caller = $function_info->{caller};
    my %current_spec = %{$function_info->{spec}};
    if ($function_info->{template}) {
        my $jstring = q{};
        my $tt = Template->new({
            INCLUDE_PATH => $template_dir,
            TRIM => 0,
            INTERPOLATE => 0,});
        print "About to process $function_info->{template} from $template_dir\n";
        use Data::Dumper;
        print Dumper $args{template_args};
        my $processed = $tt->process($function_info->{template}, $args{template_args}, \$jstring) || die $tt->error();
        $current_spec{jstring} = $jstring;
        use Data::Dumper;
        print Dumper %current_spec;
    } elsif ($args{jstring}) {
        print "Adding jstring.\n";
        $current_spec{jstring} = $args{jstring};
    } else {
        print "Neither a template nor a jstring.\n";
    }
    $current_spec{modules} = $function_info->{modules};
    if (!defined($current_spec{language})) {
        $current_spec{language} = 'bash';
    }
    if (!defined($current_spec{jname})) {
        $current_spec{jname} = $class->options->{jname};
    }
    $class->spec(\%current_spec);
}

1;
