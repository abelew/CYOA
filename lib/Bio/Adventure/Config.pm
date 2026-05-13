package Bio::Adventure::Config;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use v5.14;
use Cwd qw"abs_path fastcwd";
use Data::Printer;
use Exporter;
use File::Basename qw"basename dirname";
use FileHandle;
use File::Path qw"make_path";
use File::ShareDir qw"dist_file module_dir dist_dir";
use File::Which qw"which";
use Getopt::Long qw"GetOptionsFromArray";
use JSON;
use POSIX qw"ceil";
use String::Interpolate qw"interpolate";
use Template;
my $template_base = dist_dir('Bio-Adventure');
my $template_dir = qq"${template_base}/jobs";

our @EXPORT = qw"Get_Modules Get_Menus Get_TODOs";
our @ISA = qw"Exporter";

## 202605 survey: Invoked when reading an ini file
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

## 202605 Survey: Invoked at the beginning of BUILD to get command line arguments.
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
## 202605 Survey: called all over the place.
sub Get_Caller {
    my ($class, %args) = @_;
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
        } elsif ($sub =~ /^Bio::.*Update_State$/) {
            next COUNTER;
        }
        ## All callers which are needing modules have the form Bio::Adventure::Category::Function
        $function = $sub;
        $function =~ s/^(.*::)(\w+)$/$2/g;
        last COUNTER if (defined($function));
    }
    if (!defined($function)) {
        ## warn "The function is undefined in Config::Get_Caller().";
        $function = 'Default';
    }

    ## Now that we have the current caller, update a few pieces of the parent class.
    ## I made a decision to simultaneously store the set of all function defaults as well
    ## as the function defaults for the current caller.  I did this for two reasons:
    ## 1.  It should make it simpler in the various child functions to pull useful stuff.
    ## 2.  It should make it possible to check one function from another, I am thinking
    ##     this will be useful in the case of multi-function paths
    ## Thus, I will update the function_args slot of the class with the portion of
    ## function_defaults{} for this caller and also update spec with the portion of that
    ## which tells the default walltime/mem/etc.  Finally, I will set caller with whatever
    ## this function happens to be.
    if (!defined($function)) {
        return(undef);
    }
    ## Update the function args with this caller.
    my %all_args = %{$class->function_defaults};
    if (!defined($all_args{$function})) {
        return(undef);
    }
    my $function_defs = $all_args{$function};
    my %fn_args = %{$all_args{$function}};
    $class->function_args(\%fn_args);
    $class->caller($function);
    return($function);
}

## My thought here: Make it possible to have a set of 'global' defaults
## which may be overridden at a few levels, branching down to the indivdual sample.
## I organize my data as /type-of-project/individual-project/preprocessing/sample/
## so there are likely 4 levels at which one may most likely want to make changes.

## 202605 Survey: Invoked by BUILD to read ini-style config directives.
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

## 202605 survey: Invoked in BUILD to collect default values.
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
## 202605 Survey: Called in build to set up the menu, I want to move this to a config file.
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

## 202605 Survey: Called in BUILD to pull (currently only module) shell functions;
## The intent is to add conda/spack and whatever other future tools crop up.
## With that in mind, it kind of only makes sense in the context of module because
## (with the exception of the various condas) I expect no other tools to have multiple
## implementations -- but the existence of miniconda/anaconda/etc suggest it is worth while.
sub Get_Shell_Functions {
    my (%args) = @_;
    ## Check that the module command is available as a bash function.
    ## I was initially going to do this in BUILD(), but I think that might mess up jobs
    ## which set the jprefix parameter.
    ## This command, if run via bash, requires a controlling pty in order to find the
    ## bash function 'module'.  slurm and the various managers will not provide that
    ## so, we can use typescript to allocate one post-facto.
    my $module_test = qq!/usr/bin/script -q -c 'bash -i -c "type module"' /dev/null |!;
    my $modulecmd_handle = FileHandle->new($module_test);
    my $shell_text = qq!function get_sigterm {
  echo "A SIGTERM was sent to \$JNAME." >> $args{bash_log}
  exit 1
}
export -f get_sigterm
trap get_sigterm SIGTERM

function get_sigerr {
  echo "A ERR was sent to \$JNAME." >> $args{bash_log}
  exit 1
}
export -f get_sigerr
trap get_sigerr ERR

function module() {
!;

    my $num_modulecmd_lines = 0;
  MOD: while (my $line = <$modulecmd_handle>) {
        chomp $line;
        next MOD if ($line =~ /is a function/);
        next MOD if ($line =~ /^module/);
        if ($line =~ /^\{\s*$/) {
            next MOD;
        }
        ## In case there are weird encodings.
        $line =~ s/\s*$//g;
        $num_modulecmd_lines++;
        $shell_text .= qq"${line}\n";
    }
    $shell_text .= q"export -f module

";
    $modulecmd_handle->close();
    if (-f 'typescript') {
        print "The stupid typescript file still got created.\n";
        rm('typescript');
    }
    ## If module is not in the current environment, check for the
    ## modulecmd executable and just use that.
    if ($num_modulecmd_lines == 0) {
        my $modulecmd = which('modulecmd');
        if ($modulecmd) {
            $shell_text .= qq!${modulecmd} bash "$@"
export -f module

!;
        } else {
            $shell_text .= q!echo "The module command does not appear to exist".
}
export -f module

!;
        }
    }

    my $spack_test = qq!/usr/bin/script -q -c 'bash -i -c "type spack"' /dev/null |!;
    my $spack_handle = FileHandle->new($spack_test);
    $shell_text .= qq!function spack() {
!;
    my $num_spack_lines = 0;
  SP: while (my $line = <$spack_handle>) {
        chomp $line;
        next SP if ($line =~ /is a function/);
        next SP if ($line =~ /^spack/);
        next SP if ($line =~ /^\{\s*$/);
        ## In case there are weird encodings.
        $line =~ s/\s*$//g;
        $num_spack_lines++;
        $shell_text .= qq"${line}\n";
    }
    $shell_text .= q"export -f spack

";
    $spack_handle->close();
    if ($num_spack_lines == 0) {
        $shell_text .= q!echo "The spack command does not appear to exist".
}
export -f spack

!;
    }

    my $conda_test = which('conda');
    if (!$conda_test) {
        $shell_text .= q!
function conda() {
  echo "The conda command does not appear to exist."
}
export -f conda
!;
    }
    return($shell_text);
}

=head2 C<Get_Input_Info>

Given a file/directory name, provide a set of paths which will
hopefully prove useful for the various functions in CYOA.  I think
this function is a good candidate for replacing Check_Input() below.

=cut
## 202605 survey: Called by BUILD to create a single data structure
## describing the inputs to cyoa.
sub Get_Input_Info {
    my ($class, @inputs) = @_;
    my %ret = ();
    my $input_string = q{};
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
        $input_string .= qq"${in}:";
        my $filename = basename($in);
        my $filebase_compress = basename($in, ('.gz', '.xz', '.bz2'));
        ## FIXME: Replace this with options->{suffix_array}
        my @exts = ('.fastq', '.fasta', '.fsa', '.faa', 'fna', '.fa', '.ffn',
                    'gff', '.gff3', '.gbk', '.gbf', '.sqn', 'tbl', '.tsv', '.txt');
        my $filebase_extension = basename($filebase_compress, @exts);
        my $directory = dirname($in);
        my $dirname = basename($directory);
        ## Check that doing basename(dirname(input)) leaves us with something interesting
        if ($dirname eq '.') {
            $dirname = q{};
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

        my $suffix = $filename;
        $suffix =~ s/\.[xgb]z$//g;
        if ($suffix =~ /\..*$/) {
            $suffix =~ s/^.*\.(.*)$/$1/g;
        }

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
            input => $local_path,
            suffix => $suffix,);
        push(@outputs, \%ret);
    }
    ## Ensure the current state has the input string.
    $input_string =~ s/:$//g;
    $class->{current_state}->{input} = $input_string;
    return(\@outputs);
}

=head2 C<Get_Class_Paths>

This should be responsible for acquiring paths used by anything in
this instance of the Adventure class.

TODO: Once this is settled, take out the function-specific material
and move it to Get_Function_Paths() which will contain stuff
like output_dir.

=cut
## 202605 Survey: Called by BUILD to set up and/or find all the various database directories
sub Get_Global_Paths {
    my ($class, %args) = @_;
    my $state_ref = $class->current_state;
    my %current_state = %{$state_ref};
    my $libtype = 'genome';
    $libtype = $current_state{libtype} if (defined($current_state{libtype}));
    my $libpath_prefix = qq"$current_state{libpath}/${libtype}";
    my $libdir_prefix = qq"$current_state{libdir}/${libtype}";
    ##my $cwd_name = basename(cwd());
    my %class_paths = ();
    ## Ok, something in cwd/basename is causing a problem.
    $class_paths{cwd_name} = basename(fastcwd());
    $class_paths{blast_dir} = qq"${libpath_prefix}/blast";
    $class_paths{codon_dir} = qq"${libpath_prefix}/codon_tables";
    $class_paths{fasta_dir} = qq"${libpath_prefix}/fasta";
    $class_paths{gbk_dir} = qq"${libpath_prefix}/genbank";
    $class_paths{gff_dir} = qq"${libpath_prefix}/gff";
    $class_paths{gtf_dir} = qq"${libpath_prefix}/gtf";
    $class_paths{hmm_dir} = qq"${libpath_prefix}/hmm";
    $class_paths{index_prefix} = qq"${libpath_prefix}/indexes";
    $class_paths{index_prefix_shell} = qq"${libdir_prefix}/indexes";
    make_path($class_paths{blast_dir}) unless (-d $class_paths{blast_dir});
    make_path($class_paths{codon_dir}) unless (-d $class_paths{codon_dir});
    make_path($class_paths{fasta_dir}) unless (-d $class_paths{fasta_dir});
    make_path($class_paths{gff_dir}) unless (-d $class_paths{gff_dir});
    make_path($class_paths{gtf_dir}) unless (-d $class_paths{gtf_dir});
    make_path($class_paths{gbk_dir}) unless (-d $class_paths{gbk_dir});
    make_path($class_paths{hmm_dir}) unless (-d $class_paths{hmm_dir});
    make_path($class_paths{index_prefix}) unless (-d $class_paths{index_prefix});
    if (defined($ENV{BLAST_DB})) {
        if ($ENV{BLAST_DB} =~ /:/) {
            my @blast = split(/:/, $ENV{BLAST_DB});
            unless (grep(/$class_paths{blast_dir}/, @blast)) {
                $ENV{BLAST_DB} = qq"$ENV{BLAST_DB}:$class_paths{blast_dir}";
            }
        } else {
            unless ($class_paths{blast_dir} eq $ENV{BLAST_DB}) {
                $ENV{BLAST_DB} = qq"$ENV{BLAST_DB}:$class_paths{blast_dir}";
            }
        }
    } else {
        $ENV{BLAST_DB} = $class_paths{blast_dir};
    }
    if (defined($current_state{species})) {
        $class_paths{fasta} = qq"$class_paths{fasta_dir}/$current_state{species}.fasta";
        $class_paths{fasta_shell} = qq"${libdir_prefix}/fasta/$current_state{species}.fasta";
        $class_paths{gff} = qq"$class_paths{gff_dir}/$current_state{species}.gff";
        $class_paths{gff_shell} = qq"${libdir_prefix}/gff/$current_state{species}.gff";
        $class_paths{gtf} = qq"$class_paths{gtf_dir}/$current_state{species}.gtf";
        $class_paths{gtf_shell} = qq"${libdir_prefix}/gtf/$current_state{species}.gtf";
        $class_paths{gbk} = qq"$class_paths{gbk_dir}/$current_state{species}.gbk";
        $class_paths{gbk_shell} = qq"${libdir_prefix}/genbank/$current_state{species}.gbk";
    }
    my $set = $class->_set_global_paths(\%class_paths);
    return(\%class_paths);
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
## 202605 Survey: Called in BUILD, I want this also to move to a config file.
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

## 202605 Survey: Responsible for providing all the local variables to arbitrary functions.
sub Get_Function_Vars {
    my ($class, $funcall_args) = @_;
    my $caller = $class->Bio::Adventure::Config::Get_Caller();
    my $state_ref = $class->current_state;
    ## To an annoying degree this is effectively redoing the
    ## BUILD function from Adventure.pm
    ## I should think about how I can clean this up so that
    ## I am not doing and undoing these declarations.
    ## The caveat though is that I want to have a clean slate for
    ## every downstream function call invoked by Adventure.pm
    ## Until I think that through, keep in mind that my
    ## order of precedence is:
    ## global_defaults < instance_args < function_defaults < function_args
    my %state = %{$state_ref};
    my $global_defaults = $class->global_defaults;
    for my $k (keys %{$global_defaults}) {
        $state{$k} = $global_defaults->{$k};
    }
    my %instance_args = %{$class->instance_args};
    for my $j (keys %instance_args) {
        $state{$j} = $instance_args{$j};
    }
    my $function_state = $class->function_defaults;
    my %gl_paths = %{$class->global_paths};
    my $last_input;
    if ($funcall_args->{input}) {
        $class->{current_state}->{input} = $funcall_args->{input};
        $last_input = $class->Bio::Adventure::Config::Get_Input_Info(
            $funcall_args->{input});
        my $number_inputs = $class->Append_Inputs($last_input);
    } else {
        $last_input = $class->Get_Last_Input();
        $class->{current_state}->{input} = $last_input;
    }
    $class->current_inputs($last_input);
    ## Ensure that current_state->input matches this.

    my %all_function_args = %{$class->function_defaults};
    ## Note, I think one must make a copy of function_args here so that future invocations
    ## will not already have the variables interpolated.
    my %function_args = %{$all_function_args{$caller}};
    for my $k (keys %function_args) {
        $state{$k} = $function_args{$k};
    }
    for my $j (keys %{$funcall_args}) {
        $state{$j} = $funcall_args->{$j};
    }
    ## Define the set of accessible interpolation variables here.
    my @info_array = @{$last_input};
    my %first_file_info = %{$info_array[0]};
    my $first_input_filename = $first_file_info{filename};
    my @input_array = ();
    my @output_array = ();
    my @input_base_array = ();
    my @output_base_array = ();
    for my $input (@info_array) {
        push(@input_array, $first_file_info{filename});
        push(@input_base_array, $first_file_info{basename});
    }

    $state{output_prefix} = qq"$state{output_basedir}/$state{jprefix}";
    ## The suffix of the input file(s) excluding any compression information.
    $state{input_suffix} = $first_file_info{suffix};
    ## The directory name of the input file(s)
    $state{input_dirname} = $first_file_info{dirname};
    ## The basename of the input file(s) excluding erroneous stuff like compress/suffix/r1/etc
    $state{output_base} = $first_file_info{basename};
    ## Array of input filenames
    $state{input_array} = \@input_array;
    ## Array of the input basenames.
    $state{input_base_array} = \@input_base_array;
    ## The first input file excluding compression/suffixes
    $state{nosuffix} = qq"$first_file_info{dirname}/$first_file_info{basename}";
    ## The full output directory path, this may need some flexibility and/or to be interpolated.
    ## E.g. I want to potentially have the flexibility to add parameters/genome etc to this.
    $state{output_dir} = qq"$state{output_prefix}_$state{jname}";

    ## I want to run most/all jobs through time, so record that filename now.
    $state{time_file} = qq"$state{output_dir}/$state{jname}.time";

    ## The primary inputs/outputs are defined,
    ## now add anything from the functions.json for this function

    ## Add the input information to %available
    ## Note that the input information is a list of hash references
    ## and we are going to usually/always use the first input for stuff like basename/filename/etc.
    for my $k (keys %first_file_info) {
        $state{$k} = $first_file_info{$k};
    }
    ## output_dir is provided by fn_paths
    ## Add any arguments added via a function call.
    for my $k (keys %{$funcall_args}) {
        $state{$k} = $funcall_args->{$k};
    }
    for my $e (keys %ENV) {
        if ($e eq 'HOME' || $e eq 'PATH') {
            $state{$e} = $ENV{$e};
        }
    }
    ## Add environment variables to available

    my $first_file_base = $last_input->[0]->{filebase_extension};
    ## This depends on fn_paths
    $state{r1_base} = qq"$state{output_dir}/${first_file_base}";
    if ($last_input->[1]) {
        my $second_file_base = $last_input->[1]->{filebase_extension};
        $state{r2_base} = qq"$state{output_dir}/${second_file_base}";
    }

    ## Some specifications may provide 2 outputs but only use 1
    ## E.g. fastp can have a r1+r2 output or just r1
    ## We can check for this case by checking r2_base.
    while (my ($key) = each(%state)) {
        if ($key =~ /_r1$/) {
            if ($state{r1_base}) {
                ## print "keeping $current_state{r1_base}\n";
            } else {
                delete $state{$key};
            }
        } elsif ($key =~ /_r2$/) {
            if ($state{r2_base}) {
                ## print "Keeping $current_state{r2_base}\n";
            } else {
                ## print "DELETING $key\n";
                delete $state{$key};
            }
        }
    }

    ## Now we can create the outputs interpolation variable
    my $outputs = q{};
    STATE: for my $k (sort keys %state) {
          next STATE if ($k =~ /_comment$/);
          if ($k =~ /^output/) {
              $outputs .= qq"$state{$k}:";
          }
      }

    my $sw_added = 0;
    if ($state{modules}) {
        $sw_added++;
        $state{module_string} = qq"module add $state{modules}";
    }
    if ($state{conda}) {
        $sw_added++;
        $state{conda_string} = qq"conda load $state{conda}";
    }
    if ($state{spack}) {
        $sw_added++;
        $state{spack_string} = qq"spack load $state{spack}";
    }
    if ($sw_added == 0) {
        warn("This expects at least one software package to be added to the environment.");
    }

    $outputs =~ s/:$//g;
    $state{outputs} = $outputs;
    $function_args{inputs} = $last_input;
    my $interp_state_ref = $class->self_interp(hashref => \%state);
    $interp_state_ref = $class->self_interp(hashref => $interp_state_ref);
    ## Add this set of options to the option list in the cyoa object.
    ## Thus, we can trivially pull all of those options in future jobs.
    my @option_array = @{$class->options};
    push(@option_array, $interp_state_ref);
    $class->options(\@option_array);
    $class->current_state($interp_state_ref);
    return($interp_state_ref);
}

## For now, just overwrite keys from global->function->ini->json->cmdline
## in order to create the hash of the current state.
## Once I start getting into the dirty implementation details
## I think I will want to use those portions of the function hash to make this
## into something which sets the parameters of each step of a process.
## But while I am just experimenting with how I want to work with this, a simpler
## approach makes sense.

## 202605 survey: Called in BUILD to bring together ini/json/etc
sub Merge_Args {
    my %args = @_;
    ## This is the reverse order of precedence when overwriting parameters
    ## e.g. globals are overwritten by the ini by json by cmdline by instance
    my $global = $args{global};
    my $ini = $args{ini};
    my $json = $args{json};
    my $cmdline = $args{cmdline};
    my $instance = $args{instance};
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
    for my $k (keys %{$instance}) {
        $current{$k} = $instance->{$k};
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

## 202605 survey: Called in arbitrary functions to set the state before submitting a job.
sub Update_State {
    my ($class, %args) = @_;
    my %options;
    if (defined($args{options})) {
        %options = %{$args{options}};
    } else {
        my $function_info = Bio::Adventure::Config::Get_Function_Vars($class, \%args);
        %options = %{$function_info};
    }
    my %current_state = %{$class->current_state};

    ## FIXME:
    ## There are a couple of pieces of state set in Get_Function_Vars()
    ## in a way which means they will not get new interpolated variables
    ## freshly set variables, most notably output_prefix
    ## I think I should have output_prefix set in defaults.json
    ## and then changes to jprefix etc should carry over.
    ## Until I do that, set output_prefix here:
    $options{output_prefix} = qq"$current_state{output_basedir}/$options{jprefix}";
    $options{output_dir} = qq"$options{output_prefix}_$options{jname}";
    $options{time_file} = qq"$options{output_dir}/$options{jname}.time";
    if (-d $options{output_dir}) {
        ## print $state{output_dir} already exists.
    } else {
        make_path($options{output_dir});
    }
    ## Unfortunately, this cannot fix the problem, because the arbitrary output
    ## variables which are created all over the place have already interpolated
    ## the original jprefix containing string.
    ## I think therefore I need to have a single function responsible for
    ## performing the interpolations and making sure that it is run
    ## only in the update_state() or submit() function.
    ## In addition, any make_path() or similar secondary effect statements
    ## must but run at the same time.

    ## END FIXME:

    my $options_ref = $class->self_interp(hashref => \%options);
    if ($options{template}) {
        my $jstring = q{};
        my $templater = Template->new({
            INCLUDE_PATH => $template_dir,
            TRIM => 0,
            PRE_CHOMP => 0,
            POST_CHOMP => 0,
            INTERPOLATE => 0,});
        print "About to process $options{template} from $template_dir\n";
        ##my $processed = $tt->process($function_info->{template}, $args{template_args}, \$jstring) || die $tt->error();
        my $processed = $templater->process($options_ref->{template}, $options_ref, \$jstring);
        $current_state{jstring} = $jstring;
    } elsif ($args{jstring}) {
        print "Adding jstring.\n";
        $current_state{jstring} = $args{jstring};
    } else {
        print "Neither a template nor a jstring.\n";
    }

    if (!defined($current_state{language})) {
        $current_state{language} = 'bash';
    }
    ## I would much prefer to not need the following because
    ## if I do this, I cannot reuse the cyoa instance because this modifies
    ## the pre-interpolated variables.
    my @all_options = @{$class->options};
    push(@all_options, $options_ref);
    my $options_set = $class->options(\@all_options);
    my $set = $class->current_state(\%current_state);
    return($set);
}

1;
