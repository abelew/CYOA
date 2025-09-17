package Bio::Adventure::Defaults;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Which qw"which";

sub scalar_which {
    my $exe = $_[0];
    my $path = which($exe);
    return($path);
}

## Lets move all the default values here.
has adapter_file => (is => 'rw', default => undef);
has align_blast_format => (is => 'rw', default => 5); ## Which alignment type should we use? (5 is blastxml)
has align_jobs => (is => 'rw', default => 40); ## How many blast/fasta alignment jobs should we make when splitting alignments across nodes?
has align_parse => (is => 'rw', default => 1); ## Parse blast searches into a table?
has arbitrary => (is => 'rw', default => ''); ## Extra arbitrary arguments to pass
has array_start => (is => 'rw', default => 100);
has bamfile => (is => 'rw', default => undef); ## Default bam file for converting/reading/etc.
has basedir => (is => 'rw', default => cwd()); ## This was cwd() but I think that may cause problems.
has bash_path => (is => 'rw', default => scalar_which('bash'));
has best_only => (is => 'rw', default => 0); ## keep only the best search result when performing alignments?
has blast_args => (is => 'rw', default => ' -evalue 10 '); ## Default blast parameters
has blast_tool => (is => 'rw', default => 'blastn'); ## Default blast tool to use
has bt_default => (is => 'rw', default => '--best'); ## Default bt1 arguments.
has bt_varg => (is => 'rw', default => '-v 0');
has bt_marg => (is => 'rw', default => '-M 0');
has bt_larg => (is => 'rw', default => '-y -l 15');
has bt2_args => (is => 'rw', default => ' --very-sensitive -L 14 '); ## My favorite bowtie2 arguments
has btmulti => (is => 'rw', default => 0); ## Perform multiple bowtie searches?
has bwa_method => (is => 'rw', default => 'mem,aln'); ## Default bwa method to use.
has chosen_tag => (is => 'rw', default => 'ODDS');
has clean => (is => 'rw', default => 0); ## Cleanup after yourself?
has cluster => (is => 'rw', default => undef); ## Are we running on a cluster?
has comment => (is => 'rw', default => undef); ## Set a comment in running slurm/bash/etc scripts.
has complexity => (is => 'rw', default => 1);
has compress => (is => 'rw', default => 1); ## Compress output files?
has conda_string => (is => 'rw', default => undef);
has condition_column => (is => 'rw', default => undef);
has config => (is => 'rw', default => undef); ## Not sure
has count => (is => 'rw', default => 1); ## Quantify reads after mapping?
has correction => (is => 'rw', default => 1); ## Perform correction when using fastp?
has coverage => (is => 'rw', default => undef); ## Provide a coverage cutoff
has coverage_tag => (is => 'rw', default => 'DP');
has csv_file => (is => 'rw', default => 'all_samples.csv'); ## Default csv file to read/write.
has cutoff => (is => 'rw', default => 0.05); ## Default cutoff (looking at your vcftools, e.g. I haven't changed those yet).
has decoy => (is => 'rw', default => 1);     ## Add decoys
has debug => (is => 'rw', default => 0); ## Print debugging information.
has deduplicate => (is => 'rw', default => 1); ## Perform deduplication when using fastp
has delimiter => (is => 'rw', default => '[;:,]');
has denominator => (is => 'rw', default => undef);
has directories => (is => 'rw', default => undef); ## Apply a command to multiple input directories.
has download => (is => 'rw', default => 1);
has duplicate => (is => 'rw', default => 0);
has email => (is => 'rw', default => 'abelew@umd.edu');
has end => (is => 'rw', default => 10000);
has evalue => (is => 'rw', default => 0.001); ## Default e-value cutoff
has fasta_args => (is => 'rw', default => '-b 20 -d 20'); ## Default arguments for the fasta36 suite
has fasta_tool => (is => 'rw', default => 'ggsearch36'); ## Which fasta36 program to run?
has fastp => (is => 'rw', default => 'check'); ## Perform fastqc in a pipeline?
has fastqc => (is => 'rw', default => 'check'); ## Perform fastqc in a pipeline?
has file_column => (is => 'rw', default => undef);
has filter => (is => 'rw', default => 1); ## When performing an assembly, do a host filter?
has filtered => (is => 'rw', default => 'unfiltered'); ## Whether or not Fastqc is running on filtered data.
has freebayes => (is => 'rw', default => 0);
has fsa_input => (is => 'rw'); ## fsa genome output file for creating a genbank file
has gcode => (is => 'rw', default => '11'); ## Choose a genetic code
has genome => (is => 'rw', default => undef); ## Choose a genome to work on.
has genus => (is => 'rw', default => undef); ## Choose a genus when using prokka and potentially others like kraken
has gff => (is => 'rw', default => undef); ## Feature file to read/write
has gff_tag => (is => 'rw', default => 'ID'); ## Likely redundant with htseq_id
## Ahh I remember, htseq_type was added to facilitate performing multiple htseq counts on multiple gff files.
## Notably the interCDS for bacterial genomes.
has gff_type => (is => 'rw', default => 'gene'); ## When blank, do it on the whole gff file, otherwise use that suffix.
has gff_cds_parent_type => (is => 'rw', default => 'mRNA');
has gff_cds_type => (is => 'rw', default => 'CDS');
has help => (is => 'rw', default => undef); ## Ask for help?
has hisat_args => (is => 'rw', default => ' --sensitive ');
has htseq_args => (is => 'rw', default => ' --order=name --idattr=gene_id --minaqual=10 --type=exon --stranded=yes --mode=union '); ## Most likely htseq options
has identity => (is => 'rw', default => 70); ## Alignment specific identity cutoff
has index_file => (is => 'rw', default => 'indexes.txt'); ## File containing indexes:sampleIDs when demultiplexing samples - likely tnseq
has index_hash => (is => 'rw', default => undef);
has informat => (is => 'rw', default => '.fastq');
has initial_input => (is => 'rw', default => undef);
has input => (is => 'rw', default => undef); ## Generic input argument
has inputs => (is => 'rw', default => undef); ## Generic input argument
has input_abricate => (is => 'rw', default => 'outputs/12abricate_10prokka_09termreorder_08phageterm_07rosalind_plus/abricate_combined.tsv'); ## Used when merging annotation files into a xlsx/tbl/gbk file.
has input_aragorn => (is => 'rw', default => 'outputs/21aragorn/aragorn.txt'); ## Used when merging annotation files into a xlsx/tbl/gbk file.
has input_classifier => (is => 'rw', default => 'outputs/18classifier/ictv_filtered.tsv'); ## Similar taxa detected by tblastx
has input_genbank => (is => 'rw', default => undef); ## Existing genbank file for merging annotations.
has input_glimmer => (is => 'rw', default => 'outputs/16glimmer/glimmer.predict');
has input_fastq => (is => 'rw', default => undef);
has input_interpro => (is => 'rw', default => 'outputs/13_interproscan_10prokka_09termreorder_08phageterm_07rosalind_plus/interproscan.tsv'); ## interpro output file when merging annotations.
has input_phageterm => (is => 'rw', default => 'outputs/08phageterm_07rosalind_plus/direct-term-repeats.fasta'); ## phageterm output file when merging annotations.
has input_phanotate => (is => 'rw', default => 'outputs/16phanotate/phanotate.tsv.xz');
has input_prodigal => (is => 'rw', default => 'outputs/17prodigal/predicted_cds.gff');
has input_prokka_tsv => (is => 'rw', default => undef); ## Prokka tsv file for merging annotations.
has input_trinotate => (is => 'rw', default => '11trinotate_10prokka_09termreorder_08phageterm_07rosalind_plus/Trinotate.tsv'); ## trinotate output, used when merging annotations.
has input_umi => (is => 'rw', default => 'umi.txt');
has interactive => (is => 'rw', default => 0); ## Is this an interactive session?
has introns => (is => 'rw', default => 1); ## Is this method intron aware? (variant searching).
has iterate => (is => 'rw', default => undef);
has jobs => (is => 'rw', default => undef); ## List of currently active jobs, possibly not used right now.
has jobids => (is => 'rw', default => ''); ## A place to put running jobids, resurrected!
has jobnames => (is => 'rw', default => ''); ## A place to put running jobids, resurrected!
has jbasename => (is => 'rw', default => basename(cwd())); ## Job basename
has jcpu => (is => 'rw', default => 2); ## Number of processors to request in jobs
has jgpu => (is => 'rw', default => 0);
has jdepends => (is => 'rw', default => ''); ## Flag to start a dependency chain
has jmem => (is => 'rw', default => 20); ## Number of gigs of ram to request
has jname => (is => 'rw', default => undef); ## Job name on the cluster
has jnice => (is => 'rw', default => 10); ## Set the niceness of a job, if it starts positive, we can set a lower nice to preempt
has jpartition => (is => 'rw', default => 'dpart');
has jprefix => (is => 'rw', default => ''); ## Prefix number for the job
has jqueue => (is => 'rw', default => ''); ## What queue will jobs default to?
has jqueues => (is => 'rw', default => ''); ## Other possible queues
has jsleep => (is => 'rw', default => '0.5'); ## Set a sleep between jobs
has jstring => (is => 'rw', default => undef); ## String of the job
has jtemplate => (is => 'rw', default => undef);
has jwalltime => (is => 'rw', default => '10:00:00'); ## Default time to request
has keys => (is => 'rw', default => undef);
has kingdom => (is => 'rw', default => undef); ## Taxonomic kingdom, prokka/kraken
has language => (is => 'rw', default => 'bash'); ## What kind of script is this?
has last_job => (is => 'rw', default => ''); ## Last job in a chain.
has length => (is => 'rw', default => 17); ## kmer length, other stuff too.
has libdir => (is => 'rw', default => "\${HOME}/libraries"); ## Directory containing genomes/gff/indexes
has libpath => (is => 'rw', default => "$ENV{HOME}/libraries");
has library => (is => 'rw', default => undef); ## The library to be used for fasta36/blast searches
has libtype => (is => 'rw', default => 'genome'); ## Type of sequence to map against, genomic/rRNA/contaminants
has locus_tag => (is => 'rw', default => undef); ## Used by prokka to define gene prefixes
has logdir => (is => 'rw', default => 'outputs/logs'); ## place to dump logs
has loghost => (is => 'rw', default => 'localhost'); ## Host to which to send logs
has mapper => (is => 'rw', default => 'hisat:salmon'); ## Use this aligner if none was chosen.
has mature_fasta => (is => 'rw', default => undef); ## Database of mature miRNA sequences to search
has maximum => (is => 'rw', default => undef); ## catchall maximum threshold
has maxlength => (is => 'rw', default => 42); ## Maximum sequence length when trimming
has max_reads => (is => 'rw', default => 1000000);
has method => (is => 'rw', default => undef);
has mi_genome => (is => 'rw', default => undef); ## Set a miRbase genome to hunt for mature miRNAs
has min_depth => (is => 'rw', default => 5); ## Default use: variant searching, depth limit
has min_value => (is => 'rw', default => 0.5); ## Also variant searching.
has minimum => (is => 'rw', default => undef); ## catchall minimum threshold
has minlength => (is => 'rw', default => 8); ## Minimum length when trimming
has mirbase_data => (is => 'rw', default => undef); ## miRbase annotation dataset.
has mode => (is => 'rw', default => 'union');
has modulecmd => (is => 'rw', default => '');
has modules => (is => 'rw', default => undef); ## Environment modules to load
has module_string => (is => 'rw', default => '');
has numerator => (is => 'rw', default => undef);
has option_file => (is => 'rw', default => undef);
has orientation => (is => 'rw', default => 'start'); ## Default orientation when normalizing riboseq reads
has outformat => (is => 'rw', default => '.fasta');
has outgroup => (is => 'rw', default => undef); ## Outgroup for phylogenetic tools
has output => (is => 'rw', default => undef); ## Generic output argument
has output_base => (is => 'rw', default => 'outputs');
has outdir => (is => 'rw', default => undef);
has overlap => (is => 'rw', default => 20);
has overwrite => (is => 'rw', default => 1);
has paired => (is => 'rw', default => 1); ## Is the input paired?
has pdata => (is => 'rw', default => 'options.pdata');
has phred => (is => 'rw', default => 33); ## Minimum quality score when trimming
has postscript => (is => 'rw', default => undef); ## String to put after a cluter job.
has preprocess_dir => (is => 'rw', default => 'preprocessing');
has prescript => (is => 'rw', default => undef); ## String to put before a cluster job.
has primary_key => (is => 'rw', default => 'locus_tag'); ## Choose a keytype for merging data
has product_columns => (is => 'rw', default => 'trinity_sprot_Top_BLASTX_hit,inter_Pfam,inter_TIGRFAM'); ## When merging annotations, choose the favorites when upgrading an annotation to 'product'
has product_transmembrane => (is => 'rw', default => 'inter_TMHMM'); ## Column containing transmembrane domain data when upgrading annotations to 'product'
has product_signal => (is => 'rw', default => 'inter_signalp'); ## Column containing signal peptide domain data when upgrading annotations to 'product'
has protocol => (is => 'rw', default => 'Sassetti'); ## TNSeq protocol for TPP
has pval => (is => 'rw', default => undef); ## Pvalue cutoffs.
has qsub_args => (is => 'rw', default => '-j oe -V -m n'); ## What arguments will be passed to qsub by default?
has qsub_depends => (is => 'rw', default => 'depend=afterok:'); ## String to pass for dependencies
has qsub_dependsarray => (is => 'rw', default => 'depend=afterokarray:'); ## String to pass for an array of jobs
has qsub_path => (is => 'rw', default => scalar_which('qsub'));
has qual => (is => 'rw', default => undef); ## cutadapt quality string
has query => (is => 'rw', default => undef); ## Used for searches when input is already taken, most likely blast/fasta
has remove => (is => 'rw', default => 1);
has reference => (is => 'rw', default => undef);
has release => (is => 'rw', default => undef);
has requeue => (is => 'rw', default => 0); ## requeue failed jobs?
has restart => (is => 'rw', default => 0); ## Restart job(s) in the middle of a group
has riboanchor => (is => 'rw', default => 'start'); ## When correcting, use the start or end position as an anchor
has riboasite => (is => 'rw', default => 1); ## Count riboseq A site positions?
has ribopsite => (is => 'rw', default => 1); ## Count riboseq P site positions?
has riboesite => (is => 'rw', default => 1); ## Count riboseq E site positions?
has ribominsize => (is => 'rw', default => 24); ## Minimum size to search for ribosome positions (riboseq)
has ribomaxsize => (is => 'rw', default => 36); ## Maximum size to search for ribosome positions (riboseq)
has ribominpos => (is => 'rw', default => -30); ## Minimum position for counting (riboseq)
has ribomaxpos => (is => 'rw', default => 30); ## Maximum position for counting (riboseq)
has ribocorrect => (is => 'rw', default => 1); ## Correct ribosome positions for biases
has sizes => (is => 'rw', default => '25,26,27,28,29,30,31,32,33,34'); ## Use these sizes for riboseq reads
has runs => (is => 'rw', default => 1000); ## Number of runs for bayesian methods.
has sampleid => (is => 'rw', default => undef); ## Identifier to use for a sample.
has samtools_mapped => (is => 'rw', default => 0); ## Extract mapped reads with samtools.
has samtools_unmapped => (is => 'rw', default => 0); ## Extract unmapped reads with samtools.
has sbatch_depends => (is => 'rw', default => 'afterok:');
has sbatch_dependsarray => (is => 'rw', default => 'afterok:'); ## String to pass for an array of jobs
has sbatch_path => (is => 'rw', default => scalar_which('sbatch'));
has script_dir => (is => 'rw', default => 'scripts');
has search_string => (is => 'rw', default => 'tail');
has secondary => (is => 'rw', default => 'ignore');  ## Or ignore, use by htseq-count.
has separation => (is => 'rw', default => '10');
has seqid => (is => 'rw', default => undef);
has shell => (is => 'rw', default => '/usr/bin/bash'); ## Default qsub shell
has sort => (is => 'rw', default => 'pos'); ## Explicitly set the sorting type for htseq
has species => (is => 'rw', default => undef); ## Primarily for getting libraries to search against
has sra => (is => 'rw', default => 0);
has start => (is => 'rw', default => 1);
has starting_tree => (is => 'rw', default => undef); ## Starting tree for phylogenetic analyses
## Note 202212: Now most of the sequencing kits used by our sequencer are reverse.
has stranded => (is => 'rw', default => 'reverse'); ## Did this data come from a stranded library kit?
has suffixes => (is => 'rw', default => '.bam,.bigwig,.bw,.count,.csfasta,.csv,.fa,.faa,.fasta,.fastq,.ffn,.fna,.fsa,.gb,.gba,.gbf,.gbk,.genbank,.gff,.gz,.qual,.sam,.sf,.tbl,.tsv,.xz'); ## Suffixes to remove when invoking basename
has supplementary => (is => 'rw', default => 'ignore'); ## or ignore, used by htseq-count.
has ta_offset => (is => 'rw', default => 0); ## When counting TAs, this is either 0 or 2 depending on if the TA was removed.
has task => (is => 'rw', default => 'unknown');
has taxid => (is => 'rw', default => '353153'); ## Default taxonomy ID, unknown for now.
has test_file => (is => 'rw', default => 'direct-term-repeasts.fasta'); ## There are a few places where testing for the existence of a test file is useful.
has threshold => (is => 'rw', default => 0.05); ## A second cutoff value (looking at you, glimmer.)
has translate => (is => 'rw', default => 1);
has trim => (is => 'rw', default => 1); ## Perform trimming (rnaseq pipeline, trinity)
has type => (is => 'rw', default => undef); ## Possibly superceded by htseq_type
has varfilter => (is => 'rw', default => 1); ## use a varfilter when performing variant searches.
has verbose => (is => 'rw', default => 0); ## Print extra information while running?
has version => (is => 'rw', default => undef);
has vcf_cutoff => (is => 'rw', default => 5); ## Minimum depth cutoff for variant searches
has vcf_method => (is => 'rw', default => 'freebayes');
has vcf_minpct => (is => 'rw', default => 0.8); ## Minimum percent agreement for variant searches.
has write => (is => 'rw', default => 1);        ## Write outputs?
## A few variables which are by definition hash references and such
has methods_to_run => (is => 'rw', default => undef); ## Set of jobs to run.
has menus => (is => 'rw', default => undef); ## The menus when in an interactive session.
has term => (is => 'rw', default => undef); ## A fun Readline terminal for tab completion.
has todos => (is => 'rw', default => undef); ## Set of TODOs to perform.
has umi => (is => 'rw', default => 0);

1;
