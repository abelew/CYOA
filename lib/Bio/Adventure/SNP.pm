package Bio::Adventure::SNP;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
no warnings 'redefine';
extends 'Bio::Adventure';
use feature 'try';
no warnings 'experimental::try';

use Acme::Tools qw"btw";
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Which qw"which";
use IO::String;
use List::Util qw"sum";
use Math::Round qw":all";
use POSIX qw"floor";
use Bio::DB::SeqFeature::Store;
use Bio::Matrix::IO;
use Bio::Matrix::Scoring;
use Bio::SeqIO;
use Bio::Seq;

=head1 NAME

 Bio::Adventure::SNP - Search for variant positions given an alignment and reference genome.

=head1 SYNOPSIS

 The functions in this file handle the invocation of b/vcftools and mpileup to
 search for variant positions following an alignment.

=head1 METHODS

=head2 C<Align_SNP_Search>

 Invoke bt2, samtools, vcfutils to seek variant positions.

=cut
sub Align_SNP_Search {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],);
    my $genome = qq"$options->{libpath}/$options->{libtype}/$options->{species}.fasta";
    my $query = $options->{input};
    my $query_home = dirname(${query});
    my $query_base = basename(${query}, (".bam"));
    $query = qq"${query_home}/${query_base}";
    print "About to start Bowtie2 search against of ${query} against $options->{species}.\n";
    my $bt2_job = $class->Bio::Adventure::Map::Bowtie2(
        gff_type => 'exon',
        input => $query,
        species => $options->{species},
        jprefix => qq"$options->{jprefix}_1",);
    my $bamfile = $bt2_job->{samtools}->{output};
    print "About to start SNP search of ${bamfile} against $options->{species}\n";
    my $search = $class->Bio::Adventure::SNP::Freebayes_SNP_Search(
        input => $bamfile,
        jdepends => $bt2_job->{samtools}->{job_id},
        jprefix => qq"$options->{jprefix}_2",);
    return($search);
}

sub Dantools_RNASeq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        chosen_tag => 'PAIRED',
        coverage_tag => 'DP',
        gff_tag => 'ID',
        gff_type => 'gene',
        gff_cds_parent_type => 'mRNA',
        gff_cds_type => 'CDS',
        introns => 1,
        qual => 10,
        max_value => undef,
        min_value => 0.8,
        vcf_cutoff => 5,
        jmem => 36,
        jcpu => 4,
        jprefix => '50',
        jwalltime => '48:00:00',);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $input_fasta = $paths->{fasta};
    my $dantools_dir = $paths->{output_dir};
    my $stdout = qq"${dantools_dir}/$options->{species}.stdout";
    my $stderr = qq"${dantools_dir}/$options->{species}.stderr";

    my $dt_input = $options->{input};
    my $dt_input_args = '';
    my $input_name = $dt_input;
    if ($dt_input =~ /$options->{delimiter}/) {
        my @pair_listing = split(/$options->{delimiter}/, $dt_input);
        $dt_input_args .= qq" -1 $pair_listing[0] -2 $pair_listing[1] ";
        $input_name = $pair_listing[0];
    } else {
        $dt_input_args .= qq" --reads-u ${dt_input} ";
    }
    my $input_base = basename($dt_input, ('.xz', '.gz', '.bz2'));
    $input_base = basename($dt_input, ('.fastq'));
    my $shift_input_filename = qq"$paths->{output_dir}/outputs/${input_base}_on_$options->{species}";
    my $comment = qq!## This is a dantools search for variant against $options->{species}!;
    my $jstring = qq!
mkdir -p ${dantools_dir}
dantools pseudogen \\
  --outdir $paths->{output_dir} \\
  -r ${input_fasta} \\
  ${dt_input_args} \\
  1>>$paths->{output_dir}/pseudogen.stdout \\
  2>>$paths->{output_dir}/pseudogen.stderr
dantools shift \\
  -o $paths->{output_dir}/shifted.gff \\
  -v $shift_input_filename \\
  -f $paths->{gff} \\
  1>>$paths->{output_dir}/shift.stdout \\
  2>>$paths->{output_dir}/shift.stderr
dantools summarize-depth --outdir $paths->{output_dir} \\
  -f shifted.gff -d depth.tsv --feature gene \\
  1>>$paths->{output_dir}/summarize.stdout \\
  2>>$paths->{output_dir}/summarize.stderr
dantools label --outdir $paths->{output_dir} \\
  -v variants.vcf -f $paths->{gff} \\
  --features five_prime_UTR,CDS,three_prime_UTR \\
  1>>$paths->{output_dir}/label.stdout \\
  2>>$paths->{output_dir}/label.stderr
dantools summarize-nuc --outdir $paths->{output_dir} \\
  labeled_nucleotides.tsv \\
  1>>$paths->{output_dir}/summarize_nuc.stdout \\
  2>>$paths->{output_dir}/summarize_nuc.stderr
dantools summarize-aa --outdir $paths->{output_dir} \\
  labeled_aaleotides.tsv \\
  1>>$paths->{output_dir}/summarize_aa.stdout \\
  2>>$paths->{output_dir}/summarize_aa.stderr
!;

    my $comment_string = '## Use dantools to examine RNASeq data vs a reference.';
    my $dantools = $class->Submit(
        comment => $comment_string,
        jdepends => $options->{jdepends},
        jname => "dantools_$options->{species}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        stderr => $stderr,
        stdout => $stdout,);
    return($dantools);
}

=head2 C<Freebayes_SNP_Search>

 Invoke freebayes to create and filter a set of variant positions
 after using GATK to limit the duplicate sequence effect.
 GATK: 10.1186/s40104-019-0359-0
 freebayes: arXiv:1207.3907

 cyoa --method freebayes --species lpanamensis_v68 \
      --input outputs/hisat_lpanamensis_v68/lpanamensis-paired.bam

=cut
sub Freebayes_SNP_Search {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        chosen_tag => 'PAIRED',
        coverage_tag => 'DP',
        gff_tag => 'ID',
        gff_type => 'gene',
        gff_cds_parent_type => 'mRNA',
        gff_cds_type => 'CDS',
        introns => 1,
        qual => 10,
        max_value => undef,
        min_value => 0.8,
        vcf_cutoff => 5,
        jmem => 36,
        jcpu => 4,
        jprefix => '50',
        jwalltime => '48:00:00',);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $input_fasta = $paths->{fasta};
    my $freebayes_dir = $paths->{output_dir};
    my $output_file = qq"${freebayes_dir}/$options->{species}.vcf";
    my $output_bcf = qq"${freebayes_dir}/$options->{species}.bcf";
    my $marked = qq"${freebayes_dir}/deduplication_stats.txt";
    my $gatk_stdout = qq"${freebayes_dir}/deduplication.stdout";
    my $gatk_stderr = qq"${freebayes_dir}/deduplication.stderr";
    my $stdout = qq"${freebayes_dir}/$options->{species}.stdout";
    my $stderr = qq"${freebayes_dir}/$options->{species}.stderr";

    ## Taken from the sister mpileup function
    my $query_dir = dirname($options->{input});
    my $query_base = basename($options->{input}, ('.bam'));
    my $deduplicate_input = qq"${freebayes_dir}/${query_base}_sorted.bam";
    my $deduplicated = qq"${freebayes_dir}/${query_base}_deduplicated.bam";

    my $comment = qq!## This is a freebayes search for variant against ${input_fasta}!;
    my $jstring = qq!
mkdir -p ${freebayes_dir}
samtools sort -l 9 -@ 4 $options->{input} -o ${deduplicate_input} \\
  2>${freebayes_dir}/samtools_sort.out 1>&2
if [ "\$?" -ne "0" ]; then
    echo "samtools sort failed."
fi
gatk MarkDuplicates \\
  -I ${deduplicate_input} \\
  -O ${deduplicated} \\
  -M ${marked} --REMOVE_DUPLICATES true --COMPRESSION_LEVEL 9 \\
  2>${gatk_stderr} \\
  1>${gatk_stdout}
samtools index ${deduplicated}
echo "Finished samtools index." >> ${stdout}
freebayes -f ${input_fasta} \\
  -v ${output_file} \\
  ${deduplicated} \\
  1>>${stdout} \\
  2>>${stderr}
echo "Finished freebayes." >> ${stdout}
bcftools convert ${output_file} \\
  -Ob -o ${output_bcf} \\
  2>>${stderr} \\
  1>>${stdout}
echo "Finished bcftools convert." >> ${stdout}
bcftools index ${output_bcf} \\
  2>>${stderr} \\
  1>>${stdout}
echo "Finished bcftools index." >> ${stdout}
rm ${output_file}
!;

    my $comment_string = '## Use freebayes, bcftools, and vcfutils to get some idea about how many variant positions are in the data.';
    my $freebayes = $class->Submit(
        comment => $comment_string,
        jdepends => $options->{jdepends},
        jname => "freebayes_$options->{species}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        stderr => $stderr,
        stdout => $stdout,
        output => $output_bcf,);

    $comment_string = qq!## This little job should make unique IDs for every detected
## SNP and a ratio of snp/total for all snp positions with > 20 reads.
## Further customization may follow.
!;
    my $parse;
    if ($options->{introns}) {
        $parse = $class->Bio::Adventure::SNP::SNP_Ratio_Intron(
            chosen_tag => $options->{chosen_tag},
            coverage_tag => $options->{coverage_tag},
            gff_tag => $options->{gff_tag},
            gff_type => $options->{gff_type},
            gff_cds_parent_type => $options->{gff_cds_parent_type},
            gff_cds_type => $options->{gff_cds_type},
            input => $output_bcf,
            jcpu => 1,
            jdepends => $freebayes->{job_id},
            jname => qq"freebayes_parsenp_intron_${query_base}",
            jprefix => qq"$options->{jprefix}_1",
            min_value => $options->{min_value},
            max_value => $options->{max_value},
            qual => $options->{qual},
            species => $options->{species},
            vcf_cutoff => $options->{vcf_cutoff},
            vcf_method => 'freebayes',
            );
    } else {
        $parse = $class->Bio::Adventure::SNP::SNP_Ratio(
            chosen_tag => $options->{chosen_tag},
            coverage_tag => $options->{coverage_tag},
            gff_tag => $options->{gff_tag},
            gff_type => $options->{gff_type},
            input => $output_bcf,
            max_value => $options->{max_value},
            min_value => $options->{min_value},
            qual => $options->{qual},
            species => $options->{species},
            vcf_cutoff => $options->{vcf_cutoff},
            vcf_method => 'freebayes',
            jcpu => 1,
            jdepends => $freebayes->{job_id},
            jname => qq"freebayes_parsenp_${query_base}",
            jprefix => qq"$options->{jprefix}_1",);
    }
    $freebayes->{parse} = $parse;
    return($freebayes);
}

=head2 C<SNP_Search>

 Handle the invocation of vcfutils and such to seek high-confidence variants.

=cut
sub Mpileup_SNP_Search {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species', 'gff_tag', 'gff_type'],
        varfilter => 0,
        vcf_cutoff => 5,
        min_value => 0.8,
        qual => 10,
        jprefix => '50',);
    my $genome = qq"$options->{libpath}/$options->{libtype}/fasta/$options->{species}.fasta";
    my $query = $options->{input};
    my $query_home = dirname(${query});
    my $query_base = basename(${query}, ('.bam'));
    $query = qq"${query_home}/${query_base}";

    my $varfilter = $options->{varfilter};
    my $vcf_cutoff = $options->{vcf_cutoff};
    my $vcf_minpct = $options->{min_value};

    my $vcfutils_dir = qq"outputs/$options->{jprefix}vcfutils_$options->{species}";
    my $deduplicate_input = qq"${vcfutils_dir}/${query_base}.bam";
    my $pileup_input = qq"${vcfutils_dir}/${query_base}_deduplicated.bam";
    my $pileup_error = qq"${vcfutils_dir}/${query_base}_pileup.stderr";
    my $pileup_output = qq"${vcfutils_dir}/${query_base}_pileup.vcf";
    my $filter_error = qq"${vcfutils_dir}/${query_base}_varfilter.stderr";
    my $filter_output = qq"${vcfutils_dir}/${query_base}_filtered.vcf";
    my $call_output = qq"${vcfutils_dir}/${query_base}_summary.vcf";
    my $call_error = qq"${vcfutils_dir}/${query_base}_summary.stderr";
    my $final_output = qq"${vcfutils_dir}/${query_base}.bcf";
    my $final_error = qq"${vcfutils_dir}/${query_base}_bcf.stderr";
    my $gatk_stdout = qq"${vcfutils_dir}/deduplication.stdout";
    my $gatk_stderr = qq"${vcfutils_dir}/deduplication.stderr";
    my $marked = qq"${vcfutils_dir}/deduplication_stats.txt";
    my $jstring = qq!mkdir -p ${vcfutils_dir}
echo "Started samtools sort at \$(date)" >> ${vcfutils_dir}/vcfutils_$options->{species}.stdout
!;
    unless (-r "${pileup_input}") {
        if ($query =~ m/sorted/) {
            $jstring .= qq!if test \! -e \$(pwd)/${pileup_input}; then
  ln -sf \$(pwd)/${query}.bam \$(pwd)/${deduplicate_input}
fi
!;
        } else {
            $jstring .= qq!
samtools sort -l 9 -@ 4 ${query}.bam -o ${deduplicate_input} \\
  2>${vcfutils_dir}/samtools_sort.out 1>&2
if [ "\$?" -ne "0" ]; then
    echo "samtools sort failed."
fi
gatk MarkDuplicates \\
  -I ${deduplicate_input} \\
  -O ${pileup_input} \\
  -M ${marked} --REMOVE_DUPLICATES true --COMPRESSION_LEVEL 9 \\
  2>${gatk_stderr} \\
  1>${gatk_stdout}
samtools index ${pileup_input}
!;
        } ## End checking if the pileup input does not have sorted.
    }     ## Found the input for samtools mpileup
    $jstring .= qq!
if [ \! -r "${genome}.fai" ]; then
    samtools faidx ${genome}
fi
samtools mpileup -uvf ${genome} 2>${vcfutils_dir}/samtools_mpileup.err \\
    ${pileup_input} |\\
  bcftools call -c - 2>${vcfutils_dir}/bcftools_call.err |\\
  bcftools view -l 9 -o ${final_output} -Ob - \\
    2>${call_error}
if [ "\$?" -ne "0" ]; then
    echo "mpileup/bcftools failed."
    exit 1
fi
bcftools index ${final_output} 2>${vcfutils_dir}/bcftools_index.err
echo "Successfully finished." >> ${vcfutils_dir}/vcfutils_$options->{species}.out
!;

    my $comment_string = qq"## Use samtools, bcftools, and vcfutils to get some
## idea about how many variant positions are in the data.";
    my $pileup = $class->Submit(
        comment => $comment_string,
        input_pileup => $pileup_input,
        jdepends => $options->{jdepends},
        jcpu => 4,
        jmem => 48,
        jname => qq"mpileup_${query_base}",
        jprefix => qq"$options->{jprefix}",
        jstring => $jstring,
        jwalltime => '10:00:00',
        output_call => $call_output,
        output_filter => $filter_output,
        output_final => $final_output,
        output_pileup => $pileup_output,
        stderr => $gatk_stderr,
        stdout => $gatk_stdout,);
    $comment_string = qq!## This little job should make unique IDs for every detected
## SNP and a ratio of snp/total for all snp positions with > 20 reads.
## Further customization may follow.
!;
    my $parse = $class->Bio::Adventure::SNP::SNP_Ratio(
        gff_tag => $options->{gff_tag},
        gff_type => $options->{gff_type},
        input => ${final_output},
        qual => $options->{qual},
        species => $options->{species},
        vcf_cutoff => $vcf_cutoff,
        vcf_minpct => $vcf_minpct,
        jcpu => 1,
        jname => qq"mpileup_parsenp_${query_base}_$options->{species}",
        jdepends => $pileup->{job_id},
        jprefix => $options->{jprefix} + 1,);
    $pileup->{parse} = $parse;
    return($pileup);
}

=head2 C<CNP_Ratio>

 Given a table of variants and a genome, modify the genome to match the variants.

=cut
sub SNP_Ratio {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        introns => 0,
        qual => 10,
        vcf_method => 'freebayes',
        vcf_cutoff => 5,
        vcf_minpct => 0.8,
        jprefix => '80',
        jname => 'parsenp',
        required => ['input', 'species', 'gff_tag', 'gff_type'],);
    if ($options->{introns}) {
        my $snp_intron = Bio::Adventure::SNP::SNP_Ratio_Intron(%args);
        return($snp_intron);
    }
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $print_input = $options->{input};
    my $output_dir = dirname($print_input);
    my $print_output = qq"${output_dir}";
    my $genome = $paths->{fasta};
    my $samplename = basename(cwd());

    my $output_suffix = '';
    $output_suffix .= qq"_q-$options->{qual}" if ($options->{qual});
    $output_suffix .= qq"_c-$options->{vcf_cutoff}" if ($options->{vcf_cutoff});
    $output_suffix .= qq"_m$options->{min_value}" if ($options->{min_value});
    $output_suffix .= qq"_ctag-$options->{coverage_tag}" if ($options->{coverage_tag});
    $output_suffix .= qq"_mtag-$options->{chosen_tag}" if ($options->{chosen_tag});
    my $stdout = qq"${print_output}/stdout";
    my $stderr = qq"${print_output}/stderr";
    my $output_all = qq"${print_output}/all_tags${output_suffix}.txt.xz";
    my $output_genome = qq"${print_output}/$options->{species}-${samplename}${output_suffix}.fasta";
    my $output_by_gene = qq"${print_output}/variants_by_gene${output_suffix}.txt.xz";
    my $output_penetrance = qq"${print_output}/variants_penetrance${output_suffix}.txt.xz";
    my $output_pkm = qq"${print_output}/pkm${output_suffix}.txt.xz";
    my $output_types = qq"${print_output}/count_types${output_suffix}.txt";
    my $output_dedup = qq"${print_output}/deduplication_stats.txt";
    my $comment_string = qq!
## Parse the SNP data and generate a modified $options->{species} genome.
##  This should read the file:
## ${print_input}
##  and provide some new files:
## ${output_genome}
## ${output_by_gene}
## ${output_penetrance}
## ${output_pkm}
!;
    my $jstring = qq"
use Bio::Adventure::SNP;
my \$result = \$h->Bio::Adventure::SNP::SNP_Ratio_Worker(
  input => '$print_input',
  species => '$options->{species}',
  vcf_method => '$options->{vcf_method}',
  vcf_cutoff => '$options->{vcf_cutoff}',
  vcf_minpct => '$options->{vcf_minpct}',
  gff_tag => '$options->{gff_tag}',
  gff_type => '$options->{gff_type}',
  qual => '$options->{qual}',
  output_dir => '$output_dir',
  output => '${output_all}',
  output_dedup => '${output_dedup}',
  output_genome => '${output_genome}',
  output_by_gene => '${output_by_gene}',
  output_penetrance => '${output_penetrance}',
  output_pkm => '${output_pkm}',
  output_types => '${output_types}',
);
";
    my $parse_job = $class->Submit(
        comment => $comment_string,
        jdepends => $options->{jdepends},
        jmem => 48,
        jname => $options->{jname},
        jprefix => qq"$options->{jprefix}_1",
        jstring => $jstring,
        jwalltime => '10:00:00',
        language => 'perl',
        output_dir => $output_dir,
        output => $output_all,
        output_dedup => $output_dedup,
        output_genome => $output_genome,
        output_by_gene => $output_by_gene,
        output_penetrance => $output_penetrance,
        output_pkm => $output_pkm,
        qual => $options->{qual},
        stdout => $stdout,
        stderr => $stderr,);
    return($parse_job);
}

=head2 C<Make_SNP_Ratio>

 Given vcfutils output, make a simplified table of high-confidence variants.

 This function has been made more generic with my recent use of freebayes.
 As a result, it should probably be renamed and split apart so that the
 generation of a new haplotype/genome is separate.

=cut
sub SNP_Ratio_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        chosen_tag => 'PAIRED',
        coverage_tag => 'DP',
        gff_tag => 'ID',
        gff_type => 'gene',
        min_value => 0.8,
        max_value => undef,
        only_point => 1,
        output => 'all.txt',
        output_genome => 'new_genome.fasta',
        output_by_gene => 'counts_by_gene.txt',
        output_penetrance => 'variants_penetrance.txt',
        output_pkm => 'by_gene_length.txt',
        output_dir => 'outputs/40freebayes',
        penetrance_tag => 'SAP',
        qual => 10,
        vcf_cutoff => 5,
        vcf_method => 'freebayes',);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $species = $options->{species};
    my $out_dir = dirname($options->{output});
    my $log_file = qq"${out_dir}/snp_ratio.stdout";
    my $log = FileHandle->new(">${log_file}");
    my $subst_matrix = $class->Bio::Adventure::Align::Get_Substitution_Matrix(matrix => $options->{matrix});
    my $genome = $paths->{fasta};
    my $gff = $paths->{gff};
    print $log "Reading gff: ${gff}, extracting type: $options->{gff_type} features tagged $options->{gff_tag}.\n";
    my $in_bcf = FileHandle->new("bcftools view $options->{input} |");
    print $log "The large matrix of data will be written to: $options->{output}\n";
    my $all_out = FileHandle->new(">$options->{output}");
    my $filtered_coverage_log = qq"${out_dir}/filtered_coverage.tsv";

    ## Ok, so I want to simplify the vcf output so that I can create a pseudo
    ## count table of every potential SNP position I want to create a 2 column
    ## file where the first column is a unique ID containing 'chr_pos_ref_new'
    ## and the second column is some sort of score which I will be able to use
    ## for a PCA when I merge a bunch of these together.  Candidates include:
    ## quality score, DP*AF1 (depth * max likelihood estimate) Maybe also consider
    ## DP4 column which is comma separated: #forward-ref, #reverse-ref, #forward-alt, #reverse-alt
    ## If that, then sum all 4 and only take those greater than x (20?), then take
    ## forward-alt+reverse-alt/(sum of all) to get simple SNP ratio
    ## I am going to change my writer of the all.txt file so that it prints out the
    ## values of every observed tag.  For now I am just going to hard-code the order,
    ## but it should not be difficult to parse this out of the header lines of the
    ## bcf file.
    my @mpileup_tag_order = (
        'DP', 'ADF', 'AD', 'VDB', 'SGB', 'MQ0F', 'RPB', 'MQB', 'BQB', 'MQSB', 'ADR', 'GT',
        'ICB', 'HOB', 'AC', 'AN', 'DP4', 'MQ',);
    my @freebayes_tag_order = (
        'NS', 'DP', 'DPB', 'AC', 'AN', 'AF', 'RO', 'AO',
        'PRO', 'PAO', 'QR', 'QA', 'PQR', 'PQA', 'SRF', 'SRR', 'SAF',
        'SAR', 'SRP', 'SAP', 'AB', 'ABP', 'RUN', 'RPP', 'RPPR', 'RPL',
        'RPR', 'EPP', 'EPPR', 'DPRA', 'ODDS', 'GTI', 'TYPE', 'CIGAR',
        'NUMALT', 'LEN', 'MQM', 'MQMR', 'PAIRED', 'PAIREDR', 'MIN_DP',
        'END', 'GT', 'GQ', 'GL', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA',
        'MIN_DP');
    ## Make a data structure containing arrays for the set of observed tags in the vcf file.
    ## It should be an array which ends at length of the number of bcf entries.
    my @tag_order = ();
    if ($options->{vcf_method} eq 'freebayes') {
        @tag_order = @freebayes_tag_order;
    } else {
        @tag_order = @mpileup_tag_order;
    }
    ## Here is the global structure
    my @tag_observations = ();
    ## Use this little hash to decide what tags to write to the final tsv file.
    my %num_observed = ();
    for my $k (@tag_order) {
        $num_observed{$k} = 0;
    }
    ## Then fill it with a hash of the tags from @tag_order
    ## The following hash will be over written.

    ## There is a problem with the above, when I wrote it I assumed that all tools which
    ## boil down to a bcf file would have the same tags describing the quality of the
    ## observation.  That is a bad assumption, here are the tags from freebayes:
    ## NS,Number=1,Type=Integer,Description="Number of samples with data">
    ## DP,Number=1,Type=Integer,Description="Total read depth at the locus">
    ## DPB,Number=1,Type=Float,Description="Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype">
    ## AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
    ## AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ## AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
    ## RO,Number=1,Type=Integer,Description="Count of full observations of the reference haplotype.">
    ## AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
    ## PRO,Number=1,Type=Float,Description="Reference allele observation count, with partial observations recorded fractionally">
    ## PAO,Number=A,Type=Float,Description="Alternate allele observations, with partial observations recorded fractionally">
    ## QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">
    ## QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">
    ## PQR,Number=1,Type=Float,Description="Reference allele quality sum in phred for partial observations">
    ## PQA,Number=A,Type=Float,Description="Alternate allele quality sum in phred for partial observations">
    ## SRF,Number=1,Type=Integer,Description="Number of reference observations on the forward strand">
    ## SRR,Number=1,Type=Integer,Description="Number of reference observations on the reverse strand">
    ## SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
    ## SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
    ## SRP,Number=1,Type=Float,Description="Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding's inequality">
    ## SAP,Number=A,Type=Float,Description="Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding's inequality">
    ## AB,Number=A,Type=Float,Description="Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous">
    ## ABP,Number=A,Type=Float,Description="Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality">
    ## RUN,Number=A,Type=Integer,Description="Run length: the number of consecutive repeats of the alternate allele in the reference genome">
    ## RPP,Number=A,Type=Float,Description="Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality">
    ## RPPR,Number=1,Type=Float,Description="Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality">
    ## RPL,Number=A,Type=Float,Description="Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele">
    ## RPR,Number=A,Type=Float,Description="Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele">
    ## EPP,Number=A,Type=Float,Description="End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality">
    ## EPPR,Number=1,Type=Float,Description="End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality">
    ## DPRA,Number=A,Type=Float,Description="Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.">
    ## ODDS,Number=1,Type=Float,Description="The log odds ratio of the best genotype combination to the second-best.">
    ## GTI,Number=1,Type=Integer,Description="Number of genotyping iterations required to reach convergence or bailout.">
    ## TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
    ## CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">
    ## NUMALT,Number=1,Type=Integer,Description="Number of unique non-reference alleles in called genotypes at this position.">
    ## ,Number=A,Type=Float,Description="Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.">
    ## LEN,Number=A,Type=Integer,Description="allele length">
    ## MQM,Number=A,Type=Float,Description="Mean mapping quality of observed alternate alleles">
    ## MQMR,Number=1,Type=Float,Description="Mean mapping quality of observed reference alleles">
    ## PAIRED,Number=A,Type=Float,Description="Proportion of observed alternate alleles which are supported by properly paired read fragments">
    ## PAIREDR,Number=1,Type=Float,Description="Proportion of observed reference alleles which are supported by properly paired read fragments">
    ## MIN_DP,Number=1,Type=Integer,Description="Minimum depth in gVCF output block.">
    ## END,Number=1,Type=Integer,Description="Last position (inclusive) in gVCF output record.">
    ## GT,Number=1,Type=String,Description="Genotype">
    ## GQ,Number=1,Type=Float,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">
    ## GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
    ## DP,Number=1,Type=Integer,Description="Read Depth">
    ## AD,Number=R,Type=Integer,Description="Number of observation for each allele">
    ## RO,Number=1,Type=Integer,Description="Reference allele observation count">
    ## QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
    ## AO,Number=A,Type=Integer,Description="Alternate allele observation count">
    ## QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
    ## MIN_DP,Number=1,Type=Integer,Description="Minimum depth in gVCF output block.">

    ## As should therefore be obvious, I will need to add some logic/thought to how I handle
    ## these various tags.  The good news is that they have some nice scores referring to the
    ## likelihood of heterozygosity, which is precisely what I want to quantify in my current
    ## round of analysis

    ## In contrast, here are the tags from samtools' mpileup (btw, this is just copy/pasted
    ## from the vcf header):

    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
    ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
    ##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
    ##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
    ##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
    ##INFO=<ID=ADR,Number=R,Type=Integer,Description="Total allelic depths on the reverse strand">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
    ##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
    ##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
    ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
    ##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
    ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
    ##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">

    ## If I recall properly, the one of primary interest for looking at questions of penetrance is: DP4.

    ## Here is an arbitrary vcf line as a refresher to make sure I grabbed them all:
    ## NC_045512.2     3037    .       C       T       224     PASS    DP=244;ADF=0,121;ADR=0,121;AD=0,242;VDB=0;SGB=-0.693147;MQSB=0.976205;MQ0F=0;AC=2;AN=2;DP4=0,0,121,121;MQ=42    GT:PL:DP:SP:ADF:ADR:AD  1/1:254,255,0:242:0:0,121:0,121:0,242
    print $log "Reading bcf file.\n";
    my $count = 0; ## Use this to count the positions changed in the genome.
    my $num_variants = 0; ## Use this to count the lines read in the bcf file.
  READER: while (my $line = <$in_bcf>) {
      $num_variants++;
      next READER if ($line =~ /^#/);
      ## When using samtools mpileup with the -t LIST, we get some extra
      ## fields at the end which are called tags and tag_info here.
      my ($chr, $pos, $id, $ref, $alt, $qual, $filt, $info, $tags, $tag_info) = split(/\s+/, $line);
      ## Only accept the simplest non-indel mutations.
      if ($options->{qual} ne '' && defined($qual)) {
          if ($qual < $options->{qual}) {
              next READER;
          }
      }
      ## I think something is wrong in my complex substitution/indel logic.
      if ($options->{only_point}) {
          ## Only accept the simplest non-indel mutations.
          next READER unless ($alt eq 'A' or $alt eq 'a' or $alt eq 'T' or $alt eq 't' or
                              $alt eq 'G' or $alt eq 'g' or $alt eq 'C' or $alt eq 'c');
      }

      ## I do not know if indels are here anymore, we will see.
      my @info_list = split(/;/, $info);
      my $snp_id = qq"chr_${chr}_pos_${pos}_ref_${ref}_alt_${alt}";
      my $snp_pct = 0.0;
      my $all_sum = 0;
      my $diff_sum = 0;
      my %individual_tags = (position => $snp_id);
    TAGS: foreach my $element (@info_list) {
        my ($element_type, $element_value) = split(/=/, $element);

        ## At this point, add some post-processing for tags which are multi-element.
        if ($options->{vcf_method} eq 'mpileup' && $element_type eq 'DP4') {
            my ($same_forward, $same_reverse, $alt_forward, $alt_reverse) = split(/,/, $element_value);
            $diff_sum = $alt_forward + $alt_reverse;
            $all_sum = $same_forward + $same_reverse + $diff_sum;
            my $snp_pct = -1;
            if (defined($all_sum) && $all_sum > 0) {
                $snp_pct = ($alt_forward + $alt_reverse) / $all_sum;
            }
            $snp_pct = nearest(0.01, $snp_pct);
            $individual_tags{snp_pct} = $snp_pct;
        }
        $individual_tags{$element_type} = $element_value;
        if (!defined($num_observed{$element_type})) {
            $num_observed{$element_type} = 1;
            push(@tag_order, $element_type);
        } else {
            $num_observed{$element_type}++;
        }
    } ## End iterating over the tags in the data.
      push(@tag_observations, \%individual_tags);
  } ## End reading the bcf file.
    print $log "Finished reading ${num_variants} entries in the bcf file.\n";
    $in_bcf->close();
    ## Now we should have a big array full of little hashes
    my @used_tags = ('position', );
    ## First collect the set of tags of interest.
    my $header_line = "position\t";
    for my $k (@tag_order) {
        if ($num_observed{$k} > 0) {
            push(@used_tags, $k);
            $header_line .= "$k\t";
        }
    }
    $header_line =~ s/\t$/\n/;
    ## Write out the observations, starting with a header line comprised of the position information
    ## followed by the tags.
    print $all_out qq"${header_line}\n";
    print $log "Reading genome.\n";
    my $input_genome = $class->Bio::Adventure::Read_Genome_Fasta(
        %args, fasta => $genome,);
    my $annotations = $class->Bio::Adventure::Read_Genome_GFF(
        gff => $gff, gff_tag => $options->{gff_tag},
        gff_type => $options->{gff_type}, %args);
    my $output_by_gene = FileHandle->new(">$options->{output_by_gene}");
    print $output_by_gene qq"gene\tchromosome\tposition\tfrom_to\taa_subst\tsynonymousp\tblosum_delta\n";
    my $output_penetrance = FileHandle->new(">$options->{output_penetrance}");
    print $output_penetrance qq"chromosome\tposition\tfrom_to\tpenetrance\n";
    my $filtered_coverage = FileHandle->new(">${filtered_coverage_log}");
    my $vars_by_gene = {};
    my $all_count = 0;
    ## Make a copy of the genome so that we can compare before/after the mutation(s)
    my %original_genome;
    for my $chr_key (keys %{$input_genome}) {
        my %internal = %{$input_genome->{$chr_key}};
        $original_genome{$chr_key} = \%internal;
    }

    print $log "Filtering variant observations and writing a new genome.\n";
    my $shifters = 0;
    my $points = 0; ## Number of point mutations observed.
    my $total_delta = 0;
  SHIFTER: while (scalar(@tag_observations)) {
      my $datum = shift(@tag_observations);
      $shifters++;
      ## Now decide if we actually want to write this difference into
      ## a new genome and record it
      my $pen_tag = $options->{penetrance_tag};
      my $penetrance = $datum->{$pen_tag};
      my $chosen = $options->{chosen_tag};
      my $depth_tag = $options->{coverage_tag};
      my $vcf_cutoff = $options->{vcf_cutoff};
      if (defined($vcf_cutoff)) {
          if ($datum->{$depth_tag} < $vcf_cutoff) {
              print "Dropping $datum->{position} because depth ($datum->{$depth_tag}) is less than ${vcf_cutoff}.\n";
              print $filtered_coverage "$datum->{position}\t$datum->{$depth_tag}\t$datum->{$chosen}\n";
              next SHIFTER;
          }
      }
      if (defined($chosen) && defined($datum->{chosen})) {
          if (defined($options->{min_value}) && defined($options->{max_value})) {
              if ($datum->{chosen} > $options->{max_value}) {
                  print "$datum->{position} has exceeded max value: $options->{max_value}\n";
                  print $filtered_coverage "$datum->{position}\t$datum->{$depth_tag}\t$datum->{chosen}\n";
                  next SHIFTER;
              } elsif ($datum->{chosen} < $options->{min_value}) {
                  print "$datum->{position} is less than min value: $options->{min_value}\n";
                  print $filtered_coverage "$datum->{position}\t$datum->{$depth_tag}\t$datum->{chosen}\n";
                  next SHIFTER;
              }
          } elsif (defined($options->{min_value})) {
              if ($datum->{chosen} < $options->{min_value}) {
                  print "$datum->{position} is less than min value: $options->{min_value}\n";
                  print $filtered_coverage "$datum->{position}\t$datum->{$depth_tag}\t$datum->{chosen}\n";
                  next SHIFTER;
              }
          } elsif (defined($options->{max_value})) {
              if ($datum->{chosen} > $options->{max_value}) {
                  print "$datum->{position} has exceeded max value: $options->{max_value}\n";
                  print $filtered_coverage "$datum->{position}\t$datum->{$depth_tag}\t$datum->{chosen}\n";
                  next SHIFTER;
              }
          }
      }

      ## Now fill in the line of data for every tag that was used.
      my $datum_line = '';
      for my $j (@used_tags) {
          my $info = 'NA';
          $info = $datum->{$j} if (defined($datum->{$j}));
          $datum_line .= "${info}\t";
      }
      $datum_line =~ s/\t$/\n/;
      ## print "WRITING: $datum_line\n";
      print $all_out $datum_line;

      my $position_data = $datum->{position};
      my ($chr, $report_pos, $orig, $alt) = $position_data =~ m/^chr_(.*)_pos_(.*)_ref_(.*)_alt_(.*)$/;
      ## Guessing that comma separated means more than one allele?
      if ($alt =~ /\,/) {
          my @tmp = split(/\,/, $alt);
          $alt = $tmp[0];
      }
      my $orig_length = length($orig);
      my $replace_length = length($alt);
      my $delta_length = $replace_length - $orig_length;
      my $relative_pos = $report_pos + $total_delta;
      ## First pull the entire contig sequence.
      ## Keep in mind that we are changing the value of the genome in 'input_genome'
      ## because we aren't dereferencing it.
      my $starting_seq = $input_genome->{$chr}->{sequence};
      my $original_seq = $original_genome{$chr}{sequence};
      ## Give it a useful name, extract the length.
      my $initial_search_seq = $starting_seq;
      my $chromosome_length = length($initial_search_seq);
      ## Search a little bit of context, create a new copy of the contig to
      ## replace the nucleotide of interest.
      my $found_nt = 'X';
      if ($report_pos >= length($original_seq)) {
          warn("The report position: ${report_pos} is too large.\n");
      } else {
          $found_nt = substr($original_seq, $report_pos - 2, 5);
      }
      my $replace_seq = $starting_seq;
      ## Note that freebayes and friends also give back indels which are confusing.
      ## So, for the moment at least, only replace transitions/transversions.
      $points++;
      my $relative = $relative_pos - 1;
      my $replaced = 'X';
      if ($relative >= length($replace_seq)) {
          warn("The relative position: $relative is too large.\n");
      } else {
          $replaced = substr($replace_seq, $relative_pos - 1, $orig_length, $alt);
      }
      my $final_test_seq = $replace_seq;
      my $test_nt = 'X';
      $relative = $relative_pos - 2;
      if ($relative >= length($replace_seq)) {
          warn("The relative position: $relative is too large.\n");
      } else {
          $test_nt = substr($final_test_seq, $relative_pos - 2, 5);
      }
      $total_delta = $total_delta + $delta_length;
      ## Swap out the reference with the alt, $replace_seq gets the new data.
      ## Make a new variable so we can see the change, and
      ## replace the contig in the reference database.
      print "Original region at pos: ${report_pos} are: ${found_nt}, changing ${orig} to ${alt}.  Now they are ${test_nt}. for chr: ${chr} with $total_delta\n" if ($options->{verbose});
      $input_genome->{$chr}->{sequence} = $replace_seq;
      if (!defined($annotations->{$chr})) {
          print $log "${chr} was not defined in the annotations database.\n";
      } else {
          my %tmp = %{$annotations->{$chr}};
        FIND_GENE: foreach my $gene_id (keys %tmp) {
            my $gene = $tmp{$gene_id};
            next FIND_GENE unless (ref($gene) eq 'HASH');
            $vars_by_gene->{$gene_id}->{length} = $gene->{end} - $gene->{start};
            next FIND_GENE if ($relative_pos < $gene->{start} + $total_delta);
            next FIND_GENE if ($relative_pos > $gene->{end} + $total_delta);
            ## I am reasonably certain this variant is in a CDS.
            if (!defined($vars_by_gene->{$gene_id}->{count})) {
                $vars_by_gene->{$gene_id}->{count} = 1;
            } else {
                $vars_by_gene->{$gene_id}->{count}++;
            }
            ## Print the nucleotide changed by this position along with the amino acid substitution

            my $new_chr_string = $input_genome->{$chr}->{sequence};
            ## my $original_chr = Bio::Seq->new(-display_id => $chr, -seq => $original_genome{$chr}{sequence});
            my $original_chr = Bio::Seq->new(-display_id => $chr, -seq => $original_seq);
            ## my $new_chr = Bio::Seq->new(-display_id => $chr, -seq => $input_genome->{$chr}->{sequence});
            my $new_chr = Bio::Seq->new(-display_id => $chr, -seq => $new_chr_string);

            my $cds_relative_position;
            if ($gene->{strand} > 0) {
                $cds_relative_position = $report_pos - $gene->{start};
            } else {
                $cds_relative_position = $gene->{end} - $report_pos;
            }
            my $relative_aminos = floor($cds_relative_position / 3);

            my $original_cds = $original_chr->trunc($gene->{start}, $gene->{end});
            my $new_cds = $new_chr->trunc($gene->{start} + $total_delta, $gene->{end} + $total_delta);
            if ($gene->{strand} < 0) {
                $original_cds = $original_cds->revcom;
                $new_cds = $new_cds->revcom;
            }
            my $original_aa = $original_cds->translate;
            my $new_aa = $new_cds->translate;
            my $original_aa_string = $original_aa->seq;
            my $new_aa_string = $new_aa->seq;
            my $original_aa_position = substr($original_aa_string, $relative_aminos, 1);
            my $new_aa_position = substr($new_aa_string, $relative_aminos, 1);
            my $aa_delta_string = qq"${original_aa_position}${relative_aminos}${new_aa_position}";
            my $report_string = qq"$gene->{id}\t$chr\t${report_pos}\t${orig}_${alt}\t${aa_delta_string}\n";
            print $output_by_gene $report_string;
            last FIND_GENE;
        } ## End iterating over every gene.
      } ## End making sure the chromosome/contig has information
      ## We have looked over all of our annotations at this point, if this variant was not
      ## in any of them, assume it is intercds.
      ## Then this variant is in an inter-cds region, so lets record it as such.
      my $penetrance_string = qq"$chr\t${report_pos}\t${orig}_${alt}\t${penetrance}\n";
      print $output_penetrance $penetrance_string;

  } ## End looking at the each observation.
    $all_out->close();  ## Close out the matrix of observations.
    print $log "Iterated over ${points} attempted variant modifications.\n";
    print $log "Writing variants as (variants/kilobase gene)/megabase chromosome.\n";
    $output_by_gene->close();
    my $var_by_genelength = FileHandle->new(">$options->{output_pkm}");
    print $var_by_genelength qq"gene\tvars_by_length\n";
    foreach my $geneid (keys %{$vars_by_gene}) {
        my $gene_ratio = 0;
        my $good = 1;
        $good = 0 if (! $vars_by_gene->{$geneid}->{count});
        $good = 0 if (! $vars_by_gene->{$geneid}->{length});
        $good = 0 if (! $vars_by_gene->{$geneid}->{length});
        if ($good) {
            $gene_ratio = ($vars_by_gene->{$geneid}->{count} / $vars_by_gene->{$geneid}->{length}) * 1000.0;
        }
        print $var_by_genelength "${geneid}\t${gene_ratio}\n";
    }
    $var_by_genelength->close();
    print $log "Writing out modified genome to $options->{output_genome}.\n";
    my $output_genome = FileHandle->new(">$options->{output_genome}");
    foreach my $ch (sort keys %{$input_genome}) {
        ## my $formatted = $text->format($input_genome->{$ch}->{sequence});
        print $output_genome ">${ch}\n";
        ## Take from: https://www.biostars.org/p/70448/
        foreach my $seq_line (unpack('(a[80])*', $input_genome->{$ch}->{sequence})) {
            print $output_genome "$seq_line\n";
        }
    }
    $output_genome->close();
    $output_penetrance->close();
    print $log "Compressing matrix of all metrics.\n";
    qx"xz -9e -f $options->{output}";
    print $log "Compressing output by gene.\n";
    qx"xz -9e -f $options->{output_by_gene}";
    qx"xz -9e -f $options->{output_penetrance}";
    print $log "Compressing output pkm file.\n";
    qx"xz -9e -f $options->{output_pkm}";
    $log->close();
    $filtered_coverage->close();
    return($count);
}

sub SNP_Ratio_Intron {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        chosen_tag => 'PAIRED',
        coverage_tag => 'DP',
        filter_type => '',
        gff_tag => 'ID',
        gff_type => 'gene',
        gff_cds_parent_type => 'mRNA',
        gff_cds_type => 'CDS',
        jprefix => '80',
        jname => 'parsenp',
        min_value => 0.1,
        qual => 10,
        vcf_cutoff => 5,
        vcf_cutoff => 5,
        vcf_minpct => 0.1,);
    my $print_input = $options->{input};
    my $output_dir = dirname($print_input);
    my $print_output = qq"${output_dir}";
    my $samplename = basename(cwd());

    my $output_suffix = '';
    $output_suffix .= qq"_q-$options->{qual}" if ($options->{qual});
    $output_suffix .= qq"_c-$options->{vcf_cutoff}" if ($options->{vcf_cutoff});
    $output_suffix .= qq"_m-$options->{min_value}" if ($options->{min_value});
    $output_suffix .= qq"_ctag-$options->{coverage_tag}" if ($options->{coverage_tag});
    $output_suffix .= qq"_mtag-$options->{chosen_tag}" if ($options->{chosen_tag});
    my $genome = qq"$options->{libpath}/$options->{libtype}/$options->{species}${output_suffix}.fasta";
    my $stdout = qq"${print_output}/stdout";
    my $stderr = qq"${print_output}/stderr";
    my $output_all = qq"${print_output}/all_tags${output_suffix}.txt.xz";
    my $output_genome = qq"${print_output}/$options->{species}-${samplename}${output_suffix}.fasta";
    my $output_by_gene = qq"${print_output}/variants_by_gene${output_suffix}.txt.xz";
    my $output_penetrance = qq"${print_output}/variants_penetrance${output_suffix}.txt.xz";
    my $output_pkm = qq"${print_output}/pkm${output_suffix}.txt.xz";
    my $output_types = qq"${print_output}/count_types${output_suffix}.txt";
    my $comment_string = qq!
## Parse the SNP data and generate a modified $options->{species} genome.
##  This should read the file:
## ${print_input}
##  and provide 4 new files:
## ${print_output}/count.txt
## ${print_output}/all_tags.txt
## and a modified genome: ${output_genome}
!;
    my $jstring = qq"
use Bio::Adventure::SNP;
my \$result = \$h->Bio::Adventure::SNP::SNP_Ratio_Intron_Worker(
  input => '$print_input',
  species => '$options->{species}',
  chosen_tag => '$options->{chosen_tag}',
  vcf_cutoff => '$options->{vcf_cutoff}',
  min_value => '$options->{min_value}',
  filter_type => '$options->{filter_type}',
  gff_tag => '$options->{gff_tag}',
  gff_type => '$options->{gff_type}',
  gff_cds_parent_type => '$options->{gff_cds_parent_type}',
  gff_cds_type => '$options->{gff_cds_type}',
  qual => '$options->{qual}',
  output_dir => '$output_dir',
  output => '${output_all}',
  output_genome => '${output_genome}',
  output_by_gene => '${output_by_gene}',
  output_penetrance => '${output_penetrance}',
  output_pkm => '${output_pkm}',
  output_types => '${output_types}',
);
";
    my $parse_job = $class->Submit(
        comment => $comment_string,
        chosen_tag => $options->{chosen_tag},
        coverage_tag => $options->{coverage_tag},
        filter_type => $options->{filter_type},
        gff_tag => $options->{gff_tag},
        gff_type => $options->{gff_type},
        gff_cds_parent_type => $options->{gff_cds_parent_type},
        gff_cds_type => $options->{gff_cds_type},
        jdepends => $options->{jdepends},
        jmem => 64,
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => '10:00:00',
        language => 'perl',
        min_value => $options->{min_value},
        output_dir => $output_dir,
        qual => $options->{qual},
        output => $output_all,
        output_genome => $output_genome,
        output_by_gene => $output_by_gene,
        output_penetrance => $output_penetrance,
        output_pkm => $output_pkm,
        output_types => $output_types,
        stderr => $stderr,
        stdout => $stdout,
        vcf_cutoff => $options->{vcf_cutoff},);
    return($parse_job);
}

sub SNP_Ratio_Intron_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        chosen_tag => 'PAIRED',
        coverage_tag => 'DP',
        filter_type => '',
        gff_tag => 'ID',
        gff_type => 'gene',
        gff_exon_type => 'exon',
        gff_cds_parent_type => 'mRNA',
        gff_cds_type => 'CDS',
        max_value => undef,
        min_value => 0.1,
        output => 'all.txt',
        output_genome => 'new_genome.fasta',
        output_by_gene => 'counts_by_gene.txt',
        output_penetrance => 'counts_penetrance.txt',
        output_pkm => 'by_gene_length.txt',
        output_dir => 'outputs/40freebayes',
        output_types => 'outputs/types.txt',
        qual => 10,
        vcf_cutoff => 5,
        vcf_method => 'freebayes',);
    my $species = $options->{species};
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $output_dir = $paths->{output_dir};
    my $log_file = qq"${output_dir}/snp_ratio_intron.stdout";
    my $log = FileHandle->new(">${log_file}");
    my $subst_matrix = $class->Bio::Adventure::Align::Get_Substitution_Matrix(matrix => $options->{matrix});

    my $genome = $paths->{fasta};
    my $gff = $paths->{gff};
    print $log "Reading gff: ${gff}, extracting type: $options->{gff_type} features tagged $options->{gff_tag}.\n";
    my $in_bcf = FileHandle->new("bcftools view $options->{input} |");
    print $log "The large matrix of data will be written to: $options->{output}\n";
    my $filtered_coverage_log = qq"${output_dir}/filtered_coverage_introns.tsv";
    my $all_out = FileHandle->new(">$options->{output}");
    my $type_counter = FileHandle->new(">$options->{output_types}");
    my $filtered_coverage = FileHandle->new(">${filtered_coverage_log}");
    my %type_counts = ();

    ## Read the genome so that I may write a copy with all the nucleotides changed.
    my $input_genome = $class->Bio::Adventure::Read_Genome_Fasta(
        %args, fasta => $genome,);
    my %original_genome;
    for my $chr_key (keys %{$input_genome}) {
        my %internal = %{$input_genome->{$chr_key}};
        $original_genome{$chr_key} = \%internal;
    }
    ## Ok, so I want to simplify the vcf output so that I can create a pseudo
    ## count table of every potential SNP position I want to create a 2 column
    ## file where the first column is a unique ID containing 'chr_pos_ref_new'
    ## and the second column is some sort of score which I will be able to use
    ## for a PCA when I merge a bunch of these together.  Candidates include:
    ## quality score, DP*AF1 (depth * max likelihood estimate) Maybe also consider
    ## DP4 column which is comma separated: #forward-ref, #reverse-ref, #forward-alt, #reverse-alt
    ## If that, then sum all 4 and only take those greater than x (20?), then take
    ## forward-alt+reverse-alt/(sum of all) to get simple SNP ratio
    ## I am going to change my writer of the all.txt file so that it prints out the
    ## values of every observed tag.  For now I am just going to hard-code the order,
    ## but it should not be difficult to parse this out of the header lines of the
    ## bcf file.
    my @mpileup_tag_order = (
        'NS', 'DP', 'DPB', 'AC', 'AN', 'AF', 'RO', 'AO',
        'PRO', 'PAO', 'QR', 'QA', 'PQR', 'PQA', 'SRF', 'SRR', 'SAF',
        'SAR', 'SRP', 'SAP', 'AB', 'ABP', 'RUN', 'RPP', 'RPPR', 'RPL',
        'RPR', 'EPP', 'EPPR', 'DPRA', 'ODDS', 'GTI', 'TYPE', 'CIGAR',
        'NUMALT', 'LEN', 'MQM', 'MQMR', 'PAIRED', 'PAIREDR', 'MIN_DP',
        'END', 'GT', 'GQ', 'GL', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA',
        'MIN_DP', 'snp_pct');
    my @freebayes_tag_order = (
        'NS', 'DP', 'DPB', 'AC', 'AN', 'AF', 'RO', 'AO',
        'PRO', 'PAO', 'QR', 'QA', 'PQR', 'PQA', 'SRF', 'SRR', 'SAF',
        'SAR', 'SRP', 'SAP', 'AB', 'ABP', 'RUN', 'RPP', 'RPPR', 'RPL',
        'RPR', 'EPP', 'EPPR', 'DPRA', 'ODDS', 'GTI', 'TYPE', 'CIGAR',
        'NUMALT', 'LEN', 'MQM', 'MQMR', 'PAIRED', 'PAIREDR', 'MIN_DP',
        'END', 'GT', 'GQ', 'GL', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA',
        'MIN_DP');
    ## Make a data structure containing arrays for the set of observed tags in the vcf file.
    ## It should be an array which ends at length of the number of bcf entries.
    my @tag_order = ();
    if ($options->{vcf_method} eq 'freebayes') {
        @tag_order = @freebayes_tag_order;
    } else {
        @tag_order = @mpileup_tag_order;
    }
    ## Here is the global structure
    my @tag_observations = ();
    ## Use this little hash to decide what tags to write to the final tsv file.
    my %num_observed = ();
    for my $k (@tag_order) {
        $num_observed{$k} = 0;
    }
    ## Then fill it with a hash of the tags from @tag_order
    ## The following hash will be over written.

    print $log "Reading bcf file.\n";
    my $count = 0; ## Use this to count the positions changed in the genome.
    my $num_variants = 0; ## Use this to count the variant positions of reasonably high confidence.
    my $print_count = 0;
  READER: while (my $line = <$in_bcf>) {
      next READER if ($line =~ /^#/);
      ## When using samtools mpileup with the -t LIST, we get some extra
      ## fields at the end which are called tags and tag_info here.
      my ($chr, $pos, $id, $ref, $alt, $qual, $filt, $info, $tags, $tag_info) = split(/\s+/, $line);
      ## Only accept the simplest non-indel mutations.
      next READER unless ($alt eq 'A' or $alt eq 'a' or $alt eq 'T' or $alt eq 't' or
                          $alt eq 'G' or $alt eq 'g' or $alt eq 'C' or $alt eq 'c');
      if ($options->{qual} ne '' && $qual) {
          next READER if ($qual < $options->{qual});
      }
      ## I do not know if indels are here anymore, we will see.
      my @info_list = split(/;/, $info);
      my $snp_id = qq"chr_${chr}_pos_${pos}_ref_${ref}_alt_${alt}";
      my $snp_pct = 0.0;
      my $all_sum = 0;
      my $diff_sum = 0;
      my %individual_tags = (position => $snp_id);
    TAGS: foreach my $element (@info_list) {
        my ($element_type, $element_value) = split(/=/, $element);
        ## Add a little logic to filter for a specific variant type
        ## For the moment only accept keeping a specific type
        if ($element_type eq 'FILTER') {
            if (defined($type_counts{$element_value})) {
                $type_counts{$element_value}++;
            } else {
                $type_counts{$element_value} = 1;
            }
            next READER if ($options->{filter_type} ne '' && $element_value ne $options->{filter_type});
        }
        ## At this point, add some post-processing for tags which are multi-element.
        if ($options->{vcf_method} eq 'mpileup' && $element_type eq 'DP4') {
            my ($same_forward, $same_reverse, $alt_forward, $alt_reverse) = split(/,/, $element_value);
            $diff_sum = $alt_forward + $alt_reverse;
            $all_sum = $same_forward + $same_reverse + $diff_sum;
            $snp_pct = ($alt_forward + $alt_reverse) / $all_sum;
            $snp_pct = nearest(0.01, $snp_pct);
            $individual_tags{snp_pct} = $snp_pct;
        }
        $individual_tags{$element_type} = $element_value;
        if (!defined($num_observed{$element_type})) {
            $num_observed{$element_type} = 1;
            push(@tag_order, $element_type);
        } else {
            $num_observed{$element_type}++;
        }
    } ## End iterating over the tags in the data.
      push(@tag_observations, \%individual_tags);
  } ## End reading the bcf file.
    $in_bcf->close();
    ## Now we should have a big array full of little hashes
    my @used_tags = ('position');
    ## First collect the set of tags of interest.
    my $header_line = "position\t";
    for my $k (@tag_order) {
        if ($num_observed{$k} > 0) {
            push(@used_tags, $k);
            $header_line .= "$k\t";
        }
    }
    $header_line =~ s/\t$/\n/;
    ## Write out the observations, starting with a header line
    ## comprised of the position information followed by the tags.
    print $all_out qq"${header_line}\n";

    ## Read the data base of gff annotations, the genome, and begin
    ## hunting for variants by gene.
    print "Reading the gff into a DB::SeqFeature Store.\n";
    print $log "Reading the gff into a DB::SeqFeature Store.\n";
    my $db = Bio::DB::SeqFeature::Store->new(
        -index_subfeatures => 1,
        -adaptor => 'memory',
        -fasta => $genome,
        -gff => $gff);
    print $log "Extracting sequence data.\n";
    print "Finished the gff reader.\n";

    my %gathered_data = ();
    ## The set of @all_mrna provides the start and end points of each exon.
    ## This will need to be cross referenced against the AUG and stop provided by
    ## the first element of the CDS array in order to complete the set of
    ## the set of translated exons.  Oddly, neither of these two data structures alone
    ## are able to provide this information; I would have thought CDS sufficient.
    my @all_mrna = $db->get_features_by_type($options->{gff_cds_parent_type});
    for my $mrna (@all_mrna) {
        my %gathered_datum = ();
        my $contig = $mrna->seq_id;
        my $start = $mrna->start;
        my $end = $mrna->end;
        my $key_string = qq"${contig} ${start} ${end}";
        my $name = $mrna->load_id;
        my $strand = $mrna->strand;
        $gathered_datum{name} = $name;
        $gathered_datum{strand} = $strand;
        ## These coordinates will get us the AUG and stop coordinates.
        my @cds_features = $mrna->segments(-type => $options->{gff_cds_type});
        my $cds_feature = $cds_features[0];
        my $cds_start = $cds_feature->start;
        my $cds_end = $cds_feature->end;

        ## Now pull the exons, cross reference against the start/stop, and use this information
        my @exon_features = $mrna->segments(-type => $options->{gff_exon_type});
        my %exon_starts = ();
        my %exon_ends = ();
        for my $exon (@exon_features) {
            my $exon_start = $exon->start;
            my $exon_end = $exon->end;
            $exon_starts{$exon} = $exon_start;
            $exon_ends{$exon} = $exon_end;
        }

        my @exons_startends = ();
        my @exons_endstarts = ();
        my $found_start = 0;
        my $found_end = 0;
        my $final_sequence = '';
        my $cds_count = 0;
        my $debug = 0;
      CDS_LOOP: for my $ex (sort{$exon_starts{$a} <=> $exon_starts{$b}} keys %exon_starts) {
          my $startend = [0, 0];
          my $endstart = [0, 0];
          $cds_count++;
          my $this_start = $exon_starts{$ex};
          my $this_end = $exon_ends{$ex};
          if ($found_start == 0) {
              ## If we have not yet found the start of the CDS, and
              ## this exon ends before the start, skip it.
              if ($cds_start > $this_end) {
                  next CDS_LOOP;
              } elsif ($cds_start >= $this_start && $cds_start <= $this_end) {
                  ## Then this exon contains the start.
                  $found_start++;
                  # print "Setting the exon start to the cds start.\n";
                  $this_start = $cds_start;
              }
          }

          if ($found_end == 0) {
              if ($cds_end >= $this_start && $cds_end <= $this_end) {
                  $found_end++;
                  # print "Setting the exon end to the cds end.\n";
                  $this_end = $cds_end;
              }
          } else {
              # print "This exon is too late.\n";
              next CDS_LOOP;
          }
          # print "Extracting $this_start $this_end\n";
          my $segment_sequence = $db->fetch_sequence($contig, $this_start, $this_end);
          $startend = [$this_start, $this_end];
          $endstart = [$this_end, $this_start];
          $final_sequence .= $segment_sequence;
          push(@exons_startends, $startend);
          unshift(@exons_endstarts, $endstart);
      }

        my $cds_obj = Bio::Seq->new(-seq => $final_sequence , -alphabet => 'dna');
        if ($strand < 0) {
            $cds_obj = $cds_obj->revcom;
        }
        my $final_cds_seq = $cds_obj->seq;
        my $final_aa_seq = $cds_obj->translate->seq;
        $gathered_datum{startends} = \@exons_startends;
        $gathered_datum{endstarts} = \@exons_endstarts;
        $gathered_datum{cds} = $final_cds_seq;
        $gathered_datum{aa} = $final_aa_seq;
        $gathered_data{$key_string} = \%gathered_datum;
    } ## We should now have a data structure with the original sequences for every mRNA
    ## and their translations

    ## I saved the set of tags as @tag_observations, use that now in order to
    ## write out the modified genome/CDS/translation.
    my $output_by_gene = FileHandle->new(">$options->{output_by_gene}");
    print $output_by_gene qq"gene\tchromosome\tposition\tfrom_to\taa_subst\tsynonymousp\tblosum_delta\n";
    my $output_penetrance = FileHandle->new(">$options->{output_penetrance}");
    print $output_penetrance qq"chromosome\tposition\tfrom_to\n";
    my $vars_by_gene = {};
    my $all_count = 0;

    print $log "Filtering variant observations and writing a new genome.\n";
    my @all_mRNA_positions = keys %gathered_data;
    my $shifters = 0;
  SHIFTER: while (scalar(@tag_observations)) {
      my $datum = shift(@tag_observations);
      $shifters++;
      ## Now decide if we actually want to write this difference into
      ## a new genome and record it
      my $chosen = $options->{chosen_tag};
      my $depth_tag = $options->{coverage_tag};
      my $vcf_cutoff = $options->{vcf_cutoff};
      if (defined($vcf_cutoff)) {
          if ($datum->{$depth_tag} < $vcf_cutoff) {
              print $filtered_coverage "$datum->{position}\t$datum->{$depth_tag}\t${chosen}\n";
              next SHIFTER;
          }
      }
      if (defined($chosen) && defined($datum->{chosen})) {
          if (defined($options->{min_value}) && defined($options->{max_value})) {
              if ($datum->{chosen} > $options->{max_value} ||
                  $datum->{chosen} < $options->{min_value}) {
                  print $filtered_coverage "$datum->{position}\t$datum->{$depth_tag}\t${chosen}\n";
                  next SHIFTER;
              }
          } elsif (defined($options->{min_value})) {
              if ($datum->{chosen} < $options->{min_value}) {
                  print $filtered_coverage "$datum->{position}\t$datum->{$depth_tag}\t${chosen}\n";
                  next SHIFTER;
              }
          } elsif (defined($options->{max_value})) {
              if ($datum->{chosen} > $options->{max_value}) {
                  print $filtered_coverage "$datum->{position}\t$datum->{$depth_tag}\t${chosen}\n";
                  next SHIFTER;
              }
          }
      }

      ## Now fill in the line of data for every tag that was used.
      my $datum_line = '';
      for my $j (@used_tags) {
          my $info = 'NA';
          $info = $datum->{$j} if (defined($datum->{$j}));
          $datum_line .= "${info}\t";
      }
      $datum_line =~ s/\t$/\n/;
      ## Write out the tags for this variant and then use that
      ## information to go looking for CDS in the range specified by
      ## it.
      print $all_out $datum_line;
      my $position_data = $datum->{position};
      my ($chr, $pos, $orig, $alt) = $position_data =~ m/^chr_(.*)_pos_(.*)_ref_(.*)_alt_(.*)$/;
      ## 1 entry in one sample has a weird whitespace in it, do a little cleaning therefore.
      $chr =~ s/^\s+//g;
      $chr =~ s/\s+$//g;
      $pos =~ s/^\s+//g;
      $pos =~ s/\s+//g;

      ## Take a moment to write a copy of the chromosome with the
      ## new nucleotide which passed our filters.
      my $starting_seq = $input_genome->{$chr}->{sequence};
      ## Give it a useful name, extract the length.
      my $initial_search_seq = $starting_seq;
      my $chromosome_length = length($initial_search_seq);
      ## Search a little bit of context, create a new copy of the contig to
      ## replace the nucleotide of interest.
      my $found_nt = substr($initial_search_seq, $pos - 2, 5);
      my $replace_seq = $starting_seq;
      my $replaced = substr($replace_seq, $pos - 1, 1, $alt);
      ## Swap out the reference with the alt, $replace_seq gets the new data.
      ## Make a new variable so we can see the change, and
      ## replace the contig in the reference database.
      my $final_test_seq = $replace_seq;
      my $test_nt = substr($final_test_seq, $pos - 2, 5);
      ## print "Original region at pos: ${pos} are: ${found_nt}, changing ${orig} to ${alt}.  Now they are ${test_nt}. for chr: ${chr}\n";
      $input_genome->{$chr}->{sequence} = $replace_seq;

      ## Keep in mind that we made an array of arrays containing all
      ## the start/stop positions; so go hunting!
      my $in_cds = 0;
    MRNA: for my $mRNA_position (@all_mRNA_positions) {
        my ($mrna_chr, $mrna_start, $mrna_end) = split(/\s+/, $mRNA_position);
        next MRNA unless ($mrna_chr eq $chr);
        next MRNA if ($pos < $mrna_start);
        next MRNA if ($pos > $mrna_end);
        my %mRNA_datum = %{$gathered_data{$mRNA_position}};
        my $scan_lst = [];
        ## Keeping track of starts vs ends on two strands is too confusing for me,
        ## so I just made a separate datum for each.
        if ($mRNA_datum{strand} > 0) {
            $scan_lst = $mRNA_datum{startends};
        } else {
            $scan_lst = $mRNA_datum{endstarts};
        }
        ## The following variables will be used to define the relative
        ## position of any variants with respect to the AUG, e.g. it
        ## they provide the basis for the translation from absolute
        ## coordinates to local.
        my $relative_nt = undef;
        my $current_pos = 0;
        my $intron_dist = 0;
        my $in_exon = 0;
        my $exon_count = 0;
        my $length_of_exons = 0;
        for my $mRNA_exon (@{$scan_lst}) {
            $exon_count++;
            my $exon_length = abs($mRNA_exon->[1] - $mRNA_exon->[0]) + 1;
            ## I wish I knew btw() before, it makes this process much
            ## easier.  Thus it checks that our position is within the
            ## start/end of each exon, and is smart enough to handle
            ## both strands.
            my $found = btw($pos, $mRNA_exon->[0], $mRNA_exon->[1]);
            if ($found) {
                $in_cds++;
                ## If this position is in fact within an exon, then
                ## figure out how far in we are with respect to the
                ## AUG and switch out the nucleotide/aminoacid in
                ## order to get the new chromosome/protein sequences.
                $in_exon = $exon_count;
                my $this_distance = abs($pos - $mRNA_exon->[0]) + 1;
                my $relative_position = $this_distance + $length_of_exons;
                my $name = $mRNA_datum{name};
                my $substitution = Swap_Sequence(
                    name => $name, nt => $mRNA_datum{cds},
                    aa => $mRNA_datum{aa}, relative => $relative_position,
                    genome_from => $orig, genome_to => $alt,
                    feature_strand => $mRNA_datum{strand});
                my $synonymousp = 'S';
                if ($substitution->{aa_from} ne $substitution->{aa_to}) {
                    $synonymousp = 'NS';
                }
                my $subst_score = $subst_matrix->entry($substitution->{aa_from}, $substitution->{aa_to});
                my $aa_delta_string = $substitution->{aa};
                my $nt_delta_string = $substitution->{nt};
                my $report_string = qq"${name}\t${chr}\t${pos}\t${nt_delta_string}\t${aa_delta_string}\t${synonymousp}\t${subst_score}\n";
                print $output_by_gene $report_string;
            } else {
                ## If the current variant is not within this exon, then add its length to
                ## our exon-length counter so that it will add up when/if we eventually
                ## do get inside an exon.
                $length_of_exons = $length_of_exons + $exon_length;
            }
        } ## End iterating over the set of mRNA start/ends for each
        ## individual mRNA
    } ## End looking at all the mRNAs
      ## Then this variant is in an inter-cds region, so lets record it as such.
      my $report_string = qq"${chr}\t${pos}\t${orig}_${alt}\n";
      print $output_penetrance $report_string;
  } ## End looking at each line from bcftools
    $all_out->close();  ## Close out the matrix of observations.
    $output_by_gene->close();
    $output_penetrance->close();

    ## Finally: clean up and write the number of variants with respect to gene length.
    my $var_by_genelength = FileHandle->new(">$options->{output_pkm}");
    print $var_by_genelength qq"gene\tvars_by_length\n";
    foreach my $geneid (keys %{$vars_by_gene}) {
        my $gene_ratio = 0;
        my $good = 1;
        $good = 0 if (! $vars_by_gene->{$geneid}->{count});
        $good = 0 if (! $vars_by_gene->{$geneid}->{length});
        $good = 0 if (! $vars_by_gene->{$geneid}->{length});
        if ($good) {
            $gene_ratio = ($vars_by_gene->{$geneid}->{count} /
                           $vars_by_gene->{$geneid}->{length}) * 1000.0;
        }
        print $var_by_genelength "${geneid}\t${gene_ratio}\n";
    }
    $var_by_genelength->close();

    ## Write down how many times each variant type was observed.
    for my $vartype (sort keys %type_counts) {
        print $type_counter "${vartype}: $type_counts{$vartype}\n";
    }
    $type_counter->close();

    my $output_genome = FileHandle->new(">$options->{output_genome}");
    foreach my $ch (sort keys %{$input_genome}) {
        ## my $formatted = $text->format($input_genome->{$ch}->{sequence});
        print $output_genome ">${ch}\n";
        ## Take from: https://www.biostars.org/p/70448/
        foreach my $seq_line (unpack('(a[80])*', $input_genome->{$ch}->{sequence})) {
            print $output_genome "$seq_line\n";
        }
    }
    $output_genome->close();
    print $log "Compressing matrix of all metrics.\n";
    qx"xz -9e -f $options->{output}";
    print $log "Compressing output by gene.\n";
    qx"xz -9e -f $options->{output_by_gene}";
    qx"xz -9e -f $options->{output_penetrance}";
    print $log "Compressing output pkm file.\n";
    qx"xz -9e -f $options->{output_pkm}";
    $log->close();
    $filtered_coverage->close();
    return($count);
} ## End of the intron aware SNP Ratio worker

sub Swap_Sequence {
    my %args = @_;
    my $starting_seq = $args{nt};
    my $new_seq = $starting_seq;
    my $swap_position = $args{relative} - 1;
    my $from = $args{genome_from};
    if ($args{feature_strand} < 0) {
        $from =~ tr/AGTCagtc/TCAGtcag/;
    }
    my $to = $args{genome_to};
    if ($args{feature_strand} < 0) {
        $to =~ tr/AGTCagtc/TCAGtcag/;
    }
    my $start_seqobj = Bio::Seq->new(-display_id => $args{name}, -seq => $starting_seq,);
    my $start_nt = substr($starting_seq, $swap_position, 1);
    my $replaced_nt = substr($new_seq, $swap_position, 1, $to);
    my $new_seqobj = Bio::Seq->new(-display_id => qq"$args{name} modified", -seq => $new_seq);
    my $test_nt = substr($new_seq, $swap_position, 1);
    my $context_original = substr($starting_seq, $swap_position - 2, 5);
    my $context_new = substr($new_seq, $swap_position - 2, 5);
    my $amino_position = floor($swap_position / 3);
    my $start_aa_seq = $start_seqobj->translate->seq;
    my $new_aa_seq = $new_seqobj->translate->seq;
    my $original_aa = substr($start_aa_seq, $amino_position, 1);
    my $new_aa = substr($new_aa_seq, $amino_position, 1);
    my $aa_delta_string = qq"${original_aa}${amino_position}${new_aa}";
    my $nt_delta_string = qq"${from}${swap_position}${to}";
    my $report = {
        aa => $aa_delta_string,
        aa_from => $original_aa,
        aa_to => $new_aa,
        nt => $nt_delta_string,
        nt_from => $from,
        nt_to => $to,
    };
    return($report);
}

=head2 C<Snippy>

 Snippy provides a fast variant search tool.
 https://github.com/tseemann/snippy

=cut
sub Snippy {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],);
    my $species = $options->{species};
    my $genome = "$options->{libdir}/$options->{libtype}/${species}.fasta";
    my $query = $options->{input};
    my $query_home = dirname(${query});
    my $query_base = basename(${query}, (".fastq"));
    $query = qq"${query_home}/${query_base}";
    my $prefix_name = qq"snippy";
    my $snippy_name = qq"${prefix_name}_$options->{species}";
    my $suffix_name = $prefix_name;
    if ($options->{jname}) {
        $snippy_name .= qq"_$options->{jname}";
        $suffix_name .= qq"_$options->{jname}";
    }

    my $snippy_input = $options->{input};
    my $test_file = "";
    if ($snippy_input =~ /$options->{delimiter}/) {
        my @pair_listing = split(/$options->{delimiter}/, $snippy_input);
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        $snippy_input = qq" --R1 $pair_listing[0] --R2 $pair_listing[1] ";
        $test_file = $pair_listing[0];
    } else {
        $test_file = File::Spec->rel2abs($snippy_input);
        $snippy_input = qq" --R1 ${test_file} ";
    }

    my $snippy_dir = qq"outputs/snippy_$options->{species}";
    my $stdout = qq"${snippy_dir}/snippy.stdout";
    my $stderr = qq"${snippy_dir}/snippy.stderr";
    my $jstring = qq!mkdir -p ${snippy_dir}
echo "Started snippy at \$(date)" >> ${snippy_dir}/snippy_$options->{species}.stdout

snippy --force \\
  --outdir ${snippy_dir} \\
  --ref ${genome} \\
  ${snippy_input}
!;
    my $comment_string = '## Invoke snippy on some reads.';
    my $snippy = $class->Submit(
        comment => $comment_string,
        jdepends => $options->{jdepends},
        jname => $snippy_name,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => '10:00:00',
        jmem => 48,
        output => qq"${snippy_dir}",
        stdout => $stdout,
        stderr => $stderr,);
    return($snippy);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<samtools> L<snippy> L<vcfutils>

=cut

sub Test_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        ## required => ['species', 'input'],
        chosen_tag => 'PAIRED',
        coverage_tag => 'DP',
        gff_tag => 'ID',
        gff_type => 'gene',
        min_value => 0.5,
        max_value => undef,
        output => 'all.txt',
        output_genome => 'new_genome.fasta',
        output_gff => 'new_genome.gff',
        output_by_gene => 'counts_by_gene.txt',
        output_penetrance => 'variants_penetrance.txt',
        output_pkm => 'by_gene_length.txt',
        output_dir => 'outputs/40freebayes',
        penetrance_tag => 'SAP',
        qual => 10,
        vcf_cutoff => 1,
        vcf_method => 'freebayes',
        verbose => 1,);
    $options->{input} = 'phix.vcf';
    my $genome = 'phix.fasta';
    my $species_gff = 'phix.gff';
    $class->{species} = 'phix';
    $options->{species} = 'phix';
    my $reader = qq"<$options->{input}";
    if ($options->{input} =~ /\.bcf$/) {
        $reader = qq"bcftools view $options->{input} |";
    }
    my $out_dir = dirname($options->{output});
    my $log_file = qq"${out_dir}/snp_ratio.stdout";
    my $log = FileHandle->new(">${log_file}");
    my $gff = Bio::Adventure::Get_FH(input => $species_gff, suffix => qq"| grep -v '^#'");
    print $log "Reading gff: ${species_gff}, extracting type: $options->{gff_type} features tagged $options->{gff_tag}.\n";
    print $log "The large matrix of data will be written to: $options->{output}\n";
    my $all_out = FileHandle->new(">$options->{output}");
    my $filtered_coverage_log = qq"${out_dir}/filtered_coverage.tsv";
    my $output_by_gene = FileHandle->new(">$options->{output_by_gene}");
    print $output_by_gene qq"gene\tchromosome\tposition\tfrom_to\taa_subst\n";
    my $output_penetrance = FileHandle->new(">$options->{output_penetrance}");
    print $output_penetrance qq"chromosome\tposition\tfrom_to\tpenetrance\n";
    my $filtered_coverage = FileHandle->new(">${filtered_coverage_log}");

    ## This function returns a 3 element hash:
    ##  contigs: 1 feature / contig
    ##  sf_sequences: The seq objects for each contig (should not be needed)
    ##  sf_by_contig: The seqfeatures which effectively combine the above.
    my $species2sf = $class->Bio::Adventure::Convert::Species2SF();
    my $contig_features = $species2sf->{contigs};
    my $sf_sequences = $species2sf->{sf_sequences};
    my $seqfeatures = $species2sf->{sf_by_contig};

    ## The following funciton dumps the observed variants in a (v|b)cf file into memory.
    my $bcf_data = Extract_Variants(input => $reader, vcf_method => 'freebayes', qual => $options->{qual});
    my %observations = %{$bcf_data->{observations}};
    my @tag_observations = @{$bcf_data->{tag_observations}};
    my @tag_order = @{$bcf_data->{tag_order}};
    my %num_observed = %{$bcf_data->{num_observed}};
    print "Finished read the bcf file, observed: $observations{snp} SNPS\n";
    ## Now we should have a big array full of little hashes

    ## Use the tags from the (v|b)cf file to write the observation tsv header line.
    my @used_tags = ('position', );
    ## First collect the set of tags of interest.
    my $header_line = "position\t";
    for my $k (@tag_order) {
        if ($num_observed{$k} > 0) {
            push(@used_tags, $k);
            $header_line .= "$k\t";
        }
    }
    $header_line =~ s/\t$/\n/;
    ## Write out the observations, starting with a header line comprised of the position information
    ## followed by the tags.
    print $all_out qq"${header_line}\n";


    my $vars_by_gene = {};
    my $all_count = 0;
    print "Filtering variant observations and writing a new genome.\n";
    my $modified = Modify_Features(observations => \@tag_observations,
                                   seqfeatures => $sf_sequences,
                                   contig_features => $contig_features,
                                   penetrance_tag => $options->{penetrance_tag},
                                   chosen_tag => $options->{chosen_tag},
                                   coverage_tag => $options->{coverage_tag},
                                   vcf_cutoff => $options->{vcf_cutoff},
                                   filtered_coverage => $filtered_coverage,
                                   min_value => $options->{min_value},
                                   max_value => $options->{max_value},
                                   used_tags => \@used_tags,
                                   all_out => $all_out,
                                   output_by_gene => $output_by_gene,
                                   output_penetrance => $output_penetrance,
                                   verbose => $options->{verbose},);

    $all_out->close();  ## Close out the matrix of observations.
    $output_by_gene->close();
    my $var_by_genelength = FileHandle->new(">$options->{output_pkm}");
    print $var_by_genelength qq"gene\tvars_by_length\n";
    foreach my $geneid (keys %{$vars_by_gene}) {
        my $gene_ratio = 0;
        my $good = 1;
        $good = 0 if (! $vars_by_gene->{$geneid}->{count});
        $good = 0 if (! $vars_by_gene->{$geneid}->{length});
        $good = 0 if (! $vars_by_gene->{$geneid}->{length});
        if ($good) {
            $gene_ratio = ($vars_by_gene->{$geneid}->{count} / $vars_by_gene->{$geneid}->{length}) * 1000.0;
        }
        print $var_by_genelength "${geneid}\t${gene_ratio}\n";
    }
    $var_by_genelength->close();

    my $output_genome = new Bio::SeqIO(-file => ">$options->{output_genome}",
                                       -format => 'Fasta');
    my $string = qq!##gff-version 3
!;
    my $final_write = new FileHandle(">$options->{output_gff}");
    my $stringio = new IO::String($string);
    for my $ch (sort keys %{$sf_sequences}) {
        ## my $formatted = $text->format($input_genome->{$ch}->{sequence});
        my $seq_obj = $sf_sequences->{$ch};
        my $written = $output_genome->write_seq($seq_obj);
        my $seq_string = $seq_obj->seq;
        my $seq_id = $seq_obj->id;
        my $seq_name = $seq_obj->display_name;
        ## Take from: https://www.biostars.org/p/70448/
        #foreach my $seq_line (unpack('(a[80])*', $seq_string)) {
        #    print $output_genome "${seq_line}\n";
        #}
        my $contig_sf = $contig_features->{$ch};
        ##for my $top_feature ($contig_sf->get_SeqFeatures()) {
        ##}
        my $all_sf = $seqfeatures->{$ch};
        for my $thingie (@{$all_sf}) {
            my @tags = sort $thingie->get_all_tags();
            my $suffix = '';
          TAGS: for my $t (@tags) {
              next TAGS if ($t eq 'phase');
              next TAGS if ($t eq 'score');
              next TAGS if ($t eq 'seq_id');
              next TAGS if ($t eq 'source');
              next TAGS if ($t eq 'type');
              my @values = $thingie->get_tag_values($t);
              print "Hoping to add tag: $t with value: $values[0]\n";
              $suffix .= qq"${t}=$values[0];";
          }
            $suffix =~ s/;$/\n/;
            my $gffout = new Bio::FeatureIO(-fh => $stringio,
                                            -format => 'GFF',
                                            -version => 3);
            ## This sends the gff to my string IO, so I can add the tags to it...
            my $written = $gffout->write_feature($thingie);
            $string =~ s/$/$suffix/;
        }  ## End for every feature
    } ## End for every contig
    print $final_write $string;
    $final_write->close();
    $output_penetrance->close();

    ## Write the modified gff with the sequence features
    print $log "Compressing matrix of all metrics.\n";
    qx"xz -9e -f $options->{output}";
    print $log "Compressing output by gene.\n";
    qx"xz -9e -f $options->{output_by_gene}";
    qx"xz -9e -f $options->{output_penetrance}";
    print $log "Compressing output pkm file.\n";
    qx"xz -9e -f $options->{output_pkm}";
    $log->close();
    $filtered_coverage->close();
    return(@tag_observations);
}


sub Extract_Variants {
    my %args = @_;
    my $in_bcf = FileHandle->new($args{input});
    my @mpileup_tag_order = (
        'DP', 'ADF', 'AD', 'VDB', 'SGB', 'MQ0F', 'RPB', 'MQB', 'BQB', 'MQSB', 'ADR', 'GT',
        'ICB', 'HOB', 'AC', 'AN', 'DP4', 'MQ',);
    my @freebayes_tag_order = (
        'NS', 'DP', 'DPB', 'AC', 'AN', 'AF', 'RO', 'AO',
        'PRO', 'PAO', 'QR', 'QA', 'PQR', 'PQA', 'SRF', 'SRR', 'SAF',
        'SAR', 'SRP', 'SAP', 'AB', 'ABP', 'RUN', 'RPP', 'RPPR', 'RPL',
        'RPR', 'EPP', 'EPPR', 'DPRA', 'ODDS', 'GTI', 'TYPE', 'CIGAR',
        'NUMALT', 'LEN', 'MQM', 'MQMR', 'PAIRED', 'PAIREDR', 'MIN_DP',
        'END', 'GT', 'GQ', 'GL', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA',
        'MIN_DP');
    ## Make a data structure containing arrays for the set of observed tags in the vcf file.
    ## It should be an array which ends at length of the number of bcf entries.
    my @tag_order = ();
    if ($args{vcf_method} eq 'freebayes') {
        @tag_order = @freebayes_tag_order;
    } else {
        @tag_order = @mpileup_tag_order;
    }
    ## Here is the global structure
    my @tag_observations = ();
    ## Use this little hash to decide what tags to write to the final tsv file.
    my %num_observed = ();
    for my $k (@tag_order) {
        $num_observed{$k} = 0;
    }

    ## Here is an arbitrary vcf line as a refresher to make sure I grabbed them all:
    ## NC_045512.2     3037    .       C       T       224     PASS    DP=244;ADF=0,121;ADR=0,121;AD=0,242;VDB=0;SGB=-0.693147;MQSB=0.976205;MQ0F=0;AC=2;AN=2;DP4=0,0,121,121;MQ=42    GT:PL:DP:SP:ADF:ADR:AD  1/1:254,255,0:242:0:0,121:0,121:0,242

    ## FIXME: Redo this to use Bio::DB::HTS::VCF
    my $current_contig = '';
    my $current_position = 0;
    my $current_offset = 0;
    my %observations = (
        snp => 0,
        insertion => 0,
        deletion => 0,
        complex => 0,
        multi => 0,
        combined => 0,
        other => 0,);
    my $count = 0; ## Use this to count the positions changed in the genome.
    my $num_variants = 0; ## Use this to count the variant positions of reasonably high confidence.
  READER: while (my $line = <$in_bcf>) {
      my $difference_type = 'snp';
      next READER if ($line =~ /^#/);
      ## When using samtools mpileup with the -t LIST, we get some extra
      ## fields at the end which are called tags and tag_info here.
      my ($chr, $pos, $id, $ref, $alt, $qual, $filt, $info, $tags, $tag_info) = split(/\s+/, $line);
      ## Only accept the simplest non-indel mutations.
      if ($args{qual} ne '' && defined($qual)) {
          if ($qual < $args{qual}) {
              next READER;
          }
      }
      ## I do not know if indels are here anymore, we will see.
      my @info_list = split(/;/, $info);
      my $snp_id = qq"chr_${chr}_pos_${pos}_ref_${ref}_alt_${alt}";
      my $snp_pct = 0.0;
      my $all_sum = 0;
      my $diff_sum = 0;
      my %individual_tags = (position => $snp_id);
    TAGS: foreach my $element (@info_list) {
        my ($element_type, $element_value) = split(/=/, $element);

        ## At this point, add some post-processing for tags which are multi-element.
        if ($args{vcf_method} eq 'mpileup' && $element_type eq 'DP4') {
            my ($same_forward, $same_reverse, $alt_forward, $alt_reverse) = split(/,/, $element_value);
            $diff_sum = $alt_forward + $alt_reverse;
            $all_sum = $same_forward + $same_reverse + $diff_sum;
            $snp_pct = ($alt_forward + $alt_reverse) / $all_sum;
            $snp_pct = nearest(0.01, $snp_pct);
            $individual_tags{snp_pct} = $snp_pct;
        }
        $individual_tags{$element_type} = $element_value;
        if (!defined($num_observed{$element_type})) {
            $num_observed{$element_type} = 1;
            push(@tag_order, $element_type);
        } else {
            $num_observed{$element_type}++;
        }
    } ## End iterating over the tags in the data.
      push(@tag_observations, \%individual_tags);
      $difference_type = $individual_tags{TYPE};
      if ($individual_tags{TYPE} eq 'snp') {
          $observations{snp}++;
      } elsif ($individual_tags{TYPE} eq 'ins') {
          $observations{insertion}++;
      } elsif ($individual_tags{TYPE} eq 'del') {
          $observations{deletion}++;
      } elsif ($individual_tags{TYPE} eq 'complex') {
          $observations{complex}++;
      } elsif ($individual_tags{TYPE} eq 'mnp') {
          $observations{multi}++;
      } elsif ($individual_tags{TYPE} =~ /((ins|del|complex|mnp|snp),*)+/) {
          $observations{combined}++
      } else {
          $observations{other}++;
      }
  } ## End reading the bcf file.
    my $closed = $in_bcf->close();
    my $ret = {
        tag_observations => \@tag_observations,
        observations => \%observations,
        tag_order => \@tag_order,
        num_observed => \%num_observed,
  };
    return($ret);
}

sub Modify_Features {
    my %args = @_;
    my @observations = @{$args{observations}};
    my $seqfeatures = $args{seqfeatures};
    my $contig_features = $args{contig_features};
    my $penetrance_tag = $args{penetrance_tag};
    my $chosen_tag = $args{chosen_tag};
    my $coverage_tag = $args{coverage_tag};
    my $vcf_cutoff = $args{vcf_cutoff};
    my $filtered_coverage = $args{filtered_coverage};
    my $min_value = $args{min_value};
    my $max_value = $args{max_value};
    my @used_tags = @{$args{used_tags}};
    my $all_out = $args{all_out};
    my $output_by_gene = $args{output_by_gene};
    my $output_penetrance = $args{output_penetrance};

    ## The following SHIFTER loop is responsible for:
    ##  1. Excluding variants with insufficient depth and/or other tags
    ##  2. Writing a line to the tsv file containing all the tags/observation.
    ##  3. Adding new entries to the @{$deltas->{contig}} and @{$positions->{contig}}
    ##     Any time an observed variant is not 1:1, this should note the contig location
    ##     and the amount every following feature should be shifted.
    ##     It therefore starts at position 0 with a shift of 0.
    my $shifters = 0;
    my $points = 0; ## Number of point mutations observed.
    my $deltas = {};
    my $positions = {};
    my $total_delta = 0; ## I think I do not need this anymore?
    my $current_position = 0;
    my $current_delta = 0;
    my $var_features = {};
  SHIFTER: for my $datum (@observations) {
      my $difference_type = $datum->{TYPE};
      $shifters++;
      ## Now decide if we actually want to write this difference into
      ## a new genome and record it
      my $penetrance = $datum->{$penetrance_tag};
      my $chosen = $chosen_tag;
      if (defined($vcf_cutoff)) {
          if ($datum->{$coverage_tag} < $vcf_cutoff) {
              print "Dropping $datum->{position} because depth ($datum->{$coverage_tag}) is less than $vcf_cutoff.\n" if ($args{verbose});
              print $filtered_coverage "$datum->{position}\t$datum->{$coverage_tag}\t$datum->{chosen}\n";
              next SHIFTER;
          }
      }
      if (defined($chosen) && defined($datum->{chosen})) {
          if (defined($min_value) && defined($max_value)) {
              if ($datum->{chosen} > $max_value) {
                  print "$datum->{position} has exceeded max value: $max_value\n" if ($args{verbose});
                  print $filtered_coverage "$datum->{position}\t$datum->{$coverage_tag}\t$datum->{chosen}\n";
                  next SHIFTER;
              } elsif ($datum->{chosen} < $min_value) {
                  print "$datum->{position} is less than min value: $min_value\n" if ($args{verbose});
                  print $filtered_coverage "$datum->{position}\t$datum->{$coverage_tag}\t$datum->{chosen}\n";
                  next SHIFTER;
              }
          } elsif (defined($min_value)) {
              if ($datum->{chosen} < $min_value) {
                  print "$datum->{position} is less than min value: $min_value\n" if ($args{verbose});
                  print $filtered_coverage "$datum->{position}\t$datum->{$coverage_tag}\t$datum->{chosen}\n";
                  next SHIFTER;
              }
          } elsif (defined($max_value)) {
              if ($datum->{chosen} > $max_value) {
                  print "$datum->{position} has exceeded max value: $max_value\n" if ($args{verbose});
                  print $filtered_coverage "$datum->{position}\t$datum->{$coverage_tag}\t$datum->{chosen}\n";
                  next SHIFTER;
              }
          }
      }

      ## Now fill in the line of data for every tag that was used.
      my $datum_line = '';
      for my $j (@used_tags) {
          my $info = 'NA';
          $info = $datum->{$j} if (defined($datum->{$j}));
          $datum_line .= "${info}\t";
      }
      $datum_line =~ s/\t$/\n/;
      ## print "WRITING: $datum_line\n";
      print $all_out $datum_line;

      my $position_data = $datum->{position};
      my ($chr, $report_pos, $orig, $alt) = $position_data =~ m/^chr_(.*)_pos_(.*)_ref_(.*)_alt_(.*)$/;
      ## Guessing that comma separated means more than one allele?
      if ($alt =~ /\,/) {
          my @tmp = split(/\,/, $alt);
          $alt = $tmp[0];
      }
      my $orig_variant_length = length($orig);
      my $alt_variant_length = length($alt);
      my $delta_variant_length = $alt_variant_length - $orig_variant_length;
      $report_pos = $report_pos + $total_delta;
      if ($delta_variant_length != 0) {
          print "The delta variant length is: $delta_variant_length\n" if ($args{verbose});
          ($positions, $deltas) = Update_Delta_Position(positions => $positions, deltas => $deltas,
                                                        position => $report_pos,
                                                        delta => $delta_variant_length);
      }

      ## Create a feature object for each observation and use it later for recording the results.
      my $var_feature = Bio::SeqFeature::Generic->new(
          -start => $report_pos,
          -end => $report_pos,
          -id => $chr,
          -strand => 0,
          -display_name => $position_data,
          -tag => { ref => $orig, alt => $alt });
      my @vfs = ();
      if (defined($var_features->{$chr})) {
          @vfs = @{$var_features->{$chr}};
      }
      push(@vfs, $var_feature);
      $var_features->{$chr} = \@vfs;

      ## Now modify the contigs with the variant information.
      my $contig_sf = $contig_features->{$chr};
      my $seqio = $contig_sf->seq;
      my $seq = $seqio->seq;
      my $current_len = length($seq);

      ## The following couple of lines are a little silly, but I want to have
      ## one temporary copy of the sequence so I can check my work and another
      ## upon which to actually apply the variant, thus leaving the original
      ## sequence untouched.

      ## This first copy is named 'tmp' because I will destroy it without consideration.
      my $tmp_seq = $seq;
      my $found_nt = substr($tmp_seq, $report_pos - 2, 5);
      ## This copy is named 'rw' because it is the copy of the sequence which is actually changed.
      my $rw_copy = $seq;
      my $replaced = substr($rw_copy, $report_pos - 1, $orig_variant_length, $alt);
      my $replaced_len = length($replaced);
      ## Copy the modified sequence to tmp and extract the region of interest to compare.
      $tmp_seq = $rw_copy;
      my $test_nt = substr($tmp_seq, $report_pos - 2, 5);
      ## Swap out the reference with the alt, $replace_seq gets the new data.
      ## Make a new variable so we can see the change, and
      ## replace the contig in the reference database.
      print "Original region at pos: ${report_pos} are: ${found_nt}, changing ${orig} to ${alt}.  Now they are ${test_nt}. for chr: ${chr} with $total_delta\n" if ($args{verbose});

      ## At this point, a little reminder to self about the datastructures:
      ## The top-level structure is a hashref of contigs pointing to SeqFeatures, 1/contig.
      ## Each contig's seqfeature has a slot named (annoyingly) 'seq' which may be used
      ## to access the underlying SeqIO object containing the actual sequence.
      ## I will therefore do the following:
      ##  1: replace the seqio's sequence with my rw_copy
      ##  2: Update the total position shift (usually by 0)
      ##  3: Record the penetrance (which can probably be put somewhere else)
      ##  4: Atttach the new sequence to the contig's SeqFeature
      ##  5: Put the modified SeqFeature back on the contig hash.

      ## $replaced_seqio is either 0 or 1, $seqio gets changed.
      my $replaced_seqio = $seqio->seq($rw_copy);
      $total_delta = $total_delta + $delta_variant_length;
      my $penetrance_string = qq"${chr}\t${report_pos}\t${orig}_${alt}\t${penetrance}\n";
      print $output_penetrance $penetrance_string;

      ## I think doing sf->seq will replace the seqio
      ## portion of this seqfeature object with the new sequence.
      my $swapping_in = $contig_sf->attach_seq($seqio);
      $contig_features->{$chr} = $contig_sf;

  }  ## End iterating over the observations, now modify the features


    ## The following ~50 lines could probably be moved to a separate function, but they are making
    ## use of enough difference pieces from above they I am worried I will mess it up.
    ## Here is the idea:
    ##  1.  Iterate over the contigs and associated SeqFeatures.
    ##  2.  Iterate over the set of positions/deltas
    ##      (recall that these are arrays recording where on each contig the sequence grew or shrank).
    ##  3.  For any feature (gene/mRNA/CDS/exon/etc) which has a delta event inside it,
    ##      move the end coordinate by the appropriate amount.
    ##  4.  For any feature which is _after_ a delta event, moved the start and end coordinates.
    my @contig_ids = sort keys %{$seqfeatures};
    print "Moving features when the positions have changed.\n" if ($args{verbose});
    ## Whenever we have a change in the offset != 0, modify the start/end positions
    ## of the relevant features.
    my $changed_starts = 0;
    my $changed_ends = 0;
    ## Reset total delta back to 0 and use it again when moving features.
    $total_delta = 0;
  CHR: for my $chr (@contig_ids) {
      print "Moving gff feature positions: starting $chr\n" if ($args{verbose});
      my $chr_sf = $seqfeatures->{$chr};
      my @chr_features = $chr_sf->get_SeqFeatures();
      my $seqio = $chr_sf->seq;
      my @chr_positions = @{$positions->{$chr}};
      my @chr_deltas = @{$deltas->{$chr}};

    POS: for my $d (0 .. $#chr_positions) {
        my $this_position = $chr_positions[$d];
        my $this_delta = $chr_deltas[$d];
        $total_delta = $total_delta + $this_delta;
        next POS if ($total_delta == 0);
      MODIFYPOS: for my $rw_feat (@chr_features) {
          my $feat_seq = $rw_feat->seq;
          my $gene_id = $rw_feat->display_name;
          my $gene_end = $rw_feat->end;
          my $gene_start = $rw_feat->start;
          if ($this_position > $gene_end) {
              next MODIFYPOS;
          } elsif ($this_position > $gene_start &&
                   $this_position <= $gene_end) {
              my $new_end = $gene_end + $total_delta;
              $rw_feat->end($new_end);
              $changed_ends++;
          } elsif ($this_position < $gene_start &&
                   $this_position < $gene_end) {
              my $new_start = $gene_start + $total_delta;
              $changed_starts++;
              my $new_end = $gene_end + $total_delta;
              $changed_ends++;
              $rw_feat->start($new_start);
              $rw_feat->end($new_end);
          } else {
              carp("I do not think this should ever happen!\n");
          }
          ## Print the nucleotide changed by this position along with the amino acid substitution
      } ## End iterating over features/contig
    } ## End iterating over delta positions.
  } ## End iterating over contigs.
    print "Modified: $changed_starts and $changed_ends\n" if ($args{verbose});

    ## Once the starts/ends are set, use the new set of coordinates to find the changed codons
    ## and amino acids.
    ## In order to do this, use the variant feature created above and check if it is
    ## contained within every feature; if so pull the sequence and report it.

    ## I once again split the features to copies prefixed with rw and ro.
    ## This may not be needed.
    print "Extracting modified nucleotides/amino acids.\n" if ($args{verbose});
  CHR2: for my $chr (@contig_ids) {
      my $sfs = $seqfeatures->{$chr};
      my @seqfeats = $sfs->get_SeqFeatures;
      my $feat_count = -1;
    FEATS: for my $e (0 .. $#seqfeats) {
        my $ro_feat = $seqfeats[$e];
        my $rw_feat = $ro_feat;
        my $feat_annot = $ro_feat->annotation;
        my @gene_ids = $feat_annot->get_Annotations('gene_id');
        next FEATS unless (scalar(@gene_ids) > 0);
        my $name = $gene_ids[0]->as_text;
        my @chr_vfs = @{$var_features->{$chr}};
      VF: for my $vf (@chr_vfs) {
          my $inter = $rw_feat->contains($vf);
          next VF unless ($inter);
          my @refs = $vf->get_tag_values('ref');
          my @alts = $vf->get_tag_values('alt');

          ## Add the alt/ref to the gene feature.
          my $absolute_pos = $vf->start;
          my $rel_pos = ($absolute_pos - $rw_feat->start) + 1;
          my $ref = $refs[0];
          my $alt = $alts[0];
          my $mut_string = qq"${ref}${rel_pos}${alt}";
          my $added = $rw_feat->add_tag_value('variants', $mut_string);

          ## Now Extract the new amino acid/CDS sequences.
          ## This is why I made it ro/rw earlier, to make sure
          ## I don't accidently mess up the actual sequence objects when playing with substr()
          ## I was considering the possibility of writing out the full sequence of the
          ## CDS nucleotides or translated amino acids before/after this process
          ## Currently, the only thing done is to record the AUG/M relative position changes.
          my $report_pos = $vf->start;
          my $ro_cds_seq = $ro_feat->seq;
          my $rw_cds_seq = $rw_feat->seq;
          my $cds_relative_position;
          if ($ro_feat->strand > 0) {
              $cds_relative_position = $report_pos - $ro_feat->start;
          } else {
              $cds_relative_position = $ro_feat->end - $report_pos;
          }
          my $relative_aminos = floor($cds_relative_position / 3);
          my $ro_cds_seqstring;
          my $rw_cds_seqstring;
          if ($ro_feat->strand < 0) {
              $ro_cds_seqstring = $ro_cds_seq->revcom->seq;
              $rw_cds_seqstring = $rw_cds_seq->revcom->seq;
          } else {
              $ro_cds_seqstring = $ro_cds_seq->seq;
              $rw_cds_seqstring = $rw_cds_seq->seq;
          }

          my $ro_aa_seq = $ro_cds_seq->translate;
          my $rw_aa_seq = $rw_cds_seq->translate;
          my $ro_aa_string = $ro_aa_seq->seq;
          my $rw_aa_string = $rw_aa_seq->seq;
          my $ro_aa_position = substr($ro_aa_string, $relative_aminos, 1);
          my $rw_aa_position = substr($rw_aa_string, $relative_aminos, 1);
          my $aa_delta_string = qq"${ro_aa_position}${relative_aminos}${rw_aa_position}";
          my $report_string = qq"${name}\t${chr}\t${report_pos}\t${ref}_${alt}\t${aa_delta_string}\n";
          print $output_by_gene $report_string;
      } ## End iterating over the variant features per contig
    } ## End iterating over the sequence features per contig (e.g. genes/mRNAs/etc)
  } ## End iterating over the contigs.
    return($seqfeatures);
}

sub Update_Delta_Position {
    my %args = @_;
    my $positions = $args{positions};
    my $deltas = $args{deltas};
    my $new_position = $args{position};
    my $new_delta = $args{delta};
    my $chr = $args{chr};
    my $delta_variant_length = $args{delta_variant_length};

    my (@delts, @posits);
    if (!defined($deltas->{$chr})) {
        @delts = (0);
        @posits = (0);
    } else {
        @delts = @{$deltas->{$chr}};
        @posits = @{$positions->{$chr}};
    }
    my $current_delta = sum(@delts) + $delta_variant_length;
    print "Pushing $current_delta to the position arrays.\n";
    push(@delts, $current_delta);
    push(@posits, $new_position);
    $deltas->{$chr} = \@delts;
    $positions->{$chr} = \@posits;
    return($positions, $deltas);
}

1;
