# -*-Perl-*-
use strict;
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename dirname";
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
use Test::More qw"no_plan";

my $source_dir = dist_dir('Bio-Adventure');
my $source_phix_fasta = qq"${source_dir}/genome/fasta/phix.fasta";
my $source_phix_gff = qq"${source_dir}/genome/gff/phix.gff";
my $source_input_fastq = qq"${source_dir}/test_forward.fastq.gz";

my $start = getcwd();
my $new = 'test_output_async';
mkdir($new);
chdir($new);

ok(-r $source_phix_fasta, "Found source phix fasta at: ${source_phix_fasta}.");
ok(-r $source_phix_gff, "Found source phix gff at: ${source_phix_gff}.");
my $local_phix_fasta = 'reference/genome/fasta/phix.fasta';
my $local_phix_gff = 'reference/genome/gff/phix.gff';
my $local_input_fastq = 'test_forward.fastq.gz';
make_path('reference/genome/fasta'); ## Make a directory for the phix indexes.
make_path('reference/genome/gff'); ## Make a directory for the phix indexes.
if (-r $local_phix_fasta) {
    ok(!-z $local_phix_fasta, 'Local phix file is not null.');
} else {
    ok(cp($source_phix_fasta, $local_phix_fasta), "Copying phix fasta file from ${source_phix_fasta} to ${local_phix_fasta}.");
}

if (-r $local_phix_gff) {
    ok(!-z $local_phix_gff, 'Local phix file is not null.');
} else {
    ok(cp($source_phix_gff, $local_phix_gff), "Copying phix gff file from ${source_phix_gff} to ${local_phix_gff}.");
}

if (-r $local_input_fastq) {
    ok(!-z $local_phix_gff, 'Local fastq file is not null.');
} else {
    ok(cp($source_input_fastq, $local_input_fastq), "Copying input fastq from: ${source_input_fastq} to ${local_input_fastq}.");
}

my $cyoa = Bio::Adventure->new(
    basedir => cwd(),
    libdir => 'reference',
    species => 'phix',
    gff_tag => 'ID',
    gff_type => 'CDS',
    stranded => 'no',
    vcf_cutoff => 1,);
my $variant = $cyoa->Bio::Adventure::SNP::Align_SNP_Search(
    input => $local_input_fastq, introns => 0, jprefix => '30',);
ok($variant, 'Submit variant search.');
my $status = $cyoa->Wait(job => $variant);
ok($status->{State} eq 'COMPLETED', 'The variant search completed.');
## Some files to check out:


my $gene_variants = $variant->{parse}->{output_by_gene};
ok(-r $gene_variants, "Created the variants by gene: $gene_variants.");
my $expected = qq"gene	chromosome	position	from_to	aa_subst	synonymousp	blosum_delta
";
my $actual = qx"less ${gene_variants}";
unless(ok($expected eq $actual, 'Do we get a set of gene variants (empty).')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

my $parsed = $variant->{parse}->{output};
ok(-r $parsed, "Created the parsed variants: $parsed");
$expected = qq"position	NS	DP	DPB	AC	AN	AF	RO	AO	PRO	PAO	QR	QA	PQR	PQA	SRF	SRR	SAF	SAR	SRP	SAP	AB	ABP	RUN	RPP	RPPR	RPL	RPR	EPP	EPPR	DPRA	ODDS	GTI	TYPE	CIGAR	NUMALT	LEN	MQM	MQMR	PAIRED	PAIREDR	DP	RO	QR	AO	QA	MEANALT

";
$actual = qx"less ${parsed}";
unless(ok($expected eq $actual, 'Do we get an empty variant table?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

my $penetrance = $variant->{parse}->{output_penetrance};
ok(-r $penetrance, "Created the penetrance file: ${penetrance}.");
$expected = qq"chromosome	position	from_to
";
$actual = qx"less ${penetrance}";
unless(ok($expected eq $actual, 'Do we get an empty penetrance file?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

my $genome = $variant->{parse}->{output_genome};
ok(-r $genome, "Created the new genome: ${genome}");
$expected = qq">NC_001422
GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAG
";
$actual = qx"less ${genome} | head -n 2";
unless(ok($expected eq $actual, 'Do we get the original genome back?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}


my $types = $variant->{parse}->{output_types};
ok(-r $types, "Created the types: ${types}");
chdir($start);
