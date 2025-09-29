# -*-Perl-*-
use strict;
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename dirname";
use File::Copy qw"cp mv";
use File::Path qw"make_path remove_tree rmtree";
use File::ShareDir qw"dist_dir dist_file module_dir";
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
make_path('reference/CDS/fasta'); ## Make a directory for the phix indexes.
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
    species => 'phix',);
my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
    species => 'phix', input => $local_phix_fasta,
    gff => $local_phix_gff, libdir => 'reference',
    gff_type => 'gene', gff_tag => 'gene_id');
ok($gff2fasta, 'Submitted gff2fasta job.');
##my $status = $cyoa->Wait(job => $gff2fasta);
##ok($status->{State} eq 'COMPLETED', 'The gff2fasta completed.');
## gff2fasta does not queue itself but runs purely synchronously at this time.
ok(-f $gff2fasta->{output_nt}, "Gff2Fasta produced a nucleotide CDS file: $gff2fasta->{output_nt}.");
my $phix_transcripts = 'reference/CDS/fasta/phix.fasta';
ok(cp($gff2fasta->{output_nt}, $phix_transcripts), "Copied $gff2fasta->{output_nt} to $phix_transcripts.");
ok(-f $phix_transcripts, "Found phix transcript file: ${phix_transcripts}.");
print "Copied gff2fasta output to $phix_transcripts.\n";

my $index = $cyoa->Bio::Adventure::Index::Salmon_Index(
    input => $phix_transcripts,
    jprefix => '25',
    decoy => 0,);
ok($index, 'Submitted index job.');
my $status = $cyoa->Wait(job => $index);
ok($status->{State} eq 'COMPLETED', 'The indexer completed.');

my $salmon = $cyoa->Bio::Adventure::Map::Salmon(
    input => $local_input_fastq,
    jprefix => '25',);
ok($salmon, 'Run salmon');
$status = $cyoa->Wait(job => $salmon);
ok($status->{State} eq 'COMPLETED', 'The salmon quantifier completed.');
my $salmon_file = $salmon->{output};
ok(-f $salmon_file, 'The salmon quant.sf file was created.');

my $stats_file = $salmon->{stats}->{output};
ok(-f $stats_file, "The salmon stats file was created: ${stats_file}");
my $actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect salmon Statistics');

my $expected = qq"lib_format_counts.json,phix,44,44,44,0,0.4318181818181818";
unless(ok($expected eq $actual, 'Are the mapping stats from salmon expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

$expected = qq"Name	Length	EffectiveLength	TPM	NumReads
phiX174p01	1406	1165.522	79003.686307	7.540
chr_NC_001422_id_phiX174p01_start_1_end_136	136	4.000	0.000000	0.000
phiX174p02	890	639.000	0.000000	0.000
phiX174p03	312	62.000	287513.563515	1.460
chr_NC_001422_id_phiX174p03_start_1_end_51	51	2.000	0.000000	0.000
phiX174p04	171	6.000	0.000000	0.000
phiX174p05	261	24.000	508825.648290	1.000
phiX174p06	3618	3330.751	124657.101888	34.000
phiX174p07	634	359.172	0.000000	0.000
phiX174p08	276	32.000	0.000000	0.000
phiX174p09	117	3.000	0.000000	0.000
phiX174p10	528	277.000	0.000000	0.000
phiX174p11	987	743.819	0.000000	0.000
";

$actual = qx"less ${salmon_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}
chdir($start);
