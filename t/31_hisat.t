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
my $new = 'test_output';
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
    cluster => 0,
    basedir => cwd(),
    gff_tag => 'gene_id',
    gff_type => 'gene',
    libdir => 'reference',
    species => 'phix',
    stranded => 'no',);
my $hisat = $cyoa->Bio::Adventure::Map::Hisat2(
    input => $local_input_fastq,
    prefix => '31_1',);
ok($hisat, 'Run Hisat2 with single-end input.');
my $sam_file = $hisat->{samtools}->{main};
my $count_file = $hisat->{counter}->{main};
my $stats_file = $hisat->{stats}->{main};

ok(-f $sam_file, qq"The sorted bamfile was created to ${sam_file}.");
ok(-f $count_file, qq"The count table was created: ${count_file}.");
ok(-f $stats_file, qq"The hisat stats were recorded to ${stats_file}.");

my $actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect Hisat Statistics');

my $expected = qq"hisat_phix_genome_test.stderr,10000,49,9951,0";
unless(ok($expected eq $actual, 'Are the hisat stats as expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

## Stupid featureCounts puts a random space at the end of the first line
$expected = qq!"1"
Length
1542
1026
363
171
261
3618
634
276
117
528
987
!;

my $cmd = qq!xzcat ${count_file} | awk '{print \$6}'!;
print "Running
${cmd}\n";
$actual = qx"$cmd";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

## Repeat with some paired end input which is xz compressed.
my $input_r1 = qq"${source_dir}/r1.fastq.xz";
my $input_r2 = qq"${source_dir}/r2.fastq.xz";
my $local_r1 = basename($input_r1);
my $local_r2 = basename($input_r2);
if (!-r $local_r1) {
    print "Copying ${input_r1} to the current working directory.\n";
    cp($input_r1, $local_r1);
}
if (!-r $local_r2) {
    print "Copying ${input_r2} to the current working directory.\n";
    cp($input_r2, $local_r2);
}
my $input = qq"r1.fastq.xz:r2.fastq.xz";
my $hs_paired = $cyoa->Bio::Adventure::Map::Hisat2(
    input => $input,
    jprefix => '31_2',);

ok($hs_paired, 'Run Hisat2 with paired-end input.');
my $sam_paired = $hs_paired->{samtools}->{main};
my $count_paired = $hs_paired->{counter}->{main};
my $stats_paired = $hs_paired->{stats}->{main};
ok(-f $sam_paired, qq"The sorted bamfile was created to ${sam_file}.");
ok(-f $count_paired, qq"The count table was created: ${count_file}.");
ok(-f $stats_paired, qq"The hisat stats were recorded to ${stats_file}.");

$actual = $cyoa->Last_Stat(input => $stats_paired);
ok($actual, 'Collect Hisat Statistics for a paired-end run.');

$expected = qq"hisat_phix_genome_test.stderr,10000,49,9951,0";
unless(ok($expected eq $actual, 'Are the hisat stats as expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

## Stupid featureCounts puts a random space at the end of the first line
$expected = qq!"1"
Length
1542
1026
363
171
261
3618
634
276
117
528
987
!;

$cmd = qq!xzcat ${count_file} | awk '{print \$6}'!;
print "Running
${cmd}\n";
$actual = qx"$cmd";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}





chdir($start);
