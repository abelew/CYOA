# -*-Perl-*-
use strict;
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename";
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
    gff_tag => 'gene_id',
    gff_type => 'gene',
    libdir => 'reference',
    species => 'phix',
    stranded => 'no');

my $bt1 = $cyoa->Bio::Adventure::Map::Bowtie(
    input => $local_input_fastq,
    gff_tag => 'gene_id',
    gff_type => 'gene',
    jprefix => '23',
    species => 'phix',);
ok($bt1, 'Submitted bowtie1 jobs.');
my $status = $cyoa->Wait(job => $bt1->{htseq}->[0]);
ok($status->{State} eq 'COMPLETED', 'The bowtie jobs completed.');
my $sam_file = $bt1->{samtools}->{output};
my $htseq_file = $bt1->{htseq}->[0]->{output};
my $stats_file = $bt1->{stats}->{output};
ok(-f $sam_file, qq"The sorted bamfile was created: ${sam_file}");
ok(-f $htseq_file, qq"The count table was created: ${htseq_file}");
ok(-f $stats_file, qq"The stats file was created: ${stats_file}");

my $actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect Bowtie1 Statistics');
my $expected = qq"test,v0M1,0,10000,0,9970,0,outputs/23bowtie_phix/test_output-v0M1_rpos_sno_gene_gene_id.csv.xz";
unless(ok($expected eq $actual, "Expected output in $stats_file")) {
    my ($old, $new) = diff($expected, $actual);
    print "--Expected--\n${old}\n--Actual--\n${new}\n";
}

$expected = qq"phiX174p01\t1
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t8
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t21
__too_low_aQual\t0
__not_aligned\t9970
__alignment_not_unique\t0
";

$actual = qx"less ${htseq_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}
chdir($start);
