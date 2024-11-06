# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $phix_fasta = qq"${start_dir}/genome/phix.fasta";
my $phix_gff = qq"${start_dir}/genome/phix.gff";

my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    gff_tag => 'gene_id',
    gff_type => 'gene',
    libdir => cwd(),
    species => 'phix',
    stranded => 'no');
make_path('genome/indexes'); ## Make a directory for the phix indexes.
make_path('genome/fasta');
make_path('genome/gff');
my $paths = $cyoa->Bio::Adventure::Config::Get_Paths();

if (!-r 'test_forward.fastq.gz') {
    ok(cp($input_file, 'test_forward.fastq.gz'), 'Copying data.');
}

if (!-r 'genome/fasta/phix.fasta') {
    ok(cp($phix_fasta, $paths->{fasta}), qq"Copying phix fasta file: $paths->{fasta}.");
}

if (!-r 'genome/gff/phix.gff') {
    ok(cp($phix_gff, $paths->{gff}), qq"Copying phix gff file: $paths->{gff}.");
}

my $bt1 = $cyoa->Bio::Adventure::Map::Bowtie(
    input => qq"test_forward.fastq.gz",
    gff_tag => 'gene_id',
    gff_type => 'gene',
    jprefix => '23',
    species => 'phix',);
ok($bt1, 'Run Bowtie1');
my $sam_file = $bt1->{samtools}->{output};
my $htseq_file = $bt1->{htseq}->[0]->{output};
my $stats_file = $bt1->{stats}->{output};
ok(-f $sam_file, qq"The sorted bamfile was created: ${sam_file}");
ok(-f $htseq_file, qq"The count table was created: ${htseq_file}");
ok(-f $stats_file, qq"The stats file was created: ${stats_file}");

my $actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect Bowtie1 Statistics');
my $expected = qq"test_output-v0M1.stderr,v0M1,0,10000,0,9970,0,outputs/23bowtie_phix/test_output-v0M1_rpos_sno_gene_gene_id.count.xz";
unless(ok($expected eq $actual, 'Are the bowtie stats as expected?')) {
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
