# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename dirname";
use File::Copy qw"cp mv";
use File::Path qw"make_path remove_tree rmtree";
use File::ShareDir qw"dist_dir dist_file module_dir";
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
    stranded => 'no',);
my $paths = $cyoa->Bio::Adventure::Config::Get_Paths(subroutine => 'Hisat2');
make_path($paths->{index_dir}); ## Make a directory for the phix indexes.
ok (-d $paths->{index_dir}, qq"The index directory: $paths->{index_dir} exists\n");
my $fasta_dir = dirname($paths->{fasta});
make_path($fasta_dir);
my $gff_dir = dirname($paths->{gff});
make_path($gff_dir);
ok(-d $fasta_dir, qq"The fasta directory exists: ${fasta_dir}.");
ok(-d $gff_dir, qq"The fasta directory exists: ${gff_dir}.");
if (!-r 'test_forward.fastq.gz') {
    ok(cp($input_file, 'test_forward.fastq.gz'), 'Copying data.');
}
if (!-r $paths->{fasta}) {
    ok(cp($phix_fasta, $paths->{fasta}), qq"Copying phix fasta file to $paths->{fasta}.");
}
if (!-r $paths->{gff}) {
    ok(cp($phix_gff, $paths->{gff}), qq"Copying phix gff file to $paths->{gff}.");
}

my $hisat = $cyoa->Bio::Adventure::Map::Hisat2(
    input => 'test_forward.fastq.gz',
    jprefix => '24',);
ok($hisat, 'Run Hisat2');
my $sam_file = $hisat->{samtools}->{output};
my $htseq_file = $hisat->{htseq}->[0]->{output};
my $stats_file = $hisat->{stats}->{output};

ok(-f $sam_file, qq"The sorted bamfile was created to ${sam_file}.");
ok(-f $htseq_file, qq"The count table was created: ${htseq_file}.");
ok(-f $stats_file, qq"The hisat stats were recorded to ${stats_file}.");

my $actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect Hisat Statistics');

my $expected = qq"hisat_phix_genome_test_output.stderr,10000,49,9951,0";
unless(ok($expected eq $actual, 'Are the hisat stats as expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

$expected = qq"phiX174p01\t5
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t13
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t31
__too_low_aQual\t0
__not_aligned\t9951
__alignment_not_unique\t0
";

$actual = qx"less ${htseq_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}
chdir($start);
