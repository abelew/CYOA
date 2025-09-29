# -*-Perl-*-
use strict;
use Bio::Adventure;
use Cwd;
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
    stranded => 'no',
    gff_type => 'gene',
    gff_tag => 'gene_id',
    species => 'phix',);

my $bwa = $cyoa->Bio::Adventure::Map::BWA(input => $local_input_fastq, jprefix => '26',);
ok($bwa, 'Submitted bwa jobs.');
my $status = $cyoa->Wait(job => $bwa);
ok($status->{State} eq 'COMPLETED', 'The bwa mapping completed.');

## Some files of interest:
my $htseq_mem = $bwa->{htseq_mem}->{output};
my $htseq_aln = $bwa->{htseq_aln}->{output};
my $reporter_sam = $bwa->{reporter}->{output};
my $bwa_output_files = $bwa->{output};
my $aln_bam = $bwa->{samtools_aln}->{output};
ok(-f $htseq_mem, "htseq output from the mem alignment was created: ${htseq_mem}");
ok(-f $aln_bam, "samtools converted the sam to a bam file: ${aln_bam}");

my $expected = qq"phiX174p01\t5
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
my $actual = qx"less ${htseq_mem}";
unless(ok($expected eq $actual, "Check bwa count tables mem version.")) {
    my ($old, $new) = diff($expected, $actual);
    diag("expected:\n$old\nactual:\n$new\n");
}

$expected = qq"phiX174p01\t3
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
__ambiguous\t28
__too_low_aQual\t0
__not_aligned\t9956
__alignment_not_unique\t0
";

$actual = qx"less ${htseq_aln}";
unless(ok($expected eq $actual, "Check bwa count tables aln version.")) {
    my ($old, $new) = diff($expected, $actual);
    diag("expected:\n$old\nactual:\n$new\n");
}
