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
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $phix_fasta = qq"${start_dir}/genome/phix.fasta";
my $phix_gff = qq"${start_dir}/genome/phix.gff";

my $start = getcwd();
my $new = 'test_output_async';
mkdir($new);
chdir($new);

make_path('genome/indexes'); ## Make a directory for the phix indexes.
make_path('genome/gff');
make_path('genome/fasta');
if (!-r 'test_forward.fastq.gz') {
    ok(cp($input_file, 'test_forward.fastq.gz'), 'Copying data.');
}
if (!-r 'genome/fasta/phix.fasta') {
    ok(cp($phix_fasta, 'genome/fasta/phix.fasta'), 'Copying phix fasta file.');
    ## my $uncompressed = qx"gunzip genome/phix.fastq.gz && mv genome/phix.fasta.gz genome/phix.fasta";
}
if (!-r 'genome/gff/phix.gff') {
    ok(cp($phix_gff, 'genome/gff/phix.gff'), 'Copying phix gff file.');
    ## my $uncompressed = qx"gunzip genome/phix.gff.gz && mv genome/phix.gff.gz genome/phix.gff";
}

my $cyoa = Bio::Adventure->new(
    basedir => cwd(),
    libdir => cwd(),
    stranded => 'no',
    gff_type => 'gene',
    gff_tag => 'gene_id',
    species => 'phix',);

## Note that when using a cluster for this, one must remember to use a locally copied file for indexing/input.
my $index = $cyoa->Bio::Adventure::Index::BWA_Index(
    input => 'genome/phix.fasta',
    jprefix => '26',);
my $status = $cyoa->Wait(job => $index);
ok($status->{State} eq 'COMPLETED', 'The bwa mapping completed.');
## Check that the indexes were created:
ok(-f $index->{output_sa}, "The .sa index file was created: $index->{output_sa}");
ok(-f $index->{output_pac}, "The .sa index file was created: $index->{output_pac}");
ok(-f $index->{output_bwt}, "The .sa index file was created: $index->{output_bwt}");
ok(-f $index->{output_ann}, "The .sa index file was created: $index->{output_ann}");
ok(-f $index->{output_amb}, "The .sa index file was created: $index->{output_amb}");
ok(-f $index->{output_fa}, "The .sa index file was created: $index->{output_fa}");

my $bwa = $cyoa->Bio::Adventure::Map::BWA(
    input => 'test_forward.fastq.gz',
    jprefix => '26',);
ok($bwa, 'Submitted bwa jobs.');
$status = $cyoa->Wait(job => $bwa);
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
