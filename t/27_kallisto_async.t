# -*-Perl-*-
use strict;
use Bio::Adventure;
use Cwd;
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
    species => 'phix', input => $local_phix_fasta, gff => $local_phix_gff,
    libdir => 'reference', gff_type => 'gene', gff_tag => 'gene_id',);
ok($gff2fasta, 'Run gff2fasta.');
ok(-f $gff2fasta->{output_nt}, "Gff2Fasta produced a nucleotide CDS file: $gff2fasta->{output_nt}.");
my $phix_transcripts = 'reference/CDS/fasta/phix.fasta';
make_path('reference/CDS/fasta');
ok(cp($gff2fasta->{output_nt}, $phix_transcripts), "Copied $gff2fasta->{output_nt} to $phix_transcripts.");
ok(-f $phix_transcripts, "Found phix transcript file: ${phix_transcripts}.");

my $index = $cyoa->Bio::Adventure::Index::Kallisto_Index(
    input => $phix_transcripts, libtype => 'CDS', jprefix => '27');
my $status = $cyoa->Wait(job => $index);
ok($status->{State} eq 'COMPLETED', 'The kallisto indexer completed.');
my $kallisto = $cyoa->Bio::Adventure::Map::Kallisto(
    input => $local_input_fastq, libtype => 'CDS', jprefix => '27',);
ok($kallisto, 'Submitted kallisto job.');
$status = $cyoa->Wait(job => $kallisto);
ok($status->{State} eq 'COMPLETED', 'The bwa mapping completed.');

##my $kallisto_file = 'outputs/46kallisto_phix/abundance.tsv';
my $kallisto_file = $kallisto->{output};
ok(-f $kallisto_file, qq"The kallisto otuput file was created: ${kallisto_file}.");

my $expected = qq"target_id	length	eff_length	est_counts	tpm
phiX174p01	1406	1125.2	8.50998	189149
chr_NC_001422_id_phiX174p01_start_1_end_136	136	96.9989	0	0
phiX174p02	890	735.894	0.281695	9573.51
chr_NC_001422_id_phiX174p02_start_1_end_136	136	96.9989	0	0
phiX174p03	312	99.6032	1.20833	303402
chr_NC_001422_id_phiX174p03_start_1_end_51	51	14.352	0	0
phiX174p04	171	131.999	0	0
phiX174p05	261	179.863	1	139048
phiX174p06	3618	3648.51	29.5203	202354
phiX174p07	634	637.649	2.84075	111419
phiX174p08	276	236.999	0	0
phiX174p09	117	77.9989	0	0
phiX174p10	528	488.999	0	0
phiX174p11	987	909.775	1.63894	45054.3
";

my $actual = qx"less ${kallisto_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
