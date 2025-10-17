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



my $cyoa = Bio::Adventure->new(basedir => cwd());
my $jelly = $cyoa->Bio::Adventure::Count::Jellyfish(input => $local_phix_fasta, jprefix => 40,);
ok($jelly, 'Submit jellyfish jobs to the cluster.');
my $status = $cyoa->Wait(job => $jelly);
ok($status->{State} eq 'COMPLETED', 'The jellyfish jobs completed.');

my $first_job_outputs = $jelly->{compression}->{output};
my @output_files = split(/:/, $first_job_outputs);
for my $o (@output_files) {
    ok (-f $o, "Output file create: ${o}");
}

my $expected = qq"1 5374
";
my $mer_table = $jelly->{'13'}->{histogram_file};
my $actual = qx"less ${mer_table}";
print $actual;
unless(ok($expected eq $actual, 'Do we get an expected 9mer histogram?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}
