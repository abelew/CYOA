# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename dirname";
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
if (!-r 'genome/phix.fasta') {
    ok(cp($phix_fasta, 'genome/phix.fasta'), 'Copying phix fasta file.');
}
my $cyoa = Bio::Adventure->new(basedir => cwd());
my $jelly = $cyoa->Bio::Adventure::Count::Jellyfish(
    input => 'genome/phix.fasta',
    jprefix => 40,);
ok($jelly, 'Submit jellyfish jobs to the cluster.');
my $status = $cyoa->Wait(job => $jelly);
ok($status->{State} eq 'COMPLETED', 'The jellyfish jobs completed.');

my $first_job_outputs = $jelly->{compression}->{output};
my @output_files = split(/:/, $first_job_outputs);
for my $o (@output_files) {
    ok (-f $o, "Output file create: ${o}");
}

my $expected = qq"1 5079
2 140
3 5
4 1
";
my $actual = qx"less outputs/40jellyfish_test_output_async/test_output_async_9.hist.xz";
unless(ok($expected eq $actual, 'Do we get an expected 9mer histogram?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}
