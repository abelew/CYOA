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

my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());

my $trimmer = $cyoa->Bio::Adventure::Trim::Trimomatic_Single(
    input => $input_file,
    jprefix => '22',);
ok($trimmer, 'Run Trimomatic');

my $csv_file = $trimmer->{stats}->{output};
ok(-f $csv_file, 'The stats file got created.');
ok(my $actual = $cyoa->Last_Stat(input => $csv_file),
   'Collect Trimomatic Statistics');
my $expected = 'test_forward,10000,9316,684';
unless(ok($expected eq $actual, 'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}
chdir($start);
