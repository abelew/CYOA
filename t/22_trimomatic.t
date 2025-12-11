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

my $csv_file = $trimmer->{stats}->{outputs}->{main};
ok(-f $csv_file, 'The stats file got created.');
ok(my $actual = $cyoa->Last_Stat(input => $csv_file),
   'Collect Trimomatic Statistics');
my $expected = 'test_forward,10000,9476,524';
unless(ok($expected eq $actual, 'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("-- expected --\n${old}\n-- actual --\n${new}\n");
}

my $input = qq"${start_dir}/r1.fastq.xz:${start_dir}/r2.fastq.xz";
$trimmer = $cyoa->Bio::Adventure::Trim::Trimomatic_Pairwise(
    input => $input,
    jprefix => '22_2',);
ok($trimmer, 'Run Trimomatic Pairwise');
$csv_file = $trimmer->{stats}->{outputs}->{main};
ok(-f $csv_file, 'The stats file got created.');
ok($actual = $cyoa->Last_Stat(input => $csv_file),
   'Collect Trimomatic Statistics');
$expected = 'r1,108212,94183,10743,1074,2212';
unless(ok($expected eq $actual, 'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("-- expected --\n${old}\n-- actual --\n${new}\n");
}

chdir($start);
