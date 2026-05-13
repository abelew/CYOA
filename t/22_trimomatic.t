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

my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
my $local_file = basename($input_file);
if (-r $local_file) {
    print "The local copy: ${local_file} is readable.\n";
} else {
    print "Copying ${input_file} to the current working directory.\n";
    cp($input_file, $local_file);
}

my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
my $trimmer = $cyoa->Bio::Adventure::Trim::Trimomatic(
    input => $input_file,
    jprefix => '22_1',);
ok($trimmer, 'Run Trimomatic');

my $csv_file = $trimmer->{stats}->{main};
ok(-f $csv_file, 'The stats file got created.');
ok(my $actual = $cyoa->Last_Stat(input => $csv_file),
   'Collect Trimomatic Statistics');
my $expected = 'test,10000,9497,503';
unless(ok($expected eq $actual, 'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("-- expected --\n${old}\n-- actual --\n${new}\n");
}

my $input_r1 = qq"${start_dir}/r1.fastq.xz";
my $input_r2 = qq"${start_dir}/r2.fastq.xz";
my $local_r1 = basename($input_r1);
my $local_r2 = basename($input_r2);
if (!-r $local_r1) {
    print "Copying ${input_r1} to the current working directory.\n";
    cp($input_r1, $local_r1);
}
if (!-r $local_r2) {
    print "Copying ${input_r2} to the current working directory.\n";
    cp($input_r2, $local_r2);
}

my $input = qq"r1.fastq.xz:r2.fastq.xz";
$trimmer = $cyoa->Bio::Adventure::Trim::Trimomatic(
    input => $input,
    jprefix => '22_2',);
ok($trimmer, 'Run Trimomatic Pairwise');
$csv_file = $trimmer->{stats}->{main};
ok(-f $csv_file, 'The stats file got created.');
ok($actual = $cyoa->Last_Stat(input => $csv_file),
   'Collect Trimomatic Statistics');
$expected = 'r1,108212,94708,10366,1049,2089';
unless(ok($expected eq $actual, 'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("-- expected --\n${old}\n-- actual --\n${new}\n");
}

chdir($start);
