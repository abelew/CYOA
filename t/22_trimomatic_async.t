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
my $new = 'test_output_async';
mkdir($new);
chdir($new);
my $local_file = basename($input_file);
if (-r $local_file) {
    print "The local copy: ${local_file} is readable.\n";
} else {
    print "Copying ${input_file} to the current working directory.\n";
    cp($input_file, $local_file);
}

my $cyoa = Bio::Adventure->new(basedir => cwd());
my $trimmer = $cyoa->Bio::Adventure::Trim::Trimomatic(
    input => $local_file,
    jprefix => '22_1',);
ok($trimmer, 'Submitted trimomatic SE job.');
my $status = $cyoa->Wait();
ok($status->{State} eq 'COMPLETED', 'The trimomatic jobs completed.');
my $stats_file = $trimmer->{stats}->{main};
ok(my $actual = $cyoa->Last_Stat(input => $stats_file), 'Collect Trimomatic Statistics');
my $expected = 'test,10000,9497,503';
unless(ok($expected eq $actual, 'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

my $r1 = qq"${start_dir}/r1.fastq.xz";
my $r2 = qq"${start_dir}/r2.fastq.xz";
my $local_r1 = basename($r1);
my $local_r2 = basename($r2);
if (-r $local_r1) {
    print "The local copy: ${local_r1} is readable, assuming r2 is too.\n";
} else {
    print "Copying ${local_r1} and ${local_r2} to the current working directory.\n";
    cp($r1, $local_r1);
    cp($r2, $local_r2);
}

my $input = qq"r1.fastq.xz:r2.fastq.xz";
$trimmer = $cyoa->Bio::Adventure::Trim::Trimomatic(
    input => $input,
    jprefix => '22_2',);
ok($trimmer, 'Submitted trimomatic PE job.');
$status = $cyoa->Wait();
ok($status->{State} eq 'COMPLETED', 'Paired trimomatic completed on the cluster.');
$stats_file = $trimmer->{stats}->{main};
ok($actual = $cyoa->Last_Stat(input => $stats_file), 'Collect Trimomatic Statistics');
## Hey you, fix this, I never finished the pairwise stats collector.
$expected = 'r1,108212,94708,10366,1049,2089';
unless(ok($expected eq $actual, 'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("-- expected --\n${old}\n-- actual --\n${new}\n");
}

chdir($start);
