# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };

my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());

my $trimmer = $cyoa->Bio::Adventure::Trim::Cutadapt(
    input => $input_file,
    jprefix => '21',);
ok($trimmer, 'Run Cutadapt');

my $stats_file = $trimmer->{stats}->{output};
my $actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect cutadapt Statistics');
my $expected = 'cutst,10000,2240,11,7760,2229';
unless(ok($expected eq $actual,
          'Are the fastqc results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("$old\n$new\n");
}
## Return to the start.
chdir($start);
