# -*-Perl-*-
use strict;
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree";
use File::Basename qw"basename";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
use Test::More qw"no_plan";

my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $input_local = basename($input_file);
my $start = getcwd();
my $new = 'test_output_async';
mkdir($new);
chdir($new);
my $cyoa = Bio::Adventure->new(basedir => cwd());

if (!-r $input_local) {
    ok(cp($input_file, $input_local), 'Copying data.');
}

my $trimmer = $cyoa->Bio::Adventure::Trim::Cutadapt(
    input => $input_local,
    jprefix => '21',);
ok($trimmer, 'Submit cutadapt job.');
my $status = $cyoa->Wait(job => $trimmer->{stats});
ok($status->{State} eq 'COMPLETED', 'Cutadapt completed on cluster.');
my $stats_file = $trimmer->{stats}->{output};
my $actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect cutadapt Statistics');
my $expected = 'cutst,10000,2240,11,7760,2229';
unless(ok($expected eq $actual, 'Are the cutadapt results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("$old\n$new\n");
}
## Return to the start.
chdir($start);
