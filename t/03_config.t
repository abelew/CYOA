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

print "About to instantiate cyoa with input_file: $input_file\n";
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

my $trimmer = $cyoa->Bio::Adventure::Trim::Trimomatic_Single(
    input => $input_file,
    jprefix => '22',);
ok($trimmer, 'Run Trimomatic');
