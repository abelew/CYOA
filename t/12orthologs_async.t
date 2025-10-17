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
my $source_groupa_strep = qq"${source_dir}/genome/genbank/mgas_5005.gb.xz";
my $source_groupb_strep = qq"${source_dir}/genome/genbank/sagalactiae_cjb111.gb.xz";

my $start = getcwd();
my $new = 'test_output_async';
mkdir($new);
chdir($new);

ok(-r $source_groupa_strep, "Found source groupa strep genome at: ${source_groupa_strep}.");
ok(-r $source_groupb_strep, "Found source groupb strep genome at: ${source_groupb_strep}.");

my $cyoa = Bio::Adventure->new(
    basedir => cwd(),
    libdir => 'reference',);

## First extract the amino acid sequences of our two streptococci genbank files.
my $groupa_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => $source_groupa_strep, output_dir => '.',);
ok($groupa_convert, 'Converted the group A strep genbank to fasta/gff/etc.');
my $status = $cyoa->Wait(job => $groupa_convert);
ok($status->{State} eq 'COMPLETED', 'The group A strep conversion completed.');

ok(-r 'mgas_5005.fsa', 'Found output groupA fasta file.');
my $groupb_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => $source_groupb_strep, output_dir => '.',);
ok($groupb_convert, 'Converted the group B strep genbank to fasta/gff/etc.');
$status = $cyoa->Wait(job => $groupb_convert);
ok($status->{State} eq 'COMPLETED', 'The group B strep conversion completed.');
ok(-r 'sagalactiae_cjb111.fsa', 'Found output groupB fasta file.');

## Now invoke orthofinder
my $orthos = $cyoa->Bio::Adventure::Align::OrthoFinder(
    input => qq"$groupa_convert->{output_pep_fasta}:$groupb_convert->{output_pep_fasta}",
    jprefix => '12',);
## Note, we cannot just wait for the ortholog job, we need to wait
## for the following post-processor.
$status = $cyoa->Wait(job => $orthos->{namer});
ok($status->{State} eq 'COMPLETED', 'The orthofinder job finished.');

ok(-r $orthos->{named_out}, 'Created tsv of named orthologs.');
my $test_cmd = qq"head -n 3 $orthos->{output} | grep -v Orthogroup";
my $actual;
print "Querying tsv of named orthologs via:
$test_cmd
";
my $wtf = qx"${test_cmd}";
chomp $wtf;
$wtf =~ s/\s+$//g;
$wtf =~ s/^$//g;
## I made some changes to how I print the orthofinder outputs...
my $expected = qq!OG0000000	AAZ50854.1.231, AAZ51149.1.523, AAZ51186.1.560, AAZ51426.1.799, AAZ51473.1.840, AAZ51585.1.951, AAZ51601.1.967, AAZ51695.1.1061, AAZ51855.1.1217, AAZ51980.1.1342, AAZ52136.1.1498, AAZ52139.1.1501, AAZ52143.1.1505, AAZ52247.1.1608	QOW75830.1.1011, QOW75889.1.1073, QOW75941.1.1127, QOW75969.1.1155, QOW75973.1.1159, QOW75989.1.1178, QOW76092.1.1294, QOW76121.1.1325, QOW76148.1.1352, QOW76318.1.1536, QOW76332.1.1550, QOW76362.1.1580, QOW76375.1.1595, QOW76483.1.1706, QOW76495.1.1719, QOW76500.1.1724, QOW76676.1.1907, QOW76764.1.1997, QOW77065.1.309, QOW77115.1.360, QOW77116.1.361, QOW77166.1.414, QOW77167.1.415, QOW77187.1.440, QOW77207.1.460, QOW77578.1.860, QOW77582.1.864, QOW77612.1.894, QOW77645.1.930, QOW77719.1.32\t
OG0000001	AAZ50901.1.278, AAZ51053.1.428, AAZ51298.1.672, AAZ51422.1.795, AAZ51566.1.933, AAZ52343.1.1700	QOW75940.1.1126, QOW75947.1.1133, QOW76090.1.1292, QOW76236.1.1452, QOW76494.1.1718, QOW76683.1.1914, QOW76747.1.1980, QOW76794.1.24, QOW76903.1.138, QOW77073.1.317, QOW77338.1.596, QOW77601.1.883, QOW77610.1.892
!;
chomp $expected;
$expected =~ s/\s+$//g;
$expected =~ s/^$//g;
ok($wtf =~ /OG/, "Check the named orthologs via: $test_cmd and got $wtf");

ok(-r $orthos->{single_out}, 'Created tsv of single-named orthologs.');
print "Checking the single output: $orthos->{single_out}\n";
$test_cmd = qq"sort $orthos->{single_out} | head -n 2";
$actual = qx"${test_cmd}";
print "Invoking
${test_cmd}
to check the single output.\n";
$expected = qq!OG0000240
OG0000241
!;
unless(ok($expected eq $actual, 'Is the resulting set of named orthologs as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

ok(-r $orthos->{single_name_out}, 'Created tsv of single-named orthologs.');
print "Checking the single output: $orthos->{single_out}\n";
$test_cmd = qq"sort $orthos->{single_name_out} | head -n 2 | awk '{print \$2}'";
$actual = qx"${test_cmd}";
print "Invoking:
${test_cmd}
to check the single named output.\n";
$expected = qq!AAZ50620.1
AAZ50621.1.1
!;
unless(ok($expected eq $actual, 'Is the resulting set of named orthologs as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

chdir($start);
