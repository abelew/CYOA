# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
my $start_dir = dist_dir('Bio-Adventure');

my $groupa_strep = qq"${start_dir}/genome/mgas_5005.gb.xz";
my $groupb_strep = qq"${start_dir}/genome/sagalactiae_cjb111.gb.xz";
my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    libdir => cwd());

## First extract the amino acid sequences of our two streptococci genbank files.
my $groupa_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => $groupa_strep, output_dir => '.',);
ok($groupa_convert, 'Converted the group A strep genbank to fasta/gff/etc.');
my $groupb_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => $groupb_strep, output_dir => '.',);
ok($groupb_convert, 'Converted the group B strep genbank to fasta/gff/etc.');

## Now invoke orthofinder
my $orthos = $cyoa->Bio::Adventure::Align::OrthoFinder(
    input => qq"$groupa_convert->{output_pep_fasta}:$groupb_convert->{output_pep_fasta}",
    jprefix => '12',);

ok(-r $orthos->{named_out}, 'Created tsv of named orthologs.');
my $test_cmd = qq"head -n 2 $orthos->{named_out}";
my $actual = qx"${test_cmd}";
print "Invoking ${test_cmd} to check the named output.\n";
my $expected = qq!Orthogroup	mgas_5005	mgas_5005_name	sagalactiae_cjb111	sagalactiae_cjb111_name
"OG0000000"	"M5005_Spy0235, M5005_Spy0531, M5005_Spy0568, M5005_Spy0808, M5005_Spy0855, M5005_Spy0967, M5005_Spy0983, M5005_Spy1077, M5005_Spy1237, M5005_Spy1362, M5005_Spy1518, M5005_Spy1521, M5005_Spy1525, M5005_Spy1629"	"M5005_Spy0235"	"M5005_Spy0235|AAZ50854.1|amino acid transport ATP-binding protein|GI:71852831"	"ID870_00170, ID870_01560, ID870_01815, ID870_01820, ID870_02085, ID870_02090, ID870_02215, ID870_02315, ID870_04340, ID870_04360, ID870_04510, ID870_04690, ID870_05100, ID870_05410, ID870_05680, ID870_05820, ID870_05840, ID870_05935, ID870_06515, ID870_06675, ID870_06810, ID870_07810, ID870_07885, ID870_08035, ID870_08190, ID870_08770, ID870_08835, ID870_08860, ID870_10060, ID870_10510"	"ID870_00170"	"ID870_00170|QOW77719.1|ABC transporter ATP-binding protein"
!;
unless(ok($expected eq $actual, 'Is the resulting set of named orthologs as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

ok(-r $orthos->{single_out}, 'Created tsv of single-named orthologs.');
print "Checking the single output: $orthos->{single_out}\n";
$test_cmd = qq"sort $orthos->{single_out} | head -n 2";
$actual = qx"${test_cmd}";
print "Invoking ${test_cmd} to check the single output.\n";
$expected = qq!OG0000240
OG0000241
!;
unless(ok($expected eq $actual, 'Is the resulting set of named orthologs as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

ok(-r $orthos->{single_name_out}, 'Created tsv of single-named orthologs.');
print "Checking the single output: $orthos->{single_out}\n";
$test_cmd = qq"sort $orthos->{single_name_out} | head -n 2";
$actual = qx"${test_cmd}";
print "Invoking ${test_cmd} to check the single named output.\n";
$expected = qq!OG0000240	M5005_Spy0001	M5005_Spy0001|dnaA|AAZ50620.1|chromosomal replication initiator protein|GI:71852597	ID870_09700	ID870_09700|dnaA|QOW77779.1|chromosomal replication initiator protein DnaA
OG0000241	M5005_Spy0002	M5005_Spy0002|dnaN|AAZ50621.1|DNA polymerase III beta chain|GI:71852598	ID870_09695	ID870_09695|dnaN|QOW76609.1|DNA polymerase III subunit beta
!;
unless(ok($expected eq $actual, 'Is the resulting set of named orthologs as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}
chdir($start);
