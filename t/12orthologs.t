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
    input => qq"$groupa_convert->{pep_fasta}:$groupb_convert->{pep_fasta}",
    jprefix => '12',);

ok(-r $orthos->{named_out}, 'Created tsv of named orthologs.');
my $test_cmd = qq"head -n 2 $orthos->{named_out}";
my $actual = qx"${test_cmd}";
print "Invoking ${test_cmd} to check the named output.\n";
my $expected = qq!Orthogroup	mgas_5005	mgas_5005_name	sagalactiae_cjb111	sagalactiae_cjb111_name
"OG0000000"	"AAZ50854.1 GI_71852831 ; amino acid transport ATP-binding protein, AAZ51149.1 ftsE ; GI_71853126 ; cell division ATP-binding protein, AAZ51186.1 sagG ; GI_71853163 ; streptolysin S export ATP-binding protein, AAZ51426.1 srtF ; GI_71853403 ; lantibiotic transport ATP-binding protein, AAZ51473.1 proV ; GI_71853450 ; glycine betaine transport ATP-binding protein, AAZ51585.1 GI_71853562 ; ABC transporter ATP-binding protein, AAZ51601.1 GI_71853578 ; histidine transport ATP-binding protein, AAZ51695.1 glnQ.2 ; GI_71853672 ; glutamine transport ATP-binding protein, AAZ51855.1 artP ; GI_71853832 ; arginine transport ATP-binding protein, AAZ51980.1 GI_71853957 ; transporter, AAZ52136.1 GI_71854113 ; transporter, AAZ52139.1 GI_71854116 ; cobalt transport ATP-binding protein cbiO, AAZ52143.1 GI_71854120 ; ABC transporter ATP-binding protein, AAZ52247.1 salX ; GI_71854224 ; lantibiotic transport ATP-binding protein"	"AAZ50854.1 GI_71852831 "	""	"QOW75830.1 ATP-binding cassette domain-containing protein, QOW75889.1 ATP-binding cassette domain-containing protein, QOW75941.1 amino acid ABC transporter ATP-binding protein, QOW75969.1 ATP-binding cassette domain-containing protein, QOW75973.1 ABC transporter ATP-binding protein, QOW75989.1 ABC transporter ATP-binding protein, QOW76092.1 ABC transporter ATP-binding protein, QOW76121.1 ftsE ; cell division ATP-binding protein FtsE, QOW76148.1 amino acid ABC transporter ATP-binding protein, QOW76318.1 ABC transporter ATP-binding protein, QOW76332.1 amino acid ABC transporter ATP-binding protein, QOW76362.1 ABC transporter ATP-binding protein, QOW76375.1 betaine/proline/choline family ABC transporter ATP-binding protein, QOW76483.1 amino acid ABC transporter ATP-binding protein, QOW76495.1 ABC transporter ATP-binding protein, QOW76500.1 sugar ABC transporter ATP-binding protein, QOW76676.1 ABC transporter ATP-binding protein, QOW76764.1 ABC transporter ATP-binding protein, QOW77065.1 ABC transporter ATP-binding protein, QOW77115.1 ABC transporter ATP-binding protein, QOW77116.1 ABC transporter ATP-binding protein, QOW77166.1 ABC transporter ATP-binding protein, QOW77167.1 ABC transporter ATP-binding protein, QOW77187.1 ABC transporter ATP-binding protein, QOW77207.1 amino acid ABC transporter ATP-binding protein, QOW77578.1 ATP-binding cassette domain-containing protein, QOW77582.1 ABC transporter ATP-binding protein, QOW77612.1 ABC transporter ATP-binding protein, QOW77645.1 amino acid ABC transporter ATP-binding protein, QOW77719.1 ABC transporter ATP-binding protein"	"QOW75830.1 ATP-binding cassette domain-containing protein"	""
!;
unless(ok($expected eq $actual, 'Is the resulting set of named orthologs as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

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
AAZ50621.1
!;
unless(ok($expected eq $actual, 'Is the resulting set of named orthologs as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}
chdir($start);
