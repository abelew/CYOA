# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');
my $phix_fasta = qq"${start_dir}/genome/phix.fasta";
my $phix_gff = qq"${start_dir}/genome/phix.gff";

use Bio::Adventure;

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

ok(-r $phix_fasta, 'Have input fasta.');
ok(-r $phix_gff, 'Have input gff.');
my $phix_local = 'genome/phix.fasta';
my $gff_local = 'genome/phix.gff';
my $cds_local = 'phix_gene_gene_id_nt.fasta';
my $aa_local = 'phix_gene_gene_id_aa.fasta';
my $cds_genome = 'genome/phix_gene_gene_id_nt.fasta';
my $aa_genome = 'genome/phix_gene_gene_id_aa.fasta';

make_path('genome'); ## Make a directory for the phix indexes.
if (-r $phix_local) {
    ok(!-z $phix_local, 'Local phix file is not null.');
} else {
    ok(cp($phix_fasta, $phix_local), 'Copying phix fasta file.');
}

if (-r $gff_local) {
    ok(!-z $gff_local, 'Local gff is not null.');
} else {
    ok(cp($phix_gff, $gff_local), 'Copying phix gff file.');
}

my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
if (!-r $cds_local) {
    my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
        input => $phix_local, gff => $gff_local,
        gff_type => 'gene', gff_tag => 'gene_id',
        libpath => 'genome');
    ## In 11_blast this is:
    ## libdir => 'share');
    ok($gff2fasta, 'Run gff2fasta.');
}

ok(-r $cds_local, 'Created nucleotide fasta file.');
ok(mv($cds_local, $cds_genome), 'Moved phix cds file.');
ok(-r $aa_local, 'Created nucleotide amino acid file.');
ok(mv($aa_local, $aa_genome), 'Moved phix cds file.');

my $run_fasta = $cyoa->Bio::Adventure::Align_Fasta::Split_Align_Fasta(
    cluster => 0,
    align_jobs => 1,
    fasta_tool => 'fasta36',
    input => $cds_genome,
    library => $phix_local,
    jprefix => '10',
    parse => 1,);
ok($run_fasta, 'Run Split_Align_Fasta.');

## my $parsed_file = 'outputs/fasta_phix_cds_nt_phix/phix_cds_nt_vs_phix.parsed.txt';
my $parsed_file = $run_fasta->{output};
## Caveat: every fasta36 run will give slightly different E-values.
my $expected = qq"117
136
136
1406
171
261
276
312
3618
51
";

print "Testing ${parsed_file} for expected outputs.\n";
my $test_cmd = qq"less ${parsed_file} | awk '{print \$2}' | sort | head --lines 10";
print "Invoking: ${test_cmd}\n";
my $actual = qx"${test_cmd}";
unless(ok($expected eq $actual, 'Is the resulting table of hits expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
