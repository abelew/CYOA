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

my $source_dir = dist_dir('Bio-Adventure');
my $source_phix_fasta = qq"${source_dir}/genome/fasta/phix.fasta";
my $source_phix_gff = qq"${source_dir}/genome/gff/phix.gff";

my $start = getcwd();
my $new = 'test_output_async';
mkdir($new);
chdir($new);

ok(-r $source_phix_fasta, "Found source phix fasta at: ${source_phix_fasta}.");
ok(-r $source_phix_gff, "Found source phix gff at: ${source_phix_gff}.");
my $local_phix_fasta = 'reference/genome/fasta/phix.fasta';
my $local_phix_gff = 'reference/genome/gff/phix.gff';

my $output_phix_cds = 'phix_gene_gene_id_nt.fasta';
my $output_phix_aa = 'phix_gene_gene_id_aa.fasta';
my $local_phix_cds = 'reference/CDS/fasta/phix.fasta';
my $local_phix_aa = 'reference/amino_acid/fasta/phix.fasta';

make_path('reference/genome/fasta'); ## Make a directory for the phix indexes.
make_path('reference/genome/gff'); ## Make a directory for the phix indexes.
make_path('reference/CDS/fasta');
make_path('reference/amino_acid/fasta');

if (-r $local_phix_fasta) {
    ok(!-z $local_phix_fasta, 'Local phix file is not null.');
} else {
    ok(cp($source_phix_fasta, $local_phix_fasta), "Copying phix fasta file from ${source_phix_fasta} to ${local_phix_fasta}.");
}

if (-r $local_phix_gff) {
    ok(!-z $local_phix_gff, 'Local phix file is not null.');
} else {
    ok(cp($source_phix_gff, $local_phix_gff), "Copying phix gff file from ${source_phix_gff} to ${local_phix_gff}.");
}

my $cyoa = Bio::Adventure->new(basedir => cwd());
if (-r $output_phix_cds) {
    ok(-r $output_phix_aa, "${output_phix_cds} and ${output_phix_aa} both exist.");
} else {
    my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
        input => $local_phix_fasta, gff => $local_phix_gff,
        gff_type => 'gene', gff_tag => 'gene_id',
        libpath => 'reference');
    ## In 11_blast this is:
    ## libdir => 'share');
    ok($gff2fasta, 'Run gff2fasta.');
    my $status = $cyoa->Wait(job => $gff2fasta);
}

ok(-r $output_phix_cds, "Found nucleotide fasta file at ${output_phix_cds}.");
ok(mv($output_phix_cds, $local_phix_cds), "Moved phix cds file to ${local_phix_cds}.");
ok(-r $output_phix_aa, "Found nucleotide amino acid file at ${output_phix_aa}.");
ok(mv($output_phix_aa, $local_phix_aa), "Moved phix aa file to ${local_phix_aa}.");

print "About to run Split_Align_Fasta.\n";
my $run_fasta = $cyoa->Bio::Adventure::Align_Fasta::Split_Align_Fasta(
    align_jobs => 1, fasta_tool => 'fasta36', input => $local_phix_cds,
    jprefix => '10', library => $local_phix_fasta, parse => 1,);
ok($run_fasta, 'Run Split_Align_Fasta.');
## There is a problem in how I am testing this:
## The Worker job is allowed to finish before the parsing jobs start
## as a result it does not collect the job ID for the parser.
## It should not be difficult for me to make Wait() smart enough to find the child job IDs
## but for the moment I will simply put a sleep() here.
my $status = $cyoa->Wait(job => $run_fasta, sleep => 20);


## my $parsed_file = 'outputs/fasta_phix_cds_nt_phix/phix_cds_nt_vs_phix.parsed.txt';
my $parsed_file = $run_fasta->{output};
my $expected = qq"chr_NC_001422_id_phiX174p01_start_1_end_136
chr_NC_001422_id_phiX174p02_start_1_end_136
chr_NC_001422_id_phiX174p03_start_1_end_51
phiX174p01
phiX174p02
phiX174p03
phiX174p04
phiX174p05
phiX174p06
phiX174p07
";
print "Testing ${parsed_file} for expected outputs.\n";
my $test_cmd = qq"less ${parsed_file} | awk '{print \$1}' | sort | head --lines 10";
print "Invoking: ${test_cmd}\n";
my $actual = qx"${test_cmd}";
unless(ok($expected eq $actual, 'Is the resulting table of hits expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
