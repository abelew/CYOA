# -*-Perl-*-
use strict;
use Bio::Adventure;
use Bio::SeqIO;
use Cwd;
use File::Basename qw"basename";
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
if (!-r $output_phix_cds) {
    my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
        input => $local_phix_fasta, gff => $local_phix_gff,
        gff_tag => 'gene_id', gff_type => 'gene',
        libpath => 'reference',);
    ok($gff2fasta, 'Run gff2fasta.');
}

if (-r $local_phix_cds) {
    print "The local phix CDS input already exists at: ${local_phix_cds}.\n";
} else {
    ok(-r $output_phix_cds, "Found nucleotide fasta file at ${output_phix_cds}.");
    ok(mv($output_phix_cds, $local_phix_cds), "Moved phix cds file to ${local_phix_cds}.");
}
if (-r $local_phix_aa) {
    print "The local phix amino acid input alread exists at: ${local_phix_aa}.\n";
} else {
    ok(-r $output_phix_aa, "Found nucleotide amino acid file at ${output_phix_aa}.");
    ok(mv($output_phix_aa, $local_phix_aa), "Moved phix aa file to ${local_phix_aa}.");
}

my $run_blast = $cyoa->Bio::Adventure::Align_Blast::Split_Align_Blast(
    input => $local_phix_cds,
    jprefix => '11',
    library => $local_phix_fasta,
    number => 1,
    parse => 1,);
my $status = $cyoa->Wait(job => $run_blast, sleep => 20);
ok($status->{State} eq 'COMPLETED', 'The blast job has completed.');
my $parsed_txt = $run_blast->{parser}->{parsed_output};
my $count_txt = $run_blast->{parser}->{count_output};
my $blast_txt = $run_blast->{parser}->{input};
ok(-r $parsed_txt, "The parsed blast output file was created: ${parsed_txt}.");
ok(-r $count_txt, "The count-table blast output was created: ${count_txt}.");
my $actual_cmd = qq"head ${parsed_txt}\n";
print "Collecting outputs via:
${actual_cmd}";
my $actual = qx"${actual_cmd}";
my $expected = qq"QUERYNAME	real_name	Chromosome	Start	End	%ID	Score	Sig	CompLength	Hit_Ident	hits
Query_1	phiX174p01	NC_001422		3981	5386	100	1406	0	1406	100	1
Query_2	chr_NC_001422_id_phiX174p01_start_1_end_136	NC_001422		1	136	100	136	7.66468e-71	136	100	1
Query_3	phiX174p02	NC_001422		4497	5386	100	890	0	890	100	1
Query_4	chr_NC_001422_id_phiX174p02_start_1_end_136	NC_001422		1	136	100	136	7.66468e-71	136	100	1
Query_5	phiX174p03	NC_001422		5075	5386	100	312	2.68412e-168	312	100	1
Query_6	chr_NC_001422_id_phiX174p03_start_1_end_51	NC_001422		1	51	100	51	4.41014e-24	51	100	1
Query_7	phiX174p04	NC_001422		51	221	100	171	3.41394e-90	171	100	1
Query_8	phiX174p05	NC_001422		133	393	100	261	4.99259e-140	261	100	1
Query_9	phiX174p06	NC_001422		358	3975	100	3618	0	3618	100	1
";
unless(ok($expected eq $actual,
          'Did we parse the expected blast output?')) {
    my($old, $new) = diff($expected, $actual);
    diag("$old\n$new\n");
}
