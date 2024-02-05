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
my $start_dir = dist_dir('Bio-Adventure');
my $input = 'test_forward.fastq.gz';
my $input_file = qq"${start_dir}/${input}";
my $phix_fasta_local = 'genome/phix.fasta';
my $phix_fasta = qq"${start_dir}/${phix_fasta_local}";
my $phix_gff_local = 'genome/phix.gff';
my $phix_gff = qq"${start_dir}/${phix_gff_local}";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

make_path('genome/indexes'); ## Make a directory for the phix indexes.
if (!-r $input) {
    ok(cp($input_file, $input), 'Copying data.');
}
if (!-r $phix_fasta_local) {
    ok(cp($phix_fasta, $phix_fasta_local), 'Copying phix fasta file.');
}
if (!-r $phix_gff_local) {
    ok(cp($phix_gff, $phix_gff_local), 'Copying phix gff file.');
}

my $cyoa = Bio::Adventure->new(
    basedir => cwd(),
    libdir => cwd(),
    species => 'phix',
    gff_tag => 'ID',
    gff_type => 'CDS',
    stranded => 'no',
    vcf_cutoff => 1,
    jprefix => '30');
my $index = $cyoa->Bio::Adventure::Index::BT2_Index(input => $phix_fasta_local,);
my $status = $cyoa->Wait(job => $index);
ok($status->{State} eq 'COMPLETED', 'The bowtie2 indexing completed.');
my $variant = $cyoa->Bio::Adventure::SNP::Align_SNP_Search(
    input => $input, introns => 0);
ok($variant, 'Submitted variant search.');
$status = $cyoa->Wait(job => $variant);
ok($status->{State} eq 'COMPLETED', 'The variant search completed.');
## I think the freebayes parser runs into a problem when running on the cluster
## which I have previously associated with not properly copying the inputs
## to the current working directory.  This requires further investigation.
