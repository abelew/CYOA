# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename dirname";
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');

my $start = getcwd();
my $new = 'test_output_async';
mkdir($new);
chdir($new);

my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $phix_fasta = qq"${start_dir}/genome/phix.fasta";
my $phix_gff = qq"${start_dir}/genome/phix.gff";

my $cyoa = Bio::Adventure->new(
    basedir => cwd(),
    libdir => cwd(),
    species => 'phix',
    gff_tag => 'ID',
    gff_type => 'CDS',
    stranded => 'no',
    vcf_cutoff => 1,);
my $paths = $cyoa->Bio::Adventure::Config::Get_Paths(subroutine => 'Bowtie2');
ok (-d $paths->{index_dir}, qq"The index directory: $paths->{index_dir} exists\n");
my $fasta_dir = dirname($paths->{fasta});
make_path($fasta_dir);
my $gff_dir = dirname($paths->{gff});
make_path($gff_dir);
ok(-d $fasta_dir, qq"The fasta directory exists: ${fasta_dir}.");
ok(-d $gff_dir, qq"The fasta directory exists: ${gff_dir}.");

if (!-r 'test_forward.fastq.gz') {
    ok(cp($input_file, 'test_forward.fastq.gz'), 'Copying data.');
}
if (!-r 'genome/fasta/phix.fasta') {
    ok(cp($phix_fasta, 'genome/fasta/phix.fasta'), 'Copying phix fasta file.');
}
if (!-r 'genome/gff/phix.gff') {
    ok(cp($phix_gff, 'genome/gff/phix.gff'), 'Copying phix gff file.');
}

my $variant = $cyoa->Bio::Adventure::SNP::Align_SNP_Search(
    input => 'test_forward.fastq.gz',
    jprefix => '30',);
ok($variant, 'Submit variant search.');
my $status = $cyoa->Wait(job => $variant);
ok($status->{State} eq 'COMPLETED', 'The variant search completed.');
use Data::Dumper;
print Dumper $variant;
## Some files to check out:
## $variant->{parse}->{output_by_gene}
## $variant->{parse}->{output}  {output_penetrance}  {output_count} {output_genome} {output_types}
## Check that we got some interesting outputs:
ok(-r $variant->{parse}->{output_by_gene}, "Created the variants by gene: $variant->{parse}->{output_by_gene}.\n");
ok(-r $variant->{parse}->{output}, "Created the parsed variants: $variant->{parse}->{output}.\n");
ok(-r $variant->{parse}->{output_penetrance}, "Created the penetrance file: $variant->{parse}->{output_penetrance}.\n");
ok(-r $variant->{parse}->{output_genome}, "Created the new genome: $variant->{parse}->{output_genome}.\n");
ok(-r $variant->{parse}->{output_types}, "Created the types: $variant->{parse}->{output_types}.\n");
chdir($start);
