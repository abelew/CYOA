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
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $phix_fasta = qq"${start_dir}/genome/phix.fasta";
my $phix_gff = qq"${start_dir}/genome/phix.gff";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

make_path('genome/indexes'); ## Make a directory for the phix indexes.
if (!-r 'test_forward.fastq.gz') {
    ok(cp($input_file, 'test_forward.fastq.gz'), 'Copying data.');
}
if (!-r 'genome/phix.fasta') {
    ok(cp($phix_fasta, 'genome/phix.fasta'), 'Copying phix fasta file.');
}
if (!-r 'genome/phix.gff') {
    ok(cp($phix_gff, 'genome/phix.gff'), 'Copying phix gff file.');
}

my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    libdir => cwd(),
    species => 'phix',
    gff_tag => 'ID',
    gff_type => 'CDS',
    stranded => 'no',
    vcf_cutoff => 1,);
my $index = $cyoa->Bio::Adventure::Index::BT2_Index(input => $phix_fasta,
                                                    jprefix => '30',);
my $variant = $cyoa->Bio::Adventure::SNP::Align_SNP_Search(
    input => 'test_forward.fastq.gz',
    jprefix => '30',);
ok($variant, 'Ran variant search.');
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
