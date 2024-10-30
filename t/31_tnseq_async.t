# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename dirname";
use File::Copy qw"cp copy move mv";
use File::Path qw"make_path remove_tree rmtree";
use File::ShareDir qw"dist_dir dist_file module_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');
my $input = 'consolidated.fastq.xz';
my $input_file = qq"${start_dir}/${input}";
my $phix_fasta = qq"${start_dir}/genome/phix.fastq";
my $phix_gff = qq"${start_dir}/genome/phix.gff";

my $start = getcwd();
my $new = 'test_output_async';
mkdir($new);
chdir($new);

if (!-r $input) {
    ok(copy($input_file, $input), 'Copying data.');
}

my $genome = 'sagalactiae_cjb111';
my $cyoa = Bio::Adventure->new(
    basedir => cwd(),
    libdir => cwd(),
    species => $genome,
    gff_tag => 'locus_tag',
    gff_type => 'gene',
    stranded => 'no',);
my $paths = $cyoa->Bio::Adventure::Config::Get_Paths(subroutine => 'Transit_TPP',);
my $groupb_strep = dist_file('Bio-Adventure', qq"genome/${genome}.gb.xz");
if (!-r qq"genome/${genome}.fasta") {
    ok(copy($groupb_strep, qq"genome/${genome}.gb.xz"), 'Copying group b genome.');
}
my $groupb_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => qq"genome/${genome}.gb.xz", jprefix => '31', output_dir => 'genome',);
my $status = $cyoa->Wait(job => $groupb_convert);
ok($status->{State} eq 'COMPLETED', 'The conversion from genbank to fasta/gff completed.');
ok($groupb_convert, 'Converted the group B strep genbank to fasta/gff/etc.');
ok(-r $groupb_convert->{output_all_gff}, 'Found the newly created gff.');
ok(-r $groupb_convert->{output_fasta}, 'Found the newly created fsa.');
my $moved = move($groupb_convert->{output_fasta}, $paths->{fasta});
ok($moved, qq"Moved the ${genome}.fsa to $paths->{fasta}.");
$moved = move($groupb_convert->{output_all_gff}, $paths->{gff});
ok($moved, qq"Moved ${genome}_all.gff to $paths->{gff}.");

my $tpp = $cyoa->Bio::Adventure::TNSeq::Transit_TPP(
    input => $input,
    species => $genome,
    protocol => 'Sassetti',
    jprefix => '31',);
$status = $cyoa->Wait(job => $tpp);
ok($status->{State} eq 'COMPLETED', 'The transit preprocessing completed.');

ok(-r $tpp->{samtools}->{output}, qq"Found the samtools bam file: $tpp->{samtools}->{output}.");
ok(-r $tpp->{htseq}->[0]->{output}, qq"Found the htseq count table: $tpp->{htseq}->[0]->{output}.");
ok(-r $tpp->{output_r1}, qq"Found the r1 reads: $tpp->{output_r1}.");
ok(-r $tpp->{output_r1tr}, qq"Found the r1 trimmed reads: $tpp->{output_r1tr}.");
ok(-r $tpp->{output_count}, qq"Found the tpp counts: $tpp->{output_count}.");
ok(-r $tpp->{output_wig}, qq"Found the tpp wig: $tpp->{output_wig}.");
ok(-r $tpp->{output_bam}, qq"Found the tpp bam: $tpp->{output_bam}.");
ok(-r $tpp->{output_stats}, qq"Found the tpp stats file: $tpp->{output_stats}.");
chdir($start);
