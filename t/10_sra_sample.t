# -*-Perl-*-
use strict;
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename dirname";
use File::Copy qw"cp mv";
use File::Path qw"make_path remove_tree rmtree";
use File::ShareDir qw"dist_dir dist_file module_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
use Test::More qw"no_plan";

## Let us test the ability to download a sample from the SRA, download the genome/annotations from
## either NCBI and/or Ensembl and/or TriTrypDB,
## then perform trimming etc.

my $sra_accession = 'SRR2136722';
my $sample_accession = 'SAMN03945365';
my $species = 'lmajor_friedlin';
## I would like to make this smart enough to not need the version number...
my $assembly_accession = 'GCA_000002725.2';

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    gff_tag => 'gene_id',
    gff_type => 'gene',
    libdir => 'reference',
    species => $species,);

## Start by downloading the assembly
## I want to add a little logic to cache these files.
## Download the assembly via Bio::Adventure::Prepare::Download_NCBIDatasets_Accession

my $gff = qq"${assembly_accession}.gff";
my $genome = qq"${assembly_accession}.fsa";
my $genbank = qq"${assembly_accession}.gbff.gz";
my $cds = qq"${assembly_accession}.ffn";
my $amino_acids = qq"${assembly_accession}.faa";
my $md5_file = qq"md5sum.txt";

if (-d 'reference/genome/gff') {
    print "The gff directory exists.\n";
} else {
    make_path('reference/genome/gff');
}
if (-d 'reference/genome/fasta') {
    print "The fasta directory exists.\n";
} else {
    make_path('reference/genome/fasta');
}
if (-d 'reference/genome/genbank') {
    print "The genbank directory exists.\n";
} else {
    make_path('reference/genome/genbank');
}
if (-d 'reference/genome/amino_acids') {
    print "The amino acid directory exists.\n";
} else {
    make_path('reference/genome/amino_acids');
}
if (-d 'reference/genome/CDS') {
    print "The CDS directory exists.\n";
} else {
    make_path('reference/genome/CDS');
}
if (-d 'reference/rRNA/gff') {
    print "The ribosomal gff directory exists.\n";
} else {
    make_path('reference/rRNA/gff');
}
if (-d 'reference/rRNA/fasta') {
    print "The rRNA fasta directory exists.\n";
} else {
    make_path('reference/rRNA/fasta');
}

my $final_gff = qq"reference/genome/gff/${species}.gff";
my $final_genome = qq"reference/genome/fasta/${species}.fasta";
my $final_genbank = qq"reference/genome/genbank/${species}.gb.gz";
my $final_cds = qq"reference/genome/CDS/${species}.fasta";
my $final_amino_acids = qq"reference/genome/amino_acids/${species}.fasta";

if (-f $final_genbank) {
    print "The assembly information has been downloaded to the local reference tree.\n";
} else {
    my $assembly = $cyoa->Bio::Adventure::Prepare::Download_NCBIDatasets_Accession(
        input => $assembly_accession, delete => 0);
## Upon completion we should have a few new files:
##  GCA_000002725.2.gbff.gz, GCA_000002725.2.gff, GCA_000002725.2.faa, GCA_000002725.2.ffn,
##  GCA_000002725.2.fsa, assembly_data_report.jsonl, sequence_report.jsonl, md5sum.txt,
##  dataset_catalog.json
## So, before I am finished with this test suite, I ought to have tests for these files, their
## md5sums, and I should run my genbank converter and see that it returns information which is
## similar/identical to these files downloaded from ncbi.
## For now, I think I will just check that they exist because I am lazy.
    ok(-f $gff, qq"The L.major gff annotations were downloaded to: ${gff}.");
    ok(-f $genome, qq"The L.major genome was downloaded to: ${genome}.");
    ok(-f $genbank, qq"The L.major genbank flat file was downloaded to: ${genbank}.");
    ok(-f $cds, qq"The L.major coding sequences were downloaded to: ${cds}.");
    ok(-f $amino_acids, qq"The L.major amino acids were downloaded to: ${amino_acids}.");
    ok(-f $md5_file, qq"The checksums were downloaded to: ${md5_file}.");
    ## Move the assembly files to reference/
    my $moved_gff = mv($gff, $final_gff);
    ok($moved_gff, 'Moved the gff file to the reference tree.');
    my $moved_genome = mv($genome, $final_genome);
    ok($moved_genome, 'Moved the genome to the reference tree.');
    my $moved_genbank = mv($genbank, $final_genbank);
    ok($moved_genbank, 'Moved the genbank file to the reference tree.');
    my $moved_cds = mv($cds, $final_cds);
    ok($moved_cds, 'Moved the cds to the reference tree.');
    my $moved_aa = mv($amino_acids, $final_amino_acids);
    ok($moved_aa, 'Moved the amino acids to the reference tree.');
}

## Let us try and extract the rRNA from the assembly.
my $extracted = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => $final_genbank);
## TODO: I have been meaning for a very long time to send the outputs
## from this function to outputs/xxgb2gff/ but have not done it
## and so all the files get dumped to reference/genome/genbank
## this is bad and needs to be fixed.

ok($extracted, 'Extracted information from the genbank flat files.');
make_path('gb2gff');
my $rrna_fasta = 'reference/genome/genbank/lmajor_friedlin_rrna.fasta';
my $rrna_fasta_final = 'reference/rRNA/fasta/lmajor_friedlin.fasta';
ok(mv($rrna_fasta, $rrna_fasta_final));
my $rrna_gff = 'reference/genome/genbank/lmajor_friedlin_rrna.gff';
my $rrna_gff_final = 'reference/rRNA/gff/lmajor_friedlin.gff';
ok(mv($rrna_gff, $rrna_gff_final));
my $all_gff = 'reference/genome/genbank/lmajor_friedlin_all.gff';
my $all_gff_final = 'gb2gff/all.gff';
ok(mv($all_gff, $all_gff_final));
my $cds_gff = 'reference/genome/genbank/lmajor_friedlin_cds.gff';
my $cds_gff_final = 'gb2gff/cds.gff';
ok(mv($cds_gff, $cds_gff_final));
my $aa = 'reference/genome/genbank/lmajor_friedlin.faa';
my $aa_final = 'gb2gff/aa.fasta';
ok(mv($aa, $aa_final));
my $ffn = 'reference/genome/genbank/lmajor_friedlin.ffn';
my $ffn_final = 'gb2gff/cds.ffn';
ok(mv($ffn, $ffn_final));
my $fsa = 'reference/genome/genbank/lmajor_friedlin.fsa';
my $fsa_final = 'gb2gff/genome.fsa';
ok(mv($fsa, $fsa_final));
my $converted_gff = 'reference/genome/genbank/lmajor_friedlin.gff';
my $converted_gff_final = 'gb2gff/genome.gff';
ok(mv($converted_gff, $converted_gff_final));
my $inter_gff = 'reference/genome/genbank/lmajor_friedlin_interCDS.gff';
my $inter_final = 'gb2gff/inter.gff';
ok(mv($inter_gff, $inter_final));
my $misc = 'reference/genome/genbank/lmajor_friedlin_miscRNA.fasta';
my $misc_final = 'gb2gff/misc.fasta';
ok(mv($misc, $misc_final));
my $misc_gff = 'reference/genome/genbank/lmajor_friedlin_miscRNA.gff';
my $misc_gff_final = 'gb2gff/misc.gff';
ok(mv($misc_gff, $misc_gff_final));
my $nc = 'reference/genome/genbank/lmajor_friedlin_ncRNA.fasta';
my $nc_final = 'gb2gff/nc.fasta';
ok(mv($nc, $nc_final));
my $nc_gff = 'reference/genome/genbank/lmajor_friedlin_ncRNA.gff';
my $nc_gff_final = 'gb2gff/nc.gff';
ok(mv($nc_gff, $nc_gff_final));

## Now download the sample of interest from the sra:
#my $r1 = qq"${sample_accession}/${sample_accession}_pass_1.fastq.gz";
#my $r2 = qq"${sample_accession}/${sample_accession}_pass_2.fastq.gz";
#if (-f $r1) {
#    print "The SRA accession has already been downloaded.\n";
#} else {
#    my $sra = $cyoa->Bio::Adventure::Prepare::Fastq_Dump(input => $sample_accession);
#    ok($sra, 'Able to run fastq-dump to download the test accession.');
#    ok(-f $r1, 'Successfully downloaded read1.');
#    ok(-f $r2, 'Successfully downloaded read2.');
#}

## Trim the sample
#my $input_fastq = qq"${r1}:${r2}";
#my $trimomatic = $cyoa->Bio::Adventure::Trim::Trimomatic(input => $input_fastq,
#                                                         jprefix => '10_1');
#ok($trimomatic, qq"Run trimomatic on ${sample_accession}");
#my $fastp = $cyoa->Bio::Adventure::Trim::Fastp(input => $input_fastq);
#ok($fastp, qq"Run fastp on ${sample_accession}");

#my $fastp_output = $fastp->{main};
## Map the sample via hisat and salmon
## my $hisat = $cyoa->Bio::Adventure::Map::Hisat(input => $fastp_output,);
## my $salmon = $cyoa->Bio::Adventure::Map::Salmon(input => $fastp_output,);


#my $hisat = $cyoa->Bio::Adventure::Map::Hisat2(
#    input => $local_input_fastq,
#    prefix => '31_1',);
#ok($hisat, 'Run Hisat2 with single-end input.');
#my $sam_file = $hisat->{samtools}->{main};
#my $count_file = $hisat->{counter}->{main};
#my $stats_file = $hisat->{stats}->{main};
#ok(-f $sam_file, qq"The sorted bamfile was created to ${sam_file}.");
#ok(-f $count_file, qq"The count table was created: ${count_file}.");
#ok(-f $stats_file, qq"The hisat stats were recorded to ${stats_file}.");
#my $actual = $cyoa->Last_Stat(input => $stats_file);
#ok($actual, 'Collect Hisat Statistics');
#my $expected = qq"hisat_phix_genome_test.stderr,10000,49,9951,0";
#unless(ok($expected eq $actual, 'Are the hisat stats as expected?')) {
#    my ($old, $new) = diff($expected, $actual);
#    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
#}
## Stupid featureCounts puts a random space at the end of the first line
#$expected = qq!"1"
#Length
#1542
#1026
#363
#171
#261
#3618
#634
#276
#117
#528
#987
#!;
#my $cmd = qq!xzcat ${count_file} | awk '{print \$6}'!;
#print "Running
#${cmd}\n";
#$actual = qx"$cmd";
#unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
#    my($old, $new) = diff($expected, $actual);
#    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
#}

chdir($start);
