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
my $input_r1 = qq"${start_dir}/r1.fastq.xz";
my $input_r2 = qq"${start_dir}/r2.fastq.xz";
my $phix_fasta = qq"${start_dir}/genome/phix.fastq";
my $phix_gff = qq"${start_dir}/genome/phix.gff";
my $terminase_db = qq"${start_dir}/genome/phage_terminases.fasta";

my $start = getcwd();
my $a = 'test_output_async';
my $actual = '';
my $expected1 = '';
my $expected2 = '';
my $test_file = '';
mkdir($a);
chdir($a);

my $cyoa = Bio::Adventure->new(basedir => cwd(), libdir => cwd());
my $paths = $cyoa->Bio::Adventure::Config::Get_Paths();
## Copy the reads for running the tests.
ok(cp($input_r1, 'r1.fastq.xz'), 'Copying r1.') if (!-r 'r1.fastq.xz');
ok(cp($input_r2, 'r2.fastq.xz'), 'Copying r2.') if (!-r 'r2.fastq.xz');
ok(cp($terminase_db, "$paths->{blast_dir}/terminase.fasta"), 'Copying terminase db.') if (!-r "$paths->{blast_dir}/terminase.fasta");
## Invoke the pipeline, keep it within our test directory with basedir.

my $assemble = $cyoa->Bio::Adventure::Pipeline::Phage_Assemble(
    input => 'r1.fastq.xz:r2.fastq.xz', jprefix => '50');
my $job_id;
my $status;
## Check the trimomatic output.
my $id = '01trim';
$test_file = $assemble->{$id}->{stderr};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
my $comparison = ok(-f $test_file, qq"Checking trimomatic output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"tail -n 3 ${test_file}";
$expected1 = qq"ILLUMINACLIP: Using 1 prefix pairs, 19 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Read Pairs: 108212 Both Surviving: 94183 (87.04%) Forward Only Surviving: 10743 (9.93%) Reverse Only Surviving: 1074 (0.99%) Dropped: 2212 (2.04%)
TrimmomaticPE: Completed successfully
";
$comparison = ok($expected1 eq $actual, 'Checking trimomatic result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Look at the fastqc outputs

$id = '02fastqc';
$test_file = $assemble->{$id}->{txtfile};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking fastqc output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file}";
$expected1 = qq"PASS\tBasic Statistics\tr1-trimmed.fastq
PASS\tPer base sequence quality\tr1-trimmed.fastq
WARN\tPer tile sequence quality\tr1-trimmed.fastq
PASS\tPer sequence quality scores\tr1-trimmed.fastq
WARN\tPer base sequence content\tr1-trimmed.fastq
WARN\tPer sequence GC content\tr1-trimmed.fastq
PASS\tPer base N content\tr1-trimmed.fastq
WARN\tSequence Length Distribution\tr1-trimmed.fastq
WARN\tSequence Duplication Levels\tr1-trimmed.fastq
WARN\tOverrepresented sequences\tr1-trimmed.fastq
PASS\tAdapter Content\tr1-trimmed.fastq
";
#$comparison = ok($expected1 eq $actual, 'Checking fastqc result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected1, $actual);
#    diag("-- expected1\n${e}\n-- actual\n${a}\n");
#}

## Check out the results from RACER, this assumes the output compile flag is set.
$id = '03correction';
$test_file = $assemble->{$id}->{stdout};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking racer output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"grep changed ${test_file} | head";
$expected1 = qq"Number of changed positions\t\t58657
Number of changed positions after\t4372
Number of changed positions before\t54285
Number of changed positions\t\t423
Number of changed positions after\t418
Number of changed positions before\t5
Number of changed positions\t\t133
Number of changed positions after\t129
Number of changed positions before\t4
Number of changed positions\t\t90
";
#$comparison = ok($expected1 eq $actual, 'Checking racer result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected1, $actual);
#    diag("-- expected1\n${e}\n-- actual\n${a}\n");
#}

## Look at the standard kraken report.
$id = '04kraken_standard';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking standard kraken report: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 1 ${test_file}";
$expected1 = qq"d__Bacteria\t5003
d__Bacteria|p__Proteobacteria\t4884
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria\t4748
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales\t3523
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae\t3096
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia\t3023
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia marmotae\t3023
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Citrobacter\t72
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Citrobacter|s__Citrobacter portucalensis\t72
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Morganellaceae\t425
";
#$comparison = ok($expected1 eq $actual, 'Checking kraken standard result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected1, $actual);
#    diag("-- expected1\n${e}\n-- actual\n${a}\n");
#}

## See if the kraken-based filter worked.  The variable 'log' is already taken, so switch it out.
$id = '05host_filter';
$test_file = $assemble->{$id}->{log};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking kraken filter output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"tail -n 2 ${test_file}";
$expected1 = qq"Filtering out reads which map to GCF_002900365.1.
Symlinking final output files to outputs/05filter_kraken_host
";
$comparison = ok($expected1 eq $actual, 'Checking kraken filter result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Look at the viral kraken results.
$id = '06kraken_viral';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking kraken viral report: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected1 = qq"d__Viruses\t90498
d__Viruses|k__Heunggongvirae\t90498
d__Viruses|k__Heunggongvirae|p__Uroviricota\t90498
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes\t90498
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales\t90498
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae\t90477
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae|g__Kayfunavirus\t90464
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae|g__Kayfunavirus|s__Escherichia virus EcoDS1\t15752
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae|g__Kayfunavirus|s__Citrobacter virus CR44b\t10348
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae|g__Kayfunavirus|s__Escherichia virus ST31\t9695
";
#$comparison = ok($expected1 eq $actual, 'Checking kraken viral result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected1, $actual);
#    diag("-- expected1\n${e}\n-- actual\n${a}\n");
#}

## Check the unicycler outputs.
$id = '07assembly';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking unicycler fasta output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 1 ${test_file}";
$expected1 = qq">1 length=40082 depth=1.00x circular=true
";
$comparison = ok($expected1 eq $actual, 'Was the assembly created with expected1 size/depth?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Look at the filter_depth results.
$id = '08depth_filter';
$test_file = $assemble->{$id}->{output_log};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking depth filter output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"more ${test_file}";
$expected1 = qq"Starting depth coverage filter of outputs/57unicycler/test_output_final_assembly.fasta.
Any contigs with a coverage ratio vs. the highest coverage of < 0.2 will be dropped.
Writing filtered contigs to outputs/08filter_depth/final_assembly.fasta
The range of observed coverages is 1.00 <= x <= 1
Writing 1 with normalized coverage: 1
";
$comparison = ok($expected1 eq $actual, 'Is the depth log as expected1?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Check the phastaf results.
$id = '09phastaf';
$test_file = $assemble->{$id}->{coordinates};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking phastaf coordinate file: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file} | awk '{print \$4}'";
$expected1 = qq"PHAGE_Cronob_Dev2_NC_023558-gi|651729639|ref|YP_009005149.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729605|ref|YP_009005115.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729621|ref|YP_009005131.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729638|ref|YP_009005148.1|
PHAGE_Citrob_CR44b_NC_023576-gi|589286983|ref|YP_009007177.1|
PHAGE_Escher_vB_EcoP_GA2A_NC_031943-gi|100054|ref|YP_009324988.1|
PHAGE_Entero_K1F_NC_007456-gi|77118185|ref|YP_338107.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729630|ref|YP_009005140.1|
PHAGE_Citrob_CR44b_NC_023576-gi|589286980|ref|YP_009007173.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729625|ref|YP_009005135.1|
";
$comparison = ok($expected1 eq $actual, 'Checking phastaf result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Check the ICTV results.
$id = '10ictv';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking ICTV classifier output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 3 ${test_file} | awk '{print \$2}'";
$expected1 = qq"query_description
length=40082
length=40082
";
$comparison = ok($expected1 eq $actual, 'Checking ICTV classifier result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Check the Rosalindplus results.
$id = '11rosalindplus';
$test_file = $assemble->{$id}->{job_log};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking Rosalindplus output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 2 ${test_file}";
$expected1 = qq"Counting ORFs on the current Rosalind and Franklin strand.
If the Franklin strand is larger, flipping them.
";
$comparison = ok($expected1 eq $actual, 'Is the rosalindplus log as expected1?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Check the phageterm results.
## $test_file = 'outputs/12phageterm_11rosalindplus/contig1_nrt.txt';
$id = '12phageterm';
$test_file = $assemble->{$id}->{output_type};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking phageterm output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"more ${test_file}";
$expected1 = qq"DTR (short)";
$comparison = ok($expected1 eq $actual, 'Are the phageterm results as expected1?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Check the coverage results.
## $test_file = 'outputs/13assembly_coverage_test_output/coverage.tsv';
$id = '13coverage';
$test_file = $assemble->{$id}->{output_tsv};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking coverage output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"more ${test_file} | awk '{print \$3}'";
$expected1 = qq"Length
40082
";
$comparison = ok($expected1 eq $actual, 'Checking coverage result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Check the terminase search results.
## $test_file = 'outputs/14termreorder_12phageterm_11rosalindplus/terminase_summary.tsv';
$id = '14terminase_reorder';
$test_file = $assemble->{$id}->{output_tsv};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking terminase search output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 3 ${test_file}";
$expected1 = qq"Name\tLength\tAccession\tDescription\tScore\tSignificance\tBit\tHitStrand\tQueryStrand
AZS06569.1\t551\tAZS06569\tterminase [Mycobacterium phage JacoRen57]\t26.4\t2.6\t26.4\t0\t1
AUV61411.1\t503\tAUV61411\tlarge terminase [Pontimonas phage phiPsal1]\t32.5\t0.24\t32.5\t0\t1
";
#$comparison = ok($expected1 eq $actual, 'Checking terminase search result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected1, $actual);
#    diag("-- expected1\n${e}\n-- actual\n${a}\n");
#}

## Check the prokka run.
## $test_file = 'outputs/15prokka/test_output.ffn';
$id = '15prokka';
$test_file = $assemble->{$id}->{output_cds};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking prokka output file: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected1 = qq">test_output_async_00001 hypothetical protein
ATGGAACGTAACGCTGACGCATACTATGAGCTGCTGAATGCAACCGTTAAAGCATTTAAC
GAGCGTGTTCAGTACGACGAAATAGCTAAAGGTGATGACTACCATGATGCGCTGCATGAA
GTCGTAGACGGTCAGGTTCCGCACTATTACCACGAGATCTTCACGGTGATGGCTGCTGAT
GGTATTGACATTGAGTTTGAAGACTCTGGGCTGATGCCTGAGACCAAGGACGTAACGCGC
ATACTGCAAGCTCGCATCTATGAGGCACTGTATAACGGCGTGTCTAATAGCTCGGATGTG
GTCTGGTTTGAGGCTGAAGAGAGCGACGAAGAGGGTAAGTATTGGGTAGTTGACGCTAAA
ACGGGACTATTCGCTGAGCAAGCTATACCTCTTGAGGTCGCTATTGCATCTGCCAAAGAC
CTCTATGCGGTAGGTCATCACATGAAAGTCGAAGACATTAACGATAACGTAGTGTTCGAC
CCTGCGGCTGAAGAGGACTGCGAGTGA
";
$comparison = ok($expected1 eq $actual, 'Checking prokka result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Check the prodigal run.
## $test_file = 'outputs/16prodigal_test_output.fna/predicted_cds.fasta';
$id = '16prodigal';
$test_file = $assemble->{$id}->{output_cds};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking prodigal output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 7 ${test_file}";
$expected1 = qq">gnl|Prokka|test_output_async_1_1 # 1267 # 1773 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.489
ATGGAACGTAACGCTGACGCATACTATGAGCTGCTGAATGCAACCGTTAAAGCATTTAACGAGCGTGTTC
AGTACGACGAAATAGCTAAAGGTGATGACTACCATGATGCGCTGCATGAAGTCGTAGACGGTCAGGTTCC
GCACTATTACCACGAGATCTTCACGGTGATGGCTGCTGATGGTATTGACATTGAGTTTGAAGACTCTGGG
CTGATGCCTGAGACCAAGGACGTAACGCGCATACTGCAAGCTCGCATCTATGAGGCACTGTATAACGGCG
TGTCTAATAGCTCGGATGTGGTCTGGTTTGAGGCTGAAGAGAGCGACGAAGAGGGTAAGTATTGGGTAGT
TGACGCTAAAACGGGACTATTCGCTGAGCAAGCTATACCTCTTGAGGTCGCTATTGCATCTGCCAAAGAC
";
$comparison = ok($expected1 eq $actual, 'Checking prodigal CDS predictions:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## It appears that the x significant digits in these various outputs are going
## to be a big pita for me.  In this test, it failed due to the difference between
## 9.70 and 9.71 (orf00011)

## Check the glimmer run.
## $test_file = 'outputs/17glimmer/glimmer3.predict';
$id = '17glimmer';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking glimmer result file: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected1 = qq">gnl|Prokka|test_output_async_1
orf00003      131      265  +2     0.84
orf00005      352      305  -2     1.29
orf00010     1052     1129  +2     7.87
orf00011     1178     1225  +2     9.70
orf00013     1267     1773  +1    14.08
orf00015     1821     1928  +3     2.10
orf00016     1931     2056  +2     3.05
orf00019     2262     2411  +3     4.05
orf00020     2447     2737  +2     2.27
";
$comparison = ok($expected1 eq $actual, 'Checking glimmer CDS predictions:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## In this instance, the significant digit differences are in the range of
## 1e-5, I am thinking I really do not care.  I am thinking I will just pull
## the first column from the output and assume the rest is good enough.

## Check the phanotate run.
## $test_file = 'outputs/18phanotate/test_output_phanotate.tsv.xz';
$id = '18phanotate';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking phanotate output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file} | head | awk '{print \$1}'";
## Different versions of phanotate give slightly different outputs...
my $expected1_first = qq"#id:
#START
1
183
477
1148
1267
1773
1931
2128
";
my $expected1_second = qq"#id:
#START
<1
183
477
1148
1267
1773
1931
2128
";
## Here is the previous result, which just used head -n 5 for no good reason.
## #id:\tgnl|Prokka|test_output_1
## #START\tSTOP\tFRAME\tCONTIG\tSCORE
## 1\t117\t+\tgnl|Prokka|test_output_1\t-1.248555528686707940866691777\t
## 183\t302\t+\tgnl|Prokka|test_output_1\t-0.2175130562377134455954126775\t
## 477\t617\t+\tgnl|Prokka|test_output_1\t-0.07018008835792925643848556090\t
my $something_good = 0;
if ($expected1_first eq $actual) {
    $something_good++;
}
if ($expected1_second eq $actual) {
    $something_good++;
}
$comparison = ok($something_good, 'Checking phanotate output:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1_first, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Check the merge_cds.
## $test_file = 'outputs/19merge_cds_predictions/test_output.tsv';

## This one also changed from -5930 to -5940 (test_output_0007)

$id = '19cds_merge';
$test_file = $assemble->{$id}->{output_tsv};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking CDS merge output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected1 = qq"locus_tag	contig	type	source	start	end	strand	cds_prediciton	aa_sequence
test_output_async_0001	test_output_async_1	CDS		131	265	1	glimmer	VVVETIGWDYWLSLSLLLAAGVTAGSQWVGWVETLVCSLVSQCN
test_output_async_0002	test_output_async_1	CDS		183	302	1	phanotate, score: -0.218	LLLALLLEVSGSGGSRLSYALWSLSVINAIMVTIHERKT
test_output_async_0003	test_output_async_1	CDS		305	352	-1	glimmer	LTTVAKVSRVASAMN
test_output_async_0004	test_output_async_1	CDS		477	617	1	phanotate, score: -0.0702	LDQKFETTSHSSRTSSLPIGPLSVQTKGPTPVYHKVGPMVKTSGQR
test_output_async_0005	test_output_async_1	CDS		888	1148	-1	phanotate, score: -0.107	LLKSIPFSQRTSGRPVQCWSPPLLRCGTAYISSLLLVNYFLSSACCSYDLSGCLLNRDDPASSLSGCCRVVLTEAIKPQSRPIVNM
test_output_async_0006	test_output_async_1	CDS		1178	1225	1	glimmer	VINYRVFESTPEGPD
test_output_async_0007	test_output_async_1	CDS		1267	1773	1	phanotate, score: -5940	MERNADAYYELLNATVKAFNERVQYDEIAKGDDYHDALHEVVDGQVPHYYHEIFTVMAADGIDIEFEDSGLMPETKDVTRILQARIYEALYNGVSNSSDVVWFEAEESDEEGKYWVVDAKTGLFAEQAIPLEVAIASAKDLYAVGHHMKVEDINDNVVFDPAAEEDCE
test_output_async_0008	test_output_async_1	CDS		1773	1928	1	phanotate, score: -6.33	MVTYGLCQHHVTNARIMVKTGQLNHDATMCLLKAVYEGRKLIHNSLHAEDK
test_output_async_0009	test_output_async_1	CDS		1931	2056	1	phanotate, score: -4.91	MYQITYNSEQAFYEGCYEMMKRGACYVANHHSLTITLTGGY
";
$comparison = ok($expected1 eq $actual, 'Checking CDS merge result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## I messed up the jellyfish runs when standardizing the style.
## This test is therefore broken until I rerun.
## Check jellyfish
## $test_file = 'outputs/20jellyfish_test_output/test_output_9.hist.xz';
$id = '20jellyfish';
$test_file = $assemble->{$id}->{histogram_file};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking jellyfish output tsv: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file}";
$expected1 = qq"1 28869
2 4386
3 663
4 253
5 67
6 17
7 8
8 2
9 2
10 5
";
$expected2 = qq"# 1 37998
# 2 1153
# 3 42
# 4 177
# 5 9
# 6 1
# 7 2
# 9 2
# 10 3
";
$comparison = ok($expected1 eq $actual, 'Checking jellyfish result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Check aragorn
## $test_file = 'outputs/21aragorn/aragorn.txt';
$id = '21aragorn';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking aragorn output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"more ${test_file}";
$expected1 = qq">test_output_async_1
0 genes found
";
$comparison = ok($expected1 eq $actual, 'Checking aragorn result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## I also broke the trnascan run...
## $test_file = 'outputs/22trnascan/trnascan_relaxed.txt';
$id = '22trnascan';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking trnascan output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"tail -n 4 ${test_file}";
$expected1 = qq"number of sequences= 1
number of bases tested (one strand)=41261
number of bases tested (both strands)= 82522
number of predicted tRNA=236
";
$comparison = ok($expected1 eq $actual, 'Checking trnascan result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Trinotate
## $test_file = 'outputs/23trinotate19merge_cds_predictions/test_output.ffn.tsv';
$id = '23trinotate';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking trinotate output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 3 ${test_file} | awk '{print \$1}'";
$expected1 = qq"#gene_id
test_output_0001
test_output_0002
";
$comparison = ok($expected1 eq $actual, 'Checking trinotate results:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Abricate
## $test_file = 'outputs/24abricate_19merge_cds_predictions/abricate_summary.txt';
$id = '24abricate';
$test_file = $assemble->{$id}->{output_txt};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking abricate result: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"more ${test_file}";
$expected1 = qq"#FILE\tNUM_FOUND
outputs/74abricate_69merge_cds_predictions/abricate_argannot.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_card.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_combined.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_dbeth.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_ecoh.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_ecoli_vf.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_megares.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_mvir.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_ncbi.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_plasmidfinder.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_resfinder.tsv\t0
outputs/74abricate_69merge_cds_predictions/abricate_vfdb.tsv\t0
";
$comparison = ok($expected1 eq $actual, 'Checking abricate output:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## It turns out that every invocation of interproscan is different in pretty much every file...
## interproscan
## $test_file = 'outputs/25interproscan_19merge_cds_predictions/test_output.faa.gff3';
$id = '25interproscan';
$test_file = $assemble->{$id}->{output_tsv};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking interproscan output tsv: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"awk '{print \$1}' $test_file | sort | uniq | head -n 1";
$expected1 = qq"test_output_async_0001
";
$comparison = ok($expected1 eq $actual, 'Checking interproscan result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## merge annotations 1
## $test_file = 'outputs/26mergeannot/test_output_runlog.txt';
$id = '26merge_qualities';
$test_file = $assemble->{$id}->{output_log};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking merge_annotations output log: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected1 = qq"Merging annotations and writing new output files:
gbf: outputs/26mergeannot/test_output.gbf, tbl: outputs/26mergeannot/test_output.tbl, xlsx: outputs/26mergeannot/test_output.xlsx.
Reading tsv data from outputs/19merge_cds_predictions/test_output.tsv to start.
Checking for ICTV classification data from outputs/10classify_08filter_depth/ictv_filtered.tsv.
Wrote outputs/26mergeannot/test_output.sbt with variables filled in.
Adding trinotate annotations from outputs/23trinotate19merge_cds_predictions/test_output.tsv.
Adding interproscan annotations from outputs/25interproscan_19merge_cds_predictions/interproscan.tsv.
Adding abricate annotations from outputs/24abricate_19merge_cds_predictions/abricate_combined.tsv.
Got DTR type: DTR (short).
Adding phageterm DTRs.
";
#$comparison = ok($expected1 eq $actual, 'Did we get expected1 merge output logs:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected1, $actual);
#    diag("-- expected1\n${e}\n-- actual\n${a}\n");
#}

## merge annotations 2
$id = '27merge_unmodified';
$test_file = $assemble->{$id}->{output_log};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking merge_annotations stripped log: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected1 = qq"Merging annotations and writing new output files:
gbf: outputs/26mergeannot/test_output.gbf, tbl: outputs/26mergeannot/test_output.tbl, xlsx: outputs/26mergeannot/test_output.xlsx.
Reading tsv data from outputs/19merge_cds_predictions/test_output.tsv to start.
Checking for ICTV classification data from outputs/10classify_08filter_depth/ictv_filtered.tsv.
Wrote outputs/26mergeannot/test_output.sbt with variables filled in.
Adding trinotate annotations from outputs/23trinotate19merge_cds_predictions/test_output.tsv.
Adding interproscan annotations from outputs/25interproscan_19merge_cds_predictions/interproscan.tsv.
Adding abricate annotations from outputs/24abricate_19merge_cds_predictions/abricate_combined.tsv.
Got DTR type: DTR (short).
Adding phageterm DTRs.
";
#$comparison = ok($expected1 eq $actual, 'Checking merge_annotations run log result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected1, $actual);
#    diag("-- expected1\n${e}\n-- actual\n${a}\n");
#}

## Something something cgview...
$id = '28cgview';
$test_file = $assemble->{$id}->{output_xml};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking cgview output xml: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 2 ${test_file}";
$expected1 = qq!<?xml version="1.0" encoding="ISO-8859-1"?>
<cgview backboneRadius="4000" backboneColor="rgb(102,102,102)" backboneThickness="40" featureSlotSpacing="10" labelLineLength="450" labelPlacementQuality="better" labelLineThickness="12" rulerPadding="130" tickThickness="18" shortTickThickness="18" arrowheadLength="60" rulerFont="SansSerif, plain, 130" rulerFontColor="rgb(0,0,0)" labelFont="SansSerif, plain, 130" isLinear="true" minimumFeatureLength="1.0" sequenceLength="41261" height="10000" width="10000" globalLabel="true" moveInnerLabelsToOuter="false" featureThickness="86.54" tickLength="45" useInnerLabels="true" shortTickColor="rgb(0,51,0)" longTickColor="rgb(0,51,0)" zeroTickColor="rgb(0,51,0)" showBorder="true" borderColor="black" backgroundColor="white" tickDensity="0.5">
!;
$comparison = ok($expected1 eq $actual, 'Did we get expected1 cgview output?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## This is not working currently.
$id = '29rnafold';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking rnafold output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file} | awk '{print \$1}' | head";
$expected1 = qq!contig
test_output_async_1
test_output_async_1
test_output_async_1
test_output_async_1
test_output_async_1
test_output_async_1
test_output_async_1
test_output_async_1
test_output_async_1
!;
$comparison = ok($expected1 eq $actual, 'Checking RNAfold output:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Restriction catalog...
$id = '30research';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"The restriction endonuclease catalog exists: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"sort ${test_file} | head";
$expected1 = qq!RE	Site	Overhang	Cuts
AasI	GACNNNN^NNGTC	NN	30
AatI	AGG^CCT		1
AatII	GACGT^C	ACGT	13
AauI	T^GTACA	GTAC	9
Acc113I	AGT^ACT		1
Acc16I	TGC^GCA		10
Acc65I	G^GTACC	GTAC	1
AccB1I	G^GYRCC	GYRC	13
AccB7I	CCANNNN^NTGG	NNN	3
AccI	GT^MKAC	MK	33
!;
$comparison = ok($expected1 eq $actual, 'Checking retriction enzyme catalog output:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## This is not working currently.
$id = '31caical';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking caical output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected1 = qq!NAME	CAI
>test_output_async_0001	 0.625
>test_output_async_0002	 0.669
>test_output_async_0003	 0.602
>test_output_async_0004	 0.592
>test_output_async_0005	 0.683
>test_output_async_0006	 0.680
>test_output_async_0007	 0.692
>test_output_async_0008	 0.647
>test_output_async_0009	 0.667
!;
$expected2 = qq"# NAME	CAI
>test_output_async_0001	 0.627
>test_output_async_0002	 0.668
>test_output_async_0003	 0.609
>test_output_async_0004	 0.593
>test_output_async_0005	 0.679
>test_output_async_0006	 0.681
>test_output_async_0007	 0.695
>test_output_async_0008	 0.648
>test_output_async_0009	 0.662
";
$comparison = ok(($expected1 eq $actual || $expected2 eq $actual), 'Checking caical results:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Phagepromoter
$id = '32phagepromoter';
$test_file = $assemble->{$id}->{output_fasta};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking phagepromoter output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected1 = qq!>test_output_async_1:4118 host complement(69..98) score=0.716
TTGACCATCGGTCCAACCTTATGATAGACT
>test_output_async_1:4106 host complement(324..352) score=0.898
TTGACTACAGTAGCTAAGGTCAGTAGAGT
>test_output_async_1:4093 host complement(569..598) score=0.716
TTGACCATCGGTCCAACCTTATGATAGACT
>test_output_async_1:52 host (1036..1065) score=0.828
TTGACAAGCAGTAACGATGAGATGTAAGCT
>test_output_async_1:56 host (1134..1160) score=0.51
AATGCTCTTTAACAATCTGGATAAACT
!;
$comparison = ok($expected1 eq $actual, 'Checking phagepromoter result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Rhotermination prediction
$id = '33rhopredict';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking rhotermpredict output file: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file} | awk '{print \$2}'";
$expected1 = qq"Start
0
321
700
928
1631
2293
2854
3125
3622
";
$comparison = ok($expected1 eq $actual, 'Checking rhotermpredict result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

## Bacphlip
$id = '34bacphlip';
$test_file = $assemble->{$id}->{output};
$job_id = $assemble->{$id}->{job_id};
$status = $cyoa->Wait(job => $job_id);
$comparison = ok(-f $test_file, qq"Checking bacphlip output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file}";
$expected1 = qq"	Virulent	Temperate
0	1.0	0.0
";
$comparison = ok($expected1 eq $actual, 'Checking bacphlip result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected1, $actual);
    diag("-- expected1\n${e}\n-- actual\n${a}\n");
}

chdir($start);
