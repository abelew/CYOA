package Bio::Adventure::Structure;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use feature 'try';
no warnings 'experimental::try';

use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Capture::Tiny qw":all";
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Spec;
use File::Path qw"make_path rmtree";
use File::Which qw"which";
use File::ShareDir qw":ALL";
use IPC::Open2;
use JSON;
use Spreadsheet::Read;
use Symbol qw"gensym";

sub Check_Pair {
    my ($class, %args) = @_;
    my $first_id = $args{first};
    my $second_id = $args{second};
    my $root_dir = $args{root_dir};
    my $finished_file = qq"${root_dir}/finished.txt";
    print "Checking for ${first_id}, ${second_id} in ${finished_file}.\n";
    my $found = 0;
    if (-r $finished_file) {
        my $finished_fh = FileHandle->new("<${finished_file}");
        my $check_string = qq"${first_id},${second_id}";
      LOOP: while (my $line = <$finished_fh>) {
            chomp $line;
            if ($line eq $check_string) {
                $found++;
                last LOOP;
            }
        }
        $finished_fh->close();
    }
    else {
        print "Unable to read: $finished_file\n";
    }
    my $ret = undef;
    $ret = $found if ($found);
    return($ret);
}

sub ProteinFold {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jcpu => 8,
        jgpu => 1,
        jmem => 20,
        jwalltime => '4:00:00',
        jprefix => 79,
        jname => 'alphafold',
        libtype => 'protein',
        mode => 'separate',
        required => ['input']);
    my $output_name = basename($options->{input}, ('.gbk', '.fsa', '.fasta',));
    my $paths = $class->Bio::Adventure::Config::Get_Paths(output_name => $output_name);
    my $comment = '## Iterate over a sequence with Alphafold.';
    my $jstring = qq?
use Bio::Adventure::Structure;
\$h->Bio::Adventure::Structure::ProteinFold_Worker(
  input => '$options->{input}',
  libtype => '$options->{libtype}',
  mode => '$options->{mode}',
  output => '$paths->{output}',
  output_dir => '$paths->{output_dir}',
  jprefix => '$options->{jprefix}',
  stdout => '$paths->{stdout}',
  stderr => '$paths->{stderr}',);
?;
    my $folder = $class->Submit(
        input => $options->{input},
        output => $paths->{output},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        libtype => $options->{libtype},
        mode => $options->{mode},
        comment => $comment,
        output_dir => $paths->{output_dir},
        language => 'perl',
        stdout => $paths->{stdout},
        stderr => $paths->{stderr},);
    return($folder);
}

=head2 C<ProteinFold_PairIDs>

There are a few different things this should do:
 1. Read the txt file of IDs and extract the appropriate amino acid sequence from the database.
 2. Check that this specific pair has not already been done or is currently running.
 3. Assuminmg #2 is false, submit a pairwise job with the appropriate json entry
 4. Make an entry stating that this job is currently running.
 5. Upon completion, move that entry to the set of completed entries
 Ideally, I would want to store this queue in a SQL database which I could query from multiple
 clusters around campus, but I do not think that is currently possible due to the vpn constraints.

 These two files, hg_top3.txt and lm_top3.txt comprise the 3 genes most differentially expressed
 in hg38 and lmajor across the corpus of all samples we have which are infected/uninfected.
 Thus, I want to submit 9 jobs: hg1:lm1, hg1:lm2, hg1:lm3, hg2:lm1 ....

cyoa --method proteinpair --input hg_top3.txt:lm_top3.txt --species hg38_114:lmajor_v68

=cut
sub ProteinFold_PairIDs {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jcpu => 8,
        jgpu => 1,
        jmem => 36,
        jwalltime => '12:00:00',
        jprefix => 79,
        jname => 'proteinpair',
        libtype => 'protein',
        keys => 'transcript:gene',
        species => 'hg38_111:lmajor_v68',
        input => 'hs_top3.txt:lm_top3.txt',);
    my ($sp1, $sp2) = split(/$options->{delimiter}/, $options->{species});
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $comment = '## Iterate over pairs of seqids with proteinfold.';
    my $jstring = qq?
use Bio::Adventure::Structure;
\$h->Bio::Adventure::Structure::ProteinFold_PairIDs_Worker(
  input => '$options->{input}',
  libtype => '$options->{libtype}',
  output => '$paths->{output}',
  output_dir => '$paths->{output_dir}',
  species => '$options->{species}',
  stdout => '$paths->{stdout}',
  stderr => '$paths->{stderr}',);
?;
    my $folder = $class->Submit(
        input => $options->{input},
        output => $paths->{output},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jwalltime => $options->{jwalltime},
        jcpu => $options->{jcpu},
        jgpu => $options->{jgpu},
        jmem => $options->{jmem},
        jstring => $jstring,
        libtype => $options->{libtype},
        mode => $options->{mode},
        comment => $comment,
        output_dir => $paths->{output_dir},
        language => 'perl',
        species => $options->{species},
        stdout => $paths->{stdout},
        stderr => $paths->{stderr},);
    return($folder);
}

sub ProteinFold_PairIDs_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jname => 'proteinpair',
        jprefix => 80,
        jcpu => 8,
        jgpu => 1,
        jmem => 36,
        jwalltime => '12:00:00',
        keys => 'transcript:gene',
        species => 'hg38_111:lmajor_v68',
        input => 'hs_top3.txt:lm_top3.txt',);
    $options->{jprefix}++;
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my ($sp1, $sp2) = split(/$options->{delimiter}/, $options->{species});
    my $species1_aa = qq"$options->{libpath}/$options->{libtype}/fasta/${sp1}.fasta";
    my $species2_aa = qq"$options->{libpath}/$options->{libtype}/fasta/${sp2}.fasta";
    print "Full amino acid sequence files:
$species1_aa
$species2_aa
";
    if (!-r $species1_aa) {
        die("Could not open the first species amino acid file: ${species1_aa}.\n");
    }
    if (!-r $species2_aa) {
        die("Could not open the first species amino acid file: ${species2_aa}.\n");
    }
    my ($key1, $key2) = split(/$options->{delimiter}/, $options->{keys});
    ## Read the two amino acid sequence database into a pair of hashes.
    my $sp1_hash = $class->Bio::Adventure::Fa_to_Hash(input => $species1_aa, key => $key1);
    my $sp2_hash = $class->Bio::Adventure::Fa_to_Hash(input => $species2_aa, key => $key2);
    ## These two input files should be 1 entry/line with some ID I can find in the hash keys
    my ($idfile1, $idfile2) = split(/$options->{delimiter}/, $options->{input});
    print "Input files with sequence IDs:
${idfile1}
${idfile2}
";
    if (!-r $idfile1) {
        die("Could not open the first id file: ${idfile1}.\n");
    }
    if (!-r $idfile2) {
        die("Could not open the first id file: ${idfile2}.\n");
    }

    my $id1 = FileHandle->new("<${idfile1}");
    my $id2 = FileHandle->new("<${idfile2}");
    my @ids1 = ();
    my @ids2 = ();
    while (my $l = <$id1>) {
        chomp $l;
        push(@ids1, $l);
    }
    $id1->close();
    while (my $l = <$id2>) {
        chomp $l;
        push(@ids2, $l);
    }
    $id2->close();

    ## Now we have two hashes of sequence information and two arrays of IDs
  FIRST: for my $first (@ids1) {
      SECOND: for my $second (@ids2) {
            ## 1.  Check that this has not been finished
            my $finishedp = $class->Bio::Adventure::Structure::Check_Pair(
                first => $first, second => $second, root_dir => $paths->{output_dir},);
            if (defined($finishedp)) {
                print "This pair has been started already.\n";
                next SECOND;
            }
            print "Looking for $first and $second sequences.\n";
            my $sp1_seq = $sp1_hash->{$first};
            my $sp2_seq = $sp2_hash->{$second};
            my $result = $class->Bio::Adventure::Structure::ProteinFold_JSON_Pairwise_TwoSeq(
                options => $options, first => $sp1_seq, second => $sp2_seq, paths => $paths,);
        }
    }
    print "Finished Iterating over the pairs of IDs.\n";
}

=head2 C<ProteinFold_Worker>

 Does the actual work of submitting queries to alphafold.
 I renamed this to 'ProteinFold' because I want to generalize it to use Boltz2/Rosetta/etc.
 Which in turn suggests it should not have the prefix 'Protein' because I want to also try
 RNA:Protein folds; but for now, whatever

=cut
sub ProteinFold_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jname => 'alphafold',
        jprefix => 80,
        mode => 'separate',
        required => ['input']);
    $options->{jprefix}++;
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $output_name = basename($options->{input}, ('.fasta', '.faa', '.fsa', '.ffn'));
    my $log = qq"$paths->{output_dir}/${output_name}_runlog.txt";
    my $log_fh = FileHandle->new(">${log}");
    print $log_fh "Setting up an alphafold run using: $options->{input}.\n";
    my @inputs = split(/$options->{delimiter}/, $options->{input});
    my $num_inputs = scalar(@inputs);
    my @jobs;
    ## 'separate' mode will submit a separate job for every sequence in the input fasta file.
    if ($options->{mode} eq 'separate') {
        @jobs = $class->Bio::Adventure::Structure::ProteinFold_JSON_Separate(
            options => $options, paths => $paths);
    }
    elsif ($options->{mode} eq 'together') {
        @jobs = $class->Bio::Adventure::Structure::ProteinFold_JSON_Together(
            options => $options, paths => $paths);
    }
    elsif ($options->{mode} eq 'pairwise' && $num_inputs == 1) {
        @jobs = $class->Bio::Adventure::Structure::ProteinFold_JSON_Pairwise_OneInput(
            options => $options, paths => $paths);
    }
    elsif ($options->{mode} eq 'pairwise' && $num_inputs == 2) {
        @jobs = $class->Bio::Adventure::Structure::ProteinFold_JSON_Pairwise_TwoInput(
            options => $options, first => $inputs[0], second => $inputs[1], paths => $paths);
    }
    else {
        die("I know not this option.\n");
    }
    return(\@jobs);
}

=head2 C<ProteinFold_JSON_Separate>

  Submit a separate job for every sequence in a fasta input file.

=cut
sub ProteinFold_JSON_Separate {
    my ($class, %args) = @_;
    my $options = $args{options};
    my $paths = $args{paths};
    my $molecule_type = $options->{libtype};
    my $in = Bio::Adventure::Get_FH(input => $options->{input});
    my $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $in);
    my @jobs = ();
  SEQ: while (my $seq = $seqio->next_seq) {
        my $seqid = $seq->id;
        my $sequence = $seq->seq;
        my $json_filename = qq"$paths->{output_dir}/${seqid}.json";
        my $json_fh = FileHandle->new(">${json_filename}");
        my $peptide = {
            $molecule_type => {
                id => ['A'],
                sequence => $sequence,
            }
        };
        my @sequences = ($peptide);
        my $datum = {
            name => $seqid,
            modelSeeds => [1],
            sequences =>  \@sequences,
            dialect => 'alphafold3',
            version => 1,
        };
        my $pretty = JSON->new->pretty->encode($datum);
        print $json_fh $pretty;
        $json_fh->close();
        my $final_full = abs_path($paths->{output_dir});
        ## Note the following which applies to our cluster: CUDA
        ## Capability 7.x GPUs For all CUDA Capability 7.x GPUs
        ## (e.g. V100) the environment variable XLA_FLAGS must be
        ## changed to include
        ## --xla_disable_hlo_passes=custom-kernel-fusion-rewriter. Disabling
        ## the Tritron GEMM kernels is not necessary as they are not
        ## supported for such GPUs.  I can find the capability of our
        ## nodes via the deviceQuery which resides in the examples/
        ## directory within the cuda module tree.  For the moment,
        ## though, it appears all of the nodes are version 7, thus I
        ## am going to set the custom-kernel-fusion-rewriter globally.
        my $jstring = qq!
export TMPDIR=${final_full}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=3.2
export XLA_FLAGS="\${XLA_FLAGS} --xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_enable_triton_gemm=false"
mkdir -p $paths->{output_dir}/jax
nvcc_location=\$(command -v nvcc)
if [[ \! -z "\${nvcc_location}" ]]; then
  cuda_location=\$(dirname \$(dirname \${nvcc_location}))
  query_location=query_location="\${cuda_location}/extras/demo_suite/deviceQuery"
  if [[ -x "\${query_location}" ]]; then
    \$query_location >> $paths->{output_dir}/queryDevice.stdout
  fi
fi
/usr/bin/time -v -o $paths->{stdout}.time -a \\
  run_alphafold.py \\
    --json_path ${json_filename} \\
    --model_dir \$ALPHA_HOME/models \\
    --output_dir $paths->{output_dir} \\
    --jax_compilation_cache_dir $paths->{output_dir}/jax \\
    --flash_attention_implementation=xla \\
    1>$paths->{stdout} 2>&1
!;
        my $job = $class->Submit(
            jdepends => $options->{jdepends},
            jname => qq"$options->{jname}_${seqid}",
            jstring => $jstring,
            jprefix => $options->{jprefix},
            jmem => $options->{jmem},
            jwalltime => $options->{jwalltime},
            language => 'bash',
            stderr => $paths->{stderr},
            stdout => $paths->{stdout},
            output => $paths->{output_dir},);
        push(@jobs, $job);
    }          ## End iterating over every sequence in the input file.
    return(@jobs);
}

=head2 C<ProteinFold_JSON_Separate>

  Submit a single job comprised of every sequence in the input file(s).

=cut
sub ProteinFold_JSON_Together {
    my ($class, %args) = @_;
    my $options = $args{options};
    my $paths = $args{paths};
    my $molecule_type = $options->{libtype};
    my $in = Bio::Adventure::Get_FH(input => $options->{input});
    my $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $in);
    my @jobs = ();
    my $run_id = basename($options->{input});
    my $json_filename = qq"$paths->{output_dir}/${run_id}.json";
    my $json_fh = FileHandle->new(">${json_filename}");
    my $datum = {
        name => $run_id,
        modelSeeds => [1],
        sequences =>  [],
        dialect => 'alphafold3',
        version => 1,
    };
    my @input_sequences = ();
    my $chain_id = 'A';
  SEQ: while (my $seq = $seqio->next_seq) {
        my $seqid = $seq->id;
        my $sequence = $seq->seq;
        my $peptide = {
            $molecule_type => {
                id => [$chain_id],
                sequence => $sequence,
            }
        };
        $chain_id++;
        push(@input_sequences, $peptide);
    }                       ## End iterating over the input sequences.
    if (defined($options->{input_rna})) {
        my $in_rna = Bio::Adventure::Get_FH(input => $options->{input_rna});
        my $seqio_rna = Bio::SeqIO->new(-format => 'fasta', -fh => $in_rna);
      SEQ: while (my $seq_rna = $seqio_rna->next_seq) {
            my $seqid = $seq_rna->id;
            my $sequence = $seq_rna->seq;
            my $rna = {
                rna => {
                    id => [$chain_id],
                    sequence => $sequence,
                }
            };
            $chain_id++;
            push(@input_sequences, $rna);
        }                       ## End iterating over rna sequences
    }
    if (defined($options->{input_dna})) {
        my $in_dna = Bio::Adventure::Get_FH(input => $options->{input_dna});
        my $seqio_dna = Bio::SeqIO->new(-format => 'fasta', -fh => $in_dna);
      SEQ: while (my $seq_dna = $seqio_dna->next_seq) {
            my $seqid = $seq_dna->id;
            my $sequence = $seq_dna->seq;
            my $dna = {
                dna => {
                    id => [$chain_id],
                    sequence => $sequence,
                }
            };
            $chain_id++;
            push(@input_sequences, $dna);
        }                       ## End iterating over dna sequences
    }
    $datum->{sequences} = \@input_sequences;
    my $pretty = JSON->new->pretty->encode($datum);
    print $json_fh $pretty;
    my $final_full = abs_path($paths->{output_dir});
    my $jstring = qq!
export TMPDIR=${final_full}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=3.2
export XLA_FLAGS="\${XLA_FLAGS} --xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_enable_triton_gemm=false"
mkdir -p $paths->{output_dir}/jax
nvcc_location=\$(command -v nvcc)
if [[ \! -z "\${nvcc_location}" ]]; then
  cuda_location=\$(dirname \$(dirname \${nvcc_location}))
  query_location=query_location="\${cuda_location}/extras/demo_suite/deviceQuery"
  if [[ -x "\$query_location" ]]; then
    \$query_location >> $paths->{output_dir}/queryDevice.stdout
  fi
fi
/usr/bin/time -v -o $paths->{stdout}.time -a \\
  run_alphafold.py \\
    --json_path ${json_filename} \\
    --model_dir \$ALPHA_HOME/models \\
    --output_dir $paths->{output_dir} \\
    --jax_compilation_cache_dir $paths->{output_dir}/jax \\
    --flash_attention_implementation=xla \\
    1>$paths->{stdout} 2>&1
!;
    my $job = $class->Submit(
        jdepends => $options->{jdepends},
        jname => qq"$options->{jname}_${run_id}",
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => $options->{jmem},
        jwalltime => $options->{jwalltime},
        language => 'bash',
        stderr => $paths->{stderr},
        stdout => $paths->{stdout},
        output => $paths->{output_dir},);
    push(@jobs, $job);
    $json_fh->close();
    return(@jobs);
}

=head2 C<ProteinFold_JSON_Separate>

  Submit a single job for every pair of sequences in the input.

=cut
sub ProteinFold_JSON_Pairwise_OneInput {
    my ($class, %args) = @_;
    my $options = $args{options};
    my $paths = $options->{paths};
    my $molecule_type = $options->{libtype};
    my $in = Bio::Adventure::Get_FH(input => $options->{input});
    my $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $in);
    my @seqs = ();
  SEQ: while (my $seq = $seqio->next_seq) {
        push(@seqs, $seq);
    }
    my @jobs = ();
  FST: for my $c (0 .. ($#seqs - 1)) {
        my $first = $seqs[$c];
        next FST unless (defined($first));
        my $first_id = $first->id;
        my $first_sequence = $first->seq;
      SCD: for my $d (1 .. $#seqs) {
            my $second = $seqs[$d];
            next SCD unless (defined($second));
            my $second_id = $second->id;
            my $second_sequence = $second->seq;
            my $id_string = qq"${first_id}_${second_id}";
            my $json_filename = qq"$paths->{output_dir}/${id_string}.json";
            my $json_fh = FileHandle->new(">${json_filename}");
            my $first_peptide = {
                $molecule_type => {
                    id => ['A'],
                    sequence => $first_sequence,
                },
            };
            my $second_peptide = {
                $molecule_type => {
                    id => ['B'],
                    sequence => $second_sequence,
                },
            };
            my @sequences = ($first_peptide, $second_peptide);
            my $datum = {
                name => $id_string,
                modelSeeds => [1],
                sequences =>  \@sequences,
                dialect => 'alphafold3',
                version => 1,
            };
            my $pretty = JSON->new->pretty->encode($datum);
            print $json_fh $pretty;
            $json_fh->close();
            my $final_full = abs_path($paths->{output_dir});
            my $jstring = qq!
export TMPDIR=${final_full}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=3.2
export XLA_FLAGS="\${XLA_FLAGS} --xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_enable_triton_gemm=false"
mkdir -p $paths->{output_dir}/jax
nvcc_location=\$(command -v nvcc)
if [[ \! -z "\${nvcc_location}" ]]; then
  cuda_location=\$(dirname \$(dirname \${nvcc_location}))
  query_location=query_location="\${cuda_location}/extras/demo_suite/deviceQuery"
  if [[ -x "\$query_location" ]]; then
    \$query_location >> $paths->{output_dir}/queryDevice.stdout
  fi
fi
/usr/bin/time -v -o $paths->{stdout}.time -a \\
  run_alphafold.py \\
    --json_path ${json_filename} \\
    --model_dir \$ALPHA_HOME/models \\
    --output_dir $paths->{output_dir} \\
    --jax_compilation_cache_dir $paths->{output_dir}/jax \\
    --flash_attention_implementation=xla \\
    1>$paths->{stdout} 2>&1
!;
            my $job = $class->Submit(
                jdepends => $options->{jdepends},
                jname => qq"$options->{jname}_${id_string}",
                jstring => $jstring,
                jprefix => $options->{jprefix},
                jmem => $options->{jmem},
                jwalltime => $options->{jwalltime},
                language => 'bash',
                stderr => $paths->{stderr},
                stdout => $paths->{stdout},
                output => $paths->{output_dir},);
            push(@jobs, $job);
        }
    }
    return(@jobs);
}

sub ProteinFold_JSON_Pairwise_TwoSeq {
    my ($class, %args) = @_;
    my $options = $args{options};
    my $molecule_type = $options->{libtype};
    $molecule_type = 'protein' unless(defined($molecule_type));
    my $first = $args{first};
    my $second = $args{second};
    my $paths = $args{paths};
    return(undef) unless (defined($first) and defined($second));
    my $first_id = $first->id;
    my $first_sequence = $first->seq;
    my $second_id = $second->id;
    my $second_sequence = $second->seq;
    $first_id =~ s/\.\d{1,2}$//;
    $first_id =~ s/:.*$//;
    $second_id =~ s/\.\d{1,2}$//;
    $second_id =~ s/:.*$//;
    my $id_string = qq"${first_id}_${second_id}";
    my $final_dir = qq"$paths->{output_dir}/${id_string}";
    make_path($final_dir) unless (-d $final_dir);
    my $final_full = abs_path($final_dir);
    my $json_filename = qq"${final_dir}/${id_string}.json";
    my $json_fh = FileHandle->new(">${json_filename}");
    my $first_peptide = {
        $molecule_type => {
            id => ['A'],
            sequence => $first_sequence,
        }
    };
    my $second_peptide = {
        $molecule_type => {
            id => ['B'],
            sequence => $second_sequence,
        }
    };
    my @sequences = ($first_peptide, $second_peptide);
    my $datum = {
        name => $id_string,
        modelSeeds => [1],
        sequences =>  \@sequences,
        dialect => 'alphafold3',
        version => 1,
    };
    my $pretty = JSON->new->pretty->encode($datum);
    print $json_fh $pretty;
    $json_fh->close();
    my $xla_flag = '';
    my $needs_xla = 1;
    if ($needs_xla) {
        $xla_flag = '--flash_attention_implementation xla';
    }
    my $jstring = qq!
export TMPDIR=${final_dir}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=3.2
export XLA_FLAGS="\${XLA_FLAGS} --xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_enable_triton_gemm=false"
mkdir -p $paths->{output_dir}/jax
nvcc_location=\$(command -v nvcc)
if [[ \! -z "\${nvcc_location}" ]]; then
  cuda_location=\$(dirname \$(dirname \${nvcc_location}))
  query_location=query_location="\${cuda_location}/extras/demo_suite/deviceQuery"
  if [[ -x "\$query_location" ]]; then
    \$query_location >> $paths->{output_dir}/queryDevice.stdout
  fi
fi
/usr/bin/time -v -o $paths->{stdout}.time -a \\
  run_alphafold.py ${xla_flag} \\
    --json_path ${json_filename} \\
    --model_dir \$ALPHA_HOME/models \\
    --output_dir ${final_dir} \\
    --jax_compilation_cache_dir $paths->{output_dir}/jax \\
    --flash_attention_implementation=xla \\
    1>${final_dir}/paths->{stdout} 2>&1
echo "${first_id},${second_id}" >> $paths->{output_dir}/finished.txt
!;
    my $job = $class->Submit(
        jdepends => $options->{jdepends},
        jname => qq"proteinfold_twoseq_${id_string}",
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => $options->{jmem},
        jcpu => $options->{jcpu},
        jgpu => $options->{jgpu},
        jwalltime => $options->{jwalltime},
        language => 'bash',
        stderr => $paths->{stderr},
        stdout => $paths->{stdout},
        output_dir => $paths->{output_dir},
        output => $json_filename,);
    return($job);
}

sub ProteinFold_JSON_Pairwise_TwoInput {
    my ($class, %args) = @_;
    my $options = $args{options};
    my $molecule_type = $options->{libtype};
    my $first = $args{first};
    my $second = $args{second};
    my $paths = $args{paths};
    my $first_in = Bio::Adventure::Get_FH(input => $first);
    my $first_seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $first_in);
    my $second_in = Bio::Adventure::Get_FH(input => $second);
    my $second_seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $second_in);
    my @jobs = ();
    my @first_seqs = ();
    my @second_seqs = ();
    while (my $first_seq = $first_seqio->next_seq) {
        push(@first_seqs, $first_seq) if (defined($first_seq));
    }
    while (my $second_seq = $second_seqio->next_seq) {
        push(@second_seqs, $second_seq) if (defined($second_seq));
    }
  FIRST: for my $first_seq (@first_seqs) {
      SECOND: for my $second_seq (@second_seqs) {
            my $first_id = $first_seq->id;
            my $first_sequence = $first_seq->seq;
            my $second_id = $second_seq->id;
            my $second_sequence = $second_seq->seq;
            my $id_string = qq"${first_id}_${second_id}";
            my $json_filename = qq"$paths->{output_dir}/${id_string}.json";
            my $json_fh = FileHandle->new(">${json_filename}");
            my $first_peptide = {
                $molecule_type => {
                    id => ['A'],
                    sequence => $first_sequence,
                }
            };
            my $second_peptide = {
                $molecule_type => {
                    id => ['B'],
                    sequence => $second_sequence,
                }
            };
            my @sequences = ($first_peptide, $second_peptide);
            my $datum = {
                name => $id_string,
                modelSeeds => [1],
                sequences =>  \@sequences,
                dialect => 'alphafold3',
                version => 1,
            };
            my $pretty = JSON->new->pretty->encode($datum);
            print $json_fh $pretty;
            $json_fh->close();
            my $final_dir = abs_path($paths->{output_dir});
            my $jstring = qq!
export TMPDIR=${final_dir}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=3.2
export XLA_FLAGS="\${XLA_FLAGS} --xla_disable_hlo_passes=custom-kernel-fusion-rewriter --xla_gpu_enable_triton_gemm=false"
mkdir -p $paths->{output_dir}/jax
nvcc_location=\$(command -v nvcc)
if [[ \! -z "\${nvcc_location}" ]]; then
  cuda_location=\$(dirname \$(dirname \${nvcc_location}))
  query_location=query_location="\${cuda_location}/extras/demo_suite/deviceQuery"
  if [[ -x "\$query_location" ]]; then
    \$query_location >> $paths->{output_dir}/queryDevice.stdout
  fi
fi
/usr/bin/time -v -o $paths->{stdout}.time -a \\
  run_alphafold.py \\
    --json_path ${json_filename} \\
    --model_dir \$ALPHA_HOME/models \\
    --output_dir $paths->{output_dir} \\
    --jax_compilation_cache_dir $paths->{output_dir}/jax \\
    --flash_attention_implementation=xla \\
    1>$paths->{stdout} 2>&1
!;
            my $job = $class->Submit(
                jdepends => $options->{jdepends},
                jname => qq"$options->{jname}_${id_string}",
                jstring => $jstring,
                jprefix => $options->{jprefix},
                jmem => $options->{jmem},
                jcpu => $options->{jcpu},
                jgpu => $options->{jgpu},
                jwalltime => $options->{jwalltime},
                language => 'bash',
                stderr => $paths->{stderr},
                stdout => $paths->{stdout},
                output => $paths->{output_dir},);
            push(@jobs, $job);
        }
    }
    $first_in->close();
    $second_in->close();
    return(@jobs);
}

=head2 C<RNAFold_Windows>

 Run RNAFold on a rolling window across a sequence.
 10.1093/nar/gkn188

=cut
sub RNAFold_Windows {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jcpu => 1,
        jprefix => 80,
        jname => 'vienna',
        length => 201,
        required => ['input'],
        step => 3,);
    my $output_name = basename($options->{input}, ('.gbk', '.fsa', '.fasta',));
    my $paths = $class->Bio::Adventure::Config::Get_Paths(output_name => $output_name);
    my $comment = '## Iterate over a sequence with RNAfold.';
    my $jstring = qq?
use Bio::Adventure::Structure;
\$h->Bio::Adventure::Structure::RNAFold_Windows_Worker(
  input => '$options->{input}',
  length => '$options->{length}',
  step => '$options->{step}',
  output => '$paths->{output}',
  output_dir => '$paths->{output_dir}',
  jprefix => '$options->{jprefix}',
  stdout => '$paths->{stdout}',
  stderr => '$paths->{stderr}',);
?;
    my $folder = $class->Submit(
        input => $options->{input},
        output => $paths->{output},
        jname => 'vienna',
        jprefix => $options->{jprefix},
        length => $options->{length},
        step => $options->{step},
        jstring => $jstring,
        comment => $comment,
        output_dir => $paths->{output_dir},
        language => 'perl',
        stdout => $paths->{stdout},
        stderr => $paths->{stderr},);
    return($folder);
}

=head2 C<RNAFold_Windows_Worker>

 Does the actual work of rolling across a sequence and running rnafold.

=cut
sub RNAFold_Windows_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        length => 201,
        step => 3,
        jname => 'vienna',
        jprefix => 80,
        required => ['input', 'output'],
        modules => ['vienna']);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $input_paths = $class->Get_Path_Info($options->{output});
    ## Put the data here!  First key is location, second is structure/mfe.
    my $output = $options->{output};
    my $out_dir = dirname($output);
    my $log = FileHandle->new(">${out_dir}/rnafold.log");
    print $log "Starting rolling window fold of $options->{input} with length $options->{length} and step $options->{step}.\n";
    my $out_txt = basename($output, ('.xz'));
    my $out_txt_path = qq"${out_dir}/${out_txt}";
    make_path($out_dir);
    my $txt_writer = FileHandle->new(">${out_txt_path}");
    print $txt_writer "contig\tstart\tend\tA\tU\tG\tC\tGC\tpaired\tbp_percent\tmfe\tmfe_bp\tmfe_gc\tstructure\n";

    my $length = $options->{length};
    my $step = $options->{step};
    my $in = Bio::Adventure::Get_FH(input => $options->{input});
    my $seqio;
    if ($options->{input} =~ /\.fasta|\.fsa/) {
        $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $in);
    }
    elsif ($options->{input} =~ /\.gb|\.gbff|\.gbk|\.gbf/) {
        $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
    }
    else {
        die("I do not understand the format for $options->{input}");
    }
    my $results = {};
    my $reader = gensym();
    my $writer = gensym();
    my $pid = open2($reader, $writer, 'RNAfold -t 4 --noPS');
    my $output_line;
  CONTIGS: while (my $seq = $seqio->next_seq) {
        ## Set up the first sequence and the stepper
        my $seq_length = $seq->length;
        my $circular = $seq;
        my $id = $seq->id;
        my $make_circular = $circular->is_circular(1);
        ## Make the sequence length divisible by 3 to simplify logic later
        if (($seq_length % 3) == 1) {
            $circular = $circular->trunc(1, ($seq_length - 1));
        }
        elsif (($seq_length % 3) == 2) {
            $circular = $circular->trunc(1, ($seq_length - 2));
        }
        my $post_length = $circular->length;
        ## Start the iterator a little bit before the beginning of the sequence.
        my $start = (1 - $length) + $step;
        my $end = $start + $length;
        my $continue = 1;
        my %nt_counts = (A => 0, U => 0, G => 0, C => 0);
        my $gc_content = 0;
        my $bp_count = 0;
        my $bp_percent = 0;
        my $key = qq"${id}_${start}_${end}";
      STEP: while ($continue) {
            if ($end > $post_length) {
                $start = $start - $post_length;
                $end = $end - $post_length;
            }

            my $sub_sequence = '';
            if ($start <= 0 && $end <= 0) {
                ## This should not happen, but we can handle it easily enough
                print "Somehow got a zero on start or end, this should not happen\n";
            }
            elsif ($start <= 0) {
                my $pre_start = $post_length + $start;
                my $pre_end = $post_length;
                my $remaining = ($length + $start);
                my $post_start = 1;
                my $post_end = $post_start + $remaining;
                my $pre = $circular->subseq($pre_start, $pre_end);
                my $post = $circular->subseq($post_start, $post_end);
                $sub_sequence = $pre . $post;
            }
            elsif ($end <= 0) {
                ## This really shouldn't happen.
                print "end is less than zero, that is weird.\n";
            }
            else {
                $sub_sequence = $circular->subseq($start, $end);
            }

            ## Make sure to send a return character to RNAfold
            $sub_sequence .= "\n";
            ## Send the subsequence of interest to RNAfold.
            print $writer $sub_sequence;
            ## And read its first line of output.
            $output_line = <$reader>;
            ## If the line is comprised entirely of word characters, then it is printing the sequence.
            if ($output_line =~ /^\w+$/) {
                my $seq_line = $output_line;
                ## Count up the nucleotides.
                $nt_counts{A} = $seq_line =~ tr/A//;
                $nt_counts{U} = $seq_line =~ tr/U//;
                $nt_counts{G} = $seq_line =~ tr/G//;
                $nt_counts{C} = $seq_line =~ tr/C//;
                $gc_content = sprintf("%0.4f", ($nt_counts{G} + $nt_counts{C}) / $options->{length});
                ## Now skip down to the next line of output from RNAfold.
                $output_line = <$reader>;
            }

            ## Pick out the basepairs and calculated MFE value.
            my ($structure, $mfe) = split(/\s+/, $output_line);
            ## Count up the base pairs
            my $structure_string = $structure;
            $bp_count = $structure_string =~ tr/\.//;
            $bp_percent = sprintf("%0.4f", $bp_count / $options->{length});

            $mfe =~ s/\(|\)//g;
            my $normalized_mfe_bp = sprintf("%0.4f", $mfe / $bp_count);
            my $normalized_mfe_gc = sprintf("%0.4f", $mfe / ($nt_counts{G} + $nt_counts{C}));
            my $txt_string = qq"${id}\t${start}\t${end}\t$nt_counts{A}\t$nt_counts{U}\t$nt_counts{G}\t$nt_counts{C}\t$gc_content\t$bp_count\t$bp_percent\t${mfe}\t${normalized_mfe_bp}\t${normalized_mfe_gc}\t${structure}\n";
            print $txt_writer $txt_string;
            ## print $txt_string;
            ## Make a record of the contig/start/end so we can tell when we are finished.
            $results->{$key} = 1;

            ## Done with the step, set new start/end/key and check to see if we are done.
            $start = $start + $step;
            $end = $end + $step;
            $key = qq"${id}_${start}_${end}";
            if (defined($results->{$key})) {
                last STEP;
            }
        }
    }
    $writer->close();
    my $compressed = qx"xz -9e -f ${out_txt_path}";
    return($results);
}

1;
