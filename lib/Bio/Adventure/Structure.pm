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

sub AlphaFold {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jcpu => 1,
        jgpu => 1,
        jprefix => 80,
        jname => 'alphafold',
        libtype => 'protein',
        mode => 'separate',
        required => ['input']);
    my $output_name = basename($options->{input}, ('.gbk', '.fsa', '.fasta',));
    my $paths = $class->Bio::Adventure::Config::Get_Paths(output_name => $output_name);
    my $comment = '## Iterate over a sequence with Alphafold.';
    my $jstring = qq?
use Bio::Adventure::Structure;
\$h->Bio::Adventure::Structure::AlphaFold_Worker(
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

=head2 C<AlphaFold_Worker>

 Does the actual work of submitting queries to alphafold.

=cut
sub AlphaFold_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jname => 'alphafold',
        jprefix => 80,
        mode => 'separate',
        required => ['input']);
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
        @jobs = $class->Bio::Adventure::Structure::AlphaFold_JSON_Separate(
            options => $options, paths => $paths);
    } elsif ($options->{mode} eq 'together') {
        @jobs = $class->Bio::Adventure::Structure::AlphaFold_JSON_Together(
            options => $options, paths => $paths);
    } elsif ($options->{mode} eq 'pairwise' && $num_inputs == 1) {
        @jobs = $class->Bio::Adventure::Structure::AlphaFold_JSON_Pairwise_OneInput(
            options => $options, paths => $paths);
    } elsif ($options->{mode} eq 'pairwise' && $num_inputs == 2) {
        @jobs = $class->Bio::Adventure::Structure::AlphaFold_JSON_Pairwise_TwoInput(
            options => $options, first => $inputs[0], second => $inputs[1], paths => $paths);
    } else {
        die("I know not this option.\n");
    }
    return(\@jobs);
}

=head2 C<AlphaFold_JSON_Separate>

  Submit a separate job for every sequence in a fasta input file.

=cut
sub AlphaFold_JSON_Separate {
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
        my $jstring = qq!
run_alphafold.py \\
  --json_path ${json_filename} \\
  --model_dir \$ALPHA_HOME/models \\
  --output_dir $paths->{output_dir} \\
  2>$paths->{stderr} \\
  1>$paths->{stdout}
!;
        my $job = $class->Submit(
            jdepends => $options->{jdepends},
            jname => qq"$options->{jname}_${seqid}",
            jstring => $jstring,
            jprefix => $options->{jprefix},
            jmem => $options->{jmem},
            jwalltime => '08:00:00',
            language => 'bash',
            stderr => $paths->{stderr},
            stdout => $paths->{stdout},
            output => $paths->{output_dir},);
        push(@jobs, $job);
    } ## End iterating over every sequence in the input file.
    return(@jobs);
}

=head2 C<AlphaFold_JSON_Separate>

  Submit a single job comprised of every sequence in the input file(s).

=cut
sub AlphaFold_JSON_Together {
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
    } ## End iterating over the input sequences.
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
        } ## End iterating over rna sequences
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
        } ## End iterating over dna sequences
    }
    $datum->{sequences} = \@input_sequences;
    my $pretty = JSON->new->pretty->encode($datum);
    print $json_fh $pretty;
    my $jstring = qq!
run_alphafold.py \\
  --json_path ${json_filename} \\
  --model_dir \$ALPHA_HOME/models \\
  --output_dir $paths->{output_dir} \\
  2>$paths->{stderr} \\
  1>$paths->{stdout}
!;
    my $job = $class->Submit(
        jdepends => $options->{jdepends},
        jname => qq"$options->{jname}_${run_id}",
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => $options->{jmem},
        jwalltime => '08:00:00',
        language => 'bash',
        stderr => $paths->{stderr},
        stdout => $paths->{stdout},
        output => $paths->{output_dir},);
    push(@jobs, $job);
    $json_fh->close();
    return(@jobs);
}

=head2 C<AlphaFold_JSON_Separate>

  Submit a single job for every pair of sequences in the input.

=cut
sub AlphaFold_JSON_Pairwise_OneInput {
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
    for my $c (0 .. ($#seqs - 1)) {
        my $first = $seqs[$c];
        my $first_id = $first->id;
        my $first_sequence = $first->seq;
        for my $d (1 .. $#seqs) {
            my $second = $seqs[$d];
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
            my $jstring = qq!
run_alphafold.py \\
  --json_path ${json_filename} \\
  --model_dir \$ALPHA_HOME/models \\
  --output_dir $paths->{output_dir} \\
  2>$paths->{stderr} \\
  1>$paths->{stdout}
!;
            my $job = $class->Submit(
                jdepends => $options->{jdepends},
                jname => qq"$options->{jname}_${id_string}",
                jstring => $jstring,
                jprefix => $options->{jprefix},
                jmem => $options->{jmem},
                jwalltime => '08:00:00',
                language => 'bash',
                stderr => $paths->{stderr},
                stdout => $paths->{stdout},
                output => $paths->{output_dir},);
            push(@jobs, $job);
        }
    }
    return(@jobs);
}

sub AlphaFold_JSON_Pairwise_TwoInput {
    my ($class, %args) = @_;
    my $options = $args{options};
    my $molecule_type = $options->{libtype};
    my $first = $args{first};
    my $second = $args{second};
    my $paths = $args{paths};
    print "TESTME: First: $first second: $second\n";
    my $first_in = Bio::Adventure::Get_FH(input => $first);
    my $first_seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $first_in);
    my $second_in = Bio::Adventure::Get_FH(input => $second);
    my $second_seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $second_in);
    my @jobs = ();
    my @first_seqs = ();
    my @second_seqs = ();
    while (my $first_seq = $first_seqio->next_seq) {
        push(@first_seqs, $first_seq);
    }
    while (my $second_seq = $second_seqio->next_seq) {
        push(@second_seqs, $second_seq);
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
            my $jstring = qq!
run_alphafold.py \\
  --json_path ${json_filename} \\
  --model_dir \$ALPHA_HOME/models \\
  --output_dir $paths->{output_dir} \\
  2>$paths->{stderr} \\
  1>$paths->{stdout}
!;
            my $job = $class->Submit(
                jdepends => $options->{jdepends},
                jname => qq"$options->{jname}_${id_string}",
                jstring => $jstring,
                jprefix => $options->{jprefix},
                jmem => $options->{jmem},
                jwalltime => '08:00:00',
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
    } elsif ($options->{input} =~ /\.gb|\.gbff|\.gbk|\.gbf/) {
        $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
    } else {
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
      } elsif (($seq_length % 3) == 2) {
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
        } elsif ($start <= 0) {
            my $pre_start = $post_length + $start;
            my $pre_end = $post_length;
            my $remaining = ($length + $start);
            my $post_start = 1;
            my $post_end = $post_start + $remaining;
            my $pre = $circular->subseq($pre_start, $pre_end);
            my $post = $circular->subseq($post_start, $post_end);
            $sub_sequence = $pre . $post;
        }   elsif ($end <= 0) {
            ## This really shouldn't happen.
            print "end is less than zero, that is weird.\n";
        } else  {
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
