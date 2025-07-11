package Bio::Adventure::Align_Fasta;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::SearchIO::blast;
use Bio::SearchIO::fasta;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;
use Cwd 'abs_path';
use File::Basename qw"dirname basename";
use File::Path qw"make_path remove_tree";
use File::Which qw"which";
use POSIX qw"ceil";

=head1 NAME

 Bio::Adventure::Align_Fasta - Perform Sequence alignments with the fasta suite of tools.

=head1 SYNOPSIS

 All the functions which live here work with the fasta36 suite of
 programs.  They write the scripts, invoke them on the cluster,
 collect the results, and (optionally) parse them into simplified
 tables.

=head1 METHODS

=head2 C<Make_Fasta_Job>

 Write a job for the cluster to run a fasta36 job.

 This handles the creation of the files and directories required when
 splitting up a sequence into a bunch of pieces for a split-fasta search.

=over

=item C<Arguments>

 split(0): Split up the fasta jobs into multiple pieces?
 fasta_tool('ggsearch36'): Use this fasta36 tool.
 output_type(undef): Specify a particular output type, ideally one of
  the parseable formats.
 jdepends(): Put this in a dependency chain?
 jmem(8): Expected memory usage.
 modules('fasta'): Load this environment module.

=back

=cut
sub Make_Fasta_Job {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        align_jobs => 1,
        cluster => 'slurm',
        fasta_tool => 'ggsearch36',
        fasta_format => '9',
        jdepends => '',
        jmem => 8,
        jprefix => '90',
        modules => ['fasta', 'cyoa'],
        output_type => undef,
        split => 0,
        type => 'protein',);
    my $dep = $options->{jdepends};
    my $split = $options->{split};
    my $output_type = $options->{output_type};
    my $fasta_args = $class->Passthrough_Args(arbitrary => $options->{fasta_args});
    my $library;
    my $array_end = $options->{array_start} + $args{align_jobs};
    my $array_string = qq"$options->{array_start}-${array_end}";
    my $jstring = '';
    my $type_string = '';
    my $array_id_string = 'single_job';
    if ($options->{cluster} eq 'slurm') {
        $array_id_string = '${SLURM_ARRAY_TASK_ID}';
    } elsif ($options->{cluster} eq 'torque') {
        $array_id_string = '${PBS_ARRAYID}';
    }

    my $job_basename = $class->Get_Job_Name();
    if (defined($output_type)) {
        $type_string = "-m ${output_type}";
    }
    my $output = '';
    my $stdout = qq"$options->{workdir}/split_align.stdout";
    my $stderr = qq"$options->{workdir}/split_align.stderr";
    my $query_flag = '-n';
    if ($options->{type} eq 'protein') {
        $query_flag = '-p';
    } elsif ($options->{type} eq 'rna') {
        $query_flag = '-U';
    }

    ## Important Note:  fasta36's command line parsing fails on path names > 128 or 256 characters.
    ## Thus my usual '$options->{workdir}' will cause this to fail in many instances
    ## because it introduces many characters to the pathnames.
    if ($options->{cluster}) {
        $output = qq"$options->{workdir}/split/*/*.out";
        make_path(qq"$options->{workdir}/split") unless (-d qq"$options->{workdir}/split");
        $jstring = qq!
export IN="$options->{workdir}/split/${array_id_string}/in.fasta"
echo "This job is: ${array_id_string}"
if [[ -r "\${IN}" ]]; then
  $options->{fasta_tool} ${query_flag} -m $options->{fasta_format} \\
   $options->{fasta_args} ${type_string} -T $options->{jcpu} \\
   \${IN} \\
   $options->{library} \\
   1>$options->{workdir}/split/${array_id_string}/${array_id_string}.out \\
   2>>${stderr}
else
  echo "The input file does not exist: \${IN}"
  exit 0
fi
!;
    } else {
        $output = qq"$options->{workdir}/$options->{fasta_tool}.out";
        $jstring = qq!
$options->{fasta_tool} -m $options->{fasta_format} \\
  ${query_flag} $options->{fasta_args} ${type_string} -T $options->{jcpu} \\
  $options->{input} \\
  $options->{library} \\
  1>${output} \\
  2>>${stderr}
!;
    }
    my $comment = qq!## Running $options->{align_jobs} fasta job(s).!;

    my $fasta_jobs = $class->Submit(
        array_string => $array_string,
        cluster => $options->{cluster},
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => $job_basename,
        jstring => $jstring,
        language => 'bash',
        jprefix => qq"$options->{jprefix}_1",
        modules => $options->{modules},
        output => $output,
        stdout => $stdout,
        stderr => $stderr,);
    return($fasta_jobs);
}

=head2 C<Parse_Fasta>

 Parse a fasta36 tool result file and print a summary table.

 Given the output from one of the fasta36 programs: ggsearch36,
 fasta36, etc.  This parses it and prints a simplified table of the
 results.

=over

=item C<Arguments>

 input(required): File containing a completed fasta search output.
 best_only(0): Print only the best hit for each search sequence.
 evalue(0.0001): Filter for only hits with a better evalue than this.

=item C<Invocation>

> cyoa --task align --method parsefasta --input fasta_output.txt.xz

=back

=cut
sub Parse_Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        best_only => 0,
        evalue => 0.0001,);
    my $best_only = $options->{best_only};
    my $evalue = $options->{evalue};
    my $input = $options->{input};
    my $output = $options->{output};
    unless ($output) {
        $output = $input;
        $output =~ s/\..+//g;
        $output .= "_parsed.txt";
    }
    my $out = FileHandle->new(">${output}");
    ## This works, other attempts I've used for zipped input fail.
    my $f = $class->Get_FH(input => $input);
    my $searchio = Bio::SearchIO->new(-format => 'fasta',
                                      -fh => $f,
                                      -best => ${best_only},
                                      -signif => ${evalue});
    print $out "QueryName\tQueryLength\tLibID\tLibHitLength\tScore\tE\tIdentity\tHitMatches\tQueryStart\tQueryEnd\tLibStart\tLibEnd\tLStrand\n";
    my $results = 0;
    while (my $result = $searchio->next_result()) {
      $results++;
      my $hit_count = 0;
    HIT: while (my $hit = $result->next_hit) {
        $hit_count++;
        my $query_name = $result->query_name();
        my $query_length = $result->query_length();
        my $accession = $hit->accession();
        my $library_id = $hit->name();
        my $length = $hit->length();
        my $score = $hit->raw_score();
        my $sig = $hit->significance();
        my $ident = $hit->frac_identical();
        my $lstrand = $hit->strand('hit');
        my $hit_len;
        my $hit_matches;
        my $first_qstart = -1;
        my $first_qend = -1;
        my $first_lstart = -1;
        my $first_lend = -1;
        my $hsp_count = 0;
      HSP: while (my $hsp = $hit->next_hsp) {
          $hsp_count++;
          if ($hsp_count == 1) {
              $first_qstart = $hsp->start('query');
              $first_qend = $hsp->end('query');
              $first_lstart = $hsp->start('subject');
              $first_lend = $hsp->end('subject');
          }
          $hit_len = $hsp->length('total');
          my @matches = $hsp->matches(-seq => 'hit');
          $hit_matches = $matches[1];
      }
        print $out "${query_name}\t${query_length}\t${library_id}\t${length}\t${score}\t${sig}\t${ident}\t${hit_matches}\t${first_qstart}\t${first_qend}\t${first_lstart}\t${first_lend}\t${lstrand}\n";
        ## process the Bio::Search::Hit::HitI object
        ## while(my $hsp = $hit->next_hsp) {
        ## process the Bio::Search::HSP::HSPI object
        ## }
    }   ## End of each hit of a result.
  }       ## End of each result.
    $f->close();
    $out->close();
    return($results);
}

=head2 C<Parse_Fasta_Global>

 A separate parser for global:global searches.

 This was originally written in order to write summary tables of
 alignments of two closely related species (like L.major vs. T.cruzi)
 in order to categorize the 1:1 orthologs vs. gene families.

 Parse_Global makes a hard assumption about the structure of the hits
 in the output of the alignment program.  They should be in a
 global/global search, thus we assume 1 hit / 1 result.

 In addition, this will print summaries of hits in a few different
 files depending on how many hits were observed for each query
 sequence: singletons, doubles, triples, few (3-10), and many (10+).

=over

=item C<Arguments>

 best_only(0): Parse out only the best hit?
 fasta_format(0): Format of the search input
 evalue(0.001): Quality of hit filter
 many_cutoff(10): How many hits define 'many' hits?
 min_score(0): Provide a minimum score required for a hit?
 check_all_hits(0): I am not sure what this does.
 min_percent(0): Add a percentage identity filter

=back

=cut
sub Parse_Fasta_Global {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        best_only => 0,
        fasta_format => 'fasta',
        evalue => 0.001,
        many_cutoff => 10,
        min_score => undef,
        check_all_hits => 0,
        min_percent => undef,
        output => 'parsed.txt',
        output_counts => 'parsed_counts.txt',
        output_singles => 'parsed_singles.txt',
        output_doubles => 'parsed_doubles.txt',
        output_few => 'parsed_few.txt',
        output_many => 'parsed_many.txt',
        output_zero => 'parsed_zero.txt',
        output_all => 'parsed_all.txt',
        jmem => 8);
    my $output = $options->{input};
    my $outdir = dirname($output);
    $output = basename($output, ('.gz', '.xz'));
    $output = basename($output, ('.txt'));

    my $in = $class->Get_FH(input => $options->{input});
    my $searchio = Bio::SearchIO->new(
        -format => $options->{fasta_format},
        -fh => $in,
        -best_hit_only => $options->{best_only},
        -max_significance => $options->{evalue},
        -check_all_hits => $options->{check_all_hits},
        -min_score => $options->{min_score},);

    my $counts = FileHandle->new(qq">$options->{output_counts}");
    my $parsed = FileHandle->new(qq">$options->{output}");
    my $singles = FileHandle->new(qq">$options->{output_singles}");
    my $doubles = FileHandle->new(qq">$options->{output_doubles}");
    my $few = FileHandle->new(qq">$options->{output_few}");
    my $many = FileHandle->new(qq">$options->{output_many}");
    my $zero = FileHandle->new(qq">$options->{output_zero}");
    my $all = FileHandle->new(qq">$options->{output_all}");

    my $seq_count = 0;
    print $parsed "Query Name\tQuery length\tHit ID\tHit Length\tScore\tE\tIdentity\tHit length\tHit Matches\n";
    while (my $result = $searchio->next_result()) {
        print "Looking at result!\n";
        $seq_count++;
        my $query_name = $result->query_name();
        my $entry = qq"${query_name}\t";
        my $hit_count = 0;
      HITLOOP: while(my $hit = $result->next_hit) {
          my $query_length = $result->query_length();
          my $accession = $hit->accession();
          my $acc2 = $hit->name();
          my $length = $hit->length();
          my $score = $hit->raw_score();
          my $sig = $hit->significance();
          my $ident = $hit->frac_identical();
          my $hit_len;
          my $hit_matches;
          if (defined($args{min_percent})) {
              next HITLOOP if ($ident < $args{min_percent});
          }
          ## If we pass the min_percet threshold, increment the number of hits.
          $hit_count++;
          while (my $hsp = $hit->next_hsp) {
              $hit_len = $hsp->length('total');
              my @matches = $hsp->matches(-seq => 'hit');
              $hit_matches = $matches[1];
          }
          if ($hit_count == 1) {
              $entry .= "${acc2}\t${ident}\t${sig}";
          } else {
              $entry .= "${acc2}:${ident}:${sig} ";
          }
          print $parsed "${query_name}\t${query_length}\t${acc2}\t${length}\t${score}\t${sig}\t${ident}\t${hit_len}\t${hit_matches}\n";
      } ## End iterating over every hit for a sequence.
        $entry .= "\n";
        print $all $entry;
        if ($hit_count == 1) {
            print $singles $entry;
        } elsif ($hit_count == 2) {
            print $doubles $entry;
        } elsif ($hit_count >= 3 and $hit_count <= $options->{many_cutoff}) {
            print $few $entry;
        } elsif ($hit_count > $options->{many_cutoff}) {
            print $many $entry;
        } else {
            print $zero $entry;
        }
        print $counts "${query_name}\t${hit_count}\n";
    } ## End iterating over every result in searchio
    $counts->close();
    $parsed->close();
    $singles->close();
    $doubles->close();
    $few->close();
    $many->close();
    $zero->close();
    $all->close();
    $in->close();
    return($seq_count);
}

sub Parse_Fasta_Mismatches {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jprefix => 50,
        modules => ['fasta', 'cyoa'],
        required => ['input', 'library'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules},
                                       exe => ['fasta36']);
    my $in_base = basename($options->{input}, ('.fsa', '.faa', '.ffn', '.gbk'));
    my $lib_base = basename($options->{library}, ('.fsa', '.faa', '.ffn', '.gbk'));
    my $output_base = qq"${in_base}_vs_${lib_base}";
    my $outdir = qq"outputs/$options->{jprefix}pairwisedelta";
    make_path("${outdir}") unless(-d ${outdir});
    my $gg_output = qq"${outdir}/${output_base}.txt";
    my $gg_error = qq"${outdir}/${output_base}.err";
    my $variant_output = qq"${outdir}/${output_base}_mismatches.txt";
    my $variant_numbers = qq"${outdir}/${output_base}_mismatch_nums.txt";
    my $runner = qq"ggsearch36 -b 1 -d 1 $options->{library} $options->{input} 1>${gg_output} 2>${gg_error}";
    my $handle = IO::Handle->new;
    open($handle, "$runner |");
    while (my $line = <$handle>) {
        print "$line\n";
    }
    close($handle);
    print "Comparing: $options->{library} and $options->{input}\n";
    my $out = FileHandle->new(">${variant_output}");
    my $numbers = FileHandle->new(">${variant_numbers}");
    ## Write a header for the output tsv.
    print $out "RefID\tQueryID\tposition\ttype\treference\thit\n";
    my $results = 0;
    my $indices = {};
    my $searchio = Bio::SearchIO->new(-format => 'fasta', -file => $gg_output);
  SEARCHLOOP: while (my $result = $searchio->next_result()) {
      $results++;
      ## There is only ever 1 hit.
      my $hit = $result->next_hit;
      my $query_name = $result->query_name();
      my $query_length = $result->query_length();
      next SEARCHLOOP unless (defined($hit));
      my $accession = $hit->accession();
      my $library_id = $hit->name();
      if (!defined($accession)) {
          $accession = $library_id;
      }
      if (!defined($accession)) {
          $accession = 'unknown';
      }
      my $length = $hit->length();
      ##my $score = $hit->raw_score();
      ##my $sig = $hit->significance();
      ##my $ident = $hit->frac_identical();
      my $hstrand = $hit->strand('hit');
      my $qstrand = $hit->strand('query');
      my $index = 'unknown';
      if ($hstrand ne $qstrand) {
          ## I found a couple of reads which are messed up and reverse complementing
          ## I am not completely certain this will make them go away, but I think it will.
          print "SKIPPING: ${query_name} hit strand: $hstrand query strand: $qstrand\n";
          next SEARCHLOOP;
      }

      my $hit_len;
      my $hit_matches;
      my $first_qstart = -1;
      my $first_qend = -1;
      my $first_lstart = -1;
      my $first_lend = -1;
      my $hsp_count = 0;
      my $hits = {};

      my $hsp = $hit->next_hsp;
      ## IMPORTANT POINT OF CONFUSION:
      ## The hit is the template sequence!
      ## The query is my provided sequence.
      ##      my @template_mismatches = $hsp->seq_inds('hit', 'mismatch');
      ##      my @product_mismatches = $hsp->seq_inds('query', 'mismatch');
      ## Ok, so the various documentation strings for HSPI and friends are a bit confusing.
      ## If I want to get the non-gap non-identical indices, I need to explicitly ask for two things:
      ## First, the 'mismatches' get the positions which are completely different and denoted by ' '
      ## Then 'conserved-not-identical' provides the similar but not-identical not-gap positions
      ## which are denoted by '.'.
      my @template_mismatches = $hsp->seq_inds('hit', 'mismatch');
      my @product_mismatches = $hsp->seq_inds('query', 'mismatch');
      my @template_not_identical = $hsp->seq_inds('hit', 'conserved-not-identical');
      my @product_not_identical = $hsp->seq_inds('query', 'conserved-not-identical');
      my @template_both = sort {$a <=> $b} (@template_mismatches, @template_not_identical);
      my @product_both = sort {$a <=> $b} (@product_mismatches, @product_not_identical);

      my @template_gaps = $hsp->seq_inds('hit', 'gap');
      my @product_gaps = $hsp->seq_inds('query', 'gap');

      my $num_mis = scalar(@template_both);
      my $tem_gap = scalar(@template_gaps);
      my $pro_gap = scalar(@product_gaps);
      my $num_muts = $num_mis + $tem_gap + $pro_gap;
      print $numbers "${query_name}\t${num_muts}\n";
      my $t_seq = $hsp->hit_string;
      my @tseq = split(//, $t_seq);
      my @template = @tseq;
      my $p_seq = $hsp->query_string;
      my @pseq = split(//, $p_seq);
      my $t_raw = $t_seq;
      $t_raw =~ s/\-//g;
      $t_raw =~ s/^\s+//g;
      my @traw = split(//, $t_raw);
      my $p_raw = $p_seq;
      $p_raw =~ s/\-//g;
      $p_raw =~ s/^\s+//g;
      my @praw = split(//, $p_raw);
      ## I am increasingly confident that we cannot trust indels at the end of alignments.
      ## No matter how I mess with the gap score/penalties, I get dumb alignments at the end.
      ## The tests which look at the end of the aligned sequence attempt to address this problem
      ## but I think a more aggressive approach will be required.  I think I will have to check
      ## the position of the indel and just arbitrarily drop it if it is >= 200.
      my $print_position = 0;
      if (scalar(@template_gaps) > 0) {
          for my $i (0 .. $#template_gaps) {
              my $inc = $i + 1;
              my $total = scalar(@template_gaps);
              my $t_pos = $template_gaps[$i];
              my $position = $t_pos;
              $print_position = $position + 1;
              ## I think the next two lines are a particularly nasty hack, and
              ## I do not understand why they seem to be necessary to get
              ## correct results.
              my $mutated_to = $praw[$position + $i];
              if (!defined($mutated_to)) {
                  $mutated_to = 'undef';
              }
              my $mutated_from = $tseq[$position + $i];
              if (!defined($mutated_from)) {
                  $mutated_from = 'undef';
              }
              my $template_length = scalar(@traw);
              $hits->{$position}->{type} = 'ins';
              $hits->{$position}->{to} = $mutated_to;
              print $out qq"${query_name}\t${library_id}\t${print_position}\tInsertion\t \t${mutated_to}\n";
          }
      }
      if (scalar(@product_gaps) > 0) {
          for my $i (0 .. $#product_gaps) {
              my $p_pos = $product_gaps[$i];
              my $inc = $i + 1;
              my $total = scalar(@product_gaps);
              my $position = $p_pos;
              $print_position = $position + 1;
              my $mutated_to = $pseq[$position + $i];
              my $mutated_from = $tseq[$position + $i];
              ## Now check that the alignment did not end.
              my $product_length = scalar(@pseq);
              $hits->{$position}->{type} = 'del';
              $hits->{$position}->{from} = $mutated_from;
              print $out qq"${query_name}\t${library_id}\t${print_position}\tDeletion\t${mutated_from}\t \n";
          }
      }
      if (scalar(@template_both) > 0) {
          ## I expect the length of query_mismatches and hit_mismatches to be identical.
          for my $i (0 .. $#template_both) {
              ## Get the index position of the mismatch.
              ## This can be confusing because of indels.
              ## I want them 0 indexed, but fasta36 gives them as 1 indexed, so subtract 1.
              my $t_pos = $template_both[$i] - 1;
              my $p_pos = $product_both[$i] - 1;
              ## $mutated_from appears to be the correct one.
              my $mutated_from = $template[$t_pos];
              my $inc = $i + 1;
              my $total = scalar(@template_both);
              ## The product nucleotide is correct, except if there are indels
              ## before the mismatches...  In that case, shenanigans occur,
              ## and I absolutely cannot seem to figure out why.

              my $mutation_type = 'Mismatch';
              ##if ($template_mismatches[0] eq $template_both[$i]) {
              ##    shift @template_mismatches;
              ##    $mutation_type = 'Mismatch';
              ##} else {
              ##    shift @template_not_identical;
              ##    $mutation_type = 'Conserved';
              ##}
              my $mutated_to = $praw[$p_pos];
              my $position = $t_pos;
              $print_position = $position + 1;
              $hits->{$position}->{type} = 'mis';
              $hits->{$position}->{from} = $mutated_from;
              $hits->{$position}->{to} = $mutated_to;
              print $out qq"${query_name}\t${library_id}\t${print_position}\t${mutation_type}\t${mutated_from}\t${mutated_to}\n";
          }
      }
  } ## End of each hit of a result.
    $out->close();
    $numbers->close();
    my $unloaded = $class->Module_Reset(env => $loaded);
    return($indices);
}

=head2 C<Split_Align_Fasta>

 Split apart a set of query sequences into $args{align_jobs} pieces and align them all separately.

 This is the main callable function when one wants to perform a bunch
 of fasta36 searches in parallel.

=over

=item C<Arguments>

 input(required): File containing sequences to search _for_.
 library(required): File containing sequences to search _from_.
 align_jobs(40): How many jobs should be spawned?
 align_parse(0): Start up a parsing job upon completion?
 best_only(0): Passed to the parser, only print the best hit?

=item C<Invocation>

> cyoa --task fastasplit --input query.fasta --library library.fasta --best_only 1

=back

=cut
sub Split_Align_Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        align_jobs => 40,
        align_parse => 0,
        best_only => 0,
        fasta_tool => 'ggsearch36',
        interactive => 0,
        jmem => 8,
        jprefix => 90,
        num_dirs => 0,
        required => ['input', 'library'],
        type => 'dna',
        modules => ['fasta', 'cyoa']);
    my $loaded = $class->Module_Loader(modules => $options->{modules},
                                       exe => ['fasta36']);
    my $job_basename = $class->Get_Job_Name();
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $jname = qq"${job_basename}_fasta";
    my $comment = qq'## Running a split alignment with $options->{fasta_tool}.';
    my $jstring = qq!use Bio::Adventure::Align;
use Bio::Adventure::Align_Fasta;
my \$result = \$h->Bio::Adventure::Align_Fasta::Split_Align_Fasta_Worker(
  align_jobs => '$options->{align_jobs}',
  align_parse => '$options->{align_parse}',
  best_only => '$options->{best_only}',
  fasta_tool => '$options->{fasta_tool}',
  input => '$options->{input}',
  interactive => '$options->{interactive}',
  jname => '${jname}',
  jprefix => '$options->{jprefix}',
  library => '$options->{library}',
  num_dirs => '$options->{num_dirs}',
  output => '$paths->{output}',
  output_counts => '$paths->{output_counts}',
  output_singles => '$paths->{output_singles}',
  output_doubles => '$paths->{output_doubles}',
  output_few => '$paths->{output_few}',
  output_many => '$paths->{output_many}',
  output_zero => '$paths->{output_zero}',
  output_all => '$paths->{output_all}',
  type => '$options->{type}',);
!;
    my $align = $class->Submit(
        align_jobs => $options->{align_jobs},
        align_parse => $options->{align_parse},
        best_only => $options->{best_only},
        comment => $comment,
        fasta_tool => $options->{fasta_tool},
        input => $options->{input},
        interactive => $options->{interactive},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        library => $options->{library},
        num_dirs => $options->{num_first},
        output => $paths->{output},
        output_counts => $paths->{output_counts},
        output_singles => $paths->{output_singles},
        output_doubles => $paths->{output_doubles},
        output_few => $paths->{output_few},
        output_many => $paths->{output_many},
        output_zero => $paths->{output_zero},
        output_all => $paths->{output_all},
        type => $options->{type},
        stdout => $paths->{stdout},
        stderr => $paths->{stderr});
    return($align);
}

sub Split_Align_Fasta_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        align_jobs => 40,
        align_parse => 0,
        best_only => 0,
        fasta_tool => 'ggsearch36',
        interactive => 0,
        jmem => 8,
        jprefix => 90,
        num_dirs => 0,
        required => ['input', 'library'],
        type => 'dna',
        modules => ['fasta', 'cyoa']);
    my $loaded = $class->Module_Loader(modules => $options->{modules},
                                       exe => ['fasta36']);
    my $search_type = 'local_fasta';
    my $search_function = 'Parse_Fasta';
    if ($options->{fasta_tool} =~ /^gg/) {
        $search_type = 'global_fasta';
        $search_function = 'Parse_Fasta_Global';

    }
    my $job_basename = $class->Get_Job_Name();
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $lib = basename($options->{library}, $class->{suffixes});
    my $que = basename($options->{input}, $class->{suffixes});
    my $outdir = $paths->{output_dir};
    my $output = $paths->{output};
    my $jname = qq"${job_basename}_fasta",
    my $split_info = $class->Bio::Adventure::Align::Get_Split(
        %args, align_jobs => $options->{align_jobs});
    my $actual = $class->Bio::Adventure::Align::Make_Directories(
        workdir => $outdir,
        num_per_split => $split_info->{num_per_split},);
    print "Actually used ${actual} directories to write files.\n";
    my $alignment = $class->Bio::Adventure::Align_Fasta::Make_Fasta_Job(
        align_jobs => $actual,
        cluster => $options->{cluster},
        fasta_tool => $options->{fasta_tool},
        input => abs_path($options->{input}),
        library => abs_path($options->{library}),
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => qq"$options->{jprefix}_1",
        type => $options->{type},
        workdir => $outdir,);
    my $concat_job = $class->Bio::Adventure::Align::Concatenate_Searches(
        cluster => $options->{cluster},
        input => $alignment->{output},
        jdepends => $alignment->{job_id},
        jname => qq"${job_basename}_concatfasta",
        jprefix => qq"$options->{jprefix}_2",
        workdir => $outdir,);
    ## Make sure that the alignment job gets the final output
    $alignment->{output} = $concat_job->{output};
    my $parse_input = $alignment->{output};
    my $comment_string = qq!## I don't know if this will work.!;
    my $jstring = qq?
use Bio::Adventure;
use Bio::Adventure::Align;
use Bio::Adventure::Align_Fasta;
my \$result = \$h->Bio::Adventure::Align_Fasta::${search_function}(
  input => '$parse_input',
  best_only => '$options->{best_only}',
  output => '$paths->{output}',
  output_counts => '$paths->{output_counts}',
  output_singles => '$paths->{output_singles}',
  output_doubles => '$paths->{output_doubles}',
  output_few => '$paths->{output_few}',
  output_many => '$paths->{output_many}',
  output_zero => '$paths->{output_zero}',
  output_all => '$paths->{output_all}',
  search_type => '${search_type}',);
?;
    my $parse_job = $class->Submit(
        best_only => $options->{best_only},
        comment => $comment_string,
        cluster => $options->{cluster},
        input => $parse_input,
        jdepends => $concat_job->{job_id},
        jmem => $options->{jmem},
        jname => qq"${job_basename}_parsefasta",
        jprefix => qq"$options->{jprefix}_3",
        jstring => $jstring,
        output => $paths->{output},
        output_counts => $paths->{output_counts},
        output_singles => $paths->{output_singles},
        output_doubles => $paths->{output_doubles},
        output_few => $paths->{output_few},
        output_many => $paths->{output_many},
        output_zero => $paths->{output_zero},
        output_all => $paths->{output_all},
        stdout => $paths->{stdout},
        stderr => $paths->{stderr},
        language => 'perl',);
    $parse_job->{align} = $alignment;
    $parse_job->{concat} = $concat_job;
    return($parse_job);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<fasta36>

=cut

1;
