package Bio::Adventure::Count;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::Adventure::Config;
use Bio::DB::Sam;
use Bio::FeatureIO;
use Bio::Tools::GFF;
use Cwd;
use File::Basename;
use File::Path qw"make_path";
use File::Which qw"which";
use List::Util qw"sum";
use String::Approx qw"amatch";
use Tree::DAG_Node;

my @ISA=qw"Tree::DAG_Node";

=head1 NAME

 Bio::Adventure::Count - Perform Sequence alignments counting with HTSeq

=head1 SYNOPSIS

 These functions handle the counting of reads, primarily via htseq.

=head1 METHODS

=cut

sub Bedtools_Coverage {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jprefix => '60',
        required => ['input', 'species'],);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}bedtools_coverage_$options->{species}";
    my $genome = qq"$options->{libdir}/$options->{libtype}/$options->{species}.fasta";
    my $stderr = qq"${output_dir}/bedtools_coverage.stderr";
    my $stdout = qq"${output_dir}/bedtools_coverage_$options->{species}.bed";
    my $comment = '## Calculate coverage of an alignment vs. genome.';
    my $jstring = qq!mkdir -p ${output_dir}
bedtools genomecov -bga -ibam $options->{input} 1>${stdout} 2>${stderr}
!;
    my $bedtools = $class->Submit(
        comment => $comment,
        jcpu => 1,
        jdepends => $options->{jdepends},
        jname => qq"bedcov_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 24,
        output => $stdout,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($bedtools);
}

sub Feature_Counts {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species',],
        fractional => 1,
        minlength => 50,
        maxlength => 800,
        paired => 1,
        jcpu => 1,
        jmem => 12,
        jname => '',
        jprefix => '',);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $output_name = basename($options->{input}, ('.bam', '.sam'));
    my $output_dir = dirname($options->{input});
    my $output = qq"${output_dir}/${output_name}";
    my $sort = $options->{sort};
    my $qual = $options->{qual};
    if (!defined($qual)) {
        print "Setting minimum quality to 0.\n";
        $qual = '0';
    }
    my $stranded = $options->{stranded};
    if (defined($stranded)) {
        if ($stranded eq 'forward') {
            $stranded = 1;
        } elsif ($stranded eq 'reverse') {
            $stranded = 2;
        } else {
            $stranded = 0;
        }
    } else {
        $stranded = 0;
    }

    my $fc_invocation = qq"counted=0
gff_type=$options->{gff_type}
gff_tag=$options->{gff_tag}
max_length=$options->{maxlength}
max_reads=$options->{max_reads}
min_length=$options->{minlength}
min_quality=${qual}
mode=$options->{mode}
secondary=$options->{secondary}
stranded=${stranded}

{
  /usr/bin/time -v -o ${output}.time -a featureCounts \\
    -T 1 ";
    if ($options->{introns}) {
        $fc_invocation .= qq" -J -G $paths->{fasta}";
    }
    if ($options->{fractional}) {
        $fc_invocation .= ' -M --fraction';
    }
    $fc_invocation .= qq" -s \${stranded}";

    $fc_invocation .= qq" -Q \${min_quality}";
    if ($options->{paired}) {
        $fc_invocation .= qq" -p --countReadPairs -B -P -d \${minlength} -D \${maxlength}";
    }
    $fc_invocation .= qq" \\
    -R CORE -t \${gff_type} -g \${gff_tag} \\
    -a $paths->{gff} \\
";
    $output .= qq"_s\${stranded}_\${gff_type}_\${gff_tag}_fcounts.csv";
    my $error = basename($output, ('.csv'));
    my $stderr = qq"${output_dir}/${error}.stderr";
    my $stdout = qq"${output_dir}/${error}.stdout";
    $fc_invocation .= qq!    -o ${output} \\
    $options->{input} \\
    2>>${stderr} \\
    1>>${stdout}
  counted=\$?
} || {
  echo 'featureCounts failed.' >> ${stderr}
}
if [[ "\${counted}" = 0 ]]; then
  xz -9e -f ${output}
  xz -9e -f $options->{input}.featureCounts;
fi
!;
    $output .= '.xz';
    my $comment = qq!## Counting the number of hits in $options->{input} for each feature found in $paths->{gff} via featureCounts.
## Is this stranded? ${stranded}.
!;
    my $fc_jobname = qq"fcount_$options->{mapper}_$options->{species}_s${stranded}_$options->{gff_type}_$options->{gff_tag}";
    my $fc = $class->Submit(
        comment => $comment,
        gff_type => $options->{gff_type},
        gff_tag => $options->{gff_tag},
        input => $options->{input},
        output => $output,
        stderr => $stderr,
        stdout => $stdout,
        postscript => $args{postscript},
        prescript => $args{prescript},
        jdepends => $options->{jdepends},
        jcpu => $options->{jcpu},
        jmem => $options->{jmem},
        jname => $fc_jobname,
        jprefix => $options->{jprefix},
        jstring => $fc_invocation,);
    return($fc);
}

sub Guess_Strand {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species', 'gff_type', 'gff_tag',],
        genes => 'all',
        jmem => 12,
        jprefix => '41',
        output => 'strand_counts.txt',);
    my $job_name = $class->Get_Job_Name();
    my $output = $options->{output};
    my $stdout = dirname($output) . qq"/stranded_test.stdout";
    my $comment = '## This submits a job to figure out the strandedness of a library.';
    my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Count;
my \$result = Bio::Adventure::Count::Guess_Strand_Worker(\$h,
  input => '$options->{input}',
  genes => '$options->{genes}',
  gff_tag => '$options->{gff_tag}',
  gff_type => '$options->{gff_type}',
  output => '${output}',
  species => '$options->{species}',);
!;
    my $guess = $class->Submit(
        comment => $comment,
        genes => $options->{genes},
        gff_tag => $options->{gff_tag},
        gff_type => $options->{gff_type},
        input => $options->{input},
        jcpu => 1,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => qq"guess_strand_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        output => $output,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        species => $options->{species},);
    return($guess);
}

sub Guess_Strand_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species', 'gff_type', 'gff_tag',],
        maximum => 4000,
        coverage => 0.2,
        output_log => '',
        output => '',);
    my $log = FileHandle->new(">$options->{output}.log");
    my $gff = qq"$options->{libpath}/$options->{libtype}/$options->{species}.gff";
    my $fasta = qq"$options->{libpath}/$options->{libtype}/$options->{species}.fasta";
    my $features = $class->Read_Genome_GFF(
        gff => $gff,
        gff_type => $options->{gff_type},
        gff_tag => $options->{gff_tag},);
    my $sam = Bio::DB::Sam->new(-bam => $options->{input},
                                -fasta => ${fasta},);
    my @targets = $sam->seq_ids;
    my $num = scalar(@targets);
    my $bam = Bio::DB::Bam->open($options->{input});
    my $header = $bam->header;
    my $target_count = $header->n_targets;
    my $target_names = $header->target_name;
    my $align_count = 0;
    my $million_aligns = 0;
    my $alignstats = qx"samtools idxstats $options->{input}";
    my @alignfun = split(/\n/, $alignstats);
    my @aligns = split(/\t/, $alignfun[0]);
    my @unaligns = split(/\t/, $alignfun[1]);
    my $number_reads = $aligns[2] + $unaligns[3];
    my $output_name = qq"$options->{input}.txt";
    my $summary_string = qq"There are ${number_reads} alignments in $options->{input}.
 There are $aligns[2] aligned reads and $unaligns[3] unaligned reads.\n";
    print $log $summary_string;
    print $summary_string;
    my $forward_hit = 0;
    my $reverse_hit = 0;
    my %result = ();
    my $read_num = 0;
    my $counters = {
        counted => 0,
        paired => 0,
        proper_pair => 0,
        this_unmapped => 0,
        mate_unmapped => 0,
        this_reverse => 0,
        mate_reverse => 0,
        first_pair => 0,
        second_pair => 0,
        not_primary => 0,
        failed_sequencer => 0,
        optical_duplicate => 0,
        supplemental_align => 0,
        pair_facing_in => 0,
        pair_facing_out => 0,
        straddle_feature => 0,
        inside_feature => 0,
        outside_feature => 0,
        num_mismatches => 0,
        mean_mismatches => 0,
        total_bases => 0,
        plus_stranded => 0,
        minus_stranded => 0,
        mean_scores => 0,
        mean_scored => 0,
        grand_mean_score => 0,
    };
    print "Going to stop after $options->{maximum} reads\n";
  BAMLOOP: while (my $align = $bam->read1) {
        $read_num++;
        $counters->{counted}++;
        if ($read_num >= $options->{maximum}) {
            last BAMLOOP;
        }
        ##if ($class->{debug}) {  ## Stop after a relatively small number of reads when debugging.
        ##    last BAMLOOP if ($align_count > 200);
        ##}
        my $seqid = $target_names->[$align->tid];
        ## my $start = $align->pos + 1;
        ## my $end = $align->calend;
        my $start = $align->pos;
        my $end = $align->calend - 1;
        my $strand = $align->strand;
        my $seq = $align->query->dna;
        my $qual = $align->qual;
        my $pairedp = $align->paired;
        my $unmappedp = $align->unmapped;
        my $mate_unmappedp = $align->munmapped;
        my $reversedp = $align->reversed;
        my $mate_reversedp = $align->mreversed;
        my $mate_id = $align->mtid;
        my $mate_start = $align->mpos;
        my $properp = $align->proper_pair;
        if ($unmappedp && $mate_unmappedp) {
            next BAMLOOP;
        }

        my @seq_array = split(//, $seq);
        my @tags = $align->get_all_tags;
        ## A reminder when playing with tag values and get_all_tags()
        ## These are coming from the 2nd sam column and are the result of a decimal->binary conversion:
        ## E.g. a paired, proper, first of pair, reverse in binary is: 1100101 -> 1+2+16+64 -> 83
        ## 0x1(1) or first binary character: paired or unpaired
        ## 0x2(2) or second binary character: proper or not proper pair
        ## 0x4(4): this read is unmapped
        ## 0x8(8): the mate read is unmapped
        ## 0x10(16): this read is reverse
        ## 0x20(32): the mate is reverse
        ## 0x40(64): this is the first of a pair
        ## 0x80(128): this is the second of a pair
        ## 0x100(256): this is not the primary read
        ## 0x200(512): this failed the sequencer checks
        ## 0x400(1024): this is an optical duplicate
        ## 0x800(2048): this is a supplemental alignment.
        ## These flags get the following constant names for get_tag_values():
        ## 0x0001 => 'PAIRED', 0x0002 => 'MAP_PAIR', 0x0004 => 'UNMAPPED', 0x0008 => 'M_UNMAPPED',
        ## 0x0010 => 'REVERSED', 0x0020 => 'M_REVERSED', 0x0040 => 'FIRST_MATE',
        ## 0x0080 => 'SECOND_MATE', 0x0100 => 'NOT_PRIMARY', 0x0200 => 'QC_FAILED',
        ## 0x0400 => 'DUPLICATE', 0x0800 => 'SUPPLEMENTARY',
        $counters->{paired}++ if ($align->get_tag_values('PAIRED'));
        $counters->{proper_pair}++ if ($align->get_tag_values('MAP_PAIR'));
        $counters->{this_unmapped}++ if ($align->get_tag_values('UNMAPPED'));
        $counters->{mate_unmapped}++ if ($align->get_tag_values('M_UNMAPPED'));
        $counters->{this_reversed}++ if ($align->get_tag_values('REVERSED'));
        $counters->{mate_reversed}++ if ($align->get_tag_values('M_REVERSED'));
        $counters->{first_pair}++ if ($align->get_tag_values('FIRST_MATE'));
        $counters->{second_pair}++ if ($align->get_tag_values('SECOND_MATE'));
        $counters->{not_primary}++ if ($align->get_tag_values('NOT_PRIMARY'));
        $counters->{failed_sequencer}++ if ($align->get_tag_values('QC_FAILED'));
        $counters->{optical_duplicate}++ if ($align->get_tag_values('DUPLICATE'));
        $counters->{supplemental_align}++ if ($align->get_tag_values('SUPPLEMENTARY'));
        $counters->{total_aligned_bases} += $align->query->length;
        if ($align->strand > 0) {
            $counters->{plus_stranded}++;
        }
        else {
            $counters->{minus_stranded}++;
        }
        my $score_ref = $align->qscore;
        my @scores = @{$score_ref};
        if (scalar(@scores) > 0) {
            sub this_mean {
                return(sum (@_)/@_);
            }
            $counters->{mean_scores} += this_mean(@scores);
            $counters->{mean_scored} += scalar(@scores);
        }
            else {
                print "No score ref\n";
            }

            my $cigar = $align->cigar_str;
            if ($cigar eq '') {
                ## print "This did not align.\n";
                next BAMLOOP;
            }
            $align_count++;
            ## Look for a feature at this position...
            #my %chr_features = %{$features->{$seqid}};
            #for my $feat_id (keys %chr_features) {
            #    my %feat = %{$chr_features{$feat_id}};
            #
            #    ##  ------------------->>>>>-------<<<<<--------------------------
            #    ##                  >>>>>
            #    if ($start <= $feat{start} && $end >= $feat{start} && $end < $feat{end}) {
            #        use Data::Dumper;
            #        #print Dumper %feat;
            #    }
            #}
        }   ## End reading each bam entry
        $counters->{grand_mean_score} = $counters->{mean_scores} / $counters->{mean_scored};
        # $counters->{mean_mismatches} = $counters->{num_mismatches} / $counters->{total_bases};
        # $counters->{proportion_unmapped} = (($counters->{this_unmapped} + $counters->{mate_unmapped}) / 2) / $counters->{couted};
        print $log "
Reads counted: $counters->{counted}
Paired reads: $counters->{paired}
Properly paired: $counters->{proper_pair}
This unmapped: $counters->{this_unmapped}
Mate unmapped: $counters->{mate_unmapped}
This reverse: $counters->{this_reverse}
Mate reverse: $counters->{mate_reverse}
First pair: $counters->{first_pair}
Second pair: $counters->{second_pair}
Not primary: $counters->{not_primary}
Failed sequencer: $counters->{failed_sequencer}
Optical duplicates: $counters->{optical_duplicate}
Supplemental alignments: $counters->{supplemental_align}
Facing in: $counters->{pair_facing_in}
Facing out: $counters->{pair_facing_out}
Straddling features: $counters->{straddle_feature}
Inside features: $counters->{inside_feature}
Outside features: $counters->{outside_feature}
Total mismatches: $counters->{num_mismatches}
Mean mismatches: $counters->{mean_mismatches}
Total nucleotides: $counters->{total_bases}
Plus stranded: $counters->{plus_stranded}
Minus stranded: $counters->{minus_stranded}
Grand score_mean: $counters->{grand_mean_score}
Proportion unmapped: $counters->{proportion_unmapped}
";
        $log->close();
    }

=head2 C<HT_Multi>

 Invoke htseq multiple times with options for counting different transcript types.

=over

=item C<Arguments>

 species(required): Defines the gff/fasta files to query in {libpath}.
 input(required): The sorted/indexed .bam alignment file to count.
 htseq_stranded(required): Boolean of the stranded state of the library.
 gff_type(''): Type of feature to count.  If left alone, this will count up
   each category of feature and choose the most represented.
   (deprecated in favor of gff_type)
 gff_type('gene): Ibid.  I just have not converted everything to use this.
 gff_tag('ID'): GFF tag used to identify the genes/transcripts/etc.  In
   metazoans this is usually 'ID', in bacteria it is often
   'locus_tag', in parasites 'gene_id'.  There are lots of choices.
 libtype('genome'): Used to differentiate between genomic counts vs. rRNA
   vs. contaminants.
 modules('htseq'): List of environment modules.
 paired('1'): Is this a paired library?

=back

=cut
sub HT_Multi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my %ro_opts = %{$options};
    my $species = $options->{species};
    my $htseq_input = $options->{input};
    my $stranded = $options->{stranded};
    if ($stranded eq '1') {
        $stranded = 'yes';
    }
    elsif ($stranded eq 'forward') {
        $stranded = 'yes';
    }
    elsif ($stranded eq '0') {
        $stranded = 'no';
    }

    my $gff_type = $options->{gff_type};
    my $gff_tag = $options->{gff_tag};
    my @jobs = ();
    my $script_suffix = '';
    if ($args{suffix}) {
        $script_suffix = $args{suffix};
    }
    my $jprefix = '';
    $jprefix = $args{jprefix} if ($args{jprefix});
    my @gff_types = ('antisense', 'exon', 'fiveputr', 'interCDS',
                     'linc', 'mi', 'misc', 'nmd', 'operons', 'pseudo',
                     'rintron', 'rrna', 'sn', 'sno', 'threeputr');
    my $htseq_runs = 0;
    ## Top level directory containing the input files.
    my $top_dir = basename(getcwd());
    ## The sam/bam input basename
    my $output_name = basename($htseq_input, ('.bam', '.sam'));

    foreach my $gff_type (@gff_types) {
        my $gff = $paths->{gff};
        $gff =~ s/\.gff$/_${gff_type}\.gff/g;
        my $gtf = $gff;
        $gtf =~ s/\.gff/\.gtf/g;
        my $htseq_jobname = qq"hts_${output_name}_$options->{species}_s${stranded}_${gff_type}_${gff_tag}";
        if (-r "$gff") {
            print "Found $gff, performing htseq with it.\n";
            my $ht = $class->Bio::Adventure::Count::HTSeq(
                htseq_gff => $gff,
                input => $htseq_input,
                gff_type => $gff_type,
                gff_tag => $options->{gff_tag},
                jdepends => $options->{jdepends},
                jname => $htseq_jobname,
                jprefix => $jprefix,
                postscript => $options->{postscript},
                prescript => $options->{prescript},
                suffix => $options->{suffix},);
            push(@jobs, $ht);
            $htseq_runs++;
        } elsif (-r "$gtf") {
            print "Found ${gtf}, performing htseq with it.\n";
            my $ht = $class->Bio::Adventure::Count::HTSeq(
                gff_type => $gff_type,
                htseq_gff => $gtf,
                gff_tag => $options->{gff_tag},
                input => $htseq_input,
                jdepends => $options->{jdepends},
                jname => $htseq_jobname,
                jprefix => $options->{jprefix},
                postscript => $options->{postscript},
                prescript => $options->{prescript},
                suffix => $options->{suffix},);
            push(@jobs, $ht);
            $htseq_runs++;
        }
    } ## End foreach type
    ## Also perform a whole genome count
    my $htall_jobname = qq"htall_${output_name}_$options->{species}_s${stranded}_$ro_opts{gff_type}_$ro_opts{gff_tag}";
    if (-r "$paths->{gff}") {
        print "Found $paths->{gff}, performing htseq_all with it.\n";
        my $ht = $class->Bio::Adventure::Count::HTSeq(
            htseq_gff => $paths->{gff},
            gff_tag => $ro_opts{gff_tag},
            gff_type => $ro_opts{gff_type},
            input => $htseq_input,
            jcpu => 1,
            jdepends => $options->{jdepends},
            jname => $htall_jobname,
            jprefix => $jprefix,
            postscript => $options->{postscript},
            prescript => $options->{prescript},
            suffix => $options->{suffix},);
        push(@jobs, $ht);
    } elsif (-r "$paths->{gtf}") {
        print "Found $paths->{gtf}, performing htseq_all with it.\n";
        my $ht = $class->Bio::Adventure::Count::HTSeq(
            htseq_gff => $paths->{gff},
            gff_type => 'none',
            input => $htseq_input,
            jcpu => 1,
            jdepends => $options->{jdepends},
            jname => $htall_jobname,
            jprefix => $jprefix,
            postscript => $args{postscript},
            prescript => $args{prescript},
            suffix => $args{suffix},);
        push(@jobs, $ht);
    } else {
        print "Did not find $paths->{gff} nor $paths->{gtf}, not running htseq_all.\n";
    }
    return(\@jobs);
}

=head2 C<HT_Types>

 Guess about most appropriate flags for htseq-count.

 B
 This function reads the first 100k lines of a gff file and use that
 to guess at the most likely type of feature when invoking htseq.

=over

=item C<Arguments>

 gff_type('gene'): When set, this will just count that type.
 gff_tag('ID'): Ditto, but the GFF tag for IDs.

=back

=cut
sub HT_Types {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,);
    my $my_type = $options->{gff_type};
    my $my_id = $options->{gff_tag};
    print "Examining $options->{annotation} for type: ${my_type} and tag: ${my_id}.\n";
    $my_type = '' unless($my_type);
    $my_type = '' unless($my_id);
    my $gff_out = {};
    my $max_lines = 100000;
    my $count = 0;
    my %found_types = ();
    my %found_canonical = ();
    my $reader = new FileHandle("<$options->{annotation}");
  LOOP: while(my $line = <$reader>) {
        chomp $line;
        next LOOP if ($line =~ /^\#/);
        $count++;
        ## chr1    unprocessed_pseudogene  transcript      3054233 3054733 .       +       .       gene_id "ENSMUSG00000090025"; transcript_id "ENSMUST00000160944"; gene_name "Gm16088"; gene_source "havana"; gene_biotype "pseudogene";transcript_name "Gm16088-001"; transcript_source "havana"; tag "cds_end_NF"; tag "cds_start_NF"; tag "mRNA_end_NF"; tag "mRNA_start_NF";

        my ($chromosome, $source, $type, $start, $end, $score, $strand, $frame, $attributes) = split(/\t/, $line);
        my @attribs = split(/;\s*/, $attributes);
        my (@pairs, @names, @values) = ();
        my $splitter = '\s+';
        if ($attribs[0] =~ /=/) {
            $splitter = '=';
        }
        for my $attrib (@attribs) {
            my ($name, $value) = split(/$splitter/, $attrib);
            push(@names, $name);
            push(@values, $value);
        }
        my $canonical = $names[0];
        if ($found_types{$type}) {
            $found_types{$type}++;
        } else {
            $found_types{$type} = 1;
        }
        if ($found_canonical{$canonical}) {
            $found_canonical{$canonical}++;
        } else {
            $found_canonical{$canonical} = 1;
        }
        if ($count > $max_lines) {
            last LOOP;
        }
    }                           ## End loop
    $reader->close();
    my $found_my_type = 0;
    my $max_type = '';
    my $max = 0;
    foreach my $type (keys %found_types) {
        if ($found_types{$type} > $max) {
            $max_type = $type;
            $max = $found_types{$type};
        }
        if ($found_types{$my_type}) {
            print "The specified type: ${my_type} is in the gff file, comprising $found_types{$my_type} of the first 40,000.\n";
            $found_my_type = 1;
        }
    }                           ## End the loop

    my $max_can = 0;
    my $max_canonical = 0;
    foreach my $can (keys %found_canonical) {
        if (!defined($found_canonical{$can})) {
            $found_canonical{$can} = 0;
        }
        if ($found_canonical{$can} > $max_canonical) {
            $max_can = $can;
            $max_canonical = $found_canonical{$can};
        }
    }                           ## End the loop
    my $returned_canonical = $max_can;

    my $returned_type = "";
    if ($found_my_type == 1) {
        $returned_type = $my_type;
    } else {
        print "Did not find your specified type.  Changing it to: ${max_type} which had ${max} entries.\n";
        $returned_type = $max_type;
    }
    my $ret = [$returned_type, $returned_canonical];
    return($ret);
}

=head2 C<HTSeq>

 Run htseq-count on a sorted, indexed bam file.
 10.1093/bioinformatics/btu638

 Invoke htseq-count.  This should be able to automagically pick up
 other types of countable features and send the htseq-count results to
 a separate count file.

=over

=item C<Arguments>

 input: Sorted/indexed bam file to count.
 species: Defines the gff/fasta files to read.
 stranded: Is this library stranded?
 htseq_args(--order=name --idattr=gene_id --minaqual=10 --type exon
  --stranded=yes --mode=union): Define arbitrary htseq arguments.
 gff_type(''): Redundant with gff_type, used to choose a specific
  type to count; when left blank, HT_Types() will make a guess.
 gff_type('gene'): Deprecated, Ibid.
 gff_tag('ID'): GFF tag used to identify the counted features.
 jname(''): Job name base for the cluster.
 jprefix(''): Prefix for the job name and output directory.
 libtype('genome'): Choose the library to count,
  genomic/rRNA/contaminants/etc.
 mapper('hisat2'): Specify the mapper used so that the output file
  will make it easier to tell where it came from.
 modules('htseq'): List of environment modules to load.
 paired('1'): Is this library paired?

=back

=cut
sub HTSeq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species', 'htseq_args',],
        jcpu => 1,
        jmem => 30,
        jname => '',
        jprefix => '',);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $sort = $options->{sort};
    my $stranded = $options->{stranded};
    if ($stranded eq '1') {
        $stranded = 'yes';
    }
    elsif ($stranded eq 'forward') {
        $stranded = 'yes';
    }
    elsif ($stranded eq '0') {
        $stranded = 'no';
    }
    print "Note, this is running with stranded: $options->{stranded}\n";
    my $mode = $options->{mode};
    my $secondary = $options->{secondary};
    my $supplementary = $options->{supplementary};

    my $aqual = $options->{qual};
    if (!defined($aqual)) {
        print "Setting minimum quality to 0.\n";
        $aqual = '0';
    }

    my $gff_tag = $options->{gff_tag};
    my $htseq_input = $options->{input};
    my $gff_type = 'all';
    if (defined($options->{gff_type}) && $options->{gff_type} ne '') {
        $gff_type = $options->{gff_type};
    }
    ## Top level directory containing the input files.
    my $top_dir = basename(getcwd());
    ## The sam/bam input basename
    ## Keep in mind this input bam file likely has the genome as part of its name
    my $output_name = basename($htseq_input, ('.bam', '.sam'));

    ## And directory containing it.
    my $output_dir = dirname($htseq_input);
    my $output = qq"${output_dir}/${output_name}";
    my $gff = $paths->{gff};
    $gff = $args{htseq_gff} if ($args{htseq_gff});
    my $gtf = $gff;
    $gtf =~ s/\.gff/\.gtf/g;
    my $htseq_args = "";

    ## Set the '-t FEATURETYPE --type' argument used by htseq-count
    ## This may be provided by a series of defaults in %HT_Multi::gff_types, overridden by an argument
    ## Or finally auto-detected by HT_Types().
    ## This is imperfect to say the nicest thing possible, I need to consider more appropriate ways of handling this.
    my $gff_type_arg = '';
    my $gff_tag_arg = '';
    my $stdout = '';
    my $stderr = '';
    my $annotation = $paths->{gtf};
    if (!-r "${gtf}") {
        $annotation = $gff;
    }
    my $variable_string = '';
    if (!defined($gff_type) or $gff_type eq '' or $gff_type eq 'auto' or
        !defined($gff_tag) or $gff_tag eq '' or $gff_tag eq 'auto') {
        my $gff_type_pair = $class->Bio::Adventure::Count::HT_Types(
            annotation => $annotation,
            type => $gff_type,);
        $gff_type = $gff_type_pair->[0];
        $gff_type_arg = qq" --type ${gff_type}";
        $gff_tag_arg = qq" --idattr ${gff_tag}";
    } elsif (ref($gff_type) eq 'ARRAY') {
        $gff_type_arg = qq" --type $gff_type->[0]";
        $gff_tag_arg = qq" --idattr ${gff_tag}";
    } elsif ($gff_type eq 'none') {
        $gff_type_arg = qq'';
        $gff_tag_arg = qq'';
    } else {
        $gff_type_arg = qq" --type \${gff_type}";
        $gff_tag_arg = qq" --idattr \${gff_tag}";
    }
    if ($gff_type_arg) {
        $variable_string = qq!gff_type=${gff_type}
gff_tag=${gff_tag}
stranded=${stranded}
mode=${mode}
secondary=${secondary}
supplementary=${supplementary}
sort=${sort}
max_reads=$options->{max_reads}
!;
    }
    $output .= qq"_r${sort}_s${stranded}_${gff_type}_${gff_tag}.csv";
    if (!-r "${gff}" and !-r "${gtf}") {
        die("Unable to read ${gff} nor ${gtf}, please fix this and try again.\n");
    }
    my $error = basename($output, ('.csv'));
    $stderr = qq"${output_dir}/${error}.stderr";
    $stdout = qq"${output_dir}/${error}.stdout";

    my $htseq_jobname = qq"hts_${top_dir}_$options->{mapper}_$options->{species}_s${stranded}_${gff_type}_${gff_tag}";
    my $htseq_invocation = qq!${variable_string}
/usr/bin/time -v -o ${output}.time -a htseq-count \\
  -q -f bam \\
  -s \${stranded} -a ${aqual} -r \${sort} \\
  ${gff_type_arg} ${gff_tag_arg} --mode \${mode} \\
  --secondary-alignments \${secondary} --supplementary-alignments \${supplementary} \\
  --max-reads-in-buffer \${max_reads} \\
  --counts_output ${output} \\!;
    my $jstring = qq!${htseq_invocation}
  ${htseq_input} \\
  ${annotation} \\
  2>${stderr} \\
  1>${stdout}
xz -f -9e ${output}
!;
    $output = qq"${output}.xz";
    my $comment = qq!## Counting the number of hits in ${htseq_input} for each feature found in ${annotation}
## Is this stranded? ${stranded}.  The defaults of htseq are:
## $options->{htseq_args}
!;
    my $htseq = $class->Submit(
        comment => $comment,
        gff_type => $options->{gff_type},
        gff_tag => $options->{gff_tag},
        input => $htseq_input,
        output => $output,
        stderr => $stderr,
        stdout => $stdout,
        postscript => $args{postscript},
        prescript => $args{prescript},
        jdepends => $options->{jdepends},
        jcpu => $options->{jcpu},
        jmem => $options->{jmem},
        jname => $htseq_jobname,
        jprefix => $options->{jprefix},
        jstring => $jstring,);
    return($htseq);
}

sub Insert_Size {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        jprefix => '50',);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $output_base = basename($options->{input}, ('.bam'));
    my $output_dir = $paths->{output_dir};
    my $output_file = qq"${output_dir}/${output_base}_insert_size.txt";
    my $output_pdf = qq"${output_dir}/${output_base}_insert_size.pdf";
    my $jstring = qq!
gatk CollectInsertSizeMetrics \\
  -I $options->{input} \\
  -O ${output_file} \\
  -H ${output_pdf} \\
  -M 0.5 2>$paths->{stderr} 1>$paths->{stdout}
!;

    my $comment_string = qq"## Use gatk to collect insert size statistics.";
    my $gatk = $class->Submit(
        comment => $comment_string,
        jdepends => $options->{jdepends},
        jcpu => 4,
        jmem => 12,
        jname => qq"insertsize_${output_base}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => '2:00:00',
        output => $output_file,
        stderr => $paths->{stderr},
        stdout => $paths->{stdout},);
    return($gatk);
}

=head2 C<Jellyfish>

 Run jellyfish with multiple values of K on a fast(a|q) file(s).
 10.1093/bioinformatics/btr011

 Use jellyfish to count up all of the kmers in a set of sequence(s).
 This should also send the wacky fasta format fo counted sequences to a
 tsv of kmers and numbers, which should be a rather more tractable
 format with which to play.

=over

=item C<Arguments>

 input(required): Fast(a|q) file(s) to index.
 length('9,11,13,15'): Text list of k values to give to jellyfish.
 jprefix(18): Prefix for the job name.
 modules('jellyfish'): Module list.

=back

=cut
sub Jellyfish {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        length => '9,11,13,15',
        jprefix => 18,);
    my @kmer_array = split(/[:;,]/, $options->{length});
    my $count = 0;
    my $ret;
    if (scalar(@kmer_array) > 1) {
        my $first = shift @kmer_array;
        my $job = $class->Bio::Adventure::Count::Jellyfish(
            %{$options},
            length => $first);
      KMERARR: for my $k (@kmer_array) {
            my $job = $class->Bio::Adventure::Count::Jellyfish(
                %{$options},
                length => $k);
            if ($count == 0) {
                $ret = $job;
            }
            else {
                $ret->{$k} = $job;
            }
            $count++;
        }
        return($ret);
    }
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Path_Info($options->{input});
    my $cwd_name = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}jellyfish_${cwd_name}";
    my $input_string = $class->Get_FD(input => $options->{input});
    if ($options->{input} =~ /$options->{delimiter}/) {
        ## Then multiple files were provided.
        my @file_list = split(/$options->{delimiter}/, $options->{input});
        $input_string = '';
        for my $f (@file_list) {
            my $fd = $class->Get_FD(input => $f);
            $input_string .= qq"${fd} ";
        }
    }

    my $jelly_base = qq"${output_dir}/${cwd_name}_$options->{length}";
    my $count_file = qq"${jelly_base}.count";
    my $info_file = qq"${jelly_base}.info";
    my $histogram_file = qq"${jelly_base}.hist";
    my $count_fasta = qq"${jelly_base}_by_count.fasta";
    my $matrix_file = qq"${jelly_base}_matrix.csv";
    my $comment = '## Invoke jellyfish on some sequence!';
    my $stdout = qq"${jelly_base}.stdout";
    my $stderr = qq"${jelly_base}.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
jellyfish count -m $options->{length} \\
  -o ${count_file} \\
  -s 50000 -t 4 \\
  ${input_string} \\
  2>${stderr} 1>${stdout}
jellyfish info ${count_file} > ${info_file} \\
  2>>${stderr}
jellyfish histo ${count_file} > ${histogram_file} \\
  2>>${stderr}
jellyfish dump ${count_file} > ${count_fasta} \\
  2>>${stderr}
!;

    my $jelly = $class->Submit(
        comment => $comment,
        count_fasta => $count_fasta,
        count_file => $count_file,
        histogram_file => $histogram_file,
        info_file => $info_file,
        jdepends => $options->{jdepends},
        jname => qq"jelly_${job_name}_$options->{length}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 12,
        jcpu => 4,
        output => $matrix_file,
        output_dir => $output_dir,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stderr => $stderr,
        stdout => $stdout,);
    $comment = qq"## This should create a matrix with rows as kmers and elements
## comprised of the number of occurrences.
";
    $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Phage;
my \$result = \$h->Bio::Adventure::Count::Jellyfish_Matrix(
  input => '${count_fasta}',
  jdepends => '$jelly->{job_id}',
  jname => 'jelly_matrix',
  output_dir => '${output_dir}',  ## Defines where the pdata goes.
  output => '${matrix_file}',);
!;
    my $matrix_job = $class->Submit(
        comment => $comment,
        input => $count_fasta,
        jdepends => $jelly->{job_id},
        jname => qq"jelly_matrix_$options->{length}",
        jprefix => qq"$options->{jprefix}_1",
        jstring => $jstring,
        jcpu => 1,
        language => 'perl',
        output_dir => $output_dir,
        output => $matrix_file,
        stdout => $stdout,);
    $jelly->{matrix_job} = $matrix_job;

    my $compress_files = qq"${count_file}:${info_file}:${histogram_file}:${count_fasta}:${matrix_file}";
    my $comp = $class->Bio::Adventure::Compress::Recompress(
        comment => '## Compress the jellyfish output files.',
        input => $compress_files,
        jdepends => $matrix_job->{job_id},
        jname => qq"xzjelly_$options->{length}",
        jprefix => qq"$options->{jprefix}_2",
        jcpu => 1,
        jmem => 4,);
    $jelly->{compression} = $comp;

    $jelly->{output} = qq"${matrix_file}.xz";
    $jelly->{count_file} = qq"${count_file}.xz";
    $jelly->{info_file} = qq"${info_file}.xz";
    $jelly->{histogram_file} = qq"${histogram_file}.xz";
    $jelly->{count_fasta} = qq"${count_fasta}.xz";
    return($jelly);
}

=head2 C<Jellyfish_Matrix>

 Convert the jellyfish fasta format to a matrix-compatible tsv.

 This function is responsible for actually converting the fasta output
 from jellyfish into tsv.  It is pretty quick and dirty.

=over

=item C<Arguments>

 input(required): The peculiar fasta file produced by jellyfish.
 output('fasta_matrix.csv'): Filename to which to convert the
  jellyfish output.
 jprefix(19): Prefix for the jobname/output directory.

=back

=cut
sub Jellyfish_Matrix {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        output => 'fasta_matrix.csv',
        jprefix => 19,);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Path_Info($options->{input});
    my $fc = $class->Get_FC(input => $options->{input});
    my $in = FileHandle->new("${fc} |");
    my $counter = 1;
    my $counts = {};
    my $nmer_count = 0;
    my $nmer_identity = '';
    while (my $line = <$in>) {
        chomp($line);
        if ($counter == 1) { ## Then it is the number line, also why in the flying hell did they do that
            $nmer_count = $line;
            $nmer_count =~ s/^>//g;
            $counter++;
        }
        elsif ($counter == 2) {
            $counter--;
            $nmer_identity = $line;
            $counts->{$nmer_identity} = $nmer_count;
        }
        else {
            die('Should not get here.');
        }
    }
    $in->close();

    my $out = FileHandle->new(">$options->{output}");
    foreach my $k (sort keys %{$counts}) {
        print $out "${k}\t$counts->{$k}\n";
    }
    $out->close();
    return($nmer_count);
}

=head2 C<Kraken>

 Use kraken2 to taxonomically classify reads.

 Kraken2 is a kmer-based read classifier.  It is quite fast once its
 databases are built.

=over

=item C<Arguments>

 input(required): Set of fastq reads.
 library('viral'): Kraken2 database to classify against.
 jprefix(11): Prefix for the job name and output directory.
 modules('kraken'): Environment module to load.

=item C<Invocation>

> cyoa --task annot --method kraken --input read1.fastq.xz:read2.fastq.xz --library standard

=back

=cut
sub Kraken {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',
        jcpu => 2,
        jmem => 64,
        jprefix => '11',
        clean => 1,);
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $job_name = $class->Get_Job_Name();
    my $input_directory = basename(cwd());
    print "Classifying reads against the $options->{library} kraken database.\n";
    my $output_dir = qq"outputs/$options->{jprefix}kraken_$options->{library}";
    make_path($output_dir);
    my $input_string = "";
    if ($options->{input} =~ /$options->{delimiter}/) {
        my @in = split(/$options->{delimiter}/, $options->{input});
        my $r1_fd = $class->Get_FD(input => $in[0]);
        my $r2_fd = $class->Get_FD(input => $in[1]);
        $input_string = qq" --paired ${r1_fd} ${r2_fd} ";
        if ($in[0] =~ /\.fastq$/) {
            $input_string = qq" --paired $in[0] $in[1] ";
        }
    }
    else {
        my $r1_fd = $class->Get_FD(input => $options->{input});
        $input_string = qq"${r1_fd} ";
        if ($options->{input} =~ /\.fastq$/) {
            $input_string = qq" $options->{input} ";
        }
    }
    my $comment = qq!## This is a kraken2 submission script
!;
    my $stdout = qq"${output_dir}/kraken.stdout";
    my $stderr = qq"${output_dir}/kraken.stderr";
    my $jstring = qq!kraken2 --db "\${KRAKEN2_DB_PATH}"/$options->{library} \\
  --report ${output_dir}/kraken_report.txt --use-mpa-style \\
  --use-names ${input_string} \\
  --classified-out ${output_dir}/classified#.fastq.gz \\
  --unclassified-out ${output_dir}/unclassified#.fastq.gz \\
  2>${stderr} \\
  1>${stdout}
# shellcheck disable=SC2181
if [ "\$?" -ne "0" ]; then
  echo "Kraken returned an error."
fi
!;

    if ($options->{clean}) {
        $jstring .= qq"
## Cleaning up after running kraken.
rm -f ${stdout}
rm -f ${output_dir}/*.fastq.gz
";
    }
    my $output = qq"${output_dir}/kraken_report.txt";
    my $kraken = $class->Submit(
        comment => $comment,
        jcpu => $options->{jcpu},
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => "kraken_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $output,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stdout => $stdout,
        stderr => $stderr,);
    my $counter = $class->Bio::Adventure::Count::Kraken_to_Matrix(
        input => $output,
        jprefix => qq"$options->{jprefix}_1",
        jdepends => $kraken->{job_id},);
    $kraken->{kraken2mtrx} = $counter;
    return($kraken);
}

=head2 C<Kraken_to_Matrix>

  Given the kraken report, create a matrix of reads by genus/species/whatever.

=over

=item C<Arguments>

  input: tsv file produced by kraken.

=back

=cut
sub Kraken_to_Matrix {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => 21,);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Path_Info($options->{input});
    $inputs = $inputs->[0];
    my $output = qq"$inputs->{directory}/$inputs->{filebase_extension}_matrix.tsv";
    my $stdout = qq"$inputs->{directory}/kraken_to_matrix.stdout";
    my $stderr = qq"$inputs->{directory}/kraken_to_matrix.stderr";
    my $log = qq"$inputs->{directory}/kraken_to_matrix.log";
    my $comment = qq"## Make a matrix of kraken hits.\n";
    my $jstring = qq?
use Bio::Adventure::Count;
my \$result = \$h->Bio::Adventure::Count::Kraken_to_Matrix_Worker(
  input => '$options->{input}',);
?;
    my $matrix_job = $class->Submit(
        comment => $comment,
        input => $options->{input},
        jdepends => $options->{jdepends},
        jcpu => 1,
        jmem => 4,
        jname => 'kraken2mtrx',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        log => $log,
        output => $output,
        stdout => $stdout,
        stderr => $stderr,);
    return($matrix_job);
}

sub Kraken_to_Matrix_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],);
    my $separator = 'g__';
    my $job_log = 'kraken_matrix.log';

    my $out = FileHandle->new(">$options->{log}");
    my $mtrx = FileHandle->new(">$options->{output}");
    my $in = FileHandle->new("<$options->{input}");
    print $out "Collating reads by ${separator}.\n";
    print $out "Reading kraken report: $options->{input}.\n";
    my %genus_observed = ();
    my $line_count = 0;
    while (my $line = <$in>) {
        chomp $line;
        next unless ($line =~ m/$separator/);
        $line_count++;
        my ($entry, $number) = split(/\t/, $line);
        my ($stuff, $genus) = split(/$separator/, $entry);
        $genus =~ s/\|.*$//g;
        $genus =~ s/\s+/_/g;
        if (defined($genus_observed{$genus})) {
            $genus_observed{$genus} = $genus_observed{$genus} + $number;
        }
        else {
            $genus_observed{$genus} = $number;
        }
    }
    $in->close();
    print $mtrx "Genus\tObserved\n";
    foreach my $s (keys %genus_observed) {
        print $mtrx "${s}\t$genus_observed{$s}\n";
    }
    $mtrx->close();
    $out->close();
    return(%genus_observed);
}

=head2 C<Mash>

 Calculate distances using mash

=cut
sub Mash {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => 21,
        length => 9,
        sketch => 9,
        jcpu => 4,);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Path_Info($options->{input});
    ## Unlike most of my jobs, the input argument here is a directory, so just grab it
    my $cwd_name = basename($options->{input});
    my $output_dir = qq"outputs/$options->{jprefix}mash_${cwd_name}/dist";
    my $sketch_dir = qq"outputs/$options->{jprefix}mash_${cwd_name}/sketch";
    my $stderr = qq"outputs/$options->{jprefix}mash_${cwd_name}.stderr";
    my $stdout = qq"outputs/$options->{jprefix}mash_${cwd_name}.stdout";
    my $comment = qq"## Playing with mash";
    my $jstring = qq!mkdir -p ${output_dir}
mkdir -p ${sketch_dir}
for fasta in $options->{input}/*; do
  name="\${fasta%.*}"
  mash sketch \\
    -s $options->{sketch} -k $options->{length} \\
    $options->{input}/\${fasta} \\
    -o ${sketch_dir}/\${name} \\
    2>>${stderr} 1>>${stdout}
done

for outer in ${sketch_dir}/*; do
  name="\${outer%.*}"
  for inner in ${sketch_dir}/*; do
    mash dist \\
      ${sketch_dir}/\${outer} \\
      ${sketch_dir}/\${inner} >> \\
      ${output_dir}/pairwise_distances_s$options->{sketch}_k$options->{length}.txt \\
      2>>${stderr} 1>>${stdout}
!;

    my $mash = $class->Submit(
        comment => $comment,
        jcpu => 4,
        jdepends => $options->{jdepends},
        jmem => 12,
        jname => qq"mash_${job_name}_$options->{length}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        stderr => $stderr,
        stdout => $stdout,);
    return($mash);
}

=head2 C<Mi_Map>

 Map reads to mature/immature miRNA.

 This function has not been used in a very long time and likely will
 require some work if I wish to use it again.

=over

=item C<Arguments>

 mirbase_data(required): File containing miRNAs from the mirbase.
 mature_fasta(required): Fasta file of the mature miRNA species.
 mi_genome(required): The set of all miRNAs expected.
 bamfile(required): Input bam file to search.

=back

=cut
sub Mi_Map {
    my ($class, %args) = @_;
    eval 'use Bio::DB::Sam; 1;';
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['mirbase_data', 'mature_fasta', 'mi_genome', 'bamfile',]);
    my $bam_base = basename($options->{bamfile}, ('.bam', '.sam'));
    my $pre_map = qq"${bam_base}_mirnadb.txt";
    my $final_map = qq"${bam_base}_mature.count";

    ## Step 1:  Read the mature miRNAs to get the global IDs and sequences of the longer transcripts.
    ## This provides the IDs as 'mm_miR_xxx_3p' and short sequences ~21 nt.
    print "Starting to read miRNA sequences.\n";
    my $sequence_db = $class->Bio::Adventure::Count::Read_Mi(seqfile => $options->{mature_fasta});

    ## Step 2:  Collate those IDs against the set of miRNA_transcripts->miRNA_mature IDs using
    ## the tab delimited file downloaded from mirbase.org
    print "Starting to read miRNA mappings.\n";
    my $sequence_mappings = $class->Bio::Adventure::Count::Read_Mappings_Mi(
        mappings => $options->{mirbase_data},
        output => $pre_map,
        seqdb => $sequence_db,);
    ## At this point, we have brought together mature sequence/mature IDs to parent transcript IDs
    ## When we read the provided bam alignment, we will simultaneously read the immature miRNA database
    ## and use them to make a bridge from the mature sequence/IDs to immature miRNA IDs.
    my $final_hits = $class->Bio::Adventure::Count::Read_Bam_Mi(
        mappings => $sequence_mappings,
        mi_genome => $options->{mi_genome},
        bamfile => $options->{bamfile},);
    my $printed = $class->Bio::Adventure::Count::Final_Print_Mi(
        data => $final_hits,
        output => $final_map,);

    my $job = $printed;
    $job->{final_hits} = $final_hits;
    return($job);
}                               ## End of Mi_Map

=head2 C<Read_Mi>

 Read an miRNA database.

 This takes the fasta file from the mirbase and extracts the IDs and
 sequences.

=over

=item C<Arguments>

 seqfile: The fasta file in question.

=back

=cut
sub Read_Mi {
    my ($class, %args) = @_;
    my $fasta = Bio::SeqIO->new(-file => $args{seqfile}, -format => 'Fasta');
    my %sequences = ();
    while (my $mi_seq = $fasta->next_seq()) {
        next unless(defined($mi_seq->id));
        my $id = $mi_seq->id;
        my $length = $mi_seq->length;
        my $seq = $mi_seq->seq;
        $sequences{$id}->{sequence} = $seq;
    }
    return(\%sequences);
}

=head2 C<Mpileup>

 Run samtools mpileup.

=cut
sub Mpileup {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        jprefix => 61,);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Path_Info($options->{input});
    ## Unlike most of my jobs, the input argument here is a directory, so just grab it
    my $cwd_name = basename($options->{input});
    my $genome = qq"$options->{libpath}/$options->{libtype}/$options->{species}.fasta";
    my $output_dir = qq"outputs/$options->{jprefix}mpileup_$options->{species}";
    my $pileup_error = qq"${output_dir}/mpileup.stderr";
    my $pileup_output = qq"${output_dir}/mpileup.bcf";
    my $final_error = qq"${output_dir}/mpileup_bcf.stderr";
    my $stdout = qq"${output_dir}/mpileup.stdout";
    my $comment = qq"## Playing with samtools mpileup";
    my $jstring = qq!mkdir -p ${output_dir}
echo "Started samtools sort at \$(date)" >> ${output_dir}/mpileup_$options->{species}.stdout

samtools sort -l 9 -@ 4 $options->{input} -o ${output_dir}/sorted.bam \\
samtools index ${output_dir}/sorted.bam
if [ \! -r "${genome}.fai" ]; then
    samtools faidx ${genome}
fi
samtools mpileup -uvf ${genome} 2>${output_dir}/samtools_mpileup.stderr \\
    ${output_dir}/sorted.vcf
bcftools view -l 9 -o ${pileup_output} 2>${output_dir}/mpileup_bcftools.stderr
!;
    my $mpileup = $class->Submit(
        comment => $comment,
        jcpu => 4,
        jdepends => $options->{jdepends},
        jmem => 12,
        jname => 'mpileup',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        stderr => $final_error,
        stdout => $stdout,);
    return($mpileup);
}

=head2 C<Read_Mappings_Mi>

 Read an miRNA database and get the connections between the various IDs, mature
 sequences, and immature sequences.

=over

=item C<Arguments>

 output: Output file to write.
 seqdb: Hash of IDs to sequences from Read_Mi().
 mappings: Hash of IDs to precursors/etc.

=back

=cut
sub Read_Mappings_Mi {
    my ($class, %args) = @_;
    my $output = $args{output};
    my $seqdb = $args{seqdb};
    my $inmap = new FileHandle($args{mappings});
    my $mimap = {};
    while (my $line = <$inmap>) {
        chomp $line;
        $line =~ s/"//g;
        my ($hit_id, $ensembl, $mirbase_id, $fivep_mature, $fivep_id, $threep_mature, $threep_id) = split(/\s+/, $line);
        if (defined($fivep_id) and $fivep_id ne "") {
            $mimap->{$fivep_id}->{mirbase_id} = $mirbase_id;
            $mimap->{$fivep_id}->{hit_id} = $hit_id;
            $mimap->{$fivep_id}->{ensembl} = $ensembl;
            $mimap->{$fivep_id}->{mimat} = $fivep_id;
        }
        if (defined($threep_id) and $threep_id ne "") {
            $mimap->{$threep_id}->{mirbase_id} = $mirbase_id;
            $mimap->{$threep_id}->{hit_id} = $hit_id;
            $mimap->{$threep_id}->{ensembl} = $ensembl;
            $mimap->{$threep_id}->{mimat} = $threep_id;
        }
    }
    $inmap->close();
    my $out = FileHandle->new(">${output}");
    ## Now re-key the database of miRNAs so that there are (potentially) multiple elements to each mirbase ID
    ## This is a bit redundant, but hopefully clearer.  The idea is to first gather all miR data, then make a list of specific
    ## ID/sequences for each mature miRNA child of the mirbase_id parent entry.  Eg:
    ## mirbase_id parent ->  [ mature_id1, mature_id2, ... ]
    ## where each mature_id is a { sequence, count, id, ensembl, ... }
    my $newdb = {};
  LOOP: foreach my $id (keys %{$seqdb}) {
        next LOOP unless ($id =~ /^mmu/);
        my @hit_list = ();
        if (defined($mimap->{$id})) {
            my $new_id = $id;
            my $sequence = $seqdb->{$id}->{sequence};
            if ($mimap->{$id}->{hit_id}) {
                $new_id = $mimap->{$id}->{hit_id};
            }
            my $element = {
                count => 0,
                ensembl => $mimap->{$id}->{ensembl},
                hit_id => $mimap->{$id}->{hit_id},
                id => $id,
                mimat => $mimap->{$id}->{mimat},
                mirbase_id => $mimap->{$id}->{mirbase_id},
                sequence => $sequence,
            };
            my $ensembl_id = $element->{ensembl};
            if (defined($newdb->{$ensembl_id})) {
                ## Then there should be an existing element, append to it.
                @hit_list = @{$newdb->{$ensembl_id}};
                push(@hit_list, $element);
                $newdb->{$ensembl_id} = \@hit_list;
                print $out "${id}; $element->{mirbase_id}; $element->{ensembl}; $element->{hit_id}; $element->{mimat}; $element->{sequence}\n";
            }
            else {
                @hit_list = ();
                push(@hit_list, $element);
                $newdb->{$ensembl_id} = \@hit_list;
            }
        }
        else {
            ## print "Did not find $id\n";
        }
    }
    $out->close();
    return($newdb);
}

=head2 C<Read_Bam_Mi>

 Read a bam file and cross reference it against an miRNA database.

=over

=item C<Arguments>

 mappings: Set of database entries to query against.
 bamfile: Input alignments to search.
 mi_genome: Fasta file containing the genomic context of the miRNAs.

=back

=cut
sub Read_Bam_Mi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
    );
    my $mappings = $args{mappings};
    my $sam = Bio::DB::Sam->new(-bam => $args{bamfile}, -fasta=> $args{mi_genome},);
    my $bam = $sam->bam;
    my $header = $bam->header;
    my $target_count = $header->n_targets;
    my $target_names = $header->target_name;
    my $align_count = 0;
    print "Started reading bamfile.\n";
    ## I probably don't need to acquire most of this information, as I am only really taking
    ## the read_seq and read_seqid
  BAM: while (my $align = $bam->read1) {
        if ($options->{debug}) {
            last BAM if ($align_count > 4000);
        }
        my $read_seqid = $target_names->[$align->tid];
        ## my $read_start = $align->pos + 1;
        ## my $read_end = $align->calend;
        ## my $read_strand = $align->strand;
        ## my $read_cigar = $align->cigar_str;
        ## my @read_scores = $align->qscore;        # per-base quality scores
        ## my $read_match_qual= $align->qual;       # quality of the match
        ## my $read_length = $align->query->end;
        my $read_seq = $align->query->seq->seq; ## maybe?
        $align_count++;

        ## Check each element in the mappings to see if they are contained in the read's sequence and vice versa.
        ## If so, increment the count for that element and move on.
        $read_seqid =~ s/^chr[\d+]_//g;
        if ($mappings->{$read_seqid}) {
            my @element_list = @{$mappings->{$read_seqid}};
            my $found_element = 0;
            my $element_length = scalar(@element_list);
            ## print "Found map, checking against ${element_length} mature RNAs.\n";
            foreach my $c (0 .. $#element_list) {
                my $element_datum = $element_list[$c];
                my $element_seq = $element_list[$c]->{sequence};
                ## print "Comparing: $read_seq vs. $element_seq\n";
                my @read_vs_element = amatch($read_seq, [ 'S1' ], ($element_seq));
                my @element_vs_read = amatch($element_seq, [ 'S1' ], ($read_seq));
                if (scalar(@read_vs_element) > 0 or scalar(@element_vs_read) > 0) {
                    $mappings->{$read_seqid}->[$c]->{count}++;
                    $found_element = 1;
                }
                ## if ($read_seq =~ /$element_seq/ or $element_seq =~ /$read_seq/) {
                ##     $mappings->{$read_seqid}->[$c]->{count}++;
                ##     print "Using an exact match: ";
                ##     print "Incremented $mappings->{$read_seqid}->[$c]->{mirbase_id} to $mappings->{$read_seqid}->[$c]->{count}\n";
                ##     $found_element = 1;
                ## }
            }
        }
    }                           ## Finish iterating through every read
    return($mappings);
}

=head2 C<Final_Print_Mi>

 Print out the final counts of miRNA mappings.

=over

=item C<Arguments>

 data: Result from cross referencing the bam/genome/miRNAs.
 output: Output filename for writing a tsv of the counts.

=back

=cut
sub Final_Print_Mi {
    my ($class, %args) = @_;
    my $final = $args{data};
    my $output = FileHandle->new(">$args{output}");
    my $hits = 0;
    foreach my $immature (keys %{$final}) {
        foreach my $mature (@{$final->{$immature}}) {
            $hits = $hits++;
            print $output "$mature->{mimat}\t$mature->{count}\n";
        }
    }
    $output->close();
    return($hits);
}                               ## End of Final_Print

=head2 C<Count_Alignments>

 Compare alignments across species.

 Count alignments across a host and parasite.  This function should
 give a sense of how many reads were single-mapped to the host and
 parasite along with multi-hits across both.  This was first used to
 distinguish between T. cruzi and human hits, ergo the 'Tc' as the
 default parasite pattern.

=over

=item C<Arguments>

 input(required): Input bam file to count against.
 species(required): Set the host or parasite species to count against.
 para_patterm('^Tc'): Pattern used to differentiate parasite-mapped
  reads.
 host_pattern(''): Pattern used to differentiate host-mapped reads.
 libpath: Change the dirname of the fasta libraries.
 libtype('genome'): Change the subdirectory of the fasta libraries.

=back

=cut
sub Count_Alignments {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        para_pattern => '^Tc',
        host_pattern => '',);
    my $result = {
        mapped => 0,
        unmapped => 0,
        multi_count => 0,
        single_count => 0,
        unmapped_count => 0,
        single_para => 0,       ## Single-hit parasite
        single_host => 0,       ## Single-hit host
        multi_host => 0, ## Multi-hit host -- keep in mind htseq will not count these.
        multi_para => 0, ## Ditto parasite
        single_both => 0, ## These are the dangerzone reads, 1 hit on both. -- false positives
        single_para_multi_host => 0,
        single_host_multi_para => 0,
        multi_both => 0, ## These have multi on both and are not a problem.
        zero_both => 0,  ## This should stay zero.
        wtf => 0,
    };

    my %group = ( para => 0, host => 0,);
    my %null_group = ( para => 0, host => 0,);
    my $fasta = qq"$options->{libpath}/$options->{libtype}/$options->{species}.fasta";
    my $sam = Bio::DB::Sam->new(-bam => $options->{input},
                                -fasta => $fasta,);
    my @targets = $sam->seq_ids;
    my $num = scalar(@targets);
    my $bam = Bio::DB::Bam->open($options->{input});
    my $header = $bam->header;
    my $target_count = $header->n_targets;
    my $target_names = $header->target_name;
    my $align_count = 0;
    my $million_aligns = 0;
    my $alignstats = qx"samtools idxstats $options->{input}";
    my @alignfun = split(/\n/, $alignstats);
    my @aligns = split(/\t/, $alignfun[0]);
    my @unaligns = split(/\t/, $alignfun[1]);
    my $number_reads = $aligns[2] + $unaligns[3];
    my $output_name = qq"$options->{input}.out";
    my $out = FileHandle->new(">${output_name}");
    print $out "There are ${number_reads} alignments in $options->{input} made of $aligns[2] aligned reads and $unaligns[3] unaligned reads.\n";
    my $last_readid = "";
  BAMLOOP: while (my $align = $bam->read1) {
        $align_count++;
        if (($align_count % 1000000) == 0) {
            $million_aligns++;
            print $out "Finished $million_aligns million alignments out of ${number_reads}.\n";
        }
        my $seqid = $target_names->[$align->tid];
        my $readid = $align->qname;
        my $start = $align->pos;
        my $end = $align->calend - 1;
        my $cigar = $align->cigar_str;
        my $strand = $align->strand;
        my $seq = $align->query->dna;
        my $qual= $align->qual;
        if ($cigar eq '') {
            $result->{unmapped}++;
            next BAMLOOP;
        }
        else {
            ## Everything which follows is for a mapped read.
            $result->{mapped}++;
            my $type;
            if ($options->{para_pattern}) {
                if ($seqid =~ m/$options->{para_pattern}/) {
                    $type = 'para';
                }
                else {
                    $type = 'host';
                }
            }
            elsif ($options->{host_pattern}) {
                if ($seqid =~ m/$options->{host_pattern}/) {
                    $type = 'host';
                }
                else {
                    $type = 'para';
                }
            }
            if ($readid ne $last_readid) {
                ## Count up what is currently in the group, then reset it.
                my $reads_in_group = $group{host} + $group{para};
                if ($reads_in_group > 1) {
                    $result->{multi_count}++;
                }
                elsif ($reads_in_group == 1) {
                    $result->{single_count}++;
                }
                else {
                    $result->{unmapped_count}++;
                }
                if ($group{host} == 0) {
                    if ($group{para} == 0) {
                        $result->{zero_both}++;
                    }
                    elsif ($group{para} == 1) {
                        $result->{single_para} = $result->{single_para} + $reads_in_group;
                    }
                    elsif ($group{para} > 1) {
                        $result->{multi_para} = $result->{multi_para} + $reads_in_group;
                    }
                    else {
                        $result->{wtf}++;
                    }
                }
                elsif ($group{host} == 1) {
                    if ($group{para} == 0) {
                        $result->{single_host} = $result->{single_host} + $reads_in_group;
                    }
                    elsif ($group{para} == 1) {
                        $result->{single_both} = $result->{single_both} + $reads_in_group;
                    }
                    elsif ($group{para} > 1) {
                        $result->{single_host_multi_para} = $result->{single_host_multi_para} + $reads_in_group;
                    }
                    else {
                        $result->{wtf}++;
                    }
                }
                elsif ($group{host} > 1) {
                    if ($group{para} == 0) {
                        $result->{multi_host} = $result->{multi_host} + $reads_in_group;
                    }
                    elsif ($group{para} == 1) {
                        $result->{single_para_multi_host} = $result->{single_para_multi_host} + $reads_in_group;
                    }
                    elsif ($group{host} > 1) {
                        $result->{multi_both} = $result->{multi_both} + $reads_in_group;
                    }
                    else {
                        $result->{wtf}++;
                    }
                }
                else {
                    $result->{wtf}++;
                }
                print "Map:$result->{mapped} Un:$result->{unmapped} SingleC:$result->{single_count} MultC:$result->{multi_count} sp:$result->{single_para} sh:$result->{single_host} mh:$result->{multi_host} mp:$result->{multi_para} SB:$result->{single_both} spmh:$result->{single_para_multi_host} shmp:$result->{single_host_multi_para} bm:$result->{multi_both} bz:$result->{zero_both} wtf:$result->{wtf}\n" if ($options->{debug});
                %group = %null_group;
                ## Now set the first read for the new group to the type of this read.
                $group{$type}++;
            }
            else {
                $group{$type}++;
            }
        }
        $last_readid = $readid;
    }                           ## End reading each bam entry
    print $out "Mapped: $result->{mapped}
Unmapped: $result->{unmapped}
Multi-mapped: $result->{multi}
Single-parasite: $result->{single_para}
Single-host: $result->{single_host}
Multi-parasite, no host: $result->{multi_para}
Multi-host, no parasite: $result->{multi_host}
DANGER Single-both: $result->{single_both}
DANGER Single-parasite, multi-host: $result->{single_para_multi_host}
DANGER Single-host, multi-parasite: $result->{single_host_multi_para}
Multi-both: $result->{both_multi}\n";
    $out->close();
    return($result);
}

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

L<htseq-count> L<Bio::DB::Sam> L<Bio::SeqIO>

=cut

sub GFF_Tree {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species'],);
    my $job_name = $class->Get_Job_Name();
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $output_dir = $paths->{output_dir};
    make_path($output_dir) unless (-d $output_dir);
    my $input_gff = $paths->{gff};
    my $gff_in = Bio::Adventure::Get_FH(input => $input_gff,
                                        suffix => qq"| grep -v '^#'");
    my $gff_io = Bio::FeatureIO->new(-format => 'gff', -fh => $gff_in);
    my %id_types = ();
    my $counters = {};
    package GFFTree;
    sub new {
        my $class = shift;
        my $options = shift;
        my $self = bless $class->SUPER::new();
        $self->attributes($options);
        return $self;
    }
    sub ids {
        my $node = shift;
        my $val = shift;
        if ($val) {
            my @ids = @{$node->attributes->{ids}};
            push(@ids, $val);
            $node->attributes->{ids} = \@ids;
        }
        else {
            return $node->attributes->{ids};
        }
    }
    sub num_ids {
        my $node = shift;
        my $val = shift;
        if ($val) {
            $node->attributes->{num_ids}++;
        }
        else {
            return $node->attributes->{num_ids};
        }
    }
    sub print_num_ids {
        $_[0]->walk_down({callback=> sub {
                              my $node = shift;
                              printf "%s%.7s\t NumIDs: %d \n",
                                  "  " x $_[0]->{_depth},
                                  ## $node->name, scalar(@{$node->attributes->{ids}}),
                                  $node->name, $node->attributes->{num_ids},
                              }, _depth => 0 });
    }
    sub by_name {
        my ($self, $name) = @_;
        my @found =();
        my $retvalue = wantarray ? 1 : 0;
        $self->walk_down({callback => sub{
            if ($_[0]->name eq $name) {
                push @found, $_[0];
                return $retvalue;
            }
            1}});
        return wantarray? @found : @found ? $found[0] : undef;
    }
    package main;
    # my $gff_tree = GFFTree->new({ ids => undef });
    my $gff_tree = GFFTree->new({ num_ids => 1 });
    my $root_name = 'hg38_111';
    $gff_tree->name($root_name);
    my $last_chromosome = 'start';
    my $feat_count = 0;
  FEAT: while (my $feature = $gff_io->next_feature()) {
        $feat_count++;
        my $chr = $feature->seq_id;
        if ($chr ne $last_chromosome) {
            print "Switched chromosome to $chr\n";
            %id_types = ();
            $last_chromosome = $chr;
        }
        my $type = $feature->primary_tag;
        next FEAT unless(defined($type));
        next FEAT if ($type eq 'biological_region');
        my @ids = $feature->get_tag_values('ID');
        my @names = $feature->get_tag_values('Name');
        my @parents = $feature->get_tag_values('Parent');
        my $id = undef;
        if (scalar(@ids) > 0) {
            $id = $ids[0];
        } elsif (scalar(@names) > 0) {
            #print "Using Name instead of ID.\n";
            $id = $names[0];
        } elsif (scalar(@parents) > 0) {
            #print "Using parent as a last resort.\n";
            $id = $parents[0];
        }
        if (!defined($id)) {
            #print "This has no ID: $feat_count\n";
            $id = 'hg38_111';
        }
        unless (defined($id_types{$id})) {
            $id_types{$id} = $type;
        }
        $feat_count % 10000 == 0 and print "Tested id types: $feat_count\n";
        my $parent = $root_name;
        if (scalar(@parents) > 0) {
            $parent = $parents[0];
        }
        my $parent_type = $id_types{$parent};
        if (!defined($parent_type)) {
            #print "This entry has no parent type: $parent $id $type $feat_count\n";
            $parent_type = 'hg38_111';
        }
        #print "Searching for: $type and $id and $parent / $parent_type\n";
        my $parent_node = undef;
        if (defined($parent)) {
            $parent_node = $gff_tree->by_name($parent_type);
        }
        my $node = undef;
        if (defined($parent_node)) {
            $node = $parent_node->by_name($type);
        } else {
            $node = $gff_tree->by_name($type);
        }
        if (defined($node) && defined($parent_node)) {
            my $daughterp = $node->is_daughter_of($parent_node);
            if ($node->is_daughter_of($parent_node)) {
                ## $node->ids($id);
                $node->num_ids($id);
                #print "Parent and node are defiend, adding $id to a set of ids.\n";
            } else {
                ## $parent_node->new_daughter({ids => [$id]})->name($type);
                $parent_node->new_daughter({num_ids => 1})->name($type);
                my $parent_name = $parent_node->name;
                #print "Parent and node are defined, but node is a new child of parent.\n";
            }
        } elsif (defined($node)) {
            ## $node->ids($id);
            $node->num_ids($id);
            #print "Only node is defiend, adding $id to a set of ids.\n";
        } elsif (defined($parent_node)) {
            ## $parent_node->new_daughter({ids => [$id]})->name($type);
            $parent_node->new_daughter({num_ids => 1})->name($type);
            my $parent_name = $parent_node->name;
            #print "Creating a new child under ${parent_name}\n";
        } else {
            ## $gff_tree->new_daughter({ids => [$id]})->name($type);
            $gff_tree->new_daughter({num_ids => 1})->name($type);
            #print "Creating a new root child of type: ${type}\n";
        }
    }
    print map "$_\n", @{$gff_tree->draw_ascii_tree};
    $gff_tree->print_num_ids;
}

1;
