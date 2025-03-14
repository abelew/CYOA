package Bio::Adventure::Convert;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use feature 'try';
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::FeatureIO;
use Bio::Tools::GFF;
use Bio::Root::Exception;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Which qw"which";
use List::MoreUtils qw"uniq";
use Text::CSV_XS::TSV;

no warnings qw"experimental::try";

=head1 NAME

 Bio::Adventure::Convert - Perform conversions between various formats.


=head1 SYNOPSIS

 The functions here handle the various likely conversions one may perform.
 sam to bam, gff to fasta, genbank to fasta, etc.

=head1 METHODS

=cut

sub Any2Any {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        informat => 'fastq',
        outformat => 'fasta',
        jprefix => '78');
    my $base = basename($options->{input}, ('.gz', '.xz', '.bz2'));
    (my $base_no_extension = $base) =~ s/\.[^.]+$//;
    my $dir = dirname($options->{input});
    my $output_file = qq"${dir}/${base}.$options->{outformat}";
    my $in_fh = FileHandle->new("less $options->{input} |");
    my $input_seq = Bio::SeqIO->new(-fh => $in_fh,
                                    -format => $options->{informat},);
    my $output_seq = Bio::SeqIO->new(-file   => ">${output_file}",
                                     -format => $options->{outformat},);
    while (my $input_seq = $input_seq->next_seq) {
        $output_seq->write_seq($input_seq);
    }
    $in_fh->close();
}

=back

=head2 C<Gff2Fasta>

 It writes one file of amino acid sequence and one of nucleotides.
 Upon completion, it returns the number of entries written.

=over

=item C<Arguments>

 species: Add the .gff and .fasta suffixes to this and look in $libdir.
 gff(required): Gff file to read.
 input(required): Genome fasta file to read.
 gff_tag('gene_id'): ID tag to use when extracting features.
 gff_type(undef): gff type to extract, left undefined it will count
  them up and choose the most highly represented.

=item C<Invocation>

> cyoa --method gff2 --input lmajor.fasta \
   --gff lmajor.gff --gff_type CDS --gff_tag ID

=cut
sub Gff2Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        gff_tag => 'gene_id',
        gff_type => 'cds',);
    unless ($options->{species} || ($options->{input} && $options->{gff})) {
        die("This requires either a species or an input gff/fasta.")
    }
    my $genome;
    my $gff;
    my $species;
    if ($options->{species}) {
        my $input_base = qq"$options->{libpath}/genome/$options->{species}";
        $genome = qq"${input_base}.fasta";
        $gff = qq"${input_base}.gff";
        $species = $options->{species};
    } elsif ($options->{input} && $options->{gff}) {
        $genome = $options->{input};
        $gff = $options->{gff};
        $species = basename($genome, split(/,/, $class->{suffixes}));
        $species = basename($genome, ('.fasta', '.fna', '.faa'));
    }

    my $wanted_tag = $options->{gff_tag};
    my $wanted_type = $options->{gff_type};
    my $chromosomes = $class->Read_Genome_Fasta(genome => $genome);
    my $gff_handle = Bio::Adventure::Get_FH(input => $gff);
    my $nt_file = qq"${species}_${wanted_type}_${wanted_tag}_nt.fasta";
    my $aa_file = qq"${species}_${wanted_type}_${wanted_tag}_aa.fasta";
    my $nt_out = Bio::SeqIO->new(
        -format => 'Fasta',
        -file => qq">${nt_file}",);
    my $aa_out = Bio::SeqIO->new(
        -format => 'Fasta',
        -file => qq">${aa_file}",);
    ## Note that this and the next line might not be a good idea,
    ## HT_Types only looks at the first n (40,000) records and uses that as a heuristic
    ## to see that the wanted type is actually in the gff file.
    my $feature_type = $class->Bio::Adventure::Count::HT_Types(
        annotation => $gff,
        feature_type => $wanted_type,);
    $feature_type = $feature_type->[0];
    my $annotation_in = Bio::Tools::GFF->new(-fh => $gff_handle, -gff_version => 3);
    my $features_written = 0;
    my $features_read = 0;
  LOOP: while (my $feature = $annotation_in->next_feature()) {
      $features_read++;
      my $primary_type = $feature->primary_tag();
      unless ($primary_type eq $wanted_type) {
          next LOOP;
      }
      my $start = $feature->start();
      my $end = $feature->end();
      my $strand = $feature->strand();
      my $gff_chr = $feature->seq_id();
      my @ids;
      my $e;
      try {
          @ids = $feature->each_tag_value($wanted_tag);
      }
      catch ($e) {
          print "Did not find the tag: ${wanted_tag}, perhaps check the gff file.\n";
          next LOOP;
      }
      my $id = $ids[0];
      my $genome_obj = $chromosomes->{$gff_chr}->{obj};
      my $cds = $genome_obj->trunc($start, $end);
      if ($strand == -1) {
          $cds = $cds->revcom;
      }
      my $id_string = qq"chr_${gff_chr}_id_${id}_start_${start}_end_${end}";
      my $nt_sequence = Bio::Seq->new(
          -id => $id_string,
          -seq => $cds->seq());
      $nt_out->write_seq($nt_sequence);
      my $aa_sequence = Bio::Seq->new(
          -id => $id_string,
          -seq => $cds->translate->seq());
      $aa_out->write_seq($aa_sequence);
      $features_written++;
  } ## End LOOP
    $gff_handle->close();
    print "Wrote ${features_written} features to the aa/nt files.\n";
    my $retlist = {
        features_written => $features_written,
        output => $nt_file,
        output_nt => $nt_file,
        output_aa => $aa_file,
    };
    return($retlist);
}

=back

=head2 C<Gff2Gtf>

 Pretty much unused gff to gtf converter.

 Reads a given gff file and writes a gtf file from the features
 found therein.  It returns the number of features written.

 note that I only use gtf files for tophat, thus they must have a tag 'transcript_id'!
 This is woefully untrue for the tritrypdb gff files.  Thus I need to have a regex in this
 to make sure that web_id or somesuch is changed to transcript_id.

=cut
sub Gff2Gtf {
    my ($class, %args) = @_;
    my $input = $args{gff};

    my ($name, $path, $suffix) = fileparse($input, qr/\.gff/);
    my $out_file = $path . $name . '.gtf';

    my $in_gff = Bio::FeatureIO->new(-file => "${input}",
                                     -format => 'GFF',
                                     -version => 3);
    my $out_gtf = FileHandle->new(">${out_file}");
    local $| = 1;
    ##my $out_gtf = new FileHandle();
    ##$out_gtf->open(">$out_file");

    my $features_written = 0;
    my $in_features = 0;
  FEATURES: while (my $feature = $in_gff->next_feature()) {
      $in_features++;
      ## the output handle is reset for every file
      ## According to UCSC, GTF file contains 9 column:
      ## <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
      ## An example working gtf line: (except the double spaces are tabs...)
      ## TcChr20-S  TriTrypDB  CDS  14765  15403  .  -  0    gene_id "cds_TcCLB.397937.10-1"; transcript_id "cds_TcCLB.397937.10-1";

      my @tags = $feature->get_all_tags();
      my $seqid = $feature->seq_id();
      my $location = $feature->location();
      my $start = $feature->start();
      my $end = $feature->end();
      my $strand = $feature->strand();
      if ($strand == 1) {
          $strand = '+';
      } else {
          $strand = '-';
      }
      my $source = $feature->source_tag();
      my $primary_id = $feature->primary_tag();
      my $phase = '.';
      if ($primary_id eq 'CDS') {
          $phase = $feature->phase();
      }
      next FEATURES if ($primary_id eq 'supercontig' or
                        $primary_id eq 'region' or
                        $primary_id eq 'chromosome');
      my $string = qq"${seqid}\t${source}\t${primary_id}\t${start}\t${end}\t.\t${strand}\t${phase}\t";

      my $last_column = '';
      my $annot = $feature->{annotation};
      my $stringified;
      my $transcript_string = '';
      my $geneid_string = '';
      foreach my $key ($annot->get_all_annotation_keys()) {
          my @values = $annot->get_Annotations($key);
          foreach my $value (@values) {
              $stringified = $value->{value};
              $stringified = $value->{term}->{name} unless($stringified);
              $last_column .= qq!${key} "${stringified}"; !;
          }
          if ($key eq 'ID') {
              $transcript_string = qq!transcript_id "${stringified}"; !;
          }
          if ($key eq 'seq_id') {
              $geneid_string = qq!gene_id "${stringified}"; !;
          }
      }
      my $seq_string = $string . ${transcript_string} . ${geneid_string} . $last_column;
      $seq_string =~ s/\s+$//g;
      print $out_gtf $seq_string;
      $features_written++;
  } ## End iterating over every FEATURES
    $out_gtf->close();
    return($features_written);
}

=head2 C<Read_GFF>

 use Bio::Tools::GFF to gather information from a gff file.

 This pulls information of potential interest from a gff file.  In
 theory it could be trivially modified to be robust to different
 formats and other shenanigans.

=cut
sub Read_GFF {
    my ($class, %args) = @_;
    my $chromosomes = $args{chromosomes};
    my $gff = FileHandle->new("<$args{gff}");

    my $annotation_in = Bio::Tools::GFF->new(-fh => \$gff, -gff_version => 3);
    my $gff_out = {};
    print "Starting to read gff: $args{gff}\n";
  LOOP: while(my $feature = $annotation_in->next_feature()) {
      next LOOP unless ($feature->{_primary_tag} eq $args{gff_type});
      my $location = $feature->{_location};
      my $start = $location->start();
      my $end = $location->end();
      my $strand = $location->strand();
      my @ids = $feature->each_tag_value('ID');
      my $id = "";
      my $gff_chr = $feature->{_gsf_seq_id};
      my $gff_string = $annotation_in->gff_string($feature);
      if (!defined($chromosomes->{$gff_chr})) {
          print STDERR "Something is wrong with $gff_chr\n";
          next LOOP;
      }
      foreach my $i (@ids) {
          $i =~ s/^cds_//g;
          $i =~ s/\-\d+$//g;
          $id .= "${i} ";
      }
      $id =~ s/\s+$//g;
      my @gff_information = split(/\t+/, $gff_string);
      my $description_string = $gff_information[8];
      my $orf_chromosome = $gff_chr;
      my $annot = {
          id => $id,
          start => $start, ## Genomic coordinate of the start codon
          end => $end, ## And stop codon
          strand => $strand,
          description_string => $description_string,
          chromosome => $gff_chr,
      };
      $gff_out->{$gff_chr}->{$id} = $annot;
  } ## End looking at every gene in the gff file
    $gff->close();
    return($gff_out);
}

=head2 C<Sam2Bam>

 Sort, compress, index a sam file into a usable bam file.

 Used to invoke samtools to take the sam output from bowtie/bwa and
 convert it to an compressed-sorted-indexed bam alignment. This
 function just calls $class->Samtools(), but has a little logic
 to see if the invocation of this is for an extent .sam file or
 calling it on an existing .fastq(.gz), in which case one must
 assume it is being called on one or more files in the bowtie_out/
 directory and which start with $basename, include something like
 -trimmed-v0M1.sam.

=over

=item C<Arguments>

 input(required): Input sam file.
 species(required): Species name to hunt down the appropriate
  fasta/gff files.
 modules('samtools', 'bamtools'): Environment modules to load.

=item C<Invocation>

> cyoa --task convert --method sam2bam --input bob.sam --species hg38_91

=cut
sub Sam2Bam {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],);
    my @input_list = ();
    my $paths = $class->Get_Path_Info($options->{input});
    if ($options->{input}) {
        push(@input_list, $options->{input});
    } elsif (-r $options->{input} and $options->{input} =~ /\.sam$/) {
        push(@input_list, $options->{input});
    } elsif (-r $options->{input} and $options->{input} =~ /\.fastq$/) {
        if (-d "bowtie_out") {
            find({ wanted => sub { push(@input_list, "bowtie_out/$_") if ($_ =~ /\.sam$/); }, follow => 1 }, 'bowtie_out/');

        } else {
            foreach my $k (%{$options->{bt_args}}) {
                my $output_string = qq"bowtie_out/$paths->{jbasename}-${k}.sam";
                push(@input_list, $output_string);
            }
            my $bt = $class->Bio::Adventure::Map::Bowtie(%args);
        }
    } else {
        die("I don't know what to do without a .fastq file nor a .sam file.\n");
    }
    my $loaded = $class->Module_Load(modules => $options->{modules});
    my $sam = $class->Bio::Adventure::Convert::Samtools(%args,
        jdepends => $options->{jdepends},
        sam => \@input_list);
    return($sam);
}

=back

=head2 C<Samtools>

 Does the actual work of converting a sam to a bam file.

 $hpgl->Samtools() calls (in order): samtools view, samtools sort,
 and samtools index.  Upon completion, it invokes bamtools stats to
 see what the alignments looked like.

 It explicitly does not pipe one samtools invocation into the next,
 not for any real reason but because when I first wrote it, it
 seemed like the sorting was taking too long if I did not already
 have the alignments in a bam file.

=over

=item C<Arguments>

 input(required): Input sam file.
 species(required): Input species prefix for finding a genome.

=cut
sub Samtools {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        jmem => 30,
        jname => 'sam',
        jprefix => '',
        paired => 1,
        mismatch => 1,);
    my $input = $options->{input};
    my $output = $input;
    $output =~ s/\.sam$/\.bam/g;
    my $sorted_name = $input;
    $sorted_name =~ s/\.sam$//g;
    $sorted_name = qq"${sorted_name}-sorted";
    my $paired_name = $sorted_name;
    $paired_name =~ s/\-sorted/\-paired/g;
    my $workdir = dirname($input);
    ## Add a samtools version check because *sigh*
    my $samtools_version = qx"samtools 2>&1 | grep 'Version'";
    ## Start out assuming we will use the new samtools syntax.

    my $stderr = qq"${output}_samtools.stderr";
    my $stdout = qq"${output}_samtools.stdout";
    my $samtools_first = qq!
## If a previous sort file exists due to running out of memory,
## then we need to get rid of them first.
## hg38_100_genome-sorted.bam.tmp.0000.bam
if [[ -f "${output}.tmp.000.bam" ]]; then
  rm -f ${output}.tmp.*.bam
fi
samtools view -u -t $options->{libdir}/genome/$options->{species}.fasta \\
  -S ${input} -o ${output}  \\
  2>${stderr} \\
  1>${stdout}
!;

    my $samtools_second = qq"samtools sort -l 9 ${output} \\
  -o ${sorted_name}.bam \\
  2>>${stderr} \\
  1>>${stdout}";
    ## If there is a 0.1 in the version string, then use the old syntax.
    if ($samtools_version =~ /0\.1/) {
        $samtools_first = qq"samtools view -u \\
  -t $options->{libdir}/genome/$options->{species}.fasta \\
  -S ${input} 1>${output}";
        $samtools_second = qq"samtools sort -l 9 ${output} \\
  ${sorted_name} \\
  2>>${stderr} \\
  1>>${stdout}";
    }
    my $jstring = qq!
echo "Starting samtools"
if [[ -f "${output}" && -f "${input}" ]]; then
  echo "Both the bam and sam files exist, rerunning."
elif [[ -f "${output}" ]]; then
  echo "The output file exists, quitting."
  exit 0
elif [[ \! -f "${input}" ]]; then
  echo "Could not find the samtools input file."
  exit 1
fi
${samtools_first}
echo "First samtools command finished with \$?"
${samtools_second}
rm ${output}
rm ${input}
mv ${sorted_name}.bam ${output}
samtools index ${output} \\
  2>>${stderr} \\
  1>>${stdout}
echo "Second samtools command finished with \$?"
bamtools stats -in ${output} \\
  2>>${output}_samtools.stats 1>&2
echo "Bamtools finished with \$?"
!;
    if ($options->{paired}) {
        $jstring .= qq!
## The following will fail if this is single-ended.
samtools view -b -f 2 \\
  -o ${paired_name}.bam \\
  ${output} \\
  2>>${stderr} \\
  1>>${stdout}
samtools index ${paired_name}.bam \\
  2>>${stderr} \\
  1>>${stdout}
bamtools stats -in ${paired_name}.bam \\
  2>>${output}_samtools.stats 1>&2
!;
    } else {
        ## If it is not paired, just set the paired output name to the output name so that
        ## jobs which run this as part of a chain do not need to go looking for the paired output.
        ## I am considering this primarily in the case of gatk deduplication which demands
        ## properly paired reads or single ended; therefore if I set this here I do not need to
        ## check later for different inputs.
        $paired_name = basename($output, ('.bam'));
    }

    unless ($options->{mismatch}) {
        $jstring .= qq!
bamtools filter -tag XM:0 \\
  -in ${output} \\
  -out ${sorted_name}_nomismatch.bam \\
  2>>${output}_samtools.stats 1>&2
echo "bamtools filter finished with: \$?"
samtools index \\
  ${sorted_name}_nomismatch.bam \\
  2>>${stderr} \\
  1>>${stdout}
echo "final samtools index finished with: \$?"
!;
    }

    if ($options->{samtools_mapped}) {
        my $mapped = $input;
        $mapped =~ s/\.sam/_samtools_mapped\.fastq/g;
        $jstring .= qq"
samtools fastq ${output} -F 4 \\
  1>${mapped} \\
  2>${mapped}.stderr
xz -9e -f ${mapped}
";
    }

    if ($options->{samtools_unmapped}) {
        my $unmapped = $input;
        $unmapped =~ s/\.sam/_samtools_unmapped\.fastq/g;
        $jstring .= qq"
samtools fastq ${output} -f 4 \\
  1>${unmapped} \\
  2>${unmapped}.stderr
xz -9e -f ${unmapped}
";
    }

    my $comment = qq!## Converting the text sam to a compressed, sorted, indexed bamfile.
## Also printing alignment statistics to ${output}.stats
!;
    if ($options->{jdepends}) {
        $comment .= qq!
## This job depended on: $options->{jdepends}!;
    }
    my $jobname = qq"$options->{jname}_$options->{species}";
    my $samtools = $class->Submit(
        comment => $comment,
        depends => $options->{jdepends},
        input => $input,
        output => $output,
        paired => $options->{paired},
        paired_output => qq"${paired_name}.bam",
        stderr => $stderr,
        stdout => $stdout,
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        jcpu => 4,
        jmem => $options->{jmem},
        jname => $jobname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => '18:00:00',);
    return($samtools);
}

sub Species2SF {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species'],);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();

    ## Read the GFF
    print "Reading $paths->{gff}\n";
    my $gff_fh = Bio::Adventure::Get_FH(input => $paths->{gff}, suffix => '| grep -v "^#"');
    my $sf_by_contig = {};
    my $input_gff = new Bio::FeatureIO(-format => 'GFF', -version => 3, -fh => $gff_fh);
  INFEAT: while (my $feature = $input_gff->next_feature()) {
      my @sf_tmp_features = ();
      my $contig_id = $feature->seq_id;
      if (defined($sf_by_contig->{$contig_id})) {
          @sf_tmp_features = @{$sf_by_contig->{$contig_id}};
      } else {
          $sf_by_contig->{$contig_id} = ();
      }
      push(@sf_tmp_features, $feature);
      $sf_by_contig->{$contig_id} = \@sf_tmp_features;
  }
    $gff_fh->close();
    my @contig_ids = sort keys %{$sf_by_contig};
    my $num_contigs = scalar(@contig_ids);

    ## Read the genome
    print "Reading genome: $paths->{fasta}\n";
    my $input_genome = Bio::SeqIO->new(-file => $paths->{fasta}, -format => 'Fasta');
    my $sf_sequences = {};
    my $contig_features = {};
  GENOME: while (my $in_seq = $input_genome->next_seq()) {
      ## Create a feature containing the entire contig.
      my $id = $in_seq->id;
      my $contig_feature = new Bio::SeqFeature::Generic(
          -start => 1,
          -end => $in_seq->length,
          -id => $id,
          -strand => 0,
          -primary => 'contig',
          -display_name => $in_seq->display_name,);
      my $attached = $in_seq->add_SeqFeature($contig_feature);
      $contig_features->{$id} = $contig_feature;
      my @in_contig_features = @{$sf_by_contig->{$id}};
      my $num_gff_features = scalar(@in_contig_features);
      for my $in_feature (@in_contig_features) {
          my $in_gff_attached = $in_seq->add_SeqFeature($in_feature);
      }
      $sf_sequences->{$id} = $in_seq;

  }
    my $ret = {
        contigs => $contig_features,
        sf_sequences => $sf_sequences,
        sf_by_contig => $sf_by_contig, };
    return($ret);
}

=back

=head1 AUTHOR - atb

Email abelew@gmail.com

=head1 SEE ALSO

L<samtools> L<Bio::FeatureIO> L<Bio::Tools::GFF> L<Bio::Seq> L<Bio::SeqIO>

=cut

1;
