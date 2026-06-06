package Bio::Adventure::Index;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::Adventure::Config;
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Copy qw"cp";
use File::Path qw"make_path";
use File::Spec;
use File::Which qw"which";

=head1 NAME

 Bio::Adventure::Index - Invoke the various indexing utilities.

=head1 SYNOPSIS

=head2 C<BT1_Index>

 Create a bowtie1 index.

 This requires an 'input' argument.  This is a fasta file which will be used
 to write indexes in the libdir/libtype/indexes directory.  It assumes that the
 input fasta file is named in a way which is human-readable as the species.

 Thus, if the input fasta is 'lmajor_v46.fasta' and libtype is 'genome' it will
 be copied to $libdir/genome/lmajor_v46.fasta and indexes will have the basename
 'lmajor_v46'.

 This also uses less(1) to send an uncompressed copy of the input file to
 libdir/libtype/species.fasta.

=cut
sub BT1_Index {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input'],);
    my $species = basename($options->{input}, ('.gz', '.bz2', '.xz'));
    $species = basename($species, ('.fasta', '.fa'));
    my $output_dir = $paths->{output_dir};
    my $jstring = qq!mkdir -p ${output_dir}
bowtie-build $options->{input} \\
  $paths->{index} \\
  2>$outputs->{stderr} 1>$outputs->{stdout}
!;
    my $bt1_index = $class->Submit(
        jdepends => $options->{jdepends},
        jstring => $jstring,
        outputs => $outputs,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($bt1_index);
}

=head2 C<BT2_Index>

 Create a bowtie2 index.

 This requires an 'input' argument.  This is a fasta file which will be used
 to write indexes in the libdir/libtype/indexes directory.  It assumes that the
 input fasta file is named in a way which is human-readable as the species.

 Thus, if the input fasta is 'lmajor_v46.fasta' and libtype is 'genome' it will
 be copied to $libdir/genome/lmajor_v46.fasta and indexes will have the basename
 'lmajor_v46'.

 This also uses less(1) to send an uncompressed copy of the input file to
 libdir/libtype/species.fasta.

=cut
sub BT2_Index {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '',);
    my $libtype = $options->{libtype};
    my $libdir = $paths->{index_dir};
    my $species = basename($options->{input}, ('.gz', '.bz2', '.xz'));
    $species = basename($species, ('.fasta', '.fa'));
    my $copied_location = $paths->{fasta};
    if (!-f $copied_location) {
        my $fc = $class->Get_FC(input => $options->{input});
        my $copied = qx"${fc} > ${copied_location}";
    }
    my $jstring = qq!mkdir -p ${libdir}
bowtie2-build $options->{input} \\
  $paths->{index} \\
  2>$outputs->{stderr} \\
  1>$outputs->{stdout}
!;
    my $indexer = $class->Submit(
        jdepends => $options->{jdepends},
        jstring => $jstring,
        outputs => $outputs,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($indexer);
}

=head2 C<BWA_Index>

 Create bwa indexes.

=cut
sub BWA_Index {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input'],);
    my $species = basename($options->{input}, ('.fasta', '.fa'));
    my $copied_location = qq"$paths->{index}.fa";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $full_input = abs_path($options->{input});
    my $output_dir = qq"$options->{basedir}/outputs/$options->{jprefix}bwa_index";
    my $stdout = qq"${output_dir}/bwa_index.stdout";
    my $stderr = qq"${output_dir}/bwa_index.stderr";
    my $jstring = qq!
mkdir -p ${output_dir}
mkdir -p $paths->{index_dir}
start=\$(pwd)
cd $paths->{index_dir}
ln -sf ${full_input} ${species}.fa
bwa index ${species}.fa \\
  2>${stderr} \\
  1>${stdout}
cd \$start
!;
    my $bwa_index = $class->Submit(
        jdepends => $options->{jdepends},
        jstring => $jstring,
        outputs => $outputs,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($bwa_index);
}

=head2 C<Check_Blastdb>

 Check_Blastdb makes sure that there is an appropriately formatted blastdb for
 the library and type of blast performed.

=cut
sub Check_Blastdb {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        blast_tool => 'blastn',
        required => ['input'],
        type => 'prot',);
    my $libname;
    if (defined($options->{library})) {
        $libname = $options->{library};
    } else {
        $libname = $options->{input};
    }
    $libname = basename($libname, ('.fasta', '.fa', '.faa', '.fsa', '.fna'));

    ## First check for the relevant library in $ENV{BLASTDB}
    ## If it isn't there, make one in basedir/blastdb/
    my $foundlib = 0;
    my $libtype = $options->{type};
    if ($libtype ne 'prot' && $libtype ne 'nucl') {
        die(qw"makeblastdb requires either a type of 'prot' or 'nucl'.");
    }
    if ($options->{blast_tool} eq 'blastn') {
        $libtype = 'nucl';
    } elsif ($options->{blast_tool} eq 'blastp') {
        $libtype = 'prot';
    } else {
        $libtype = 'nucl';
    }

    my $test_in = Bio::SeqIO->new(-file => $options->{input}, -format => 'Fasta');
    my $guesses = {
        prot => 0,
        nucl => 0,
    };
    my $sum = 0;
    while (my $seq = $test_in->next_seq) {
        $sum++;
        my $guess = $seq->alphabet;
        if ($guess eq 'protein') {
            $guesses->{prot}++;
        } else {
            $guesses->{nucl}++;
        }
    }
    my $guessed = 'nucl';
    if ($guesses->{prot} > $guesses->{nucl}) {
        $guessed = 'prot';
    }
    print "Out of ${sum} sequences, $guesses->{prot} were guessed prot ";
    print "and $guesses->{nucl} were guessed nucleotide.\n";
    if ($guessed eq $libtype) {
        print "The guessed type and provided type agree!\n";
    } else {
        print "The guessed type and provided type disagree, going with: $guessed.\n";
        $libtype = $guessed;
    }

    my $checklib = qq"${libname}.psq";
    my $checklib_zero = qq"${libname}.00.psq";
    if ($libtype eq 'nucl') {
        $checklib_zero = qq"${libname}.00.nsq";
        $checklib = qq"${libname}.nsq";
    }
    my %modules = Get_Modules(caller => 1);
    my $loaded = $class->Module_Loader(%modules);
    my $db_directory = $ENV{BLASTDB};
    my $unloaded = $class->Module_Reset(env => $loaded);
    my $lib = '';
    my $relative_directory = 'blastdb';
    if (!defined($db_directory)) {
        $ENV{BLASTDB} = "$options->{basedir}/blastdb";
        $db_directory = "$options->{basedir}/blastdb";
    } else {
        $relative_directory = qq"$ENV{BLASTDB}";
    }

    print "Looking for ${checklib} / ${checklib_zero} in either $ENV{BLASTDB} or $options->{basedir}/blastdb.\n";
    ## Start with BLASTDB
    if (-f "${db_directory}/${checklib}" or
        -f "${db_directory}/${checklib_zero}") {
        $foundlib++;
        print "Found an existing blast database at ${libname}.\n";
        return($libname);
    } else {
        print "Did not find an existing blast database.\n";
    }

    ## If we do not find the blast database, create it in the basedir.
    if (!$foundlib) {
        if (!-d $db_directory) {
            make_path($db_directory);
        }
        my $formatdb_command = qq"
export BLASTDB=$ENV{BLASTDB}
makeblastdb \\
  -in $options->{input} \\
  -dbtype ${libtype} \\
  -out ${db_directory}/${libname} \\
  2>$outputs->{stderr} \\
  1>$outputs->{stdout}";
        print "The makeblastdb command run is: ${formatdb_command}\n";
        my $formatdb_ret = qx"${formatdb_command}";
    }
    my $final_directory = qq"${db_directory}/${libname}";
    return($final_directory);
}

=head2 C<Extend_Kraken_DB>

 Add some more sequences to an existing kraken2 database.

=cut
sub Extend_Kraken_DB {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',);
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $jstring = qq!
kraken2-build --download-taxonomy --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>$outputs->{stderr} \\
              1>$outputs->{stdout}
kraken2-build --add-to-library $options->{input} --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>$outputs->{stderr} \\
              1>>$outputs->{stdout}
kraken2-build --download-library $options->{library} \\
              --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>$outputs->{stderr} \\
              1>>$outputs->{stdout}
kraken2-build --build --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>$outputs->{stderr} \\
              1>>$outputs->{stdout}
!;
    my $kraken = $class->Submit(
        jdepends => $options->{jdepends},
        jstring => $jstring,
        outputs => $outputs,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($kraken);
}

=head2 C<Hisat2_Index>

 Create a hisat2 index using ${species}.fasta and leaves it in the indexes/
 directory.

=cut
sub Hisat2_Index {
    my ($class, %args) = @_;
    my $function_info = $class->Bio::Adventure::Config::Get_Function_Vars(\%args);
    my %options = %{$function_info};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];
    my $index_directory = dirname($options{index});
    if (!-d $index_directory) {
        make_path($index_directory);
    }
    my $updated = $class->Bio::Adventure::Config::Update_State(
        options => \%options,);
    my $index = $class->Submit();
    return($index);
}

=head2 C<Kallisto_Index>

 Use kallisto and an annotated_CDS fasta sequence library to create an index.

=cut
sub Kallisto_Index {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        jprefix => '15',
        required => ['input'],);
    my $cds = basename($options->{input}, ('.fasta', '.fa'));
    my $species = $cds;
    $species =~ s/_cds//g;
    $species =~ s/_nt//g;
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${cds}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $libtype = $options->{libtype};
    my $input = File::Spec->rel2abs($options->{input});
    my $jstring = qq!kallisto index -i $paths->{index_file} \\
  ${input} \\
  2>$outputs->{stderr} 1>$outputs->{stdout}
!;
    my $ka_index = $class->Submit(
        jdepends => $options->{jdepends},
        jstring => $jstring,
        outputs => $outputs,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($ka_index);
}

=head2 C<Make_Codon_Table>

 Assist caical by creating a codon table for an arbitrary species in
 the expected format.

 Given a reference assembly, print out a codon table for use with a
 codon adapatation calculator.  I wrote a little piece of code which
 works from a fasta/gff pair, but currently it is dumb.  I ought to
 improve its logging and have it return a useful/interesting data
 structure.

=cut
sub Make_Codon_Table {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['species'],);
    ## I have a few suffixes for writing genbank files.
    my @potential_suffixes = ('gbff', 'gbk', 'gbf', 'gb', 'genbank');
    my @compressions = ('gz', 'xz', 'bz2');
    my $in_gbff = '';
  POTENTIAL: for my $potential (@potential_suffixes) {
        my $start = qq"$paths->{gbk_dir}/$options->{species}.${potential}";
        print "SEARCHING FOR $start\n";
        if (-r $start) {
            $in_gbff = $start;
            last POTENTIAL;
        }
        for my $c (@compressions) {
            my $comp_test = qq"${start}.${c}";
            if (-r $comp_test) {
                $in_gbff = $comp_test;
                last POTENTIAL;
            }
        }
    }
    ## The table already exists, move on -- or maybe I should have it overwrite?
    print "Checking for the existence of a codon table: $paths->{codon_table}\n";
    if (-r $paths->{codon_table}) {
        return(1);
    }

    ## This is a little dumb, but an easy way to think through writing out the table.
    ## E.g. I will do a for(for(for())) over these to get the 64 codons.
    my @first = ('T', 'C', 'A', 'G');
    my @second = ('T', 'C', 'A', 'G');
    my @third = ('T', 'C', 'A', 'G');
    ## Set up the pieces which will hold the data of interest.
    my $total_codons = 0;
    my %codon_counts = ();
    print "Getting a FH for ${in_gbff}\n";
    my $in = $class->Get_FH(input => $in_gbff);
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
    my $seq_count = 0;
  SEQ: while (my $seq = $seqio->next_seq) {
        $seq_count++;
        ## print "Starting sequence $seq_count\n";
        my @feature_list = $seq->get_SeqFeatures();
        my $f_count = 0;
      FEAT: for my $f (@feature_list) {
            next FEAT unless ($f->primary_tag eq 'CDS');
            $f_count++;
            ## print "Feature: $f_count\n";
            my $sequence_string = $f->seq->seq;
            my $trans = $f->seq->translate->seq;
            my @seq_str = split(//, $sequence_string);
            my @tr_str = split(//, $trans);
          CHOMP: while (scalar(@tr_str) > 0) {
                my $amino = shift @tr_str;
                my $nt1 = shift @seq_str;
                my $nt2 = shift @seq_str;
                my $nt3 = shift @seq_str;
                my $codon = qq"${nt1}${nt2}${nt3}";
                if (defined($codon_counts{$codon})) {
                    $codon_counts{$codon}++;
                } else {
                    $codon_counts{$codon} = 1;
                }
                $total_codons++;
            } ## Done pulling apart the sequence arrays.
        } ## Iterating over the features in this sequence.
    } ## End going through the sequences of the assembly.
    $in->close();

    if ($total_codons == 0) {
        print "This failed: $options->{species}\n";
        return(undef);
    }
    my $divisor = 1000.0 / $total_codons;
    my $table = FileHandle->new(">$paths->{codon_table}");
    for my $f (@first) {
        for my $s (@second) {
            my $string = '';
            for my $t (@third) {
                my $codon = qq"${f}${s}${t}";
                my $per_thousand = sprintf("%.1f", $codon_counts{$codon} * $divisor);
                $string .= qq"${codon} ${per_thousand}($codon_counts{$codon})   ";

            }
            print $table qq"${string}\n";
        }
        print $table "\n";
    }
    $table->close();
}

=over

=item C<RSEM_Index>

 Use RSEM and an annotated_CDS fasta sequence library to create a transcript index.

=cut
sub RSEM_Index {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input'],);
    my $species = basename($options->{input}, ('.fasta', '.fa'));
    $species =~ s/_cds//g;
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $jstring = qq!mkdir -p $paths->{output_dir}
rsem-prepare-reference --bowtie2 $options->{input} ${species} \\
  2>$outputs->{stderr} 1>$outputs->{stdout}
!;
    my $jobid = $class->Submit(
        jdepends => $options->{jdepends},
        jstring => $jstring,
        outputs => $outputs,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($jobid);
}

=back

=head2 C<Salmon_Index>

 Invoke salmon with an annotated_CDS fasta sequence library to create a
 transcript index.  Note that newer version of salmon would like a set
 of decoys, which may be performed in a couple of ways, the second of
 which I am copy/pasting from the documentation.

 The second is to use the entire genome of the organism as the decoy
 sequence. This can be done by concatenating the genome to the end of
 the transcriptome you want to index and populating the decoys.txt
 file with the chromosome names. Detailed instructions on how to
 prepare this type of decoy sequence is available here. This scheme
 provides a more comprehensive set of decoys, but, obviously, requires
 considerably more memory to build the index

=cut
sub Salmon_Index {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        decoy => 1,);
    my $genome = File::Spec->rel2abs($options->{input});
    my $cds = basename($options->{input}, ('.fasta', '.fa'));
    my $cds_dir = dirname($options->{input});
    my $species = $cds;
    my $log = FileHandle->new(">$paths->{log}");
    ## Drop the suffixes which might be annoying.
    $species =~ s/_cds//g;
    $species =~ s/_nt//g;
    my $species_file = qq"${cds_dir}/${species}.fasta";
    my $index_basedir = qq"$options->{libpath}/$options->{libtype}";
    my $species_location = qq"${index_basedir}/${species}.fasta";
    my $decoy_copy_string = qq'';
    my $jstring = qq'';
    my $index_string = qq!
salmon index -t $options->{input} \\
  -i $paths->{index_dir} !;
    if ($options->{decoy}) {
        if (!-f $species_location) {
            cp($species_file, $species_location);
        }
        my $decoy_location = qq"$options->{libpath}/$options->{libtype}/fasta/$options->{species}_decoys.fasta";
        my $decoy_txt_name = basename($decoy_location, ('.fasta'));
        my $decoy_txt_dir = dirname($decoy_location);
        my $decoy_text_file = qq"${decoy_txt_dir}/${decoy_txt_name}.txt";
        my $genome_location = qq"$options->{libpath}/genome/fasta/$options->{species}.fasta";
        ## Write out the required concatenation of the transcripts+genome.
        print $log "Writing transcripts and decoy genome to: ${decoy_location}\n";
        my $decoys = FileHandle->new(">${decoy_location}");
        my $tx = FileHandle->new("<$options->{input}");
        while (my $line = <$tx>) {
            print $decoys $line;
        }
        $tx->close();
        my $genome = FileHandle->new("<${genome_location}");
        print $log "Writing modified fasta headers to the decoy text file: ${decoy_text_file}\n";
        my $decoy_txt = FileHandle->new(">${decoy_text_file}");
        while (my $line = <$genome>) {
            chomp $line;
            if ($line =~ /^>/) {
                $line =~ s/\s+|:|\|//g;
                my $decoy_line = $line;
                $decoy_line =~ s/^\>//g;
                print $decoy_txt "$decoy_line\n";
            }
            print $decoys "$line\n";
        }
        $genome->close();
        $decoy_txt->close();
        print "Writing indexes to: $paths->{index_dir}\n";
        $jstring = qq!salmon index -t ${decoy_location} \\
  -i $paths->{index_dir} \\
  --decoys ${decoy_text_file} \\
  2>$outputs->{stderr} \\
  1>$outputs->{stdout}
if [[ \$? -eq "0" ]]; then
  echo "Salmon indexing completed successfully, deleting the decoys fasta and text file."
  rm ${decoy_text_file}
  rm ${decoy_location}
else
  echo "Salmon indexing completed with \$?, not deleting decoys fasta and text."
fi
!;
    } else {
        warn("This function would prefer to make a decoy aware index set which requires the full genome.");
        say("Waiting 3 seconds to see if you want to quit and gather that genome,
otherwise a decoy-less index will be generated.");
        sleep(3);
        $jstring = qq!${index_string} \\
  2>$outputs->{stderr} \\
  1>$outputs->{stdout}
!;
    }
    $log->close();
    my $jobid = $class->Submit(
        jdepends => $options->{jdepends},
        jstring => $jstring,
        outputs => $outputs,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($jobid);
}

=head2 C<STAR_Index>

 Create indexes appropriate for STAR.

=cut
sub STAR_Index {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['input',],);
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = basename($options->{input}, ('.fasta', '.fa'));
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $jstring = qq!STAR \\
  --runMode genomeGenerate \\
  --runThreadN 12 \\
  --genomeDir $paths->{index_dir} \\
  --genomeFastaFiles $options->{libdir}/${libtype}/${species}.fasta \\
  --sjdbGTFfile $options->{libdir}/${libtype}/${species}.gtf \\
  --limitGenomeGenerateRAM 160000000000 \\
  2>$outputs->{stderr} 1>$outputs->{stdout}

!;
    my $jobid = $class->Submit(
        jdepends => $options->{depends},
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        outputs => $outputs,);
    return($jobid);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
