package Bio::Adventure::Parsers;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use File::Basename;

## Adding this as a place to put functions which read the various formats we run across.

## Some candidates to move here soon:
## Read_Genome_Fasta, Read_Genome_GFF (though these two can likely just be removed)

sub Count_Extract_GFF {
    my ($class, %args) = @_;
    my $gff = $args{gff};
    my $log = $args{log};
    my $gff_type = $args{gff_type};
    my $gff_tag = $args{gff_tag};
    my $gff_in = Bio::Adventure::Get_FH(input => $gff, suffix => qq"| grep -v '^#'");
    ## The input gff file and parser:
    ##my $gff_in = FileHandle->new("${fc} | grep -v '^#' |");
    my $gff_io = Bio::FeatureIO->new(-format => 'gff', -fh => $gff_in);
    my $counters = {
        total_features => 0,
        id_features => 0,
        plus_features => 0,
        minus_features => 0,
        unstranded_features => 0,
        unnamed_features => 0,
        nogene_contigs => 0,
        contigs => {},
        reads => 0,
        quantified_reads => 0,
        unstranded_reads => 0,
        plus_reads => 0,
        minus_reads => 0,
        notplusminus_reads => 0,
        primary_reads => 0,
        supplemental_reads => 0,
        secondary_reads => 0,
        unmapped_reads => 0,
        unmapped_pairs => 0,
        same_strand_in_gene => 0,
        opposite_strand => 0,
        unannotated_read => 0,
    };
    my $observed = {};
    my $feature_count = 0;
  FEAT: while (my $feature = $gff_io->next_feature()) {
        my $contig_id = $feature->seq_id;
        ## print "TESTME: Contig: $contig_id\n";
        $counters->{total_features}++;
        my $tag = $feature->primary_tag;
        if (defined($gff_type)) {
            if ($feature->primary_tag ne $gff_type) {
                ## print "Skipping feature because of type mismatch.\n";
                next FEAT;
            }
        }
        my @seqnames = $feature->get_tag_values($gff_tag);
        my $gene_name = $seqnames[0];
        if (!defined($gene_name)) {
            $counters->{unnamed_features}++;
            ## print "Skipping: This feature does not have a name.\n";
            next FEAT;
        }
        my $str = $feature->strand;
        if (!defined($str)) {
            ## print "This strand is undefined for this feature: ${gene_name}.\n";
            $counters->{unstranded_features}++;
            next FEAT;
        }
        $str = '1' if ($str eq '+1');
        if ($str ne '-1' && $str ne '1') {
            print "The strand is neither +1 nor -1 for this feature: ${gene_name}, ${str}.\n";
            $counters->{unstranded_features}++;
            next FEAT;
        }

        my $info = {};
        my $contig = $feature->seq_id;
        $info->{contig} = $contig;
        if (!defined($observed->{$contig})) {
            ## print "TESTME: STarting a new contig\n";
            $observed->{$contig} = [];
            $counters->{contigs}->{$contig} = 0;
            $feature_count = -1;
        } else {
            ## print "TESTME: In existing contig: $contig\n";
        }
        $feature_count++;
        $counters->{contigs}->{$contig}++;
        $info->{primary_tag} = $feature->primary_tag;
        $info->{start} = $feature->start;
        $info->{end} = $feature->end;
        $info->{strand} = $str;
        $info->{feature} = $feature;
        $info->{most_upstream} = 0;
        $info->{observed} = {};
        $info->{next_distance} = 0;
        $info->{next_start} = 0;
        $info->{next_end} = 0;
        $info->{next_strand} = 0;
        $info->{previous_distance} = 0;
        $info->{previous_end} = 0;
        $info->{previous_strand} = 0;
        $info->{gene_name} = $gene_name;
        my @recorded_features = @{$observed->{$contig}};
        push(@recorded_features, $info);
        $observed->{$contig} = \@recorded_features;
        my $current_num = scalar(@recorded_features) - 1;
        my $previous_num = $current_num - 1;
        ## print "Checking contig num: $previous_num\n";
        if ($previous_num >= 0) {
            my $prev = $observed->{$contig}->[$previous_num];
            my $next_distance = abs($feature->start - $prev->{end});
            ## Set the next_distance for the final contig in a round-about fashion
            ## $observed->{$contig}[$contig_count]->{next_distance} = $contig_length - $feature->end;
            $observed->{$contig}->[$previous_num]->{next_distance} = $next_distance;
            $observed->{$contig}->[$previous_num]->{next_strand} = $str;
            $observed->{$contig}->[$previous_num]->{next_start} = $feature->start;
            $observed->{$contig}->[$previous_num]->{next_end} = $feature->end;
            $observed->{$contig}->[$previous_num]->{next_name} = $gene_name;
            ##print "Set next_name to: $gene_name for $prev->{gene_name}.\n";
            $info->{previous_distance} = $feature->start - $prev->{end};
            $info->{previous_end} = $prev->{end};
            $info->{previous_strand} = $prev->{strand};
            $info->{previous_name} = $prev->{gene_name};
            ## print "In previous gene: $prev->{gene_name}, set previous to $prev->{previous_name} and next to $prev->{next_name}.\n";
        } else {
            ## If there is no previous element of the contig, this is the first.
            $info->{previous_distance} = $feature->start;
        }
        ## print "Finished reading $feature_count features of the contig.\n";
    }
    print "Finished reading all features of the gff.\n";
    for my $k (keys %{$observed}) {
        my @features = @{$observed->{$k}};
        my @sorted = sort { $a->{start} <=> $b->{start} } @features;
        $observed->{$k} = \@sorted;
    }
    $gff_in->close();
    return($counters, $observed);
}

1;
