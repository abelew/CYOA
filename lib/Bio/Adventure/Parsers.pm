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
    my $contig_count = 0;
  FEAT: while (my $feature = $gff_io->next_feature()) {
        $counters->{total_features}++;
        my $tag = $feature->primary_tag;
        if (defined($gff_type)) {
            next FEAT unless ($feature->primary_tag eq $gff_type);
        }
        my @seqnames = $feature->get_tag_values($gff_tag);
        my $gene_name = $seqnames[0];
        if (!defined($gene_name)) {
            $counters->{unnamed_features}++;
            next FEAT;
        }
        my $str = $feature->strand;
        if (!defined($str)) {
            print $log "This strand is undefined for this feature: ${gene_name}.\n";
            $counters->{unstranded_features}++;
            next FEAT;
        }
        $str = 1 if ($str eq '+1');
        if ($str ne '-1' && $str ne '1') {
            print $log "The strand is neither +1 nor -1 for this feature: ${gene_name}, ${str}.\n";
            $counters->{unstranded_features}++;
            next FEAT;
        }

        my $info = {};
        my $contig = $feature->seq_id;
        $info->{contig} = $contig;
        if (!defined($observed->{$contig})) {
            $observed->{$contig} = [];
            $counters->{contigs}->{$contig} = 0;
            $contig_count = -1;
        }
        $contig_count++;
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
        $info->{next_strand} = $str;
        $info->{gene_name} = $gene_name;

        push(@{$observed->{$contig}}, $info);
        if (defined($observed->{$contig}[$contig_count - 1])) {
            my $prev = $observed->{$contig}[$contig_count - 1];
            my $next_distance = abs($feature->start - $prev->{end});
            $observed->{$contig}[$contig_count - 1]->{next_distance} = $next_distance;
            $observed->{$contig}[$contig_count - 1]->{next_strand} = $str;
            $observed->{$contig}[$contig_count - 1]->{next_start} = $feature->start;
        }
    }
    print "Finished reading $contig_count contigs of gff.\n";
    for my $k (keys %{$observed}) {
        my @features = @{$observed->{$k}};
        my @sorted = sort { $a->{start} <=> $b->{start} } @features;
        $observed->{$k} = \@sorted;
    }
    $gff_in->close();
    return($counters, $observed);
}

1;
