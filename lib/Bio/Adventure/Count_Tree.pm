package Bio::Adventure::Count_Tree;
use strict;
use Tree::DAG_Node;
our @ISA = qw"Tree::DAG_Node";

=head1 NAME

 Bio::Adventure::Count_Tree - Use Tree::DAG_Node to count up heirarchical data

 This is taken almost verbatim from: https://www.perlmonks.org/?node_id=153259
 It seeks to provide an easy way to count up the nodes of a gff file by type.

=head2 C<new>

 The constructor uses the @our ISA above to create a new instance of
 Tree::DAG_Node and set whatever options we want with it.
 Given that starting point, the functions which follow will be used to
 do the actual counting and ensuring that we get names which make sense.

=cut
sub new {
    my ($class, $options) = @_;
    my $super = $class->SUPER::new();
    my $self = bless $super;
    $self->attributes($options);
    return $self;
}

=head2 C<ids>

 This extracts the IDs from the nodes in our tree.

=cut
sub ids {
    my ($node, $val) = @_;
    if ($val) {
        my @ids = @{$node->attributes->{ids}};
        push(@ids, $val);
        $node->attributes->{ids} = \@ids;
    } else {
        return $node->attributes->{ids};
    }
}

=head2 C<num_ids>

 This counts up IDs by type in our tree.

=cut
sub num_ids {
    my ($node, $val) = @_;
    if ($val) {
        $node->attributes->{num_ids}++;
    } else {
        return $node->attributes->{num_ids};
    }
}

=head2 C<print_num_ids>

 This prints out the IDs at each level of the tree.

=cut
sub print_num_ids {
    $_[0]->walk_down({
        callback=> sub {
            my $node = shift;
            printf "%s%.7s\t NumIDs: %d \n",
                "  " x $_[0]->{_depth},
                ## $node->name, scalar(@{$node->attributes->{ids}}),
                $node->name, $node->attributes->{num_ids},
            }, _depth => 0 });
}

=head2 C<by_name>

 This performs the walk by name.

=cut
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

1;
