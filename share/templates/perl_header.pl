#!/usr/bin/env perl
[% IF cluster == 'slurm' %]
  [%- INCLUDE slurm_perl_header.pl -%]
[% END %]

use strict;
use FileHandle;
use Bio::Adventure;
my $out = FileHandle->new(">>outputs/log.txt");
my $d = qx'date';
chdir("[% workdir %]");
my $h = Bio::Adventure->new();

[% IF pdata %]
if (-r "[% pdata %]") {
  ## Pull options from the option file: outputs/01trimomatic/kraken2mtrxoULZ.pdata
  my $loaded = $h->Load_Vars(input => '[% pdata %]');
}
[% END %]
