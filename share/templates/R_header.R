#!/usr/bin/env Rscript
[% IF cluster == 'slurm' %]
  [%- INCLUDE slurm_bash_header.sh -%]
[% END %]

