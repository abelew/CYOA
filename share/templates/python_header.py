#!/usr/bin/env python
[% IF cluster == 'slurm' %]
  [%- INCLUDE slurm_bash_header.sh -%]
[% END %]

