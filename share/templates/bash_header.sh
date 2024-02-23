#!/usr/bin/env bash
[% IF cluster == 'slurm' %]
  [%- INCLUDE slurm_bash_header.sh -%]
[% END %]
cd [% basedir %]
set -o errexit
set -o errtrace
set -o pipefail
export LESS='--buffers 0 -B'
[% module_string %]
[% conda_string %]
