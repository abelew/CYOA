package Bio::Adventure::Defaults;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Which qw"which";
