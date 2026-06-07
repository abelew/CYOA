package Bio::Adventure::Software;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use v5.14;
use feature qw"try";
no warnings qw"experimental::try";;
use Capture::Tiny qw":all";
use Cwd qw"abs_path fastcwd";
use Data::Printer;
use Env::Modulecmd;
use Exporter;
use File::Basename qw"basename dirname";
use FileHandle;
use File::Path qw"make_path";
use File::ShareDir qw"dist_file module_dir dist_dir";
use File::Which qw"which";
use Getopt::Long qw"GetOptionsFromArray";
use JSON;
use POSIX qw"ceil";
use String::Interpolate qw"interpolate";
use Template;
my $template_base = dist_dir('Bio-Adventure');
my $template_dir = qq"${template_base}/jobs";

our @EXPORT = qw"Module_Loader Module_Reset Spack_Loader Spack_Reset";
our @ISA = qw"Exporter";

sub Conda_Loader {
    print "I am a stub for working with conda.\n";
}

sub Conda_Reset {
    print "I am a stub for resetting conda environments.\n";
}

=head2 C<Get_Condacmd>

  Figure out what command this computer is using to invoke conda, likely candidates
  include: 'conda', 'miniconda', and a bunch of others which I am not inclined to
  remember right now.  This function should figure out which one(s) the user environment
  currently has working.

=cut
sub Get_Condacmd {
    my $shell_text = q{};
    my $conda_test = which('conda');
    if (!$conda_test) {
        $shell_text .= q!
function conda() {
  echo "The conda command does not appear to exist."
}
export -f conda
!;
    }
    my $conda = {
        found => 0,
        shell_text => $shell_text,
    };
    return($conda);
}

=head2 C<Get_Modulecmd>

  There are also multiple different implementations of the module command.
  This function should find whatever is in the current shell.

=cut
sub Get_Modulecmd {
    my %args = @_;
    my $module_test = qq!/usr/bin/script -q -c 'bash -i -c "type module"' /dev/null |!;
    my $modulecmd_handle = FileHandle->new($module_test);
    my $found = 0;
    my $shell_text = qq!function get_sigterm {
  echo "A SIGTERM was sent to \$JNAME." >> $args{bash_log}
  exit 1
}
export -f get_sigterm
trap get_sigterm SIGTERM

function get_sigerr {
  echo "A ERR was sent to \$JNAME." >> $args{bash_log}
  exit 1
}
export -f get_sigerr
trap get_sigerr ERR

function module() {
!;

    my $num_modulecmd_lines = 0;
  MOD: while (my $line = <$modulecmd_handle>) {
        chomp $line;
        next MOD if ($line =~ /is a function/);
        next MOD if ($line =~ /^module/);
        if ($line =~ /^\{\s*$/) {
            next MOD;
        }
        ## In case there are weird encodings.
        $line =~ s/\s*$//g;
        $num_modulecmd_lines++;
        $shell_text .= qq"${line}\n";
    }
    $shell_text .= q"export -f module

";
    $modulecmd_handle->close();
    if (-f 'typescript') {
        rm('typescript');
    }
    ## If module is not in the current environment, check for the
    ## modulecmd executable and just use that.
    if ($num_modulecmd_lines == 0) {
        my $modulecmd = which('modulecmd');
        if ($modulecmd) {
            $found++;
            $shell_text .= qq!${modulecmd} bash "$@"
export -f module

!;
        } else {
            $shell_text .= q!echo "The module command does not appear to exist".
}
export -f module

!;
        }
    } else {
        $found++;
    }
    my $module = {
        found => $found,
        shell_text => $shell_text,
    };
    return($module);
}

=head2 C<Get_Spackarch>

  Discover the string spack uses to define the current cpu architecture.

=cut
sub Get_Spackarch {
    my $arch_string = qx"spack arch";
    if ($arch_string =~ /^The spack command does not appear/) {
        return(undef);
    }
    return($arch_string);
}

=head2 C<Get_Spackcmd>

  At this point, I sort of just assume spack will have multiple different invocations/versions.
  So, though for now there is only one, figure it out.

=cut
sub Get_Spackcmd {
    my %args = @_;
    my $shell_text = q{};
    my $spack_test = qq!/usr/bin/script -q -c 'bash -i -c "type spack"' /dev/null |!;
    my $spack_handle = FileHandle->new($spack_test);
    my $found = 0;
    $shell_text .= qq!function spack() {
!;
    my $num_spack_lines = 0;
  SP: while (my $line = <$spack_handle>) {
        chomp $line;
        next SP if ($line =~ /is a function/);
        next SP if ($line =~ /^spack/);
        next SP if ($line =~ /not found/);
        next SP if ($line =~ /^\{\s*$/);
        ## In case there are weird encodings.
        $line =~ s/\s*$//g;
        $num_spack_lines++;
        $shell_text .= qq"${line}\n";
    }
    $shell_text .= q"export -f spack

";
    $spack_handle->close();
    if (-f 'typescript') {
        rm('typescript');
    }
    if ($num_spack_lines == 0) {
        $shell_text .= q!echo "The spack command does not appear to exist".
}
export -f spack

!;
    } else {
        print "Software::Get_Spack: There are: $num_spack_lines \n";
        $found = 1;
    }
    my $spack_cmd = {
        found => $found,
        shell_text => $shell_text,
    };
    return($spack_cmd);
}

=head2 C<Module_Catalog>

  Create a data structure containing all of the available modules, their versions, and
  which one(s) are the default.  We may want to have some logic to exclude spack-derived
  modules and use the Spack_Catalog to find them.

=cut
sub Module_Catalog {
    my %args = @_;
    my $skip_spack = 1;
    if (defined($args{skip_spack})) {
        $skip_spack = $args{skip_spack};
    }
    my $module_avail = qq!bash -i -c "module avail -j" 2>&1 |!;
    my $module_cmd = FileHandle->new($module_avail);
    my $module_string = q{};
    while (my $line = <$module_cmd>) {
        $module_string .= $line;
    }
    $module_cmd->close();
    my $json = JSON->new->allow_nonref;
    $json = $json->relaxed([1]);
    my $modules_ref = $json->decode($module_string);
    ## This results in a hash of hashes of hashes, where the outer key is the module directory,
    ## the middle key is the package_name/package_version and the inner keys include:
    ## name, pathname, symbols (default), tags, type, and via.
    ## I am going therefore to rekey this so that it is: (outer->inner)
    ## name -> version -> { default }
    ## with the option to add some more stuff later.
    my $pkg_info = {};
  DIRS: for my $module_dir (keys %{$modules_ref}) {
        if ($skip_spack) {
            next DIRS if ($module_dir =~ /\/spack\//);
        }
        my $module_dir_ref = {};
        my %module_dir_info = %{$modules_ref->{$module_dir}};
        my $found_default = 0;
        my $previous_name = q{};
        my $counter = 0;
        my $version = q{};
        my $name = q{};
        for my $pkg (sort keys %module_dir_info) {
            $counter++;
            ($name, $version) = split(/\//, $pkg, 2);
            if (!defined($pkg_info->{$name})) {
                $pkg_info->{$name} = {};
                ## If we have moved to a new package from an old one, and we never found
                ## a default version, then the last version observed is the default.
                if (!$found_default && $counter > 0) {
                    $pkg_info->{$previous_name}->{$version}->{default} = 1;
                } else {
                    $found_default = 0;
                }
            }
            $pkg_info->{$name}->{$version} = {
                mod_dir => $module_dir,
                default => 0,
            };
            if (defined($module_dir_info{$pkg}->{symbols}->[0])) {
                print "HEY: Symbols are defined!\n";
                p $module_dir_info{$pkg}->{symbols}->[0];
                $pkg_info->{$name}->{$version}->{default} = 1;
                $found_default = 1;
            }
            $previous_name = $name;
        } ## End iterating over packages in an individual module directory.
        if (!$found_default) {
            $pkg_info->{$previous_name}->{$version}->{default} = 1;
        }
    }
    p $pkg_info;
    return($pkg_info);
}

sub Module_Catalog_bad {
    my $modules_ref = {};
    my $module_dir = q{};
    my $module_dir_ref = {};
    my $avail_test = qq!bash -i -c "module avail --no-pager --redirect --color=never -t -C" |!;
    my $avail_handle = FileHandle->new($avail_test);
    my $previous_module = q{};
    my $found_default = 0;
    my $counter = 0;
    my $default = 0;
  MOD: while (my $line = <$avail_handle>) {
        chomp $line;
        $line =~ s/\s+$//g;
        $counter++;
        $default = 0;
        my $version = q{};
        my $name = q{};
        if ($line =~ /:$/) {
            if ($counter > 1) {
                $modules_ref->{$module_dir} = $module_dir_ref;
                print "Adding a new hash to modules\n";
            }
            $module_dir = $line;
            $module_dir =~ s/:$//;
            $module_dir_ref = {};
            $previous_module = q{};
            $found_default = 0;
        } else {
            $line =~ s/^(\S+)\s+.*$/$1/g;
            ($name, $version) = split(/\//, $line, 2);
            next MOD if ($name eq '');
            if ($name ne $previous_module) {
                if ($found_default) {
                    $found_default = 0;
                } else {
                    ## If we switched to a new module, and never found a default
                    ## for the previous, then the last one printed is the default.
                    $module_dir_ref->{$previous_module}->{$version} = 1;
                }
            }
            if ($version =~ /\(default\)$/) {
                $version =~ s/^(.+)\(.*$/$1/g;
                $default = 1;
                $found_default = 1;
            }
            $module_dir_ref->{$name}->{$version} = $default;
            $previous_module = $name;
        } ## End checking for a module directory.
    }
    $modules_ref->{$module_dir} = $module_dir_ref;
    return($modules_ref);
}

=head2 C<Module_Loader>

This loads environment-modules into the current ENV.

=cut
sub Module_Loader {
    my ($class, %args) = @_;
    my $action = 'load';
    if ($args{action}) {
        $action = $args{action};
    }

    ## Start with the current ENV and an empty set of things to load
    ## and test.
    my $module_string = q"";
    my %old_env = %ENV;
    my @mod_lst = ();
    my $mod_class = ref($args{modules});
    my @test_lst = ();

    my $test_class = ref($args{exe});
    if ($mod_class eq 'SCALAR') {
        push(@mod_lst, $args{modules});
    } elsif ($mod_class eq 'ARRAY') {
        @mod_lst = @{$args{modules}};
    } elsif ($mod_class eq 'HASH') {
        my %mods = %{$args{modules}};
        @mod_lst = keys %mods;
        for my $k (@mod_lst) {
            push(@test_lst, $mods{$k});
        }
    } elsif (!$mod_class) {
        push(@mod_lst, $args{modules});
    } else {
        print "I do not know this class: $mod_class, $args{modules}\n";
    }

    if ($test_class eq 'SCALAR') {
        push(@test_lst, $args{exe});
    } elsif ($test_class eq 'ARRAY') {
        for my $e (@{$args{exe}}) {
            push(@test_lst, $e);
        }
    } elsif (!$test_class) {
        push(@test_lst, $args{exe});
    } else {
        print "I do not know this class: $test_class, $args{exe}\n";
    }

    my $count = 0;
    my $mod;
    if ($action eq 'unload') {
        @mod_lst = reverse(@mod_lst);
    }
    my @messages = ();
    my @loads = ();
    my $loaded = 0;
  LOADER: for $mod (@mod_lst) {
        ## Skip one's self.
        my ($stdout, $stderr, @result) = capture {
            try {
                my $test;
                push(@loads, $mod);
                if ($action eq 'load') {
                    $test = Env::Modulecmd::load($mod);
                    $loaded++;
                    if ($loaded == 1) {
                        $module_string = q'module add ';
                    }
                    $module_string .= qq"$mod ";
                } else {
                    $test = Env::Modulecmd::unload($mod);
                }
            }
            catch ($e) {
                ## print "There was an error loading ${mod}, ${e}\n";
                my $message = qq"There was an error loading ${mod}: ${e}";
                push(@messages, $message);
            }
            ;
        };
        $count++;
    } ## End iterating over every mod in the list

    ## If we got a set of test programs, check that they are in the PATH:
    if (scalar(@test_lst) > 0) {
      EXELOOP: for my $exe (@test_lst) {
            next EXELOOP unless($exe);
            my $check = which($exe);
            die(qq"Could not find ${exe} in the PATH.") unless ($check);
        }
    }
    my $set = $class->modules(\@mod_lst);
    ## Note, this module_string does not appear to be getting added to the function-local state.
    my $state_ref = $class->current_state;
    $state_ref->{module_string} = $module_string;
    $state_ref->{spack_string} = q{};
    $state_ref->{conda_string} = q{};
    my $state_set = $class->current_state($state_ref);
    return(\%old_env);
}

sub Module_Reset {
    my ($class, %args) = @_;
    my $orig_env = $args{env};
    for my $k (keys %{$orig_env}) {
        $ENV{$k} = $orig_env->{$k};
    }
}

=head2 C<Spack_Catalog>

  Create a catalog of the available spack versions/hashes/etc

=cut
sub Spack_Catalog {
    my $spack_by_name = {};
    my $spack_string = qx"spack find --json";
    my $json = JSON->new->allow_nonref;
    $json = $json->relaxed([1]);
    my $spack_ref = $json->decode($spack_string);
    ## Each element of spack_ref is a hash containing a wealth of useful stuff.
    ## I want to rekey this information so that I can easily get all
    ## installed spack specs which match a tool I want.
    for my $spec (@{$spack_ref}) {
        my $current_name = $spec->{name};
        if (defined($spack_by_name->{$current_name})) {
            my @current_specs = @{$spack_by_name->{$current_name}};
            push(@current_specs, $spec);
            $spack_by_name->{$current_name} = \@current_specs;
        } else {
            $spack_by_name->{$current_name} = [ $spec ];
        }
    }
    return($spack_by_name);
}

sub Spack_Loader {
    print "I am a stub for loading spack packages/envs.\n";
}

sub Spack_Reset {
    print "I am a stub for resetting spack packages/envs.\n";
}

1;
