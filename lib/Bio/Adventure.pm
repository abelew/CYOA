package Bio::Adventure;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
use vars qw"$VERSION";
use feature qw"try";
no warnings qw"experimental::try";
use Exporter;
our @EXPORT = qw"Submit Wait Extra_Options";
our @ISA = qw"Exporter";

use AppConfig qw":argcount :expand";
use Bio::SeqIO;
use Bio::DB::Universal;
use Bio::Root::RootI;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Capture::Tiny qw":all";
use Carp qw"croak carp confess cluck longmess";
use Cwd qw"abs_path getcwd cwd";
use Data::Printer;
use Env qw"COMPRESSION XZ_OPTS XZ_DEFAULTS HOME";
use File::Basename;
use File::Find;
use File::Spec;
use File::Which qw"which";
use File::Path qw"make_path remove_tree";
use File::Temp;
use FileHandle;
use Getopt::Long qw"GetOptionsFromArray";
use IO::String;
use List::MoreUtils qw"uniq";
use List::Util qw"any";
use Log::Log4perl;
use PerlIO;
use Pod::Usage;
use Signal::StackTrace qw"TTIN";
use Storable qw"lock_store lock_retrieve";
use String::Interpolate qw"interpolate";
use Term::ReadLine;
use Term::UI;

use Bio::Adventure::Config;
use Bio::Adventure::Languages;
use Bio::Adventure::Software;
use Bio::Adventure::Slurm;
use Bio::Adventure::Bash;

use Bio::Adventure::Align;
use Bio::Adventure::Align_Blast;
use Bio::Adventure::Align_Fasta;
use Bio::Adventure::Annotation;
use Bio::Adventure::Annotation_Genbank;
use Bio::Adventure::Assembly;
use Bio::Adventure::Cleanup;
use Bio::Adventure::Count;
use Bio::Adventure::Compress;
use Bio::Adventure::Convert;
use Bio::Adventure::Feature_Prediction;
use Bio::Adventure::Index;
use Bio::Adventure::Map;
use Bio::Adventure::Metadata;
use Bio::Adventure::Parsers;
use Bio::Adventure::Phage;
use Bio::Adventure::Phylogeny;
use Bio::Adventure::PopGen;
use Bio::Adventure::Prepare;
use Bio::Adventure::Resistance;
use Bio::Adventure::Riboseq;
use Bio::Adventure::QA;
use Bio::Adventure::SeqMisc;
use Bio::Adventure::Splicing;
use Bio::Adventure::Structure;
use Bio::Adventure::SNP;
use Bio::Adventure::TNSeq;
use Bio::Adventure::Torque;
use Bio::Adventure::Trim;
use Bio::Adventure::Visualization;
use Bio::Adventure::Pipeline;

## Available software packages.
has available_sw => (is => 'rw', default => sub{ {} }, writer => '_set_available_sw');
## The current function.
has caller => (is => 'rw', default => undef);
## Arguments provided by the commandline.
has cmdline_args => (is => 'ro', writer => '_set_cmdline_args', default => sub { {} });
## Keep a reference to the input information for the current function here.
has current_inputs => (is => 'rw', default => sub{ {} });
## The current merged state of the various arguments.
has current_state => (is => 'rw', default => sub{ {} });
## Args passed from function to function.
has function_args => (is => 'rw', default => sub { {} });
## Default file from which to acquire function-specfic defaults
has function_file => (is => 'ro', writer => '_set_function_file', default => 'functions.json');
## Arguments provided in the functions.json.
has function_defaults => (is => 'ro', writer => '_set_function_defaults', default => sub { {} });
## Default file from which to acquire global defaults.
## If the file is not found, check in the dist-share directory.
has function_paths => (is => 'rw', default => sub{ {} });
has global_file => (is => 'ro', writer => '_set_global_file', default => 'defaults.json');
## Global default arguments.
has global_defaults => (is => 'ro', writer => '_set_global_defaults', default => sub { {} });
## Arguments provided by an ini file.
has global_paths => (is => 'ro', writer => '_set_global_paths', default => sub{ {} });
has ini_file => (is => 'ro', writer => '_set_ini_file', default => 'cyoa.conf');
has ini_args => (is => 'ro', writer => '_set_ini_args', default => sub{ {} });
has instance_args => (is => 'rw', default => sub{ {} });
## Information about the inputs used over time by cyoa.  I am changing
## this to an array ref for which every element is the result
## of Get_Input_Info (which is also an array ref); but this way
## we can track all inputs from beginning->end of a pipeline.
has inputs => (is => 'rw', default => sub{ [] });
## Information about each submitted job.
has jobs => (is => 'rw', default => sub{ [] });
## Array of job IDs
has job_ids => (is => 'rw', default => sub{ [] });
## Arbitrary json file containing arguments
has json_file => (is => 'ro', writer => '_set_json_file', default => undef);
## Arguments provided by a json file.
has json_args => (is => 'rw', writer => '_set_json_args', default => sub{ {} });
## A menu of possible adventures.
has menu => (is => 'ro', writer => '_set_menu', default => sub{ {} });
## The set of currently loaded modules.
has modules => (is => 'rw', default => sub{ [] });
## Similar to inputs, have options be a list of the option state
## during each invoked job.
has options => (is => 'rw', default => sub{ [] });
## Information about the current outputs.
has outputs => (is => 'rw', default => sub{ {} });
## json file containing path information
has path_file => (is => 'rw', default => 'paths.json');
## Useful paths for the current job.
has paths => (is => 'rw', default => sub{ {} });
## Set of things to do.
has todos => (is => 'ro', writer => '_set_todos', default => sub{ []} );

=head1 NAME

    Bio::Adventure - A perl library to make preprocessing high-throughput data easier.

=head1 SYNOPSIS

    Bio::Adventure.pm  Methods to simplify creating job scripts and submitting
 them to a
    computing cluster.  Bio::Adventure makes heavy use of command line options and as a result
    accepts many options.  The simplest way to access these methods is via the cyoa script, but
    everything may be called directly in perl.

    use Bio::Adventure;
    my $hpgl = Bio::Adventure->new(cluster => 0);  ## Force it to run locally.
    $hpgl->{species} = 'mmusculus';  ## There had better be a mmusculus.[fasta&gff] in $hpgl->{libdir}
    my $hpgl = Bio::Adventure->new(species => 'scerevisiae', libdir => '/home/bob/libraries');
    $hpgl->{input} = 'test.fastq';  ## Also via with cyoa --input=test.fastq

    ## Run Fastqc on untrimmed sequence.
    my $bp = $hpgl->Fastqc()
    ## Generates test-trimmed.fastq from test.fastq
    my $trim = $hpgl->Trimomatic();

    Conversely, just run 'cyoa' from the command line and follow the prompts and/or provide
    all the parameters with the appropriate command line options:

    > cyoa --task map --method hisat --species hg38_100 --input r1-trimmed.fastq.xz:r2-trimmed.fastq.xz
    > cyoa --task align --method blast --library nr --input transcripts.fasta

=head1 DESCRIPTION

    This library should write out PBS/slurm compatible job files and
    submit them to the appropriate computing cluster.  It should also
    collect the outputs and clean up the mess.

=head2 Options

  Most options for cyoa are defined in Adventure.pm, there are a small number are used
  only for specific functions defined within them.  The document L<Bio::Adventure/Arguments.pm>
  will seek to explain how to access them and their various purposes.
  the most likely purpose for each parameter follows:

=cut
our $AUTOLOAD;
##our @EXPORT_OK = qw"";
$VERSION = '202601';
$COMPRESSION = 'xz';
$XZ_OPTS = '-9e';
$XZ_DEFAULTS = '-9e';
$ENV{PERL_USE_UNSAFE_INC} = 0;
if (!defined($ENV{TMPDIR})) {
    $ENV{TMPDIR} = '/tmp';
}
## Added to test Bio:DB::SeqFeature::Store
if (!defined($ENV{PERL_INLINE_DIRECTORY})) {
    my $filename = File::Temp::tempnam($ENV{TMPDIR}, "inline_$ENV{USER}");
    $ENV{PERL_INLINE_DIRECTORY} = $filename;
    my $made = make_path($filename);
}

sub Append_Inputs {
    my ($class, $input) = @_;
    my @current = @{$class->inputs};
    push(@current, $input);
    my $current_num = scalar(@current);
    $class->inputs(\@current);
    return($current_num);
}

sub Get_Last_Input {
    my ($class, %args) = @_;
    my @in = @{$class->inputs};
    my $last = $in[$#in];
    $class->{current_state}->{input} = $last;
    return($last);
}

sub Extra_Options {
    my ($class, %args) = @_;
    my %opt = %{$args{options}};
    my $extras = $args{extras};
    for my $k (keys %{$extras}) {
        $opt{$k} = $extras->{$k};
    }
    return(%opt);
}

=head2 C<Help>

    Help() always gives 0.
    Before it returns, it will hopefully print some useful information
    regarding ways to invoke Bio::Adventure.pm.

=cut
sub Help {
    my ($class, $args) = @_;
    my $fh = \*STDOUT;
    ##    my $usage = pod2usage(-output => $fh, -verbose => 99, -sections => "SYNOPSIS");
    use Pod::Find qw"pod_where";
    pod2usage(-input => pod_where({-inc => 1}, __PACKAGE__),
              -verbose => 2, -output => $fh,
              -sections => "NAME|SYNOPSIS|DESCRIPTION|VERSION",
              -exitval => 'NOEXIT');
    return(0);
}

=head2 C<BUILD>

This is the main CYOA constructor.  It takes up the variables at the
top of Adventure.pm along with any arguments pulled in from @ARGV and
incorporates them into the class.

Most of the logic in it is working toward setting up GetOpt::Long
options.  Thus, if there is confusion about what is happening, that is
the first place to look.

=cut
sub BUILD {
    my ($class, $args) = @_;
    ## I made some pretty big changes to the logic of cyoa, I now want
    ## to gather all of my options/arguments/etc into a few hash/array
    ## references and pass them about from function to function rather
    ## than have the morass of individual variables.

    ## I therefore created the file defaults.json in share/ from which
    ## I would like to acquire the global defaults.  I will then
    ## gather function-specific parameters from functions.json, etc.

    ## Read the global defaults json filend extract variables from it.
    if (!%{$class->global_defaults}) {
        my $defaults = Bio::Adventure::Config::Get_JSON(
            input_file => $class->global_file,);
        my $set = $class->_set_global_defaults($defaults);
    }

    ## Pull anything from ARGV into the cmdline_args slot.
    if (!%{$class->cmdline_args}) {
        my $cmdline_args = Bio::Adventure::Config::Get_ARGV();
        my $set = $class->_set_cmdline_args($cmdline_args);
    }

    ## Read the json file of function-specific default values.
    ## Ideally I should be pulling quite a lot out of the global
    ## defaults into this.
    if (!%{$class->function_defaults}) {
        ## print "The function-specific defaults are undefined.\n";
        my $defaults = Bio::Adventure::Config::Get_JSON(
            input_file => $class->function_file);
        ## Fill in some elements of the defaults if they are not provided in the json file.
        ## Notably, each function should have a stderr/stdout/main/module/jname/template.
        ## TODO: I want a function which collects all software packages from spack/conda/module
        ## and will provide them by name.
        for my $fn (keys %{$defaults}) {
            my $lc_fn = lc($fn);
            my %inner = %{$defaults->{$fn}};
            if (!defined($inner{stderr})) {
                $defaults->{$fn}->{stderr} = q"${output_dir}/stderr";
            }
            if (!defined($inner{stdout})) {
                $defaults->{$fn}->{stdout} = q"${output_dir}/stdout";
            }
            if (!defined($inner{main})) {
                $defaults->{$fn}->{main} = q"${output_dir}/main";
            }
            if (!defined($inner{jname})) {
                $defaults->{$fn}->{jname} = $lc_fn;
            }
            if (!defined($inner{modules})) {
                $defaults->{$fn}->{modules} = [ $lc_fn ];
            }
            if (!defined($inner{template})) {
                $defaults->{$fn}->{template} = qq"${lc_fn}.tt";
            }
            if (!defined($inner{jprefix})) {
                $defaults->{$fn}->{jprefix} = '00';
            }
        } ## Finished checking expected keys of each function.
        my $set = $class->_set_function_defaults($defaults);
    }

    ## If there is an ini file in the current-working-tree
    ## Use variables from it.
    if (!%{$class->ini_args}) {
        ## print "Arguments from an ini file are undefined.\n";
        my $ini = Bio::Adventure::Config::Get_Ini(
            input_file => $class->ini_file);
        my $set = $class->_set_ini_args($ini);
    }

    ## Similarly to the ini file, pull json options.
    if (!%{$class->json_args}) {
        my $json = Bio::Adventure::Config::Get_JSON(
            input_file => $class->json_file);
        my $set = $class->_set_json_args($json);
    }

    ## In theory, this should give a hash of arguments passed to new()
    ## It needs testing.
    my $set = $class->instance_args($args);
    ## Create the current state of cyoa from the various
    ## arguments/ini/etc above.
    my $current_state = Bio::Adventure::Config::Merge_Args(
        instance => $class->instance_args,
        global => $class->global_defaults,
        cmdline => $class->cmdline_args,
        ini => $class->ini_args,
        json => $class->json_args);

    ## Coerce a few values from strings here
    my @suffixes = ('.gz', '.xz', '.bz2');
    if (defined($current_state->{suffixes})) {
        @suffixes = split(/[;:,]/, $current_state->{suffixes});
    }
    $current_state->{suffix_array} = \@suffixes;
    my $path_agrees = Check_Libpath(libdir => $current_state->{libdir},
                                    libpath => $current_state->{libpath});
    $current_state->{libpath} = $path_agrees->{libpath};
    $current_state->{libdir} = $path_agrees->{libdir};
    ## Write out a few functions to assist the shell: error traps and module particularly.
    my $shell_functions = Bio::Adventure::Config::Get_Shell_Functions(
        bash_log => $current_state->{bash_log});
    my $exists = make_path(dirname($current_state->{bash_functions}));
    if (!-f $current_state->{bash_functions}) {
        my $source_fh = FileHandle->new(">$current_state->{bash_functions}");
        print $source_fh $shell_functions->{shell_text};
        $source_fh->close();
    }
    $current_state->{bash_functions} = abs_path($current_state->{bash_functions});
    ## FIXME: There is a better place to set the absolute path of the basedir
    $current_state->{abs_basedir} = abs_path($current_state->{basedir});
    $current_state->{bash_source} = qq"source $current_state->{bash_functions}";
    if ($current_state->{modules}) {
        $current_state->{module_string} = qq"module add $current_state->{modules}";
    }

    if (!%{$class->available_sw}) {
        my $sw_module = {};
        if ($shell_functions->{module_found}) {
            $sw_module = Bio::Adventure::Software::Module_Catalog();
        }
        my $sw_spack = {};
        if ($shell_functions->{spack_found}) {
            $sw_spack = Bio::Adventure::Software::Spack_Catalog();
        }
        my $sw_conda = {};
        if ($shell_functions->{conda_found}) {
            $sw_conda = Bio::Adventure::Software::Conda_Catalog();
        }
        my $sw_container = {};
        my $sw = {
            module => $sw_module,
            spack => $sw_spack,
            conda => $sw_conda,
            containers => $sw_container,
        };
        my $set = $class->_set_available_sw($sw);
    }

    $class->current_state($current_state);
    ## Get the relatively small number of paths which should prove universal.
    my $global_paths = $class->Bio::Adventure::Config::Get_Global_Paths();

    my @available_clusters = Get_Clusters();
    $current_state->{available_clusters} = \@available_clusters;
    ## Now that we have collected all the various potential options, set the state with them.
    $class->current_state($current_state);
    ## Add the first element to inputs
    my @inputs = @{$class->inputs};
    if ($current_state->{input}) {
        my $first = $class->Bio::Adventure::Config::Get_Input_Info(
            $current_state->{input});
        push(@inputs, $first);
    }
    else {
        ## If no input was provided, send a blank.
        push(@inputs, []);
    }
    $class->inputs(\@inputs);

    ## If it turns out we want a tui menu, get them here.
    my $menus = {};
    if ($current_state->{interactive}) {
        $menus = Bio::Adventure::Config::Get_Menus();
    }
    $class->_set_menu($menus);

    ## Finally, use the command line options to figure out what tasks should be performed
    my $todos = Bio::Adventure::Config::Get_TODOs(%{$class->cmdline_args});
    $class->_set_todos($todos);
    return($args);
}

sub Add_Config_File {
    my ($class, %args) = @_;
    my $defaults = $args{defaults};
    my $steps = $args{steps};
    my $input = $args{input};
    my $config = AppConfig->new({ CREATE => '.*' ,
                                  GLOBAL => { ARGCOUNT => AppConfig::ARGCOUNT_ONE } });
    my $info = $config->file($input);
    my $variable = q{};
    my $step = '0';
    my %data = $config->varlist(".*");
    for my $config_option (keys %data) {
        my $step = 'defaults';
        $variable = $config_option;
        my $config_value = $data{$config_option};
        if ($config_option =~ /^defaults_/) {
            $variable =~ s/^defaults_//g;
        }
        if ($config_option =~ /^(\d+)_(.+)$/) {
            $step = $1;
            $variable = $2;
            if (!defined($steps->{$step})) {
                $steps->{$step} = {};
            }
            $steps->{$step}->{$variable} = $config_value;
        }
        else {
            $defaults->{$variable} = $config_value;
        }
    }
    return($defaults, $steps);
}

=head2 C<Check_Input>

Given a set of inputs, do a little checking to try to ensure that they
are usable.

This function is not used in many places and may need culling.  Either
that or it should be improved a little and propagated.

=cut
sub Check_Input {
    my ($class, %args) = @_;
    my $file_list;
    if (ref($args{files}) eq 'SCALAR' || ref($args{files}) eq q{}) {
        if ($args{files} =~ /[:,;]/) {
            $args{files} =~ s/[:,;]$//g;
            my @tmp = split(/[:,;]/, $args{files});
            $file_list = \@tmp;
        }
        else {
            $file_list->[0] = $args{files};
        }
    }
    elsif (ref($args{files}) eq 'ARRAY') {
        $file_list = $args{files};
    }
    else {
        my $unknown_class = ref($args{files});
        warn("I do not know type: ${unknown_class}.");
    }
    my $found = {};
    foreach my $file (@{$file_list}) {
        $found->{$file} = 0;
        my $first_test = $file;
        my $second_test = basename($file, $class->{suffixes});
        my $third_test = basename($second_test, $class->{suffixes});
        $found->{$file} = $found->{$file} + 1 if (-r $first_test);
        $found->{$file} = $found->{$file} + 1 if (-r $second_test);
        $found->{$file} = $found->{$file} + 1 if (-r $third_test);
    }
    for my $f (keys %{$found}) {
        if ($found->{$f} == 0) {
            die("Unable to find a file corresponding to $f.");
        }
    }
    return($found);
}

sub Check_Job {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args);
    my $result;
    if ($options->{cluster} eq 'slurm') {
        $result = $class->Bio::Adventure::Slurm::Check_Job(%args);
    }
    else {
        print "This runner doesn't yet have a Check_Job.\n";
    }
    return($result);
}

sub Check_Libpath {
    my %args = @_;
    my $provided_libdir = 0;
    if (defined($args{libdir})) {
        $provided_libdir = 1;
    }
    else {
        $args{libdir} = '\$HOME/libraries';
    }
    my $provided_libpath = 0;
    if (defined($args{libpath})) {
        $provided_libpath = 1;
    }
    else {
        $args{libpath} = "$ENV{HOME}/libraries";
    }
    ## Make sure that the libdir and libpath agree with one another.
    ## I am intending to use libdir often(always?) as the shell variable $HOME.
    ## But this may change from host to host and I want to be able to
    ## perform operations on stuff in that directory while still having some
    ## degree of flexibility.
    ## Thus, this will check if libdir starts with a \$, if so, extract that variable
    ## from the environment, and check if the libpath agrees with it.

    my $interpolated_libdir = $args{libdir};
    if ($args{libdir} =~ /\$/) {
        if ($interpolated_libdir =~ /^(.*)?\$\{*([A-Z]+)\}*(.*)$/) {
            my $prefix = $1;
            if (!defined($prefix)) {
                $prefix = q{};
            }
            my $varname = $2;
            if (!defined($varname)) {
                $varname = q{};
            }
            my $suffix = $3;
            if (!defined($suffix)) {
                $suffix = q{};
            }
            $interpolated_libdir = qq"${prefix}$ENV{$varname}${suffix}";
        }
        else {
            $interpolated_libdir = $args{libdir};
        }
    }

    my $ret = {
        libdir => $args{libdir},
        libpath => $interpolated_libdir,
    };

    if ($provided_libdir && $provided_libpath) {
        ## Both were passed to the constructure, check that they agree, if not
        ## Then warn and go with libdir.
    }

    return($ret);
}

sub Fa_to_Hash {
    my ($class, %args) = @_;
    my $in = Bio::Adventure::Get_FH(input => $args{input});
    my $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $in);
    my $key = $args{key};
    my $output = {};
  SEQ: while (my $seq = $seqio->next_seq) {
        my $seqid;
        if (defined($key)) {
            $seqid = $seq->desc;
        }
        else {
            $seqid = $seq->id;
        }
        my $short_id = $seqid;
        if ($seqid =~ m/${key}\s*[:|=](\S+)/) {
            $short_id = $1;
            $short_id =~ s/\.\d{1,2}$//; ## Get rid of terminal ID numbers (assume it doesn't pass 99)
        }
        my $sequence = $seq->seq;
        $output->{$short_id} = $seq;
    }
    $in->close();
    return($output);
}

sub Get_Clusters {
    my %checks = (
        bash => ['bash'],
        slurm => ['sacctmgr', 'srun', 'sbatch'],
        torque => ['qsub'],
        perl => ['perl'],
    );
    my @available = ();
  TYPE: for my $check (keys %checks) {
        my @tests = @{$checks{$check}};
        my $avail = 1;
      EXE: for my $exe (@tests) {
            my $check = which($exe);
            if (!defined($check)) {
                $avail = 0;
                last EXE;
            }
        }
        if ($avail) {
            push(@available, $check);
        }
    }
    return(@available);
}

=head2 C<Get_FC>

  Returns an arbitrary command used to send a file to STDOUT.  Note that if a program requires
  O_SEEK, then this will not work.

=cut
sub Get_FC {
    my ($class, %args);
    if ((scalar(@_) % 2) == 1) {
        ($class, %args) = @_;
    }
    else {
        %args = @_;
    }
    if (defined($args{needs_seek})) {
        my $opener = $class->Get_FC_Seekable(%args);
        return($opener);
    }
    my $ignored = q{};
    $ignored = $args{ignored} if ($args{ignored});
    my $opener = $args{input};
    my $input_dir = dirname($args{input});
    if ($args{input} =~ /\.xz$/ && $ignored !~ /xz/) {
        $opener = qq"xz -dc $args{input}";
    }
    elsif ($args{input} =~ /\.bz2$/ && $ignored !~ /bz2/) {
        $opener = qq"bunzip2 -c $args{input}";
    }
    elsif ($args{input} =~ /\.gz$/ && $ignored !~ /gz/) {
        $opener = qq"gunzip -c $args{input}";
    }
    else {
        $opener = qq"less $args{input}";
    }
    return($opener);
}

=head2 C<Get_FC_Seekable>

  Returns an arbitrary command used to send a file to STDOUT.  Use this if the downstream
  command requires the ability to SEEK the input; this will provide a command to preprocess
  the file (decompress) to its seekable format, a new filename resulting from that process, and
  a command to cleanup the resulting mess.

=cut
sub Get_FC_Seekable {
    my %args = @_;
    my $opener = $args{input};
    my $input_dir = dirname($args{input});
    if ($args{input} =~ /\.xz$/) {
        my $new_input = qq"${input_dir}/" . basename($args{input}, ('.xz'));
        $opener = {
            decompress => qq"xz -d $args{input}",
            input => $new_input,
            cleanup => qq"rm $args{input}"
        };
    }
    elsif ($args{input} =~ /\.bz2$/) {
        my $new_input = qq"${input_dir}/" . basename($args{input}, ('.bz2'));
        $opener = {
            decompress => qq"bunzip2 $args{input}",
            input => $new_input,
            cleanup => qq"rm $args{input}"
        };
    }
    elsif ($args{input} =~ /\.gz$/) {
        my $new_input = qq"${input_dir}/" . basename($args{input}, ('.gz'));
        $opener = {
            decompress => qq"gunzip $args{input}",
            input => $new_input,
            cleanup => qq"rm $args{input}"
        };
    }
    else {
        my $new_input = $args{input};
        $opener = {
            input => $new_input,
        };
    }
    return($opener);
}

=head2 C<Get_FD>

  Returns an arbitrary command in a process substitution which may be used instead of a filename.
  Unfortunately, if the command requires  O_SEEK, then this will not work.

=cut

sub Get_FD {
    my ($class, %args);
    if ((scalar(@_) % 2) == 1) {
        ($class, %args) = @_;
    }
    else {
        %args = @_;
    }
    my $ignored = q{};
    $ignored = $args{ignored} if ($args{ignored});
    my $opener = $args{input};
    if ($args{input} =~ /\.xz$/ && $ignored !~ /xz/) {
        $opener = qq"<(xz -dc $args{input})";
    }
    elsif ($args{input} =~ /\.bz2$/ && $ignored !~ /bz2/) {
        $opener = qq"<(bunzip2 -c $args{input})";
    }
    elsif ($args{input} =~ /\.gz$/ && $ignored !~ /gz/) {
        $opener = qq"<(gunzip -c $args{input})";
    }
    return($opener);
}

sub Get_FH {
    my ($class, %args);
    if ((scalar(@_) % 2) == 1) {
        ($class, %args) = @_;
    }
    else {
        %args = @_;
    }
    my $opener = undef;
    if ($args{input} =~ /\.xz$/) {
        $opener = 'xz -dc';
    }
    elsif ($args{input} =~ /\.bz2$/) {
        $opener = 'bunzip2 -c';
    }
    elsif ($args{input} =~ /\.gz$/) {
        $opener = 'gunzip -c';
    }
    my $fh;
    my $cmd;
    if (defined($opener)) {
        $cmd = $opener;
        if ($args{suffix}) {
            $cmd .= qq" $args{input} $args{suffix} |";
        }
        else {
            $cmd .= qq" $args{input} |";
        }
    }
    else {
        if ($args{suffix}) {
            $cmd = qq"cat $args{input} $args{suffix} |";
        }
        else {
            $cmd = qq"<$args{input}";
        }
        $fh = FileHandle->new($cmd);
    }
    $fh = FileHandle->new($cmd);
    return($fh);
}

sub Get_Script_Filename {
    my ($class, %args) = @_;
    my $state_ref = $class->current_state;
    my %state = %{$state_ref};
    my $script_file = qq"$state{basedir}/$state{script_dir}/$state{jprefix}$state{jname}.sh";
    my $exists_test = 0;
    $exists_test = 1 if (-r $script_file);
    my $exists_count = 0;
    while ($exists_test) {
        $exists_count++;
        my $suffix_base = basename($state{basedir});
        my $suffix = qq"_${suffix_base}_${exists_count}_";
        $script_file = qq"$state{basedir}/$state{script_dir}/$state{jprefix}${suffix}$state{jname}.sh";
        unless (-r $script_file) {
            $exists_test = 0;
            print "Writing script to: ${script_file}\n";
        }
    }    ## Add a counter to the suffix until we have a free filename.
    return($script_file);
}

=head2 C<Get_Term>

Used to get a useful interactive terminal when more options are required.
Ideally, this should get us a terminal with tab completion.

=cut
sub Get_Term {
    my $term = undef;
    if (-t STDOUT) {
        $term = Term::ReadLine->new('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $OUT = $term->OUT || \*STDOUT;
        $Term::UI::VERBOSE = 0;
        $term->ornaments(0);
    }
    else {
        print "This process does not have a controlling terminal.\n";
    }
    return($term);
}

=head2 C<Get_Job_Name>

This attempts to make a reasonable job name given a filename or
cwd().  In practice it just strips off the suffixes from the input
file(s) and uses that.

=cut
sub Get_Job_Name {
    my ($class, $args) = @_;
    my $name = 'unknown';
    $name = $args->{input} if ($args->{input});
    if ($name =~ /[:;,\s+]/) {
        my @namelst = split(/[:;,\s+]/, $name);
        $name = $namelst[0];
    }
    $name = basename($name, split(/,/, $class->{suffixes}));
    $name = basename($name, split(/,/, $class->{suffixes}));
    $name =~ s/\-trimmed_*//g;
    $name =~ s/_*forward_*//g;
    $name =~ s/_*[R|r]1_*//g;
    return($name);
}

=head2 C<Get_Input>

    Get_Input() attempts to standardize the inputs passed to Bio::Adventure.
    It returns a stringified and standard representation of the likely
    input(s) to the script.

    There are a few problems with how I send input to these scripts:
    Sometimes I put in --hpgl hpgl0415 when I mean --hpgl HPGL0415,
    Sometimes I put in --input hpgl0415.fastq  when I mean hpgl0415.fastq.(gz|xz),
    Sometimes I put in -i hpgl0415.fastq.(gz|xz) when I mean hpgl0415.fastq,
    Sometimes I put in -i hpgl0415.fastq when I mean hpgl0415-trimmed.fastq(.gz|.xz)
    Sometimes I put in -i hpgl0415_forward.fastq:hpgl0415_reverse.fastq

    So, this function should make this unambiguous and consistent no matter what I type

=cut
sub Get_Input {
    my $class = shift;
    my $input = $class->{input};
    my $id = $class->{hpgl};
    my $actual = "";

    my @input_list = undef;
    if ($input =~ /:/) {
        @input_list = split(/:/, $input);
    }
    else {
        push(@input_list, $input);
    }

    foreach my $in (@input_list) {
        my $low = $in;
        $low =~ tr/[A-Z]/[a-z]/;
        unless ($in eq $low) {
            qx(qq"mv ${in}.fastq ${low}.fastq");
            if (-r "${in}_forward.fastq") {
                qx(qq"mv ${in}_forward.fastq ${low}_forward.fastq");
            }
            if (-r "${in}_reverse.fastq") {
                qx(qq"mv ${in}_reverse.fastq ${low}_reverse.fastq");
            }
        }

    }
    $input =~ tr/[A-Z]/[a-z]/;
    return($input);
}

=head2 C<Get_Vars>

  Handle the peculiar mix of instance options held in $class->{options},
  the set of arguments passed to %args, a list of required and potentially
  missing options, and whatever default values one wishes to set.

  my $options = $class->Get_Vars(args => \%args, required => ['input', 'genome'], bob => 'jane');

=cut

sub Juggle {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(args => \%args);
    ## I really ought to make a ForkManager type, that could be super-fun to play with
    my $runner;
    if ($options->{cluster} eq 'slurm') {
        $runner = Bio::Adventure::Slurm->new();
    }
    elsif ($options->{cluster} eq 'torque') {
        $runner = Bio::Adventure::Torque->new();
    }
    elsif ($options->{cluster} eq 'bash') {
        ## I should probably have something to handle gracefully bash jobs.
        $runner = Bio::Adventure::Bash->new();
    }
    elsif ($options->{cluster} eq '0') {
        ## On occasion I set cluster to 0 which is bash.
        $runner = Bio::Adventure::Bash->new();
    }
    else {
        carp("Could not find sbatch, qsub, nor bash.");
        print "Assuming this is running on a local shell.\n";
        $runner = Bio::Adventure::Bash->new();
    }
    my $result = $runner->Juggle($class, %args);
}

sub Load_Vars {
    my ($class, %args) = @_;
    my $num_changed = 0;
    for my $varname (keys %args) {
        if ($varname eq 'input') {
            local $Storable::Eval = 1;
            my $input_options = lock_retrieve($args{$varname});
            for my $opt (keys %{$input_options}) {
                $num_changed ++;
                $class->{$opt} = $input_options->{$opt};
            }
        }
        else {
            $num_changed++;
            $class->{$varname} = $args{$varname};
        }
    }
    return($num_changed);
}

=head2 C<Reset_Vars>

  Make sure that the environment is returned to a useful default state when perturbed.

=cut
sub Reset_Vars {
    my ($class, %args) = @_;
    my %original_values = (
        array_string => undef,
        jbasename => undef,
        jdepends => undef,
        jmem => undef,
        jname => undef,
        jname_prefix => undef,
        jprefix => undef,
        jstring => undef,
        jwalltime => undef,
        language => 'bash',
        prescript => q{},
        postscript => q{},
        subjob => undef,
        shell => '/usr/bin/env bash',);
    $class->{last_prefix} = $class->{jprefix};
    for my $k (keys %original_values) {
        $class->{$k} = $original_values{$k};
        $class->{variable_current_state}->{$k} = $original_values{$k};
    }
    return($class);
}

sub Last_Stat {
    my ($class, %args) = @_;
    if (!defined($args{input})) {
        print "Input for Last_Stat was undefined.\n";
        return(undef);
    }
    my $input = FileHandle->new("<$args{input}");
    if (!defined($input)) {
        print "The input file for Last_Stat appears to be empty.\n";
        return(undef);
    }
    my ($line, $last);
    while ($line = <$input>) {
        chomp $line;
        ## I am doing this due to terminal newlines on the file, there is probably a much smrtr way.
        $last = $line unless ($line =~ /^$/);
    }
    $input->close();
    return($last);
}

sub Process_Input_String {
    my ($class, %args) = @_;
    my $input_datum = $args{input};
    my $input_string = $input_datum->{input};
    my $interpreter = undef;
    if ($input_string =~ /^(\w+\s*:)?\s*(.+)$/) {
        $interpreter = $1;
        my $input_command = $2;
        $interpreter =~ s/:$//g;
        my $inputter = FileHandle->new("${input_command} |");
        my $in_string = q{};
      IN: while (my $line = <$inputter>) {
            chomp $line;
            $in_string .= qq"$line:";
        }
        $inputter->close();
        $input_string = $in_string;
    }
    return($input_string);
}

=head2 C<Passthrough_Args>

Add some flexibility in how arbitrary arguments get passed through to
various tools invoked by cyoa.  Using this function you may have any
number of extra args passed through to the final script by separating
them.  If it sees a single letter, it will prefix a single dash.  If
there are multiple letters it will double-dash it.  Conversely, if
the arbitrary string starts with a dash, it will leave it alone.

Thus:

--arbitrary ',--funkytown=bob;L 10,very-sensitive'

Will get the following arguments passed to your downstream tool:

'--funkytown=bob -L 10 --very-sensitive'.

 Note, since trimomatic uses colons, I excluded them from the default delimiter.

=cut
sub Passthrough_Args {
    my ($class, %args) = @_;
    my $argstring = $args{arbitrary};
    my $new_string = q{};
    return($new_string) unless (defined($argstring));
    my $splitter = qr/[:;,]/;
    for my $arg (split($splitter, $argstring)) {
        if ($arg =~ /^\w{1}$|^\w{1}\W+/) {
            ## print "Single letter arg passthrough\n";
            $arg = qq" -${arg} ";
        }
        elsif ($arg =~ /^\w{2}/) {
            ## print "Multi letter arg passthrough\n";
            $arg = qq" --${arg} ";
        }
        elsif ($arg =~ /^\-/) {
            $arg= qq" ${arg} ";
        }
        $new_string .= $arg;
    }
    return($new_string);
}

=head2 C<Read_Genome_Fasta>

Read a fasta file and return the chromosomes.

=cut
sub Read_Genome_Fasta {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        extra => 0);
    my $chromosomes = {};
    ## FIXME: This should be removed and I need to
    ## ensure that all calls use either fasta or genome.
    my $input_file = q{};
    if (defined($options->{genome})) {
        $input_file = $options->{genome};
    }
    else {
        if (defined($options->{input})) {
            $input_file = $options->{input};
        }
        elsif (defined($options->{fasta})) {
            $input_file = $options->{fasta};
        }
    }
    my $fasta_name = basename($input_file, ('.gz', '.xz', '.bz2'));
    $fasta_name = basename($fasta_name, ('.fasta', '.fa', '.ffn', '.fsa',));
    my $genome_file = qq"$options->{basedir}/$options->{script_dir}/pdata/${fasta_name}.pdata";
    ##if (-r $genome_file) {
    ##    $chromosomes = lock_retrieve($genome_file);
    ##} else {
    my $fh = Bio::Adventure::Get_FH(input => $input_file);
    my $input_genome = Bio::SeqIO->new(-fh => $fh, -format => 'Fasta');
    while (my $genome_seq = $input_genome->next_seq()) {
        next unless(defined($genome_seq->id));
        my $id = $genome_seq->id;
        my $sequence = $genome_seq->seq;
        my $length = $genome_seq->length;
        my $empty_forward = [];
        my $empty_reverse = [];
        if ($options->{extra}) {
            for my $c (0 .. $length) {
                $empty_forward->[$c] = 0;
                $empty_reverse->[$c] = 0;
            }
        }
        $chromosomes->{$id}->{forward} = $empty_forward;
        $chromosomes->{$id}->{sequence} = $sequence;
        $chromosomes->{$id}->{reverse} = $empty_reverse;
        $chromosomes->{$id}->{obj} = $genome_seq;
    }          ## End reading each chromosome
    ## my $stored = lock_store($chromosomes, $genome_file);
    ##} ## End checking for a .pdata file
    return($chromosomes);
}

=head2 C<Read_Genome_GFF>

Read a GFF file and extract the annotation information from it.
Now that I am more comfortable with the various gff parsers, this seems dumb.

=cut
sub Read_Genome_GFF {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(
        args => \%args,
        required => ['gff'],
        gff_type => 'gene',
        gff_tag => 'locus_tag');
    my $annotation_in = Bio::Tools::GFF->new(-file => "$options->{gff}", -gff_version => 3);
    my $gff_out = {
        stats => {
            chromosomes => [],
            lengths => [],
            feature_names => [],
            cds_starts => [],
            cds_ends => [],
            inter_starts => [],
            inter_ends => [],
        },
    };
    my $gff_name = basename($options->{gff}, ['.gff']);
    ##my $genome_file = qq"$options->{basedir}/${gff_name}.pdata";
    print "Reading $options->{gff}, seeking features tagged (last column ID) $options->{gff_tag},
 type (3rd column): $options->{gff_type}.\n";
    ##if (-r $genome_file && !$options->{debug}) {
    ##    $gff_out = lock_retrieve($genome_file);
    ##} else {
    print "Starting to read gff: $options->{gff}\n" if ($class->{debug});
    my $hits = 0;
    my @chromosome_list = @{$gff_out->{stats}->{chromosomes}};
    my @feature_names = @{$gff_out->{stats}->{feature_names}};
    my @cds_starts = @{$gff_out->{stats}->{cds_starts}};
    my @inter_starts = @{$gff_out->{stats}->{inter_starts}};
    my @cds_ends = @{$gff_out->{stats}->{cds_ends}};
    my @inter_ends = @{$gff_out->{stats}->{inter_ends}};
    my $start = 1;
    my $old_start = 1;
    my $end = 1;
    my $old_end = 1;
  LOOP: while(my $feature = $annotation_in->next_feature()) {
        next LOOP unless ($feature->primary_tag eq $options->{gff_type});
        $hits++;
        my $location = $feature->location;
        $old_start = $start;
        $old_end = $end;
        $start = $feature->start();
        $end = $feature->end();
        my $strand = $location->strand();
        my @ids = $feature->each_tag_value($options->{gff_tag});
        my $id = "";
        my $gff_chr = $feature->{_gsf_seq_id};
        ## my $gff_chr = $feature->display_name;
        my $gff_string = $annotation_in->gff_string($feature);
        foreach my $i (@ids) {
            $i =~ s/^cds_//g;
            $i =~ s/\-\d+$//g;
            $id .= "$i ";
        }
        $id =~ s/\s+$//g;
        my @gff_information = split(/\t+/, $gff_string);
        my $description_string = $gff_information[8];
        my $orf_chromosome = $gff_chr;
        ## Add the chromosome to the list of chromosomes if it is not there already.
        push(@chromosome_list, $orf_chromosome) if ($orf_chromosome !~~ @chromosome_list);
        push(@feature_names, $id);
        push(@cds_starts, $start);
        push(@inter_starts, $old_start + 1);
        push(@cds_ends, $end);
        ##push(@inter_ends, $old_start-1);
        push(@inter_ends, $start - 1);
        my $annot = {
            id => $id,
            start => $start,  ## Genomic coordinate of the start codon
            end => $end,      ## And stop codon
            strand => $strand,
            description_string => $description_string,
            chromosome => $gff_chr,
        };
        $gff_out->{$gff_chr}->{$id} = $annot;
    }                     ## End looking at every gene in the gff file
    $gff_out->{stats}->{chromosomes} = \@chromosome_list;
    $gff_out->{stats}->{feature_names} = \@feature_names;
    $gff_out->{stats}->{cds_starts} = \@cds_starts;
    $gff_out->{stats}->{cds_ends} = \@cds_ends;
    $gff_out->{stats}->{inter_starts} = \@inter_starts;
    $gff_out->{stats}->{inter_ends} = \@inter_ends;
    print STDERR "Not many hits were observed, do you have the right feature type?  It is: $options->{gff_type}\n" if ($hits < 1000);
    ## if (-f $genome_file) {
    ##         unlink($genome_file);
    ##     }
    ##     my $stored = lock_store($gff_out, $genome_file);
    ## } ## End looking for the gff data file
    return($gff_out);
}

=head2 C<Record_Options>

 After submitting a job to a local shell/cluster/etc, write down some of the options used.

 My thinking is to be able to, in future invocations, skip passing some options like input.

=cut
sub Record_Options {
    my ($class, %args) = @_;
    my $runner = $args{runner};
    my $result = $args{result};
    my $arg_fh = FileHandle->new(qq">>$class->{output_basedir}/args.tsv");
    my $method = Bio::Adventure::Config::Get_Caller();
    my @written = ();
    if ($result->{stdout}) {
        print $arg_fh "${method}\tstdout\t$result->{stdout}\n";
        push(@written, $result->{stdout});
    }
    for my $k (keys %{$result}) {
        if ($k =~ /main|output/) {
            if (!defined($result->{$k})) {
                warn("When writing options, ${k} has an undefined result.");
                print $arg_fh "${method}\t${k}\tundefined_result\n";
            }
            else {
                print $arg_fh "${method}\t${k}\t$result->{$k}\n";
            }
            push(@written, $k);
        }
    }
    $arg_fh->close();
    return(@written);
}

=head2 C<Reorder_Fasta>

 Some tools which pass back a set of chromosomes do not always return
 them in the most sensible order.  This attempts to fix that, and make
 them prettier.

=cut
sub Reorder_Fasta {
    my ($class, %args) = @_;
    my ($options, $paths, $inputs, $outputs) = $class->Get_Vars(args => \%args);
    my $genome = $class->Read_Genome_Fasta(%args);

    my @chromosomes = sort { $a cmp $b } keys %$genome;
    my $out = FileHandle->new(">$options->{output}");
    foreach my $c (@chromosomes) {
        my $start = $genome->{c};
        my $end = join("\n", ($start =~ m/.{1,80}/g));
        print $out ">$c
$end
";
    }
    $out->close();
}

sub scalar_which {
    my $exe = $_[0];
    my $path = which($exe);
    return($path);
}

sub self_interp {
    my ($class, %args) = @_;
    my %hash = %{$args{hashref}};
    for my $k (keys %hash) {
        my $value = $hash{$k};
        next unless (defined($value));
        if ($value =~ /\$/) {
            my $matched = $value =~ m/^(.*)?\$\{*(?<varname>[a-zA-Z_]+)\}*(.*)$/;
            my $varname = $+{varname};
            next unless (defined($varname));
            if ($hash{$varname}) {
                my $new = interpolate($hash{$k}, \%hash);
                $hash{$k} = $new;
            }
        }
    }
    return(\%hash);
}

=head2 C<Submit>

 Pass a job to slurm or torque after doing a little sanity checking.  I
 am increasinbly certain that keeping a copy of the options in a .pdata
 file is not needed.

=cut
sub Submit {
    my ($class, %args) = @_;
    my $state = $class->current_state;
    my %state_deref = %{$state};
    my $all_options = $class->options;
    my $last_optref = $all_options->[$#$all_options];
    my %last_options = %{$last_optref};
    ## Bring together the current options and state.
    ##OPTS: for my $k (keys %last_options) {
    ##}

    ## my $modules = $state->{modules};
    my $modules = $state_deref{modules};
    ## my $avail_clusters = $state->{available_clusters};
    my $avail_clusters = $state_deref{available_clusters};
    ## Before creating the runner, add the required modules.
    my $orig_env = $class->Module_Loader(action => 'load',
                                         modules => $modules);
    my $runner;
    ## if ($state->{cluster} eq 'slurm') {
    if ($state_deref{cluster} eq 'slurm') {
        if (any {$_ eq 'slurm'} @{$avail_clusters}) {
            $runner = Bio::Adventure::Slurm->new();
        } else {
            $runner = Bio::Adventure::Bash->new();
        }
        ## } elsif ($state->{cluster} eq 'torque') {
    } elsif ($state_deref{cluster} eq 'torque') {
        if (any {$_ eq 'torque'} @${avail_clusters}) {
            $runner = Bio::Adventure::Torque->new();
        } else {
            $runner = Bio::Adventure::Bash->new();
        }
        ## } elsif ($state->{cluster} eq 'bash') {
    } elsif ($state_deref{cluster} eq 'bash') {
        ## I should probably have something to handle gracefully bash jobs.
        $runner = Bio::Adventure::Bash->new();
        ## } elsif ($state->{cluster} eq '0') {
    } elsif ($state_deref{cluster} eq '0') {
        ## On occasion I set cluster to 0 which is bash.
        $runner = Bio::Adventure::Bash->new();
    } else {
        carp("Could not find sbatch, qsub, nor bash.");
        print "Assuming this is running on a local shell.\n";
        $runner = Bio::Adventure::Bash->new();
    }
    ## Add the current options to the runner:
    my $result = $runner->Submit($class, %args);

    ## After the runner is finished, unload the modules.
    my $reset = $class->Module_Reset(env => $orig_env);
    ## We should also unset a few pieces of the state
    ## so that they may be in future jobs:
    ## In order to do this with maximum cleanliness, I should record
    ## the set of things set in get_gfv and clean them out here.

    ## When keys from modified_state evaluate to a defined value, then
    ## They came from something before the current function call
    ## and should be returned to that state.  If they are not defined,
    ## just delete them.
    print "TESTME: About to check modified_state: $state_deref{main}\n";
    my %modified_state = %{$state_deref{modified_state}};
  RESETTERS: for my $reset_key (keys %modified_state) {
        my $old_value = $modified_state{$reset_key};
        ## The config file uses the explicit string 'undef'
        ## for placeholders, so skip them.
        if (defined($old_value)) {
            ## next RESETTERS if ($old_value eq 'undef');
            $state_deref{$reset_key} = $old_value;
        }
        else {
            my $deleted = delete($state_deref{$reset_key});
        }
    }
    %modified_state = ();
    $state_deref{modified_state} = \%modified_state;
    my $reset_state = $class->current_state(\%state_deref);
    ## Now the class' state should be similar/identical
    ## to what it was before we messed with it.

    my @jobs = @{$class->jobs};
    push(@jobs, $result);
    $class->jobs(\@jobs);
    my @job_ids = @{$class->job_ids};
    push(@job_ids, $result->{job_id});
    $class->job_ids(\@job_ids);
    return($result);
}

sub Test_Job {
    my ($class, %args) = @_;
    my $test = $class->Bio::Adventure::SNP::Test_Worker(%args);
}

=head2 C<Wait>

 Wait on a job, when possible collect information about it.

=cut
sub Wait {
    my ($class, %args) = @_;
    my %options = %{$class->current_state};
    my $all_inputs = $class->inputs;
    my $inputs = $all_inputs->[$#$all_inputs];
    my $status = undef;
    if ($options{cluster} eq 'slurm') {
        $status = $class->Bio::Adventure::Slurm::Wait(%args);
    }
    elsif ($options{cluster} eq 'torque') {
        $status = $class->Bio::Adventure::Torque::Wait(%args);
    }
    elsif ($options{cluster} eq 'bash') {
        ## I should probably have something to handle gracefully bash jobs.
        $status = $class->Bio::Adventure::Bash::Wait(%args);
    }
    elsif ($options{cluster} eq '0') {
        ## On occasion I set cluster to 0 which is bash.
        $status = $class->Bio::Adventure::Bash::Wait(%args);
    }
    else {
        carp("Could not find sbatch, qsub, nor bash.");
        print "Assuming this is running on a local shell.\n";
        return(undef);
    }
    return($status);
}

sub Adventure_Help {
    my ($class, %args) = @_;
    my %methods = %{$class->{methods}};
    print qq"The command line program 'cyoa' has a series of shortcuts intended to make it easy to use and flexible.
The following comprises the set of strings you may feed it as 'methods':\n";
    my $c = 0;
    for my $k (sort keys %methods) {
        my $sep = '\t';
        if (($c % 3) == 0) {
            $sep = '\n';
        }
        else {
            $sep = '\t';
        }
        $c++;
        print "${k}${sep}";
    }
    print "\n";
    print qq"cyoa uses GetOptions, so you can shortcut all the 'methods', so:
'cyoa --task ri --method cut --input test.fastq' calls the cutadapt with
  options suitable for ribosome profiling data.  Conversely:
'cyoa --task rna --method top --input test.fastq' calls tophat assuming
  rnaseq data.
'cyoa --method blastsplit --query test.fasta --library nr --blast_tool blastp'
  Splits the test.fasta into a bunch of pieces (settable with --number), calls
  blastp on them, merges the outputs, and parses the result into a table of hits.
";
    $class->Help();
    return(0);
}

=head1 AUTHOR - atb

Email abelew@gmail.com

=head1 SEE ALSO

    L<Bio::Seq> L<autodie> L<local::lib> L<Bio::Adventure::Assembly> L<Bio::Adventure::Count>
    L<Bio::Adventure::QA> L<Bio::Adventure::Trim> L<Bio::Adventure::Align_Blast>
    L<Bio::Adventure::Align_Fasta> L<Bio::Adventure::Align> L<Bio::Adventure::Compress>
    L<Bio::Adventure::Convert> L<Bio::Adventure::Prepare> L<Bio::Adventure::Riboseq>
    L<Bio::Adventure::SeqMisc> L<Bio::Adventure::SNP> L<Bio::Adventure::TNSeq>

=cut

1;
