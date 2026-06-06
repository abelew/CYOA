package Bio::Adventure::Config;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use v5.14;
use Cwd qw"abs_path fastcwd";
use Data::Printer;
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

our @EXPORT = qw"Get_Menus Get_TODOs";
our @ISA = qw"Exporter";

## 202605 survey: Invoked when reading an ini file
sub Add_Config_File {
    my %args = @_;
    my $input = $args{input};
    my %current = %{$args{values}};
    my $config = AppConfig->new({ CREATE => '.*' , GLOBAL => {
        ARGCOUNT => AppConfig::ARGCOUNT_ONE } });
    my $info = $config->file($input);
    my $variable = q{};
    my $steps = {};
    my $defaults = {};
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
        } else {
            $defaults->{$variable} = $config_value;
        }
    }
    my $new = Merge_Current_Ini(current => %current, steps => $steps, defaults => $defaults);
    return($new);
}

## 202605 Survey: Invoked at the beginning of BUILD to get command line arguments.
sub Get_ARGV {
    my $cmdline_args = {};
    if ($#ARGV <= 0) {
        return $cmdline_args;
    }
    my $arg_string = q{};
    my @spec_strings = ();
    foreach my $a (@ARGV) {
        $arg_string .= "$a ";
        if ($a =~ /^\-/) {
            my $spec = $a;
            $spec =~ s/^\-*//g;
            push(@spec_strings, $spec);
        }
    }

    my %argv_information = ();
    my %conf_specification = ();
    for my $spec (@spec_strings) {
        $argv_information{$spec} = undef;
        my $def_string = qq"${spec}:s";
        $conf_specification{$def_string} = \$argv_information{$spec};
    }
    my $argv_result = GetOptions(%conf_specification);
    $cmdline_args = \%argv_information;
    $cmdline_args->{argv_string} = $arg_string;
    return($cmdline_args);
}

=head2 C<Get_Caller>

Figure out the calling function so that I may fill in some configuration information.

=cut
## 202605 Survey: called all over the place.
sub Get_Caller {
    my ($class, %args) = @_;
    my $caller_number = 2;
    $caller_number = $args{caller} if ($args{caller});
    my $attempts = 4;
    my $found_sub = undef;
    my $function;
    my $last_attempt = $caller_number + $attempts;
  COUNTER: for my $attempt ($caller_number .. $last_attempt) {
        my ($package, $filename, $line, $sub, $hasargs,
            $wantarray, $evaltext, $is_require, $hints, $bitmask, $hinthash) = caller($attempt);
        last COUNTER unless (defined($sub));
        if ($sub =~ /^Bio::.*Submit$/) {
            next COUNTER;
        } elsif ($sub =~ /^Bio::.*Get_Vars$/) {
            next COUNTER;
        } elsif ($sub =~ /^Bio::.*Get_Caller$/) {
            next COUNTER;
        } elsif ($sub =~ /^Bio::.*Update_State$/) {
            next COUNTER;
        }
        ## All callers which are needing modules have the form Bio::Adventure::Category::Function
        $function = $sub;
        $function =~ s/^(.*::)(\w+)$/$2/g;
        last COUNTER if (defined($function));
    }
    if (!defined($function)) {
        ## warn "The function is undefined in Config::Get_Caller().";
        $function = 'Default';
    }

    ## Now that we have the current caller, update a few pieces of the parent class.
    ## I made a decision to simultaneously store the set of all function defaults as well
    ## as the function defaults for the current caller.  I did this for two reasons:
    ## 1.  It should make it simpler in the various child functions to pull useful stuff.
    ## 2.  It should make it possible to check one function from another, I am thinking
    ##     this will be useful in the case of multi-function paths
    ## Thus, I will update the function_args slot of the class with the portion of
    ## function_defaults{} for this caller and also update spec with the portion of that
    ## which tells the default walltime/mem/etc.  Finally, I will set caller with whatever
    ## this function happens to be.
    if (!defined($function)) {
        return(undef);
    }
    ## Update the function args with this caller.
    my %all_args = %{$class->function_defaults};
    if (!defined($all_args{$function})) {
        return(undef);
    }
    my $function_defs = $all_args{$function};
    my %fn_args = %{$all_args{$function}};
    $class->function_args(\%fn_args);
    $class->caller($function);
    return($function);
}

## My thought here: Make it possible to have a set of 'global' defaults
## which may be overridden at a few levels, branching down to the indivdual sample.
## I organize my data as /type-of-project/individual-project/preprocessing/sample/
## so there are likely 4 levels at which one may most likely want to make changes.

## 202605 Survey: Invoked by BUILD to read ini-style config directives.
sub Get_Ini {
    my %args = @_;
    my $filename = basename($args{input_file});
    my $ini_values = {};
    ## At the super-project level, e.g. all the RNASeq analyses may share defaults
    ## which are different from the TNSeq etc.
    my $super_super_super_file = qq"../../../${filename}";
    if (-f $super_super_super_file) {
        $ini_values = Add_Config_File(
            input => $super_super_super_file,
            values => $ini_values,
        );
    }
    ## At the project directory, which is likely redundant to the preprocessing directory
    ## but I can imagine having a separate directory for other types of analysis
    my $super_super_file = qq"../../${filename}";
    if (-f $super_super_file) {
        $ini_values = Add_Config_File(
            input => $super_super_file,
            values => $ini_values,
        );
    }
    ## At the preprocessing directory
    my $super_file = qq"../${filename}";
    if (-f $super_file) {
        $ini_values = Add_Config_File(
            input => $super_file,
            values => $ini_values,
        );
    }
    ## Change defaults at the sample level
    my $file = qq"${filename}";
    if (-f $file) {
        $ini_values = Add_Config_File(
            input => $file,
            values => $ini_values,
        );
    }
    return($ini_values);
}

## 202605 survey: Invoked in BUILD to collect default values.
sub Get_JSON {
    my %args = @_;
    if (!defined($args{input_file})) {
        return({});
    }
    my $default_base = basename($args{input_file});
    if (!defined($default_base)) {
        return({});
    }
    my $input_file;
    my $dist_dir = dist_dir('Bio-Adventure');
    my $dist_file = qq"${dist_dir}/${default_base}";
    if (-r $args{input_file}) {
        $input_file = $args{input_file};
    } elsif (-f $dist_file) {
        $input_file = $dist_file;
    } else {
        say "Unable to find global default file, returning an empty hash ref.";
        return({});
    }
    my $json_string = q{};
    my $in = FileHandle->new(qq"<${input_file}");
    while (my $line = <$in>) {
        $json_string .= $line;
    }
    $in->close();
    my $json = JSON->new->allow_nonref;
    $json = $json->relaxed([1]);
    my $information = $json->decode($json_string);
    return($information);
}

=head2 C<Get_Menus>

This is responsible for printing out the cyoa menus and connecting the
various options to the functions in the package.

The top-level of the hash is the set of tasks, so stuff like sequence alignment,
mapping, indexing, conversion, counting, etc.

There are a few keys inside each of these; the names are the option which may be used
with --task to differentiate between contexts (e.g. if one wishes to use mapping
options specific for tnseq vs. snps, for example).  The message is just silly.
Finally, the choices are a set of strings which point to a function.
Choosing the number associated with the string will therefore lead to the
invocation of that function.

=cut
## 202605 Survey: Called in build to set up the menu, I want to move this to a config file.
sub Get_Menus {
    my $start = Bio::Adventure::Config::Get_JSON(
        input_file => 'menus.json',);
    my $menus = {};
    for my $k (keys %{$start}) {
        my %menu_elements = ();
        $menu_elements{name} = $start->{$k}{name};
        $menu_elements{message} = $start->{$k}{message};
        my %choice_elements = ();
        my @inner_json = @{$start->{$k}{choices}};
        for my $inner_element (@inner_json) {
            my $inner_key = qq"($inner_element->{name}) $inner_element->{description}";
            my $fn_name = $inner_element->{function};
            my $fn_deref = \&${fn_name};
            $choice_elements{$inner_key} = $fn_deref;
        }
        $menu_elements{choices} = \%choice_elements;
        $menus->{$k} = \%menu_elements;
    }
    return($menus);
}

## 202605 Survey: Called in BUILD to pull (currently only module) shell functions;
## The intent is to add conda/spack and whatever other future tools crop up.
## With that in mind, it kind of only makes sense in the context of module because
## (with the exception of the various condas) I expect no other tools to have multiple
## implementations -- but the existence of miniconda/anaconda/etc suggest it is worth while.
sub Get_Shell_Functions {
    my (%args) = @_;
    ## Check that the module command is available as a bash function.
    ## I was initially going to do this in BUILD(), but I think that might mess up jobs
    ## which set the jprefix parameter.
    ## This command, if run via bash, requires a controlling pty in order to find the
    ## bash function 'module'.  slurm and the various managers will not provide that
    ## so, we can use typescript to allocate one post-facto.
    my $shell_text = q{};
    my $module_string = Bio::Adventure::Software::Get_Modulecmd(@_);
    $shell_text .= $module_string;
    my $spack_string = Bio::Adventure::Software::Get_Spackcmd(@_);
    $shell_text .= $spack_string;
    my $conda_string = Bio::Adventure::Software::Get_Condacmd(@_);
    $shell_text .= $conda_string;
    return($shell_text);
}

=head2 C<Get_Input_Info>

Given a file/directory name, provide a set of paths which will
hopefully prove useful for the various functions in CYOA.  I think
this function is a good candidate for replacing Check_Input() below.

=cut
## 202605 survey: Called by BUILD to create a single data structure
## describing the inputs to cyoa.
sub Get_Input_Info {
    my ($class, @inputs) = @_;
    my %ret = ();
    my $input_string = q{};
    return(undef) unless($inputs[0]);
    if ($inputs[0] =~ /[:;,\s]+/) {
        @inputs = split(/[:;,\s]+/, $inputs[0]);
    }
    my $num_inputs = scalar(@inputs);
    if ($num_inputs == 0) {
        die("This requires an input filename.");
    }
    my @outputs = ();
    for my $in (@inputs) {
        $input_string .= qq"${in}:";
        my $filename = basename($in);
        my $filebase_compress = basename($in, ('.gz', '.xz', '.bz2'));
        ## FIXME: Replace this with options->{suffix_array}
        my @exts = ('.bam', '.fastq', '.fasta', '.fsa', '.faa', 'fna', '.fa', '.ffn',
                    'gff', '.gff3', '.gbk', '.gbf', '.sam', '.sqn', 'tbl', '.tsv', '.txt');
        my $filebase_extension = basename($filebase_compress, @exts);
        my $directory = dirname($in);
        my $dirname = basename($directory);
        ## Check that doing basename(dirname(input)) leaves us with something interesting
        if ($dirname eq '.') {
            $dirname = q{};
        }
        ## This test/path build is because abs_path only works on stuff which exists.
        if (! -e $directory) {
            make_path($directory);
        }
        my $full_path = abs_path($directory);
        $full_path .= "/${filename}";
        my $local_path = qq"${directory}/${filename}";

        ## Make a new job_basename which doesn't have some cruft
        my $jbasename = $filebase_extension;
        $jbasename =~ s/_forward.*//g;
        $jbasename =~ s/_R1.*//g;
        $jbasename =~ s/-trimmed.*//g;
        $jbasename =~ s/_trimmed.*//g;

        my $suffix = $filename;
        $suffix =~ s/\.[xgb]z$//g;
        if ($suffix =~ /\..*$/) {
            $suffix =~ s/^.*\.(.*)$/$1/g;
        }

        my %ret = (
            filename => $filename,
            filebase_compress => $filebase_compress,
            filebase_extension => $filebase_extension,
            directory => $directory,
            dirname => $dirname,
            localpath => $local_path,
            fullpath => $full_path,
            jbasename => $jbasename,
            basename => $jbasename,
            input => $local_path,
            suffix => $suffix,);
        push(@outputs, \%ret);
    }
    ## Ensure the current state has the input string.
    $input_string =~ s/:$//g;
    $class->{current_state}->{input} = $input_string;
    return(\@outputs);
}

=head2 C<Get_Class_Paths>

This should be responsible for acquiring paths used by anything in
this instance of the Adventure class.

TODO: Once this is settled, take out the function-specific material
and move it to Get_Function_Paths() which will contain stuff
like output_dir.

=cut
## 202605 Survey: Called by BUILD to set up and/or find all the various database directories
sub Get_Global_Paths {
    my ($class, %args) = @_;
    my $state_ref = $class->current_state;
    my %current_state = %{$state_ref};
    my $libtype = 'genome';
    $libtype = $current_state{libtype} if (defined($current_state{libtype}));
    my $libpath_prefix = qq"$current_state{libpath}/${libtype}";
    my $libdir_prefix = qq"$current_state{libdir}/${libtype}";
    ##my $cwd_name = basename(cwd());
    my %class_paths = ();
    ## Ok, something in cwd/basename is causing a problem.
    $class_paths{cwd_name} = basename(fastcwd());
    $class_paths{blast_dir} = qq"${libpath_prefix}/blast";
    $class_paths{codon_dir} = qq"${libpath_prefix}/codon_tables";
    $class_paths{fasta_dir} = qq"${libpath_prefix}/fasta";
    $class_paths{gbk_dir} = qq"${libpath_prefix}/genbank";
    $class_paths{gff_dir} = qq"${libpath_prefix}/gff";
    $class_paths{gtf_dir} = qq"${libpath_prefix}/gtf";
    $class_paths{hmm_dir} = qq"${libpath_prefix}/hmm";
    $class_paths{index_prefix} = qq"${libpath_prefix}/indexes";
    $class_paths{index_prefix_shell} = qq"${libdir_prefix}/indexes";
    make_path($class_paths{blast_dir}) unless (-d $class_paths{blast_dir});
    make_path($class_paths{codon_dir}) unless (-d $class_paths{codon_dir});
    make_path($class_paths{fasta_dir}) unless (-d $class_paths{fasta_dir});
    make_path($class_paths{gff_dir}) unless (-d $class_paths{gff_dir});
    make_path($class_paths{gtf_dir}) unless (-d $class_paths{gtf_dir});
    make_path($class_paths{gbk_dir}) unless (-d $class_paths{gbk_dir});
    make_path($class_paths{hmm_dir}) unless (-d $class_paths{hmm_dir});
    make_path($class_paths{index_prefix}) unless (-d $class_paths{index_prefix});
    if (defined($ENV{BLAST_DB})) {
        if ($ENV{BLAST_DB} =~ /:/) {
            my @blast = split(/:/, $ENV{BLAST_DB});
            unless (grep(/$class_paths{blast_dir}/, @blast)) {
                $ENV{BLAST_DB} = qq"$ENV{BLAST_DB}:$class_paths{blast_dir}";
            }
        } else {
            unless ($class_paths{blast_dir} eq $ENV{BLAST_DB}) {
                $ENV{BLAST_DB} = qq"$ENV{BLAST_DB}:$class_paths{blast_dir}";
            }
        }
    } else {
        $ENV{BLAST_DB} = $class_paths{blast_dir};
    }
    if (defined($current_state{species})) {
        $class_paths{fasta} = qq"$class_paths{fasta_dir}/$current_state{species}.fasta";
        $class_paths{fasta_shell} = qq"${libdir_prefix}/fasta/$current_state{species}.fasta";
        $class_paths{gff} = qq"$class_paths{gff_dir}/$current_state{species}.gff";
        $class_paths{gff_shell} = qq"${libdir_prefix}/gff/$current_state{species}.gff";
        $class_paths{gtf} = qq"$class_paths{gtf_dir}/$current_state{species}.gtf";
        $class_paths{gtf_shell} = qq"${libdir_prefix}/gtf/$current_state{species}.gtf";
        $class_paths{gbk} = qq"$class_paths{gbk_dir}/$current_state{species}.gbk";
        $class_paths{gbk_shell} = qq"${libdir_prefix}/genbank/$current_state{species}.gbk";
    }
    my $set = $class->_set_global_paths(\%class_paths);
    return(\%class_paths);
}

=head2 C<Get_TODOs>

The keys of possible_todos provide the GetOpt::Long options required
to figure out if the user is requesting the various functions.  When
using GetOpt::Long, if an options has a '+' suffix, its value is
incremented by one when it is set without a following argument.  Thus
doing --abricate will set the
todo_list{todo}{Bio::Adventure::Resistance::Abricate} value from 0 to 1.
The cyoa script will then iterate over the todo hash and look for
anything with a positive value and run its associated function.

=cut
## 202605 Survey: Called in BUILD
sub Get_TODOs {
    my %args = @_;
    my $todo_list = ();
    my $methods = Bio::Adventure::Config::Get_JSON(
        input_file => 'methods.json',);
    my %possible_todos = ();
    for my $k (keys %{$methods}) {
        my $new_key = qq"${k}+";
        my $value = $methods->{$k};
        $possible_todos{$new_key} = \$todo_list->{todo}{$value};
    }
    if (defined($args{method})) {
        my $m = '--' . $args{method};
        my @method = ($m,);
        my $array_result = GetOptionsFromArray(\@method, %possible_todos);
    }

    my @methods_to_run = ();
    my $set_of_todos = $todo_list->{todo};
    foreach my $job (keys %{$set_of_todos}) {
        my $thing = $set_of_todos->{$job};
        if ($thing) {
            if (!defined($args{task})) {
                print "Going on an adventure to ${job}.\n";
            } else {
                print "Going on an adventure to ${job} in the $args{task} context.\n";
            }
            my $final = \&${job};
            push(@methods_to_run, $final);
        }
    }
    return(\@methods_to_run);
}

## 202605 Survey: Responsible for providing all the local variables to arbitrary functions.
sub Get_Function_Vars {
    my ($class, $funcall_args) = @_;
    my $caller = $class->Bio::Adventure::Config::Get_Caller();
    my $state_ref = $class->current_state;

    ## Store the changed state elements here.
    ## They should get reset at the end of Submit()
    ## We will leave global defaults and instance args alone
    my %modified_state = ();
    ## To an annoying degree this is effectively redoing the
    ## BUILD function from Adventure.pm
    ## I should think about how I can clean this up so that
    ## I am not doing and undoing these declarations.
    ## The caveat though is that I want to have a clean slate for
    ## every downstream function call invoked by Adventure.pm
    ## Until I think that through, keep in mind that my
    ## order of precedence is:
    ## global_defaults < instance_args < function_defaults < function_args
    my %state = %{$state_ref};
    $state{caller} = lc($caller);
    my $global_defaults = $class->global_defaults;
    for my $k (keys %{$global_defaults}) {
        $state{$k} = $global_defaults->{$k};
    }
    my %instance_args = %{$class->instance_args};
    for my $j (keys %instance_args) {
        $state{$j} = $instance_args{$j};
    }
    my $function_state = $class->function_defaults;
    my %gl_paths = %{$class->global_paths};
    my $last_input;
    if ($funcall_args->{input}) {
        $class->{current_state}->{input} = $funcall_args->{input};
        $last_input = $class->Bio::Adventure::Config::Get_Input_Info(
            $funcall_args->{input});
        my $number_inputs = $class->Append_Inputs($last_input);
    } else {
        $last_input = $class->Get_Last_Input();
        $class->{current_state}->{input} = $last_input;
    }
    $class->current_inputs($last_input);
    ## Ensure that current_state->input matches this.

    my %all_function_args = %{$class->function_defaults};
    ## Note, I think one must make a copy of function_args here so that future invocations
    ## will not already have the variables interpolated.
    my %function_args = %{$all_function_args{$caller}};
    for my $k (keys %function_args) {
        if (defined($state{$k})) {
            $modified_state{$k} = $state{$k};
        } else {
            $modified_state{$k} = undef;
        }
        $state{$k} = $function_args{$k};
    }

    MODIFIED: for my $j (keys %{$funcall_args}) {
          if (defined($state{$j})) {
              ## I set aside the string 'undef' for
              ## key/value pairs which should be reserved
              ## but are left initially undefined, so
              ## have a specific check for them here.
              next MODIFIED if ($state{$j} eq 'undef');
                  $modified_state{$j} = $state{$j};
          } else {
              $modified_state{$j} = undef;
          }
          $state{$j} = $funcall_args->{$j};
      }

    ## Add the global paths
    for my $i (keys %gl_paths) {
        $state{$i} = $gl_paths{$i};
    }

    ## Define the set of accessible interpolation variables here.
    my @info_array = @{$last_input};
    my %first_file_info = %{$info_array[0]};
    my $first_input_filename = $first_file_info{filename};
    my @input_array = ();
    my @output_array = ();
    my @input_base_array = ();
    my @output_base_array = ();
    for my $input (@info_array) {
        push(@input_array, $first_file_info{filename});
        push(@input_base_array, $first_file_info{basename});
    }

    $state{output_prefix} = qq"$state{output_basedir}/$state{jprefix}";
    ## The suffix of the input file(s) excluding any compression information.
    $state{input_suffix} = $first_file_info{suffix};
    ## The directory name of the input file(s)
    $state{input_dirname} = $first_file_info{dirname};
    ## The basename of the input file(s) excluding erroneous stuff like compress/suffix/r1/etc
    $state{output_base} = $first_file_info{basename};
    ## Array of input filenames
    $state{input_array} = \@input_array;
    ## Array of the input basenames.
    $state{input_base_array} = \@input_base_array;
    ## The first input file excluding compression/suffixes
    $state{nosuffix} = qq"$first_file_info{dirname}/$first_file_info{basename}";
    ## The full output directory path, this may need some flexibility and/or to be interpolated.
    ## E.g. I want to potentially have the flexibility to add parameters/genome etc to this.
    $state{output_dir} = qq"$state{output_prefix}_$state{jname}";

    ## I want to run most/all jobs through time, so record that filename now.
    $state{time_file} = qq"$state{output_dir}/$state{jname}.time";

    ## The primary inputs/outputs are defined,
    ## now add anything from the functions.json for this function

    ## Add the input information to %available
    ## Note that the input information is a list of hash references
    ## and we are going to usually/always use the first input for stuff like basename/filename/etc.
    for my $k (keys %first_file_info) {
        $state{$k} = $first_file_info{$k};
    }
    ## output_dir is provided by fn_paths
    ## Add any arguments added via a function call.
    for my $k (keys %{$funcall_args}) {
        $state{$k} = $funcall_args->{$k};
    }
    for my $e (keys %ENV) {
        if ($e eq 'HOME' || $e eq 'PATH') {
            $state{$e} = $ENV{$e};
        }
    }
    ## Add environment variables to available

    my $first_file_base = $last_input->[0]->{filebase_extension};
    ## This depends on fn_paths
    $state{r1_base} = qq"$state{output_dir}/${first_file_base}";
    if ($last_input->[1]) {
        my $second_file_base = $last_input->[1]->{filebase_extension};
        $state{r2_base} = qq"$state{output_dir}/${second_file_base}";
    }

    ## Some specifications may provide 2 outputs but only use 1
    ## E.g. fastp can have a r1+r2 output or just r1
    ## We can check for this case by checking r2_base.
    while (my ($key) = each(%state)) {
        if ($key =~ /_r1$/) {
            if ($state{r1_base}) {
                ## print "keeping $current_state{r1_base}\n";
            } else {
                delete $state{$key};
            }
        } elsif ($key =~ /_r2$/) {
            if ($state{r2_base}) {
                ## print "Keeping $current_state{r2_base}\n";
            } else {
                ## print "DELETING $key\n";
                delete $state{$key};
            }
        }
    }

    ## Now we can create the outputs interpolation variable
    my $outputs = q{};
    STATE: for my $k (sort keys %state) {
          next STATE if ($k =~ /_comment$/);
          if ($k =~ /^output/) {
              $outputs .= qq"$state{$k}:";
          }
      }

    my $sw_added = 0;
    if ($state{modules}) {
        $sw_added++;
        $state{module_string} = qq"module add $state{modules}";
    }
    if ($state{conda}) {
        $sw_added++;
        $state{conda_string} = qq"conda load $state{conda}";
    }
    if ($state{spack}) {
        $sw_added++;
        $state{spack_string} = qq"spack load $state{spack}";
    }
    if ($sw_added == 0) {
        warn("This expects at least one software package to be added to the environment.");
    }

    $outputs =~ s/:$//g;
    $state{outputs} = $outputs;
    $function_args{inputs} = $last_input;
    my $interp_state_ref = $class->self_interp(hashref => \%state);
    $interp_state_ref = $class->self_interp(hashref => $interp_state_ref);
    ## Add this set of options to the option list in the cyoa object.
    ## Thus, we can trivially pull all of those options in future jobs.
    my @option_array = @{$class->options};
    push(@option_array, $interp_state_ref);
    $interp_state_ref->{modified_state} = \%modified_state;
    $class->options(\@option_array);
    $class->current_state($interp_state_ref);
    return($interp_state_ref);
}

## For now, just overwrite keys from global->function->ini->json->cmdline
## in order to create the hash of the current state.
## Once I start getting into the dirty implementation details
## I think I will want to use those portions of the function hash to make this
## into something which sets the parameters of each step of a process.
## But while I am just experimenting with how I want to work with this, a simpler
## approach makes sense.

## 202605 survey: Called in BUILD to bring together ini/json/etc
sub Merge_Args {
    my %args = @_;
    ## This is the reverse order of precedence when overwriting parameters
    ## e.g. globals are overwritten by the ini by json by cmdline by instance
    my $global = $args{global};
    my $ini = $args{ini};
    my $json = $args{json};
    my $cmdline = $args{cmdline};
    my $instance = $args{instance};
    my %current;
    for my $k (keys %{$global}) {
        $current{$k} = $global->{$k};
    }
    for my $k (keys %{$ini}) {
        $current{$k} = $ini->{$k};
    }
    for my $k (keys %{$json}) {
        $current{$k} = $json->{$k};
    }
    for my $k (keys %{$cmdline}) {
        $current{$k} = $cmdline->{$k};
    }
    for my $k (keys %{$instance}) {
        $current{$k} = $instance->{$k};
    }
    return(\%current);
}

sub Merge_Current_Ini {
    my %args = @_;
    my %current = $args{current};
    my $steps = $args{steps};
    my $defaults = $args{defaults};
    for my $def_key (keys %{$defaults}) {
        my $current_key = $def_key;
        while (defined($current{$current_key})) {
            $current_key .= "_";
            say "${def_key} is already defined in the current state, adding an _.";
        }
        $current{$current_key} = $defaults->{$def_key};
    }
    for my $step_key (keys %{$steps}) {
        my $current_key = $step_key;
        while (defined($current{$current_key})) {
            $current_key .= "_";
            say "${step_key} is already defined in the current state, adding an _.";
        }
        $current{$current_key} = $steps->{$step_key};
    }
    return(\%current);
}

## 202605 survey: Called in arbitrary functions to set the state before submitting a job.
sub Update_State {
    my ($class, %args) = @_;
    my %options;
    if (defined($args{options})) {
        %options = %{$args{options}};
    } else {
        my $function_info = Bio::Adventure::Config::Get_Function_Vars($class, \%args);
        %options = %{$function_info};
    }
    my %current_state = %{$class->current_state};

    ## FIXME:
    ## There are a couple of pieces of state set in Get_Function_Vars()
    ## in a way which means they will not get new interpolated variables
    ## freshly set variables, most notably output_prefix
    ## I think I should have output_prefix set in defaults.json
    ## and then changes to jprefix etc should carry over.
    ## Until I do that, set output_prefix here:
    $options{output_prefix} = qq"$current_state{output_basedir}/$options{jprefix}";
    $options{output_dir} = qq"$options{output_prefix}_$options{jname}";
    $options{time_file} = qq"$options{output_dir}/$options{jname}.time";
    if (-d $options{output_dir}) {
        ## print $state{output_dir} already exists.
    } else {
        make_path($options{output_dir});
    }

    ## Set stdout/stderr if they are not provided by the functions.json
    if (!$options{stdout}) {
        $options{stdout} = qq"$options{output_dir}/stdout";
    }
    if (!$options{stderr}) {
        $options{stderr} = qq"$options{output_dir}/stderr";
    }
    if ($options{main}) {
        print "TESTME: MAIN, main is set to $options{main}\n";
    }

    ## If prefix is not provided, set it to the number of inputs obseved so far.
    if (!$options{jprefix}) {
        my $num_inputs = scalar(@{$class->inputs});
        $options{jprefix} = $num_inputs;
    }
    ## Unfortunately, this cannot fix the problem, because the arbitrary output
    ## variables which are created all over the place have already interpolated
    ## the original jprefix containing string.
    ## I think therefore I need to have a single function responsible for
    ## performing the interpolations and making sure that it is run
    ## only in the update_state() or submit() function.
    ## In addition, any make_path() or similar secondary effect statements
    ## must but run at the same time.

    ## END FIXME:

    my $options_ref = $class->self_interp(hashref => \%options);
    ## Add portions of this to the current state
  OPTS: for my $k (keys %{$options_ref}) {
        next OPTS unless (defined($options_ref->{$k}));
        next OPTS if (defined($current_state{$k}));
        my $opt_class = ref($options_ref->{$k});
        next OPTS unless ($opt_class eq 'SCALAR');
        next OPTS if ($options_ref->{$k} eq 'undef');
        $current_state{$k} = $options_ref->{$k};
    }

    if ($options_ref->{template}) {
        my $jstring = q{};
        my $templater = Template->new({
            INCLUDE_PATH => $template_dir,
            TRIM => 0,
            PRE_CHOMP => 0,
            POST_CHOMP => 0,
            INTERPOLATE => 0,});
        print "TESTME: About to process template: $options_ref->{template}\n";

        ##my $processed = $tt->process($function_info->{template}, $args{template_args}, \$jstring) || die $tt->error();
        my $processed = $templater->process($options_ref->{template}, $options_ref, \$jstring) or die $templater->error();
        $current_state{jstring} = $jstring;
        print "TESTME: current_state{jstring}: $current_state{jstring}\n";
    } elsif ($args{jstring}) {
        print "Adding jstring.\n";
        $current_state{jstring} = $args{jstring};
    } else {
        print "Neither a template nor a jstring.\n";
    }

    if (!defined($current_state{language})) {
        $current_state{language} = 'bash';
    }

    ## I would much prefer to not need the following because
    ## if I do this, I cannot reuse the cyoa instance because this modifies
    ## the pre-interpolated variables.
    my @all_options = @{$class->options};
    push(@all_options, $options_ref);
    my $options_set = $class->options(\@all_options);
    my $set = $class->current_state(\%current_state);
    return($set);
}

1;
