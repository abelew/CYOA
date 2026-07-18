package Bio::Adventure::Languages;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use FileHandle;
use File::Path qw"make_path remove_tree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use JSON qw"to_json from_json";
my $template_base = dist_dir('Bio-Adventure');
my $template_dir = qq"${template_base}/jobs";

sub Get_JSON {
    my ($class, $option_file) = @_;
    my $input_string = '';
    my $in = FileHandle->new(qq"<${option_file}");
    while (my $line = <$in>) {
        $input_string .= $line;
    }
    $in->close();
    my $information = from_json($input_string);
    my $class_opts = $information->{class};
    for my $k (keys %{$class_opts}) {
        $class->{$k} = $class_opts->{$k};
    }
    $class->{options} = $information->{options};
    $class->{paths} = $information->{paths};
    $class->{inputs} = $information->{inputs};
    $class->{outputs} = $information->{outputs};
    return($information);
}

sub Make_Comment {
    my %args = @_;
    my $comment_string = '';
    if (!defined($args{jcomment})) {
        ## I am moving these comments to the job templates
        return(q{});
    }
    if ($args{language} =~ /lisp/) {
        $comment_string = qq";; $args{jcomment}";
    } elsif ($args{language} eq 'c' || $args{language} eq 'c++' ||
             $args{language} eq 'java') {
        $comment_string = qq"/*
     $args{jcomment}
 */";
    } else {
        ## Assume bash-ish languages
        $comment_string = qq"## $args{jcomment}";
    }
    return($comment_string);
}

=head2 C<Make_Sub_Job>

 Create a sub-job: This is a secondary script run by an arbitrary scripting language
 and submitted as a job.

 As of the time of this writing, I have only played with perl and R in this context;
 but with the recent changes to cyoa, it should be pretty easy to define jobs which
 run in other programming languages and get them to run with whatever
 parameters one wants.




=cut
sub Make_Sub_Job {
    my ($class, %args) = @_;
    my $state_ref = $class->current_state;
    my %state = %{$state_ref};
    print "TESTME: CREATING A SUBJIB with $state{perl_job}\n";
    my $primary_job = $state{jstring};
    if (defined($state{perl_job})) {
        ## perl_job is the location of the .pl which will be run by the primary job.
        ## The following template process() requires that a
        ## perl_jstring be defined; when I start creating ruby/R/python/julia/etc
        ## the same should be true for them and be comprised of the functions
        ## run by the appropriate runtime.
        $state{perl_stderr} = qq"$state{output_dir}/$state{jname}.stderr";
        $state{perl_stdout} = qq"$state{output_dir}/$state{jname}.stdout";

        ## Ensure a few assumed variables exist.
        if (!defined($state{log})) {
            $state{log} = qq"$state{basedir}/outputs/log.txt";
        }
        if (!defined($state{perl_shebang})) {
            $state{perl_shebang} = qq"#!/usr/bin/env perl";
        }

        ## The submitting function is responsible for providing a perl_job as above and
        ## a 'worker' which is the name of the function which will be called; there should
        ## be a corresponding file in the jobs/ tree which templates
        ## the various parameters which this job will use.
        my $worker_name = lc($state{worker});
        my $worker_file = qq"${worker_name}.tt";
        my $templater = Template->new({
            INCLUDE_PATH => $template_dir,
            TRIM => 0,
            INTERPOLATE => 0,});
        my $perl_jstring = q{};
        my $perl_jstring_processed = $templater->process($worker_file, \%state, \$perl_jstring);
        $state{perl_jstring} = $perl_jstring;
        print "\n\nTESTME In subjob: $perl_jstring\n\n";
        sleep(3);
        my $processed = $templater->process('perl.tt', \%state, $state{perl_job}) or die $templater->error;
        say "Writing job script to $state{perl_job}.";
        chmod(0755, $state{perl_job});
    } ## End checking for subjob.
    return($primary_job);
}

sub Save_JSON {
    my ($class, $options, $paths, $inputs, $outputs) = @_;
    my $template_dir = qq"$options->{basedir}/$options->{script_dir}/json";
    make_path($template_dir) unless (-r $template_dir);
    my $option_file = File::Temp->new(
        TEMPLATE => qq"$options->{jname}XXXX",
        DIR => $template_dir,
        UNLINK => 0,
        SUFFIX => '.json',);

    my $class_args = {};
    CLASS_ARGS: for my $k (keys %{$class}) {
        $class_args->{k} = $class->{$k};
    }
    my $option_args = {};
  OPT_ARGS: for my $k (keys %{$options}) {
        next OPT_ARGS if ($k eq 'variable_getopt_overrides');
        $option_args->{$k} = $options->{$k};
    }
    my $information = {
        class => $class_args,
        options => $option_args,
        paths => $paths,
        inputs => $inputs,
        outputs => $outputs,
    };
    my $json_text = to_json($information, {utf8 => 1, pretty => 1});
    ## my $json_text = to_json($information, {utf8 => 1});
    my $out = FileHandle->new(qq">${option_file}");
    print $out $json_text;
    $out->close();
    $class->{option_file} = $option_file;
    return($option_file);
}

1;
