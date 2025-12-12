package Bio::Adventure::Languages;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use FileHandle;
use File::Path qw"make_path remove_tree";
use JSON qw"to_json from_json";

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

sub Make_Sub_Job {
    my ($class, %args) = @_;
    my $options = $class->options;
    my $spec = $class->spec;

    my $primary_job = $spec->{jstring};
    if (defined($spec->{subjob})) {
        if ($spec->{subjob} eq 'perl') {
            my $perl_file = qq"$options->{basedir}/$options->{script_dir}/$options->{jname}.pl";
            my $perl_stderr = qq"$options->{output_dir}/$options->{jname}.stderr";
            my $perl_stdout = qq"$options->{output_dir}/$options->{jname}.stdout";
            my $bash_log = qq"$options->{basedir}/outputs/log.txt";
            $options->{prescript} = '' unless (defined($options->{prescript}));
            $options->{postscript} = '' unless (defined($options->{postscript}));
            my $subjob = qq?#!/usr/bin/env perl
use strict;
use Modern::Perl;
use FileHandle;
use Bio::Adventure;
my \$out = FileHandle->new(">>${bash_log}");
my \$d = qx'date';
print \$out "###Started $options->{jname} at \${d}";
chdir("$options->{basedir}");
my \$h = Bio::Adventure->new();
my \$information = \$h->Bio::Adventure::Languages::Get_JSON('$class->{option_file}');
## $options->{jcomment}
$options->{prescript}
$options->{jstring}
$options->{postscript}
## The following lines give status codes and some logging
my \$jobid = "";
\$jobid = \$ENV{SLURM_JOBID} if (\$ENV{SLURM_JOBID});
my \$end_d = qx'date';
print \$out "####Finished \${jobid} $options->{jname} at \${end_d}.";
close(\$out);
?;
            my $perl_script = FileHandle->new(">${perl_file}");
            print $perl_script $subjob;
            $perl_script->close();
            chmod(0755, $perl_file);
            $primary_job = qq?
${perl_file} \\
  2>${perl_stderr} \\
  2>${perl_stdout}
?;
        } else {
            print "I currently only know how to set up a perl job.\n";
        }
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
