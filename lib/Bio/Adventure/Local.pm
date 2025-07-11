package Bio::Adventure::Local;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use feature 'try';
no warnings 'experimental::try';

use Cwd;
use File::Basename qw "basename dirname";
use File::Path qw"make_path remove_tree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use File::Which qw"which";
use IO::Handle;
use Template;
my $template_base = dist_dir('Bio-Adventure');
my $template_dir = qq"${template_base}/templates";

=head1 NAME

Bio::Adventure::Local - Write out jobs using scripts for running jobs on the local computer.

=head1 SYNOPSIS

use Bio::Adventure;
my $hpgl = new Bio::Adventure::Local;
$hpgl->Submit(%job_args);

This should invoke scripts with appropriate parameters for various jobs on the local host.

=head1 Methods

=head2 C<Submit>

This is responsible for the final writeup and submission of a local script.

=cut
sub Submit {
    my ($class, $parent, %args) = @_;
    my $options = $parent->Get_Vars(
        args => \%args,
        jname => 'unknown',);
    ## For arguments to bash, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    my $bash_log = qq"$options->{basedir}/outputs/log.txt";
    if (!defined($options->{jname})) {
        $options->{jname} = $class->Get_Job_Name();
    }
    my $jname = qq"$options->{jprefix}$options->{jname}";
    if (!defined($jname)) {
        print "In Local::Submit, jname was not defined, this should not happen.\n";
        $jname = 'undefined';
    }
    my $finished_file = qq"$options->{basedir}/outputs/logs/${jname}.finished";
    my $job = {};
    foreach my $k (keys %args) {
        next if ($k eq 'jstring');
        next if ($k eq 'comment');
        $job->{$k} = $args{$k};
    }
    my @wanted_vars = ('basedir', 'depends_string', 'input',
                       'jcpu', 'jmem', 'jname', 'jqueue', 'jwalltime', 'output');
    foreach my $w (@wanted_vars) {
        $job->{$w} = $options->{$w} if (!defined($job->{$w}));
    }
    if ($options->{restart} && -e $finished_file) {
        print "The restart option is on, and this job appears to have finished.\n";
        return($job);
    }
    my $script_file = qq"$options->{basedir}/$options->{script_dir}/$options->{jprefix}$options->{jname}.sh";
    my $mycwd = getcwd();
    my $out_dir;
    if (!defined($options->{stdout}) && !defined($options->{output})) {
        die("Every job must have either output or stdout defined.");
    } elsif (!defined($options->{stdout})) {
        my $subroutine = Bio::Adventure::Config::Get_Caller();
        $out_dir = dirname($options->{output});
        warn("Every job should have a stdout defined, " .
             "setting it to the output directory: ${out_dir} for ${subroutine}.");
        $options->{stdout} = $out_dir;
        sleep(5);
    } elsif (!defined($options->{output})) {
        my $subroutine = Bio::Adventure::Config::Get_Caller();
        warn("Every job should have an output defined, " .
             "setting it to $options->{stdout} for ${subroutine}.");
        $options->{output} = $options->{stdout};
        $out_dir = dirname($options->{stdout});
        sleep(5);
    } else {
        $out_dir = dirname($options->{stdout});
    }

    make_path($out_dir, {verbose => 0}) unless (-r $out_dir);
    make_path("$options->{logdir}", {verbose => 0}) unless (-r qq"$options->{logdir}");
    unless (-r qq"$options->{basedir}/scripts") {
        make_path("$options->{basedir}/$options->{script_dir}/pdata", {verbose => 0});
    }
    my $script_base = basename($script_file);

    ## Remove the need for two functions that do the same thing except one for perl and one for bash
    if ($options->{language} eq 'perl') {
        my $perl_file = qq"$options->{basedir}/$options->{script_dir}/$options->{jprefix}$options->{jname}.pl";
        my $perl_stderr = qq"${out_dir}/$options->{jprefix}$options->{jname}.stderr";
        my $perl_stdout = qq"${out_dir}/$options->{jprefix}$options->{jname}.stdout";
        my $perl_start = qq?#!/usr/bin/env perl
use strict;
use FileHandle;
use Bio::Adventure;
my \$out = FileHandle->new(">>${bash_log}");
my \$d = qx'date';
print \$out "###Started ${script_base} at \${d}";
chdir("$options->{basedir}");
my \$h = Bio::Adventure->new();
?;
        if (defined($parent->{option_file})) {
            $perl_start .= qq!
if (-r "$parent->{option_file}") {
  ## Pull options from the option file: $parent->{option_file}
  my \$loaded = \$h->Load_Vars(input => '$parent->{option_file}');
}
!;
        }
        my $perl_end = qq!## The following lines give status codes and some logging
my \$jobid = "";
\$jobid = \$ENV{SLURM_JOBID} if (\$ENV{SLURM_JOBID});
my \$end_d = qx'date';
print \$out "####Finished \${jobid} ${script_base} at \${end_d}.";
close(\$out);
!;
        my $total_perl_string = qq"${perl_start}\n";
        $total_perl_string .= qq"$options->{comment}\n" if ($options->{comment});
        $total_perl_string .= qq"$options->{prescript}\n" if ($options->{prescript});
        $total_perl_string .= qq"$options->{jstring}\n" if ($options->{jstring});
        $total_perl_string .= qq"${perl_end}\n";

        my $perl_script = FileHandle->new(">${perl_file}");
        print $perl_script $total_perl_string;
        $perl_script->close();
        chmod(0755, $perl_file);
        $options->{jstring} = qq"
${perl_file} \\
  2>${perl_stderr} \\
  1>${perl_stdout}
";
    } ## End extra processing for submission of a perl script (perhaps not needed for slurm?

    if ($options->{jtemplate}) {
        print "In local, processing with dir: ${template_dir}\n";
        print "Template is outputting to: ${script_file}\n";
        my $tt = Template->new({
            INCLUDE_PATH => $template_dir,
            OUTPUT => $script_file,
            TRIM => 1,
            INTERPOLATE => 1,});
        $tt->process($options->{jtemplate}, $options) || die $tt->error();
    } else {
        my $script_start = qq?#!$options->{shell}
?;
        $script_start .= qq"set -x\n" if ($options->{debug});
        $script_start .= qq?cd $options->{basedir}
set -o errexit
set -o errtrace
set -o pipefail
echo "## Started ${script_base} at \$(date) on \$(hostname)." >> ${bash_log}
function get_sigterm {
  echo "A SIGTERM was sent to ${jname}." >> ${bash_log}
  exit 1
}
trap get_sigterm SIGTERM
function get_sigerr {
  echo "A ERR was sent to ${jname}." >> ${bash_log}
  exit 1
}
trap get_sigerr ERR
?;
        $script_start .= $options->{module_string} if ($options->{module_string});
        $script_start .= $options->{conda_string} if ($options->{conda_string});
        my $script_end = qq!
## The following lines give status codes and some logging
echo "Job status:\$?" >> ${bash_log}
echo " \$(hostname) Finished ${script_base} at \$(date), it took \$(( SECONDS / 60 )) minutes." >> ${bash_log}
!;

        my $total_script_string = '';
        $total_script_string .= qq"${script_start}\n";
        $total_script_string .= qq"$options->{comment}\n" if ($options->{comment});
        $total_script_string .= qq"$options->{prescript}\n" if ($options->{prescript});
        $total_script_string .= qq"$options->{jstring}\n" if ($options->{jstring});
        $total_script_string .= qq"$options->{postscript}\n" if ($options->{postscript});
        $total_script_string .= qq"${script_end}\n";
        my $script = FileHandle->new(">${script_file}");
        if (!defined($script)) {
            die("Could not write the script: ${script_file}, check its permissions.")
        }
        print $script $total_script_string;
        $script->close();
    }
    chmod(0755, $script_file);
    my $job_text = '';
    my $handle;
    my $bash_pid = open($handle, '-|', ${script_file}) or
        die("The script: ${script_file}
failed with error: $!.\n");
    print "Starting a bash job: ${bash_pid} $options->{jname}";
    if ($options->{jdepends}) {
        print ", depending on $options->{jdepends}.";
    }
    print "\n";
    while (my $line = <$handle>) {
        $job_text = $job_text . $line;
        print $line if ($options->{debug});
    }
    my $closed;
    try {
        $closed = close($handle);
    } catch ($e) {
        warn "Unable to close script filehandle return: $? error: $!\n";
    }
    print "Finished running, outputs should be in $options->{output}.\n\n";
    if (!defined($bash_pid)) {
        print "It should be impossible for bash_pid to be undefined, yet here we are.\n";
        $bash_pid = 'undefined';
    }

    $job->{log} = $bash_log;
    $job->{job_id} = $bash_pid;
    $job->{script_file} = $script_file;
    my $job_logger = FileHandle->new(">>outputs/logs/jobs.txt");
    print $job_logger "${jname}\tbash\t${bash_pid}\n";
    $job_logger->close();
    ## Reset the environment in case we left any cruft behind
    my $reset = Bio::Adventure::Reset_Vars($class);
    return($job);
}

sub Wait {
    return(undef);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Adventure::Slurm> L<Bio::Adventure::Torque>

=cut

1;
