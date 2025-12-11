package Bio::Adventure::Local;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
use feature 'try';
no warnings 'experimental::try';

use Bio::Adventure::Languages;
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
    my ($runner, $class, %args) = @_;
    my $spec = $class->spec;
    my $state = $class->current_state;
    my $caller = $class->caller;
    my $all_options = $class->options;
    my $job_options = $all_options->[$#$all_options];
    my $primary_jstring = $class->Bio::Adventure::Languages::Make_Sub_Job();

    my $comment_string = Bio::Adventure::Languages::Make_Comment(
        language => $spec->{language}, jcomment => $spec->{jcomment});
    ## For arguments to bash, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    my $bash_log = qq"$state->{basedir}/outputs/log.txt";
    if (!defined($spec->{jname})) {
        print "In Local::Submit, jname was not defined, this should not happen.\n";
        $spec->{jname} = 'undefined';
    }
    my $finished_file = qq"$state->{basedir}/outputs/logs/$spec->{jname}.finished";
    my $job = {};
    for my $k (keys %args) {
        next if ($k eq 'jstring');
        next if ($k eq 'comment');
        $job->{$k} = $args{$k};
    }

    if ($spec->{restart} && -e $finished_file) {
        print "The restart option is on, and this job appears to have finished.\n";
        return($job);
    }
    my $script_file = qq"$state->{basedir}/$state->{script_dir}/$spec->{jname}.sh";
    my $mycwd = getcwd();
    make_path("$state->{logdir}", {verbose => 0}) unless (-r qq"$state->{logdir}");
    unless (-r qq"$state->{basedir}/$state->{script_dir}/json") {
        make_path("$state->{basedir}/$state->{script_dir}/json", {verbose => 0});
    }
    print "Writing job script to $script_file\n";
    my $script_base = basename($script_file);
    ## Remove the need for two functions that do the same thing except one for perl and one for bash
    if ($spec->{jtemplate}) {
        print "In local, processing with dir: ${template_dir}\n";
        print "Template is outputting to: ${script_file}\n";
        my $tt = Template->new({
            INCLUDE_PATH => $template_dir,
            OUTPUT => $script_file,
            TRIM => 1,
            INTERPOLATE => 1,});
        $tt->process($spec->{jtemplate}, $state) || die $tt->error();
    } else {
        my $script_start = qq?$state->{shebang}
?;
        $script_start .= qq"set -x\n" if ($state->{debug});
        $script_start .= qq?cd $state->{basedir} || exit
set -o errexit
set -o errtrace
set -o pipefail
echo "## Started ${script_base} at \$(date) on \$(hostname)." >> ${bash_log}
function get_sigterm {
  echo "A SIGTERM was sent to $spec->{jname}." >> ${bash_log}
  exit 1
}
trap get_sigterm SIGTERM
function get_sigerr {
  echo "A ERR was sent to $spec->{jname}." >> ${bash_log}
  exit 1
}
trap get_sigerr ERR
?;
        $script_start .= $state->{module_string} if ($state->{module_string});
        $script_start .= $state->{conda_string} if ($state->{conda_string});
        my $script_end = qq!
## The following lines give status codes and some logging
echo "Job status:\$?" >> ${bash_log}
echo " \$(hostname) Finished ${script_base} at \$(date), it took \$(( SECONDS / 60 )) minutes." >> ${bash_log}
!;
        if (!defined($comment_string)) {
            print "The comment string for ${caller} is undefined.\n";
            $comment_string = '## currently undefined';
        }
        ## This and the following check on postscript should not be needed.
        if (!defined($spec->{prescript})) {
            ##print "The prescript for ${caller} is undefined.\n";
            $spec->{prescript} = '';
        }
        if (!defined($spec->{postscript})) {
            ##print "The postscript for ${caller} is undefined.\n";
            $spec->{postscript} = '';
        }
        my $total_script_string = qq?${script_start}
${comment_string}
$spec->{prescript}
${primary_jstring}
$spec->{postscript}
${script_end}
?;
        my $script_fh = FileHandle->new(">${script_file}");
        if (!defined($script_fh)) {
            die("Could not write the script: ${script_file}, check its permissions.")
        }
        print $script_fh $total_script_string;
        $script_fh->flush();
        $script_fh->close();
    }
    chmod(0755, $script_file);
    my $job_text = '';
    my $handle;
    my $bash_pid = open($handle, '-|', ${script_file}) or
        die("The script: ${script_file}
failed with error: $!.\n");
    print "Starting a bash job: ${bash_pid} $spec->{jname}";
    if ($spec->{jdepends}) {
        print ", depending on $spec->{jdepends}.";
    }
    print "\n";
    while (my $line = <$handle>) {
        $job_text = $job_text . $line;
        print $line if ($state->{debug});
    }
    my $closed;
    try {
        $closed = close($handle);
    } catch ($e) {
        warn "Unable to close script filehandle return: $? error: $!\n";
    }
    if (!defined($bash_pid)) {
        print "It should be impossible for bash_pid to be undefined, yet here we are.\n";
        $bash_pid = 'undefined';
    }

    $job->{log} = $bash_log;
    $job->{job_id} = $bash_pid;
    $job->{script_file} = $script_file;
    $class->{language} = 'bash';
    my $logging_dir = qq"$state->{basedir}/outputs/logs";
    make_path($logging_dir);
    my $logging_file = qq"${logging_dir}/jobs.txt";
    my $job_logger = FileHandle->new(">>${logging_file}");
    print $job_logger "$spec->{jname}\tbash\t${bash_pid}\n";
    $job_logger->close();
    ## Reset the environment in case we left any cruft behind
    ## my $reset = Bio::Adventure::Reset_Vars($class);
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
