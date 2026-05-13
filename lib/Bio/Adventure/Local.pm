package Bio::Adventure::Bash;
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
my $template_dir = qq"${template_base}/jobs";

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
    my $state_ref = $class->current_state;
    my %state = %{$state_ref};
    my $caller = $class->caller;
    my $primary_jstring = $class->Bio::Adventure::Languages::Make_Sub_Job();

    ##########################
    ## TODO TODO TODO
    ##########################
    ## Replace this foolishness with a template ASAP.

    my $comment_string = Bio::Adventure::Languages::Make_Comment(
        language => $state{language}, jcomment => $state{jcomment});
    ## For arguments to bash, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    $state{bash_log} = qq"$state{basedir}/outputs/log.txt";
    if (!defined($state{jname})) {
        $state{jname} = 'undefined';
    }
    my $finished_file = qq"$state{basedir}/outputs/logs/$state{jname}.finished";
    my $job = {};
    for my $k (keys %args) {
        next if ($k eq 'jstring');
        next if ($k eq 'comment');
        $job->{$k} = $args{$k};
    }
    my $script_file = $class->Get_Script_Filename();
    my $mycwd = getcwd();
    make_path("$state{logdir}", {verbose => 0}) unless (-r qq"$state{logdir}");
    unless (-r qq"$state{basedir}/$state{script_dir}/json") {
        make_path("$state{basedir}/$state{script_dir}/json", {verbose => 0});
    }
    $state{debug_string} = q{};
    if ($state{debug}) {
        $state{debug_string} = 'set -x';
    }
    $state{prescript} = q{};
    if (defined($state{prescript})) {
        $state{prescript} = $state{prescript};
    }
    $state{postscript} = q{};
    if (defined($state{postscript})) {
        $state{postscript} = $state{postscript};
    }
    $state{primary_jstring} = $primary_jstring;
    $state{script_base} = basename($script_file);
    my $templater = Template->new({
        INCLUDE_PATH => $template_dir,
        TRIM => 0,
        INTERPOLATE => 0,});
    my $processed = $templater->process('local.tt', \%state, $script_file);
    print "Writing job script to $script_file\n";
    chmod(0755, $script_file);
    my $job_text = '';
    my $handle;
    my $bash_pid = open($handle, '-|', ${script_file}) or
        die("The script: ${script_file}
failed with error: $!.\n");
    print "Starting a bash job: ${bash_pid} $state{jname}.\n";
    while (my $line = <$handle>) {
        $job_text = $job_text . $line;
        print $line if ($state{debug});
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
    $job->{log} = $state{bash_log};
    $job->{job_id} = $bash_pid;
    $job->{script_file} = $script_file;
    for my $k (keys %state) {
        $job->{$k} = $state{$k};
    }
    my $logging_dir = qq"$state{basedir}/outputs/logs";
    make_path($logging_dir);
    my $logging_file = qq"${logging_dir}/jobs.txt";
    my $job_logger = FileHandle->new(">>${logging_file}");
    print $job_logger "$state{jname}\tbash\t${bash_pid}\n";
    $job_logger->close();
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
