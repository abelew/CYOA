package Bio::Adventure::Splicing;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use feature 'try';
no warnings 'experimental::try';

use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Capture::Tiny qw":all";
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Spec;
use File::Path qw"make_path rmtree";
use File::Which qw"which";
use File::ShareDir qw":ALL";
use IPC::Open2;
use Spreadsheet::Read;
use Symbol qw"gensym";


=head2 C<Suppa>

  Set up and invoke Suppa, a differential transcript and splicing
  event calculator.  10.1186/s13059-018-1417-1

  I am putting this function in this file because it is relevant to
  the gene structure, which is admittedly a bit of a stretch.

  This function should be able to invoke the tool which creates the
  catalog of splicing events (suppa.py generateEvents), create the
  table of TPM values by transcript.  I am thinking to use the sample
  sheet as the input variable with a parameter containing the
  salmon/etc filename and have this generate the associated table.
  Given that information, this should be able to run the psiPerIsoform
  function and/or psiPerEvent.  Note, it turns out that suppa provides
  a helper function for creating the expression table 'joinFiles'; so
  that will simplify things quite a bit.

=cut
sub Suppa {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        type => 'SE',
        condition_column => 'drug',
        file_column => 'hg38100salmonfile',
        jprefix => '90',);
    my $suppa_dir = qq"outputs/$options->{jprefix}suppa_$options->{species}";
    my $gff_file = qq"$options->{libpath}/$options->{libtype}/$options->{species}.gtf";
    my $reader = Spreadsheet::Read->new($options->{input});
    my $sheet = $reader->sheet(1);
    ## Strangely, Spreadsheet::Read only uses numeric values for columns, so grab the first row
    ## and get the number of my column from it...
    my @column_names = $sheet->cellrow(1);
    print "TESTME: @column_names\n";
    my $file_number = undef;
    my $condition_number = undef;
    my $count = 0;
  COLUMNS: for my $name (@column_names) {
      $count++;
      if ($name eq $options->{file_column}) {
          print "Found $options->{file_column} as number: $count\n";
          $file_number = $count;
      }
      if ($name eq $options->{condition_column}) {
          print "Found $options->{condition_column} as number: $count\n";
          $condition_number = $count;
      }
  }
    my @filenames = $sheet->cellcolumn($file_number);
    my @conditions = $sheet->cellcolumn($condition_number);
    my $files_by_condition = {};
    my $cond_count = -1;
    for my $cond (@conditions) {
        print "Working on ${cond}\n";
        $cond_count++;
        my @cond_filenames;
        if (defined($files_by_condition->{$cond})) {
            @cond_filenames = @{$files_by_condition->{$cond}};
        } else {
            @cond_filenames = ();
        }
        my $this_file = $filenames[$cond_count];
        push(@cond_filenames, $this_file);
        print "TESTME: @cond_filenames\n";
        $files_by_condition->{$cond} = \@cond_filenames;
    }
    my $join_string = '';
    my $psi_string = '';
    my $psi_local_string = '';
    print "Creating joinfiles, psiperisoform, and psiperevent strings.\n";
    for my $cond (sort keys %{$files_by_condition}) {
        my $file_string = '';
        for my $f (@{$files_by_condition->{$cond}}) {
            $file_string .= qq" ${f} ";
        }
        $join_string .= qq"suppa.py joinFiles -f tpm -i ${file_string} -o ${cond}.tpm
";
        print "TESTME: $join_string\n";
        $psi_string .= qq"suppa.py psiPerIsoform -g ${gff_file} -e ${cond}.tpm -o ${cond}.psi
";
        $psi_local_string .= qq"suppa.py psiPerEvent --ioe-file local_as_events.txt --expression-file ${cond}.tpm -o ${cond}_local.psi
";
    }

    my $file_string = '';
    for my $f (@filenames) {
        $file_string .= " ${f} ";
    }
    my $jstring = qq!mkdir -p ${suppa_dir}
suppa.py generateEvents -i ${gff_file} -f ioe -e SE -o local_as_events_SE.txt
suppa.py generateEvents -i ${gff_file} -f ioe -e SS -o local_as_events_SS.txt
suppa.py generateEvents -i ${gff_file} -f ioe -e MX -o local_as_events_MX.txt
suppa.py generateEvents -i ${gff_file} -f ioe -e RI -o local_as_events_RI.txt
suppa.py generateEvents -i ${gff_file} -f ioe -e FL -o local_as_events_FL.txt
suppa.py generateEvents -i ${gff_file} -f ioi -o transcript_events.txt
${join_string}
${psi_string}
${psi_local_string}
!;
    my $suppa_job = $class->Submit(
        jstring => $jstring,
        );
    return($suppa_job);
}


1;
