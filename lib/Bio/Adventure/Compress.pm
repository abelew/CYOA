package Bio::Adventure::Compress;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename qw"basename dirname";
use File::Path qw"make_path";
use File::Which qw"which";

=head1 NAME

Bio::Adventure::Compress - Handle the (de/re)compression of data files.

=head1 SYNOPSIS

Some tools can handle compressed input, some cannot.  Bio::Adventure makes heavy use
of bash subshells <() to automagically handle any format, but at times
one must still decompress/recompress some data.

=head1 METHODS

=head2 C<Recompress>

Invoke xz to recompress a given input file.

=cut
sub Compress {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        comment => '## Compressing files.',
        jmem => 12,
        jname => 'xz',
        jprefix => '',
        jwalltime => '36:00:00',
        modules => undef,);
    my $paths = $class->Bio::Adventure::Config::Get_Paths();
    my $output_prefix = $paths->{output_prefix};
    my $output_dir = qq"${output_prefix}$options->{jname}";
    make_path($output_dir);
    my $stderr = qq"${output_dir}/xz.stderr";
    my $stdout = qq"${output_dir}/xz.stdout";
    my $input_paths = $class->Get_Path_Info($options->{input});
    my $check = which('xz');
    die('Could not find xz in your PATH.') unless($check);
    my $jstring = "";
    my $output_string = '';
  INPUTS: for my $in (@{$input_paths}) {
        my $in_dir = $in->{directory};
        my $in_base = $in->{filebase_compress};
        my $in_full = $in->{fullpath};
        my $output_file = qq"${in_dir}/${in_base}.xz";
        if (-r $in_full) {
            my $size = (stat $in_full)[7];
            ## Skip files less than 1Mb.
            if ($size < 1000000) {
                print "The file: ${in_full} is only ${size} bytes, skipping it.\n";
                next INPUTS;
            }
            if (-l $in_full) {
                print "The file: ${in_full} is a symlink, skipping it.\n";
                next INPUTS;
            }
        }
        $output_string .= qq"${output_file}:";
        $jstring .= qq!
## Compressing ${in_full}
echo "Compressing ${in_full}" 1>>${stdout}
compressed=1
if [ -f "${in_full}" ]; then
  {
    /usr/bin/du -h ${in_full} 1>>${stdout} 2>>${stderr}
    /usr/bin/time -a -o ${stdout} xz -9e -f ${in_full} 1>>${stdout} 2>>${stderr}
  } || {
    echo "Compression of ${in_full} failed." >> ${stdout}
    compressed=0
  }
else
  echo "The input: ${in_full} does not exist." 1>>${stderr}
fi
!;
    } ## End looping over input files.
    $output_string =~ s/:$//g;

    my $compression = $class->Submit(
        comment => $options->{comment},
        jdepends => $options->{jdepends},
        input => $options->{input},
        jcpu => 1,
        jmem => $options->{jmem},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => $options->{jwalltime},
        output => $output_string,);
    return($compression);
}

sub Recompress {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        comment => '## Recompressing files.',
        jname => 'xz',
        jmem => 8,
        jprefix => '99',
        jwalltime => '12:00:00',);
    my $input_paths = $class->Get_Path_Info($options->{input});

    my $jstring = "";
    my $output_string = '';
    for my $in (@{$input_paths}) {
        my $in_dir = $in->{directory};
        my $in_base = $in->{filebase_compress};
        my $in_full = $in->{fullpath};
        my $output_file = qq"${in_dir}/${in_base}.xz";
        my $fc = $class->Get_FC(input => $in_full);
        $output_string .= qq"${output_file}:";
        if ($in_full =~ /\.gz|\.bz2|\.zip/) {
            $jstring .= qq!
${fc} | xz -9e -f > ${in_dir}/${in_base}.xz
!;
        } else {
            $jstring .= qq!
xz -9e -f ${in_full}
!;
        }
    }
    $output_string =~ s/:$//g;

    my $compression = $class->Submit(
        comment => $options->{comment},
        jdepends => $options->{jdepends},
        input => $options->{input},
        jcpu => 1,
        jmem => $options->{jmem},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => $options->{jwalltime},
        output => $output_string,);
    return($compression);
}

=head2 C<Spring>

 Using spring's lossless mode to recompress fastq archives.

 https://doi.org/10.1093/bioinformatics/bty1015

 I checked out spring recently and was pleased with its speed.  I am not sure how I feel
 about lossy compression of the data; but the compression ratios were undeniably nice.

=cut
sub Spring {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        comment => '## Recompressing files.',
        jname => 'spring',
        jmem => 12,
        jprefix => '99',
        jwalltime => '12:00:00',);
    my $input_paths = $class->Get_Path_Info($options->{input});
    my $paths = $class->Bio::Adventure::Config::Get_Paths();

    my @compression_jobs = ();
    my $number_inputs = scalar(@{$input_paths});
    print "TESTME: Input number: $number_inputs\n";
    my $jstring = '';
    if ($number_inputs == 0) {
        die("No input provided.");
    } elsif ($number_inputs == 1) {
        print "Assuming a single, single-ended sample.\n";
        my $in = $input_paths->[0];
        my $in_dir = $in->{directory};
        my $in_base = $in->{filebase_compress};
        $in_base = basename($in_base, ('.fastq'));
        my $jname = qq"$options->{jname}_${in_base}";
        my $in_full = $in->{fullpath};
        if (-l $in_full) {
            print "The file: ${in_full} is a symlink, skipping it.\n";
            return undef;
        }
        my $size = (stat $in_full)[7];
        unless(defined($size)) {
            print "The size appears undefined for ${in_full}, skipping it.\n";
            return(undef);
        }
        ## Skip files less than 1Mb.
        if ($size < 1000000) {
            print "The file: ${in_full} is only ${size} bytes, skipping it.\n";
            return(undef);
        }
        my $output_file = qq"$paths->{output_dir}/${in_base}.spring";
        my $fc = '';
        my $gz_flag = '';
        if ($in_full =~ /.gz$/) {
            $gz_flag = ' -g';
            $fc = { input => $in_full };
        } else {
            $fc = $class->Get_FC(input => $in_full, needs_seek => 1);
        }
        if ($fc->{decompress}) {
            $jstring .= qq!$fc->{decompress}
!;
        }
        $jstring .= qq!
start=\$(du -sb ${in_full} | awk '{print \$1}')
/usr/bin/time -v -o $paths->{stdout}.time -a \\
  spring -c${gz_flag} -i $fc->{input} \\
    -o ${output_file} \\
    2>>$paths->{stderr} 1>>$paths->{stdout}
final=\$(du -sb ${output_file} | awk '{print \$1}')
ratio=\$(perl -e "print \${final} / \${start}")
echo "" >> $paths->{stdout}
echo "The input size is \${start}." >> $paths->{stdout}
echo "The output size of ${in_base} is: \${final}." >> $paths->{stdout}
echo "The ratio is: \${ratio}." >> $paths->{stdout}
mv ${output_file} ${in_dir}/
!;
        if ($fc->{cleanup}) {
            $jstring .= qq"$fc->{cleanup}
";
        } else {
            $jstring .= qq"rm ${in_full}
";
        }
        my $compression = $class->Submit(
            comment => $options->{comment},
            jdepends => $options->{jdepends},
            input => $options->{input},
            jcpu => 1,
            jmem => $options->{jmem},
            jname => $jname,
            jprefix => $options->{jprefix},
            jstring => $jstring,
            jwalltime => $options->{jwalltime},
            output => $output_file,);
        push(@compression_jobs, $compression);
    } elsif ($number_inputs == 2) {
        print "Assuming a single, paired-ended sample.\n";
        my $r1 = $input_paths->[0];
        my $r2 = $input_paths->[1];
        my $in_dir = $r1->{directory};
        my $in_base = $r1->{filebase_compress};
        $in_base = basename($in_base, ('.fastq'));
        my $r1_full = $r1->{fullpath};
        my $r2_full = $r2->{fullpath};
        if (-l $r1_full || -l $r2_full) {
            print "The file: ${r1_full} or ${r2_full} is a symlink, skipping it.\n";
            return undef;
        }
        my $r1_size = (stat $r1_full)[7];
        my $r2_size = (stat $r2_full)[7];
        ## Skip files less than 1Mb.
        if (($r1_size + $r2_size) < 2000000) {
            print "The file: ${r1_full} is only ${r1_size} bytes, skipping it.\n";
            next INPUTS;
        }
        my $jname = qq"$options->{jname}_${in_base}";
        my $output_file = qq"$paths->{output_dir}/${in_base}.spring";
        my $r1_fc = '';
        my $r2_fc = '';
        my $gz_flag = '';
        if ($r1_full =~ /.gz$/) {
            $gz_flag = ' -g';
            $r1_fc = { input => $r1_full };
            $r2_fc = { input => $r2_full };
        } else {
            $r1_fc = $class->Get_FC(input => $r1_full, needs_seek => 1);
            $r2_fc = $class->Get_FC(input => $r2_full, needs_seek => 1);
        }

        if ($r1_fc->{decompress}) {
            $jstring .= qq!$r1_fc->{decompress}
$r2_fc->{decompress}
!;
        }
        $jstring .= qq!
start=\$(( \$(du -sb ${r1_full} | awk '{print \$1}' | perl -pe 's/[a-zA-Z]\$//g') + \$(du -sb ${r1_full} | awk '{print \$1}' | perl -pe 's/[a-zA-Z]\$//g') ))
/usr/bin/time -v -o $paths->{stdout}.time -a \\
  spring -c${gz_flag} -i $r1_fc->{input} $r2_fc->{input} \\
    -o ${output_file} \\
    2>>$paths->{stderr} 1>>$paths->{stdout}
final=\$(du -sb ${output_file} | awk '{print \$1}')
echo "" >> $paths->{stdout}
echo "The input size is \${start}." >> $paths->{stdout}
echo "The output is: \${final}." >> $paths->{stdout}
echo "The ratio is: \$(perl -e \\"print \${final} / \${start}\\")." >> $paths->{stdout}"
mv ${output_file} ${in_dir}/
!;
        if ($r1_fc->{cleanup}) {
            $jstring .= qq"$r1_fc->{cleanup}
$r2_fc->{cleanup}
";
        } else {
            $jstring .= qq"rm ${r1_full}
rm ${r2_full}
";
        }
        my $compression = $class->Submit(
            comment => $options->{comment},
            jdepends => $options->{jdepends},
            input => $options->{input},
            jcpu => 1,
            jmem => $options->{jmem},
            jname => $jname,
            jprefix => $options->{jprefix},
            jstring => $jstring,
            jwalltime => $options->{jwalltime},
            output => $output_file,);
         push(@compression_jobs, $compression);
    } else {
        print "Assuming multiple single-ended samples.\n";
      INPUTS: for my $in (@{$input_paths}) {
            my $in_dir = $in->{directory};
            my $in_base = $in->{filebase_compress};
            $in_base = basename($in_base, ('.fastq'));
            my $in_full = $in->{fullpath};
            my $size = (stat $in_full)[7];
            ## Skip files less than 1Mb.
            if ($size < 1000000) {
                print "The file: ${in_full} is only ${size} bytes, skipping it.\n";
                next INPUTS;
            }
            if (-l $in_full) {
                print "The file: ${in_full} is a symlink, skipping it.\n";
                next INPUTS;
            }
            $in_base = basename($in_base, ('.fastq'));
            my $jname = qq"$options->{jname}_${in_base}";
            my $output_file = qq"$paths->{output_dir}/${in_base}.spring";
            my $fc = '';
            my $gz_flag = '';
            if ($in_full =~ /.gz$/) {
                $gz_flag = ' -g';
                $fc = { input => $in_full };
            } else {
                $fc = $class->Get_FC(input => $in_full, needs_seek => 1);
            }
            if ($fc->{decompress}) {
                $jstring .= qq!$fc->{decompress}
!;
            }
            $jstring .= qq!
start=\$(du -sb ${in_full} | awk '{print \$1}')
/usr/bin/time -v -o $paths->{stdout}.time -a \\
  spring -c${gz_flag} -i $fc->{input} \\
    -o ${output_file} \\
    2>>$paths->{stderr} 1>>$paths->{stdout}
final=\$(du -sb ${output_file} | awk '{print \$1}')
echo "" >> $paths->{stdout}
echo "The input size is \${start}." >> $paths->{stdout}
echo "The output size is \${final}." >> $paths->{stdout}
echo "The ratio is: \$(perl -e \\"print \${final} / \${start}\\")." >> $paths->{stdout}
mv ${output_file} ${in_dir}/
!;
            if ($fc->{cleanup}) {
                $jstring .= qq"$fc->{cleanup}
";
            } else {
                $jstring .= qq"rm ${in_full}
";
            }
            my $compression = $class->Submit(
                comment => $options->{comment},
                jdepends => $options->{jdepends},
                input => $options->{input},
                jcpu => 1,
                jmem => $options->{jmem},
                jname => $jname,
                jprefix => $options->{jprefix},
                jstring => $jstring,
                jwalltime => $options->{jwalltime},
                output => $output_file,);
             push(@compression_jobs, $compression);
        } ## End iterating over input files.
    }
    return(\@compression_jobs);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<xz> L<gzip> L<bash>

=cut

1;
