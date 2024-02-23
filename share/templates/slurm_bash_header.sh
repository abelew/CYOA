#SBATCH --export=ALL --requeue --mail-type=NONE --open-mode=append
#SBATCH --chdir=[% basedir %] --nodes=1 --ntasks=1 --output=outputs/log.txt.sbatch
[% IF jname.defined %]
#SBATCH --job-name=[% jname %]
[% END %]
[% IF jnice.defined %]
#SBATCH --nice=[% jnice %]
[% ELSE %]
#SBATCH --nice=0
[% END %]
[% IF jcpu.defined %]
#SBATCH --cpus-per-task=[% jcpu %]
[% END %]
[% IF jmem.defined %]
#SBATCH --mem=[% jmem %]G
[% END %]
[% IF walltime.defined %]
#SBATCH --time=[% jwalltime %]
[% END %]
