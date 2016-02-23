#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# all_run_PBS.pl
######################################

use strict;
use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;


my ($cmd, $helpFlag, $noGapFlag, $minRate_indelReplace);
my ($lib1, $lib2, $refFn, $contigFn, $contigFn2, $outHeader, $step, $iter_num);

$step = 1;
$iter_num = 5;
GetOptions(
	"h|?|help"	=> \$helpFlag,

	"out=s"	=> \$outHeader,
	"ref=s"	=> \$refFn,
	"c=s"	=> \$contigFn,
	"c2=s"	=> \$contigFn2,

	"lib1=s"	=> \$lib1,
	"lib2=s"	=> \$lib2,
	"num=s" => \$iter_num,

	"step=s"	=> \$step,
) || help(1);


help(0) if defined $helpFlag;


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 -o <header of output>  -r <reference FASTA file> -c <contig FASTA file> [-c2 <secondcontig FASTA file>] -lib1 <fastq files> [-lib2 <fastq files>]  [-n <number of iteration>] [-step <1,2,3,4,5,6 or 7>]\n\n";
	print STDERR "  ex 1) $0 -o myHeader -r ref.fa -c contigs.fa -lib1 'lib1_1.fq lib1_2.fq'\n";
	print STDERR "  ex 2) $0 -o myHeader -r ref.fa -c contigs.fa -lib1 lib1.fq -lib2 'lib2_1.fq lib2_2.fq'\n";
	exit($return);
}

if(!defined $outHeader){ print STDERR "\nNeed a header of files!!\n\n"; help(1); }
if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }
if(! -e $refFn){ print STDERR "\n[Error] Could not find the reference file '$refFn'.\n"; exit(1); }
	
if(!defined $contigFn){ print STDERR "\nNeed a contig sequence file!!\n\n"; help(1); }
if(! -e $contigFn){ print STDERR "\n[Error] Could not find the contig sequence file '$contigFn'.\n"; exit(1); }
if($contigFn2 && ! -e $contigFn){ print STDERR "\n[Error] Could not find the contig sequence file '$contigFn2'.\n"; exit(1); }

if(!defined $lib1 && !defined $lib2){ print STDERR "\nNeed fastq sequence files!!\n\n"; help(1); }

if(!defined $lib1){ 
	$lib1 = $lib2;
	$lib2 = undef;
}


my $EXE_DIR= dirname($0);

my $qsubFn = "qsub.$outHeader.sh";
my $out;
open($out, "> $qsubFn") || die "Could not create a script file : $!\n";

print $out "#!/bin/bash\n";
print $out "#PBS -q sandybridge_q\n#PBS -Wgroup_list=sfx\n";
print $out "#PBS -N $outHeader\n";
print $out "#PBS -j oe\n";
print $out "#PBS -l walltime=150:00:00 \n"; 
print $out "#PBS -l nodes=1:ppn=4\n"; ## any 1 node and 5 processes
print $out "#PBS -d $ENV{'PWD'}\n";
print $out "echo start $outHeader `date`\n";

print $out "module load samtools\n";
print $out "module load bwa\n";

$cmd = "$EXE_DIR/reviseq_pipeline.pl -o $outHeader -n $iter_num -r $refFn -c $contigFn" . ($contigFn2 ? " -c2 $contigFn2" : ""). " -lib1 '$lib1'" . ($lib2?" -lib2 '$lib2'":"") . " -step $step";
print $out "echo #### $cmd\n";
print $out "$cmd\n";


print $out "echo end $outHeader `date`\n";
close($out);

system("qsub $qsubFn");
#unlink($qsubFn);
