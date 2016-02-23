#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# reviseq_pipeline.pl
#   It takes an output header, a reference file and FASTQ files containing read sequences as inputs, and runs BWA. 
#   And it generates "[header].sam", "[header].bam", "[header].pileup", "[header].align.coverage", "[header].align.gap", "[header].align.var" and "[header].align.seq" files.
#
# Copyright   : GPL V3 (http://www.gnu.org/licenses/gpl-3.0.html)
######################################

use strict;
use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use IO::Handle;
STDOUT->autoflush(1);
STDERR->autoflush(1);


my ($cmd, $helpFlag, $noGapFlag, $minRate_indelReplace);
my ($lib1, $lib2, $refFn, $contigFn, $contigFn2, $outHeader);

$cmd = "$0 @ARGV";  ### $0 : command..  store all inputs(@ARGV) as a string in $cmd..

my $step = 1;
my $iter_num = 5;

GetOptions(
	"h|?|help"	=> \$helpFlag,

	"out=s"	=> \$outHeader,
	"ref=s"	=> \$refFn,
	"c=s"	=> \$contigFn,
	"c2=s"	=> \$contigFn2,

	"lib1=s"	=> \$lib1,
	"lib2=s"	=> \$lib2,
	"num=s"	=> \$iter_num,

	"step=s"	=> \$step,
) || help(1);


help(0) if defined $helpFlag;


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 -o <header of output> -r <reference FASTA file> -c <contig FASTA file> [-c2 <second contig FASTA file>] -lib1 <fastq files> [-lib2 <fastq files>] [-n <number of iteration>]\n\n";
	print STDERR "  ex 1) $0 -o myHeader -r ref.fa -c contigs.fa -lib1 'lib1_1.fq lib1_2.fq'\n";
	print STDERR "  ex 2) $0 -o myHeader -r ref.fa -c contigs.fa -lib1 lib1.fq -lib2 'lib2_1.fq lib2_2.fq'\n\n";
	exit($return);
}

if(!defined $outHeader){ print STDERR "\nNeed an output header!!\n\n"; help(1); }
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

######################## print starting time
my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "### $cmd\n[Start] at $Month[$mon] $mday $hour:$min:$sec\n";
######################## print starting time

my $EXE_DIR= dirname($0);

runProcess("$EXE_DIR/run_bwa.pl -o $outHeader.p1 -r $refFn -lib1 '$lib1'" . ($lib2?" -lib2 '$lib2'":"")) if($step == 1);

if(isSameSeq($refFn, "$outHeader.p1.align.seq")){
	print "There is no difference between the reference and consensus sequences.\n\n" ;
	exit(0);
}

runProcess("$EXE_DIR/run_bwasw_contigs.pl -o $outHeader.p2.contigs.sam -r $refFn -m $outHeader.p1.align.seq -c $contigFn" . ($contigFn2 ? " -c2 $contigFn2" : "")) if($step <= 2);
runProcess("$EXE_DIR/run_bwa.pl -o $outHeader.p3 -nopileup -r $outHeader.p2.contigs.align.seq -lib1 '$lib1'" . ($lib2?" -lib2 '$lib2'":"")) if($step <= 3);
runProcess("$EXE_DIR/localAssembly.pl -o $outHeader.p3.new.fa -r $outHeader.p2.contigs.align.seq -s $outHeader.p3.sam") if($step <= 4);

#system("rm -f $outHeader.p[123]*.sam $outHeader.p[123]*.pileup");

if(isSameSeq($refFn, "$outHeader.p3.new.fa")){
	print "There is no difference between the reference and consensus sequences.\n\n" ;
	exit(0);
}

runProcess("$EXE_DIR/run_bwa.pl -o $outHeader.p4 -r $outHeader.p3.new.fa -nogap -lib1 '$lib1'" . ($lib2?" -lib2 '$lib2'":"") ." -indel 0.2") if($step <= 5);
#system("rm -f $outHeader.p4.sam $outHeader.p4.pileup");

runProcess("$EXE_DIR/run_iteration.pl -o $outHeader -r $outHeader.p4.align.seq  -i 1 -n $iter_num -lib1 '$lib1'" . ($lib2?" -lib2 '$lib2'":"")) if($step <= 6);


my $in = openInput("$outHeader.it.flag"); ### this file is generated by 'p3.iter.pl' 
my $i = <$in>;   ### the last index of the iteration in 'p3.iter.pl' 
$i =~ s/\n//g;
close($in);

my $seqFn = "$outHeader.final.genome.fa";
system("mv -f $outHeader.it_$i.align.seq $seqFn");

runProcess("$EXE_DIR/twoSeq_comparison.pl -o $outHeader.final -r $refFn -i $seqFn");

runProcess("$EXE_DIR/run_bwa.pl -o $outHeader.final.readAlign -r $seqFn -nopileup  -lib1 '$lib1'" . ($lib2?" -lib2 '$lib2'":""));
runProcess("$EXE_DIR/readCountsAtLoci.pl -ori $refFn -g $outHeader.final.sam -s $outHeader.final.readAlign.sam  -o $outHeader.final.readAlign.var");

#system("rm -f $outHeader.it_$i.sam $outHeader.it_$i.pileup");
#system("rm -f $outHeader.p[2-4]*.[bs]am $outHeader.p[1-4]*.seq.* $outHeader.p[1-4]*.fa*");
for(my $j=0;$j < $i;$j++){
	#system("rm -f $outHeader.it_$j.*");
}

($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "[End $0] at $Month[$mon] $mday $hour:$min:$sec\n";

sub runProcess {
	my ($cmd, $log) = @_;
	my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
	my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
	print "### $cmd\n[Start] at $Month[$mon] $mday $hour:$min:$sec\n";
	print "---[Process start at $Month[$mon] $mday $hour:$min:$sec] -----------------\n$cmd" . (defined $log ? " > $log" : ""). "\n";
	if(defined $log) {$cmd .= " > $log";}
	elsif($cmd !~ />/) {$cmd .= " 2>&1";}
	
	my $res = system($cmd);

	($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
	print "---[Process end at $Month[$mon] $mday $hour:$min:$sec] -------------------\n\n";
	if($res != 0){
		print STDERR "[Error] Stop processing. " .(defined $log ? "See the $log file for detail." : ""). "\n";
		exit(1);
	}
}



sub isSameSeq {
	my ($f1, $f2) = @_;
	my $in1 = openInput($f1);
	my $in2 = openInput($f2);

	$/ = ">";

	my ($str1, $str2, $seq1, $seq2, $name);
	$str1 = <$in1>;
	$str2 = <$in2>;

	my $bSame = 1;

	while($str1 && $str2){
		while($str1){
			if(!(($name, $seq1) = $str1 =~ /(.*?)\n(.*)/s)){
				$str1 = <$in1>;
				next
			}
			last;
		}
		while($str2){
			if(!(($name, $seq2) = $str2 =~ /(.*?)\n(.*)/s)){
				$str2 = <$in2>;
				next
			}
			last;
		}
		
		if($seq1 ne $seq2){
			$bSame = 0;
		}

		$str1 = <$in1>;
		$str2 = <$in2>;
	}

	close($in1);
	close($in2);

	$/ = "\n";
	
	return $bSame;
}




sub openInput
{
	my ($fn) = @_;

	my ($fd);
	if($fn =~ /\.bam/) { open($fd, "samtools view $fn -h |");  }
	else { open($fd, $fn =~ /\.gz/ ? "zcat $fn|" : ($fn =~ /\.bz2/ ? "bunzip2 -c $fn|" : $fn)) || die "Could not open '$fn' : $!\n"; }
	return $fd;
}

sub openOutput
{
	my ($fn) = @_;

	return *STDOUT unless defined $fn;

	my ($fd);
	open($fd, $fn =~ /.gz$/ ? "| gzip -c > $fn" : ($fn =~ /\.bz2/ ? "| bzip2 -c > $fn" : ">$fn")) || die "Could not write '$fn' : $!\n";
	return $fd;
}