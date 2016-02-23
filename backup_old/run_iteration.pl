#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w


#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# run_iteration.pl
#  It takes an output header, a reference file, FASTQ files containing read sequences, 
#   an iteration starting index i and a number of iteration n as inputs, and runs BWA 'n' times. 
#  Each iteration step has an index, "i + [iteration count]" 
#   and generates "[header].it_[index].sam", "[header].it_[index].bam", "[header].it_[index].pileup", "[header].it_[index].out.coverage", "[header].it_[index].out.gap", "[header].it_[index].out.var" and "[header].it_[index].out.seq" files.
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

my ($cmd, $helpFlag, $iterStart, $iterNum, $outHeader);
my ($lib1, $lib2, $refFn);


my $EXE_DIR= dirname($0);

GetOptions(
	"h|?|help"	=> \$helpFlag,
	
	"ref=s"	=> \$refFn,
	"out=s"	=> \$outHeader,

	"iter=s"	=> \$iterStart,
	"num=s"	=> \$iterNum,

	"lib1=s"	=> \$lib1,
	"lib2=s"	=> \$lib2,
) || help(1);


help(0) if defined $helpFlag;


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 -o <header of output> -r <reference FASTA file> -i <iteration starting index>  -n <number of iteration> -lib1 <fastq files> [-lib2 <fastq files>] \n\n";
	print STDERR "  ex 1) $0 -o myHeader -r ref.fa -i 1 -n 5 -lib1 'lib1_1.fq lib1_2.fq'\n";
	print STDERR "  ex 2) $0 -o myHeader -r ref.fa -i 1 -n 5 -lib1 lib1.fq -lib2 'lib2_1.fq lib2_2.fq'\n";
	exit($return);
}

if(!defined $outHeader){ print STDERR "\nNeed an output header!!\n\n"; help(1); }
if(!defined $iterStart){ print STDERR "\nNeed an iteration starting index!!\n\n"; help(1); }
if(!defined $iterNum){ print STDERR "\nNeed the number of iteration!!\n\n"; help(1); }

if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }
if(! -e $refFn){ print STDERR "\n[Error] Could not find the reference file '$refFn'.\n"; exit(1); }
	
if(!defined $lib1 && !defined $lib2){ print STDERR "\nNeed fastq sequence files!!\n\n"; help(1); }

if(!defined $lib1){ 
	$lib1 = $lib2;
	$lib2 = undef;
}


my $i;
my $end = $iterStart + $iterNum - 1;
my $bNoChange = 0;
for($i = $iterStart; $i <= $end; $i++){
	runProcess("$EXE_DIR/run_bwa.pl -o $outHeader.it_$i -r $refFn -nogap -lib1 '$lib1'" . ($lib2?" -lib2 '$lib2'":""));
	if(isSameSeq($refFn, "$outHeader.it_$i.align.seq") ) { 
		if($iterStart < $i && $i < $end){			
			runProcess("$EXE_DIR/localAssembly.pl -o $outHeader.it_$i.align.seq -r $refFn -s $outHeader.it_$i.sam");
			if(isSameSeq($refFn, "$outHeader.it_$i.align.seq") ) { 
				$bNoChange = 1; last;
			}
		}
		else{
			$bNoChange = 1; last;
		}
	}
	last if($i == $end);
	
	#$refFn = "$outHeader.it_$i.align.seq";
	$refFn = "$outHeader.it_".($i+1).".fa";
	system("cat $outHeader.it_$i.align.seq | sed 's/block/ref/' > $refFn");	
}
my $out = openOutput("$outHeader.it.flag");
print $out "$i\n$bNoChange\n";
close($out);

runProcess("$EXE_DIR/getProblematicReads.sh $refFn $outHeader.it_$i.sam");

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


sub runProcess {
	my ($cmd, $log) = @_;
	my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
	my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
	print "---[Process start at $Month[$mon] $mday $hour:$min:$sec] -----------------\n###$cmd" . (defined $log ? " > $log" : ""). "\n";
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
