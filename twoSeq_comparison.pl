#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# twoSeq_comparison.pl
#   It takes an output header, a reference file and a revised genome sequence as inputs, 
#     and runs BWASW to align the revised genome sequence to the reference sequence. 
#   And it generates the "[header].sam", "[header].bam", "[header].pileup", "[header].out.coverage", "[header].out.gap", "[header].out.var" and "[header].out.seq" files.
#
# Copyright   : GPL V3 (http://www.gnu.org/licenses/gpl-3.0.html)
######################################

use strict;
use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use IO::Handle;
STDOUT->autoflush(1);
STDERR->autoflush(1);

my ($cmd, $helpFlag);
my ($seqFn, $refFn, $outHeader);


my $EXE_DIR= dirname($0);

#my $SAMTOOLS="samtools"; ## use this line if you use 'samtools' from in 'module'
my $SAMTOOLS="$EXE_DIR/samtools_0.1.16";

#my $BWA="bwa";
my $BWA="$EXE_DIR/bwa_0.6.1";

`$BWA 2>&1`;
if($? == -1){
	print STDERR "\nCould not execute the '$BWA' program.\n\n";  exit(1);
}

`$SAMTOOLS 2>&1`;
if($? == -1){
	print STDERR "\nCould not execute the '$SAMTOOLS' program.\n\n";  exit(1);
}

GetOptions(
	"h|?|help"	=> \$helpFlag,

	"out=s"	=> \$outHeader,
	"ref=s"	=> \$refFn,

	"input=s"	=> \$seqFn,
) || help(1);


help(0) if defined $helpFlag;


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 -r <reference FASTA file>  -i <input sequence file>  -o <header of output>\n\n";
	print STDERR "  ex) $0 -r ref.fa -i new.out.seq -o new.final\n\n";
	exit($return);
}

if(!defined $outHeader){ print STDERR "\nNeed an output header!!\n\n"; help(1); }

if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }
if(! -e $refFn){ print STDERR "\n[Error] Could not find the reference file '$refFn'.\n"; exit(1); }

if(!defined $seqFn){ print STDERR "\nNeed a consensus sequence file!!\n\n"; help(1); }
if(! -e $seqFn){ print STDERR "\n[Error] Could not find the consensus file '$seqFn'.\n"; exit(1); }


######################## print starting time
$cmd = "$0 @ARGV";  ### $0 : command..  store all inputs(@ARGV) as a string in $cmd..
my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "### $cmd\n[Start] at $Month[$mon] $mday $hour:$min:$sec\n";
######################## print starting time


my $samFn = "$outHeader.sam";

if (! -e "$refFn.bwt"){	
	runProcess("$BWA index $refFn");
}

runProcess("$BWA bwasw -w 200 $refFn $seqFn -f $samFn");
my $bamFn = "$outHeader.bam";
runProcess("$SAMTOOLS view -hS $samFn -b | $SAMTOOLS sort - $outHeader");
runProcess("$SAMTOOLS view -h $bamFn -o $samFn");
my $log= "$outHeader.log";
my $fixed="$outHeader.fixed.sam";
runProcess("$EXE_DIR/fix_bwasw_sam.pl -r $refFn $samFn -o $fixed -v", $log);

runProcess("$SAMTOOLS view -hS $fixed -b | $SAMTOOLS sort - $outHeader");
runProcess("$SAMTOOLS view -h $bamFn -o $samFn");

my $pileupFn = $outHeader.".pileup";
runProcess("$SAMTOOLS pileup -f $refFn -S $samFn > $pileupFn");
runProcess("$EXE_DIR/getConsensusFromPileup.pl $pileupFn $samFn");

unlink($fixed);
unlink($bamFn);

($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "[End $0] at $Month[$mon] $mday $hour:$min:$sec\n";



sub runProcess {
	my ($cmd, $log) = @_;
	($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
	print "---[Process start at $Month[$mon] $mday $hour:$min:$sec] -----------------\n###$cmd" . (defined $log ? " > $log" : ""). "\n";
	if(defined $log) {$cmd .= " > $log";}
	elsif($cmd !~ />/) {$cmd .= " 2>&1";}
	
	my $res = system($cmd);
	if($res != 0){
		print STDERR "[Error] Stop processing. " .(defined $log ? "See the $log file for detail." : ""). "\n";
		exit(1);
	}
	($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
	print "---[Process end at $Month[$mon] $mday $hour:$min:$sec] -------------------\n\n";
}


sub openInput
{
	my ($fn) = @_;

	return *STDIN unless defined $fn;

	my ($fd);
	open($fd, $fn =~ /\.gz/ ? "zcat $fn|" : 
		($fn =~ /\.bz2/ ? "bunzip2 -c $fn|" : 
		($fn =~ /\.bam/ ? "samtools view -h $fn|" : $fn))) || die "Could not open '$fn' : $!\n";
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

