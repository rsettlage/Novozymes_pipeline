#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# run_bwa.pl
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


my ($cmd, $helpFlag, $noGapFlag, $noPileupFlag, $minRate_indelReplace);
my ($lib1, $lib2, $refFn, $outHeader);


my $EXE_DIR= dirname($0);

my $SAMTOOLS="$EXE_DIR/samtools_0.1.16";

#my $BWA="bwa" ## use this line if you use 'bwa' installed at the machine.
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
	"nogap"		=> \$noGapFlag,     #### for getConsensusFromPileup.pl
	"nopileup"		=> \$noPileupFlag,

	"indel"		=> \$minRate_indelReplace,   #### for getConsensusFromPileup.pl

	"out=s"	=> \$outHeader,
	"ref=s"	=> \$refFn,

	"lib1=s"	=> \$lib1,
	"lib2=s"	=> \$lib2,
) || help(1);


help(0) if defined $helpFlag;


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 -o <header of output>  -r <reference FASTA file> -lib1 <fastq files> [-lib2 <fastq files>] \n\n";
	print STDERR "  ex 1) $0 -o myHeader -r ref.fa -lib1 'lib1_1.fq lib1_2.fq'\n";
	print STDERR "  ex 2) $0 -o myHeader -r ref.fa -lib1 lib1.fq -lib2 'lib2_1.fq lib2_2.fq'\n";
	exit($return);
}

if(!defined $outHeader){ print STDERR "\nNeed an output header!!\n\n"; help(1); }
if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }
if(! -e $refFn){ print STDERR "\n[Error] Could not find the reference file '$refFn'.\n"; exit(1); }
	
if(!defined $lib1 && !defined $lib2){ print STDERR "\nNeed fastq sequence files!!\n\n"; help(1); }

if(!defined $lib1){ 
	$lib1 = $lib2;
	$lib2 = undef;
}

my $dirName = dirname($outHeader);
if($dirName ne "."){
	mkdir($dirName);
}



######################## print starting time
$cmd = "$0 @ARGV";  ### $0 : command..  store all inputs(@ARGV) as a string in $cmd..
my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "### $cmd\n[Start] at $Month[$mon] $mday $hour:$min:$sec\n";
######################## print starting time



my (@lib1_fn, @lib2_fn);

$lib1 =~ s/^[\s,]+//; $lib1 =~ s/[\s,]+$//; #### removing spaces at the left and right ends of the string
@lib1_fn = split /[\s,]+/, $lib1;
if($#lib1_fn == -1){ print STDERR "\nNeed fastq sequence files!!\n\n"; help(1); }
foreach my $fn (@lib1_fn) {
	if(! -e $fn){ print STDERR "\n[Error] Could not find the '$fn' file.\n"; exit(1); }
}
$#lib1_fn = 1 if($#lib1_fn > 1);  ## only take the first two files

if($lib2){
	$lib2 =~ s/^[\s,]+//; $lib2 =~ s/[\s,]+$//;
	@lib2_fn = split /[\s,]+/, $lib1;
	foreach my $fn (@lib2_fn) {
		if(! -e $fn){ print STDERR "\n[Error] Could not find the '$fn' file.\n"; exit(1); }
	}
	$#lib2_fn = 1 if($#lib2_fn > 1);
}


my ($i, $t, @thread, @sa1_Fn, @sa2_Fn);

############################################################## mapping
runProcess("$BWA index $refFn");

for($i = 0, $t = 0; $i <= $#lib1_fn; $i++, $t++){
	$sa1_Fn[$i] = $outHeader.".lib1.$i.sai";
	$thread[$t] = new Thread(\&run_aln, $refFn, $lib1_fn[$i], $sa1_Fn[$i]);
}
for($i = 0; $i <= $#lib2_fn; $i++, $t++){
	$sa2_Fn[$i] = $outHeader.".lib2.$i.sai";	
	$thread[$t] = new Thread(\&run_aln, $refFn, $lib2_fn[$i], $sa2_Fn[$i]);
}
my $tNum = $t;

for($t = 0; $t < $tNum; $t++) {
	$thread[$t]->join;
}
############################################################## mapping


############################################################# generating sam file
my $samFn = $outHeader . ".sam";
$thread[0] = new Thread(\&run_sam, $refFn, \@lib1_fn, \@sa1_Fn, $samFn, $#lib1_fn);
$tNum = 1;
if($lib2){	
	$thread[1] = new Thread(\&run_sam, $refFn, \@lib2_fn, \@sa2_Fn, "$outHeader.lib2.sam", $#lib2_fn);
	$tNum = 2;
}

for($t = 0; $t < $tNum; $t++) {
	$thread[$t]->join;
}
############################################################# generating sam file



############################################################# merging sam files
if($lib2){ 	system("grep '^[A-Z]' $outHeader.lib2.sam >> $samFn"); }


############################################################# sorting
runProcess("$SAMTOOLS view -hS $samFn -b | $SAMTOOLS sort - $outHeader");
my $bamFn = $outHeader.".bam";
runProcess("$SAMTOOLS view $bamFn -h -o $samFn");

$outHeader =~ /(\w+)$/;
my $blockHeader = $1;

if(!$noPileupFlag){
	my $pileupFn = $outHeader.".pileup";
	unlink("$refFn.fai");
	runProcess("$SAMTOOLS pileup -f $refFn $bamFn > $pileupFn");

	runProcess("$EXE_DIR/getConsensusFromPileup.pl $pileupFn $samFn -pos -header $blockHeader" .
		($noGapFlag?" -nogap -ref $refFn":"") . ($minRate_indelReplace?" -indel $minRate_indelReplace":""));
}

system("rm $outHeader.lib* -f");

($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "[End $0] at $Month[$mon] $mday $hour:$min:$sec\n";




sub run_aln
{
	my ($refFn, $fq, $sai) = @_;
	runProcess("$BWA aln $refFn $fq -f $sai");
}


sub run_sam
{
	my ($refFn, $fqFiles, $saiFiles, $samFn, $pe) = @_;
	if($pe){ ### paired end mapping
		runProcess("$BWA sampe $refFn $saiFiles->[0] $saiFiles->[1] $fqFiles->[0] $fqFiles->[1] -f $samFn");
	}
	else{ ### single end mapping
		runProcess("$BWA samse $refFn $saiFiles->[0] $fqFiles->[0] -f $samFn");
	}
}


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

	#return *STDIN unless defined $fn;

	my ($fd);
	if($fn =~ /\.bam/) { open($fd, "$SAMTOOLS view $fn -h |");  }
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
