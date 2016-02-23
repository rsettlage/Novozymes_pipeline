#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# run_bwasw_contigs.pl
#   It takes an output SAM file name, a reference file, a consensus sequence file and contig sequence files as inputs, 
#     and runs BWASW to align the consensus and contig sequences to the reference sequence.
#   And it generates the output SAM file, "[header].bam", "[header].pileup", "[header].out.coverage", "[header].out.gap", "[header].out.var" and "[header].out.seq" files. 
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
my ($contigFn, $contigFn2, $consensusFn, $refFn, $samOutFn);


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

	"out=s"	=> \$samOutFn,
	"ref=s"	=> \$refFn,

	"map=s"	=> \$consensusFn,
	"c=s"	=> \$contigFn,
	"c2=s"	=> \$contigFn2,
) || help(1);


help(0) if defined $helpFlag;


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 -o <output sam file>  -r <reference FASTA file> -m <consensus file from mapping>  -c <contig file from assembly> [-c2 <second contig FASTA file>]\n\n";
	print STDERR "  ex) $0 -o output.sam -r ref.fa -m p1.align.seq -c contigs.fa\n";
	exit($return);
}

if(!defined $samOutFn){ print STDERR "\nNeed a name of an output sam file!!\n\n"; help(1); }

if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }
if(! -e $refFn){ print STDERR "\n[Error] Could not find the reference file '$refFn'.\n"; exit(1); }

if(!defined $consensusFn){ print STDERR "\nNeed a consensus sequence file!!\n\n"; help(1); }
if(! -e $consensusFn){ print STDERR "\n[Error] Could not find the consensus file '$consensusFn'.\n"; exit(1); }

if(!defined $contigFn){ print STDERR "\nNeed a contig sequence file!!\n\n"; help(1); }
if(! -e $contigFn){ print STDERR "\n[Error] Could not find the contig file '$contigFn'.\n"; exit(1); }
if($contigFn2 && ! -e $contigFn){ print STDERR "\n[Error] Could not find the contig sequence file '$contigFn2'.\n"; exit(1); }



######################## print starting time
$cmd = "$0 @ARGV";  ### $0 : command..  store all inputs(@ARGV) as a string in $cmd..
my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "### $cmd\n[Start] at $Month[$mon] $mday $hour:$min:$sec\n";
######################## print starting time



my $cSam=basename($consensusFn);
$cSam =~ s/.(fa|fasta|seq|fna|txt)$/.sam/;
$cSam = dirname($samOutFn) . "/$cSam";

if (! -e "$refFn.bwt"){	
	runProcess("$BWA index $refFn");
}

runProcess("$BWA bwasw $refFn $consensusFn -f $cSam");
my $outHeader=$cSam;
$outHeader =~ s/.sam$//;
my $bamFn = "$outHeader.bam";
runProcess("$SAMTOOLS view -hS $cSam -b | $SAMTOOLS sort - $outHeader");
runProcess("$SAMTOOLS view -h $bamFn -o $cSam");
my $log= "$outHeader.log";
my $fixed="$outHeader.fixed.sam";
runProcess("$EXE_DIR/fix_bwasw_sam.pl -r $refFn $cSam -o $fixed -v", $log);

system("mv -f $fixed $cSam");
unlink($bamFn);

runProcess("$BWA bwasw -w 100 $refFn $contigFn -f $samOutFn");
$outHeader=$samOutFn;
$outHeader =~ s/.sam$//;
$bamFn = "$outHeader.bam";
runProcess("$SAMTOOLS view -hS $samOutFn -b | $SAMTOOLS sort - $outHeader");
runProcess("$SAMTOOLS view -h $bamFn -o $samOutFn");
$log= "$outHeader.log";
$fixed="$outHeader.fixed.sam";
runProcess("$EXE_DIR/fix_bwasw_sam.pl -r $refFn $samOutFn -o $fixed -v", $log);
unlink($bamFn);

if($contigFn2){
	my $c2SamHeader = $samOutFn;
	$c2SamHeader =~ s/.sam$//;
	$c2SamHeader .= ".c2";
	my $c2Sam = "$c2SamHeader.sam";
	$bamFn = "$c2SamHeader.bam";
	
	runProcess("$BWA bwasw -w 100 $refFn $contigFn2 -f $c2Sam");	
	runProcess("$SAMTOOLS view -hS $c2Sam -b | $SAMTOOLS sort - $c2SamHeader");
	runProcess("$SAMTOOLS view -h $bamFn -o $c2Sam");
	$log= "$c2SamHeader.log";
	my $c2fixed="$c2SamHeader.fixed.sam";
	runProcess("$EXE_DIR/fix_bwasw_sam.pl -r $refFn $c2Sam -o $c2fixed -v", $log);
	
	system("echo \"cat $c2fixed | perl -ne 'print if(not /^@/)' >> $fixed\"");
	system("cat $c2fixed | perl -ne 'print if(not /^@/)' >> $fixed");
	system("echo \"cat $cSam | perl -ne 'print if(not /^@/)' >> $fixed\"");
	system("cat $cSam | perl -ne 'print if(not /^@/)' >> $fixed");
	unlink($c2fixed);
	unlink($c2Sam);
	unlink($bamFn);
}

system("echo \"cat $cSam | perl -ne 'print if(not /^@/)' >> $fixed\"");
system("cat $cSam | perl -ne 'print if(not /^@/)' >> $fixed");
system("echo \"cat $cSam | perl -ne 'print if(not /^@/)' >> $fixed\"");
system("cat $cSam | perl -ne 'print if(not /^@/)' >> $fixed");
unlink($cSam);

#exit;
############################################################# sorting
runProcess("$SAMTOOLS view -hS $fixed -b | $SAMTOOLS sort - $outHeader");
$bamFn = "$outHeader.bam";

#exit;
runProcess("$SAMTOOLS view -h $bamFn -o $samOutFn");
unlink($fixed);

my $pileupFn = $outHeader.".pileup";
runProcess("$SAMTOOLS pileup -f $refFn $bamFn > $pileupFn");

$outHeader =~ /(\w+)$/;
my $blockHeader = $1;

runProcess("$EXE_DIR/getConsensusFromPileup.pl $pileupFn $samOutFn -header $blockHeader -indel 0.3 -nogap");

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

