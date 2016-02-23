#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# clippedNum2Wig.pl
######################################

use strict;
use warnings "all";
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, $verboseFlag, $outFn, $binSize);
my ($samFn);
$cmd = "$0 @ARGV";  ### $0 : command..  store all inputs(@ARGV) as a string in $cmd..
my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
print "### $cmd\n[Start] at $Month[$mon] $mday $hour:$min:$sec\n";

$binSize = 1;
GetOptions(
	"h|?|help"	=> \$helpFlag,
	"verbose"		=> \$verboseFlag,

	"bin=i" => \$binSize,
	"sam=s"	=> \$samFn,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;


$samFn = shift if(!defined $samFn);
$outFn = shift if(!defined $outFn);


if(!defined $samFn){ print STDERR "\nNeed a sam file!!\n\n"; help(1); }
if(!-e $samFn){ print STDERR "\nCould not find the sam file $samFn!!\n\n"; help(1); }
if(!defined $outFn){ print STDERR "\nNeed a output file!!\n\n"; help(1); }


if(defined $outFn && ($samFn eq $outFn))
{ print STDERR " Error) Input and output files are same \n"; exit(1); }


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  [-s] <sam file>  [-o] <output wig file>\n";
	print STDERR "  ex) $0 -s test.sam -o test.wig\n\n";
	exit($return);
}
my ($in, $out, @arr, $i, $refName, @cov);
$refName = "";

$in = openInput($samFn);
$out = openOutput($outFn);
while(<$in>){
	s/[\r\n]+//g;
	next if(/^@/ || /^\s*$/);

	## 0: ID, 1: flag, 2: reference, 3: position, 4: mapping score, 5: cigar, 6/7/8: mate info, 9: seq, 10: quality, 11..: tag.
	my @arr = split /\t/;
	next if($arr[5] !~ /S/);
	if($refName ne "" && $refName ne $arr[2]){
		print $out "variableStep\tchrom=$refName\n";
		for($i = 0; $i <= $#cov; $i++){
			next if(!$cov[$i] || $cov[$i] <3);
			print $out sprintf("%d\t%.2f\n", ($i*$binSize+1), $cov[$i]/$binSize);
		}
		@cov = ();
	}
	$refName = $arr[2];

	my $pos = $arr[3]-1;
	if($arr[5] =~ /^(\d+)S/ && $1 > 3){
		$cov[int($pos/$binSize)]++;
	}
	if($arr[5] =~ /(\d+)S$/ &&  $1 > 3){
		while($arr[5] =~ /(\d+)[MD]/g){
			$pos += $1;
		}
		$cov[int($pos/$binSize)]++;
	}	
}

close($in);

if($refName ne ""){
	print $out "variableStep\tchrom=$refName\n";
	for($i = 0; $i <= $#cov; $i++){
		next if(!$cov[$i] || $cov[$i] <3);
		print $out sprintf("%d\t%.2f\n", ($i*$binSize+1), $cov[$i]/$binSize);
	}
}

($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "[End $0] at $Month[$mon] $mday $hour:$min:$sec\n";




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


sub rc #### return reverse complement
{
	my ($seq, $type) = @_;

	my $rc = reverse $seq;
	if(defined $type && ($type eq "rna" || $type eq "RNA")) # RNA
	{   $rc =~ y/acgtuACGTU/ugcaaUGCAA/;  }
	else ## DNA
	{   $rc =~ y/acgtuACGTU/tgcaaTGCAA/;  }

	return $rc;
}

