#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2010
#
# pileup2wig_wBin.pl
# reads pileup, gets coverages within bin size and writes in a WIG format.
# Copyright   : GPL V3 (http://www.gnu.org/licenses/gpl-3.0.html)
######################################

use strict;
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, $chr, $verbose, $binSize, $pileupFn, $targetFastaFn, $headerFn, $outFn);
$cmd = "$0 @ARGV";
$binSize = 10;

GetOptions(
	"h|?|help"		=> \$helpFlag,
	"verbose"		=> \$verbose,
	"bin=i"	=> \$binSize,
	"chr=s"	=> \$chr,
	"sam=s"	=> \$headerFn,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;

if(!defined $headerFn){ print STDERR "\nNeed a sam file for reference lengths!!\n\n"; help(1); }
if($#ARGV == -1){ print STDERR "\nNeed pileup files!!\n\n"; help(1); }
if(defined $outFn){
	foreach my $fn (@ARGV) {
		if($fn eq $outFn){ print STDERR " Error) Input and output files are same \n"; exit; }
	}
}


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 -bin  <bin size>  -s <sam header or sam file>  pileup1   [-o <out file>] \n";
	print STDERR "  ex) $0 -bin $binSize -s s_3.B.suis_1330.sam s_3.B.suis_1330.pileup -o s_3.B.suis_1330.$binSize.wig\n";
	exit($return);
}

my ($in, $out, @arr, $i, $name);


my (%refLen);
$in = openInput($headerFn);
while(<$in>){
	s/[\r\n]+//g;
	next if(/^#/ || /^\s*$/);
	last if(!/^@/);
	if(/^\@SQ\tSN:([\w_+\-:\.\|#]+)\tLN:(\d+)/){
		$refLen{$1} = $2;		
	}	
}
close($in);

my (@allCov, $fNum, %maxI, $j, @refNames);

$fNum = $#ARGV+1;
$i = 0;
#for($i = 0; $i < $fNum; $i++) {
	getCoverageFromPileup($ARGV[$i]);

	$ARGV[$i] = `basename $ARGV[$i]`;	
	$ARGV[$i] =~ s/[\r\n]//g;
	$ARGV[$i] =~ s/\.pileup$//;
#}


if(!defined $outFn){
	$outFn = $ARGV[0] . ".wig";
}

$out = openOutput($outFn);
print " Writing $outFn..\n";


foreach my $ref (@refNames) { 
	print $out "variableStep\tchrom=$ref\n";
	for($i = 0; $i <= $maxI{$ref}; $i++){
		print $out sprintf("%d\t%.2f\n", ($i*$binSize+1), (defined $allCov[0]{$ref}[$i] ? $allCov[0]{$ref}[$i] : 0) );
	}
}



close($out) if(defined $outFn);



sub getCoverageFromPileup{
	my ($pileupFn) = @_;

	my($i, $ref, $readCnt, %cov, $found);
	$in = openInput($pileupFn);
	$found = 0;

	while(<$in>){
		s/[\r\n]+//g;
		next if(/^#/ || /^\s*$/);

		@arr = split /\t/;
#		if($arr[0] =~ /^[0-9]+$/ || $arr[0] eq "X" || $arr[0] eq "Y"){
#			$arr[0] = "chr" . $arr[0];
#		}

		last if(defined $chr && $found == 1 && $chr ne $arr[0]);
		next if(defined $chr && $chr ne $arr[0]);

		$refNames[$#refNames + 1] = $arr[0] if($#refNames == -1 || $refNames[$#refNames] ne $arr[0]);

		$i = int($arr[1]/$binSize);
		$cov{$arr[0]}[$i] = 0 if(!defined $cov{$arr[0]}[$i]);
		$readCnt = (defined $arr[7] ? $arr[7] : $arr[3]);
		if(!defined $cov{$arr[0]}[$i] || !defined $readCnt){
			print "$ARGV[0], i ; $i, arr[0] : $arr[0], readCnt $readCnt, $cov{$arr[0]}[$i] += $readCnt;\n$_\n"; exit;
		}
		$cov{$arr[0]}[$i] += $readCnt;
		$found = 1;
	}
	close($in);
	
	foreach $ref (@refNames) {
		for($i = 0; $i <= $#{$cov{$ref}}; $i++){
			if(!defined $cov{$ref}[$i]) { $cov{$ref}[$i] = 0; }
			else {$cov{$ref}[$i] = sprintf("%.2f", $cov{$ref}[$i]/$binSize); }
		}
	
		$maxI{$ref} = $#{$cov{$ref}} if(!defined $maxI{$ref} || $maxI{$ref} < $#{$cov{$ref}});
	}
	$allCov[$#allCov +1] = \%cov;
}

sub openInput
{
	my ($fn) = @_;

	return *STDIN unless defined $fn;

	my ($fd);
	open($fd, $fn =~ /\.gz/ ? "zcat $fn|" : ($fn =~ /\.bz2/ ? "bunzip2 -c $fn|" : $fn)) || die "Could not open '$fn' : $!\n";
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


sub rc
{
	my ($seq, $type) = @_;

	my $rc = reverse $seq;
	if(defined $type && ($type eq "rna" || $type eq "RNA")) # RNA
	{   $rc =~ y/acgtuACGTU/ugcaaUGCAA/;  }
	else ## DNA
	{   $rc =~ y/acgtuACGTU/tgcaaTGCAA/;  }

	return $rc;
}


