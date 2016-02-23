#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# getLocalAssemblyPoints.pl
# searches local assembly points from a SAM file..
######################################

use strict;
use warnings "all";
use Getopt::Long qw(:config no_ignore_case);
use IO::Handle;
STDOUT->autoflush(1);
STDERR->autoflush(1);


my ($cmd, $helpFlag, $verboseFlag, $outFn, $samFn, $refFn);
my ($mergeGap, @motifRepeatCnt, $numClipped, $maxTargetLen, $maxMotifLen);

$maxTargetLen = 200;
$maxMotifLen = 8;
$mergeGap = 5;
$numClipped = 5;

$cmd = "$0 @ARGV";  ### $0 : command..  store all inputs(@ARGV) as a string in $cmd..
my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
print "### $cmd\n[Start] at $Month[$mon] $mday $hour:$min:$sec\n";


GetOptions(
	"h|?|help"	=> \$helpFlag,
	"verbose"		=> \$verboseFlag,

	"num=s"	=> \$numClipped,
	"sam=s"	=> \$samFn,
	"ref=s"	=> \$refFn,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;

$samFn = shift if(!defined $samFn);
$refFn = shift if(!defined $refFn);
$outFn = shift if(!defined $outFn);

if(!defined $samFn){ print STDERR "\nNeed a sam file!!\n\n"; help(1); }
if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }
if(!defined $outFn){ print STDERR "\nNeed a output file!!\n\n"; help(1); }


if(defined $outFn && $samFn eq $outFn)
{ print STDERR " Error) Input and output files are same \n"; exit(1); }


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  [-s] <sam file>  -r <reference file> -o <out file> [-n number of clipped reads to search (default : $numClipped)]\n";
	print STDERR "  ex) $0 -s test.sam  -r ref.fa -o test.out.txt\n\n";
	exit($return);
}
my ($in, $out, $i);


################################### set the number of repeat for each motif length.
$motifRepeatCnt[1] = 10;
$motifRepeatCnt[2] = 6;
$motifRepeatCnt[3] = 5;
$motifRepeatCnt[4] = 4;
for($i = 5; $i <= $maxMotifLen; $i++){
	$motifRepeatCnt[$i] = 3;
}
####################################


my ($rBuf, $rin);
my $readLen = 0;

my (@leftClipped, @rightClipped, $tandemRepeatList, $refName, $refSeq, %refLen);

$rin = openInput($refFn);
$in = openInput($samFn);
$out = openOutput($outFn);

while(<$in>){
	if(/^@/){
		if(/SN:(.*)\tLN:(\d+)/){
			$refLen{$1} = $2;
		}
	}

	my @arr = split /\t/;
	next if($#arr < 10 || ($arr[1]&0x4) != 0);

	my $len = length($arr[9]);
	$readLen = $len if($readLen < $len);

	if(!$refName){
		readRefSeq($arr[2]);
		$tandemRepeatList = searchRepeat(\$refSeq);
	}
	if($refName && $refName ne $arr[2]){
		findAssemblyPoints();
		
		readRefSeq($arr[2]);
		$tandemRepeatList = searchRepeat(\$refSeq);

		$#leftClipped = -1;
		$#rightClipped = -1;
	}

	if($arr[5] =~ /^(\d+)S/){
		$leftClipped[$arr[3]]++;		
	}
	elsif($arr[5] =~ /(\d+)S$/){
		my $tEnd = $arr[3] - 1;
		while($arr[5] =~ /(\d+)([IMD])/g){
			my ($len, $type) = ($1, $2);	
			$tEnd += $len if($type ne "I");
		}
		$rightClipped[$tEnd]++;
	}	
	$refName = $arr[2];
}
findAssemblyPoints() if($refName);
close($out);

($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "[End $0] at $Month[$mon] $mday $hour:$min:$sec\n";


sub findAssemblyPoints{
	print "testing $refName\n";
	my ($i, $j, $prevJ, $r, $rNum);

	$r = 0;
	$rNum = $#$tandemRepeatList+1;  ### repeat num;
	
	for($i = 1; $i <= $refLen{$refName}; $i++){
		if($leftClipped[$i] && $leftClipped[$i] >= $numClipped){
			while($r < $rNum && $tandemRepeatList->[$r]{end} < $i){
				$r++;
			}
			print "$i, left : $leftClipped[$i]\n";
			$prevJ = 0;
			for($j = $i; $j <= $refLen{$refName} && $j < $i+$readLen-10; $j++){
				if($rightClipped[$j] && $rightClipped[$j] >= $numClipped){
					print "$j, right : $rightClipped[$j]\n";
					$prevJ = $j;
				}
			}
			if($prevJ != 0){
				my ($start, $end) = ($i, $prevJ);
				if($r < $rNum && $tandemRepeatList->[$r]{end} >= $i && $tandemRepeatList->[$r]{start} <= $prevJ){
					$start = $tandemRepeatList->[$r]{start} if($start > $tandemRepeatList->[$r]{start});
					$end = $tandemRepeatList->[$r]{end} if($end < $tandemRepeatList->[$r]{end});
				}  ### overlap to tandemRepeat..
				print $out "$refName\t$start\t$end\n";
				print "$refName\t$start:$i($leftClipped[$i])\t$end:$prevJ($rightClipped[$prevJ])\n";
				$i = $end+1;
			}
			print "i : $i, j : $j\n";
		}
	}
}




sub searchRepeat
{	
	my ($tSeq) = @_; ### $tSeq is a scalar.. not simple variable.
	my ($i, $j, $nr, @arrRepeat, @arr); # store info of all SSRs from each sequence

	$nr = 0;
	for ($i=1; $i <= $maxMotifLen; $i++) # scalar(@typ); $i++) #check each motif class
	{
		my $motiflen = $i; #$typ[$i];
		my $minreps = $motifRepeatCnt[$i]-1; #$typrep{$typ[$i]} - 1;		
		my $search = "(([ATGC]{$motiflen})\\2{$minreps,})";
		while ( $$tSeq =~ /$search/ig ) #scan whole sequence for that class
		{
			my $motif = uc $2;
			my $redundant; #reject false type motifs [e.g. (TT)6 or (ACAC)5]
			for ($j = $motiflen - 1; $j > 0; $j--)
			{
				my $redmotif = "([ACGT]{$j})\\1{".($motiflen/$j-1)."}";				
				$redundant = 1 if ( $motif =~ /$redmotif/ )
			};

			next if $redundant;
			$arrRepeat[$nr]{motif} = $motif;

			my $ssr = $1;			
			$arrRepeat[$nr]{repeat} = length($ssr) / $motiflen;
			$arrRepeat[$nr]{end} = pos($$tSeq);
			$arrRepeat[$nr]{start} = $arrRepeat[$nr]{end} - length($ssr) + 1;
			$nr++;
			
		}
	}
	@arr = sort {$a->{start} <=> $b->{start}} @arrRepeat; #put SSRs in right order

	@arrRepeat = ();

	for($i = 0, $j = -1; $i < scalar(@arr); $i++){		
		if($i == 0 || $arr[$i-1]{end} < $arr[$i]{start} - $mergeGap){
			$j++;
			$arrRepeat[$j]{start} = $arr[$i]{start};
			$arrRepeat[$j]{end} = $arr[$i]{end};
			$arrRepeat[$j]{lastMotif} = $arr[$i]{motif};
			$arrRepeat[$j]{motif} = $arr[$i]{motif};
		}
		else{
			$arrRepeat[$j]{end} = $arr[$i]{end};
			$arrRepeat[$j]{lastMotif} = $arr[$i]{motif};
			$arrRepeat[$j]{motif} .= "|" .$arr[$i]{motif};
		}
		#print "$arrRepeat[$j]{start}, $arrRepeat[$j]{end}, $arrRepeat[$j]{motif}\n";
	}

	for($i = 0; $i < scalar(@arrRepeat); $i++){
		my $motif = $arrRepeat[$i]{lastMotif};
		my $extended = 0;
		my $seq = substr($$tSeq, $arrRepeat[$i]{start}-1, ($arrRepeat[$i]{end} - $arrRepeat[$i]{start} + 1));
		for($j = 0; $j < length($motif); $j++){
			my $c = substr($motif, $j, 1);			
			if($c eq substr($$tSeq, $arrRepeat[$i]{end} + $j, 1)) {			
				$seq .= $c;
				$extended++;
			}
			else{
				last;
			}
		}
		$arrRepeat[$i]{end} += $extended;
	}
	return \@arrRepeat;	
}


sub readRefSeq{
	my $rName = shift;

	my ($name);
	print STDERR "reading '$rName' sequence...\n";
	$rBuf = <$rin> if(!$rBuf); 
	while(defined $rBuf){
		if($rBuf =~ /^>/){			
			$name = $';
			$name =~ s/^\s*//g; $name =~ s/\s.*//g;
			if($name eq $rName){
				$rBuf = <$rin>;
				last;
			}
		}
		$rBuf = <$rin>;
	}
	
	$refSeq = "";
	
	while(defined $rBuf){
		last if($rBuf =~ /^>/);
		$rBuf =~ s/[\r\n]+//g;
		$refSeq .= uc $rBuf;
		$rBuf = <$rin>;
	}

	if($refSeq eq ""){
		print STDERR "[Error] Could not read sequences for '$rName'. The order of references in '$refFn' and '$samFn' may be differ.\n";
		exit(1);
	}
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

