#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# localAssembly.pl
# searchs local assembly points and assemble sequences at each point.
#
# Copyright   : GPL V3 (http://www.gnu.org/licenses/gpl-3.0.html)
######################################

use strict;
use Getopt::Long qw(:config no_ignore_case);
use IO::Handle;
STDOUT->autoflush(1);
STDERR->autoflush(1);

use Class::Struct;

struct AlignRes => {
	seqA => '$',
	seqAMatchStart => '$',
	seqAMatchEnd => '$',
	seqB => '$',
	seqBMatchStart => '$',
	seqBMatchEnd => '$',
	gap => '$',
	mismatch => '$',
	variations => '@',
};


my ($helpFlag, $verboseFlag, $samFn, $queryFn, $refFn, 
	$outFn, $cmd);
$cmd = "$0 @ARGV";

my ($mergeGap, @motifRepeatCnt, $numClipped, $maxTargetLen, $maxMotifLen);

$maxTargetLen = 200;
$maxMotifLen = 8;
$mergeGap = 5;
$numClipped = 5;

my $totNumNewDB = 0;

GetOptions(
	"h|?|help"		=> \$helpFlag,
	"verbose"		=> \$verboseFlag,

	"num=s"	=> \$numClipped,
	"s|sam=s"	=> \$samFn,
	"r|reference=s"	=> \$refFn,
	"output=s"	=> \$outFn
) || help(1);

help(0) if defined $helpFlag;


if(!defined $samFn){ print STDERR "\nNeed sam files from BWA!!\n\n"; help(1); }
if(!defined $refFn){ print STDERR "\nNeed a reference file!!\n\n"; help(1); }

if(! -e $samFn){ print STDERR "Could not find the '$samFn' file!!\n\n"; help(1); }
if(! -e $refFn){ print STDERR "Could not find the '$refFn' file!!\n\n"; help(1); }

if(defined $outFn && ($refFn eq $outFn || $samFn eq $outFn))
{ print STDERR " Error) Input and output files are same.  Please choose a different output name.\n"; exit(1); }

sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  -r <reference file>  -s <input SAM file> -o <out file>  [-n number of clipped reads to trigger assembly (default : $numClipped)]\n\n";
	print STDERR "  ex) $0 -o ref.revised.fa  -r ref.fa -s input.sam \n\n";
	exit($return);
}


if(!defined $outFn){
	$outFn = $refFn;
	$outFn =~ s/\.fa$./.new.fa/;
	$outFn =~ s/\.seq$/.new.seq/;
	$outFn =~ s/\.fasta$/.new.fasta/;
	$outFn =~ s/\.fna$/.new.fna/;
	$outFn .= ".new.fa" if($outFn eq $refFn);
}

my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
print "### $cmd\n[Start] at $Month[$mon] $mday $hour:$min:$sec\n";

my ($in, $out, $i, $j, $inR, $ref, $refName, $revRefSeq, $refLen, $refSeq, $rBuf, %allPointsInRef);


########################################## find assembly point....
print "###### searching assembly points......\n";
################################### set the number of repeat for each motif length.
$motifRepeatCnt[1] = 10;
$motifRepeatCnt[2] = 6;
$motifRepeatCnt[3] = 5;
$motifRepeatCnt[4] = 4;
for($i = 5; $i <= $maxMotifLen; $i++){
	$motifRepeatCnt[$i] = 3;
}
####################################

my (@leftClipped, @rightClipped, $tandemRepeatList, %hashRefLen);
my $readLen = 0;

$inR = openInput($refFn);
$in = openInput($samFn);

while(<$in>){
	if(/^@/){
		if(/SN:(.*)\tLN:(\d+)/){
			$hashRefLen{$1} = $2;
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
	

	if($arr[5] =~ /^(\d+)[SI]/){
		$leftClipped[$arr[3]]++;		
	}
	elsif($arr[5] =~ /(\d+)[SI]$/){
		my $tEnd = $arr[3] - 1;
		while($arr[5] =~ /(\d+)([MD])/g){
			my ($len, $type) = ($1, $2);	
			$tEnd += $len;
		}
		$rightClipped[$tEnd]++;
	}	
	$refName = $arr[2];
}
findAssemblyPoints() if($refName);
close($refFn);
close($samFn);
########################################## find assembly point....





print "\n###### assembling targets......\n";

my ($pI, $prevTStart, $pointNum, $qualChar,
	$arrPoints, $refAddStart, $newRefSeq, $printedRefID);

my ($ri, @arrRefIDs, %refIDs);
$ri = 0;

my ($totReads, $totDiffReads) = (0,0);
$inR = openInput($refFn);


my ($F, $R) = (0,1);
$pI = 0;
$printedRefID = -1;
$prevTStart = 0;
$refName = "";
$qualChar = "H";
$refAddStart = 1;
$out = openOutput($outFn);

open($in, $samFn =~ /.bam$/ ? "samtools view -h $samFn |"  : $samFn);	

$refName = "";
$prevTStart = 0;

my $line = 0;
my $prevUnmapped = 0;
while(<$in>){
	if(/^@/ && /SN:([\w\-\_:\|\.\/\+]+)/ && !defined $refIDs{$1}) {
		$refIDs{$1} = $ri;
		$arrRefIDs[$ri] = $1;
		$ri++;
		next;
	}
	next if(/^@/);

	$line++;

	my @arr = split /\t/;
	if($arr[2] eq "*"){ 
		if($line < 100){ $prevUnmapped = 1; next; }
		elsif($arr[2] eq "*") {last;}  ### if printed all mapped reads, stop the loop
	}

	if($prevUnmapped){
		exitMsg(__LINE__, "Reads in the '$samFn' file is out of order.\n");
	}
	processRead(\@arr) if(($arr[1]&0x04) == 0);
}

printNewRefSeq();
for($i = $printedRefID+1; $i <= $#arrRefIDs;$i++){
	readRefSeq($arrRefIDs[$i]);
	printFormattedSeq($arrRefIDs[$i], \$refSeq, 80);
}

close($in);
close($inR);
close($out);


($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "[End] at $Month[$mon] $mday $hour:$min:$sec\n";



sub processRead{
	my ($arr) = @_;


	my ($qLen, $seq, $rvSeq, $xa, @ar,
		$strand, $tStart, $cigar, $clippedHead, $clippedTail,
		$gapNum, $isMic, $qStart, $qEnd,
		$len, $blockLen, $rPos, $motifLen, $aln, $templateSeq, $qSeq, $qPos,
		$matchLen,);

	my %read;
	$read{arr} = $arr;	
	$read{len} = length($arr->[9]);
	
	$tStart = $arr->[3];
	$cigar = $arr->[5];
	$clippedHead = 0;
	$clippedTail = 0;

	$cigar =~ s/^(\d+)I/$1S/;
	$cigar =~ s/(\d+)I$/$1S/;

	if($cigar =~ /^(\d+)[HS]/){
		$clippedHead = $1;
	}
	if($cigar =~ /(\d+)[HS]$/){
		$clippedTail = $1;
	}	

	$read{clippedHead} = $clippedHead;
	$read{clippedTail} = $clippedTail;

	if(($clippedHead != 0 && $clippedTail != 0) || $arr->[5] eq "*") { ### discard reads whose both ends were clipped.		
		return;
	}
	

	$seq = $arr->[9];
	$rvSeq = rc($arr->[9]);
	$qLen = length($seq);

	my $c = substr($arr->[10], 0, 1);
	$qualChar = $c if($c ne "*" && (!defined $qualChar || ord($qualChar) < ord($c)));

	if(!defined $refName || $refName eq ""){		
		getTargetRefSeq($arr->[2]);
		$refAddStart = 1;
		$newRefSeq = "";
		$arrPoints = $allPointsInRef{$arr->[2]};
		$pointNum = $#{$arrPoints} + 1;	
		$pI = 0;
	}
	elsif($refName ne $arr->[2]){
		while(defined $arrPoints->[$pI]){ ### to print all remaining assmbly points.
			findNewRefOnPoints($pI);
			$pI++;
		}
		
		printNewRefSeq();

		delete $allPointsInRef{$refName};

		$pI = 0;
		$prevTStart = 0;
		
		getTargetRefSeq($arr->[2]);
		$refAddStart = 1;
		$newRefSeq = "";
		$arrPoints = $allPointsInRef{$arr->[2]};
		$pointNum = $#{$arrPoints} + 1;
	}

	$refName = $arr->[2];
	if($tStart < $prevTStart){
		exitMsg(__LINE__, "Reads in the sam file is out of order.($arr->[0] start : $tStart < previous start : $prevTStart) They should be sorted.\n");
	}
	$prevTStart = $tStart;

	#### before reference start..		
	### the reads which are in least 5 bases distance from point. : $arrPoints->[$pI]{end}+5  < $tStart
	### ->
	### the reads should overlap at least 5 bases of point. : $arrPoints->[$pI]{end}-5  < $tStart
	while($pI < $pointNum && $arrPoints->[$pI]{end} < $tStart){ 
		if(!defined $arrPoints->[$pI]{skip}){				
			findNewRefOnPoints($pI);
		}
		$pI++;
	}
	return if($pI >= $pointNum);
	
	$qPos = $clippedHead + 1;
	$qEnd = $clippedHead;

	$rPos = $tStart;
	my $tEnd = $tStart -1;
	my (@arrBlock, $gapLen);
	$i = 0;
	while($cigar =~ /(\d+)(\w)/g){		
		$blockLen = $1;
		next if($2 eq "S" || $2 eq "H");

		if($2 eq "M"){
			$tEnd += $blockLen;
			$qEnd += $blockLen;

			$arrBlock[$i]{type} = $2;
			$arrBlock[$i]{tStart} = $rPos;
			$arrBlock[$i]{tEnd} = $tEnd;
			$arrBlock[$i]{qStart} = $qPos;
			$arrBlock[$i]{qEnd} = $qEnd;
			$arrBlock[$i]{len} = $blockLen;

			$i++;
			$rPos = $tEnd + 1;
			$qPos = $qEnd + 1;
		}
		elsif($2 eq "D" || $2 eq "N"){			
			$tEnd += $blockLen;			
			
			$arrBlock[$i]{type} = $2;
			$arrBlock[$i]{tStart} = $rPos;
			$arrBlock[$i]{tEnd} = $tEnd;
			$arrBlock[$i]{qStart} = $qPos;
			$arrBlock[$i]{qEnd} = $qPos;
			$arrBlock[$i]{len} = $blockLen;
			$i++;
			$rPos = $tEnd + 1;
		}
		elsif($2 eq "I"){
			$qEnd += $blockLen;

			if(!defined $rPos || !defined $arrPoints->[$pI]{end} || !defined $arrPoints->[$pI]{start}){
				exitMsg(__LINE__, "undefined .... if($rPos <= $arrPoints->[$pI]{end} && $rPos >= $arrPoints->[$pI]{start}) \n");
			}

			$arrBlock[$i]{type} = $2;
			$arrBlock[$i]{tStart} = $rPos;
			$arrBlock[$i]{tEnd} = $rPos;
			$arrBlock[$i]{qStart} = $qPos;
			$arrBlock[$i]{qEnd} = $qEnd;
			$arrBlock[$i]{len} = $blockLen;
			$i++;
			$qPos = $qEnd + 1;
		}
	}	
	$read{block} = \@arrBlock;
	$read{blockCnt} = $i;
	
	### the reads which are in least 5 bases distance from AP. : $arrPoints->[$pI]{start}-5 > $tEnd
	#->
	### the reads should overlap at least 5 bases of AP. : $arrPoints->[$pI]{start}+5 > $tEnd
	if(!defined $arrPoints->[$pI] || $arrPoints->[$pI]{start} > $tEnd){
		return;
	}


	$read{tStart} = $tStart;
	$read{tEnd} = $tEnd;
	$rPos = $tStart;
	
	### check whether another assembly point exist in the same read.  
	my $tmpI = $pI;
	while($tmpI < $pointNum){
		last if(!defined $tmpI || !defined $arrPoints->[$tmpI] || $arrPoints->[$tmpI]{start} > $tEnd);
				
		#### this process is to find a different INDEL length from the reference.
		$gapLen = 0;
		foreach my $block (@arrBlock) {
			if($block->{type} eq "D" && 
				$block->{tStart} <= $arrPoints->[$tmpI]{end} && 
				$block->{tEnd} >= $arrPoints->[$tmpI]{start})
			{
				$gapLen -= $block->{len};
			}
			elsif($block->{type} eq "I" && 
				$block->{tStart} <= $arrPoints->[$tmpI]{end} && 
				$block->{tEnd} >= $arrPoints->[$tmpI]{start})
			{
				$gapLen += $block->{len};
			}
		}			
		if($gapLen == 0 && $tStart < $arrPoints->[$tmpI]{start} -5 && $tEnd > $arrPoints->[$tmpI]{end} +5){
			$tmpI++;
			next;
		}
		###############
		
		$read{pointIDs}[ $#{$read{pointIDs}} + 1 ] = $tmpI;
		
		$arrPoints->[$tmpI]{reads}[$#{$arrPoints->[$tmpI]{reads}}+1] = \%read;
		$arrPoints->[$tmpI]{headClippedReads}[ $#{$arrPoints->[$tmpI]{headClippedReads}}+1 ] = \%read if($clippedHead != 0);
		$arrPoints->[$tmpI]{tailClippedReads}[ $#{$arrPoints->[$tmpI]{tailClippedReads}}+1 ] = \%read if($clippedTail != 0);
		$arrPoints->[$tmpI]{longGapRead_cnt}++ if(abs($gapLen) > 2); 
		$arrPoints->[$tmpI]{headClipped_cnt}++ if($clippedHead != 0);
		$arrPoints->[$tmpI]{tailClipped_cnt}++ if($clippedTail != 0);
	
		$tmpI++;
	}
}


sub getTargetRefSeq
{
	my $rName = shift;
	my $i;
	my $rID = $refIDs{$rName};
	for($i = $printedRefID+1; $i < $rID;$i++){
		readRefSeq($arrRefIDs[$i]);
		printFormattedSeq($arrRefIDs[$i], \$refSeq, 80);
	}
	readRefSeq($rName);
}

sub printNewRefSeq
{
	if($refAddStart == 1){
		printFormattedSeq($refName, \$refSeq, 80);
	}
	else{
		print "adding sequence, $refName : $refAddStart-" .length($refSeq). "\n";
		$newRefSeq .= substr($refSeq, $refAddStart - 1);
		print "length " . length($newRefSeq). "\n";
		printFormattedSeq($refName, \$newRefSeq, 80);
	}
	$printedRefID = $refIDs{$refName};
}

sub printFormattedSeq
{
	my ($rName, $seq, $lenLine) = @_;  ### $seq is sclar.. to reduce time to copy the sequence..

	my $i;
	my $len = length($$seq);
	print "printing $rName (length : $len NTs)..\n";
	print $out ">$rName\n";
	for($i = 0; $i < $len; $i+= $lenLine){
		print $out substr($$seq, $i, $lenLine) . "\n";
	}
}



#### $allPointsInRef{$ref}[$i]{start, end, len, motif, gid, loc, str}
sub findNewRefOnPoints{
	my ($pI) = @_;

	print "> $refName : $arrPoints->[$pI]{start}-$arrPoints->[$pI]{end}, point ID : $pI,  read count : ".($#{$arrPoints->[$pI]{reads}}+1) . 
		", head clipped : " . (!defined $arrPoints->[$pI]{headClipped_cnt} ? 0 : $arrPoints->[$pI]{headClipped_cnt}) .
		", tail clipped : " . (!defined $arrPoints->[$pI]{tailClipped_cnt} ? 0 : $arrPoints->[$pI]{tailClipped_cnt}) .
		", long gap reads : " . (!defined $arrPoints->[$pI]{longGapRead_cnt} ? 0 : $arrPoints->[$pI]{longGapRead_cnt}) .
		"\n" if($verboseFlag);
	
	### the first read of the point.
	my $read;	
	
	### if there is no clipped reads, then print it.
	if( 
		(!defined $arrPoints->[$pI]{headClipped_cnt} || $arrPoints->[$pI]{headClipped_cnt} < 2 || 
			!defined $arrPoints->[$pI]{tailClipped_cnt} || $arrPoints->[$pI]{tailClipped_cnt} < 2 ) &&
		(!defined $arrPoints->[$pI]{cnt_blatRead} || $arrPoints->[$pI]{cnt_blatRead} < 2 )&&
		(!defined $arrPoints->[$pI]{longGapRead_cnt} || $arrPoints->[$pI]{longGapRead_cnt} < 2 )
	){ ## number of reads >= 2 for cliped ends 
		print " No significant error read\n" if($verboseFlag);
		
		deleteAPReads($arrPoints->[$pI]);		
		return;
	}

	$arrPoints->[$pI]{readCnt} = $#{$arrPoints->[$pI]{reads}};

	my @arrReads = sort { $a->{tStart} <=> $b->{tStart} || $b->{clippedHead} <=> $a->{clippedHead}} @{$arrPoints->[$pI]{reads}};
	
	print "Aassembling a sequence from ".($#arrReads+1) . " reads..\n" if($verboseFlag);
	my $gapLen = assembleTarget($pI, \@arrReads);
	
	deleteAPReads($arrPoints->[$pI]);
	print "\n" if($verboseFlag);
}

sub assembleTarget
{
	my ($pI, $arrReads) = @_;

	## pointLen_1
	#######
	my ($tStart, $tEnd, $arr, $i, $j, $k,	$overLen, $len, $read);	
	my ($startRead, $endRead);
	my $seedLen = 30; ### min. required overlapping length between two reads to make a link/graph

	## read->{arr}[]
	## arr 0: ID, 1: flag, 2: reference, 3: rStart, 4: mapping score, 5: cigar, 6/7/8: mate info, 9: seq, 10: quality, 11..: tag.,

	makegraph($arrReads, 30);

	my ($sameCnt, $m, $n, $refStart, $refEnd);
	my $mergedSeq = "";		
	

	$j = int($#$arrReads/2);
	$k = 0;
	my $max = 0;
	## find a starting read	
	for($i = 0; $i < $j && ($i < 50 || $max ==0); $i++) {
		
		if($max < $arrReads->[$i]{rightScore}) {
			$max = $arrReads->[$i]{rightScore};
			$k = $i;
		}
	}
	$i = $k;
	

	if($max < 3){		
		print "Read sequences are not consistent. They may be from mis-mapping, chimeric DNAs or sequencing errors.\n" if($verboseFlag);#exit(1);
		return 0;
	}

	my @arrLinkeReadID;
	$arrLinkeReadID[0] = $i;
	$startRead = $arrReads->[$i];
	$arrReads->[$i]{mergeStart} = 1;
	$arrReads->[$i]{mergeEnd} = $arrReads->[$i]{len}; 
	#print STDERR "@{$arrReads->[$i]{arr}}, $arrReads->[$i]{len}\n"; exit;

	my ($isFirst) = (0);
	$refStart = $arrReads->[$i]{arr}[3];

	while($i < $#$arrReads) {
		$max = 0;
		
		$m = 0;
		for($j = 0;$j < $arrReads->[$i]{rightLinkNum} && $j < 10; $j++) {
			$k = $arrReads->[$i]{rightLink}[$j]{id};
			next if($arrReads->[$i]{tEnd} > $arrReads->[$k]{tEnd});
			if($max < $arrReads->[$k]{rightScore}){
				$max = $arrReads->[$k]{rightScore};
				$m = $j;
			}
		}
		$j = $m;
		$k = $arrReads->[$i]{rightLink}[$j]{id};	
				
		#### is it possible??
		if($j != 0 && $j == $arrReads->[$i]{rightLinkNum}){
			print "#### is it possible??\n";
			$mergedSeq .= $arrReads->[$k]{arr}[9];
			last;
		}

		$arrLinkeReadID[$#arrLinkeReadID+1] = $k;
		
		$arrReads->[$i]{mergeStart} = length($mergedSeq)+1;
		$arrReads->[$i]{mergeEnd} = $arrReads->[$i]{mergeStart}+$arrReads->[$i]{len}-1;
		$mergedSeq .= substr($arrReads->[$i]{arr}[9], 0, $arrReads->[$i]{rightLink}[$j]{lPos})
			if($arrReads->[$i]{rightLink}[$j]{lPos} != 0);
		
		
		### if it is the last read.
		if($arrReads->[$k]{rightLinkNum} == 0){		
			##### set the last read in the link as the end read.
			$arrReads->[$k]{mergeStart} = length($mergedSeq)+1;
			$arrReads->[$k]{mergeEnd} = $arrReads->[$k]{mergeStart}+$arrReads->[$k]{len}-1;
			$mergedSeq .= $arrReads->[$k]{arr}[9];

			##### or set the second last read in the link as the end read, because the last read may have a sequencing error at the right end
			last;
		}
		$i = $k;		

	}
	$refEnd = $arrReads->[$k]{tEnd};
	##### set the last read in the link as the end read.
	$endRead = $arrReads->[$k];

	if($startRead->{arr}[4] <= 5 && $endRead->{arr}[4] <= 5){ ### if the mapping quality of both end reads are low
		print "The mapping qualities of reads were too low.\n" if($verboseFlag);#exit(1);
		return 0;
	}
	print "Total number of reads : ".($#$arrReads+1).", number of reads in a link : ".($#arrLinkeReadID+1)."\n" if($verboseFlag);

	my $bStop = 0;

	my ($pointStart, $pointEnd, $block, $pointTStart, $rightDel);
	for($i = 0; $i < $startRead->{blockCnt}; $i++){ #each $block ($startRead->{block}) 
		$block = $startRead->{block}[$i];
		
		if($block->{tStart} <=  $arrPoints->[$pI]{start} && $block->{tEnd} >=  $arrPoints->[$pI]{start}) 
		{
			if($block->{type} eq "D"){
				$pointStart = $block->{qStart};
				$pointTStart = $block->{tStart};
			}
			else{
				$pointStart = $block->{qStart} + ($arrPoints->[$pI]{start}-$block->{tStart}); ### qPos of point at the read
				$pointTStart = $arrPoints->[$pI]{start};
			}

			## $pointStart is the point starting postion in a read.

			last;
		}
	}

	if(!defined $pointStart || $pointStart-1 < 20){ ### $pointStart-1 < 20 : at least 20 bases at the left flanking ...
		if(($arrPoints->[$pI]{headClipped_cnt}+$arrPoints->[$pI]{headClipped_cnt}) > ($#$arrReads+1)*20 && searchInsertedSequence($pI) ){
			return 1;
		}
		else{
			print "The merged sequence does not cover the whole target sequence enough.\n" if($verboseFlag);#exit(1);
			return 0;
		}
	}

	my $pointTEnd; ### in case the point end position is adjusted..
	for($i = $endRead->{blockCnt} -1; $i >= 0; $i--){ #each $block ($endRead->{block}) 
		$block = $endRead->{block}[$i];
		
		if($block->{tStart} <=  $arrPoints->[$pI]{end} && $block->{tEnd} >=  $arrPoints->[$pI]{end}) 
		{
			if($block->{type} eq "D"){
				$pointEnd = $endRead->{block}[$i-1]->{qEnd};
				$pointTEnd = $block->{tEnd};			
			}
			elsif($block->{type} eq "I"){	
				next;				
			}
			else
			{
				if($i < $endRead->{blockCnt}-1 && $block->{tEnd} == $arrPoints->[$pI]{end} && $endRead->{block}[$i+1]{type} eq "I"){
					$pointEnd = $endRead->{block}[$i+1]{qEnd};
				}
				else{
					$pointEnd = $block->{qStart} + ($arrPoints->[$pI]{end}-$block->{tStart}); ### qPos at $arrPoints->[$pI]{start}
				}
				$pointTEnd = $arrPoints->[$pI]{end};
			}
			last;
		}
	}
		
	
	my $gap = 0;
	### current $pointEnd is the end position of the point in the endRead.
	if(!defined $pointEnd || ($endRead->{len}-$pointEnd) < 20){ ### $endRead->{len}-$pointEnd < 20 : at least 20 bases at the right flanking ...
		if(($arrPoints->[$pI]{headClipped_cnt}+$arrPoints->[$pI]{headClipped_cnt}) > ($#$arrReads+1)*20 && searchInsertedSequence($pI) ){
			return 1;
		}
		else{
			print "The merged sequence does not cover the flanking sequences of the target enough.\n"  if($verboseFlag);#exit(1);		
			return 0;
		}
	}
	
	### convert the end position of the point in the endRead to the position in the mergedSeq
	$pointEnd = length($mergedSeq) -($endRead->{len}-$pointEnd);
		
	my $assembleLen = ($pointEnd - $pointStart +1);

	$gap = $assembleLen - ($pointTEnd-$pointTStart+1);

	if($assembleLen <= 0){
		$gap = searchInsertedSequence($pI);
		#print "\nThe merged sequence is too much different from the reference.\n" if($verboseFlag);
		print "\nThe merged sequence is invalid.\n" if($verboseFlag && $gap == 0);
	}
	elsif(abs($gap) <= 2){
		$gap = searchInsertedSequence($pI);
		print "Gap length of the assembled sequence is $gap bases. Too small difference.\n" if($verboseFlag && $gap == 0);
	}
	else{
		print "original sequence : " . substr($refSeq, $pointTStart-1, ($pointTEnd-$pointTStart+1)) . "\n" if($verboseFlag);
		my $newSeq = substr($mergedSeq, $pointStart -1, $assembleLen);
		$newRefSeq .= substr($refSeq, $refAddStart - 1, ($pointTStart-$refAddStart)) . $newSeq;
		print "replacing sequence, $refName : $refAddStart-$pointTStart + $newSeq\n" if($verboseFlag);
		$refAddStart = $pointTEnd+1;
		
	}

	if($gap == 0 && ($arrPoints->[$pI]{headClipped_cnt}+$arrPoints->[$pI]{headClipped_cnt}) > ($#$arrReads+1)*20 ){
		$gap = searchInsertedSequence($pI);
	}
	
	#if($verboseFlag){
	#	for($i = 0; $i <= $#$arrReads; $i++) {
	#		print "$i : $arrReads->[$i]{arr}[0], $arrReads->[$i]{arr}[9], t:$arrReads->[$i]{tStart}-$arrReads->[$i]{tEnd}, ". 
	#			(defined $arrReads->[$i]{mergeStart} ? "m:$arrReads->[$i]{mergeStart}-$arrReads->[$i]{mergeEnd}, " : "") .
	#			"(left $arrReads->[$i]{leftLinkNum}, right $arrReads->[$i]{rightLinkNum})\n";
	#	}
	#	print "link : "; for($i = 0; $i <= $#arrLinkeReadID; $i++){ print "$arrLinkeReadID[$i], ";}	print "\n";
	#}
	
	deleteReadLink($arrReads);

	return $gap;
}


### $seedLen = 30; ### min. required overlapping length between two reads to make a link/graph
sub makegraph
{
	my ($arrReads, $seedLen, $isRev) = @_;  ### if $isRev == 1, make score from right to left

	my ($short, $continueNoHit, $pos, $i, $j, $read);
	

	deleteReadLink($arrReads);

	$arrReads->[0]{leftLinkNum} = 0;
	$arrReads->[0]{rightLinkNum} = 0;
	$continueNoHit = 0; 
	my $refName = $arrReads->[0]{arr}[2];
	
	for($i = 1; $i <= $#$arrReads; $i++) {
		$read = $arrReads->[$i];
		$read->{leftLinkNum} = 0 if(!defined $read->{leftLinkNum}) ;
		$read->{rightLinkNum} = 0 if(!defined $read->{rightLinkNum});
		next if($read->{len} <= $seedLen);
	
		$continueNoHit = 0;
		$short = substr($read->{arr}[9], 0, $seedLen);
		for($j = $i - 1; $j >= 0; $j--){
			next if($arrReads->[$j]{len} <= $seedLen);
			$pos = -1;
						
			while(($pos = index($arrReads->[$j]{arr}[9], $short, $pos+1)) > -1){					
				if(!defined $arrReads->[$j]{arr}[9] || length($arrReads->[$j]{arr}[9]) < $pos+$seedLen || !defined $read->{arr}[9] || $arrReads->[$j]{len} < $pos){
					print "short $short, $arrReads->[$j]{arr}[0], $arrReads->[$j]{arr}[9], $pos+$seedLen,  $read->{arr}[0], $read->{arr}[9] \n"; exit;
				}
				last if(substr($arrReads->[$j]{arr}[9], $pos+$seedLen) eq substr($read->{arr}[9], $seedLen, $arrReads->[$j]{len} - $pos - $seedLen));
			}

			if($pos != -1){				
				$read->{leftLink}[$read->{leftLinkNum}]{id} = $j;
				$read->{leftLink}[$read->{leftLinkNum}]{lPos} = $pos; ### position of the right side read in the left side read
				$arrReads->[$j]{rightLink}[$arrReads->[$j]{rightLinkNum}]{id} = $i;
				$arrReads->[$j]{rightLink}[$arrReads->[$j]{rightLinkNum}]{lPos} = $pos;
				$read->{leftLinkNum}++;
				$arrReads->[$j]{rightLinkNum}++;
				$continueNoHit = 0;
			}
			else{ ### if number of continous no hits is over 20, no connection..
				$continueNoHit++;
				last if($continueNoHit >= 100)
			}
		}			
	}
	
	if($isRev){
		for($i = 0; $i <= $#$arrReads; $i++) {
			$arrReads->[$i]{leftScore} = $arrReads->[$i]{leftLinkNum} if(!defined $arrReads->[$i]{leftScore});
			for($j = 0;$j < $arrReads->[$i]{rightLinkNum}; $j++) {
				$arrReads->[ $arrReads->[$i]{rightLink}[$j]{id} ]{leftScore} += $arrReads->[$i]{leftScore} + $arrReads->[$i]{leftLinkNum};
			}
		}
	}
	else{
		for($i = $#$arrReads; $i >= 0; $i--) {
			$arrReads->[$i]{rightScore} = $arrReads->[$i]{rightLinkNum} if(!defined $arrReads->[$i]{rightScore});
			for($j = 0;$j < $arrReads->[$i]{leftLinkNum}; $j++) {
				$arrReads->[ $arrReads->[$i]{leftLink}[$j]{id} ]{rightScore} += $arrReads->[$i]{rightScore} + $arrReads->[$i]{rightLinkNum};
			}
		}
	}
}

sub searchInsertedSequence
{
	my $pI = shift;

	return 0 if($arrPoints->[$pI]{headClipped_cnt} < $arrPoints->[$pI]{readCnt}*0.1 || $arrPoints->[$pI]{tailClipped_cnt} < $arrPoints->[$pI]{readCnt}*0.1);

	my @arrReads = sort { $a->{tStart} <=> $b->{tStart}} @{$arrPoints->[$pI]{tailClippedReads}};
	makegraph(\@arrReads, 40, 1);

	my ($i, $j, $k);
	$j = int($#arrReads/2);
	$k = 0;
	my $max = 0;
	## find a starting read	
	for($i = $#arrReads; $i > $j && ($i > $#arrReads-50 || $max ==0); $i--) {
		
		if($max < $arrReads[$i]{leftScore}) {
			$max = $arrReads[$i]{leftScore};
			$k = $i;
		}
	}
	$i = $k;

	return 0 if(!defined $arrReads[$i]{leftLink}[0] || ($#{$arrReads[$i]{leftLink}}+1) < ($#arrReads+1)*0.3);
	
	my $rightRead = $arrReads[$arrReads[$i]{leftLink}[0]{id}];
	print "Total tail clipped reads : ".($#arrReads+1).
		", Left link : ".($#{$arrReads[$i]{leftLink}}+1)."\n" if($verboseFlag);
	for($j = $#{$arrReads[$i]{leftLink}}; $j >= 0; $j--){
		$k = $arrReads[$i]{leftLink}[$j]{id};
		#print "$arrReads[$k]{arr}[0], $arrReads[$k]{arr}[3], $arrReads[$k]{arr}[5], $arrReads[$k]{arr}[9]\n" if($verboseFlag);
		if($rightRead->{clippedTail} < $arrReads[$k]->{clippedTail}){
			$rightRead = $arrReads[$k];
		}
	}

	my $insert5prime; ## 5' of insertion..
	
	return 0 if($rightRead->{clippedTail} < 20);
	$insert5prime = substr($rightRead->{arr}[9], length($rightRead->{arr}[9]) - $rightRead->{clippedTail});

	

	@arrReads = sort { $a->{tStart} <=> $b->{tStart} || $b->{clippedHead} <=> $a->{clippedHead}} @{$arrPoints->[$pI]{headClippedReads}};
	makegraph(\@arrReads, 40);

	$j = int($#arrReads/2);
	$k = 0;
	$max = 0;
	## find a starting read	
	for($i = 0; $i < $j && ($i < 50 || $max ==0); $i++) {
		
		if($max < $arrReads[$i]{rightScore}) {
			$max = $arrReads[$i]{rightScore};
			$k = $i;
		}
	}
	$i = $k;
	return 0 if(!defined $arrReads[$i]{rightLink}[0] || ($#{$arrReads[$i]{rightLink}}+1) < ($#arrReads+1)*0.3);

	my $leftRead = $arrReads[$arrReads[$i]{rightLink}[0]{id}]; ### considering that $arrReads[$i] read may have errors at the end, choose its right linked read.
	print "\nTotal head clipped reads : ".($#arrReads+1).
		", Right link : ".($#{$arrReads[$i]{rightLink}}+1)."\n" if($verboseFlag);
	for($j = 0; $j <= $#{$arrReads[$i]{rightLink}}; $j++){
		$k = $arrReads[$i]{rightLink}[$j]{id};
		#print "$arrReads[$k]{arr}[0], $arrReads[$k]{arr}[3], $arrReads[$k]{arr}[5], $arrReads[$k]{arr}[9]\n" if($verboseFlag);
		if($leftRead->{clippedHead} < $arrReads[$k]->{clippedHead}){
			$leftRead = $arrReads[$k];
		}
	}

	my $insert3prime; ## 3' of insertion..
	return 0 if($leftRead->{clippedHead} < 20);
	$insert3prime = substr($leftRead->{arr}[9], 0, $leftRead->{clippedHead});

	print "searching sequence : $insert5prime  ....  $insert3prime\n"  if($verboseFlag);

	
	my ($selStart, $selEnd, $start, $end) = (0,0,0,0);
	$i = index($refSeq, $insert5prime);
	$j = index($refSeq, $insert3prime, $i+1);

	while($i != -1 && $j != -1){
		$start = $i;
		$end = $j;
		if($selStart == 0 || ($end-$start) < ($selEnd-$selStart) ){
			($selEnd, $selStart) = ($end, $start);
		}
		$i = index($refSeq, $insert5prime, $i+1);
		$j = index($refSeq, $insert3prime, $i+1) if($i > $j);		
	}

	my $bRev = 0;
	$i = index($revRefSeq, $insert5prime);
	$j = index($revRefSeq, $insert3prime, $i+1);
	while($i != -1 && $j != -1){
		$start = $i;
		$end = $j;
		if($selStart == 0 || ($end-$start) < ($selEnd-$selStart) ){
			($selEnd, $selStart) = ($end, $start);
			$bRev = 1;
		}
		$i = index($revRefSeq, $insert5prime, $i+1);
		$j = index($revRefSeq, $insert3prime, $i+1) if($i > $j);		
	}

	my ($insSeq, $insLen);
	$insLen = ($selEnd-$selStart+$leftRead->{clippedHead});
	if($selStart == 0 || $insLen > 2000) ## minimum size of the sequence to insert is 2000;
	{
		print "No proper sequence to be inserted\n"  if($verboseFlag);		
		return -1;
	}

	if($bRev){
		$insSeq = substr($revRefSeq, $selStart, $insLen);		
	}
	else{
		$insSeq = substr($refSeq, $selStart, $insLen);
		print "".($selStart+1). "-" . ($selEnd+1+$leftRead->{clippedHead}) . "\n";
	}	
	
	$newRefSeq .= substr($refSeq, $refAddStart - 1, ($rightRead->{tEnd}+1-$refAddStart)) . $insSeq;
	$refAddStart = $leftRead->{arr}[3];
	print "inserted sequence : $insSeq\n" if($verboseFlag);


	return $insLen;
}

sub deleteReadLink
{
	my $arrReads = shift;
	
	for(my $i = 0; $i <= $#$arrReads; $i++){
		delete $arrReads->[$i]{leftLink};
		delete $arrReads->[$i]{rightLink};
		delete $arrReads->[$i]{leftLinkNum};
		delete $arrReads->[$i]{rightLinkNum};
		delete $arrReads->[$i]{leftScore};
		delete $arrReads->[$i]{rightScore};
		delete $arrReads->[$i]{mergeStart};
		delete $arrReads->[$i]{mergeEnd};
	}
}


sub deleteAPReads {
	my $point = shift;

	deleteReadLink($point->{reads});
	delete $point->{reads};
	delete $point->{headClippedReads};
	delete $point->{tailClippedReads};
}

sub readRefSeq{
	my $refName = shift;

	my ($name);
	print "reading '$refName' sequence...\n";
	$rBuf = <$inR> if(!$rBuf); 
	while(defined $rBuf){
		if($rBuf =~ /^>/){			
			$name = $';
			$name =~ s/^\s*//g; $name =~ s/\s.*//g;
			if($name eq $refName){
				$rBuf = <$inR>;
				last;
			}
		}
		$rBuf = <$inR>;
	}
	
	$refSeq = "";
	
	while(defined $rBuf){
		last if($rBuf =~ /^>/);
		$rBuf =~ s/[\r\n]+//g;
		$refSeq .= uc $rBuf;
		$rBuf = <$inR>;
	}

	if($refSeq eq ""){
		exitMsg(__LINE__, "Could not read sequences for '$refName'. The order of references in '$refFn' and the sam may be differ.\n");
	}
	$revRefSeq = rc($refSeq);
	$refLen = length($refSeq);
}




sub findAssemblyPoints{
	print "seraching assembly points from $refName\n" if($verboseFlag);
	my ($i, $j, $prevJ, $r, $k, $rNum);

	$k = 0; ### index of an assembly point
	$r = 0;
	$rNum = $#$tandemRepeatList+1;  ### repeat num;
	
	for($i = 2; $i < $hashRefLen{$refName}; $i++){
		if($leftClipped[$i] && $leftClipped[$i] >= $numClipped){
			while($r < $rNum && $tandemRepeatList->[$r]{end} < $i){
				$r++;
			}
			print "$i, left clipped reads : $leftClipped[$i]\n" if($verboseFlag);
			$prevJ = 0;
			for($j = $i; $j < $hashRefLen{$refName} && $j < $i+$readLen-10; $j++){
				if($rightClipped[$j] && $rightClipped[$j] >= $numClipped){
					print "$j, right clipped reads : $rightClipped[$j]\n" if($verboseFlag);
					$prevJ = $j;
				}
			}
			if($prevJ != 0){
				my ($start, $end) = ($i, $prevJ);
				if($r < $rNum && $tandemRepeatList->[$r]{end} >= $i && $tandemRepeatList->[$r]{start} <= $prevJ){
					$start = $tandemRepeatList->[$r]{start} if($start > $tandemRepeatList->[$r]{start});
					$end = $tandemRepeatList->[$r]{end} if($end < $tandemRepeatList->[$r]{end});
				}  ### overlap to tandemRepeat..
				
				$allPointsInRef{$refName}[$k]{start} = $start;
				$allPointsInRef{$refName}[$k]{end} = $end;
				$allPointsInRef{$refName}[$k]{len} = $end - $start + 1;
				
				print "$refName\t$start:$i($leftClipped[$i])\t$end:$prevJ($rightClipped[$prevJ])\n" if($verboseFlag);
				$i = $end+1;

				$k++;
			}
			#print "i : $i, j : $j\n";
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


sub exitMsg
{
	my ($line, $msg) = @_;
	print STDERR "[Error] at the line $line : $msg\n" if($msg);
	exit(1);
}

sub openInput
{
	my ($fn) = @_;

	print "reading '$fn'\n";
	#return *STDIN unless defined $fn;

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




sub SW_align_NoEndsFree{
	my ($seqA, $seqB) = @_;	
	
	$seqA = uc $seqA;
	$seqB = uc $seqB;
	
	my ($lenA, $lenB);
	$lenA = length($seqA);
	$lenB = length($seqB);
	
	## the first column and the first row of matrix will have 0  
	## and next column and row will start with the first result..
	## so ...the result of position $A and $B will always be stored in $matrix[$B+1][$A+1]..	 
	my @matrix;
	  
	## it will store the information for the directions..
	## 0 from diagonal, 1 from above, 2 from left 
	my @moveFrom; 
	
	
	## @sa, @sb : arraies of sequences.
	## $A, $B : positions of sa and sb
	## $i, $j : matrix rows..
	## $val : temporary value
	my(@sa, @sb, $A, $B, $i, $j, $val, $matchScore, $misPenalty, $gapPenalty, $gapExtPenalty);
	@sa = split(//, $seqA);
	@sb = split(//, $seqB);
	
	$matchScore = 3;
	$misPenalty = -6;#-6;
	$gapPenalty = -8;#-8;
	$gapExtPenalty = -2; #-4;
	
	$matrix[0][0][0] = 0;
	$matrix[0][0][1] = 0;
	$matrix[0][0][2] = 0;

	for($B = 1; $B <= $lenB; $B++){ 
		$matrix[$B][0][0] = -1000000; 
		$matrix[$B][0][1] = -1000000;
		$matrix[$B][0][2] = -1000000;
	}
	for($A = 1; $A <= $lenA; $A++){
		$matrix[0][$A][0] = -1000000;
		$matrix[0][$A][1] = -1000000;	
		$matrix[0][$A][2] = -1000000; 
	}
	
	my @path = qw("\" "^" "<"); #### used to show the path..
	my ($max, $maxI, $penalty);
	
	$j = 0; ### matrix row of previous row of B. 
	for($B = 0; $B < $lenB; $B++){
				
		for($A = 0; $A < $lenA; $A++){
						
			# for the previous score....
			## At the first, the value from diagonal direction is calculated..	
			if($sa[$A] eq $sb[$B]) { 
				$penalty = $matchScore;
			}
			else{
				$penalty = $misPenalty;
			}

			$max = -1000000;
			for($i = 0; $i < 3; $i++){
				next if(!defined $matrix[$B][$A][$i]);
				#if(!defined $matrix[$B][$A][$i]){
				#	print "\$val = \$matrix[$B][$A][$i] + \$penalty =  $matrix[$B][$A][$i] + $penalty; \n"; exit(1);
				#}
				$val = $matrix[$B][$A][$i] + $penalty; 
				
				if($max < $val){
					$max = $val;
					$maxI = $i;
				}
			}		
			$matrix[$B+1][$A+1][0] = $max;
			$moveFrom[$B][$A][0] = $maxI;
			
			if(($A==0 && $B==0) || ($A==$lenA-1 && $B==$lenB-1)){ ## do not allow end free.				
				next;
			}

			## value from above direction	## one more base at seqB..
			$max = -1000000;
			for($i = 0; $i < 3; $i++){
				next if(!defined $matrix[$B][$A+1][$i]);
				$val = $matrix[$B][$A+1][$i] + ($B>1 && $i==1?$gapExtPenalty:$gapPenalty);
				$val -= 100 if($i == 2);
				if(($A > 0 && $sa[$A-1] eq $sa[$A] && $max <= $val) || $max < $val) {
				#if($max <= $val) {
					$max = $val;
					$maxI = $i;
				}
			}		
			$matrix[$B+1][$A+1][1] = $max;
			$moveFrom[$B][$A][1] = $maxI;
			
						
			## value from left direction  #### applying ends free alignment for B			
			$max = -1000000;
			for($i = 0; $i < 3; $i++){
				next if(!defined $matrix[$B+1][$A][$i]);
				$val = $matrix[$B+1][$A][$i] + ($A>1 && $i==2?$gapExtPenalty:$gapPenalty);
				$val -= 100 if($i == 1);
				if(($B > 0 && $sb[$B-1] eq $sb[$B] && $max <= $val) || $max < $val){
				#if($max <= $val){
					$max = $val;
					$maxI = $i;
				}
			}		
			$matrix[$B+1][$A+1][2] = $max;
			$moveFrom[$B][$A][2] = $maxI;
		}
	}
	
	my (@ssa, @ssb); # temporary sequence storages.
	
	my $res = new AlignRes();
	my (@variations, $insNum, $delNum, $delBases);
	my ($mismatchNum, $gapNum) = (0,0);
	($insNum, $delNum, $delBases) = (0,0, "");
	
	$max = -100000;
	$B = $lenB;
	$A = $lenA;
	for($i = 0; $i < 3; $i++){
		next if(!defined $matrix[$B][$A][$i]);
		if($max < $matrix[$B][$A][$i]){			
			$max = $matrix[$B][$A][$i];
			$maxI = $i;
		}
	}
	
	### for the gap position.....
	$i = 0; 
	$j = 0;
	
	### record from the last position to the first position
	for($B = $lenB-1, $A = $lenA -1; $B >=0 || $A >= 0; $i++, $j++){
		
		#print "$A, $B -- \n" if(!defined $moveFrom[$B][$A]);
		if($A >= 0 && $B >= 0 && $maxI == 0 ){ #from diagnol direction
			$ssa[$i] = $sa[$A];
			$ssb[$j] = $sb[$B];
		
			$res->{seqAMatchEnd} = $A+1 if(!defined $res->{seqAMatchEnd});  ### record only once... 
			$res->{seqBMatchEnd} = $B+1 if(!defined $res->{seqBMatchEnd});
			$res->{seqAMatchStart} = $A+1;
			$res->{seqBMatchStart} = $B+1;
			
			#print "$A, $B, $ssa[$i], $ssb[$j], $moveFrom[$B][$A], $matrix[$B+1][$A+1]\n"; 
			if($insNum != 0){
				my %node;
				$node{type} = "I";
				$node{pos} = $B+2;
				$node{aPos} = $A+2;
				$node{len} = $insNum;
				push(@variations, \%node);
			}
			elsif($delNum != 0){
				my %node;
				$node{type} = "D";
				$node{pos} = $B+2;
				$node{aPos} = $A+2;
				$node{len} = $delNum;
				$node{bases} = $delBases;
				push(@variations, \%node);
			}

			if($sa[$A] ne $sb[$B]){
				$mismatchNum++ ;
				my %node;
				$node{type} = "M";
				$node{pos} = $B+1;
				$node{aPos} = $A+1;
				$node{len} = 1;
				$node{bases} = $sa[$A];
				push(@variations, \%node);
			}
			$maxI = $moveFrom[$B][$A][$maxI];
			$A--;
			$B--; 
			
			($insNum, $delNum, $delBases) = (0,0, "");
		}
		elsif($B >=0 && ($A < 0 || $maxI == 1)){ #from above direction
		
			$ssa[$i] = "-";
			$ssb[$j] = $sb[$B];
			#print "$A, $B, $ssa[$i], $ssb[$j], $moveFrom[$B][$A], $matrix[$B+1][$A+1]\n"; 
			if($A >=0 && $A != $lenA-1){
				$gapNum++;
				$insNum++;
			}

			$maxI = ($A == -1 ? 1 : $moveFrom[$B][$A][$maxI]);
			$B--; 
		}
		else { #from left direction
		
			$ssa[$i] = $sa[$A];
			$ssb[$j] = "-";
			#print "$A, $B, $ssa[$i], $ssb[$j], $moveFrom[$B][$A], $matrix[$B+1][$A+1]\n"; 
			if($B >=0 && $B != $lenB-1){
				$gapNum++;
				$delNum++;
				$delBases = $sa[$A] . $delBases;
			}

			$maxI = ($B == -1 ? 2 : $moveFrom[$B][$A][$maxI]);
			$A--; 
		}	
			
	}
	
	$res->{mismatch} = $mismatchNum;
	$res->{gap} = $gapNum;

	@variations = reverse(@variations);	
	#@variations = sort {$a->{pos} <=> $b->{pos}} @variations;
	
	$res->{variations} = \@variations;

	### convert...
	for($A = $#ssa; $A >= 0; $A--){
		$sa[$#ssa-$A] = $ssa[$A];
	}
	for($B = $#ssb; $B >=0 ; $B-- ){
		$sb[$#ssb-$B] = $ssb[$B];
	}
	$res->{seqA} = join("", @sa);
	$res->{seqB} = join("", @sb);
	
	return $res;
}


