#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2011
#
# fix_bwasw_sam.p
# fix a alignment at the first base of a reference and a splitted read of which parts are overlapping each other.
######################################

use strict;
use warnings "all";
use Getopt::Long qw(:config no_ignore_case);


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


my ($cmd, $helpFlag, $verboseFlag, $outFn);

my ($samFn, $refFn);
$cmd = "$0 @ARGV";  ### $0 : command..  store all inputs(@ARGV) as a string in $cmd..

GetOptions(
	"h|?|help"	=> \$helpFlag,
	"verbose"		=> \$verboseFlag,
	"sam=s"	=> \$samFn,
	"ref=s"	=> \$refFn,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;


$samFn = shift if(!defined $samFn);
$refFn = shift if(!defined $refFn);
$outFn = shift if(!defined $outFn);

if(!defined $samFn){ print STDERR "\nNeed a sam or bam file!!\n\n"; help(1); }
if(!-e $samFn){ print STDERR "\nCould not find $samFn!!\n\n"; help(1); }
if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }
if(!-e $refFn){ print STDERR "\nCould not find $refFn!!\n\n"; help(1); }
#if(!defined $outFn){ print STDERR "\nNeed a output file!!\n\n"; help(1); }


if(defined $outFn && ($samFn eq $outFn || $refFn eq $outFn))
{ print STDERR " Error) Input and output files are same \n"; exit 1; }


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  [-s] <sam file>  [-r] <reference file>  [[-o] <out file>]\n";
	print STDERR "  ex) $0 -s test.sam -r sequence.fa -o test.fixed.sam\n\n";
	exit($return);
}
my ($in, $out, @arr, $i);

if(!defined $outFn){
	$outFn = $samFn;
	$outFn =~ s/(sam|bam)$/fixed.sam/;
}

my ($rBuf, $refSeq, $rin);
$rin = openInput($refFn);

my ($qName, %refLen, $len, $refName, $strand, %readAln, %usedReads);
open($in, $samFn =~ /.bam$/ ? "samtools view -h $samFn |"  : $samFn);
$refName = "";

my (@arrAlign, $emptyI, $isNewAdd, $flag);

$out = openOutput($outFn);
while(<$in>){
	s/[\r\n]+//g;
	if(/^@/ || /^\s*$/){ print $out "$_\n"; next; }

	## 0: ID, 1: flag, 2: reference, 3: position, 4: mapping score, 5: cigar, 6/7/8: mate info, 9: seq, 10: quality, 11..: tag.
	my @arr = split /\t/;
	my $ci = 5; ### cigar I..

	
	$len = 0;
	while($arr[$ci] =~ /(\d+)(\w)/g){				
		my $gapLen = $1;
		my $flag = $2;
#		if($flag eq "M" || $flag eq "D"){
#			$rPos += $gapLen;
#		}
		
		if($flag eq "M" || $flag eq "I"){
			$len += $gapLen;
		}
	}
	#print "cigar len : $len, seq len : " . length($arr[9]) . "\n"; exit;

	$strand = (($arr[1]&0x0010) == 0 ? "+" : "-"); 

	$i= (($arr[1]&0x0080) == 0 ? 1 : 2); ####($arr[1]&0x0040)/64 + ($arr[1]&0x0080)/64;		
	$qName = $arr[0] . ($arr[0] =~ /#0$/ || $arr[1]&0x0001 ? "/$i" : "");	

	if($arr[2] ne $refName){
		printAllReads();
		last if($arr[2] eq "*");
		readRefSeq($arr[2]);
		$refName = $arr[2];
		#if($refName =~ /^(\d+)$/ || $refName =~ /^[MXY]$/){	$refName = "chr" . $refName; }
	}

	if($arr[3] == 2 && $arr[$ci] =~ /(\d+)S(\d+)M/){
		my ($clipped, $match) = ($1, $2);

		my ($md, $nm, $str);
		$clipped--;
		$match++;
		$str = ($clipped == 0 ? $match . "M" : $clipped . "S$match" . "M"); 
		$arr[$ci] =~ s/\d+S\d+M/$str/;
		$arr[3] = 1;
	}

	my ($qStart, $qEnd, $rStart, $rEnd, $gapLen);
	$qStart = 1;
	$qEnd = 0;
	$rStart = $arr[3];
	$rEnd = $rStart -1;

	while($arr[$ci] =~ /(\d+)(\w)/g){		
		$gapLen = $1; $flag = $2;
		if($2 eq "M"){
			$rEnd += $gapLen;			
			$qEnd += $gapLen;
		}
		elsif($2 eq "D" || $2 eq "N"){
			$rEnd += $gapLen;			
		}
		elsif($2 eq "I"){
			$qEnd += $gapLen;
		}
		elsif($2 eq "S" && $qEnd == 0){
			$qStart = $gapLen + 1;
			$qEnd = $gapLen;
		}
	}

	#foreach my $name (keys %read) {
	#	if($readAln{$name}{end} < 
	#}

	if(!defined $readAln{$qName}){
		$readAln{$qName}{seq} = $arr[9];
		$readAln{$qName}{qual} = $arr[10];
	}
	
	$readAln{$qName}{arr}[$#{$readAln{$qName}{arr}}+1] = \@arr;	
	
	
	if($arr[$ci] =~ /M$/ || ($arr[$ci] =~ /(\d+)S$/ && $1 < 50))   ### it may be the last segment of the read
	{ 
		if($#{$readAln{$qName}{arr}} == 0){
			print "Printing $qName\n" if($verboseFlag);
			print $out join("\t", @arr). "\n"; # exit;
		}
		else{			
			print "$qName has " . ($#{$readAln{$qName}{arr}}+1) . " fragments.\n" if($verboseFlag);
			mergeAlignment(\%{$readAln{$qName}}, $qName);
		}
		delete $readAln{$qName};
		$usedReads{$qName} = 1;
	}
	else{		
		$arr[9] = "";
		$arr[10] = "";
		$readAln{$qName}{qEnd} = $qEnd;
		$readAln{$qName}{rEnd} = $rEnd;
	}	
}

printAllReads();

close($in);
close($out);


sub printAllReads{
	foreach $qName (keys %readAln) {		
		mergeAlignment(\%{$readAln{$qName}}, $qName);
		delete $readAln{$qName};
	}	
}

sub mergeAlignment{
	my ($readAll, $qName) = @_;
		
	my (%prev, $strand,  $i, $c,
		$qStart, $qEnd, $rStart, $rEnd, $gapLen, $ci, $prevCigar, $postCigar);
	$ci = 5; #cigar 

	$i = 0;

	print "Printing $qName\n" if($verboseFlag);

	
	my $ii = 0;
	my $qLen = length($readAll->{seq});
	### sam format
	## 0: ID, 1: flag, 2: reference, 3: position, 4: mapping score, 5: cigar, 6/7/8: mate info, 9: seq, 10: quality, 11..: tag.
	foreach my $aln (@{$readAll->{arr}}) {
		#print "$ii , rPos : $aln->[3]\n"; $ii++;

		$strand = (($aln->[1]&0x0010) == 0 ? "+" : "-"); 

		$qStart = 1;
		$qEnd = 0;
		$rStart = $aln->[3];
		$rEnd = $rStart -1;

		while($aln->[$ci] =~ /(\d+)(\w)/g){
			$gapLen = $1; $flag = $2;
			if($2 eq "M"){
				$rEnd += $gapLen;			
				$qEnd += $gapLen;
			}
			elsif($2 eq "D" || $2 eq "N"){
				$rEnd += $gapLen;			
			}
			elsif($2 eq "I"){
				$qEnd += $gapLen;
			}
			elsif($2 eq "S" && $qEnd == 0){
				$qStart = $gapLen + 1;
				$qEnd = $gapLen;
			}
		}		

		if(defined $prev{strand} && $prev{strand} ne $strand){
			if($rEnd <= $prev{rEnd}){
				next;
			}
			my @arrCigar;

			cigar2Array($prev{arr}[$ci], \@arrCigar);
			$prev{arr}[$ci] = array2Cigar(\@arrCigar, $prev{rStart}, \$refSeq, 1, \$readAln{$qName}{seq});
			print $out join("\t", @{$prev{arr}}) . "\n";
			%prev = ();
		}

		
=head		
		print "query $qStart - $qEnd, ref $rStart - $rEnd\n"; 
		if(defined $prev{strand} && $qStart - $prev{qEnd} > 0 && $rStart - $prev{rEnd} > 0 ){
			print "-  " .($qStart - $prev{qEnd}). " -- " .($rStart - $prev{rEnd}). " -\n";
		}
		$aln->[9] = $readAll->{seq};
		$aln->[10] = $readAll->{qual};
		$prev{strand} = $strand;
		$prev{arr} = $aln;
		$prev{qStart} = $qStart;
		$prev{qEnd} = $qEnd; 
		$prev{rStart} = $rStart; 
		$prev{rEnd} = $rEnd;			
		next;
=cut

		my $bGoNext = 0;

		if(!defined $prev{strand}){
			$bGoNext = 1;
		}
		elsif($prev{rEnd} < $rStart && $prev{qEnd} == $qStart-1){
			$gapLen = $rStart - $prev{rEnd} -1;
			$prev{arr}[$ci] =~ s/\d+S$//;
			$aln->[$ci] =~ s/^\d+S//;
			$aln->[$ci] = $prev{arr}[$ci] . ($gapLen == 0 ? "" : $gapLen . "D") . $aln->[$ci];
			$aln->[3] = $prev{rStart};
			$rStart = $prev{rStart};
			$bGoNext = 1;
			#print "+++ 1\n";
		}
		elsif($prev{rEnd} == $rStart -1 && $prev{qEnd} < $qStart){
			$gapLen = $qStart - $prev{qEnd} -1;
			$prev{arr}[$ci] =~ s/\d+S$//;
			$aln->[$ci] =~ s/^\d+S//;
			$aln->[$ci] = $prev{arr}[$ci] . $gapLen . "I" . $aln->[$ci];
			$aln->[3] = $prev{rStart};
			$rStart = $prev{rStart};
			$bGoNext = 1;
			#print "+++ 2\n";
		}
		elsif($prev{rEnd} < $rStart && $prev{qEnd} == $qStart){			
			$gapLen = $rStart - $prev{rEnd};
			$prev{arr}[$ci] =~ s/\d+S$//;
			$aln->[$ci] =~ s/^\d+S//;
			$aln->[$ci] =~ /^(\d+)M/;
			$c = ($1-1) . "M";
			$aln->[$ci] =~ s/^\d+M/$c/;

			$aln->[$ci] = $prev{arr}[$ci] . ($gapLen == 0 ? "" : $gapLen . "D") . $aln->[$ci];
			$aln->[3] = $prev{rStart};
			$rStart = $prev{rStart};
			$bGoNext = 1;
			#print "+++ 3\n";
		}
		elsif($prev{rEnd} == $rStart && $prev{qEnd} < $qStart){
			$gapLen = $qStart - $prev{qEnd};
			$prev{arr}[$ci] =~ s/\d+S$//;
			$aln->[$ci] =~ s/^\d+S//;
			$aln->[$ci] =~ /^(\d+)M/;
			$c = ($1-1) . "M";
			$aln->[$ci] =~ s/^\d+M/$c/;

			$aln->[$ci] = $prev{arr}[$ci] . ($gapLen == 0 ? "" : $gapLen . "I") . $aln->[$ci];
			$aln->[3] = $prev{rStart};
			$rStart = $prev{rStart};
			$bGoNext = 1;
			#print "+++ 4\n";
		}
		elsif($rEnd <= $prev{rEnd}){
			next;
		}
		elsif(($rStart - $prev{rEnd} > 20 && $qStart - $prev{qEnd} > 20) || ($qStart < $prev{qStart})
			|| $rStart - $prev{rEnd} > 100000 || $qStart - $prev{qEnd} > 100000){
			my @arrCigar;
			cigar2Array($prev{arr}[$ci], \@arrCigar);
			$prev{arr}[$ci] = array2Cigar(\@arrCigar, $prev{rStart}, \$refSeq, 1, \$readAln{$qName}{seq});
			print $out join("\t", @{$prev{arr}}) . "\n";
			$bGoNext = 1;
		}

		if($bGoNext){			
			$aln->[9] = $readAll->{seq};
			$aln->[10] = $readAll->{qual};
			$prev{strand} = $strand;
			$prev{arr} = $aln;
			$prev{qStart} = $qStart;
			$prev{qEnd} = $qEnd; 
			$prev{rStart} = $rStart; 
			$prev{rEnd} = $rEnd;	
			
			#print "--next $rStart; \n";
			next;
		}

		my ($alnQStart, $alnQEnd, $alnRStart, $alnREnd, $qPos, $rPos) = (0,0,0,0,0,0);
		my $diff = 0;
		
		if($prev{rEnd} < $rStart && $prev{qEnd} < $qStart){
			($alnQStart, $alnQEnd, $alnRStart, $alnREnd) = 
				($prev{qEnd}+1, $qStart -1, $prev{rEnd}+1, $rStart-1);
			$prevCigar = $prev{arr}[$ci];
			$prevCigar =~ s/\d+S$//;
			#$prevCigar =~ /(\d+)M$/;
			#$c = ($1-1)."M";
			#$prevCigar =~ s/(\d+)M$/$c/;
			$postCigar = $aln->[$ci];
			$postCigar =~ s/^\d+S//;
			
			#print "if($prev{rEnd} < $rStart && $prev{qEnd} < $qStart)\n($prev{qEnd}+1  = $alnQStart, $alnQEnd, $alnRStart, $alnREnd)\n-----------------------\n";
		}
		else { #if($prev{rEnd} > $rStart){# && $prev{qEnd} <= $qStart){
			print "$qName, q prev($prev{qStart},$prev{qEnd}), current($qStart,$qEnd), r prev($prev{rStart},$prev{rEnd}), current($rStart,$rEnd)\n"  if($verboseFlag);

			my ($qLeft, $rLeft) = ($qStart, $rStart);
			#$qLeft = $prev{qEnd} - 100 if($prev{qEnd} - $qStart > 100);
			#$rLeft = $prev{rEnd} - 100 if($prev{rEnd} - $rStart > 100);

			#print "-- prePos $prev{rStart} $prev{arr}[$ci]\n";
			##### for the previous segment
			$qPos = 1;
			$rPos = $prev{rStart};
			$prevCigar = "";
			while($prev{arr}[$ci] =~ /(\d+)(\w)/g){
				$gapLen = $1; $flag = $2;
				if($flag eq "M"){					
					$rPos += $gapLen;
					$qPos += $gapLen;
					$alnRStart = $rPos;
					$alnQStart = $qPos; #print "alnQStart $alnQStart -- \n";

					if(($rPos - $rLeft) > ($qPos - $qLeft) && $rPos >= $rLeft){
						$diff = ($rPos - $rLeft);
						$prevCigar .= ($gapLen-$diff) . $flag;						
						$alnRStart -= $diff;
						$alnQStart -= $diff;
						#print "alnQStart $alnQStart, diff $diff -\n";
						last;
					}
					elsif($qPos >= $qLeft){
						$diff = ($qPos - $qLeft);
						$prevCigar .= ($gapLen-$diff) . $flag;						
						$alnRStart -= $diff;
						$alnQStart -= $diff;
						#print "alnQStart $alnQStart, diff $diff +\n";
						last;
					}
				}
				elsif($flag eq "D" || $flag eq "N"){
					$rPos += $gapLen;
					#print "last $rPos >= $rLeft ----\n";
					last if($rPos >= $rLeft);
				}
				elsif($flag eq "I"){
					$qPos += $gapLen;
					#print "last $gapLen$flag, $prev{arr}[$ci], $alnQStart, $qPos >= $qLeft ++++\n";
					last if($qPos >= $qLeft);
				}
				elsif($flag eq "S" && $qPos == 1){
					$qPos = $gapLen + 1;
				}
				$prevCigar .= $gapLen . $flag;
			}

			##### for the current segment
			my ($qRight, $rRight) = ($prev{qEnd}, $prev{rEnd});
			$qRight = $qStart + 100 if($prev{qEnd} - $qStart > 100);
			$rRight = $rStart + 100 if($prev{rEnd} - $rStart > 100);

			my $bGetPos = 0;
			$qPos = 1;
			$rPos = $rStart;
			$postCigar = "";
			while($aln->[$ci] =~ /(\d+)(\w)/g){
				$gapLen = $1; $flag = $2;
				if($flag eq "M"){		
					if($bGetPos == 0){ 
						#print "if($rRight < $rPos && $qRight < $qPos){\n";
						if($rRight < $rPos && $qRight < $qPos){
							$alnREnd = $rPos - 1;
							$alnQEnd = $qPos - 1;  #####################
							#print " alnQEnd $alnQEnd\n";
							$bGetPos = 1;
						}
						elsif($rRight < $rPos + $gapLen && $qRight < $qPos + $gapLen ){
							if($rRight - $rPos > $qRight - $qPos){
								$diff = ($rRight - $rPos + 1);
							}
							else{
								$diff = ($qRight - $qPos + 1);
							}
							$postCigar .= ($gapLen-$diff -1) . $flag;	######################
							$alnREnd = $rPos + $diff;
							$alnQEnd = $qPos + $diff;  ######################
							$bGetPos = 1;
							next;
						}
					}
					$rPos += $gapLen;
					$qPos += $gapLen;

				}
				elsif($flag eq "D" || $flag eq "N"){
					$rPos += $gapLen;
				}
				elsif($flag eq "I"){
					$qPos += $gapLen;
				}
				elsif($flag eq "S" && $qPos == 1){
					$qPos = $gapLen + 1;
				}
				$postCigar .= $gapLen . $flag if($bGetPos); #if($prev{rEnd} < $rPos && $prev{qEnd} < $qPos);				
			}
		}
		#print "0000 1 $alnQEnd"."M$postCigar\n";
		
		print "(alnQStart $alnQStart, alnQEnd $alnQEnd, alnRStart $alnRStart, alnREnd $alnREnd)\n"  if($verboseFlag);
		#print "$prevCigar, ".($alnQEnd-$alnQStart+1)."M, $postCigar\n";
#		my $newCigar = ($alnQEnd-$alnQStart+1)."M";
#=head
#		exit;
		my $newCigar = "";
		my ($lenA, $lenB, $res, $seqA, $seqB);
		
		

		$seqA = substr($refSeq, $alnRStart-1, $alnREnd-$alnRStart+1);
		$seqB = substr($readAll->{seq}, $alnQStart-1, $alnQEnd-$alnQStart+1);

		$res = SW_align($seqA, $seqB);
		$lenA = length($seqA);
		$lenB = length($seqB);
		
		#print " newCigar before : $newCigar,  variantions count : $#{$res->{variations}}+1\n";
		#print "$res->{seqA}\n$res->{seqB}\n";
		
		my $var;
		my $prevEnd = 0;
		$i = 0;
		$qPos = 0;
		
		if($#{$res->{variations}} > -1){			
			for(; $i <= $#{$res->{variations}} ; $i++) {
				$var = $res->{variations}[$i];
				
				if($var->{type} ne "M"){
					$newCigar .= ($var->{aPos}-$prevEnd - 1) . "M" if($var->{aPos} != 1);
					$newCigar .= $var->{len} . $var->{type};
					
					$prevEnd = $var->{aPos} - 1 + ($var->{type} eq "D" ? $var->{len} : 0);	
					$qPos = $var->{pos} - 1 + ($var->{type} eq "I" ? $var->{len} : 0);	
				}
			}
		}
		
		$newCigar .= ($lenB - $qPos) . "M" if($lenB > $qPos);		
		#print "newCigar $newCigar\n";

#=cut

		#print "----- rStart $rStart\n";
		$aln->[9] = $readAln{$qName}{seq};
		$aln->[10] = $readAln{$qName}{qual};

		my @arrCigar;
		cigar2Array($prevCigar.$newCigar.$postCigar, \@arrCigar);
		$aln->[$ci] = array2Cigar(\@arrCigar, $prev{rStart}, \$refSeq, $prev{qStart}, \$readAln{$qName}{seq});
		$aln->[3] = $prev{rStart};

		$prev{strand} = $strand;
		$prev{arr} = $aln;
		#$prev{qStart} = $qStart;
		$prev{qEnd} = $qEnd; 
		#$prev{rStart} = $rStart; 
		$prev{rEnd} = $rEnd;		

		#last if($ii++>1);
	}

	my @arrCigar;
	cigar2Array($prev{arr}[$ci], \@arrCigar);
	$prev{arr}[$ci] = array2Cigar(\@arrCigar, $prev{rStart}, \$refSeq, 1, \$readAln{$qName}{seq});
	print $out join("\t", @{$prev{arr}}) . "\n";
}
###########################



sub cigar2Array
{
	my ($cigar, $arr) = @_;
	
	while($cigar =~ /(\d+)([A-Z])/g){	
		$arr->[$#$arr+1]{len} = $1;
		$arr->[$#$arr]{type} = $2;
	}
}


sub array2Cigar
{
	my ($arr, $tPos, $refSeq, $qPos, $qSeq) = @_;

	$tPos -= 1; ### 1 based position to 0 based position
	$qPos -= 1;

	my $num = $#$arr + 1;
	my $cigar = "";
	my @cArr;
	my ($c, $diff, $i) = (-1, 0);	

	for($i = 0; $i < $num; $i++){
		if($i == $num - 1){
			$c++;
			$cArr[$c]{len} = $arr->[$i]{len};
			$cArr[$c]{type} = $arr->[$i]{type};
			last;
		}

		if($arr->[$i]{type} eq $arr->[$i+1]{type}){
			$arr->[$i+1]{len} += $arr->[$i]{len};
		}
		elsif($arr->[$i]{type} eq "M" || $arr->[$i+1]{type} eq "M")
		{
			$c++;
			$cArr[$c]{len} = $arr->[$i]{len};
			$cArr[$c]{type} = $arr->[$i]{type};
			$qPos += $arr->[$i]{len} if($arr->[$i]{type} ne "D");
			$tPos += $arr->[$i]{len} if($arr->[$i]{type} ne "I");
			#print " ** array2Cigar, tPos $tPos, qPos $qPos\n";
		}
		else #### $arr->[$i]{type} is I or D, and $arr->[$i+1]{type} D or I
		{
			if($arr->[$i]{len} == $arr->[$i+1]{len} && $arr->[$i]{type} ne $arr->[$i+1]{type}){
				if($c >= 0) { $cArr[$c]{len} += $arr->[$i]{len}; }			### $cArr[$c]{type} is 'M'. 
				elsif($arr->[$i+2]{type} eq "M") { $arr->[$i+2]{len} += $arr->[$i]{len}; }
				else{
					$c++;
					$cArr[$c]{len} = $arr->[$i]{len};
					$cArr[$c]{type} = "M";					
				}
				$i++;
			}
			elsif($arr->[$i]{len} < $arr->[$i+1]{len}){
				$diff = $arr->[$i+1]{len} - $arr->[$i]{len};				
				if(defined $arr->[$i+2] && $arr->[$i+2]{type} eq "M" &&
					(
					($arr->[$i]{type} eq "D" && substr($$refSeq, $tPos, $arr->[$i]{len}) eq  substr($$qSeq, $qPos+$diff, $arr->[$i]{len})) ||
					($arr->[$i]{type} eq "I" && substr($$qSeq, $qPos, $arr->[$i]{len}) eq  substr($$refSeq, $tPos+$diff, $arr->[$i]{len})) ## if two sequence are aligned at the right end.. 
					)
				)
				{
					$arr->[$i+2]{len} += $arr->[$i]{len};
				}
				elsif($c == -1 || $cArr[$c]{type} ne "M"){
					$c++;
					$cArr[$c]{len} = $arr->[$i]{len};
					$cArr[$c]{type} = "M";
					#print " ***** strainge at line " . __LINE__ ." of $0\n"; 
					#return undef;
				}
				else{
					$cArr[$c]{len} += $arr->[$i]{len};
				}
				
				$arr->[$i+1]{len} = $diff;
			}
			else ##### if($arr->[$i]{len} > $arr->[$i+1]{len})
			{	
				$diff = $arr->[$i]{len} - $arr->[$i+1]{len};				
				if(defined $arr->[$i+2] && $arr->[$i+2]{type} eq "M" &&
					(
					($arr->[$i]{type} eq "D" && substr($$qSeq, $qPos, $arr->[$i+1]{len}) eq  substr($$refSeq, $tPos+$diff, $arr->[$i+1]{len})) ||
					($arr->[$i]{type} eq "I" && substr($$refSeq, $tPos, $arr->[$i+1]{len}) eq  substr($$qSeq, $qPos+$diff, $arr->[$i+1]{len})) ## if two sequence are aligned at the right end.. 
					)
				)
				{
					$arr->[$i+2]{len} += $arr->[$i+1]{len};
				}
				elsif($c == -1 || $cArr[$c]{type} ne "M"){
					$c++;
					$cArr[$c]{len} = $arr->[$i+1]{len};
					$cArr[$c]{type} = "M";
					#print " ***** strainge at line " . __LINE__ ." of $0\n"; 					
				}
				else{
					$cArr[$c]{len} += $arr->[$i+1]{len};
				}

				$arr->[$i+1]{len} = $diff;
				$arr->[$i+1]{type} = $arr->[$i]{type};
			}
		}
	}

	foreach my $ci (@cArr) {
		$cigar .= $ci->{len} . $ci->{type};
	}
	return $cigar;
}



sub readRefSeq{
	my $refName = shift;

	my ($name);
	print STDERR "reading '$refName' sequence...\n";
	$rBuf = <$rin> if(!$rBuf); 
	while(defined $rBuf){
		if($rBuf =~ /^>/){			
			$name = $';
			$name =~ s/^\s*//g; $name =~ s/\s.*//g;
			if($name eq $refName){
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
		print STDERR "[Error] Could not read sequences for '$refName'. The order of references in '$refFn' and '$samFn' may be differ.\n";
		exit 1;
	}
}

sub min
{
	my ($i, $min);
	$min = $_[0];
	for($i = 1; $i <= $#_; $i++){
		$min = $_[$i]  if($min > $_[$i]);
	}
}


sub openInput
{
	my ($fn) = @_;

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



sub SW_align{
	my ($seqA, $seqB, $alignFlag) = @_;	  ### $alignFlag : align to right or left.. adjust gap open penalty...
	
	$seqA = uc $seqA;
	$seqB = uc $seqB;
	
	my ($lenA, $lenB);
	$lenA = length($seqA);
	$lenB = length($seqB);
	
	## the first column and the first row of matrix will have 0  
	## and next column and row will start with the first result..
	## so ...the result of position $A and $B will always be stored in $matrix[$B+1][$A+1]..	 
	my (@matrix);
	  
	## it will store the information for the directions..
	## 0 from diagonal, 1 from above, 2 from left 
	my (@moveFrom);
	
	

	## @sa, @sb : arraies of sequences.
	## $A, $B : positions of sa and sb
	## $i, $j : matrix rows..
	## $val : temporary value
	my(@sa, @sb, $A, $B, $i, $j, $val, $matchScore, $misPenalty, $gapPenalty, $gapExtPenalty);
	@sa = split(//, $seqA);
	@sb = split(//, $seqB);
	
	$matchScore = 3;
	$misPenalty = -6;
	$gapPenalty = -11;
	$gapExtPenalty = -1; #-4;
	
	$matrix[0][0][0] = 0;
	$matrix[0][0][1] = 0;
	$matrix[0][0][2] = 0;

	for($B = 1; $B <= $lenB; $B++){ 
		$matrix[$B][0][1] = $gapPenalty + ($gapExtPenalty+1)*($B-1);  ### from above. 
		$matrix[$B][0][0] = -1000000;
		$matrix[$B][0][2] = -1000000;
		$moveFrom[$B][0][1]= 1;  ### from above. 
	}
	for($A = 1; $A <= $lenA; $A++){	
		$matrix[0][$A][2] = $gapPenalty + ($gapExtPenalty+1)*($A-1);
		
		$matrix[0][$A][0] = -1000000;
		$matrix[0][$A][1] = -1000000;
		$moveFrom[0][$A][2]= 2; 
	}
	
	my @path = qw("\" "^" "<"); #### used to show the path..
	
	my ($max, $maxI, $penalty);

	
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
				$val = $matrix[$B][$A][$i] + $penalty; 				
				if($max < $val) {
					$max = $val;
					$maxI = $i;
				}
			}		
			$matrix[$B+1][$A+1][0] = $max;
			$moveFrom[$B][$A][0] = $maxI;
			
			
			## value from above direction	
			$max = -1000000;
			for($i = 0; $i < 3; $i++){
				$val = $matrix[$B][$A+1][$i] + ($B>1 && $i==1?($A==$lenA-1?$gapExtPenalty+1:$gapExtPenalty):$gapPenalty);
				$val -= 1 if($i == 2);
				if($max < $val){
					$max = $val;
					$maxI = $i;
				}
			}		
			$matrix[$B+1][$A+1][1] = $max;
			$moveFrom[$B][$A][1] = $maxI;
			
						
			## value from left direction  #### applying ends free alignment for B			
			$max = -1000000;
			for($i = 0; $i < 3; $i++){				
				$val = $matrix[$B+1][$A][$i] + ($A>1 && $i==2?($B==$lenB-1?$gapExtPenalty+1:$gapExtPenalty):$gapPenalty);
				$val -= 1 if($i == 1);
				if($max < $val){				
					$max = $val;
					$maxI = $i;
				}
			}		
			$matrix[$B+1][$A+1][2] = $max;
			$moveFrom[$B][$A][2] = $maxI;
			
			#print "$matrix[$B+1][$A+1]($path[$moveFrom[$B][$A]]\t";# $errMatrix[$B+1][$A+1])\t";
		}
		#print "\n";
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
			
			$gapNum++;
			$insNum++;
			
			$maxI = ($A == -1 ? 1 : $moveFrom[$B][$A][$maxI]);
			$B--; 
		}
		else { #from left direction
			$ssa[$i] = $sa[$A];
			$ssb[$j] = "-";
			#print "$A, $B, $ssa[$i], $ssb[$j], $moveFrom[$B][$A], $matrix[$B+1][$A+1]\n"; 
			
			$gapNum++;
			$delNum++;
			$delBases = $sa[$A] . $delBases;
			
			#print "\$maxI = \$moveFrom[$B][$A][$maxI]; $maxI = $moveFrom[$B][$A][$maxI];\n";
			#exit if(!defined $maxI);
			$maxI = ($B == -1 ? 2 : $moveFrom[$B][$A][$maxI]);
			$A--; 
		}
	}

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



