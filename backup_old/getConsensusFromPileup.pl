#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2010
#
# getCoverageAtTargetInPileup.pl
# Copyright   : GPL V3 (http://www.gnu.org/licenses/gpl-3.0.html)
######################################

use strict;
use Getopt::Long qw(:config no_ignore_case);
use IO::Handle;
STDOUT->autoflush(1);
STDERR->autoflush(1);

my ($cmd, $helpFlag, $noGapFlag, $header, $pileupFn, $samFn, $conPosFlag, $outFn);
$cmd = "$0 @ARGV";
my $minRate_varPrint = 0.2;  ## min. rate of reads to print the variants supported by the reads
my $minRate_indelReplace = 0.4;  ## min. rate of reads to replace reference bases with INDELs supported by the reads
### INDELs will be printed if they are supported by reads in MIN($minRate_varPrint, $minRate_indelReplace) rate
my $minDepth = 1;

$header = "align";

GetOptions(
	"h|?|help"		=> \$helpFlag,
	"nogap"		=> \$noGapFlag,
	"pos"		=> \$conPosFlag,
	#"contig"		=> \$contigFlag,
	"var=f"		=> \$minRate_varPrint,
	"indel=f"		=> \$minRate_indelReplace, 
	"depth=i"		=> \$minDepth,
	"sam=s"	=> \$samFn,
	"head|header=s"	=> \$header,
) || help(1);

help(0) if defined $helpFlag;

$pileupFn = shift;
if(defined $pileupFn && ($pileupFn =~ /sam$/ || $pileupFn =~ /bam$/)){
	$samFn = $pileupFn;
	$pileupFn = shift;
}

$samFn = shift if(!defined $samFn);
if(!defined $pileupFn){ print STDERR "\nNeed a pileup file!!\n\n"; help(1); }
if(!defined $samFn){ print STDERR "\nNeed a SAM file!!\n\n"; help(1); }

#if(defined $outFn && ($pileupFn eq $outFn || $samFn eq $outFn))
#{ print STDERR " Error) Input and output files are same \n"; exit 1; }

sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 <pileup file>  <sam file> \n"
		."\t[-header <string to be used for sequence names. default : '$header'>\n"
		."\t[-var <min. rate of reads to print variants supported by the reads, default : $minRate_varPrint>]\n"
		."\t[-d <min. depth to print variants, default : $minDepth>]\n"
		."\t[-indel <min. rate of reads to replace references with INDELs supported by the reads, default : $minRate_indelReplace>]\n"
		."\t[-pos (a flag to print positions of variants on the consensus)\n"
		."\t[-nogap (a flag to concatenate all blocks on a chromosome)\n"
		#."\t[-contig (a flag to use contigs. A consensus includes any INDELs supported by reads in rates higher than 0.2)\n"
		;
	print STDERR " \n ex) $0 contigs.merged.pileup contigs.merged.sam -var 0.3 -pos\n\n";
	exit($return);
}

my ($in, $out, @arr, $i);

my (%base2ID, @id2Base);
my $str = "BDEFHIJKLMNOPQRSUVWXYZ";
for($i = 0; $i < length($str); $i++){
	$base2ID{substr($str, $i, 1)} = 4;
}

$base2ID{'a'} = 0;
$base2ID{'A'} = 0;
$base2ID{'t'} = 1;
$base2ID{'T'} = 1;
$base2ID{'g'} = 2;
$base2ID{'G'} = 2;
$base2ID{'c'} = 3;
$base2ID{'C'} = 3;
$base2ID{'N'} = 4;
$id2Base[0] = 'A';
$id2Base[1] = 'T';
$id2Base[2] = 'G';
$id2Base[3] = 'C';
$id2Base[4] = 'N';


my (%refLen, $ref, $len, $start, $end, $readLen);
if(! -e $samFn){ print STDERR "Could not find $samFn.\n\n"; exit; }

open($in, $samFn =~ /.bam$/ ? "samtools view -h $samFn |"  : $samFn);
$i = 0;
$readLen = 0;
while(<$in>){
	s/[\r\n]+//g;
	next if(/^#/ || /^\s*$/);
	last if($i > 1000);  ### check the lengths of the first 1000 reads.
	if(/^\@SQ\tSN:([\w_+\-:\.\|#]+)\tLN:(\d+)/){
		$ref = $1; 
		$len = $2;
		if($ref =~ /^(\d+)$/ || $ref eq "X" || $ref eq "M" || $ref eq "Y"){
			$ref = "chr" . $ref;
		}
		$refLen{$ref} = $len;		
	}	
	elsif(!/^@/){
		@arr=split /\t/;
		my $len = length($arr[9]);
		$readLen = $len if($readLen < $len); 
		$i++;
	}
}
close($in);

my($totCov, $cov, $totCovLen, $totLen, $coveredLen, $prePos, $blockNum, $rPos, $qPos);

$in = openInput($pileupFn);

$outFn = $pileupFn;

$outFn =~ s/.mpileup$/.align.seq/;
$outFn =~ s/.mpileup.gz$/.align.seq/;
$outFn =~ s/.pileup$/.align.seq/;
$outFn =~ s/.pileup.gz$/.align.seq/;
if($outFn eq $pileupFn) { $outFn = $pileupFn . ".align.seq"; }

print "Writing\n$outFn\n";
my $outSeq = openOutput($outFn);

$outFn =~ s/.seq$/.coverage/;
print "$outFn\n";
my $outCov = openOutput($outFn);

$outFn =~ s/.coverage$/.gap/;
print "$outFn\n";
my $outGap = openOutput($outFn);

$outFn =~ s/.gap$/.var/;
print "$outFn\n";
my $outV = openOutput($outFn);
print $outCov "#[Reference]\t[Ref. Length]\t[Average coverage]\t[Covered bases (percentage)]\t[Base level identity for covered]\t[Base level identity for whole]\n";
print $outV "#[Reference]\t[Pos in ref.]\t[Ref. base]\t" . (defined $conPosFlag ? "[Block]\t[Pos in block]\t" : "") . 
	"[Var.]\t[# of reads in Var./# of total reads]\n";

my $prevEndCnt = 0;
my $bRightPrinted = 0;
my ($varNum, $totVarNum) = (0,0);
$ref = undef;
$cov = 0;
$totCov = 0;
$totCovLen = 0;
$totLen = 0;
$coveredLen = 0;
$prePos = 0;
$blockNum = 1;
$qPos = 1;
my $isGapPrinted = 0;

my ($readCnt_all, $pileupType, $baseStr, $rBase, $conBase);

my ($lenPrintRight, $skipForSeq) = (1000, 0);
my ($conSeq) = ("");

my $bINDEL = 0;
while(<$in>){
	s/[\r\n]+//g;
	next if(/^#/ || /^\s*$/);
	
	####0:reference, 1:position in ref., 2:reference base 3:Consensus base, 4: ,5: 6: , 7:number of reads, 8: match or mismatch, 9:quality,
	@arr = split /\t/;
	
	if(!defined $pileupType){
		if($#arr == 5){
			$pileupType = "normal";
		}
		elsif($#arr == 9){
			$pileupType = "maq";
		}
		else{
			print STDERR "Unknown format : $pileupFn\n"; exit(1);
		}
	}

	if($arr[1] != $prePos && $skipForSeq > 0){
		#if($bINDEL) {print "skipForSeq $skipForSeq, rPos $arr[1], qPos : $qPos\n"; }

		$skipForSeq--;
		$prePos = $arr[1];
		next;
	}


	if($arr[0] =~ /^(\d+)$/ || $arr[0] eq "X" || $arr[0] eq "M" || $arr[0] eq "Y"){
		$arr[0] = "chr" . $arr[0];
	}

	if(!defined $refLen{$arr[0]}){		
		print STDERR "$samFn does not have length information for '$arr[0]'\n"; exit 1;
	}
	
	if(!defined $ref || $ref ne $arr[0]){
		printBlock() if($prePos != 0);
		printCoverage() if(defined $ref);		
		printTailGap();
		print $outGap "\n" if(defined $ref);		
		print $outGap ">$arr[0]\n";
		$cov = 0;
		$ref = $arr[0];
		$coveredLen = 0;
		$prePos = 0;
		$conSeq = "";
		$bRightPrinted = 0;
		$varNum = 0;
		$isGapPrinted = 0;
	}

	
	$rPos = $arr[1];
	$rBase = $arr[2];
	if($pileupType eq "normal"){
		$readCnt_all = $arr[3];
		$baseStr = $arr[4];
	}
	else{
		$readCnt_all = $arr[7];
		$baseStr = $arr[8];
	}

	$cov += $readCnt_all;
	

	if(!$noGapFlag && $rPos - 1 > $prePos){	#### if there is a gap...
		printBlock() if($prePos != 0);

		print $outGap "- gap start : " . ($prePos+1). ", end : " . ($rPos-1). "\n";
		my $len = length($conSeq);
		if($len < 100) { $i = 0;}
		else { $i = $len - 100; $len = 100; }
		print $outGap "  left : " . substr($conSeq, $i, $len) . "\n"; 
		$conSeq = "";
		$lenPrintRight = 0;
		$isGapPrinted = 1;
	}
	if(!$noGapFlag && $isGapPrinted == 1 && $lenPrintRight == 101){
		print $outGap "  right : " . substr($conSeq, 0, 100) . "\n\n";
		$bRightPrinted = 1;
		$lenPrintRight++;
	}

	if($prePos != $rPos){ 
		$coveredLen++;
		
		my (%indelReadCnt, $indelStr, @indels);
		$prevEndCnt = 0;
		$i = 0;

		my $chCnt = length($baseStr);
		my $depth = $chCnt;
		my (@cnt, $c, $j);
		for($i = 0; $i < $chCnt; $i++) {
			$c = substr($baseStr, $i, 1);
			if($c eq "^") { $i++; next; }
			next if($c eq "^" || $c eq "~" || $c eq "*" || $c eq "!" || $c eq "F" || $c eq "]" );
			if($c eq "\$") { $prevEndCnt++; next; }

			if($c eq "+" || $c eq "-"){
				$indelStr = $c;
				$i++;

				$str = substr($baseStr, $i);
				$str =~ /^(\d+)/;								
				$i += length($1) + $1 -1;
				$indelStr .= uc substr($str, length($1),$1);
				$depth -= (length($1)+$1);

				$indelReadCnt{$indelStr}{type} = $c;
				$indelReadCnt{$indelStr}{cnt}++;
				
				next;
			}
				
			if($c eq "." || $c eq ","){ $cnt[$base2ID{$rBase}]++; }
			elsif(defined $base2ID{$c}){
				$cnt[$base2ID{$c}]++;
			}
		}

		my ($max, $mID, $seMax, $seMID) = (0, -1, 0, -1);
		for($i = 0; $i <= $#cnt; $i++){
			if(defined $cnt[$i] && $max <= $cnt[$i] ){
				$seMax = $max; $seMID = $mID; $max = $cnt[$i]; $mID = $i;
			}				
		}		

		if($mID != -1 && $seMID != -1 && $max == $seMax && $id2Base[$seMID] eq $rBase){
			$i = $mID; $mID = $seMID; $seMID = $i;
		}
		$conBase = $id2Base[$mID];

		if($mID != -1 && $id2Base[$mID] ne 'N'){
			#$arr[3] = $id2Base[$mID] if($max > $readCnt_all*0.5 || $arr[3] !~ /[ATGC]/);
			
			
			if(!defined $base2ID{$rBase} || !defined $mID  || !defined $max ){
				print "Error : (\$base2ID{$rBase}, $base2ID{$rBase} != $mID  && $max > $readCnt_all*$minRate_varPrint )\n"; exit 1;
			}

			
			if($depth >= $minDepth){
				if($base2ID{$rBase} != $mID  && $max > $readCnt_all* $minRate_varPrint){
					if(defined $conPosFlag) { print $outV "$ref\t$rPos\t$rBase\t$header.block_$blockNum\t$qPos\t$conBase\t$max/$readCnt_all\n"; }
					else { print $outV "$ref\t$rPos\t$rBase\t$conBase\t$max/$readCnt_all\n"; }

					$varNum++;			
				}
				elsif($base2ID{$rBase} != $seMID  && 
					$seMID > -1 && $base2ID{$rBase} != $seMID  && $seMax > $readCnt_all*$minRate_varPrint )
				{
					if(defined $conPosFlag) { print $outV "$ref\t$rPos\t$rBase\t$header.block_$blockNum\t$qPos\t$id2Base[$seMID]\t$seMax/$readCnt_all\n"; }
					else { print $outV "$ref\t$rPos\t$rBase\t$id2Base[$seMID]\t$seMax/$readCnt_all\n"; }
					$varNum++;
				}
			}
		}
		
		$conSeq .= $conBase;
		if( !defined $conBase ) {print "$_\n"; exit 1;}
		if($lenPrintRight < 101){
			$lenPrintRight++;
		}
		$qPos++;

		################################## for INDELs
		if($rPos > 1 && $rPos < $refLen{$ref}){
			@indels = sort {$indelReadCnt{$b}{cnt} <=> $indelReadCnt{$a}{cnt}} keys %indelReadCnt;
			#if($rPos == 736306){
			#	print "---\n";
			#	foreach my $st (@indels) {
			#		print "$st, $indelReadCnt{$st}{cnt} \n";
			#	}
			#}

			### if there are several different length of INDELs, then take the longest of the first two INDELs.
			if($#indels > 0 && $indelReadCnt{$indels[0]}{type} eq $indelReadCnt{$indels[1]}{type} && 
					$indelReadCnt{$indels[0]}{cnt} < 2*$indelReadCnt{$indels[1]}{cnt}){ ### if # of reads for the second highest read indel is more than half of ..
				$indelReadCnt{$indels[0]}{cnt} += $indelReadCnt{$indels[1]}{cnt}; ### count as one INDELs
				if(length($indels[0]) <length($indels[1])) { ## take the longest..
					$indels[0] = $indels[1];
				}
			}
			if($#indels > -1){
				$readCnt_all -= $prevEndCnt;

				$indelStr = substr($indels[0], 1);
				my $indelLen = length($indelStr);				

				my ($indelType, $indelReadCnt) =  ($indelReadCnt{$indels[0]}{type}, $indelReadCnt{$indels[0]}{cnt});
				
				my $cutOffCnt = $readCnt_all * ($minRate_varPrint<$minRate_indelReplace ? $minRate_varPrint:$minRate_indelReplace) 
					#* (($readLen - $indelLen + 1)/$readLen);
					* ($indelType eq "+" ? ($readLen - $indelLen + 1)/$readLen : 1);
				
				#print "$rPos, depth : $depth,  readAll : $readCnt_all, $indelReadCnt >= $cutOffCnt, \n" if($rPos == 736269);
				$cutOffCnt = 1 if($cutOffCnt < 1);
				
				if($indelReadCnt >= $cutOffCnt)
				{				
					if($depth >= $minDepth){
						if(defined $conPosFlag) { print $outV "$ref\t" . ($rPos+1) . "\t*\t$header.block_$blockNum\t$qPos\t$indels[0]\t$indelReadCnt/$readCnt_all\n"; }
						else { print $outV "$ref\t" . ($rPos+1) . "\t*\t$indels[0]\t$indelReadCnt/$readCnt_all\n"; }

						$varNum += $indelLen;
					}				
					
					if($indelReadCnt >= $readCnt_all * $minRate_indelReplace){ ### if it is a contigs, it takes any gaps supported by at least 20% of reads.
						if($indelType eq "+"){			
							$conSeq .= $indelStr;
							if($lenPrintRight < 100){
								$lenPrintRight += $indelLen;
							}
							$qPos += $indelLen;
						}
						elsif($indelType eq "-"){
							$skipForSeq = $indelLen;
						}
						
						$bINDEL = 1;			
					}
				}
			}
		}
		################################## for INDELs		
	}
	$prePos = $rPos;
}

printBlock() if($prePos != 0);
printTailGap();
printCoverage() if(defined $ref);

print $outCov "\tTotal length $totLen bases, average coverage : ". sprintf("%.2f", $totCov/$totLen). "\n";
print $outCov sprintf("\tTotal Base level identity for covered : %.4f %%, ", ($totCovLen-$totVarNum)*100/$totCovLen) . 
						sprintf("Total Base level identity for whole : %.4f %%\n", ($totCovLen-$totVarNum)*100/$totLen );

close($in);
if(defined $outFn){ 
	close($outSeq); close($outCov); close($outV); 
}

sub processINDELs{
	my ($pos, $seq, $type, $readCnt_indel, $readCnt_all);

	print $outV "$ref\t" . ($pos+1) . "\t*\t$header.block_$blockNum\t" . ($qPos+1) . "\t$type$seq\t$readCnt_indel/$readCnt_all\n";
				
	my $indelLen = length($seq);
	$varNum += $indelLen;
	if($readCnt_indel >= $readCnt_all/2){
		if($type eq "+"){			
			$conSeq .= $seq;
			if($lenPrintRight < 100){
				$lenPrintRight += $indelLen;
			}
			$qPos += $indelLen;
		}
		elsif($type eq "-"){
			$skipForSeq = $indelLen;
		}
	}
}

sub printTailGap{	
	if(defined $ref && $refLen{$ref} - $prePos > 1){		
		print $outGap "- gap start : " . ($prePos+1). ", end : $refLen{$ref}\n";
		my $len = length($conSeq);
		if($len < 100) { $i = 0; $len -= 1; }
		else { $i = $len - 101; $len = 100; }
		print $outGap "  left : " . substr($conSeq, $i, $len) . "\n"; $conSeq = "";
	}
}


sub printBlock{
	if($conSeq ne "" && $isGapPrinted == 1 && $bRightPrinted == 0){ 
		print $outGap "  right : " . substr($conSeq, 0, 100) . "\n\n";
	}
	
	print $outSeq ">$header.block_$blockNum\n";
	print $outGap "[$header.block_$blockNum]\n\n";
	printSeqByMultiLine($outSeq, \$conSeq, 80);
	print $outSeq "\n";
	$blockNum++;
	$qPos = 1;
}

sub printCoverage{
	print $outCov "$ref\t" .(defined $refLen{$ref} ? $refLen{$ref} : 0). "\t". (defined $refLen{$ref} ? sprintf("%.2f", $cov/$refLen{$ref}) : 0) . 
		" Covered $coveredLen (" . sprintf("%.3f", $coveredLen/$refLen{$ref} * 100) . "%)";

	print $outCov sprintf("\t%.4f%%\t%.4f%%\n", ($coveredLen-$varNum)*100/$coveredLen, ($coveredLen-$varNum)*100/$refLen{$ref});
	$totVarNum += $varNum;
	$totCov += $cov;
	$totCovLen += $coveredLen;
	$totLen += $refLen{$ref};

}


sub printSeqByMultiLine
{
	my ($out, $seq, $lineLen) = @_;
	my ($i, $len);
	$len = length($$seq);

	for($i = 0; $i < $len; $i += $lineLen){		
		print $out substr($$seq, $i, $lineLen) . "\n";
	}
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


