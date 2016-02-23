#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# readCountsAtLoci.pl
#   It shows read counts supporting a reference allele at each locous.
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
my ($oriFn, $newFn, $genomeSamFn, $samFn, $outFn);

GetOptions(
	"h|?|help"	=> \$helpFlag,

	"ori=s"	=> \$oriFn,
	"new=s"	=> \$newFn,
	"o|out=s"	=> \$outFn,
	"genome=s"	=> \$genomeSamFn,
	"sam=s"	=> \$samFn,
) || help(1);


help(0) if defined $helpFlag;


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 -ori <sequence file of an original genome seq.>   -g <SAM|BAM file containing alignments between an original and a new genome seq.>  -s <SAM|BAM file containing reads mapped to a new genome>  -o <output>\n\n";
	print STDERR "  ex) $0 -ori ori.seq  -g final.genome.sam  -s final.genome.reads.sam  -t final.align.var -o out.txt\n\n";
	exit($return);
}
# -new <sequence file of a new genome seq.> ,  -new final.genome.fa 


if(!defined $oriFn){ print STDERR "\nNeed a sequence file of the original genome sequence!!\n\n"; help(1); }
if(! -e $oriFn){ print STDERR "\n[Error] Could not find '$oriFn'.\n"; exit(1); }

#if(!defined $newFn){ print STDERR "\nNeed a sequence file of the original genome sequence!!\n\n"; help(1); }
#if(! -e $newFn){ print STDERR "\n[Error] Could not find '$newFn'.\n"; exit(1); }

if(!defined $samFn){ print STDERR "\nNeed a SAM|BAM file!!\n\n"; help(1); }
if(! -e $samFn){ print STDERR "\n[Error] Could not find '$samFn'.\n"; exit(1); }

if(!defined $genomeSamFn){ print STDERR "\nNeed a target loci file!!\n\n"; help(1); }
if(! -e $genomeSamFn){ print STDERR "\n[Error] Could not find '$genomeSamFn'.\n"; exit(1); }


######################## print starting time
$cmd = "$0 @ARGV";  ### $0 : command..  store all inputs(@ARGV) as a string in $cmd..
my @Month = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec");
my ($sec,$min,$hour,$mday,$mon,$year, $wday,$yday,$isdst) = localtime time;
print "### $cmd\n[Start] at $Month[$mon] $mday $hour:$min:$sec\n";
######################## print starting time



my ($i, %list_loci, @arr, $refName, $newRef, $str, %refIDs, @arrRef, $l, 
	$tStart, $tEnd, $tPos, $nPos);

my ($refFn, $rBuf, $inR, $refSeq, $refLen,);
$refFn = $oriFn;
$refName = "";
$newRef = "";

my $in = openInput($genomeSamFn);
print "$genomeSamFn\n";
while(<$in>){
	next if(/^@/);
	chomp;
	@arr = split /\t/;

	my $loci;
	if($refName ne $arr[2]){
		$refName = $arr[2];
		readRefSeq($refName);
	}
	if($newRef ne $arr[0]){
		$newRef = $arr[0];
	}
	
	my $arrLoci = getMismatches($arr[3], $arr[5], \$arr[9]);
	push(@{$list_loci{$newRef}}, @$arrLoci);
}

=head
my @loci = sort {$a->{start} <=> $b->{start}} @{$list_loci{$newRef}};
my ($cntL, $l) = ($#loci + 1, 0);
while($l < $cntL){
	print "$loci[$l]{tStart}\t$loci[$l]{ori}\t$loci[$l]{start}\t$loci[$l]{new}\n";
	$l++;
}
exit;
=cut

sub getMismatches {
	my ($refPos, $cigar, $qSeq) = @_;
	
	my ($r, $q) = ($refPos-1, 0);
	my ($rC, $qC);

	my (@loci, $l);

	
	if($cigar =~ /^(\d+)S/){
		$q = $1;
		$cigar =~ s/^(\d+)S//;
	}

	$l = 0;
	while($cigar =~ /(\d+)([IMD])/g){
		my ($len, $type) = ($1, $2);		
		
		if($type eq "M"){
			for($i = 0; $i < $len; $i++, $r++, $q++){
				if($r >= $refLen) { print STDERR "$_\nref len : $refLen <= $r\n"; exit(1);}
				$rC = substr($refSeq, $r, 1);
				$qC = substr($$qSeq, $q, 1);
				if(!$rC || !$qC) { print STDERR "$_\ncigar : $cigar, $refName, ref len : $refLen, r($r, $rC), qLen : ".length($$qSeq).", q($q, $qC)\n"; exit(1);}
				if($rC ne $qC){
					$loci[$l]{start} = $q+1;
					$loci[$l]{end} = $q+1;

					$loci[$l]{tStart} = $r+1;
					$loci[$l]{tEnd} = $r+1;
					
					$loci[$l]{type} = "M";
					$loci[$l]{ori} = $rC;
					$loci[$l]{new} = $qC;
					$loci[$l]{sup} = 0;
					$loci[$l]{tot} = 0;
					$l++;
				}
			}
		}
		elsif($type eq "I"){
			$loci[$l]{type} = "I";
			$loci[$l]{ori} = "*";
			$loci[$l]{new} = "+" . substr($$qSeq, $q, $len);
			$loci[$l]{sup} = 0;
			$loci[$l]{tot} = 0;

			$loci[$l]{start} = $q+1;
			$q += $len;
			$loci[$l]{end} = $q;
			$loci[$l]{tStart} = $r+1;			
			$loci[$l]{tEnd} = $r+1;
			$l++;			
		}
		elsif($type eq "D"){			
			$loci[$l]{type} = "D";
			$loci[$l]{ori} = "*";
			$loci[$l]{new} = "-" . substr($refSeq, $r, $len);
			$loci[$l]{sup} = 0;
			$loci[$l]{tot} = 0;


			$loci[$l]{start} = $q+1;
			$loci[$l]{end} = $q+1;
			$loci[$l]{tStart} = $r+1;
			$r += $len;
			$loci[$l]{tEnd} = $r;
			$l++;
		}
	}

	return \@loci;
}

$in = openInput($samFn);
my $out = openOutput($outFn);

$refName = "";  ### it will contains a name of the new genome sequence
my (@loci, $cntL, );
($cntL, $l) = (0,0);

while(<$in>){
	next if(/^@/);
	chomp;
	@arr = split /\t/;

	last if($arr[2] eq "*");

	if($refName ne $arr[2]){
		#print " +++ $arr[2]\n";
		while($l < $cntL){ 
			print $out "$refName\t$loci[$l]{start}\t$loci[$l]{ori}\t$loci[$l]{new}\t$loci[$l]{sup}/$loci[$l]{tot}\n";
			$l++;
		}

=head
		if($refName ne "" && $refIDs{$refName} + 1 != $refIDs{$arr[2]}){
			for($i = $refIDs{$refName} + 1; $i < $refIDs{$arr[2]}; $i++){
				@loci = sort {$a->{start} <=> $b->{start}} @{$list_loci{$arrRef[$i]}};
				($cntL, $l, $p) = ($#loci + 1, 0);
				while($l < $cntL){ 
					print $out "$loci[$l]{ori}\t$loci[$l]{sup}/$loci[$l]{tot}\n";
					$l++;
				}
			}
		}
=cut

		if($list_loci{$arr[2]}){
			@loci = sort {$a->{start} <=> $b->{start}} @{$list_loci{$arr[2]}};
			($cntL, $l) = ($#loci + 1, 0);
		}
		else{
			($cntL, $l) = (0, 0);
		}		
		$refName = $arr[2];
	}

	
	$tStart = $tEnd = $arr[3];

	while($l < $cntL && $tStart > $loci[$l]{start}){ 
		print $out "$refName\t$loci[$l]{start}\t$loci[$l]{ori}\t$loci[$l]{new}\t$loci[$l]{sup}/$loci[$l]{tot}\n";
		$l++;
	} 
	next if($l >= $cntL);

	$tEnd--;
	my (@arrMatch, $len, $flag);
	while($arr[5] =~ /(\d+)([A-Z])/g){
		($len, $flag) = ($1, $2);
		if($flag eq "S" || $flag eq "H"){			
		}
		elsif($flag eq "I"){
		}
		elsif($flag eq "D"){
			$tEnd += $len;
		}
		elsif($flag eq "M"){
			$i = $#arrMatch+1;
			$arrMatch[$i]{start} = $tEnd + 1;
			$tEnd += $len;
			$arrMatch[$i]{end} = $tEnd;
		}
	}

	my $tmpL = $l;
	while($tmpL < $cntL && $tEnd > $loci[$tmpL]{start}){
		if($tEnd < $loci[$tmpL]{end}){ $tmpL++; next; } 
		$loci[$tmpL]{tot}++;


		my ($md, $mdi);
		for($i = 11; $i <= $#arr; $i++)
		{
			if($arr[$i] =~ /^MD:Z:/){
				$md = $';
				$mdi = $i;
				last
			}
		}

		my (@arrMismatch, $str, $j);
		$md =~ /^(\d+)/;
		my $pos = $arr[3] + $1;#+1;
		$md =~ s/^\d+//;
		while($md =~ /([A-Z\^]+)(\d+|)/g) {
			($str, $len) = ($1, $2);
			#print "$pos : $1\n";
			
			if($1 =~ /\^/){
				$pos += length($str)-1 + $len;
			}
			else{
				$arrMismatch[$#arrMismatch+1] = $pos;#+1;
				$pos += 1 + $len;
			}
		}

		if($loci[$tmpL]{type} eq "I" || $loci[$tmpL]{type} eq "D"){
			for($i = 0; $i <= $#arrMatch; $i++){
				if($arrMatch[$i]->{start} <= $loci[$tmpL]{start} && $loci[$tmpL]{end} <= $arrMatch[$i]->{end}){
					$loci[$tmpL]{sup}++;
					last;
				}
				last if($arrMatch[$i]->{start} >= $loci[$tmpL]{start});
			}
		}
		else{ ### mismatch
			for($i = 0; $i <= $#arrMatch; $i++){
				if($arrMatch[$i]->{start} <= $loci[$tmpL]{start} && $loci[$tmpL]{end} <= $arrMatch[$i]->{end}){
					my $bMismatch = 0;
					for($j = 0; $j <= $#arrMismatch; $j++){
						#if($arr[0]=~/:1111:12462:9874/){ print "mismatch : $arrMismatch[$j], loci : $loci[$tmpL]{start} \n"; }
						if($arrMismatch[$j] == $loci[$tmpL]{start}) {
							$bMismatch = 1;
						}
						last if($arrMismatch[$j] >= $loci[$tmpL]{start});
					}
					$loci[$tmpL]{sup}++ if($bMismatch == 0);
					last;
				}
				last if($arrMatch[$i]->{start} >= $loci[$tmpL]{start});
			}
		}
		$tmpL++;
	}

}

if($refName ne ""){
	while($l < $cntL){ 
		print $out "$refName\t$loci[$l]{start}\t$loci[$l]{ori}\t$loci[$l]{new}\t$loci[$l]{sup}/$loci[$l]{tot}\n";
		$l++;
	}
}

close($in);
close($out);





sub readRefSeq{
	my $rName = shift;

	if($rName eq ""){
		print STDERR "Wrong reference name.. '$rName'.\n$_\n";
	}

	if(!defined $inR){
		$inR = openInput($refFn);
	}
	
	my ($name);
	print STDERR "reading '$rName' sequence...\n";
	$rBuf = <$inR> if(!$rBuf); 
	while(defined $rBuf){
		if($rBuf =~ /^>/){			
			$name = $';
			$name =~ s/^\s*//g; $name =~ s/\s.*//g;
			if($name eq $rName){
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
		print STDERR "Could not read sequences for '$rName'. The order of references in '$refFn' and '$genomeSamFn' may be different.\n";
	}
	$refLen = length($refSeq);
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

