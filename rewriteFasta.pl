#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# rewriteFasta.pl
# reformats a FASTA file.
# Copyright   : GPL V3 (http://www.gnu.org/licenses/gpl-3.0.html)
######################################

use strict;
use warnings "all";
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, $seqFn, $lenLine, $outFn);

$lenLine = 80;
GetOptions(
	"h|?|help"		=> \$helpFlag,
	"input=s"	=> \$seqFn,	
	"length=s"	=> \$lenLine,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;

sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 <fasta file> -l <length per line (default:$lenLine)> [-o <out file>]\n";
	print STDERR "  ex) $0 sequence.fa -l 80 \n";
	print STDERR "  ex) $0 sequence.fa -l 60 -o new.fa\n\n";
	exit($return);
}

$seqFn = shift if(!defined $seqFn);
if(!defined $seqFn){ print STDERR "\nNeed a sequence file!!\n\n"; help(1); }
$outFn = shift if(!defined $outFn);


my ($in, $out);

my ($name, $seq, $len, $i, $tmpFn);
$/ = ">"; 


$in = openInput($seqFn);
if($outFn){
	$out = openOutput($outFn);
}
else{
	$tmpFn = $seqFn."__";
	$out = openOutput($tmpFn);
}


while(<$in>)
{	
	next unless (($name,$seq) = /(.*?)\n(.*)/s);
	$seq =~ s/[\d\s>]//g; #remove digits, spaces, line breaks,...
	$name =~ s/^\s*//g; $name =~ s/\s.*//g;
	print "Processing $name.....\n";
	$len = length($seq);

	print $out ">$name\n";
	for($i = 0; $i < $len; $i+= $lenLine){
		print $out substr($seq, $i, $lenLine) . "\n";
	}
	
}

close($in);
close($out);
if(!$outFn){
	system("mv -f $tmpFn $seqFn");
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

