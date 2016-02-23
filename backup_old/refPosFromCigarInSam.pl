#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# template.pl
# This file will be used as a template for other perl scripts.
######################################

use strict;
use warnings "all";
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, $verboseFlag, $outFn);

if($#ARGV < 1){
	print STDERR "Usage : $0 <start pos> <cigar>\n\n";  exit 1;
}

my $i = $ARGV[0]-1;
my $len = 0;
my ($qStart, $qEnd) = (1, 0);
while($ARGV[1] =~ /(\d+)([SIDNME])/g){
	$i += $1 if($2 eq "M" || $2 eq "D");
	$len += $1 if($2 eq "M" || $2 eq "I" || $2 eq "S");
	$qEnd += $1 if($2 eq "M" || $2 eq "I");
	if($2 eq "S" && $qEnd < 2){ $qStart += $1; $qEnd += $1; }
}

print "length $len, qStart : $qStart, qEnd : $qEnd, rStart : $ARGV[0], rEnd : $i \n";

