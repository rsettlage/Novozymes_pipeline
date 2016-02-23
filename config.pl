#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2012
#
# config.pl
# Copyright   : GPL V3 (http://www.gnu.org/licenses/gpl-3.0.html)
######################################

use File::Basename;
my $EXE_DIR= dirname($0);


## change only this part ##########################################
$ENV{SAMTOOLS} = "$EXE_DIR/samtools_0.1.16"; #samtools should support 'pileup'

$ENV{BWA} = "$EXE_DIR/bwa_0.6.1";
#$ENV{BWA} = "bwa"; ## use this line if you use 'bwa' installed at the machine.
###################################################################


`$ENV{BWA} 2>&1`;
if($? == -1){
	print STDERR "\nCould not execute the '$ENV{BWA}' program.\n\n";  exit(1);
}

`$ENV{SAMTOOLS} 2>&1`;
if($? == -1){
	print STDERR "\nCould not execute the '$ENV{SAMTOOLS}' program.\n\n";  exit(1);
}

$ENV{REVISEQ_VER} = "0.1.1";

1;

