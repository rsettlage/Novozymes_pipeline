#!/bin/bash

if [ -z $2 ]
then
        echo $0 [reference file] [sam file]
        exit -1
fi

ref=$1
sam=$2

EXE_DIR=`dirname $0`

#SAMTOOLS="samtools" ## use this line if you use 'samtools' from in 'module'
SAMTOOLS="$EXE_DIR/samtools_0.1.16"



echo `date` start getProblematicReads.sh $sam
head=`echo $sam | sed 's/.sam$//'`

out=$head".long.sam"
perl -e 'while(<>){if(/^@/){print;next;} @arr=split /\t/; print if(abs($arr[8]) > 1000);}' $sam > $out
pile=$head".long.pileup"
cmd="$SAMTOOLS pileup -f $ref -S $out > $pile"
echo;echo $cmd
if !(eval $cmd); then echo "exit 1"; exit 1; fi;

wig=$head".long.wig"
cmd="$EXE_DIR/pileup2wig_wBin.pl -bin 10 -s $sam $pile -o $wig"
echo;echo $cmd
if !(eval $cmd); then echo "exit 1"; exit 1; fi;


out=$head".clipped.sam"
perl -e 'while(<>){if(/^@/){print;next;} @arr=split /\t/; print if($arr[5] =~ /(\d+)S/ && $1 > 3);}' $sam > $out
pile=$head".clipped.pileup"
cmd="$SAMTOOLS pileup -f $ref -S $out > $pile"
echo;echo $cmd
if !(eval $cmd); then echo "exit 1"; exit 1; fi;

wig=$head".clipped.wig"
cmd="$EXE_DIR/pileup2wig_wBin.pl -bin 10 -s $sam $pile -o $wig"
echo;echo $cmd
if !(eval $cmd); then echo "exit 1"; exit 1; fi;


out=$head".pairUnmapped.sam"
perl -e 'while(<>){if(/^@/){print;next;} @arr=split /\t/; print if($arr[1]&0x8);}' $sam > $out
pile=$head".pairUnmapped.pileup"
cmd="$SAMTOOLS pileup -f $ref -S $out > $pile"
echo;echo $cmd
if !(eval $cmd); then echo "exit 1"; exit 1; fi;

wig=$head".pairUnmapped.wig"
cmd="$EXE_DIR/pileup2wig_wBin.pl -bin 10 -s $sam $pile -o $wig"
echo;echo $cmd
if !(eval $cmd); then echo "exit 1"; exit 1; fi;



echo `date` end getProblematicReads.sh $sam
