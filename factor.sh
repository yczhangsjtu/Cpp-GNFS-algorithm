#!/bin/bash

path=/tmp/
nfile=${path}n.txt
polyfile=${path}poly.txt
basefile=${path}factorbase.txt
pairfile1=${path}abpairs1.txt
pairfile2=${path}abpairs2.txt
factorfile=${path}factor.txt

ccred=$(echo -e "\033[0;31m")
ccgreen=$(echo -e "\033[0;32m")
ccend=$(echo -e "\033[0m")
info="[${ccgreen}Info${ccend}]: "
error="[${ccred}Error${ccend}]: "

if [ -z "$1" ]
then
	n=$(python -c 'print 2**50+1')
else
	n=$1
fi
echo "${info}Factoring $n..."
echo $n > $nfile

echo "${info}Selecting polynomial..."
./polyselect $nfile $polyfile
if [ $? != "0" ]
then
	echo "${error}polyselect failed!"
	exit 1
fi
echo "${info}Forming factor bases..."
./factorbase $polyfile $basefile
echo "${info}Sieving..."
./latticesieve $basefile $pairfile1
echo "${info}Solving linear system..."
./linear $pairfile1 $pairfile2
echo "${info}Sqrting..."
./sqrt $pairfile2 $factorfile
echo "${info}Result: $(cat $factorfile)"
