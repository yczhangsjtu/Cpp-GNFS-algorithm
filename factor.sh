#!/bin/bash

path=/tmp/
nfile=${path}n.txt
polyfile=${path}poly.txt
basefile=${path}factorbase.txt
pairfile1=${path}abpairs1.txt
pairfile2=${path}abpairs2.txt
factorfile=${path}factor.txt

sieve=latticesieve

ccred=$(echo -e "\033[0;31m")
ccgreen=$(echo -e "\033[0;32m")
ccend=$(echo -e "\033[0m")
info="[${ccgreen}Info${ccend}]: "
error="[${ccred}Error${ccend}]: "

if [[  -z "$1" || (! "$(grep "^[ [:digit:] ]*$" <<< $1)") ]]; then
	n=$(python -c 'print 2**50+1')
else
	n=$1
	shift
fi

while [ -n "$1" ]; do
	if [ "$1" == "-lattice" ]; then
		sieve=latticesieve
		shift
	fi
	if [ "$1" == "-linear" ]; then
		sieve=sieve
		shift
	fi
	if [ "$1" == "-b" ]; then
		abratio="-b $2"
		shift 2
	fi
	if [ "$1" == "-a" ]; then
		Afactor="-a $2"
		shift 2
	fi
done

echo "${info}Factoring $n..."
echo $n > $nfile

echo "${info}Selecting polynomial..."
./polyselect $nfile $polyfile
if [ $? != "0" ]; then
	echo "${error}polyselect failed!"
	exit 1
fi

echo "${info}Forming factor bases..."
./factorbase $polyfile $basefile
if [ $? != "0" ]; then
	echo "${error}factorbase failed!"
	exit 1
fi

echo "${info}Sieving..."
/usr/bin/time -f %e -o /tmp/timefile ./$sieve $basefile $pairfile1 $Afactor $abratio
if [ $? != "0" ]; then
	echo "${error}$sieve failed!"
	exit 1
fi
echo "${info}Time consumed by $sieve: $(head -n 1 /tmp/timefile)s"

echo "${info}Checking if there is redundancy in the sieved pairs..."
./checkredundancy $pairfile1
if [ $? != "0" ]; then
	echo "${error}checkredundancy failed!"
	exit 1
fi
echo "${info}No redundant."

echo "${info}Solving linear system..."
./linear $pairfile1 $pairfile2
if [ $? != "0" ]; then
	echo "${error}linear failed!"
	exit 1
fi

echo "${info}Sqrting..."
./sqrt $pairfile2 $factorfile
if [ $? != "0" ]; then
	echo "${error}sqrt failed!"
	exit 1
fi

echo "${info}Result: $(cat $factorfile)"
