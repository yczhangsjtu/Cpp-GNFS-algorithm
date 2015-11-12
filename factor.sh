#!/bin/bash

path=/tmp/
nfile=${path}n.txt
polyfile=${path}poly.txt
basefile=${path}factorbase.txt
pairfile1=${path}abpairs1.txt
pairfile2=${path}abpairs2.txt
factorfile=${path}factor.txt

sieve=./sieve

ccred=$(echo -e "\033[0;31m")
ccgreen=$(echo -e "\033[0;32m")
ccend=$(echo -e "\033[0m")
info="[${ccgreen}Info${ccend}]: "
error="[${ccred}Error${ccend}]: "

if [[  -z "$1" || (! "$(grep "^[ [:digit:] ]*$" <<< $1)") ]]; then
	n=$(python -c 'print 2**92+1')
else
	n=$1
	shift
fi

sievearg=""
np="-np 4"

while [ -n "$1" ]; do
	if [ "$1" == "-lattice" ]; then
		sieve="latticesieve"
		shift
	fi
	if [ "$1" == "-linear" ]; then
		sieve="sieve"
		shift
	fi
	if [ "$1" == "-b" ]; then
		sievearg="$sievearg -b $2"
		shift 2
	fi
	if [ "$1" == "-a" ]; then
		sievearg="$sievearg -a $2"
		shift 2
	fi
	if [ "$1" == "-s" ]; then
		sievearg="$sievearg -s $2"
		shift 2
	fi
	if [ "$1" == "-np" ]; then
		np="-np $2"
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
printf "${info}"
/usr/bin/time -f %e -o /tmp/timefile mpirun $np ./$sieve $basefile $pairfile1 $sievearg
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
/usr/bin/time -f %e -o /tmp/timefile ./linear $pairfile1 $pairfile2
if [ $? != "0" ]; then
	echo "${error}linear failed!"
	exit 1
fi
echo "${info}Time consumed by linear: $(head -n 1 /tmp/timefile)s"

echo "${info}Sqrting..."
./sqrt $pairfile2 $factorfile
if [ $? != "0" ]; then
	echo "${error}sqrt failed!"
	exit 1
fi

echo "${info}Result: $(cat $factorfile)"
