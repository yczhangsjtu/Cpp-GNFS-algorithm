flags=-lflint -lgmp

.PHONY: all
.DEFAULT: all
all: polyselect factorbase sieve latticesieve linear checkredundancy sqrt

polyselect: polyselect.cpp polyselect.h GNFS.h poly.o util.o mypair.h
	g++ $< poly.o util.o $(flags) -o $@

factorbase: factorbase.cpp factorbase.h GNFS.h poly.o util.o mypair.h
	g++ $< poly.o util.o $(flags) -o $@

sieve: sieve.cpp sieve.h GNFS.h rational.o algebraic.o util.o poly.o
	mpic++ $< poly.o util.o rational.o algebraic.o $(flags) -o $@

latticesieve: latticesieve.cpp latticesieve.h mypair.h GNFS.h poly.o util.o latticeutil.o
	mpic++ $< poly.o util.o latticeutil.o $(flags) -o $@

linear: linear.cpp linear.h GNFS.h poly.o util.o mypair.h
	g++ $< poly.o util.o $(flags) -o $@

checkredundancy: checkredundancy.cpp GNFS.h mypair.h
	g++ $< $(flags) -o $@

sqrt: sqrt.cpp sqrt.h GNFS.h poly.o util.o mypair.h
	g++ $< poly.o util.o $(flags) -o $@

algebraic.o: algebraic.cpp algebraic.h mypair.h GNFS.h
	g++ -c $< -o $@

poly.o: poly.cpp poly.h mypair.h GNFS.h
	g++ -c $< -o $@

rational.o: rational.cpp rational.h mypair.h
	g++ -c $< -o $@

util.o: util.cpp util.h mypair.h
	g++ -c $< -o $@

latticeutil.o: latticeutil.cpp latticeutil.h mypair.h
	g++ -c $< -o $@
