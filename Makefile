flags=-lflint -lgmp

.SECONDEXPANSION:

.PHONY: all
.DEFAULT: all
all: polyselect factorbase sieve latticesieve linear sqrt

polyselect: polyselect.cpp polyselect.h GNFS.h poly.o util.o
	g++ $< poly.o util.o $(flags) -o $@

factorbase: factorbase.cpp factorbase.h GNFS.h poly.o util.o
	g++ $< poly.o util.o $(flags) -o $@

sieve: sieve.cpp sieve.h GNFS.h rational.o algebraic.o
	g++ $< poly.o util.o rational.o algebraic.o $(flags) -o $@

latticesieve: latticesieve.cpp latticesieve.h mypair.h GNFS.h GNFS-lattice.h poly.o util.o latticeutil.o
	g++ $< poly.o util.o latticeutil.o $(flags) -o $@

linear: linear.cpp linear.h GNFS.h poly.o util.o
	g++ $< poly.o util.o $(flags) -o $@

sqrt: sqrt.cpp sqrt.h GNFS.h poly.o util.o
	g++ $< poly.o util.o $(flags) -o $@

algebraic.o: algebraic.cpp algebraic.h mypair.h GNFS.h
	g++ -c $< -o $@

factorbase.o: factorbase.cpp factorbase.h mypair.h GNFS.h
	g++ -c $< -o $@

linear.o: linear.cpp linear.h mypair.h GNFS.h
	g++ -c $< -o $@

poly.o: poly.cpp poly.h mypair.h GNFS.h
	g++ -c $< -o $@

polyselect.o: polyselect.cpp polyselect.h mypair.h GNFS.h
	g++ -c $< -o $@

rational.o: rational.cpp rational.h mypair.h GNFS.h
	g++ -c $< -o $@

sieve.o: sieve.cpp sieve.h mypair.h GNFS.h
	g++ -c $< -o $@

sqrt.o: sqrt.cpp sqrt.h mypair.h GNFS.h
	g++ -c $< -o $@

util.o: util.cpp GNFS.h
	g++ -c $< -o $@

latticeutil.o: latticeutil.cpp latticeutil.h mypair.h GNFS.h GNFS-lattice.h
	g++ -c $< -o $@
