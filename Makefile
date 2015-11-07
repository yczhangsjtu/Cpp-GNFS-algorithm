flags=-lflint -lgmp

listf=algebraic.cpp\
	factorbase.cpp\
	linear.cpp\
	NFS.cpp\
	poly.cpp\
	polyselect.cpp\
	rational.cpp\
	sieve.cpp\
	sqrt.cpp\
	util.cpp

listl=algebraic.cpp\
	factorbase.cpp\
	linear.cpp\
	poly.cpp\
	polyselect.cpp\
	rational.cpp\
	sieve.cpp\
	sqrt.cpp\
	util.cpp\
	latticeutil.cpp\
	latticesieve.cpp

.SECONDEXPANSION:

listof=$(subst cpp,o,$(listf))
listol=$(subst cpp,o,$(listl))

flint: GNFS-flint.cpp $(listof) GNFS.h
	g++ $< $(listof) $(flags) -o $@

lattice: GNFS-lattice.cpp $(listol) GNFS.h GNFS-lattice.h
	g++ $< $(listol) $(flags) -o $@

algebraic.o: algebraic.cpp algebraic.h mypair.h GNFS.h
	g++ -c $< -o $@

factorbase.o: factorbase.cpp factorbase.h mypair.h GNFS.h
	g++ -c $< -o $@

linear.o: linear.cpp linear.h mypair.h GNFS.h
	g++ -c $< -o $@

NFS.o: NFS.cpp GNFS.h
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

latticesieve.o: latticesieve.cpp latticesieve.h mypair.h GNFS.h GNFS-lattice.h
	g++ -c $< -o $@
