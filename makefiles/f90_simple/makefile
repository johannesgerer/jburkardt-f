#
#  On Chili and Kachina, this is where the F90 compiler is:
#
F90 = /usr/bin/fort

f90_simple.out: f90_simple
	./f90_simple > f90_simple.out
f90_simple: f90_simple.o midpoint.o f.o
	$(F90) -o f90_simple f90_simple.o midpoint.o f.o
f90_simple.o: f90_simple.f90
	$(F90) -c f90_simple.f90
midpoint.o: midpoint.f90
	$(F90) -c midpoint.f90
f.o: f.f90
	$(F90) -c f.f90
clean:
	/bin/rm f90_simple.out f90_simple f90_simple.o midpoint.o f.o



