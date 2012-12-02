#
#  On Chili and Kachina, this is where the F90 compiler is:
#
F90 = /usr/bin/fort

sge_mod_prb.out: sge_mod_prb
	./sge_mod_prb > sge_mod_prb.out
sge_mod_prb: sge_mod_prb.o r_swap.o rvec_print.o sge_create.o sge_delete.o\
  sge_det.o sge_fa.o sge_sl.o sge_slt.o timestamp.o sge_module.o
	$(F90) -o sge_mod_prb *.o
sge_mod_prb.o: sge_mod_prb.f90 sge_module.o
	$(F90) -c sge_mod_prb.f90
r_swap.o: r_swap.f90
	$(F90) -c r_swap.f90
rvec_print.o: rvec_print.f90
	$(F90) -c rvec_print.f90
sge_create.o: sge_create.f90 sge_module.o
	$(F90) -c sge_create.f90
sge_delete.o: sge_delete.f90 sge_module.o
	$(F90) -c sge_delete.f90
sge_det.o: sge_det.f90 sge_module.o
	$(F90) -c sge_det.f90
sge_fa.o: sge_fa.f90 sge_module.o
	$(F90) -c sge_fa.f90
sge_sl.o: sge_sl.f90 sge_module.o
	$(F90) -c sge_sl.f90
sge_slt.o: sge_slt.f90 sge_module.o
	$(F90) -c rvec_print.f90
timestamp.o: timestamp.f90
	$(F90) -c timestamp.f90
sge_module.o: sge_module.f90
	$(F90) -c sge_module.f90
clean:
	/bin/rm sge_mod_prb.out sge_mod_prb *.o sge_module.mod



