.SUFFIXES:
FC := gfortran
FCFLAGS := $(FCFLAGS) -O5 -Wall -fPIC
LDFLAGS = -L/usr/local/include
LDLIBS = -lopenblas

OS := $(shell uname -s)

DEPSDIR := ./dependencies

DEPSSRC := $(shell find $(DEPSDIR)/ -name \*.f)
XSSRC := xsbbrd.f xsdskb.f
SCOSRC := sco_global.f90 \
  sco_arrays.f90 \
  sco_simpson.f90 \
  sco_mppinv.f90 \
  sco_band_integration.f90 \
  sco_par.f90 \
  sco_model.f90 \
  sco_model_LOG.f90 \
  sco_model_LOGbb.f90 \
  sco_model_LOG_dskb.f90

DEPSOBJ := $(DEPSSRC:.f=.o)
SCOOBJ := $(SCOSRC:.f90=.o)
XSOBJ := $(XSSRC:.f=.o)
#SCOMOD := $(SCOSRC:.f90=.mod)

SRC_LIN := sco_program.f90
SRC_LOG := sco_programLOG.f90
SRC_DISKBB := sco_programDSKB.f90
SRC_BB := sco_program_BB.f90

$(DEPSDIR)/%.o: $(DEPSDIR)/%.f
	$(FC) $(FCFLAGS) -c $< -o $@

./%.o: ./%.f
	$(FC) $(FCFLAGS) -c $< -o $@

./%.o: ./%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

#./%.mod: ./%.f90
#	$(FC) $(FCFLAGS) -c $<


all: loader version vkompth_lin vkompth_log vkompth_bb vkompth_dk makefile_libs vkompthbb vkompthdk vkdualbb vkdualdk vkddka pyvkompth

model: vkompth_lin vkompth_log vkompth_bb vkompth_dk

wrappers: loader version makefile_libs vkompthbb vkompthdk vkdualbb vkdualdk vkddka pyvkompth

loader:
	@echo "\nApply PATHTO to load_vkompth.xcm...\n"
	sed "s|/PATHTO|$$(pwd)|g" load_vkompth.xcm > load_vkompth.tmp; mv load_vkompth.tmp load_vkompth.xcm
	@echo "\n   ... Done.\n\n"

version:
	@echo "\nApply VERSION number $$(cat VERSION) to XSPEC model wrappers...\n"
	for f90 in $(shell ls */*.f90); do echo $$f90; sed s/VERSION/$$(cat VERSION)/g $$f90 > $$f90.tmp; mv $$f90.tmp $$f90; done
	@echo "\n   ... Done.\n\n"

makefile_libs:
	@echo "\nApply LDFLAGS and LDLIBS variables to wrappers Makefile_libs files...\n"
	cd vkompthbb; \
	sed "5s,.*,LDFLAGS=$(LDFLAGS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	sed "6s,.*,LDLIBS=$(LDLIBS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	cd ../vkompthdk; \
	sed "5s,.*,LDFLAGS=$(LDFLAGS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	sed "6s,.*,LDLIBS=$(LDLIBS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	cd ../vkdualbb; \
	sed "5s,.*,LDFLAGS=$(LDFLAGS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	sed "6s,.*,LDLIBS=$(LDLIBS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	cd ../vkdualdk; \
	sed "5s,.*,LDFLAGS=$(LDFLAGS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	sed "6s,.*,LDLIBS=$(LDLIBS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	cd ../vkddka ; \
	sed "5s,.*,LDFLAGS=$(LDFLAGS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	sed "6s,.*,LDLIBS=$(LDLIBS)," Makefile_libs > Makefile_libs.tmp; mv Makefile_libs.tmp Makefile_libs; \
	cd .. ;
	@echo "\n   ... Done.\n\n"

vkompth_lin: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_LIN) $(LDLIBS) -o $@

vkompth_log: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_LOG) $(LDLIBS) -o $@

vkompth_bb: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_BB) $(LDLIBS) -o $@

vkompth_dk: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_DISKBB) $(LDLIBS) -o $@

vkompthbb:
	cd vkompthbb; initpackage vkompthbb lmod_vkompthbb.dat .; cp Makefile_libs Makefile;	hmake; cd ..

vkompthdk:
	cd vkompthdk; initpackage vkompthdk lmod_vkompthdk.dat .; cp Makefile_libs Makefile;	hmake; cd ..

vkdualbb:
	cd vkdualbb; initpackage vkdualbb lmod_vkdualbb.dat .; cp Makefile_libs Makefile;	hmake; cd ..

vkdualdk:
	cd vkdualdk; initpackage vkdualdk lmod_vkdualdk.dat .; cp Makefile_libs Makefile;	hmake; cd ..

vkddka:
	cd vkddka; initpackage vkddka lmod_vkddka.dat .; cp Makefile_libs Makefile;	hmake; cd ..

pyvkompth:
	cd pyvkompth; \
	python -m numpy.f2py $(LDFLAGS) $(LDLIBS) -c pyvkompthbb.pyf pyvkompthbb.f90 \
	   ../sco_arrays.f90 ../sco_global.f90 ../sco_band_integration.f90 \
	   ../sco_mppinv.f90 ../sco_simpson.f90 ../sco_par.f90 \
	   ../sco_model_LOGbb.f90 ../dependencies/*.f ../xsbbrd.f; \
	python -m numpy.f2py $(LDFLAGS) $(LDLIBS) -c pyvkompthdk.pyf pyvkompthdk.f90 \
	   ../sco_arrays.f90 ../sco_global.f90 ../sco_band_integration.f90 \
	   ../sco_mppinv.f90 ../sco_simpson.f90 ../sco_par.f90 \
	   ../sco_model_LOG_dskb.f90 ../dependencies/*.f ../xsdskb.f; cd ..

.PHONY: clean vkompthbb vkompthdk vkdualbb vkdualdk vkddka model wrappers version pyvkompth
clean:
	-rm -f $(DEPSDIR)/*.o *.o *.mod *.dat vkompth_lin vkompth_log vkompth_bb vkompth_dk; \
	cd vkompthbb; hmake clean; cd ..; cd vkompthdk; hmake clean; cd ..; \
	cd vkdualbb; hmake clean;	cd ..; cd vkdualdk; hmake clean; cd ..;  \
	cd vkddka; hmake clean; cd ..; \
	cd pyvkompth; rm pyvkompth*cpython*; cd ..; \
	echo '\n\n   FINISHED make clean all. Ready to re-compile using make. \n'
