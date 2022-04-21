.SUFFIXES:
FC := gfortran
FCFLAGS := $(FCFLAGS) -O5 -Wall -fPIC
LDFLAGS := $(LDFLAGS) -lopenblas

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

SRC_LIN := sco_program.f90
SRC_LOG := sco_programLOG.f90
SRC_DISKBB := sco_programDSKB.f90
SRC_BB :=	sco_program_BB.f90

$(DEPSDIR)/%.o: $(DEPSDIR)/%.f
	$(FC) $(FCFLAGS) -c $< -o $@

./%.o: ./%.f
	$(FC) $(FCFLAGS) -c $< -o $@

./%.o: ./%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

#./%.mod: ./%.f90
#	$(FC) $(FCFLAGS) -c $<

all: vkompth_lin vkompth_log vkompth_bb vkompth_dk

sco_lin: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(TARGET_ARCH) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_LIN) $(LDFLAGS) -o $@

sco_log: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(TARGET_ARCH) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_LOG) $(LDFLAGS) -o $@

sco_log_BB: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(TARGET_ARCH) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_BB) $(LDFLAGS) -o $@

sco_log_diskBB: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(TARGET_ARCH) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_DISKBB) $(LDFLAGS) -o $@

.PHONY: clean
clean:
	-rm -f $(DEPSDIR)/*.o *.o *.mod vkompth_lin vkompth_log vkompth_bb vkompth_dk
