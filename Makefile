FVERSION := $(shell fversion 2>/dev/null || echo FAIL)
ifeq ($(FVERSION),FAIL)
  $(error 'fversion' not found. You should activate HEASOFT first)
endif

HEADASVERSION := $(shell echo $(FVERSION) | sed -E 's/.*V6\.([0-9]+).*/\1/')

VK_VERSION := $(shell cat VERSION 2>/dev/null || echo "1.3.0")

ifeq ($(shell test $(HEADASVERSION) -ge 36 && echo yes),yes)
  EXTRA_FFLAGS := -cpp -DHEASOFT_GE_36 -DVKOMPTH_VER='"$(VK_VERSION)"'
else
  EXTRA_FFLAGS := -cpp -DVKOMPTH_VER='"$(VK_VERSION)"'
endif

.SUFFIXES:
FC := gfortran
FCFLAGS := $(FCFLAGS) -O5 -Wall -fPIC $(EXTRA_FFLAGS)
PYTHON ?= python3
LDFLAGS = -L/usr/local/include
LDLIBS = -lopenblas

OS := $(shell uname -s)

DEPSDIR := ./dependencies
DEPSSRC := $(shell find $(DEPSDIR)/ -name \*.f)
XSSRC   := xsbbrd.f xsdskb.f
SCOSRC  := sco_global.f90 \
  sco_arrays.f90 \
  sco_simpson.f90 \
  sco_mppinv.f90 \
  sco_band_integration.f90 \
  sco_par.f90 \
  sco_model.f90 \
  sco_model_LOG.f90 \
  sco_model_LOGbb.f90 \
  sco_model_LOG_dskb.f90 \
  sco_model_LOG_dskb_dL.f90

DEPSOBJ := $(DEPSSRC:.f=.o)
SCOOBJ  := $(SCOSRC:.f90=.o)
XSOBJ   := $(XSSRC:.f=.o)

SRC_LIN    := sco_program.f90
SRC_LOG    := sco_programLOG.f90
SRC_DISKBB := sco_programDSKB.f90
SRC_BB     := sco_program_BB.f90
SRC_DL     := sco_programDSKB_dL.f90

WRAPPERS := vkompthbb vkompthdk vkompthdl vkdualbb vkdualdk vkdualdl vkddka

$(DEPSDIR)/%.o: $(DEPSDIR)/%.f
	$(FC) $(FCFLAGS) -c $< -o $@

./%.o: ./%.f
	$(FC) $(FCFLAGS) -c $< -o $@

./%.o: ./%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

all: fversion loader version vkompth_lin vkompth_log vkompth_bb vkompth_dk vkompth_dl makefile_libs $(WRAPPERS) pyvkompth

model: vkompth_lin vkompth_log vkompth_bb vkompth_dk vkompth_dl

wrappers: loader version makefile_libs $(WRAPPERS) pyvkompth

fversion:
	@echo "Running HEASOFT version $(HEADASVERSION) (vKompth v$(VK_VERSION))..."

loader:
	@echo "\nGenerating load_vkompth.xcm from load_vkompth.xcm.in..."
	@sed "s|/PATHTO|$$(pwd)|g" load_vkompth.xcm.in > load_vkompth.xcm
	@echo "   ... Done.\n"

version:
	@echo "\nvKompth target version: $(VK_VERSION)\n"

makefile_libs:
	@echo "\nPreparing build configurations for XSPEC wrappers..."
	@for dir in $(WRAPPERS); do \
		cd $$dir; \
		if [ $(HEADASVERSION) -ge 37 ]; then \
			sed -e "s|^LDFLAGS=.*|LDFLAGS=$(LDFLAGS)|" \
			    -e "s|^LDLIBS=.*|LDLIBS=$(LDLIBS)|" \
			    -e "s| -DINITPACKAGE| -DINITPACKAGE -cpp -DHEASOFT_GE_36 -DVKOMPTH_VER='\"$(VK_VERSION)\"'|" \
			    -e "s/ -lXS / /g" \
			    Makefile_libs36 > Makefile.temp; \
		elif [ $(HEADASVERSION) -eq 36 ]; then \
			sed -e "s|^LDFLAGS=.*|LDFLAGS=$(LDFLAGS)|" \
			    -e "s|^LDLIBS=.*|LDLIBS=$(LDLIBS)|" \
			    -e "s| -DINITPACKAGE| -DINITPACKAGE -cpp -DHEASOFT_GE_36 -DVKOMPTH_VER='\"$(VK_VERSION)\"'|" \
			    Makefile_libs36 > Makefile.temp; \
		else \
			sed -e "s|^LDFLAGS=.*|LDFLAGS=$(LDFLAGS)|" \
			    -e "s|^LDLIBS=.*|LDLIBS=$(LDLIBS)|" \
			    -e "s| -DINITPACKAGE| -DINITPACKAGE -cpp -DVKOMPTH_VER='\"$(VK_VERSION)\"'|" \
			    Makefile_libs > Makefile.temp; \
		fi; \
		cd ..; \
	done
	@echo "   ... Done.\n"

vkompth_lin: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_LIN) $(LDLIBS) -o $@

vkompth_log: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_LOG) $(LDLIBS) -o $@

vkompth_bb: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_BB) $(LDLIBS) -o $@

vkompth_dk: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_DISKBB) $(LDLIBS) -o $@

vkompth_dl: $(DEPSOBJ) $(SCOOBJ) $(XSOBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(DEPSOBJ) $(SCOOBJ) $(XSOBJ) $(SRC_DISKBB) $(LDLIBS) -o $@

vkompthbb:
	cd vkompthbb && initpackage vkompthbb lmod_vkompthbb.dat . && mv Makefile.temp Makefile && hmake && cd ..

vkompthdk:
	cd vkompthdk && initpackage vkompthdk lmod_vkompthdk.dat . && mv Makefile.temp Makefile && hmake && cd ..

vkompthdl:
	cd vkompthdl && initpackage vkompthdl lmod_vkompthdl.dat . && mv Makefile.temp Makefile && hmake && cd ..

vkdualbb:
	cd vkdualbb && initpackage vkdualbb lmod_vkdualbb.dat . && mv Makefile.temp Makefile && hmake && cd ..

vkdualdk:
	cd vkdualdk && initpackage vkdualdk lmod_vkdualdk.dat . && mv Makefile.temp Makefile && hmake && cd ..

vkdualdl:
	cd vkdualdl && initpackage vkdualdl lmod_vkdualdl.dat . && mv Makefile.temp Makefile && hmake && cd ..

vkddka:
	cd vkddka && initpackage vkddka lmod_vkddka.dat . && mv Makefile.temp Makefile && hmake && cd ..

pyvkompth:
	cd pyvkompth; \
	$(PYTHON) -m numpy.f2py $(LDFLAGS) $(LDLIBS) --f90flags="$(EXTRA_FFLAGS)" -c pyvkompthbb.pyf pyvkompthbb.f90 \
	   ../sco_arrays.f90 ../sco_global.f90 ../sco_band_integration.f90 \
	   ../sco_mppinv.f90 ../sco_simpson.f90 ../sco_par.f90 \
	   ../sco_model_LOGbb.f90 ../dependencies/*.f ../xsbbrd.f; \
	$(PYTHON) -m numpy.f2py $(LDFLAGS) $(LDLIBS) --f90flags="$(EXTRA_FFLAGS)" -c pyvkompthdk.pyf pyvkompthdk.f90 \
	   ../sco_arrays.f90 ../sco_global.f90 ../sco_band_integration.f90 \
	   ../sco_mppinv.f90 ../sco_simpson.f90 ../sco_par.f90 \
	   ../sco_model_LOG_dskb.f90 ../dependencies/*.f ../xsdskb.f; \
	$(PYTHON) -m numpy.f2py $(LDFLAGS) $(LDLIBS) --f90flags="$(EXTRA_FFLAGS)" -c pyvkompthdl.pyf pyvkompthdl.f90 \
	   ../sco_arrays.f90 ../sco_global.f90 ../sco_band_integration.f90 \
	   ../sco_mppinv.f90 ../sco_simpson.f90 ../sco_par.f90 \
	   ../sco_model_LOG_dskb_dL.f90 ../dependencies/*.f ../xsdskb.f; cd ..

.PHONY: clean vkompthbb vkompthdk vkdualbb vkdualdk vkdualdl vkompthdl vkddka model wrappers version pyvkompth fversion loader makefile_libs

clean:
	@rm -f $(DEPSDIR)/*.o *.o *.mod *.dat vkompth_lin vkompth_log vkompth_bb vkompth_dk vkompth_dl load_vkompth.xcm
	@for dir in $(WRAPPERS); do \
		if [ -d "$$dir" ]; then \
			cd $$dir && (hmake clean >/dev/null 2>&1 || true) && \
			rm -f Makefile Makefile.temp *.so *.o *.cxx *.h *.tcl *.mod *.bck && \
			cd ..; \
		fi; \
	done
	@if [ -d "pyvkompth" ]; then \
		cd pyvkompth && rm -f *.so *.o *.cpython* *.pyc && cd ..; \
	fi
	@echo "\n   FINISHED make clean. Working tree clean. Ready to re-compile using make.\n"
