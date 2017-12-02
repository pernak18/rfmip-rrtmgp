RRTMGP_ROOT = /Users/robert/Codes/RRTMGP/trunk
BUILD_DIR = $(RRTMGP_ROOT)/build
include $(BUILD_DIR)/Makefile.conf
include $(BUILD_DIR)/Makefile.rules

# Override default rules
# %: %.o
# 	$(ECOMPILE) -o ../$@ $^ $(LDFLAGS) $(LIBS)
%.o: %.f
	$(FC) $(F77FLAGS) $(FCINCLUDE) -c  $<

VPATH = ./:$(BUILD_DIR):$(RRTMGP_ROOT)/extensions::$(RRTMGP_ROOT)/test/util/src/io

#
# RRTMGP library, module files
#
LDFLAGS   += -L$(BUILD_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(BUILD_DIR)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS = mo_rfmip_io.o mo_load_coefficients.o

F77FLAGS += -O2

rrtmgp_rfmip_lw: rrtmgp_rfmip_lw.o $(ADDITIONS)

rrtmgp_rfmip_lw.o: rrtmgp_rfmip_lw.F90 $(ADDITIONS)

clean:
	-rm ../shdompp_solve ../change_albedo_emis_sza *.o *.mod *.optrpt ../*.optrpt
