include ../make.inc

MODULES = nulib.F90 \
	  input_parser.F90

SOURCES = fermi.F90 \
	  helpers.F90 \
	  gauss_laguerre_helpers.F90 \
	  gauss_legendre_helpers.F90 \
	  electron_positron_annihilation.F90 \
	  absorption_crosssections.F90 \
	  emissivities.F90 \
	  scattering.F90 \
	  weak_magnetism_correction.F90 \
	  inelastic_electron_scattering.F90 \
	  eos_interface.F90

NT_SOURCES=nulibtable.F90 \
	   nulibtable_reader.F90 \
	   linterp_many_mod.F90


ifeq ($(NUCLEI_HEMPEL),1)
	MODULES += nuclei_distribution_helpers.F90
	F_SOURCES = extra_code_and_tables/sfho_frdm_comp/sfho_frdm_composition_module.f
	F_OBJECTS=$(F_SOURCES:.f=.o )
	MODINC += -I./extra_code_and_tables/sfho_frdm_comp
	DEFS += -DNUCLEI_HEMPEL 
endif

F2PY := $(shell which f2py)
F90WRAP := $(shell which f90wrap)
ifeq ($(PYTHON_BINDINGS),1)
ifeq ($(F2PY),)
$(error "Build python bindings flag is set, but f2py is not available on this system.")
endif
ifeq ($(F90WRAP),)
$(error "Build python bindings flag is set, but f90wrap is not available on this system.")
endif
endif

#the weak_rates module requires M. Hempel's
#EOS dependent NSE distributions
ifeq ($(WEAK_RATES),1)
ifeq ($(NUCLEI_HEMPEL),0)
$(error NUCLEI_HEMPEL must be 1 for weakrates)
endif

	MODULES += weakrates/ratetable.F90 \
		   weakrates/rateapprox.F90 \
		   weakrates/ratelibrary.F90 \
		   weakrates/interface.F90 \

	DEFS += -DWEAK_RATES
	SOURCES += pyinterface.F90

endif


#dual openmp+mpi jobs may work, but are not currently supported
ifeq ($(OPENMP),1)
ifeq ($(MPI),1)
$(error MPI and OPENMP cannot be used simultaneously)
endif
endif

ifeq ($(OPENMP),1)
ifeq ($(F90),gfortran)
	F90FLAGS += -fopenmp
else
	F90FLAGS += -openmp
endif	
endif

ifeq ($(MPI),1)
	F90 = mpif90
	DEFS += -D__MPI__
endif

CFLAGS=-Wall -O3 -fPIC
MOD_OBJECTS=$(MODULES:.F90=.o )
OBJECTS=$(SOURCES:.F90=.o )
NT_OBJECTS=$(NT_SOURCES:.F90=.o )
F90FLAGS += -fPIC
LDFLAGS=-shared

#take care of EOS dependences etc
EXTRAINCS = -I./nuc_eos
EXTRADEPS = nuc_eos/nuc_eos.a
EXTRAOBJECTS = nuc_eos/nuc_eos.a

#extra C library dependence
EXTRACDEPS = ./weakrates/Clibs.a
CLIB_SRC = ./weakrates/calc_delta_e.c
CLIB_OBJECTS = $(CLIB_SRC:.c=.o)

ifeq ($(HELMHOLTZ_EOS),1)
	DEFS += -DHELMHOLTZ_EOS
	f90_SOURCES = other_eos/helmholtz/helmholtz.f90
	f90_OBJECTS=$(f90_SOURCES:.f90=.o )
else
	DEFS += -DHAVE_NUC_EOS
	f90_SOURCES =
	f90_OBJECTS=
endif
EXTRADEPS += requested_interactions.inc constants.inc ../make.inc Makefile
EXTRAINCS += $(HDF5INCS)
EXTRAOBJECTS += $(HDF5LIBS) -Wl,-rpath,$(HDF5DIR)/lib

ifeq ($(WEAK_RATES),1)
all: nulib.a ../nulibtable_driver ../point_example ../make_table_weakrates ../libNuLib.so PyNuLib
else
all: nulib.a ../nulibtable_driver ../point_example ../make_table_example ../libNuLib.so
endif

../nulibtable_driver:  $(EXTRADEPS) $(NT_OBJECTS) $(F_OBJECTS) $(f90_OBJECTS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINC) -o $@ nulibtable_driver.F90 $(NT_OBJECTS) $(EXTRAOBJECTS) $(f90_OBJECTS)

../point_example:  $(EXTRADEPS) $(F_OBJECTS) $(MOD_OBJECTS) $(EXTRACDEPS) $(OBJECTS)  $(f90_OBJECTS)point_example.F90
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -o $@ point_example.F90 $(MOD_OBJECTS) $(EXTRACDEPS) $(OBJECTS) $(F_OBJECTS) $(EXTRAOBJECTS) $(f90_OBJECTS)

../make_table_example:  $(EXTRADEPS) $(F_OBJECTS) $(MOD_OBJECTS) $(EXTRACDEPS) $(OBJECTS)  $(f90_OBJECTS)make_table_example.F90
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -o $@ make_table_example.F90 $(MOD_OBJECTS) $(EXTRACDEPS) $(OBJECTS) $(F_OBJECTS) $(EXTRAOBJECTS) $(f90_OBJECTS)

../make_table_weakrates:  $(EXTRADEPS) $(F_OBJECTS) $(MOD_OBJECTS) $(EXTRACDEPS) $(OBJECTS)  $(f90_OBJECTS)make_table_weakrates.F90
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -o $@ make_table_weakrates.F90 $(MOD_OBJECTS) $(EXTRACDEPS) $(OBJECTS) $(F_OBJECTS) $(EXTRAOBJECTS) $(f90_OBJECTS)

nulib.a: $(EXTRADEPS) $(F_OBJECTS) $(MOD_OBJECTS) $(OBJECTS) $(NT_OBJECTS) $(CLIB_OBJECTS)
	ar -r nulib.a *.o nuc_eos/*.o

../libNuLib.so: $(EXTRADEPS) $(F_OBJECTS) $(MOD_OBJECTS) $(EXTRACDEPS) $(OBJECTS) $(NT_OBJECTS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) $(LDFLAGS) -o $@ $(MOD_OBJECTS) $(EXTRACDEPS) $(OBJECTS) $(F_OBJECTS) $(EXTRAOBJECTS)

./weakrates/Clibs.a: $(CLIB_OBJECTS)
	ar -rvs ./weakrates/Clibs.a $(CLIB_OBJECTS)
$(CLIB_OBJECTS): %.o: %.c
	gcc $(CFLAGS) -c $< -o $@

PyNuLib:
ifeq ($(PYTHON_BINDINGS),1)
	$(F90WRAP) --abort-func "abort; print *," -m $@ pyinterface.F90
	$(F90) $(F90FLAGS) $(EXTRAINCS) $(MODINC) f90wrap_*.f90 -c
	$(F2PY) -m _$@ -h $@.pyf f90wrap_*.f90
	LDFLAGS="$(LDFLAGS) -Wl,-rpath,$(HDF5DIR)/lib -Wl,-rpath,\\\\\$$ORIGIN" $(F2PY)  --f90flags="-ffree-line-length-none" --f77flags="-ffree-line-length-none -ffixed-line-length-none" -c $@.pyf pyinterface.o f90wrap_pyinterface.o $(HDF5LIBS) -L.. -lNuLib
	mv _$@.so ..
	mv $@.py ..
endif

test:	../libNuLib.so  test.F90
	$(F90) $(F90FLAGS) $(DEFS) -o test $^ -Wl,-rpath,src

$(OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

$(MOD_OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

$(F_OBJECTS): %.o: %.f 
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@	

$(f90_OBJECTS): %.o: %.f90 
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@	

$(NT_OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

nuc_eos/nuc_eos.a: nuc_eos/*.F90 nuc_eos/*.f
	$(MAKE) -C nuc_eos

clean:
	rm -rf ../make_table_example
	rm -rf ../make_table_weakrates
	rm -rf ../point_example
	rm -rf ../nulibtable_driver
	rm -rf test
	rm -rf *.o
	rm -rf *.so
	rm -rf *~
	rm -rf *.mod
	rm -rf *.a
	rm -rf *.pyf
	rm -rf *.pyc
	rm -rf ../*.pyf
	rm -rf ../*.pyc
	rm -rf ../__pycache__
	rm -rf ../PyNuLib.py
	rm -rf ../*.so
	rm -rf f90wrap*
	$(MAKE) -C nuc_eos clean
	rm -rf ./extra_code_and_tables/sfho_frdm_comp/sfho_frdm_composition_module.o
	rm -rf ./other_eos/helmholtz/*.o
	rm -rf ./weakrates/*.o
	rm -rf ./weakrates/*.a
