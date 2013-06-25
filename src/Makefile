include ../make.inc

SOURCES=nulib.F90 \
	fermi.F90 \
	helpers.F90 \
	gauss_laguerre_helpers.F90 \
	gauss_legendre_helpers.F90 \
	electron_positron_annihilation.F90 \
	absorption_crosssections.F90 \
	emissivities.F90 \
	scattering.F90 \
	weak_magnetism_correction.F90 \
	inelastic_electron_scattering.F90 \
	nuclei_distribution_helpers.F90

NT_SOURCES=nulibtable.F90 \
	   nulibtable_reader.F90 \
	   linterp_many_mod.F90

ifeq ($(NUCLEI_HEMPEL),1)
	F_SOURCES = extra_code_and_tables/sfho_frdm_comp/sfho_frdm_composition_module.f
	F_OBJECTS=$(F_SOURCES:.f=.o )
	MODINC += -I./extra_code_and_tables/sfho_frdm_comp
	DEFS += -DNUCLEI_HEMPEL
endif

EXTRADEPS = requested_interactions.inc constants.inc

OBJECTS=$(SOURCES:.F90=.o )
NT_OBJECTS=$(NT_SOURCES:.F90=.o )

#take care of EOS dependences etc
DEFS += -DHAVE_NUC_EOS
EXTRAINCS = $(HDF5INCS) -I./nuc_eos
EXTRADEPS += nuc_eos/nuc_eos.a
EXTRAOBJECTS = nuc_eos/nuc_eos.a $(HDF5LIBS)

all: nulibtable_driver point_example make_table_example

nulibtable_driver:  $(EXTRADEPS) $(NT_OBJECTS) $(F_OBJECTS) 
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINC) -o ../nulibtable_driver nulibtable_driver.F90 $(NT_OBJECTS) $(EXTRAOBJECTS)

point_example:  $(EXTRADEPS) $(F_OBJECTS) $(OBJECTS) point_example.F90
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -o ../point_example point_example.F90 $(OBJECTS) $(F_OBJECTS) $(EXTRAOBJECTS)

make_table_example:  $(EXTRADEPS) $(F_OBJECTS) $(OBJECTS) make_table_example.F90
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -o ../make_table_example make_table_example.F90 $(OBJECTS) $(F_OBJECTS) $(EXTRAOBJECTS)

$(OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

$(F_OBJECTS): %.o: %.f 
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@	

$(NT_OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

nuc_eos/nuc_eos.a: nuc_eos/*.F90 nuc_eos/*.f
	$(MAKE) -C nuc_eos

clean:
	rm -rf ../make_table_example
	rm -rf ../point_example
	rm -rf ../nulibtable_driver
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.a
	$(MAKE) -C nuc_eos clean	
	rm -rf ./extra_code_and_tables/sfho_frdm_comp/sfho_frdm_composition_module.o
