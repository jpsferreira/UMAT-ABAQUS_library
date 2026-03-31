# Makefile for UMAT Builder modular library
#
# Usage:
#   make              Build all object files for linking with ABAQUS
#   make clean        Remove build artifacts
#
# For ABAQUS (single-file submission):
#   make abaqus_single   Concatenate all sources into one .f90 for abaqus user=
#
# For ABAQUS (pre-compiled):
#   abaqus job=mymodel user=build/umat_all.f90

FC       = gfortran
FFLAGS   = -O2 -ffree-form -fPIC
BUILDDIR = build

# Module compilation order (dependency chain)
MODULES = \
	$(BUILDDIR)/mod_constants.o \
	$(BUILDDIR)/mod_tensor.o \
	$(BUILDDIR)/mod_kinematics.o \
	$(BUILDDIR)/mod_continuum.o \
	$(BUILDDIR)/mod_hyperelastic.o \
	$(BUILDDIR)/mod_anisotropic.o \
	$(BUILDDIR)/mod_network.o \
	$(BUILDDIR)/mod_damage.o \
	$(BUILDDIR)/mod_viscosity.o

UMAT_OBJ = $(BUILDDIR)/umat_builder.o $(BUILDDIR)/uexternaldb.o

.PHONY: all clean abaqus_single

all: $(BUILDDIR) $(MODULES) $(UMAT_OBJ)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Module dependencies (order matters)
$(BUILDDIR)/mod_constants.o: mod_constants.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILDDIR)

$(BUILDDIR)/mod_tensor.o: mod_tensor.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_kinematics.o: mod_kinematics.f90 $(BUILDDIR)/mod_tensor.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_continuum.o: mod_continuum.f90 $(BUILDDIR)/mod_tensor.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_hyperelastic.o: mod_hyperelastic.f90 $(BUILDDIR)/mod_tensor.o $(BUILDDIR)/mod_kinematics.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_anisotropic.o: mod_anisotropic.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_network.o: mod_network.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_damage.o: mod_damage.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/mod_viscosity.o: mod_viscosity.f90 $(BUILDDIR)/mod_constants.o
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/umat_builder.o: umat_builder.f90 $(MODULES)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/uexternaldb.o: uexternaldb.f90
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

# Concatenate all sources into a single file for ABAQUS user= submission
abaqus_single: $(BUILDDIR)
	cat mod_constants.f90 mod_tensor.f90 mod_kinematics.f90 \
	    mod_continuum.f90 mod_hyperelastic.f90 mod_anisotropic.f90 \
	    mod_network.f90 mod_damage.f90 mod_viscosity.f90 \
	    umat_builder.f90 uexternaldb.f90 > $(BUILDDIR)/umat_all.f90
	@echo "Created $(BUILDDIR)/umat_all.f90 — submit to ABAQUS with: abaqus job=X user=$(BUILDDIR)/umat_all.f90"

clean:
	rm -rf $(BUILDDIR)
