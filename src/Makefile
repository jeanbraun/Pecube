# ═══════════════════════════════════════════════════════════════════════════════
# Pecube Makefile — cross-platform (Linux, macOS, Windows)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Usage:
#   make all          build all binaries (serial + MPI)
#   make Pecube       serial Pecube only
#   make MPI          MPI-enabled PecubeMPI only
#   make clean        remove object files and binaries
#
# Compiler selection (in order of priority):
#   1. Command-line:   make FC=ifort MPIFC=mpiifort all
#   2. Environment:    set automatically by conda-forge build environments
#                      and HPC environment modules (module load gcc openmpi)
#   3. Defaults below: used only when make's built-in defaults are in effect
# ═══════════════════════════════════════════════════════════════════════════════

# ── Compilers ──────────────────────────────────────────────────────────────────
# $(origin VAR) returns 'default' when make's built-in value is in effect
# (e.g. FC=f77, CC=cc). We replace these with sensible defaults while
# respecting any value set by the environment (conda-forge, module system)
# or on the command line.

ifeq ($(origin FC), default)
    FC := gfortran
endif
ifeq ($(origin F77), default)
    F77 := $(FC)
endif
ifeq ($(origin CC), default)
    CC := gcc
endif

MPIFC  ?= mpif90
MPIF77 ?= mpif77
MPICC  ?= mpicc

# Internal aliases (keep targets readable)
FF90     = $(FC)
FF77     = $(F77)
LINK     = $(FC)
FF90MPI  = $(MPIFC)
FF77MPI  = $(MPIF77)
CCMPI    = $(MPICC)

# ── Compiler flags ─────────────────────────────────────────────────────────────
# Flags are selected based on the compiler family, not the platform.
#
# F77FLAGS:    used for serial .f/.f77 compilation
#              -fd-lines-as-comments: treats lines starting with 'd' as comments,
#              which disables the MPI calls in the serial NA code (naScalar.f)
#
# F77MPIFLAGS: used for MPI .f compilation (naMPI.f only)
#              -fd-lines-as-code: treats 'd' lines as active code,
#              enabling the MPI calls in naMPI.f
#
# This 'd'-line mechanism is how the same source selectively includes/excludes
# MPI calls for serial vs parallel builds.

OPT     ?= -O
LDFLAGS ?=
CFLAGS  ?=

# FFLAGS is also passed at link time so that driver flags set by conda-forge
# activation scripts (e.g. -fuse-ld=lld, -fms-runtime-lib=dll on Windows)
# take effect during linking.

ifneq (,$(findstring gfortran, $(FC)))
    FFLAGS      ?= -fbounds-check -fallow-invalid-boz
    F77FLAGS    := $(FFLAGS) -fd-lines-as-comments
    F77MPIFLAGS := $(FFLAGS) -fd-lines-as-code
else ifneq (,$(findstring ifort, $(FC)))
    FFLAGS      ?= -check bounds
    F77FLAGS    := $(FFLAGS)
    F77MPIFLAGS := $(FFLAGS)
else
    # flang or unknown compiler: no special flags
    # (conda-forge sets FFLAGS via activation scripts, respected via ?=)
    FFLAGS      ?=
    F77FLAGS    := $(FFLAGS)
    F77MPIFLAGS := $(FFLAGS)
    # Use clang for C when using flang — same LLVM toolchain, MSVC-ABI compatible
    CC := clang
    # Locate and link the flang Fortran runtime (provides dtime_ and other intrinsics)
    FLANG_RT := $(firstword $(wildcard \
        $(subst \,/,$(CONDA_PREFIX))/Library/lib/clang/*/lib/x86_64-pc-windows-msvc/flang_rt.runtime.dynamic.lib))
    ifneq ($(FLANG_RT),)
        LDFLAGS += -Wl,-defaultlib:$(FLANG_RT)
    endif
endif

# ── Object lists ───────────────────────────────────────────────────────────────
# OBJECTS_MPI: same as OBJECTS but without naScalar.o — the MPI build
# replaces it with naMPI.o (compiled separately in the MPI target).

# Compiler-specific compatibility shims
ifneq (,$(findstring flang, $(FC)))
    EXTRA_OBJS = dtime_compat.o
else
    EXTRA_OBJS =
endif

OBJECTS = \
    module_Pecube.o \
    Pecube_utils.o \
    Trapped_charge.o \
    Mad_He.o \
    Mad_Trax.o \
    Mad_Trax_Zircon.o \
    Pecube.o \
    calculate_fault_parameters.o \
    calculate_ages.o \
    create_pecube_in.o \
    find_distance.o \
    find_dt.o \
    find_element.o \
    find_neighbours.o \
    find_temp.o \
    find_triangle.o \
    find_velo.o \
    find_velocity.o \
    four1.o \
    heat_production.o \
    hinterpolate.o \
    interpol3d.o \
    interpolate.o \
    isostatic_rebound.o \
    make_matrix.o \
    move_fault.o \
    naScalar.o \
    read_in_fault_parameters.o \
    realft.o \
    screen_counter.o \
    sinft.o \
    solve_iterative.o \
    thermal_diffusivity.o \
    tridag.o \
    KETCH.o \
    calculate_Ttpath.o \
    extract_Ttpath.o \
    Diffusion1D.o \
    element_surface.o \
    read_input_file.o \
    He_diffusion.o \
    scanfile.o \
    read_data_files.o \
    read_table.o \
    Ketch07.o \
    RD_Models.o \
    ExtractDEM.o \
    initialize_random_seed.o \
    calculate_TL.o \
    calculate_OSL.o \
    calculate_ESR.o \
    smooth_topo.o \
    Topography.o

OBJECTS_MPI = $(filter-out naScalar.o, $(OBJECTS))

OBJECTS_TEST = \
    testFault.o \
    calculate_fault_parameters.o \
    find_distance.o \
    find_velocity.o \
    move_fault.o \
    read_in_fault_parameters.o \
    read_input_file.o \
    read_data_files.o \
    read_table.o \
    scanfile.o \
    testTopo.o \
    ExtractDEM.o \
    Pecube_utils.o \
    smooth_topo.o \
    Topography.o

OBJECTS_VTK = \
    Pecube_utils.o \
    vtk.o

OBJECTS_PECUBE324 = \
    Pecube324.o

OBJECT_AGE = \
    CalcAge.o \
    He_diffusion.o \
    Mad_He.o \
    interpolate.o \
    Ketch07.o \
    KETCH.o \
    RD_Models.o \
    Pecube_utils.o \
    tridag.o

# ── Suffix rules ───────────────────────────────────────────────────────────────
.SUFFIXES: .o .f90 .F90 .f95 .f .f77 .c

.f90.o:
	$(FF90) -c $(OPT) $(FFLAGS) $*.f90

.F90.o:
	$(FF90) -c $(OPT) $(FFLAGS) $*.F90

.f95.o:
	$(FF90) -c $(OPT) $(FFLAGS) $*.f95

.f.o:
	$(FF77) -c $(OPT) $(F77FLAGS) $*.f

.f77.o:
	$(FF77) -c $(OPT) $(F77FLAGS) $*.f77

.c.o:
	$(CC) -c $(OPT) $(CFLAGS) $*.c

# ── Targets ────────────────────────────────────────────────────────────────────
.PHONY: all Pecube MPI Test Vtk Pecube324 Age_specific clean

# Build MPI target only if the MPI Fortran wrapper is available
HAS_MPI := $(shell where $(MPIF77) 2>nul || which $(MPIF77) 2>/dev/null)
ifneq ($(HAS_MPI),)
all: Pecube Test Vtk Pecube324 Age_specific MPI
else
all: Pecube Test Vtk Pecube324 Age_specific
    $(info Note: $(MPIF77) not found — skipping MPI build. Run 'make MPI' explicitly when MPI is available.)
endif

Pecube: $(OBJECTS) $(EXTRA_OBJS)
	$(LINK) $(FFLAGS) -o ../bin/Pecube $(OBJECTS) $(EXTRA_OBJS) $(LDFLAGS)

# naMPI.f compiled with the MPI F77 wrapper using F77MPIFLAGS so that
# 'd'-prefixed lines are treated as active code (enabling MPI calls).
# All other objects are reused from the serial build.
MPI: $(OBJECTS_MPI) $(EXTRA_OBJS)
	$(FF77MPI) -c $(OPT) $(F77MPIFLAGS) naMPI.f
	$(FF90MPI) $(FFLAGS) -o ../bin/PecubeMPI $(OBJECTS_MPI) $(EXTRA_OBJS) naMPI.o $(LDFLAGS)

Test: $(OBJECTS_TEST) $(EXTRA_OBJS)
	$(LINK) $(FFLAGS) -o ../bin/Test $(OBJECTS_TEST) $(EXTRA_OBJS) $(LDFLAGS)

Vtk: $(OBJECTS_VTK) $(EXTRA_OBJS)
	$(LINK) $(FFLAGS) -o ../bin/Vtk $(OBJECTS_VTK) $(EXTRA_OBJS) $(LDFLAGS)

Pecube324: $(OBJECTS_PECUBE324) $(EXTRA_OBJS)
	$(LINK) $(FFLAGS) -o ../bin/Pecube324 $(OBJECTS_PECUBE324) $(EXTRA_OBJS) $(LDFLAGS)

Age_specific: $(OBJECT_AGE) $(EXTRA_OBJS)
	$(LINK) $(FFLAGS) -o ../bin/CalcAge $(OBJECT_AGE) $(EXTRA_OBJS) $(LDFLAGS)

# ── Clean ──────────────────────────────────────────────────────────────────────
clean:
	-$(RM) *.o
	-$(RM) ../bin/Pecube ../bin/PecubeMPI ../bin/Test ../bin/Vtk \
	       ../bin/Pecube324 ../bin/CalcAge
