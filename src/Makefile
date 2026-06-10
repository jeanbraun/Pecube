.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .f95 .f77 /mpif90.e .r .y .yr .ye .l .p .sh .csh .h

# Check system
ifeq ($(OS), Windows_NT)
    SYSTEM := Windows
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S), Linux)
        SYSTEM := Linux
    else ifeq ($(UNAME_S), Darwin)
        SYSTEM := MacOS
    else
        SYSTEM := Unknown
    endif
endif

# Check if windows or unix system

# use these lines to produce a scalar version using the gfortran compiler
# this removes the mpi commands in NA

FLAGS = -c -O
FLAGS = -c -O -g
#FF90 = gfortran -fsanitize=address -fbounds-check -fcheck=all -fbacktrace -Wextra
FF90 = gfortran -fbounds-check
FF77 = gfortran -fd-lines-as-comments -fbounds-check
CC = gcc -fbounds-check -static-libgcc -static-libstdc++ -lquadmath
LINK = gfortran -fbounds-check -static-libgfortran -static-libgcc -static-libstdc++

# use these lines to produce a parallel (mpi) version to parallelize the inversion
# Note that some compilers will not allow to have an "include" command preceded by a "d" as an indication
#   of a commented line; if this is the case for your compiler, the routine na.f will not compile properly
#   You will need to edit the na.f file and remove the "d" in the first column of three lines that
#   start with include 'mpif.h'

FLAGSMPI = -c -O
ifeq ($(SYSTEM), Windows)
	FF90MPI := mpif90 -fallow-invalid-boz
	FF77MPI = mpif77 -fd-lines-as-code 
else
	FF90MPI := mpif90
	FF77MPI = mpif77 -fd-lines-as-code
endif
CCMPI = mpicc

#FLAGS = -c -O
#FF90 = ifort
#FF77 = ifort
#CC = icc

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


OBJECTS_MPI = \
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
calculate_TL.o\
calculate_OSL.o\
calculate_ESR.o \
smooth_topo.o \
Topography.o

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
ExtractDEM.o\
Pecube_utils.o \
smooth_topo.o\
Topography.o

OBJECTS_VTK = \
Pecube_utils.o \
vtk.o

OBJECTS_PECUBE324 = \
Pecube324.o

.f95.o:
	$(FF90) $(FLAGS) $*.f95

.f90.o:
	$(FF90) $(FLAGS) $*.f90

.F90.o:
	$(FF90) $(FLAGS) $*.F90

.f.o:
	$(FF77) $(FLAGS) $*.f

.f77.o:
	$(F77) $(FLAGS) $*.f77

.c.o:
	$(CC) $(FLAGS) $*.c

Pecube:	$(OBJECTS)
	$(LINK) -O $(OBJECTS) -o ../bin/Pecube

MPI: $(OBJECTS_MPI)
	mpif77 -c -O -fallow-invalid-boz naMPI.f 
	$(FF90MPI) -O $(OBJECTS_MPI) naMPI.o -o ../bin/PecubeMPI

Test: $(OBJECTS_TEST)
	$(LINK) -O $(OBJECTS_TEST) -o ../bin/Test

Vtk: $(OBJECTS_VTK)
	$(LINK) -O $(OBJECTS_VTK) -o ../bin/Vtk

Pecube324:	$(OBJECTS_PECUBE324)
	$(LINK) -O $(OBJECTS_PECUBE324) -o ../bin/Pecube324
	
Age_specific: $(OBJECT_AGE)
	$(LINK) -o ../bin/CalcAge $(OBJECT_AGE)

all: $(OBJECTS) $(OBJECTS_TEST) $(OBJECTS_VTK) $(OBJECTS_PECUBE324) $(OBJECT_AGE) $(OBJECTS_MPI)
	make Pecube
	make Test
	make Vtk
	make Pecube324
	make Age_specific
	make MPI
	
clean:
	rm *.o
	rm *.exe
