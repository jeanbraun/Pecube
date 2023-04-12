.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .f95 .f77 /mpif90.e .r .y .yr .ye .l .p .sh .csh .h

# use these lines to produce a scalar version using the gfortran compiler
# this removes the mpi commands in NA

FLAGS = -c -O
FF90 = gfortran #-fbounds-check
FF77 = gfortran -fd-lines-as-comments #-fbounds-check
CC = gcc #-fbounds-check
LINK = gfortran -static-libgfortran -static-libgcc -static-libstdc++ #-fbounds-check 

# use these lines to produce a parallel (mpi) version to parallelize the inversion
# Note that some compilers will not allow to have an "include" command preceded by a "d" as an indication
#   of a commented line; if this is the case for your compiler, the routine na.f will not compile properly
#   You will need to edit the na.f file and remove the "d" in the first column of three lines that
#   start with include 'mpif.h'

FLAGSMPI = -c -O
FF90MPI = mpif90
FF77MPI = mpif77 -fd-lines-as-code
CCMPI = mpicc

#FLAGS = -c -O
#FF90 = ifort
#FF77 = ifort
#CC = icc

OBJECTS = \
module_Pecube.o \
TLModel.o \
OSLModel.o \
ESRModel.o \
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
scanfile.o \
read_data_files.o \
read_table.o \
ExtractDEM.o \
initialize_random_seed.o \
calculate_TL.o \
calculate_OSL.o \
calculate_ESR.o \
smooth_topo.o

OBJECTS_MPI = \
module_Pecube.o \
TLModel.o \
OSLModel.o \
ESRModel.o \
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
scanfile.o \
read_data_files.o \
read_table.o \
ExtractDEM.o \
initialize_random_seed.o \
calculate_TL.o\
calculate_OSL.o\
calculate_ESR.o \
smooth_topo.o

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
smooth_topo.o

OBJECTS_VTK = \
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
	mpif77 -c -O naMPI.f
	$(FF90MPI) -O $(OBJECTS_MPI) naMPI.o -o ../bin/PecubeMPI

Test: $(OBJECTS_TEST)
	$(LINK) -O $(OBJECTS_TEST) -o ../bin/Test

Vtk: $(OBJECTS_VTK)
	$(LINK) -O $(OBJECTS_VTK) -o ../bin/Vtk

Pecube324:	$(OBJECTS_PECUBE324)
	$(LINK) -O $(OBJECTS_PECUBE324) -o ../bin/Pecube324

all: $(OBJECTS) $(OBJECTS_TEST)  $(OBJECTS_VTK)
	make Pecube
	make Test
	make Vtk
	make Pecube324

clean:
	rm *.o
