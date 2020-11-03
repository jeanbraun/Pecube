!--------------------------------------------------------------------------------------------

module Pecube

! Key module for Pecube

!--------------------------------------------------------------------------------------------

! "version" type that contains information about Pecube version

  type version

  character*5 :: str = "4.2.1"
  integer :: major = 4
  integer :: minor = 2
  integer :: patch = 0

  end type version

!--------------------------------------------------------------------------------------------

! "param" type that contains all input parameters

  type parameters

  character*5 :: run_name = "RUN00"
  character*128 :: run_name_desc = "Run and corresponding folder name"

  character*128 :: topo_file_name = "Nil"
  character*128 :: topo_file_name_desc = "DEM file name (stored in RUN)"

  integer :: nx = 31
  character*128 :: nx_desc = "Grid resolution in the x-direction"

  integer :: ny = 31
  character*128 :: ny_desc = "Grid resolution in the y-direction"

  integer :: nz = 21
  character*128 :: nz_desc = "Grid resolution in the z-direction"

  double precision :: lon0 = 0.
  character*128 :: lon0_desc = "x-axis origin (in degrees longitude)"

  double precision :: lat0 = 0.
  character*128 :: lat0_desc = "y-axis origin (in degrees latitude)"

  double precision :: dlon = 0.0083333
  character*128 :: dlon_desc = "Grid spacing in x-direction (in degrees longitude)"

  double precision :: dlat = 0.0083333
  character*128 :: dlat_desc = "Grid spacing in y-direction (in degrees latitude)"

  integer :: nskip = 1
  character*128 :: nskip_desc = "Skipping factor (integer >= 1)"

  double precision :: thickness = 35.d0
  character*128 :: thickness_desc = "Depth to bottom of the model (in m)"

  integer :: ntime = 1
  character*128 :: ntime_desc = "Number of time steps used to describe the topographic evolution"

  double precision :: erosional_time_scale = 0.d0
  character*128 :: erosional_time_scale_desc = "Erosional time scale (in yr) characterizing transition between two time steps"

  double precision, dimension(:), allocatable :: time_topo
  character*128 :: time_topo_desc = "Topographic time stamp (in yr in past)"

  double precision, dimension(:), allocatable :: amplification
  character*128 :: amplification_desc = "Topographic amplification factor"

  double precision, dimension(:), allocatable :: offset
  character*128 :: offset_desc = "Topographic offset factor (in m)"

  integer, dimension(:), allocatable :: output
  character*128 :: output_desc = "Ouput flag (0/1)"

  integer :: isostasy = 0
  character*128 :: isostasy_desc = "Isostatic flag (0/1)"

  double precision :: rho_crust = 2400.d0
  character*128 :: rho_crust_desc = "Crustal rock density (in kg/m^3)"

  double precision :: rho_asthenosphere = 3150.d0
  character*128 :: rho_asthenosphere_desc = "Asthenospheric rock density (in kg/m^3)"

  double precision :: young_modulus = 1.d11
  character*128 :: young_modulus_desc = "Young's modulus (in Pa)"

  double precision :: poisson_ratio = 0.25d0
  character*128 :: poisson_ratio_desc = "Poisson's ratio"

  double precision :: EET = 20.d0
  character*128 :: eet_desc = "Effective elastic plate thickness (in m)"

  integer :: nx_isostasy = 1024
  character*128 :: nx_isostasy_desc = "FFT resolution for flexural isostasy solution in x-direction"

  integer :: ny_isostasy = 1024
  character*128 :: ny_isostasy_desc = "FFT resolution for flexural isostasy solution in y-direction"

  double precision :: thermal_diffusivity = 25.d0
  character*128 :: thermal_diffusivity_desc = "Thermal diffusivity (in m^2/yr)"

  double precision :: basal_temperature = 700.d0
  character*128 :: basal_temperature_desc = "Basal temperature (in C)"

  double precision :: sea_level_temperature = 0.d0
  character*128 :: sea_level_temperature_desc = "temperature at sea level (in C)"

  double precision :: lapse_rate = 0.d0
  character*128 :: lapse_rate_desc = "Atmospheric lapse rate (in C/m)"

  double precision :: heat_production = 0.d0
  character*128 :: heat_production_desc = "Crustal heat production (in C/yr)"

  character*128 :: data_folder = "Nil"
  character*128 :: data_folder_desc = "Data folder name, containing all ther,ochronological data"

  double precision :: default_age = -1.d0
  character*128 :: default_age_desc = "Default age for non reset rocks (in yr)"

  integer :: FT_code_flag = 0
  character*128 :: FT_code_flag_desc = "Flag for choice of FT routine (0=vanderBeek/1=Ketcham)"

  integer :: fault_advect_flag = 0
  character*128 :: fault_advect_flag_desc = "Flag for choice of fault advection (0=off/1=on)"

  double precision :: shear_heating = 0.d0
  character*128 ::  shear_heating_desc = "Shear heating friction coefficient"

  integer :: age_TL_flag = 0
  character*128 :: age_TL_flag_desc = "Flag for computing TL ages (0=no/1=yes)"

  integer :: age_AHe_flag = 1
  character*128 :: age_AHe_flag_desc = "Flag for computing helium ages in apatite (0=no/1=yes)"

  integer :: age_ZHe_flag = 0
  character*128 :: age_ZHe_flag_desc = "Flag for computing helium ages in zircon (0=no/1=yes)"

  integer :: age_AFT_flag = 0
  character*128 :: age_AFT_flag_desc = "Flag for computing FT ages in apatite (0=no/1=yes)"

  integer :: age_ZFT_flag = 0
  character*128 :: age_ZFT_flag_desc = "Flag for computing FT ages in zircon (0=no/1=yes)"

  integer :: age_FTL_flag = 0
  character*128 :: age_FTL_flag_desc = "Flag for computing FT lengths (0=no/1=yes)"

  integer :: age_KAr_flag = 0
  character*128 :: age_KAr_flag_desc = "Flag for computing Argon ages in K-feldspar (0=no/1=yes)"

  integer :: age_BAr_flag = 0
  character*128 :: age_BAr_flag_desc = "Flag for computing Argon ages in biotite (0=no/1=yes)"

  integer :: age_MAr_flag = 0
  character*128 :: age_MAr_flag_desc = "Flag for computing Argon ages in muscovite (0=no/1=yes)"

  integer :: age_HAr_flag = 0
  character*128 :: age_HAr_flag_desc = "Flag for computing Argon ages in hornblend (0=no/1=yes)"

  integer :: nfault = 0
  character*128 :: nfault_desc = "Number of faults"

  double precision :: x1 = 0.d0
  character*128 :: x1_desc = "Longitude of first fault defining fault trace (in degrees longitude)"

  double precision :: y1 = 0.d0
  character*128 :: y1_desc = "Latitude of first fault defining fault trace (in degrees latitude)"

  double precision :: x2 = 0.d0
  character*128 :: x2_desc = "Longitude of second fault defining fault trace (in degrees longitude)"

  double precision :: y2 = 0.d0
  character*128 :: y2_desc = "Latitude of second fault defining fault trace (in degrees latitude)"

  double precision :: bottom_left = 1.d0
  character*128 :: bottom_left_desc = "Scaling value for uplift function applied at bottom left corner of grid"

  double precision :: bottom_right = 1.d0
  character*128 :: bottom_right_desc = "Scaling value for uplift function applied at bottom right corner of grid"

  double precision :: upper_right = 1.d0
  character*128 :: upper_right_desc = "Scaling value for uplift function applied at top right corner of grid"

  double precision :: upper_left = 1.d0
  character*128 :: upper_left_desc = "Scaling value for uplift function applied at top left corner of grid"

  integer, dimension(:), allocatable :: npoint
  character*128 :: npoint_desc = "number of points used to describe each fault geometry"

  double precision, dimension(:,:), allocatable :: r
  character*128 :: r_desc = "r-coordinates of points used to describe each fault geometry (in m)"

  double precision, dimension(:,:), allocatable :: s
  character*128 :: s_desc = "r-coordinates of points used to describe each fault geometry (in m)"

  integer, dimension(:), allocatable :: nstep
  character*128 :: nstep_desc = "Number of time steps used to describe the motion on each fault"

  integer, dimension(:), allocatable :: static
  character*128 :: static_desc = "Flag to activate velocity field (static=0) or not (static=1)"

  double precision, dimension(:,:), allocatable :: time_start
  character*128 :: time_start_desc = "Starts of time steps used to describe the motion on each fault (in yr)"

  double precision, dimension(:,:), allocatable :: time_end
  character*128 :: time_end_desc = "Ends of time steps used to describe the motion on each fault (in yr)"

  double precision, dimension(:,:), allocatable :: velo
  character*128 :: velo_desc = "Velocities used to describe the motion on each fault (in m/yr)"

  integer :: logarithmic_velocity = 0
  character*128 :: logarithmic_velocity_desc = "Flag to enable logarithmic velocity specification (mostly for inversion purposes)"

  integer :: debug = 0
  character*128 :: debug_desc = "Flag to enable debugging mode"

  integer :: save_PTT_paths = 0
  character*128 :: save_PTT_paths_desc = "Flag to save PTT paths at observation points"

  integer :: save_ages_inversion = 0
  character*128 :: save_ages_inversion_desc = "Flag to save predicted ages for each model runs during NA inversion"

  integer :: save_eroded_volume = 0
  character*128 :: save_eroded_volume_desc = "Flag to save eroded volume history"

  integer :: echo_input_file =0
  character*128 :: echo_input_file_desc = "Flag to echo the reading of the input file to the screen/terminal"

  double precision :: misfit_weight_AGE = 1.d0
  character*128 :: misfit_weight_AGE_desc = "Weight that multiplies the age part of the misfit"

  double precision :: misfit_weight_FTLD = 1.d0
  character*128 :: misfit_weight_FTLD_desc = "Weight that multiplies the Fission Track Length Distribution part of the misfit"

  double precision :: misfit_weight_TH = 1.d0
  character*128 :: misfit_weight_TH_desc = "Weight that multiplies the thermal hostories part of the misfit"

  double precision :: misfit_weight_43He = 1.d0
  character*128 :: misfit_weight_43He_desc = "Weight that multiplies the 43Helium part of the misfit"

  double precision :: misfit_weight_TL = 1.d0
  character*128 :: misfit_weight_TL_desc = "Weight that multiplies the Thermoluminescence part of the misfit"

  integer :: misfit_corrected = 0
  character*128 :: misfit_corrected_desc = "Flag to divide the misfit function by (N-nd-1) if set to 1"

  integer :: misfit_slope = 0
  character*128 :: misfit_slope_desc = "Flag for choice of misfit function (0=ages/1=slopes)"

  integer :: maximum_number_of_iterations = 4
  character*128 :: maximum_number_of_iterations_desc = "Maximum number of NA iterations "//&
          "(not including the first one)"

  integer :: sample_size_for_first_iteration = 8
  character*128 :: sample_size_for_first_iteration_desc = "Number of model runs (samples) to be performed"//&
          " by NA during the first iteration"

  integer :: sample_size_for_all_other_iterations = 8
  character*128 :: sample_size_for_all_other_iterations_desc = "Number of model runs (samples) "//&
          "to be performed by NA during all subsequent iterations"

  integer :: number_of_cells_to_resample = 4
  character*128 :: number_of_cells_to_resample_desc = "Number of Voronoi cells to be resampled "//&
          "by NA at each iteration"

  double precision :: TL_doser = 5.d0
  character*128 :: TL_doser_desc = "Dose rate"

  double precision :: TL_D0 = 800.d0
  character*128 :: TL_D0_desc = "Onset of dose saturation"

  double precision :: TL_a = 1.8d0
  character*128 :: TL_a_desc = "Kinetic orders of trapping"

  double precision :: TL_b = 1.8d0
  character*128 :: TL_b_desc = "Kinetic orders of detrapping"

  double precision :: TL_Et = 1.4d0
  character*128 :: TL_Et_desc = "Activation energy"

  double precision :: TL_logs = 12.d0
  character*128 :: TL_logs_desc = "Logarithm of thermal frequency factor"

  double precision :: TL_logrho = -5.5d0
  character*128 :: TL_logrho_desc = "Logarithm of dimensionless recombination center density"

  end type parameters

!--------------------------------------------------------------------------------------------

! "faulttype" type that is used to define faults

  type faulttype

! A fault vertical cross-sectional (2D) geometry is defined by a set of points (x,y)
! Its 3D location is defined by two points in the horizontal plane (x1,y1),(x2,y2)
! Those points limit its lateral extent
! In the (x,y) coordinate system (in which the faul vertical geometry is defined)
! x lies to the right of the (x1,y1),(x2,y2) line and z is vertical upwards
! and zero at the surface

! n is the number of points used to define the fault
! x(n) and y(n) are the coordinates of the points
! (x1,y1) and (x2,y2) are the coordinates of the trace of the fault at the surface
! nstep is the number of time intevals used to define the motion story of the fault
! per interval: timestart,timeend and velo
! timestart and timeend are the start and end time of the interval
! velo is the velocity across the fault (along the plane of the fault) for the interval

! By convention as one goes along the fault from point 1 to point n
! it is the block to the right that moves

! The sign of velo determines whether the sign of the x-component of velocity

! from that information we compute xs(n-1),ys(n-1), the direction
! of each of the (n-1) segments connecting the n points

! as well as xn,yn the normal to the (x1,y1)-(x2,y2) line in the z=0 plane

  integer n,nstep,static
  double precision,dimension(:),pointer::x,y
  double precision x1,y1,x2,y2
  double precision,dimension(:),pointer::timestart,timeend,velo
  double precision,dimension(:),pointer::xs,ys
  double precision xn,yn

  end type faulttype

!--------------------------------------------------------------------------------------------

!type "edge"

! this type is to store edges in a trianglulation
! it is used to update (in a generalized Delaunay sense)
! the triangulation of the 3D points on the surfaces
! for each edge:
! n1, n2 are the node numbers defining the edge
! t1, t2 are the triangle numbers on either side of the edge
! going from n1 to n2, t1 is to the left and t2 is to the right
! m1, m2 are the node numbers of the two other nodes making t1 and t2

  type edge
  integer n1,n2,m1,m2,t1,t2
  end type edge

!--------------------------------------------------------------------------------------------

! following is a general interface to read stuff from a file

  interface scanfile

    subroutine iscanfile (fnme,text,res,res_desc,ires,vocal,nd,range,par)
    character*(*) fnme,text,res_desc
    integer,intent(out)::res
    integer,intent(out)::ires
    integer,intent(in)::vocal
    integer,intent(inout)::nd
    real*4,intent(inout)::range(2,1024),par(1024)
    end subroutine iscanfile

    subroutine dscanfile (fnme,text,res,res_desc,ires,vocal,nd,range,par)
    character*(*) fnme,text,res_desc
    double precision,intent(out)::res
    integer,intent(out)::ires
    integer,intent(in)::vocal
    integer,intent(inout)::nd
    real*4,intent(inout)::range(2,1024),par(1024)
    end subroutine dscanfile

    subroutine cscanfile (fnme,text,res,res_desc,ires,vocal,nd,range,par)
    character*(*) fnme,text,res_desc
    character*(*),intent(out)::res
    integer,intent(out)::ires
    integer,intent(in)::vocal
    integer,intent(inout)::nd
    real*4,intent(inout)::range(2,1024),par(1024)
    end subroutine cscanfile

  end interface scanfile

!--------------------------------------------------------------------------------------------

! following is the interface for the ucase function so that one does not need to assume a length for its
! input/output

  interface ucase

    function ucase(in) result (out)

! function that forces a string to be in upper case

    implicit none

    character (*), intent(in)  :: in
    character(:), allocatable  :: out

    end function ucase

  end interface ucase

!--------------------------------------------------------------------------------------------

end module Pecube

!--------------------------------------------------------------------------------------------

module read_data_module

!--------------------------------------------------------------------------------------------

interface read_data_files

subroutine read_data_files (fnme,t,nx,ny,f,s)

character*(*) :: fnme
character*100, dimension(:,:), allocatable :: t
character*100, dimension(:), allocatable :: f,s
integer :: nx,ny

end subroutine read_data_files

end interface

!--------------------------------------------------------------------------------------------

interface ucase

function ucase(in) result (out)

character (*), intent(in)  :: in
character(:), allocatable  :: out

end function ucase

end interface

!--------------------------------------------------------------------------------------------

end module read_data_module

!--------------------------------------------------------------------------------------------

module read_string_module

!--------------------------------------------------------------------------------------------

interface read_string

!--------------------------------------------------------------------------------------------

function read_string (unit, istring, jstring) result (out)

integer, intent(in)  :: unit, istring, jstring
character(:), allocatable  :: out

end function read_string

end interface

!--------------------------------------------------------------------------------------------

end module read_string_module

!--------------------------------------------------------------------------------------------

module DEM

!--------------------------------------------------------------------------------------------

interface ExtractDEM

subroutine ExtractDEM (lonmin, latmin, nx, ny, demout, PecubeFnme, DEMFnme)

double precision, intent(in) :: lonmin, latmin
integer, intent(in) :: nx, ny
double precision, dimension(:,:), pointer, optional, intent(out) :: demout
character*(*), optional, intent(in) :: PecubeFnme, DEMFnme

end subroutine ExtractDEM

end interface ExtractDEM

!--------------------------------------------------------------------------------------------

end module DEM

!--------------------------------------------------------------------------------------------
