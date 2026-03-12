!--------------------------------------------------------------------------------------------

module Pecube

! Key module for Pecube

!--------------------------------------------------------------------------------------------

! "version" type that contains information about Pecube version

  type version

  character*5 :: str = "4.3.6"
  integer :: major = 4
  integer :: minor = 3
  integer :: patch = 6

  end type version

!--------------------------------------------------------------------------------------------

! "param" type that contains all input parameters

  type parameters

  character*5 :: run_name = "RUN00"
  character*128 :: run_name_desc = "Run and corresponding folder name"

  !!!!!!!!!!!!!! Topography !!!!!!!!!!!!!!!!!
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
  
  integer :: ntime = 1
  character*128 :: ntime_desc = "Number of time steps used to describe the topographic evolution"
  
  double precision :: erosional_time_scale = 0.d0
  character*128 :: erosional_time_scale_desc = "Erosional time scale (in yr) characterizing transition between two time steps"
  
  double precision :: topo_wavelength = 0.d0
  character*128 :: topo_wavelength_desc = "Topography wavelength (km). Used to build synusoidal topography."
  
  double precision :: topo_offset = 0.d0
  character*128 :: topo_offset_desc = "Vertical offset for topography (km). Used to build synusoidal topography."
  
  double precision :: topo_amp = 1.d0
  character*128 :: topo_amp_desc = "Relief amplification for present-day synthetitc sinusoidal topography."
  
  double precision :: topo_phase = 0.d0
  character*128 :: topo_phase_desc = "Phase shift for present-day synthetitc sinusoidal topography (km)."

  double precision, dimension(:), allocatable :: time_topo
  character*128 :: time_topo_desc = "Topographic time stamp (in yr in past)"

  double precision, dimension(:), allocatable :: amplification
  character*128 :: amplification_desc = "Topographic amplification factor"

  double precision, dimension(:), allocatable :: offset
  character*128 :: offset_desc = "Topographic offset factor (in km)"

  integer, dimension(:), allocatable :: output
  character*128 :: output_desc = "Ouput flag (0/1)"
  
  double precision, dimension(:), allocatable :: phase
  character*128 :: phase_desc = "Topographic lateral phase shift (in km)"
  
  integer :: topo_ref = 0 !Maxime
  character*128 :: topo_ref_desc = "Flag to set the reference elevation for the computation of topography evolution."
  
  double precision :: topo_ref_custom = 0.d0
  character*128 :: topo_ref_custom_desc = "Custom reference elevation (m)."
  
  ! This is for simulating headward propagation of erosion (test)
  integer :: do_Headward = 0
  character*128 :: do_Headward_desc = "Flag to enable headward propagation of erosion."
  
  double precision :: tauH = 0.d0
  character*128 :: tauH_desc = "scaling factor for headward propagation rate."
  
  double precision :: depth_incision = 0.d0
  character*128 :: depth_incision_desc = "Relative depth of incision for wave propagation (%)."
  
  double precision :: time_incision_start = 0.d0
  character*128 :: time_incision_start_desc = "Onset time of headward propagation (Ma)."
  
  double precision :: time_incision_stop = 0.d0
  character*128 :: time_incision_stop_desc = "End time of headward propagation (Ma)"
  

  !!!!!!!!!!!!!! Thermal !!!!!!!!!!!!!!!!!!!!!!!!!
  double precision :: thickness = 35.d0
  character*128 :: thickness_desc = "Depth to bottom of the model (in km)"
  
  double precision :: thermal_diffusivity = 25.d0
  character*128 :: thermal_diffusivity_desc = "Thermal diffusivity (in m^2/yr)"

  double precision :: basal_temperature = 700.d0
  character*128 :: basal_temperature_desc = "Basal temperature (in C)"
  
  double precision :: D_heat = 20.d0
  character*128 :: D_heat_desc = "e-folding depth for exponential decrease in heat production."

  double precision :: sea_level_temperature = 0.d0
  character*128 :: sea_level_temperature_desc = "temperature at sea level (in C)"

  double precision :: lapse_rate = 0.d0
  character*128 :: lapse_rate_desc = "Atmospheric lapse rate (in C/m)"

  integer :: heatprod_flag = 0
  character*128 :: heatprod_flag_desc = "Flag for heat production: 0 = uniform HP, 1 = vertically non-uniform (e-folding)"

  double precision :: heat_production = 0.d0
  character*128 :: heat_production_desc = "Crustal heat production (in C/yr)"
  
  double precision :: min_dt = 1
  character*128 :: min_dt_desc = "Minimum time step to use in Pecube (Myr)"


  !!!!!!!!!!!!!!!!!! Isostasy !!!!!!!!!!!!!!!!!!!!!
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
  character*128 :: eet_desc = "Effective elastic plate thickness (in km)"

  integer :: nx_isostasy = 1024
  character*128 :: nx_isostasy_desc = "FFT resolution for flexural isostasy solution in x-direction"

  integer :: ny_isostasy = 1024
  character*128 :: ny_isostasy_desc = "FFT resolution for flexural isostasy solution in y-direction"

  
  !!!!!!!!!!!!!!!!!! Data/Thermochronometers !!!!!!!!!!!!!!!!!!!!!!
  character*128 :: data_folder = "Nil"
  character*128 :: data_folder_desc = "Data folder name, containing all ther,ochronological data"

  double precision :: default_age = -1.d0
  character*128 :: default_age_desc = "Default age for non reset rocks (in Ma)"

  integer :: FT_code_flag = 2
  character*128 :: FT_code_flag_desc = "Flag for choice of FT routine (0=vanderBeek/1=Ketcham (1999)/2=Ketcham (2007))"

  integer :: age_OSL_flag = 0
  character*128 :: age_OSL_flag_desc = "Flag for computing OSL ages (0=no/1=yes)"
  
  integer :: age_TL_flag = 0
  character*128 :: age_TL_flag_desc = "Flag for computing TL ages (0=no/1=yes)"

  integer :: age_ESR_flag = 0
  character*128 :: age_ESR_flag_desc = "Flag for computing ESR ages (0=no/1=yes)"

  integer :: age_AHe_flag = 1
  character*128 :: age_AHe_flag_desc = "Flag for computing helium ages in apatite (0=no/1=yes)"

  integer :: age_ZHe_flag = 0
  character*128 :: age_ZHe_flag_desc = "Flag for computing helium ages in zircon (0=no/1=yes)"

  integer :: age_AFT_flag = 0
  character*128 :: age_AFT_flag_desc = "Flag for computing FT ages in apatite (0=no/1=yes)"
  
  integer :: Init_FTL_Model = 0
  character*128 :: Init_FTL_Model_desc = "Calculate initial track length from  0 = Single value,&
                    & 1 = kinetic parameter."
  
  double precision :: Init_FTL_value = 16.3
  character*128 :: Init_FTL_value_desc = "Value for initial track length (µm)."
  
  integer :: Kinetic_FTL_Parameter = 4
  character*128 :: Kinetic_FTL_Parameter_desc = "Kinetic parameter for fission track annealing, 0 = Dpar, &
                    & 1= Cl (apfu), 2 = OH (apfu), 3 = Cl (wt%), 4 = RMR0."
  
  double precision :: Kinetic_FTL_Parameter_value_AFT = 0.79
  character*128 :: Kinetic_FTL_Parameter_value_AFT_desc = "Value for the kinetic parameter for fission track annealing."
  
  double precision :: MFTL_error = 0.4
  character*128 :: MFTL_error_desc = "Error on the MFTL (µm) to apply in the misfit calculation (inversion)."
  
  double precision :: MFTL_std_error = 0.3
  character*128 :: MFTL_std_error_desc = "Error on the std of MFTL (µm) to apply in the misfit calculation (inversion)."

  integer :: age_ZFT_flag = 0
  character*128 :: age_ZFT_flag_desc = "Flag for computing FT ages in zircon (0=no/1=yes)"

  integer :: age_FTL_flag = 0
  character*128 :: age_FTL_flag_desc = "Flag for computing FT lengths (0=no/1=MFTL/2=FTL distribution)"

  integer :: age_KAr_flag = 0
  character*128 :: age_KAr_flag_desc = "Flag for computing Argon ages in K-feldspar (0=no/1=yes)"

  integer :: age_BAr_flag = 0
  character*128 :: age_BAr_flag_desc = "Flag for computing Argon ages in biotite (0=no/1=yes)"

  integer :: age_MAr_flag = 0
  character*128 :: age_MAr_flag_desc = "Flag for computing Argon ages in muscovite (0=no/1=yes)"

  integer :: age_HAr_flag = 0
  character*128 :: age_HAr_flag_desc = "Flag for computing Argon ages in hornblend (0=no/1=yes)"
  
  integer :: age_computation = 0 !Maxime
  character*128 :: age_computation_desc = "Flag for computing ages as usual or sample-specific &
   &(0=no ages/1=for all nodes/2=sample-specific)"
  
  ! Parameters for AHe age computation - Maxime
  double precision :: rhoST = 0.893
  character*128 :: rhoST_desc = "Track density reduction in age standard."
  
  double precision :: grainradius = 60.d0
  character*128 :: grainradius_desc = "Sphere equivalent radius of the grain."
  
  double precision :: ASize = 60.d0
  character*128 :: ASize_desc = "Sphere equivalent radius of the apatite grain."
  
  ! Helium thermochronometry
  double precision :: AUPPM = 26.2d0
  character*128 :: AUppm_desc = "Uranium concentration in ppm."
  
  double precision :: AThPPM = 8.d0
  character*128 :: AThppm_desc = "Thorium concentration in ppm."
  
  integer :: RDmodel = 0
  character*128 :: RDmodel_desc = "Flag for the diffusion model for apatite."
  
  double precision :: D0_AHe = 50.d0
  character*128 :: D0_AHe_desc = "Apatite diffusivity for infinite temperature (cm²/s)."
  
  double precision :: Ea_AHe = 137.6536
  character*128 :: Ea_AHe_desc = "Apatite activation energy (kJ/mol)."

  double precision :: Kinetic_FTL_Parameter_value_AHe = 0.79
  character*128 :: Kinetic_FTL_Parameter_value_AHe_desc = "Value for the kinetic parameter for fission track annealing for&
                    &radiation damage model in AHe."
  
  double precision :: D0_ZHe = 0.46
  character*128 :: D0_ZHe_desc = "Zircon diffusivity for infinite temperature (cm²/s)."
  
  double precision :: Ea_ZHe = 169.0336
  character*128 :: Ea_ZHe_desc = "Zircon activation energy (kJ/mol)."
  
  integer :: RDmodelz = 5
  character*128 :: RDmodelz_desc = "Flag for the diffusion model for zircon."
  
  double precision :: ZUppm = 200.2d0
  character*128 :: ZUppm_desc = "Zircon uranium concentration in ppm."
  
  double precision :: ZThppm = 150.d0
  character*128 :: ZThppm_desc = "Zircon thorium concnetration in ppm."
  
  double precision :: ZSize = 60.d0
  character*128 :: ZSize_desc = "Sphere equivalent radius of the Zircon grain."
  
  integer :: Alpha_Ejec_Flag = 2
  character*128 :: Alpha_Ejec_Flag_desc = "Flag for alpha stopping distances."
  
  ! 4He/3He thermochronometer
  integer :: He43_flag = 0
  character*128 :: He43_flag_desc = "Flag to compute 4He/3He spectra (edge age)."
  
  integer :: n43step = 20
  character*128 :: n43step_desc = "Number of steps in 4He/3He heating schedule."
  
  integer :: He43_inv_flag = 0
  character*128 :: He43_inv_flag_desc = "Flag to choose what 4He/3He data to invert for."
  
  double precision, dimension(20) :: duration = (/0.38,0.38,0.51,0.66,0.66,0.46,0.45,0.48,&
                                          &0.66,0.53,0.48,0.5,0.56,0.63,0.5,0.5,0.5,0.5,0.5,0.5/)
  double precision, dimension(20) :: theating =(/200.0,270.0,290.0,300.0,310.0,330.0,340.0,&
                                          &350.0,350.0,370.0,400.0,410.0,420.0,440.0,460.0,475.0,&
                                          &500.0,600.0,700.0,900.0/)
                                          
  ! Feldspar K-Ar
  double precision :: D0_KAr = 0.014
  character*128 :: D0_KAr_desc = "Feldspar diffusivity for infinite temperature (cm²/s)."
  
  double precision :: Ea_KAr = 120.0
  character*128 :: Ea_KAr_desc = "Feldspar activation energy (kJ/mol)."
  
  ! Biotite K-Ar
  double precision :: D0_BAr = 0.4
  character*128 :: D0_BAr_desc = "Biotite diffusivity for infinite temperature (cm²/s)."
  
  double precision :: Ea_BAr = 211.09
  character*128 :: Ea_Bar_desc = "Biotite activation energy (kJ/mol)."
  
  ! Muscovite K-Ar
  double precision :: D0_MAr = 0.04
  character*128 :: D0_Mar_desc = "Muscovite diffusivity for infinite temperature (cm²/s)."
  
  double precision :: Ea_MAr = 217.36
  character*128 :: Ea_MAr_desc = "Muscovite activation energy (kJ/mol)."
  
  ! Hornblende K-Ar
  double precision :: D0_HAr = 0.06
  character*128 :: D0_Har_desc = "Hornblende diffusivity for infinite temperature (cm²/s)."
  
  double precision :: Ea_HAr = 276.0
  character*128 :: Ea_HAr_desc = "Hornblende activation energy (kJ/mol)."
  
  
  ! Luminescence
  double precision :: TL_doser = 5.d0
  character*128 :: TL_doser_desc = "Dose rate (Gy/Kyr)"

  double precision :: TL_D0 = 800.d0
  character*128 :: TL_D0_desc = "Characteristic dose (Gy)"

  double precision :: TL_a = 1.d0
  character*128 :: TL_a_desc = "Kinetic orders of trapping"

  double precision :: TL_b = 1.d0
  character*128 :: TL_b_desc = "Kinetic orders of detrapping"

  double precision :: TL_Et = 1.4d0
  character*128 :: TL_Et_desc = "Activation energy (eV)"

  double precision :: TL_logs = 12.d0
  character*128 :: TL_logs_desc = "Logarithm of thermal frequency factor (1/s)"

  double precision :: TL_logrho = -5.5d0
  character*128 :: TL_logrho_desc = "Logarithm of dimensionless recombination center density"
  
  integer :: TL_Model = 0
  character*128 :: TL_Model_desc = "Flag to choose the TL model."

   double precision :: OSL_doser = 5.d0
  character*128 :: OSL_doser_desc = "Dose rate (Gy/Kyr)"

  double precision :: OSL_D0 = 800.d0
  character*128 :: OSL_D0_desc = "Characteristic dose (Gy)"

  double precision :: OSL_Et = 1.4d0
  character*128 :: OSL_Et_desc = "Activation energy"

  double precision :: OSL_Eu = 12.d0
  character*128 :: OSL_Eu_desc = "Logarithm of thermal frequency factor"

  double precision :: OSL_logs = 12.d0
  character*128 :: OSL_logs_desc = "Logarithm of thermal frequency factor"

  double precision :: OSL_logrho = -5.5d0
  character*128 :: OSL_logrho_desc = "Logarithm of dimensionless recombination center density"
  
  double precision :: OSL_Lmax = 1.d0
  character*128 :: OSL_Lmax_desc = "Maximum luminescence signal (Lx/Tx)."
  
  integer :: OSL_Model = 0
  character*128 :: OSL_Model_desc = "Flag to choose the OSL model. 0: SSE-BTS-FAD, &
        &1: GOK-Gauss-FAD, 2: GOK-FAD."
  
  double precision :: OSL_a = 1.0
  character*128 :: OSL_a_desc = "Kinetic order for electron trapping."

  double precision :: OSL_b = 1.0
  character*128 :: OSL_b_desc = "Kinetic order for electron detrapping."

  double precision :: ESR_doser = 4.265d0
  character*128 :: ESR_doser_desc = "Dose rate"

  double precision :: ESR_D0 = 5186.8d0
  character*128 :: ESR_D0_desc = "Onset of dose saturation"

  double precision :: ESR_logs = 14.58d0
  character*128 :: ESR_logs_desc = "Logarithm of frequency factor"

  double precision :: ESR_Et = 1.763d0
  character*128 :: ESR_Et_desc = "Center of Gaussian distribution of activation energies"

  double precision :: ESR_sigmaEt = 0.096d0
  character*128 :: ESR_sigmaEt_desc = "Gaussian width of activation energy distribution"
  
  double precision :: ESR_a = 1.d0
  character*128 :: ESR_a_desc = "Kinetic orders of trapping for ESR"
  
  double precision :: ESR_b = 1.d0
  character*128 :: ESR_b_desc = "Kinetic orders of detrapping for ESR"
  
  double precision :: ESR_Lmax = 1.d0
  character*128 :: ESR_Lmax_desc = "Maximum ESR signal (Lx/Tx)."
  
  integer :: ESR_Model = 0
  character*128 :: ESR_Model_desc = "Flag to choose the ESR model."
  
                                  
  

  !!!!!!!!!!!!!!!!!!!!!! Tectonic !!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: fault_advect_flag = 0
  character*128 :: fault_advect_flag_desc = "Flag for choice of fault advection (0=off/1=on)"

  double precision :: shear_heating = 0.d0
  character*128 ::  shear_heating_desc = "Shear heating friction coefficient"

  integer :: nfault = 0
  character*128 :: nfault_desc = "Number of faults"

  double precision :: x1 = 0.d0
  character*128 :: x1_desc = "Longitude of first point defining fault trace (in degrees longitude)"

  double precision :: y1 = 0.d0
  character*128 :: y1_desc = "Latitude of first point defining fault trace (in degrees latitude)"

  double precision :: x2 = 0.d0
  character*128 :: x2_desc = "Longitude of second point defining fault trace (in degrees longitude)"

  double precision :: y2 = 0.d0
  character*128 :: y2_desc = "Latitude of second point defining fault trace (in degrees latitude)"

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
  character*128 :: r_desc = "r-coordinates of points used to describe each fault geometry (in km)"

  double precision, dimension(:,:), allocatable :: s
  character*128 :: s_desc = "s-coordinates of points used to describe each fault geometry (in km)"

  integer, dimension(:), allocatable :: nstep
  character*128 :: nstep_desc = "Number of time steps used to describe the motion on each fault"

  integer, dimension(:), allocatable :: static
  character*128 :: static_desc = "Flag to activate velocity field (static=0) or not (static=1)" 
  
  double precision, dimension(:,:), allocatable :: time_start
  character*128 :: time_start_desc = "Starts of time steps used to describe the motion on each fault (in Ma)"

  double precision, dimension(:,:), allocatable :: time_end
  character*128 :: time_end_desc = "Ends of time steps used to describe the motion on each fault (in Ma)"

  double precision, dimension(:,:), allocatable :: velo
  character*128 :: velo_desc = "Velocities used to describe the motion on each fault (in km/Myr)"

  integer :: logarithmic_velocity = 0
  character*128 :: logarithmic_velocity_desc = "Flag to enable logarithmic velocity specification (mostly for inversion purposes)"


  !!!!!!!!!!!!!!!!!!  Outputs !!!!!!!!!!!!!!!!!!!!!!!!
  integer :: debug = 0
  character*128 :: debug_desc = "Flag to enable debugging mode"

  integer :: save_PTT_paths = 0
  character*128 :: save_PTT_paths_desc = "Flag to save PTT paths at observation points"
  integer :: save_ages_inversion = 0
  character*128 :: save_ages_inversion_desc = "Flag to save predicted ages for each model runs during NA inversion"                         
                                                                                                                   
  integer :: save_cooling_rates = 0 !Maxime
  character*128 :: save_cooling_rates_desc = "Flag to save PTT paths for all nodes"
  
  integer :: save_eroded_volume = 0
  character*128 :: save_eroded_volume_desc = "Flag to save eroded volume history"

  integer :: echo_input_file =0
  character*128 :: echo_input_file_desc = "Flag to echo the reading of the input file to the screen/terminal"

  
  !!!!!!!!!!!!!!!!!!!!!!!! Inversion !!!!!!!!!!!!!!!!!!!!!!!
  double precision :: error_predictions = 0.d0
  character*128 :: error_predictions_desc = "Consider error on age predictions (%)."
  
  double precision :: misfit_weight_AHE = 1.d0
  character*128 :: misfit_weight_AHE_desc = "Weight that multiplies the AHE part of the misfit"
  
  double precision :: misfit_weight_ZHE = 1.d0
  character*128 :: misfit_weight_ZHE_desc = "Weight that multiplies the ZHE part of the misfit"
  
  double precision :: misfit_weight_AFT = 1.d0
  character*128 :: misfit_weight_AFT_desc = "Weight that multiplies the AFT part of the misfit"
  
  double precision :: misfit_weight_ZFT = 1.d0
  character*128 :: misfit_weight_ZFT_desc = "Weight that multiplies the ZFT part of the misfit"
  
  double precision :: misfit_weight_BAR = 1.d0
  character*128 :: misfit_weight_BAR_desc = "Weight that multiplies the BAR part of the misfit"
  
  double precision :: misfit_weight_MAR = 1.d0
  character*128 :: misfit_weight_MAR_desc = "Weight that multiplies the MAR part of the misfit"
  
  double precision :: misfit_weight_HAR = 1.d0
  character*128 :: misfit_weight_HAR_desc = "Weight that multiplies the HAR part of the misfit"
  
  double precision :: misfit_weight_KAR = 1.d0
  character*128 :: misfit_weight_KAR_desc = "Weight that multiplies the KAR part of the misfit"

  double precision :: misfit_weight_FTLD = 1.d0
  character*128 :: misfit_weight_FTLD_desc = "Weight that multiplies the Fission Track Length Distribution part of the misfit"

  double precision :: misfit_weight_TH = 1.d0
  character*128 :: misfit_weight_TH_desc = "Weight that multiplies the thermal histories part of the misfit"

  double precision :: misfit_weight_43He = 1.d0
  character*128 :: misfit_weight_43He_desc = "Weight that multiplies the 43Helium part of the misfit"

  double precision :: misfit_weight_TL = 1.d0
  character*128 :: misfit_weight_TL_desc = "Weight that multiplies the Thermoluminescence part of the misfit"

  double precision :: misfit_weight_OSL = 1.d0
  character*128 :: misfit_weight_OSL_desc = "Weight that multiplies the OSL part of the misfit"

  double precision :: misfit_weight_ESR = 1.d0
  character*128 :: misfit_weight_ESR_desc = "Weight that multiplies the ESR part of the misfit"

  integer :: ESR_misfit_Target = 0
  character*128 :: ESR_misfit_Target_desc = "Flag to calculate misfit on ESR age or n/N values. (0= Age, 1= n/N)"
  
  integer :: OSL_misfit_Target = 0
  character*128 :: OSL_misfit_Target_desc = "Flag to calculate misfit on OSL age or n/N values. (0= Age, 1= n/N)"

  integer :: misfit_slope = 0
  character*128 :: misfit_slope_desc = "Flag for choice of misfit function (0=ages/1=slopes)"
  
  integer :: misfit_function = 1
  character*145 :: misfit_function_desc = "Flag for choice of misfit function&
       &(1:Chi-squared,2: Reduced chi-squared, 3: L2-norm, 4:L1-norm,&
       &5: Reduced L1-norm, 6: Log-scale misfit)"

  integer :: maximum_number_of_iterations = 4
  character*128 :: maximum_number_of_iterations_desc = "Maximum number of NA iterations "//&
          &"(not including the first one)"

  integer :: sample_size_for_first_iteration = 8
  character*128 :: sample_size_for_first_iteration_desc = "Number of model runs (samples) to be performed"//&
          &" by NA during the first iteration"

  integer :: sample_size_for_all_other_iterations = 8
  character*128 :: sample_size_for_all_other_iterations_desc = "Number of model runs (samples) "//&
          &"to be performed by NA during all subsequent iterations"

  integer :: number_of_cells_to_resample = 4
  character*128 :: number_of_cells_to_resample_desc = "Number of Voronoi cells to be resampled "//&
          &"by NA at each iteration"
          
  integer :: score_with_data = 0
  character*128 :: score_with_data_desc = "Score model by the number of data that fit within uncertainty"
  
  integer :: inversion_mode = 0
  character*128 :: inversion_mode_desc = "Select the inversion mode: 0 = NA, 1= NA from previous na.nad file, 2 = MC (batch)."

  
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
