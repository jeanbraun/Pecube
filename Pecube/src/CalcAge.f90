
program calculate_ages_specific
!-------------------------------------------------------------------------------
! % This script computes (U-Th)/He ages and FT on apatite and zircon, and simulates degassing experiments to predict
! % 4He/3He profiles. It can be called after a Pecube run to calculate ages, if the thermal paths
! % are unchanged. The routine to calculate thermochronometer ages are shared with Pecube.
! %
! % It includes alpha ejection models and alpha stopping distances from:
! %  Farley et al. (1996)
! %  Ketcham et al. (2011)
! % 
! % It also includes radiation damage models from
! %  Flowers et al. (2009)
! %  Gautheron et al. (2009)
! %  Willett et al. (2017)
! %
! %

! % Last updates: 12/03/2026
! % Author: Maxime Bernard

!################### Initiate variables ####################

implicit none

character*5 run
integer :: io,numarg,ngrains,i,nsamples,TotGrains,j,nstep,istep,ntemp,rr
character*9 arg
double precision,dimension(:),allocatable::time
double precision time_step_array(1000),temp_step_array(1000)
integer, dimension(:), allocatable:: ngrains_array
! 4He/3He thermochronometer parameters
double precision,dimension(:),allocatable::He_ratio,TotHe3_released,He_ratioEjec
double precision,dimension(:,:),allocatable::TotHe3_released_array,He_ratioEjec_array,He_ratio_array
double precision,dimension(:,:),allocatable::temperature_samples,Heating_schedule,temperature_samples_temp
double precision, dimension(:,:), allocatable :: Step_duration,Step_temperature,RsRb,dRsRb,Sum3He,dSum3He 
integer :: nstep_He43,nsteps_He43_array(50),col_Heating,row_Heating,He43_flag, nb_He43, nb_step_max
integer :: khe, cnt_43
double precision,dimension(:),allocatable:: He43obs,SF3Heobs
real, dimension(:), allocatable :: UPPM43,THPPM43,RMR043,age43,dage43,size43, age43_obs
character(len=20), dimension(:), allocatable :: sample_list_He43
! Other parameters
integer :: ntime,kg,ks,cnt,ntime_temp,Graincounter
integer :: AHe_flag, AFT_flag, ZHe_flag,KAr_flag,BAr_flag,MAr_flag,HAr_flag,nb_ADAM_data
integer :: THL_flag, OSL_flag, ESR_flag
real,dimension(:,:),allocatable::ageAHe_array,ageAFT_array,ageZHe_array
real,dimension(:,:),allocatable::ageKAr_array,ageBAr_array,ageMAr_array,ageHAr_array
real,dimension(:,:),allocatable::MFTL_array, DMFTL_array
real,dimension(:),allocatable::ageTime
double precision,dimension(:,:),allocatable::rho_r
double precision data_ADAM_array(3,1000000)
real *4 Ug, Thg
real *4 gsize,grainsize,dt0
real *8 D0,Ea
real *8 D0z,Eaz !Zircon
real *8 D0k,Eak !Feldspar
real *8 D0b,Eab !Biotite
real *8 D0m,Eam !Muscovite
real *8 D0h,Eah !Hornblende
character*4 c(9999)
character*20 sample_name(300)
character*3 grain_ID
character*12 cdummy
integer sample_temp_flag, irec
character*20 dataFolder
real,dimension(:),allocatable::AHEOBS,THLOBS,OSLOBS,ESROBS
real,dimension(:),allocatable::AFTOBS,ZHEOBS,ZFTOBS
real,dimension(:),allocatable::KAROBS,BAROBS,MAROBS
real,dimension(:),allocatable::HAROBS,FT01OBS,FT02OBS
real,dimension(:),allocatable::FT03OBS,FT04OBS,FT05OBS
real,dimension(:),allocatable::FT06OBS,FT07OBS,FT08OBS
real,dimension(:),allocatable::FT09OBS,FT10OBS,FT11OBS,FT12OBS
real,dimension(:),allocatable::FT13OBS,FT14OBS,FT15OBS,FT16OBS
real,dimension(:),allocatable::FT17OBS,FT18OBS,FT19OBS,FT20OBS
real,dimension(:),allocatable::gsize_array,gsizez_array,Ug_array,Thg_array,Ugz_array,Thgz_array
double precision, dimension(:), allocatable:: KINAFT, KINAHE
integer,dimension(:),allocatable::AFT_Obs, AHE_Obs, ZHE_Obs, ThL_Obs, OSL_Obs, ESR_Obs
integer,dimension(:),allocatable::KAR_Obs, BAR_Obs, MAR_Obs, HAR_Obs
real,dimension(:),allocatable::rho_red
double precision, dimension(:), allocatable:: dtf
double precision,dimension(:),allocatable::time_interp,temp_interp
! AFT
double precision oldest_age,final_age, fmeanp,dfmeanp,Init_FTL_value
double precision fdist(200)
integer Init_FTL_Model, kinFTLID
double precision ts, rhoST
logical is_unix

!##################### Input Parameters ####################
! The input parameters are read from the samples_settings.txt file created by 
! PecubeGUI
 
! Alpha_Flag = To account for alpha ejection (0:no computation, 1: Ejection, 2:Redistribution (not yet available))
! Alpha_Ejec_Flag = alpha ejection stopping distances to consider (0:no computation, 1:Farley, 2: Ketcham)
! RDmodel = Radiation Damage model (0:no RD model, 1:Flowers 2009, 2: Willett 2017) 
! ztime = time vector of a thermal history (Ma)
! ztemp = temperature vector of a thermal history (°C)
! ztime_heating = duration vector of a heating schedule (hours) - 4He/3He 
! ztemp_heating = temperature vector of a heating schedule (°C) - 4He/3He 
 
! Flags
integer Alpha_flag
integer RDmodel, Alpha_Ejec_flag, dummy
integer DiffzModel,AnnModel
real dummyReal
 
! Thermal Scenarios
double precision,dimension(:),allocatable::ztime_double,ztemp_double
real, dimension(:),allocatable:: ztemp_real,ztime_real
double precision,dimension(:),allocatable:: ztime_heating,ztemp_heating
real*4,dimension(:),allocatable::ztime,ztemp


! Ages
real age
 
write (*,*) '-------------------------------------------------------------'
write (*,*) '--------------- Computing Samples specific ------------------'
write (*,*) '-------------------------------------------------------------'

! Get the folder name (Pecube project)
numarg=command_argument_count()
if (numarg.eq.0) then
write (*,*) 'You need to specify a run directory (i.e. RUN00 for example)'
stop 'End of run'
else
   call getarg (1,arg)
   run = arg
endif

! Initiate writing file number 
do i = 1, 9999
write (c(i),'(i4)') i
if (i.lt.10) c(i)(1:3)='000'
if (i.lt.100) c(i)(1:2)='00'
if (i.lt.1000) c(i)(1:1)='0'
enddo
      
! Check Os system
call GetOperatingSystem(is_unix)
!---------------------------------------------
!--------- Read samples specific file --------
!---------------------------------------------
! To retrieve the flags and parameters for age computation
! ngrains_array = array with the number of grain per sample
! gsize_array = grain size for each grain (µm)
! Ug_array = Uranium concentration for each grain (ppm) - assumed to be uniform
! Thg_array = Thorium concentration for each grain (ppm) - assumed to be uniform

! Open and read the input file Samples_settings.txt
if (is_unix) then
    open (52,file=run//'/data/Samples_settings.txt',status='old',action='read',err=991)
    goto 992
    991 write (*,*) 'Cant open '//run//'/data/Samples_settings.txt'
    stop
else 
    open (52,file=run//'\data\Samples_settings.txt',status='old',action='read',err=999)
    goto 992
    999 write (*,*) 'Cant open '//run//'\data\Samples_settings.txt'
    stop
endif
992 continue

read (52,*) dataFolder
read (52,*) nsamples
read (52,*) TotGrains
!---------------------------------------------
!----- Read observation from data file -----
!---------------------------------------------
allocate  (AHEOBS(TotGrains), THLOBS(TotGrains), OSLOBS(TotGrains), ESROBS(TotGrains),&
    AFTOBS(TotGrains),ZHEOBS(TotGrains),&
    ZFTOBS(TotGrains),KAROBS(TotGrains),BAROBS(TotGrains),&
    MAROBS(TotGrains),HAROBS(TotGrains),FT01OBS(TotGrains),&
    FT02OBS(TotGrains),FT03OBS(TotGrains),FT04OBS(TotGrains),&
    FT05OBS(TotGrains),FT06OBS(TotGrains),&
    FT07OBS(TotGrains),FT08OBS(TotGrains),FT09OBS(TotGrains),&
    FT10OBS(TotGrains),FT11OBS(TotGrains),&
    FT12OBS(TotGrains),FT13OBS(TotGrains),FT14OBS(TotGrains),&
    FT15OBS(TotGrains),FT16OBS(TotGrains),&
    FT17OBS(TotGrains),FT18OBS(TotGrains),FT19OBS(TotGrains),&
    FT20OBS(TotGrains),gsize_array(TotGrains),Ug_array(TotGrains),Thg_array(TotGrains),&
    gsizez_array(TotGrains),Ugz_array(TotGrains),Thgz_array(TotGrains))
allocate (ThL_Obs(TotGrains),OSL_Obs(TotGrains),ESR_Obs(TotGrains))
allocate (AFT_Obs(TotGrains),AHE_Obs(TotGrains),ZHE_Obs(TotGrains))
allocate (KAR_Obs(TotGrains),BAR_Obs(TotGrains),MAR_Obs(TotGrains),HAR_Obs(TotGrains))
allocate (KINAFT(TotGrains), KINAHE(TotGrains))

ThL_Obs = -9999
OSL_Obs = -9999
ESR_Obs = -9999
AFT_Obs = -9999
AHE_Obs = -9999
ZHE_Obs = -9999
KAR_Obs = -9999
BAR_Obs = -9999
MAR_Obs = -9999
HAR_Obs = -9999
KINAHE = -9999
KINAFT = -9999
call read_observations(run,dataFolder,TotGrains,nsamples,&
        gsize_array,Ug_array,Thg_array,KINAHE,KINAFT,gsizez_array,Ugz_array,Thgz_array,&
        AFT_Obs,AHE_Obs,ZHE_Obs,KAR_Obs,BAR_Obs,MAR_Obs,HAR_Obs,is_unix)
        
!-----------------------------------------------
read (52,*) AHe_flag,ZHe_flag,AFT_flag,KAr_flag,BAr_flag,MAr_flag,HAr_flag,THL_flag,&
            OSL_flag, ESR_flag
read (52,*) He43_flag
 ! Zonation not considered (yet)
read (52,*) dummy
! sample names
do i=1,nsamples
    read (52,*) sample_name(i),dummyReal,dummyReal,dummyReal
enddo

! Read the number of grain for each sample
allocate (ngrains_array(nsamples))
do i=1,nsamples
    read (52,*) ngrains_array(i)
enddo

! Read other flags
read (52,*) Alpha_flag,Alpha_Ejec_Flag,RDModel,D0,Ea ! Apatite He
read (52,*) DiffzModel,D0z,Eaz ! Zircon He
read (52,*) AnnModel, rhoST ! Apatite Fission track
read (52,*) D0k, Eak ! KAr
read (52,*) D0b, Eab ! BAr
read (52,*) D0m, Eam ! MAr
read (52,*) D0h, Eah ! HAr

if (He43_flag.eq.1) then
    !Read number of 4He/3He profiles
    call get_number_43He (dataFolder, nb_He43, is_unix, run, nsteps_He43_array)
    print *, 'Number of 4He/3He spectrum: ', nb_He43
    print *, 'Number of heating steps: ', nsteps_He43_array(1:nb_He43)

    ! allocate arrays
    allocate (sample_list_He43(nb_He43))
    allocate (size43(nb_He43),age43(nb_He43),dage43(nb_He43),UPPM43(nb_He43),age43_obs(nb_He43))
    allocate (THPPM43(nb_He43), RMR043(nb_He43))

    ! Read 4He/3He observations
    nb_step_max = maxval(nsteps_He43_array) ! take max number of heating step across spectra
    allocate (Step_duration(nb_step_max,nb_He43),Step_temperature(nb_step_max,nb_He43),RsRb(nb_step_max,nb_He43))
    allocate (dRsRb(nb_step_max,nb_He43),Sum3He(nb_step_max,nb_He43),dSum3He(nb_step_max,nb_He43))

    call read_data_files_for_43He_age_specific (dataFolder, nb_He43, is_unix,run, sample_list_He43,&
    size43,age43_obs,dage43,UPPM43, THPPM43, RMR043, nb_step_max, Step_temperature,Step_duration,&
    Sum3He,dSum3He,RsRb,dRsRb)

endif
close (52)

!Print within the console
write(*,*) 'Number of samples =' ,nsamples
write(*,*) 'Number of grains =', TotGrains
write(*,*) 'Thermochronometers to predict:'
if (AHe_flag.eq.1) then
    write(*,*) 'Apatite He: '
    if (Alpha_Ejec_Flag.eq.1) then
        write(*,*) '    Alpha stopping distances: ', 'Farley et al. (1996)'
    elseif (Alpha_Ejec_Flag.eq.2) then 
        write(*,*) '    Alpha stopping distances: ', 'Ketcham et al. (2011)'
    else
        write(*,*) '    Alpha stopping distances: ', 'no'
    endif

    if (RDModel.eq.0) then
        write(*,*) '    Diffusion model: ', 'Farley et al. (2000)'
    elseif (RDModel.eq.1) then
        write(*,*) '    Diffusion model: ', 'Shuster et al. (2006)'
    elseif (RDModel.eq.2) then
        write(*,*) '    Diffusion model: ', 'Gautheron et al. (2009)'
    elseif (RDModel.eq.3) then
        write(*,*) '    Diffusion model: ', 'Flowers et al. (2009) - RDAAM'
    elseif (RDModel.eq.4) then
        write(*,*) '    Diffusion model: ', 'Willett et al. (2017) - ADAM'
    endif

    write(*,*) '    D0 = ', D0, ' cm2/s'
    write(*,*) '    Ea = ', Ea, ' kJ/mol'
endif

if (He43_flag.eq.0) then
    write(*,*) 'Predict 4He/3He profiles: ', 'no'
else
    write(*,*) 'Predict 4He/3He profiles: ', 'yes'
endif
if (ZHe_flag.eq.1) then
  write(*,*) 'Zircon He: '
   write(*,*) '     D0 = ', D0z,  ' cm2/s'
   write(*,*) '     Ea = ', Eaz, ' kJ/mol'
   if (DiffzModel.eq.5) then
       write(*,*) '    Diffusion model: ', 'Reiners et al. (2004)'
   elseif (DiffzModel.eq.6) then
       write(*,*) '    Diffusion model: ',  'Guenthner et al. (2013)'
   endif
endif
if (KAr_flag.eq.1) then
  write(*,*) 'Feldspar Ar: '
   write(*,*) '     D0 = ', D0k,  ' cm2/s'
   write(*,*) '     Ea = ', Eak,  ' kJ/mol'
endif
if (BAr_flag.eq.1) then
  write(*,*) 'Biotite Ar: '
   write(*,*) '     D0 = ', D0b,  ' cm2/s'
   write(*,*) '     Ea = ', Eab,  ' kJ/mol'
endif
if (MAr_flag.eq.1) then
  write(*,*) 'Muscovite Ar: '
   write(*,*) '     D0 = ', D0m,  ' cm2/s'
   write(*,*) '     Ea = ', Eam,  ' kJ/mol'
endif
if (HAr_flag.eq.1) then
  write(*,*) 'Hornblende Ar: '
   write(*,*) '     D0 = ', D0h,  ' cm2/s'
   write(*,*) '     Ea = ', Eah,  ' kJ/mol'
endif
if (THL_flag.eq.1) then
  write(*,*) 'ThermoLuminescence: '
endif
if (OSL_flag.eq.1) then
  write(*,*) 'OSL: '
endif
if (ESR_flag.eq.1) then
  write(*,*) 'ESR: '
endif
write(*,*) '-------------------------------------------------'

! !!! ro read from file !!!
Init_FTL_Model = 0
kinFTLID = 4
Init_FTL_value = 16.3

! name of directory where data are stored 
dataFolder = trim(dataFolder)

!---------------------------------------------
!----- Read Tt paths of specific samples -----
!---------------------------------------------
! First read file to have the number of ntime
if (is_unix) then
    open (77,file=run//'/output/TimeTemperaturePaths.csv',status='old',action='read',IOSTAT=io)
else
    open (77,file=run//'\output\TimeTemperaturePaths.csv',status='old',action='read',IOSTAT=io)
endif
ntime = 0
do
    read (77,*,END=10)
    if (io.gt.0) then
        write (*,*) "Error happened when reading 'TimeTemperaturePaths.csv' - in CalcAge.f90"
    elseif (io.lt.0) then
        exit
    else
        ntime = ntime + 1
    endif
enddo
rewind (77)
10 close(77)

!----------------------------------------------------------------------------------
!Get number of Tt paths (to handle nsamples vs ngrains)
ntime = ntime -1 !one line counts for headers,another for blank line and the end of the file
allocate (temperature_samples_temp(ntime,nsamples))
allocate (time(ntime),temperature_samples(ntime,TotGrains))

sample_temp_flag = 1
if (is_unix) then
    open (77,file=run//'/output/TimeTemperaturePaths.csv',status='old',action='read',IOSTAT=io)
else
    open (77,file=run//'\output\TimeTemperaturePaths.csv',status='old',action='read',IOSTAT=io)
endif
rewind (77)
read (77,*) !skip first line
do i=1,ntime
    read (77,*, iostat=io) time(i), temperature_samples_temp(i,:)
    if (io.lt.0.or.io.gt.0) then
        sample_temp_flag = 0
        exit
    endif
enddo
rewind (77)
11 close(77)

if (sample_temp_flag.eq.1) then !The number of Tt paths is equal to the number of samples not grains
    deallocate (temperature_samples)
    allocate (temperature_samples(ntime,nsamples))
    temperature_samples = temperature_samples_temp
    deallocate (temperature_samples_temp)
else
    call read_Pecube_file (run,TotGrains,ntime,time,temperature_samples)
endif

!----------------------------------------------------------------------------------
ntime_temp = ntime
! ! Open file to write ages
! if (AHe_flag.eq.1) then
!   open (82,file=run//'/output/TimeAgeAHe.csv',status='unknown',action='write')
!   write (82, '(a)',advance='no') "Time"
!   do ks=1,nsamples !loop through samples
!       ngrains = ngrains_array(ks)
!       do i=1,ngrains
!           write (82,'(",",a)',advance='no') "Age"//trim(sample_name(ks))//"_"//c(i)
!       enddo
!   enddo
! write (82,*) 
! endif
allocate (ageAHe_array(ntime,TotGrains+1),ageAFT_array(ntime,TotGrains+1),&
            ageZHe_array(ntime,TotGrains+1),ageTime(ntime),rho_red(ntime),&
            ageKAr_array(ntime,TotGrains+1),ageBAr_array(ntime,TotGrains+1),&
            ageMAr_array(ntime,TotGrains+1),ageHAr_array(ntime,TotGrains+1),&
            MFTL_array(ntime,TotGrains+1),DMFTL_array(ntime,TotGrains+1))
ageAHe_array = -9999
ageAFT_array = -9999
ageZHe_array = -9999
ageKAr_array = -9999
ageBAr_array = -9999
ageMAr_array = -9999
ageHAr_array = -9999
MFTL_array = -9999
DMFTL_array = -9999
ageAHe_array(:,1) = time
ageAFT_array(:,1) = time
ageZHe_array(:,1) = time
ageKAr_array(:,1) = time
ageBAr_array(:,1) = time
ageMAr_array(:,1) = time
ageHAr_array(:,1) = time
MFTL_array(:,1) = time
DMFTL_array(:,1) = time
ageTime(:) = 0.0

! ! Open file for AFT
! if (AFT_flag.eq.1) then
! open (83,file=run//'/output/TimeAgeAFT.csv',status='unknown',action='write')
! write (83, '(a)',advance='no') "Time"
! do ks=1,nsamples !loop through samples
!       ngrains = ngrains_array(ks)
!       do i=1,ngrains
!           write (83,'(",",a)',advance='no') "Age"//trim(sample_name(ks))//"_"//c(i)
!       enddo
!   enddo
! write (83,*) 
! endif

! ! Open file for ZHe
! if (ZHe_flag.eq.1) then
! open (84,file=run//'/output/TimeAgeZHe.csv',status='unknown',action='write')
! write (84, '(a)',advance='no') "Time"
! do ks=1,nsamples !loop through samples
!       ngrains = ngrains_array(ks)
!       do i=1,ngrains
!           write (84,'(",",a)',advance='no') "Age"//trim(sample_name(ks))//"_"//c(i)
!       enddo
!   enddo
! write (84,*) 
! endif

! Do we want to compute 4He/3He predictions?
allocate (He_ratio(nb_step_max),TotHe3_released_array(nb_step_max,Totgrains),&
    He_ratioEjec(nb_step_max),TotHe3_released(nb_step_max),&
    He_ratio_array(nb_step_max,TotGrains),He_ratioEjec_array(nb_step_max,TotGrains))
TotHe3_released = 0.
TotHe3_released_array = 0.
He_ratio = 0.
He_ratio_array = 0.
He_ratioEjec = 0.
He_ratioEjec_array = 0.

! For Radiation damage model if you of ADAM model, read data once
if (RDModel.eq.4) then
    open (17, file='src/Willett_data.txt',status='old',action='read')
    !First row = EDD; Second row = Ea, third row = ln(D0/a²)
    read (17,*) data_ADAM_array
    close (17)
endif


!----------------------------------------------
!---- Loop through grains characteristics -----
!----------------------------------------------
! Here we compute for each grain the (U-Th)/He ages for apatite, and if
! specified we predict 4He/3He profile.
! If we want to predict 4He/3He profiles, the geological model (i.e., AHe age computation)
! is first run, and the 4He profile is stored in a vector at the end of the computation. Then,
! The diffusion model is run starting from the 4He concentration of the geological model, and uses
! the heating schedule provided by the user. 

cnt = 1
cnt_43 = 1
GrainCounter = 1
do ks=1,nsamples !loop through samples
    ngrains = ngrains_array(ks)
    
    !#################### Thermal Scenarios #############################
    allocate (ztime(ntime), ztemp(ntime),ztime_double(ntime),ztemp_double(ntime))
    allocate (ztime_real(ntime), ztemp_real(ntime))
    ztime = time
    ztime_double = time
 
    !Iterate through grains in samples
    do kg=1,ngrains
        
        allocate (ztime_heating(10),ztemp_heating(10))
        allocate (He43obs(10),SF3Heobs(10))

        print '(/,"Doing grain ",i2,1x,"of",1x,i2)',GrainCounter,TotGrains
        !Works iwe change the number of grain without running Pecube again and if all
        !samples have the same number of grains
        if (sample_temp_flag.eq.1) then
            ztemp = temperature_samples(:,ks)
        else
            ztemp = temperature_samples(:,GrainCounter)
        endif
        ztemp_double = ztemp
        gsize = gsize_array(GrainCounter)
        Ug = Ug_array(GrainCounter)
        Thg = Thg_array(GrainCounter)
         
         
        !------------ Compute AFT age? ----------------------------------------
        if (AFT_Obs(GrainCounter).eq.1.and.AnnModel.eq.2) then
            ! reverse time and temperature array (from 0 to x Ma)
            do irec=1,ntime
                ztime_real(irec)=ztime_double(ntime-irec+1) 
                if (ztemp_double(ntime-irec+1).gt.500) ztemp_double(ntime-irec+1) = 500
                ztemp_real(irec)=ztemp_double(ntime-irec+1) 
            enddo
          
            final_age = 0.0
            fmeanp = 0.0
            dfmeanp = 0.0
            oldest_age  = 0.0
            if (KINAFT(i).ne.-9999) then
                call ketch_main(ntime,ztime_real,ztemp_real,AnnModel,rhoST,Init_FTL_Model,KINAFT(GrainCounter),&
                kinFTLID,Init_FTL_value,final_age,oldest_age,fmeanp,dfmeanp,fdist)

                print *, 'Age AFT: ', final_age, fmeanp, dfmeanp
                if (final_age.le.1e-4) then
                    final_age = 1e-4
                endif
        
                ageAFT_array(:,cnt+1)=real(final_age,4)
                MFTL_array(:,cnt+1) = fmeanp
                DMFTL_array(:,cnt+1) = dfmeanp
            endif
            ! call FissionTrackAge(ztime_double,ztemp_double,ntime,rho_red,dummyReal,0,rhoST,FTLKIN(GrainCounter))
            ! do rr=2,ntime !See Ketcham (2005)
            !   ageAFT_array(rr,cnt+1) = rho_red(rr)/rhoST
            ! enddo
                
        endif
        !--------------------- End AFT ----------------------------------------

        age = 0.
        !-------------- Compute AHe age? --------------------------------------
        if (AHE_Obs(GrainCounter).eq.1) then

             ! Do we want 4He/3He profile ?
             if (He43_flag.eq.1) then 
                 ! find corresponding 4He/3He spectra by comparing sample names
                if (kg.le.9) then
                    write (grain_ID, '(I1)') kg
                elseif (kg.gt.9) then
                    write (grain_ID, '(I2)') kg
                else
                    write (grain_ID, '(I3)') kg
                endif
                 do khe = 1,nb_He43
                    print *, trim(sample_name(ks))//"_"//trim(grain_ID), trim(sample_list_He43(khe)), nb_He43
                    if (trim(sample_name(ks))//"_"//trim(grain_ID).eq.trim(sample_list_He43(khe))) then

                     print *, 'Predict 4He/3He spectrum for sample: ', trim(sample_name(ks))//"_"//trim(grain_ID)  

                     deallocate (ztime_heating,ztemp_heating)
                     deallocate (He43obs,SF3Heobs)

                     allocate (ztime_heating(nsteps_He43_array(khe)),ztemp_heating(nsteps_He43_array(khe)))
                     allocate (He43obs(nsteps_He43_array(khe)),SF3Heobs(nsteps_He43_array(khe)))

                     ! Get Heating schedule
                     ztime_heating = 0.d0
                     ztemp_heating = 0.d0
                     ztime_heating = Step_duration(1:nsteps_He43_array(ks),ks)
                     ztemp_heating = Step_temperature(1:nsteps_He43_array(ks),ks)

                     ! Compute age
                     call Hediff(ztime_double, ztemp_double, ntime, ztime_heating,ztemp_heating,&
                     nsteps_He43_array(ks), age43(khe), size43(khe),Alpha_Ejec_Flag, RDmodel,Uppm43(khe),THPPM43(khe),&
                     D0,Ea,KINAHE(GrainCounter),ageTime,&
                     He43_flag,He_ratio_array(1:nsteps_He43_array(ks),cnt_43),&
                     TotHe3_released_array(1:nsteps_He43_array(ks),cnt_43),He_ratioEjec_array(1:nsteps_He43_array(ks),cnt_43),&
                     RsRb,Sum3He,0,data_ADAM_array,1,0)

                     cnt_43 = cnt_43 +1 ! count 4He/3He profiles
                    !  go to 907

                    endif
                 enddo

             else
                ! Compute age
                print*, "No 4He/3He calculation"
                call Hediff(ztime_double, ztemp_double, ntime, ztime_heating,ztemp_heating,nstep_He43, age, gsize,&
                Alpha_Ejec_Flag, RDmodel,Ug,Thg,D0,Ea,KINAHE(GrainCounter),ageTime,0,He_ratio(1:nstep_He43),&
                TotHe3_released(1:nstep_He43),He_ratioEjec,He43obs,SF3Heobs,0,data_ADAM_array,0,0)
             endif
            
            ! 907 continue 
            
            ageAHe_array(:,cnt+1) = ageTime(:)
            
        endif
        !--------------------- End AHe ----------------------------------------
       
        ntime = ntime_temp

        !--------------------- Compute ZHe age? -------------------------------
        age = 0.
        if (ZHE_Obs(GrainCounter).eq.1) then
            call Hediff(ztime_double, ztemp_double, ntime,ztime_heating,ztemp_heating,nstep_He43, age, gsizez_array(GrainCounter),&
            Alpha_Ejec_Flag, DiffzModel, Ugz_array(GrainCounter),Thgz_array(GrainCounter),D0z,Eaz,KINAHE(GrainCounter),ageTime,&
            0,He_ratio(1:nstep_He43),TotHe3_released(1:nstep_He43),He_ratioEjec,He43obs,SF3Heobs,0,&
            data_ADAM_array,0,1)
            ageZHe_array(:,cnt+1) =  ageTime(:)
        endif
        !--------------------- End ZHe ----------------------------------------

        !--------------------- Compute KAr age? -------------------------------
        age = 0.
        if (KAR_Obs(GrainCounter).eq.1) then
            call Mad_He (ztime,ztemp,ntime,age,3,grainsize,D0k,Eak)
            ageKAr_array(ntime,cnt+1) = age
        endif
        !--------------------- End KAr ----------------------------------------

        !--------------------- Compute BAr age? -------------------------------
        age = 0.
        if (BAR_Obs(GrainCounter).eq.1) then
            call Mad_He (ztime,ztemp,ntime,age,4,grainsize,D0b,Eab)
            ageBAr_array(ntime,cnt+1) = age
        endif
        !--------------------- End BAr ----------------------------------------

        !--------------------- Compute MAr age? -------------------------------
        age = 0.
        if (MAR_Obs(GrainCounter).eq.1) then
            call Mad_He (ztime,ztemp,ntime,age,5,grainsize,D0m,Eam)
            ageMAr_array(ntime,cnt+1) = age
        endif
        !--------------------- End MAr ----------------------------------------

        !--------------------- Compute HAr age? -------------------------------
        age = 0.
        if (HAR_Obs(GrainCounter).eq.1) then
            call Mad_He (ztime,ztemp,ntime,age,6,grainsize,D0h,Eah)
            ageHAr_array(ntime,cnt+1) = age
        endif
        !--------------------- End HAr ----------------------------------------
        
        cnt = cnt +1
        GrainCounter = GrainCounter + 1
        
        
        deallocate (ztime_heating,ztemp_heating)
        deallocate (He43obs,SF3Heobs)
        
    enddo ! End of grains iteration
    deallocate (ztime, ztemp,ztime_double,ztemp_double,ztime_real,ztemp_real)
enddo ! End of samples iteration

! Write AHe ages in file
call write_ages(run,ntime,TotGrains,ageAHe_array(ntime,2:),ageAFT_array(ntime,2:),&
    ageZHe_array(ntime,2:),AHe_flag,AFT_flag,KAr_flag,BAr_flag,MAr_flag,HAr_flag,&
    AnnModel,ZHe_flag,ageKAr_array(ntime,2:),ageBAr_array(ntime,2:),&
    ageMAr_array(ntime,2:),ageHAr_array(ntime,2:),MFTL_array(ntime,2:),DMFTL_array(ntime,2:),is_unix)

! Write 4He/3He prediction in file 'Compare43HE.csv'
if (He43_flag.eq.1) then 
    call write_files_for_43He_age_specific (nb_He43, is_unix, run, sample_list_He43,&
    age43_obs, age43, nsteps_He43_array, nb_step_max,RsRb,He_ratio_array,Sum3He,&
    TotHe3_released_array,He_ratioEjec_array)
endif

 
write(*,*) '------------------------------------------------'
write(*,*) '------- End of Compute specific ages -----------'
write(*,*) '------------------------------------------------'
 
deallocate (time,temperature_samples,He_ratio,He_ratioEjec,TotHe3_released)
deallocate (ageAHe_array,ageZHe_array,ageAFT_array,ageKAr_array,ageBAr_array)
deallocate (ageMAr_array, MFTL_array,DMFTL_array)
deallocate (gsizez_array,Ugz_array,Thgz_array,ageTime,TotHe3_released_array)
deallocate (He_ratioEjec_array,KINAHE,KINAFT,rho_red)
deallocate (AFT_Obs,AHE_Obs,ZHE_Obs,KAR_Obs,BAR_Obs,MAR_Obs,HAR_obs)
if (He43_flag.eq.1) then
  deallocate (sample_list_He43)
  deallocate (size43,age43,dage43,UPPM43,THPPM43, RMR043,age43_obs)
  deallocate (Step_duration,Step_temperature,RsRb)
  deallocate (dRsRb,Sum3He,dSum3He)
endif
deallocate (AHEOBS,AFTOBS,ZHEOBS,ZFTOBS,KAROBS,BAROBS,&
    MAROBS,HAROBS,FT01OBS,FT02OBS,FT03OBS,FT04OBS,&
    FT05OBS,FT06OBS,FT07OBS,FT08OBS,FT09OBS,&
    FT10OBS,FT11OBS,FT12OBS,FT13OBS,FT14OBS,&
    FT15OBS,FT16OBS,FT17OBS,FT18OBS,FT19OBS,FT20OBS,&
    gsize_array,Ug_array,Thg_array)
    
end program calculate_ages_specific



!--------------------------------------------------------------------------
subroutine read_Pecube_file (run,TotGrains,ntime,time,temperature_samples)

! This subroutine read the TimeTemperature.csv file provided by Pecube,
! and stores the time and temperature vectors.
implicit none

character*(*) :: run
integer i,TotGrains,ntime,io
double precision time(ntime),temperature_samples(ntime,TotGrains)

open (77,file=run//'/output/TimeTemperaturePaths.csv',status='old',action='read')
rewind (77)
read (77,*) 
do i=1,ntime
    read (77,*,iostat=io) time(i), temperature_samples(i,:)
enddo
close (77)

end subroutine



!--------------------------------------------------------------------------

function read_string (unit, istring, jstring) result (out)

! Returns the string located at location istring,jstring in a csv file
! istring is column number and jstring is line or row number

implicit none

integer, intent(in)  :: unit, istring, jstring
character(:), allocatable  :: out

character*10000 line
character*1 del
integer :: i, eof, start, end, comma

out='NotFound'

rewind (unit)
  do  i = 1, jstring
  read (unit, '(a)', iostat = eof) line
  if (eof.ne.0) return
  enddo

del = ','
start = 1

  do i = 1, istring-1
  comma = index(line(start:), del)
  if (comma.eq.0) return
  start = start + comma
  enddo
comma = index(line(start:), del)
out = ''
if (comma.eq.1) return
end  = start + comma - 2
if (end.lt.start) end = len(trim(line))

out = line(start:end)

end function read_string



!--------------------------------------------------------------------------
function ucase(in) result (out)

! returns the uppercase version of the string in

implicit none

character (*), intent(in)  :: in
character(:), allocatable  :: out

integer                    :: i, j

out = trim(in)
do i = 1, len_trim(out)
j = iachar(out(i:i))
if ((j-97)*(j-122).le.0) out(i:i) = achar(j-32)
end do

end function ucase



!--------------------------------------------------------------------------

subroutine find_string (unit, word, iword, jword)

! Returns the position (istring, jstring) of a string from a csv file (unit=unit)
! jstring is line (row) number and istring is column number


implicit none

!----------------------------

interface ucase

function ucase(in) result (out)

character (*), intent(in)  :: in
character(:), allocatable  :: out

end function ucase

end interface
!--------------------------

character*(*) :: word
integer :: iword, jword, unit

character*1 del
character*10000 line
integer eof, pos, comma, start, end, ends

rewind (unit)

del = ','

eof = 0

jword=0

pos = 0

  do while (eof.ne.-1.and.pos.eq.0)
  read (unit, '(a)', iostat = eof) line
  jword = jword + 1
  pos = index (ucase(line), ucase(word))
    if (pos.ne.0) then
    comma = -1
    iword = 0
    end = pos - 1
    start = 1
      do while (comma.ne.0)
      comma = index(line(start:end), del)
      iword = iword + 1
      start = start + comma
      enddo
    ends = index(line(pos:), del)
      if (ends.eq.0) then
      ends = len_trim(line)
      else
      ends = pos + ends - 2
      endif
    if (ucase(line(start:ends)).eq.ucase(word)) return
    endif
  enddo

iword = 0
jword = 0

return

end subroutine find_string


!--------------------------------------------------------------------------
subroutine read_observations (run,dataFolder,TotGrains,nsamples,&
        gsize_array,Ug_array,Thg_array,KINAHE,KINAFT,gsizez_array,Ugz_array,Thgz_array,&
        AFT_Obs,AHE_Obs,ZHE_Obs,KAR_Obs,BAR_Obs,MAR_Obs,HAR_Obs,is_unix)

! This subroutine read the TimeTemperature.csv file provided by Pecube,
! and stores the time and temperature vectors.

implicit none
interface read_string

function read_string (unit, istring, jstring) result (out)

integer, intent(in)  :: unit, istring, jstring
character(:), allocatable  :: out

end function read_string

end interface

character*(*) :: run,dataFolder
integer :: nsamples, i, TotGrains,io,cnt, loc
integer :: isample,jsample,izuppm,jzuppm,izthppm,jzthppm,izsize,jzsize
integer :: iauppm,jauppm,iathppm,jathppm,iasize,jasize
integer :: irmr0,jrmr0,ikinAFT,jkinAFT,iahe,jahe,izhe,jzhe,iaft,jaft
integer :: ikar,ibar,imar,ihar,jkar,jbar,jmar,jhar
real gsize_array(TotGrains),Ug_array(TotGrains),Thg_array(TotGrains)
double precision KINAHE(TotGrains), KINAFT(TotGrains)
real gsizez_array(TotGrains),Ugz_array(TotGrains),Thgz_array(TotGrains)
integer AFT_Obs(TotGrains),AHE_Obs(TotGrains),ZHE_Obs(TotGrains)
integer KAR_Obs(TotGrains),BAR_Obs(TotGrains),MAR_Obs(TotGrains)
integer HAR_Obs(TotGrains)
! real,dimension(:),allocatable::AFTOBS,ZHEOBS,ZFTOBS
! real,dimension(:),allocatable::KAROBS,BAROBS,MAROBS
! real,dimension(:),allocatable::HAROBS,FT01OBS,FT02OBS
! real,dimension(:),allocatable::FT03OBS,FT04OBS,FT05OBS
! real,dimension(:),allocatable::FT06OBS,FT07OBS,FT08OBS
! real,dimension(:),allocatable::FT09OBS,FT10OBS,FT11OBS,FT12OBS
! real,dimension(:),allocatable::FT13OBS,FT14OBS,FT15OBS,FT16OBS
! real,dimension(:),allocatable::FT17OBS,FT18OBS,FT19OBS,FT20OBS
logical is_unix

character*100 :: sample, field


! Open input AGE file. it should have been written by Pecube
if (is_unix) then
    open (12,file=run//'/data/'//trim(dataFolder)//'/'//trim(dataFolder)//'.csv',status='old')
else
    open (12,file=run//'\data\'//trim(dataFolder)//'\'//trim(dataFolder)//'.csv',status='old')
endif

!###### Read into the file to store parameters ######
! First line are the headers
! find column of specific label
jsample = 0
call find_string (12, 'SAMPLE', isample, jsample)

loc = 0
if (isample.ne.0) then
    sample=''
    cnt = 0
    do while (trim(sample).ne.'NotFound')
        jsample = jsample + 1
        cnt = cnt + 1
        sample = read_string(12, isample, jsample)
        if (sample.ne.'') then
            if (sample.ne.'NotFound') then
                ! Find grain size
                call find_string (12, 'ASIZE', iasize, jasize)
                ! Find U concentration
                call find_string (12, 'AUPPM', iauppm, jauppm)
                ! Find Th concentration
                call find_string (12, 'ATHPPM', iathppm, jathppm)
                ! Find rmr0 AHE
                call find_string (12, 'KINAHE', irmr0, jrmr0)
                ! Find rmr0 AFT
                call find_string (12, 'KINAFT', ikinAFT, jkinAFT)
                ! Find grain size
                call find_string (12, 'ZSIZE', izsize, jzsize)
                ! Find U concentration
                call find_string (12, 'ZUPPM', izuppm, jzuppm)
                ! Find Th concentration
                call find_string (12, 'ZTHPPM', izthppm, jzthppm)
                ! Find AHe obs
                call find_string (12, 'AHE', iahe, jahe)
                ! Find ZHe obs
                call find_string (12, 'ZHE', izhe, jzhe)
                ! Find AFT obs
                call find_string (12, 'AFT', iaft, jaft)
                ! Find KAr obs
                call find_string (12, 'KAR', ikar, jkar)
                ! Find BAR obs
                call find_string (12, 'BAR', ibar, jbar)
                ! Find MAR obs
                call find_string (12, 'MAR', imar, jmar)
                ! Find HAR obs
                call find_string (12, 'HAR', ihar, jhar)
                
                !  Read zircon grains size
                if (iasize.ne.0) then
                    loc = 1
                    field = read_string(12, iasize, jsample)
                    if (field.eq.'NotFound') goto 110
                    if (field.ne.'') then
                        read(field,*,err=899) gsize_array(cnt)
                    else
                        gsize_array(cnt) = -9999
                    endif
                endif
                110 continue
                ! Read uppm
                if (iauppm.ne.0) then
                    loc = 2
                    field = read_string(12, iauppm, jsample)
                    if (field.eq.'NotFound') goto 111
                    if (field.ne.'') then
                        read(field,*,err=899) Ug_array(cnt)
                    else
                        Ug_array(cnt) = -9999
                    endif
                endif
                111 continue
                ! Read thppm 
                if (iathppm.ne.0) then
                    loc = 3
                    field = read_string(12, iathppm, jsample)
                    if (field.eq.'NotFound') goto 112
                    if (field.ne.'') then
                         read(field,*,err=899) Thg_array(cnt)
                    else
                         Thg_array(cnt) = -9999
                    endif
                endif
                112 continue
                ! Read rmr0 AHE
                if (irmr0.ne.0) then
                    loc = 4
                    field = read_string(12, irmr0, jsample)
                    if (field.eq.'NotFound') goto 113
                    if (field.ne.'') then
                        read(field,*,err=899) KINAHE(cnt)
                    else
                     KINAHE(cnt) = -9999
                    endif
                endif
                113 continue
                ! Read rmr0 AFT
                if (ikinAFT.ne.0) then
                    loc = 5
                    field = read_string(12, ikinAFT, jsample)
                    if (field.eq.'NotFound') goto 114
                    if (field.ne.'') then
                        read(field,*,err=899) KINAFT(cnt)
                    else
                     KINAFT(cnt) = -9999
                    endif
                endif
                114 continue
                ! Read zircon grains size
                if (izsize.ne.0) then
                    loc = 6
                    field = read_string(12, izsize, jsample)
                    if (field.eq.'NotFound') goto 115
                    if (field.ne.'') then
                        read(field,*,err=899) gsizez_array(cnt)
                    else
                        gsizez_array(cnt) = -9999
                    endif
                endif
                115 continue
                ! Read uppm
                if (izuppm.ne.0) then
                    loc = 7
                    field = read_string(12, izuppm, jsample)
                    if (field.eq.'NotFound') goto 116
                    if (field.ne.'') then
                        read(field,*,err=899) Ugz_array(cnt)
                    else
                        Ugz_array(cnt) = -9999
                    endif
                endif
                116 continue
                ! Read thppm 
                if (izthppm.ne.0) then
                    loc = 8
                    field = read_string(12, izthppm, jsample)
                    if (field.eq.'NotFound') goto 117
                     if (field.ne.'') then
                         read(field,*,err=899) Thgz_array(cnt)
                     else
                         Thgz_array(cnt) = -9999
                     endif
                endif
                117 continue
               ! is AHE obs ? 
               if (iahe.ne.0) then
                   loc = 9
                   field = read_string(12, iahe, jsample)
                   if (field.eq.'NotFound') goto 118
                    if (field.ne.'') then
                        AHE_Obs(cnt) = 1
                    else
                        AHE_Obs(cnt) = 0
                    endif
               endif
               118 continue
               ! is ZHE obs ? 
               if (izhe.ne.0) then
                   loc = 10
                   field = read_string(12, izhe, jsample)
                   if (field.eq.'NotFound') goto 119
                    if (field.ne.'') then
                        ZHE_Obs(cnt) = 1
                    else
                        ZHE_Obs(cnt) = 0
                    endif
                endif
                119 continue
               ! is AFT obs ? 
               if (iaft.ne.0) then
                   loc = 11
                   field = read_string(12, iaft, jsample)
                   if (field.eq.'NotFound') goto 120
                    if (field.ne.'') then
                        AFT_Obs(cnt) = 1
                    else
                        AFT_Obs(cnt) = 0
                    endif
               endif
               120 continue
               ! is KAr obs ? 
               if (ikar.ne.0) then
                   loc = 12
                   field = read_string(12, ikar, jsample)
                   if (field.eq.'NotFound') goto 121
                    if (field.ne.'') then
                        KAR_Obs(cnt) = 1
                    else
                        KAR_Obs(cnt) = 0
                    endif
                endif
                121 continue
               ! is BAr obs ? 
               if (ibar.ne.0) then
                   loc = 13
                   field = read_string(12, ibar, jsample)
                   if (field.eq.'NotFound') goto 122
                    if (field.ne.'') then
                        BAR_Obs(cnt) = 1
                    else
                        BAR_Obs(cnt) = 0
                    endif
               endif
               122 continue
               ! is MAr obs ? 
               if (imar.ne.0) then
                   loc = 14
                   field = read_string(12, imar, jsample)
                    if (field.ne.'') then
                        MAR_Obs(cnt) = 1
                    else
                        MAR_Obs(cnt) = 0
                    endif
               endif
               ! is HAr obs ? 
               if (ihar.ne.0) then
                   loc = 15
                   field = read_string(12, ihar, jsample)
                   if (field.eq.'NotFound') goto 123
                    if (field.ne.'') then
                        HAR_Obs(cnt) = 1
                    else
                        HAR_Obs(cnt) = 0
                    endif
               endif
               123 continue
            endif
        endif
    enddo
endif
close(12)

return
    
899 continue
print*,'Error in reading observations at line ',jsample, ' loc = ', loc, ' cnt = ', cnt
stop 'Wrong format'


    
    
end subroutine read_observations



!--------------------------------------------------------------------------
subroutine write_ages (run,ntime,nsamples,ageAHe_array,ageAFT_array,ageZHe_array,&
        AHe_flag,AFT_flag,KAr_flag,BAr_flag,MAr_flag,HAr_flag,AnnModel,ZHe_flag,&
        ageKAr_array,ageBAr_array,ageMAr_array,ageHAr_array,MFTL_array,DMFTL_array,is_unix)

! This subroutine read the TimeTemperature.csv file provided by Pecube,
! and stores the time and temperature vectors.
implicit none
character*(*) :: run
integer nsamples, i, ntime, TotGrains,io
integer AHe_flag,AFT_flag,ZHe_flag,AnnModel
integer KAr_flag,BAr_flag,MAr_flag,HAr_flag
real ageAHe_array(nsamples),ageAFT_array(nsamples),ageZHe_array(nsamples)
real ageKAr_array(nsamples),ageBAr_array(nsamples),ageMAr_array(nsamples)
real ageHAr_array(nsamples), MFTL_array(nsamples), DMFTL_array(nsamples)
character*50,dimension(:),allocatable::SAMPLE
double precision,dimension(:),allocatable::LON,LAT,HeightObs,HeightPred,AHEOBS,AHEPred
double precision,dimension(:),allocatable::AFTOBS,AFTPRED,ZHEOBS,ZHEPRED,ZFTOBS,ZFTPRED
double precision,dimension(:),allocatable::KAROBS,KARPRED,BAROBS,BARPRED,MAROBS,MARPRED
double precision,dimension(:),allocatable::HAROBS,HARPRED,MFTLOBS,MFTLPRED,DMFTLOBS,DMFTLPRED
double precision,dimension(:),allocatable::FT01OBS,FT01PRED,FT02OBS,FT02PRED
double precision,dimension(:),allocatable::FT03OBS,FT03PRED,FT04OBS,FT04PRED,FT05OBS,FT05PRED
double precision,dimension(:),allocatable::FT06OBS,FT06PRED,FT07OBS,FT07PRED,FT08OBS,FT08PRED
double precision,dimension(:),allocatable::FT09OBS,FT09PRED,FT10OBS,FT10PRED,FT11OBS,FT11PRED,FT12OBS,FT12PRED
double precision,dimension(:),allocatable::FT13OBS,FT13PRED,FT14OBS,FT14PRED,FT15OBS,FT15PRED,FT16OBS,FT16PRED
double precision,dimension(:),allocatable::FT17OBS,FT17PRED,FT18OBS,FT18PRED,FT19OBS,FT19PRED,FT20OBS,FT20PRED
logical is_unix

! Open CompareAGE file. it should have been written by Pecube
open (13,file=run//'/output/CompareAGE.csv',status='old')

! Allocate arrays
allocate  (SAMPLE(nsamples),LON(nsamples),LAT(nsamples),HeightObs(nsamples),HeightPred(nsamples),AHEOBS(nsamples),&
    AHEPred(nsamples),AFTOBS(nsamples),AFTPRED(nsamples),ZHEOBS(nsamples),ZHEPRED(nsamples),&
    ZFTOBS(nsamples),ZFTPRED(nsamples),KAROBS(nsamples),KARPRED(nsamples),BAROBS(nsamples),BARPRED(nsamples),&
    MAROBS(nsamples),MARPRED(nsamples),HAROBS(nsamples),HARPRED(nsamples),MFTLOBS(nsamples),MFTLPRED(nsamples),&
    DMFTLOBS(nsamples),DMFTLPRED(nsamples),&
    FT01OBS(nsamples),FT01PRED(nsamples),&
    FT02OBS(nsamples),FT02PRED(nsamples),FT03OBS(nsamples),FT03PRED(nsamples),FT04OBS(nsamples),&
    FT04PRED(nsamples),FT05OBS(nsamples),FT05PRED(nsamples),FT06OBS(nsamples),FT06PRED(nsamples),&
    FT07OBS(nsamples),FT07PRED(nsamples),FT08OBS(nsamples),FT08PRED(nsamples),FT09OBS(nsamples),&
    FT09PRED(nsamples),FT10OBS(nsamples),FT10PRED(nsamples),FT11OBS(nsamples),FT11PRED(nsamples),&
    FT12OBS(nsamples),FT12PRED(nsamples),FT13OBS(nsamples),FT13PRED(nsamples),FT14OBS(nsamples),&
    FT14PRED(nsamples),FT15OBS(nsamples),FT15PRED(nsamples),FT16OBS(nsamples),FT16PRED(nsamples),&
    FT17OBS(nsamples),FT17PRED(nsamples),FT18OBS(nsamples),FT18PRED(nsamples),FT19OBS(nsamples),&
    FT19PRED(nsamples),FT20OBS(nsamples),FT20PRED(nsamples))
    
AHEOBS = -9999
AHEPRED = -9999
AFTOBS = -9999
AFTPRED = -9999
ZHEOBS = -9999
ZHEPRED = -9999
ZFTOBS = -9999
ZFTPRED = -9999
KAROBS = -9999
KARPRED = -9999
BAROBS = -9999
BARPRED = -9999
MAROBS = -9999
MARPRED = -9999
HAROBS = -9999
HARPRED = -9999
MFTLOBS = -9999
MFTLPRED = -9999
DMFTLOBS = -9999
DMFTLPRED = -9999

!###### Read into the file to store parameters ######
! First line are the headers
! Sample, 'LON','LAT','HEIGHTOBS','HEIGHTPRED', &
!     'AHEOBS','AHEPRED','AFTOBS','AFTPRED','ZHEOBS','ZHEPRED','ZFTOBS','ZFTPRED', &
!     'KAROBS','KARPRED','BAROBS','BARPRED','MAROBS','MARPRED','HAROBS','HARPRED', &
!     'FT01OBS','FT01PRED','FT02OBS','FT02PRED','FT03OBS','FT03PRED','FT04OBS','FT04PRED', &
!     'FT05OBS','FT05PRED','FT06OBS','FT06PRED','FT07OBS','FT07PRED','FT08OBS','FT08PRED', &
!     'FT09OBS','FT09PRED','FT10OBS','FT10PRED','FT11OBS','FT11PRED','FT12OBS','FT12PRED', &
!     'FT13OBS','FT13PRED','FT14OBS','FT14PRED','FT15OBS','FT15PRED','FT16OBS','FT16PRED', &
!     'FT17OBS','FT17PRED','FT18OBS','FT18PRED','FT19OBS','FT19PRED','FT20OBS','FT20PRED'
read (13,*) 
do i=1,nsamples
    read (13,*,iostat=io) SAMPLE(i),LON(i),LAT(i),HeightObs(i),HeightPred(i),AHEOBS(i),&
        AHEPRED(i),AFTOBS(i),AFTPRED(i),ZHEOBS(i),ZHEPRED(i),&
        ZFTOBS(i),ZFTPRED(i),KAROBS(i),KARPRED(i),BAROBS(i),BARPRED(i),&
        MAROBS(i),MARPRED(i),HAROBS(i),HARPRED(i),MFTLOBS(i),MFTLPRED(i),&
        DMFTLOBS(i),DMFTLPRED(i),FT01OBS(i),FT01PRED(i),&
        FT02OBS(i),FT02PRED(i),FT03OBS(i),FT03PRED(i),FT04OBS(i),&
        FT04PRED(i),FT05OBS(i),FT05PRED(i),FT06OBS(i),FT06PRED(i),&
        FT07OBS(i),FT07PRED(i),FT08OBS(i),FT08PRED(i),FT09OBS(i),&
        FT09PRED(i),FT10OBS(i),FT10PRED(i),FT11OBS(i),FT11PRED(i),&
        FT12OBS(i),FT12PRED(i),FT13OBS(i),FT13PRED(i),FT14OBS(i),&
        FT14PRED(i),FT15OBS(i),FT15PRED(i),FT16OBS(i),FT16PRED(i),&
        FT17OBS(i),FT17PRED(i),FT18OBS(i),FT18PRED(i),FT19OBS(i),&
        FT19PRED(i),FT20OBS(i),FT20PRED(i)
enddo
close(13)

! Write updates in CompareAge.csv
if (is_unix) then
    open (13,file=run//'/output/CompareAGE.csv',status='unknown')
else 
    call CheckFileExistence(run//'\output\CompareAGE.csv')
    open (13,file=run//'\output\CompareAGE.csv',status='unknown')
endif
! first write Headers
write (13,'(a,64(",",a))') 'SAMPLE','LON','LAT','HEIGHTOBS','HEIGHTPRED', &
            'AHEOBS','AHEPRED','AFTOBS','AFTPRED','ZHEOBS','ZHEPRED','ZFTOBS','ZFTPRED', &
            'KAROBS','KARPRED','BAROBS','BARPRED','MAROBS','MARPRED','HAROBS','HARPRED', &
            'MFTLOBS','MFTLPRED','DMFTLOBS','DMFTLPRED',&
            'FT01OBS','FT01PRED','FT02OBS','FT02PRED','FT03OBS','FT03PRED','FT04OBS','FT04PRED', &
            'FT05OBS','FT05PRED','FT06OBS','FT06PRED','FT07OBS','FT07PRED','FT08OBS','FT08PRED', &
            'FT09OBS','FT09PRED','FT10OBS','FT10PRED','FT11OBS','FT11PRED','FT12OBS','FT12PRED', &
            'FT13OBS','FT13PRED','FT14OBS','FT14PRED','FT15OBS','FT15PRED','FT16OBS','FT16PRED', &
            'FT17OBS','FT17PRED','FT18OBS','FT18PRED','FT19OBS','FT19PRED','FT20OBS','FT20PRED'

! Then write new predictions
! first check thermochronometers flag
! if no prediction, get the default value read from CompareAge.csv file
if (AHe_flag.eq.0) ageAHe_array = AHEPRED
print *, 'AHe pred 2: ', AHEPRED
if (AFT_flag.eq.0) ageAFT_array = AFTPRED
if (ZHe_flag.eq.0) ageZHe_array = ZHEPRED
if (BAr_flag.eq.0) ageBAr_array = BARPRED
if (KAr_flag.eq.0) ageKAr_array = KARPRED
if (MAr_flag.eq.0) ageMAr_array = MARPRED
if (HAr_flag.eq.0) ageHAr_array = HARPRED
MFTL_array = MFTLOBS
DMFTL_array = DMFTLOBS

do i=1,nsamples 

    write (13,'(A,",",g12.6,63(",",g12.6))') trim(SAMPLE(i)),LON(i),LAT(i),HeightObs(i),HeightPred(i),AHEOBS(i),&
        ageAHe_array(i),AFTOBS(i),ageAFT_array(i),ZHEOBS(i),ageZHe_array(i),&
        ZFTOBS(i),ZFTPRED(i),KAROBS(i),ageKAr_array(i),BAROBS(i),ageBAr_array(i),&
        MAROBS(i),ageMAr_array(i),HAROBS(i),ageHAr_array(i),MFTL_array(i),MFTLPRED(i),&
        DMFTL_array(i),DMFTLPRED(i),FT01OBS(i),FT01PRED(i),&
        FT02OBS(i),FT02PRED(i),FT03OBS(i),FT03PRED(i),FT04OBS(i),&
        FT04PRED(i),FT05OBS(i),FT05PRED(i),FT06OBS(i),FT06PRED(i),&
        FT07OBS(i),FT07PRED(i),FT08OBS(i),FT08PRED(i),FT09OBS(i),&
        FT09PRED(i),FT10OBS(i),FT10PRED(i),FT11OBS(i),FT11PRED(i),&
        FT12OBS(i),FT12PRED(i),FT13OBS(i),FT13PRED(i),FT14OBS(i),&
        FT14PRED(i),FT15OBS(i),FT15PRED(i),FT16OBS(i),FT16PRED(i),&
        FT17OBS(i),FT17PRED(i),FT18OBS(i),FT18PRED(i),FT19OBS(i),&
        FT19PRED(i),FT20OBS(i),FT20PRED(i)
enddo
close(13)

end subroutine write_ages


!---------------------------------------------------------------------------------------------------
subroutine get_number_43He (dataFolder, nobs, is_unix, run, nhist)

use read_string_module

implicit none

character*(*) :: dataFolder,run
character(len=1024) :: fnme
integer :: isample, jsample, idurh, jdurh
integer :: nobs
integer, dimension(50) :: nhist
character*100 :: sample,sample_prev,field
character(len=1024) :: folderOS
logical is_unix

! first count the number of 43 samples in all files

nobs = 0
nhist = 0 ! Number os heating step for each spectrum

! Open input 4He/3He file. it should have been written by Pecube
if (is_unix) then
    open (9,file=run//'/data/'//trim(dataFolder)//'/43_data.csv',status='old')
else
    open (9,file=run//'\data\'//trim(dataFolder)//'\43_data.csv',status='old')
endif

! Find the location of SAMPLE field
call find_string (9, 'SAMPLE', isample, jsample)
call find_string (9, 'DUR43', idurh, jdurh)

if (isample.ne.0) then ! if found
    sample=''
    sample_prev = ''
    ! Read sample names
    do while (trim(sample).ne.'NotFound')
     jsample = jsample + 1
     ! Read sample ID
     sample = read_string(9, isample, jsample)
     if (trim(sample).eq.trim(sample_prev)) goto 888
     sample_prev = sample
     ! If reading is a success add sample to count
     if (sample.ne.''.and.sample.ne.'NotFound') then
        nobs = nobs + 1
     endif 
     888   continue
     ! Read duration
     field = read_string(9, idurh, jsample)
     if (field.ne.''.and.field.ne.'NotFound') then
        nhist(nobs) = nhist(nobs) + 1
     endif
    enddo ! end do while
endif

close (9) ! close 43He input file

if (nobs.eq.0) return ! No 4He/3He data found

end subroutine get_number_43He



!---------------------------------------------------------------------------------------------------
subroutine read_data_files_for_43He_age_specific (folder, nobs, is_unix, run, sample_list,&
    size43,age43,dage43,UPPM43, THPPM43, RMR043, nb_step_max, temph,durh,rel,drel,arel,darel)

! this subroutine finds and reads the thermal histories from all files contained in folder
! and writes it to unit 111 (the data input file for Pecube)
! temph: temperature at steps (°C)
! durh: duration of step (h)
! rel: cumsum3He ratio od steps
! drel: uncertainty on rel
! arel: Rs/Rb ratio od steps
! darel: uncertainty on arel

use read_string_module

implicit none

character*(*) :: folder,run
character(len=1024) :: fnme
integer :: idurh, jdurh, itemph, jtemph, irel, jrel, idrel, jdrel, iarel, jarel, idarel, jdarel
integer :: ilat, jlat, ilon, jlon, iheight, jheight, iage, jage, idage, jdage, isize, jsize, isample, jsample
integer :: nobs, nhist, iuppm, juppm, ithppm, jthppm, irmr0, jrmr0, nb_step_max
character*100 :: sample, field,sample_prev
double precision :: lat, lon
double precision, dimension(nb_step_max,nobs) :: temph,durh,rel,drel,arel,darel
real UPPM43(nobs), THPPM43(nobs), RMR043(nobs),age43(nobs),dage43(nobs),size43(nobs)
integer :: i,ksample, countSample
character(len=1024) :: folderOS
logical is_unix
character(len=20), dimension(nobs) :: sample_list

countSample = 0
fnme = '43_data.csv'

! Open input 4He/3He file. it should have been written by Pecube
if (is_unix) then
    open (9,file=run//'/data/'//trim(folder)//'/'//trim(fnme),status='old')
else
    open (9,file=run//'\data\'//trim(folder)//'\'//trim(fnme),status='old')
endif

! find sample line, heating schedule, and gas released data
call find_string (9, 'SAMPLE', isample, jsample)
call find_string (9, 'DUR43', idurh, jdurh)
call find_string (9, 'TEMP43', itemph, jtemph)
call find_string (9, 'rel', irel, jrel) ! Maxime
call find_string (9, 'drel', idrel, jdrel) !Maxime
call find_string (9, 'arel', iarel, jarel)
call find_string (9, 'darel', idarel, jdarel)

if (idurh.ne.0) then
 call find_string (9, 'SAMPLE', isample, jsample)

 if (itemph.eq.0) then
    print*,'Found DUR43 field but not TEMP43 field in file ',trim(folder)
    stop
 endif

 if (isample.ne.0) then
    sample=''
    sample_prev = ''
    nhist=0
    ! Loop until the end of file
    do while (trim(sample).ne.'NotFound')
     jsample = jsample + 1
     sample = read_string(9, isample, jsample)
     if (trim(sample).ne.trim(sample_prev).and.sample.ne.'NotFound') then
        countSample = countSample +1
        ! Write data in temporary file
        sample_list(countSample) = trim(sample)
        ! if (nhist.ne.0) then
        
        ! endif
        nhist = 0 
     endif
     sample_prev = sample
     if (trim(sample).ne.''.and.trim(sample).ne.'NotFound') then
      ! index Latitude
      call find_string (9, 'LAT', ilat, jlat)
      if (ilat.eq.0) then
        print*,'Need LAT to locate sample ',sample,' in file ',trim(folderOS)
        stop
      endif
      ! index Longitude
      call find_string (9, 'LON', ilon, jlon)
      if (ilat.eq.0) then
        print*,'Need LON to locate sample ',sample,' in file ',trim(folderOS)
        stop
      endif
      ! index grain size
      call find_string (9, 'SIZE', isize, jsize)
      if (isize.eq.0) then
        print*,'Need GRAIN SIZE for 43He data in sample ',sample,' in file ',trim(folderOS)
        stop
      endif
      ! index Age
      call find_string (9, 'AGE43', iage, jage)
      if (iage.eq.0) then
        print*,'Need AGE43 for 43He data in sample ',sample,' in file ',trim(folderOS)
        stop
      endif
      ! Maxime - Find U concentration
      call find_string (9, 'UPPM', iuppm, juppm)
      if (iuppm.eq.0) then
        print*,'Need UPPM for 43He data in sample ',sample,' in file ',trim(folderOS)
        stop
      endif
      ! Maxime - Find Th concentration
      call find_string (9, 'THPPM', ithppm, jthppm)
      if (ithppm.eq.0) then
        print*,'Need THPPM for 43He data in sample ',sample,' in file ',trim(folderOS)
        stop
      endif
      ! Maxime - Find RMR0 concentration
      call find_string (9, 'KINAHE', irmr0, jrmr0)
      if (irmr0.eq.0) then
        print*,'Need KINAHE for 43He data in sample ',sample,' in file ',trim(folderOS)
        stop
      endif
      ! index Error age
      call find_string (9, 'DAGE43', idage, jdage)
      ! Read latitude
      field = read_string(9, ilat, jsample)
      read(field,*) lat
      ! Read longitude
      field = read_string(9, ilon, jsample)
      read(field,*) lon
      ! Read grian size
      field = read_string(9, isize, jsample)
      read(field,*,err=899) size43(countSample)
      ! Read age
      field = read_string(9, iage, jsample)
      read(field,*,err=898) age43(countSample)
      ! Maxime - read UPPM
      field = read_string(9, iuppm, jsample)
      read(field,*,err=900) UPPM43(countSample)
      ! Maxime - read THPPM
      field = read_string(9, ithppm, jsample)
      read(field,*,err=901) THPPM43(countSample)
      ! Maxime - read RMR0
      field = read_string(9, irmr0, jsample)
      read(field,*,err=902) RMR043(countSample)
      ! Read error age
      dage43 = age43*0.1d0
      if (idage.ne.0) then
        field = read_string(9, idage, jsample)
        read(field,*,err=897) dage43(countSample)
      endif
      ! Read duration
      field = read_string(9, idurh, jsample)
      if (field.ne.''.and.field.ne.'NotFound') then
        nhist = nhist + 1
        field = read_string(9, idurh, jsample)
        read (field,*,err=896) durh(nhist,countSample)
        ! Read temperature
        field = read_string(9, itemph, jsample)
        read (field,*,err=895) temph(nhist,countSample)
        ! Read release 3He
        field = read_string(9, irel, jsample)
        read (field,*,err=894) rel(nhist,countSample)
        ! Read error release 3He
        if (idrel.ne.0) then
            field = read_string(9, idrel, jsample)
            read (field,*,err=893) drel(nhist,countSample)
        else
            drel(nhist,countSample) = rel(nhist,countSample)*0.1d0
        endif
        ! Read release 4He
        field = read_string(9, iarel, jsample)
        read (field,*,err=892) arel(nhist,countSample)
        ! Read error release 4He
        if (idarel.ne.0) then
            field = read_string(9, idarel, jsample)
            read (field,*,err=891) darel(nhist,countSample)
        else
            darel(nhist,countSample) = arel(nhist,countSample)*0.1d0
        endif
      endif
     endif

    enddo

 endif

endif

close (9)

994 continue

return

902 continue
print*,'Error in reading RMR0 concentration in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

901 continue
print*,'Error in reading Th concentration in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

900 continue
print*,'Error in reading U concentration in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

899 continue
print*,'Error in reading grain size in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

898 continue
print*,'Error in reading 43He age in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

897 continue
print*,'Error in reading uncertainty in 43He age in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

896 continue
print*,'Error in reading 43He heating step duration in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

895 continue
print*,'Error in reading 43He heating step temperature in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

894 continue
print*,'Error in reading percent gas released in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

893 continue
print*,'Error in reading uncertainty in percent gas released in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

892 continue
print*,'Error in reading age gas released in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

891 continue
print*,'Error in reading uncertainty in age gas released in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

end subroutine read_data_files_for_43He_age_specific


!---------------------------------------------------------------------------------------------------
subroutine write_files_for_43He_age_specific ( nobs, is_unix, run, sample_list,&
    age43_obs,age43_pred, nsteps_list, nstep_max,RsRb_obs,RsRb_pred,Sum3He_obs,&
    Sum3He_pred,RsRb_ejec)

! Write output file 'Compare43HE.csv' with predicted 4He/3He

implicit none

logical is_unix
character*(*) :: run
integer :: nsteps_list(50), nstep_max, i, j, nobs
real age43_obs(nobs), age43_pred(nobs)
character(len=20), dimension(nobs) :: sample_list
double precision :: RsRb_obs(nstep_max,nobs),RsRb_pred(nstep_max,nobs),RsRb_ejec(nstep_max,nobs)
double precision :: Sum3He_obs(nstep_max,nobs),Sum3He_pred(nstep_max,nobs)

! Write updates in CompareAge.csv
if (is_unix) then
    open (13,file=run//'/output/Compare43HE.csv',status='unknown')
else 
    call CheckFileExistence(run//'\output\Compare43HE.csv')
    open (13,file=run//'\output\Compare43HE.csv',status='unknown')
endif

! first write Headers
write (13,'(a,8(",",a))') 'SAMPLE','AGE43OBS','AGE43PRED','RELEASED',&
    'RELEASEDPRED','AGERELEASED','AGERELEASEDPRED','ARELEJEC'

! Write 4He/3He profiles
do i=1,nobs ! loop through sample
    do j=1,nsteps_list(i) ! loop through heating steps
    write (13,'(a,7(",",g12.6))')trim(sample_list(i)),age43_obs(i),age43_pred(i),&
        Sum3He_obs(j,i),Sum3He_pred(j,i),RsRb_obs(j,i),RsRb_pred(j,i),RsRb_ejec(j,i)
    enddo
enddo
close(13)


end subroutine write_files_for_43He_age_specific