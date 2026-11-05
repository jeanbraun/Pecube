!----------------------------------------------------------------------------------------------------
subroutine RD_model(diffusivities,rhov,RDmodel,ntime,D0,Ea,rho_r,temperature,time,He_flag,&
                    Uppm_array,Thppm_array,nrad,data_ADAM_array,grainsize,eU)

! This subroutine computes diffusivities from diffusion/radiation damage models.
! Radiation model <5 for apatite, >= 5 for zircon
! time in years

! Set variables
integer nrad, He_flag,RDmodel
real*4 lambda235, lambda238, lambda232, lambdaf, R
real,dimension(:),allocatable:: annealedDam
real*4 ETrap, psi, omega, etaq, Letch
real*4 t1, t2, trapDiff, temp, qfac, cfac, Eb, eU
real diffusivities(ntime), rhov(ntime)
real Uppm_array(nrad), Thppm_array(nrad)
double precision time(ntime),temperature(ntime)
real *8 D0, Ea
double precision rho_r(ntime,ntime)
double precision,dimension(3,1000000)::data_ADAM_array
real *4 grainsize
real U235atom(nrad),U238atom(nrad),TH232atom(nrad), alphae(ntime)
double precision fa, Ba, DI, lint, tort,SV,Dtort,dtf(ntime)
real*4 D01, E1,DN17,DN17a2
real lint_lattice,ts

! Decay constants
lambda232 = 4.948e-11
lambda238 = 1.511e-10
lambda235 = 9.849e-10
lambdaf = 8.46e-17
! Parameters
psi = 1.0e-13
omega = 1e-22
ETrap = 34 !Activation Energy for traps (kJ/mol)
etaq = 0.91 !For Durango
c = 3e-7 !(ppm/Ma)
Eb = 31.4 !kJ/mol
Letch = 8.15e-4 !Half of total etchable FT length
R = 0.008314472 ! kj/mol.k
! Zircon
Ba = 5.48e-19 ! g amorphized per alpha event
SV = 1.669 ! nm^-1 track surface to volume ratio
lint_lattice = 45920 ! nm, extrapolated to a zircon with 1e14 alphas/g
D01 = 193188 ! cm^2/s
E1 = 165 !KJ/mol
cfac = 3e-7 !1/(ppm.Ma)

! debug
!print *, 'Parameters RD_models: ', RDmodel,ntime,D0,Ea,He_flag,Uppm_array(nrad),Thppm_array(nrad),nrad,grainsize,eU

! Get diffusivity
! Compute effective track density
allocate(annealedDam(ntime))
annealedDam = 0.
alphae = 0.

if (RDmodel.eq.3) then !Flowers - Apatite --------------------------------------------------------
    do i=2,ntime
        annealedDam(i) = sum(rho_r(2:i,i)*(lambdaf/lambda238)*rhov(i)*etaq*Letch)
    enddo

    ! Get diffusivities
    do i=1,ntime-1
        temp = (temperature(i)+273.15 + temperature(i+1)+273.15) / 2 ! average temperature
        trapDiff = psi * annealedDam(i+1) + omega * annealedDam(i+1)**3
        diffusivities(i) = (D0*exp(-Ea/(R*temp)))/(1+ ( trapDiff*exp(ETrap/(R*temp)) ))*1e8 !µm²/s
    enddo

elseif (RDmodel.eq.2) then ! means Gautheron et al.(2009) ----------------------------------------
    
    ! Get diffusivities
    do i=2,ntime
        annealedDam(i) = sum(rho_r(2:i,i))*(time(i-1)-time(i))/(maxval(time)-time(i))
    enddo
    
    do i=1,ntime-1
        temp = (temperature(i)+273.15 + temperature(i+1)+273.15) / 2
        qfac = cfac * eU * annealedDam(i+1) * (maxval(time)-time(i))
        diffusivities(i) = (D0*exp(-Ea/(R*temp)))/(1+qfac*exp(Eb/(R*temp)))*1e8 !µm²/s
    enddo
    
elseif (RDmodel.eq.4) then !Willett et al. (2017) ------------------------------------------------
    call RD17(time,temperature,ntime,nrad,Uppm_array,Thppm_array,grainsize*1e-4,&
    D0/(grainsize*1e-4)**2,He_flag,diffusivities,rhov(ntime),data_ADAM_array)
    diffusivities(:) = (/diffusivities(2:size(diffusivities)),0./)

elseif (RDmodel.eq.6) then ! Guenthner et al. 2013 - Zircon --------------------------------------

    rhor_r = 0.d0
    dtf(2:ntime) = abs(time(2:ntime)-time(1:ntime-1))
    ts = 1e6*365.25*24*3600 !time in s

    call Guenthner13(ntime,time,temperature,rho_r,dtf*ts)
    ! Compute alpha dose - parent concentration in atoms/g
    U235atom = ((Uppm_array/1e6)*0.007204/238.02891)*6.02214179e23
    U238atom = ((Uppm_array/1e6)*0.992745/238.02891)*6.02214179e23
    TH232atom = ((Thppm_array/1e6)/232.03805)*6.02214179e23
    do i=1,ntime-1
      alphae(i) = 8*U238atom(1)*(exp(lambda238*time(i)*1e6)-exp(lambda238*time(i+1)*1e6)) &
        + 7*U235atom(1)*(exp(lambda235*time(i)*1e6)-exp(lambda235*time(i+1)*1e6)) &
        + 6*TH232atom(1)*(exp(lambda232*time(i)*1e6)-exp(lambda232*time(i+1)*1e6)) 
    enddo
    ! Equivalent alpha dose
    do i=1,ntime-1
        annealedDam(i) = annealedDam(i) + alphae(i)*sum(rho_r(2:i+1,i+1))
    enddo
    do i=1,ntime-1
        temp = (temperature(i)+273.15 + temperature(i+1)+273.15) / 2
        fa = 1 - exp(-Ba * annealedDam(i))
        DI = 1 - exp(-Ba * annealedDam(i)*3)
        lint = 4.2 / (fa*SV) - 2.5
        tort = (lint_lattice/lint)**2
        Dtort = (1/tort)*D01*exp(-E1/(R*temp)) ! cm^2/s
        Dtorta2 = Dtort / (grainsize*(1/1e4)*(1-DI))**2 ! 1/s
        DN17 = D0 * exp(-Ea/(R*temp)) !cm^2/s
        DN17a2 = DN17/(grainsize*(1/1e4)*DI)**2 ! 1/s
        diffusivities(i) = (DI/DN17a2+(1-DI)/Dtorta2)**(-1) ! cm^2/s
        diffusivities(i) = diffusivities(i)*(grainsize)**2! um^2/s
    enddo
  
else ! No radiation damage model ------------------------------------------------------------------
      do i=1,ntime-1
          temp = (temperature(i)+273.15 + temperature(i)+273.15) /2
          if (He_flag.eq.1) then
              temp = temperature(i)+273.15
          endif
          diffusivities(i) = D0*exp(-Ea/(R*temp))*1e8; !µm²/s
      enddo
endif
  

deallocate (annealedDam)

end subroutine RD_model


!-------------------------------------------------------------------------------------------
subroutine ftdW15(Uppm_bulk, Thppm_bulk, time_points, npts, erhos_W)
! Calculates an effective damage density at one time step for apatite He diffusion
! time must have the same endpoint as time_points
! Return the value for EDD at current time step

! Adapted from Greg Balco -- Berkeley Geochronology Center -- Feb 2009
! with previous modification by C. Willett and M. Fox -- 2015/2016
! Adapted by Maxime Bernard (Potsdam University) 

implicit none

!Variables
integer npts
double precision Uppm_bulk,Thppm_bulk,lambda238,lambda235,lambda232,lambdaf
double precision maxt,dty,LL,nq,Nu,N238,N235,N232,rho,rhov_B,erhos_W
double precision D238,D235,D232
double precision time_points(2),tago_W15(2),rho_v(2)

maxt = maxval(time_points)
dty = abs(time_points(npts) - time_points(npts-1))

!Set constants
lambda232 = 4.948e-11
lambda238 = 1.511e-10
lambda235 = 9.849e-10
lambdaf = 8.46e-17
 
!Track density eaquation constants 
LL = 8.1e-4
nq = 0.91
rho = 3.2 !Density for apatite

!----------- Compute atom number density concentrations of U, Th -------------
!Atoms/g concentrations of U, Th
NU = Uppm_bulk*1.e-6*6.02e23/238.02891 !g/g * (g/mol)^-1 = mol/g * ats/mol = ats/g
N238 = NU*(137.88/(137.88+1.))
N235 = NU*(1./(137.88+1))
N232 = Thppm_bulk*1.e-6*6.02e23/232.03806

!Convert to atoms/cm3
D238 = N238*rho !ats/g * g/cm3 = ats/cm3
D235 = N235*rho
D232 = N232*rho

!--------- Calculate track production in given time step ---------------------
!What time is it?
tago_W15 = maxt - time_points

!Do track generation calculation
!rho_v in tracks/cm3
rho_v = D238*(exp(lambda238*(tago_W15+dty))-exp(lambda238*tago_W15)) + &
        (7./8.)*D235*(exp(lambda235*(tago_W15+dty))-exp(lambda235*tago_W15)) + &
        (3./4.)*D232*(exp(lambda232*(tago_W15+dty))-exp(lambda232*tago_W15))
        
!Step 1 is a dummy step; set 0 to avoid confusion
rhov_B = rho_v(2)

!Effective damage density
!number of decays times lambda_fission/lambda_alpha times etching efficiency
!times etchable range of single fission fragment
erhos_W = rhov_B*(lambdaf/lambda238)*nq*LL

end subroutine ftdW15


!-------------------------------------------------------------------------------------------
subroutine RD17(time, temperature, ntime, nrad, Uppm_array, Thppm_array, rad, D0a2, He_flag, diffusivities, rhov_stored,&
            data_ADAM_array)

! This subroutine compute the diffusivity according to the amount of radiation damages (ADAM)
! This script has been provided by Chelsea Willett (Willett et al., 2017)
! Authors: Chelsea Willett (University of California, Berkeley) and Matthew Fox (University College London)
! Adapted by Maxime Bernard (Potsdam University)

implicit none
integer ntime,nrad,i,j,nstep,lower_t_fine,upper_t_fine,ind_EDD,aa,He_flag
double precision time(ntime), temperature(ntime)
real diffusivities(ntime)
real Uppm_array(nrad), Thppm_array(nrad)
real *4 rad
real *8 D0a2
double precision Uppm_W, Thppm_W
double precision Etrap, D0a2_W
real omega, psi, rhov_stored

!Decay constants
double precision lambda235, lambda238, lambda232

double precision R
real dtf(ntime)

!Variables
double precision EDD_prev,N238,N235,temp_step,u_fine,erhos_W,dt_small,EDD_now,EL_W
double precision lnD0a2_interp,Ea_interp,EDD_lower,EDD_higher,Ea_lower,Ea_higher
double precision lnD0a2_lower,lnD0a2_higher,weight,check_interp,c0_Ea,c1_Ea,c2_Ea,c3_Ea,c3_Ea_interp
double precision c0_d0,c1_d0,c2_d0,alpha0,alpha_Ea,delta_Ea,Ea_new,c3_d0,c3_d0_interp
double precision alpha_d0,delta_LnD0a2,lnD0a2_new,EDD_out_Ea_interp,Ea_new_lower
double precision Ea_new_higher,check_weight,EDD_out_lnD_interp,lnD_new_lower,lnD_new_higher
double precision EDD_new,diff_calc_EDD_Ea,D_temp,time_years(ntime),tempK(ntime)
double precision NTh_vector(ntime),N238_vector(ntime),N235_vector(ntime),NU_vector(ntime)
integer ind_Ea,ind_lnD,nb_data
double precision,dimension(:),allocatable::EDD_data,Ea_data,lnD0a2_data,Tago
double precision,dimension(:),allocatable::K_step
double precision,dimension(:),allocatable::time_points_fine,temp_points_fine,Tt_fine,temp_fine
double precision,dimension(:),allocatable::time_points_small
double precision,dimension(3,1000000)::data_ADAM_array

Etrap = 34*1.e3 !J/mol
omega = 1.e-22
psi= 10.e-13
nb_data = 1.e6
D0a2_W = D0a2

!No zonation
Uppm_W = Uppm_array(1)
Thppm_W = Thppm_array(1)

!Decay constants
lambda232 = 4.948e-11
lambda238 = 1.511e-10
lambda235 = 9.849e-10

R = 8.3144598 !J/(K*mol)
if (He_flag.eq.1) then
    dtf(2:ntime) = abs(time(2:ntime)-time(1:ntime-1)) !a minute time step 
    time_years = time !time is already in years
else
    dtf(:) = 1.e5 !smaller time step for diffusion
    time_years = (maxval(time)-time)*1e6
endif

tempK = temperature + 273.15 !temperature in Kelvin

EDD_prev = 0. !initialize to be cumulatively summed in solver loop

!----------- Initialisation ------------
!Open up the curves that relate Ea and lnD0a2 to damage
!Flowers et al (2009) supplementary figures
! allocate (data_array(3,nb_data))
allocate (EDD_data(nb_data),Ea_data(nb_data), lnD0a2_data(nb_data))

! open (17, file='src/Willett_data.txt',status='old',action='read')
! !First row = EDD; Second row = Ea, third row = ln(D0/a²)
! read (17,*) data_array
! close (17)

!Assign data to variables
EDD_data(:) = data_ADAM_array(1,:)
Ea_data(:) = data_ADAM_array(2,:)
lnD0a2_data(:) = data_ADAM_array(3,:)
EL_W = minval(Ea_data(:))*1000.

!Calculate how much parent you have at each time step going backwards
!from today's observed concentrations, no zonation
!Calculate [Th] going back
allocate (Tago(ntime))
Tago = maxval(time_years) - time_years
NTh_vector = Thppm_W*exp(lambda232*Tago)

!Calculate [U] going back, account for two isotopes
N238 = Uppm_W*(137.88/(137.88+1))
N238_vector = N238*exp(lambda238*Tago)

N235 = Uppm_W*(1/(137.88+1))
N235_vector = N235*exp(lambda235*Tago)
NU_vector = N238_vector + N235_vector

!Define parameters - from Willett et al. 2017 EPSL
c0_Ea = 58.6;
c1_Ea = 1.;
c2_Ea = -21820.;
c0_d0 = 58.4; ! fit done in Millions-of-years
c1_d0 = 1.; ! fit done in Millions-of-years
c2_d0 = -21700.; ! fit done in Millions-of-years

!Initialize vector to store diffusivities
allocate (K_step(ntime))
K_step = 0.

!------------- Compute diffusivities -------------------
allocate (time_points_fine(2), temp_points_fine(2))
allocate (time_points_small(2))
do i=2, ntime !step 1 is all zeros, start at step 2

    !Compute temperature
    temp_step = (tempK(i) + tempK(i-1))/2
    
    !Generate a miniature Tt path for subloop
    time_points_fine(:) = (/time_years(i-1), time_years(i)/)
    temp_points_fine(:) = (/tempK(i-1), tempK(i)/)
    nstep = (time_points_fine(2)-time_points_fine(1))/dtf(i)
    
    allocate(Tt_fine(nstep+1),temp_fine(nstep+1))
    Tt_fine = (/(time_points_fine(1)+(j-1)*dtf(i),j=1,nstep+1)/)

    !Redefine temperature on a regular vector
    temp_fine(:) = 0.
    temp_fine(1) = temp_points_fine(1)
    lower_t_fine = 1 !index of time points below current value of Tt
    upper_t_fine = 2 !index of time points above current value of Tt
    do j=2,nstep
        u_fine = (time_points_fine(lower_t_fine)-Tt_fine(j))&
            / (time_points_fine(lower_t_fine)-time_points_fine(upper_t_fine))
        temp_fine(j) = (1-u_fine)*temp_points_fine(lower_t_fine) + &
            u_fine*temp_points_fine(upper_t_fine)
        if (Tt_fine(j+1).gt.time_points_fine(upper_t_fine)) then
            lower_t_fine = lower_t_fine + 1
            upper_t_fine = upper_t_fine + 1
        endif
    enddo
    temp_fine(nstep+1) = temp_points_fine(2)

    !Now go within the detailed time-history
    do aa=2,size(temp_fine)
        dt_small = Tt_fine(aa) - Tt_fine(aa-1)
        
        !Compute EDD
        erhos_W = 0.
        time_points_small = (/Tt_fine(aa-1), Tt_fine(aa)/)
        
        if (He_Flag.eq.1.and.aa.eq.2.and.i.eq.2) then !start from erhos_W from the end of the geological model
            EDD_prev = rhov_stored
        endif
        
        call ftdW15(NU_vector(i),NTh_vector(i),time_points_small,2,erhos_W)
        EDD_now = EDD_prev + erhos_W

        !Read Ea off of Flowers 2009 Damage-Ea curves
        ind_EDD = 0.
        do j=1,nb_data
            if (EDD_now.le.EDD_data(j)) then
                ind_EDD = j
                exit
            endif
        enddo

        !interpolate between this selected neighbor and the other neighbor
        !that's just higher than the target value
        if (ind_EDD.eq.1) then
            Ea_interp = Ea_data(1)
            lnD0a2_interp = lnD0a2_data(1)
        else
            EDD_lower = EDD_data(ind_EDD-1)
            EDD_higher = EDD_data(ind_EDD)
            Ea_lower = Ea_data(ind_EDD-1)
            Ea_higher = Ea_data(ind_EDD)
            lnD0a2_lower = lnD0a2_data(ind_EDD-1)
            lnD0a2_higher = lnD0a2_data(ind_EDD)
            
            weight = (EDD_higher-EDD_now)/(EDD_higher-EDD_lower)
            check_interp = EDD_lower*weight+EDD_higher*(1-weight)
            
            Ea_interp = Ea_lower*weight+Ea_higher*(1-weight)
            lnD0a2_interp = lnD0a2_lower*weight+lnD0a2_higher*(1-weight)
        endif

        !Function to determine delta_Ea from ADAM paper
        !Load up the relationship between Ea and EDD from Flowers (2009)
        !Electronic Annex Figure 2
        if (ind_EDD.eq.1) then
            c3_Ea = Ea_data(ind_EDD) - EL_W/1000.
        else
            c3_Ea_interp = Ea_lower*weight+Ea_higher*(1-weight)
            c3_Ea = c3_Ea_interp - EL_W/1000.
        endif

        if (c3_Ea.lt.0) then
            print *, 'In RD17.f90: c3_Ea < 0 !'
            stop
            c3_Ea = 0
        endif
        
        !Equation (1) - Willett et al. 2017
        alpha0 = exp(c0_Ea)
        alpha_Ea = alpha0*exp(c2_Ea/temp_step)
        delta_Ea = c3_Ea*(exp(-alpha_Ea*(dt_small/1e6)**c1_Ea)-1) !for fit in Myrs
        !Add delta_Ea to previous Ea (in kJ)
        Ea_new = delta_Ea + Ea_interp

        !Function to determine delta_lnD0a2 from ADAM paper
        !Use relationship between lnD0a2 and EDD from Flowers et
        !al 2009 Electronic Annex Figure 2
        if (ind_EDD.eq.1) then
            c3_d0 = lnD0a2_data(ind_EDD) - log(D0a2_W)
        else
            c3_d0_interp = lnD0a2_lower*weight+lnD0a2_higher*(1-weight)
            c3_d0 = c3_d0_interp - log(D0a2_W)
        endif
        
        !Equation (2) - Willett et al. 2017
        alpha0 = exp(c0_d0)
        alpha_d0 = alpha0*exp(c2_d0/temp_step)
        delta_LnD0a2 = c3_d0*(exp(-alpha_d0*(dt_small/1e6)**c1_d0)-1) !for fit in Myrs
        !Add delta_Ea to previous Ea (in kJ)
        lnD0a2_new = lnD0a2_interp + delta_lnD0a2
        
        !Use new Ea and lnd0a2 to get new EDD
        ind_Ea = 0.
        do j=1,nb_data
            if (Ea_new.le.Ea_data(j)) then
                ind_Ea = j
                exit
            endif
        enddo
        ind_lnD = 0.
        do j=1,nb_data
            if (lnD0a2_new.le.lnD0a2_data(j)) then
                ind_lnD = j
                exit
            endif
        enddo
        
        !EDD from new Ea value
        if (ind_Ea.eq.1) then
            EDD_out_Ea_interp = EDD_data(1)
        else
            Ea_new_lower = Ea_data(ind_Ea-1)
            Ea_new_higher = Ea_data(ind_Ea)
            weight = (Ea_new_higher-Ea_new)/(Ea_new_higher-Ea_new_lower)
            check_weight = Ea_new_lower*weight+Ea_new_higher*(1-weight)
            if (abs(check_weight-Ea_new).gt.1.e-10) then
                print *,i, abs(check_weight-Ea_new)
                stop 'In RD17.f90: Rounding error for new EDD value'
            endif
            EDD_lower = EDD_data(ind_Ea-1)
            EDD_higher = EDD_data(ind_Ea)
            EDD_out_Ea_interp = EDD_lower*weight+EDD_higher*(1-weight)
        endif
        
        !EDD from new lnD0a2 value
        if (ind_lnD.eq.1) then
            EDD_out_lnD_interp = EDD_data(1)
        else
            lnD_new_lower = lnD0a2_data(ind_lnD-1)
            lnD_new_higher = lnD0a2_data(ind_lnD)
            weight = (lnD_new_higher-lnD0a2_new)/(lnD_new_higher-lnD_new_lower)
            check_weight = lnD_new_lower*weight+lnD_new_higher*(1-weight)
            EDD_lower = EDD_data(ind_lnD-1)
            EDD_higher = EDD_data(ind_lnD)
            EDD_out_lnD_interp = EDD_lower*weight+EDD_higher*(1-weight)
        endif
        EDD_new = EDD_out_Ea_interp
        EDD_prev = EDD_out_Ea_interp
        
        !Use damage and temperature to calculate the diffusivity 
        !according to RDAAM
        diff_calc_EDD_Ea = (D0a2_W*exp(-EL_W/(R*temp_step)))/ &
                        (1+(psi*EDD_new+omega*(EDD_new**3)) &
                        * exp(Etrap/(R*temp_step)))
                        
        !Convert to K for use in age equation and He diffusion
        D_temp = diff_calc_EDD_Ea*(rad**2) !Now in cm2/s
        
    enddo
    K_step(i) = D_temp
    
    deallocate (Tt_fine,temp_fine)
    enddo
    !store EDD of the geological model to be used later
    !for 4He/3He predicitons
    if (He_Flag.eq.0) then
        rhov_stored = Edd_prev
    endif
    diffusivities(:) = K_step(:)*1e8 !µm²/s

deallocate (EDD_data, Ea_data, lnD0a2_data)
deallocate (Tago)
deallocate (K_step,time_points_fine,temp_points_fine,time_points_small)

end subroutine RD17



!-------------------------------------------------------------------------------------------
subroutine Guenthner13(ntime,time,temperature,damage,dtf)

  ! Damage annealing model for Zircon - Guenthner et al. 2013
  ! provided by William Guenthner - University of Illinois
  
  implicit none
  
  integer ntime, i, j
  double precision time(ntime),temperature(ntime)
  double precision damage(ntime,ntime),temp,teq,rcb2(ntime,ntime),t
  double precision c0,c1,c2,c3,alpha,eqTotAnnLen,last_rs,dtf(ntime)
  
  ! Fitting coefficients of Yamada data
  c0 = 6.24534
  c1 = -0.11977
  c2 = -314.937
  c3 = -14.2868
  alpha = - 0.05721
  
  ! Create damage array to hold tracks
  rcb2 = 0.
  do i=2,ntime
     temp = (temperature(i-1)+temperature(i))/2 +273.15
     rcb2(i,i) = 1./( (c0+c1*(log(dtf(i))-c2) / (log(1/temp)-c3))**(1/alpha)+1)
  enddo

  do i=3,ntime
      temp = (temperature(i-1)+temperature(i))/2+273.15
      do j=2,i-1
          last_rs = rcb2(j,i-1)
          teq = exp(c2 + ((log(1/temp)-c3)/c1)*((((1/last_rs)-1)**alpha)-c0))
          rcb2(j,i) = 1/( (c0 + c1* (log(dtf(j)+teq)-c2) / ((log(1./temp))-c3) )**(1./alpha) + 1 )
      enddo
   enddo

  ! Volume conversion
  damage(:,:) = 1.25*(rcb2(:,:) - 0.2)
  
  ! Suppress situations where damage < equivalent total annealing length
  eqTotAnnLen = 0.36/1.25 + 0.2
  do i=1,ntime
     do j=1,ntime
         if (damage(i,j).lt.eqTotAnnLen) then
             damage(i,j) = 0.
         endif
     enddo
   enddo
   !print *, 'damages = ', damage

    
end subroutine Guenthner13
