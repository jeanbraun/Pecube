!----------------------------------------------------------------------------------------------------
subroutine Hediff(time, temperature, ntime, duration, temp_heat, ndur, Apatite_age, grainsize,&
 Alpha_Ejec_flag,RDmodel,Uppm,Thppm,D0,Ea,rmr0,ageTime,He_flag,He_ratio,&
 TotHe3_released,He_ratioEjec,He43obs,SF3Heobs,He43_inv_flag,data_ADAM_array,Pecube,mineral)

!-------------------------------------------------------------------------------
! % This script computes (U-Th)/He ages on apatite and zircon, and simulates degassing
! % experiments to predict 4He/3He profiles.
! %
! % It includes alpha ejection models and alpha stopping distances from:
! %  Farley et al. (1996)
! %  Ketcham et al. (2011)
! % 
! % It also includes radiation damage models from
! %  Flowers et al. (2009)
! %  Gautheron et al. (2009)
! %  Willett et al. (2017)
! %  (see file 'RD_Model.f90')
! %
! %
! % This script is based on the Matlab script provided by Guenther et al.
! % (2020)

!

! % Last updates: 16/11/2022
! % Author: Maxime Bernard - Potsdam University (maxime.bernard@uni-potsdam.de)

 implicit none
 ! Counters
 integer i, tstep,counter,Pecube,stepCounter,nstep,mineral
 
 ! Decay constants
 real*4 lambda235, lambda238, lambda232
 real*4 lambdaf, lambdaD 
 integer Alpha_Ejec_flag, RDmodel, kk
 
 real*4 dr
 
 ! Input parameters
 double precision time(ntime),temperature(ntime),time_temp(ntime)
 double precision time_reversed(ntime),temperature_reversed(ntime)
 double precision,dimension(:),allocatable::time_out,temp_out,time_highres,temp_highres
 real *4 grainsize, Apatite_age
 integer ntime,He_flag, ntime_highres, status
 integer MAX_NUM_TIME_STEPS
 
 ! Parent
 real *4,dimension(:),allocatable::U238_array, U235_array, Th232_array, Uppm_array, Thppm_array
 real *4 Uppm, Thppm, eU, Uc, dti
 
 ! Grain radius
 integer nrad !number of point in radius
 real,dimension(:),allocatable::radius
 
 ! Alpha stopping distances
 real*4 s238, s235, s232, X0, Xs,D,Fr
 real*4 nd232, nd235, nd238, rad_pos, total232, total235, total238
 real*4 d232, d235, d238, ft232, ft235, ft238, innerVol, outerVol, vol
 real,dimension(:),allocatable::aE238_rad,aE235_rad,aE232_rad
 real,dimension(:),allocatable::Eprod_U238,Eprod_U235,Eprod_Th232
 
 ! Radiation damage model
 real*8 D0, Ea
 double precision rmr0
 real*4 t1, t2, trapDiff, temp
 
 ! Age computation
 real,dimension(:),allocatable::diffusivities, rhov, rhov_stored
 double precision,dimension(:,:),allocatable::rho_r

 ! Diffusion model
 real, dimension(:), allocatable::diag, upper, lower, HeSphere, f
 real, dimension(:), allocatable::heProfile, integrand
 real*4 dt0, DT, beta, He_prod, minHe, res, totalHe, leftSUm, loAge, hiAge
 real*4 midAge, midVal, heModelAge,Temp_crit
 integer ageConv
 real, parameter::pi = 3.1415926
 double precision HeSphere_stored(513),HeSphere_prev(513)
 double precision,dimension(3,1000000)::data_ADAM_array
 double precision, dimension(:), allocatable:: dtf
 real ageTime(ntime), ts
 double precision ageTime_dp(ntime)
 double precision, dimension(:), allocatable::ageTime_highres
 
 ! 4He/3He
 integer counterStep, ndur,istep,nstep_He43,He43_inv_flag
 double precision duration(ndur), temp_heat(ndur), time_heat(ndur)
 double precision duration2(ndur+1), temp_heat2(ndur+1)
 double precision, dimension(:), allocatable::time_step_array,temp_step_array
 real He4_frac,He3_frac,res3,totalHe3,He4_frac_next,He3_frac_next,He3_frac_step,He4_fracEjec_next
 real He3_frac_cum,He4_frac_cum,He4Ejec_frac_cum,He4_frac_Ejec,fstep
 real f3(513),He3Sphere(513),integrand3(513),he3Profile(513),He3Sphere_init,He3Sphere_prev(513)
 real HeSphereEjec_prev(513),HeSphereEjec(513),fejec(513),integrandEjec(513),heEjecProfile(513)
 double precision TotHe3_released_temp(5001),bulk, He_ratio_temp(5001),He3_released(5001)
 double precision He_ratioEjec(ndur),He_ratioEjec_temp(5001),HeSphereEjec_stored(513)
 double precision He_ratio(ndur), TotHe3_released(ndur)
 double precision He43obs(ndur),SF3Heobs(ndur)
 
 !##########################################################################
 !################# set parameters of one grain ############################
 !##########################################################################
 nrad = 513 !number of points through the grain radius
 dr = grainsize/nrad !distance between each point
 
 ! Decay constants
 lambda232 = 4.948e-11
 lambda238 = 1.511e-10
 lambda235 = 9.849e-10
 lambdaf = 8.46e-17
 lambdaD = lambda238
 He_ratio_temp = 0.0
 He3_released = 0.0
 MAX_NUM_TIME_STEPS = 5000   ! <-- adjust as needed
 
 ! Handle zeros value U and Th
 Uppm = Uppm + 1e-10
 Thppm = Thppm + 1e-10
 
 ! Parent concentration of a grain 
 eU = Uppm + 0.235*Thppm
 Uc = Uppm * 1e-6 * 6.02e23 / 238.02891 !Total uranium atomes/g
 if (He_flag.eq.0.and.Pecube.eq.0) then
    print '(/,"     eU (ppm) = ",f7.2)', eU
 elseif (He_flag.gt.1) then
    print *, 'He43_flag = ', He_flag
    stop 'option not supported in Hediff.f90 - 4He/3He flag'
 endif
 allocate (radius(nrad), U238_array(nrad), U235_array(nrad), Th232_array(nrad), Uppm_array(nrad), Thppm_array(nrad))

 ! Fill radius array
 do i=1,nrad
    radius(i) = (i-1)*dr
 enddo
 ! Parent concentration of one grain
 ! Fill U and Th arrays along radius - assumed to be uniform
 Uppm_array = 0.0
 Thppm_array= 0.0
 do i=1,nrad
    U238_array(i) = Uc * (137.88/(137.88+1))
    U235_array(i) = Uc * (1/(137.88+1))
    Th232_array(i) = Thppm * 1e-6 * 6.02e23 / 232.0306
    Uppm_array(i) = Uppm
    Thppm_array(i) = Thppm
 enddo

 ! Debug
!   if (mineral.eq.1) then ! 1: zircon, 0: apatite
!     print *, 'Zircon'
!     print *, 'Size = ', grainsize
!     print *, 'Uppm = ', Uppm
!     print *, 'Thppm = ', Thppm
!     print *, 'rmr0 = ', rmr0
!     ! print *, 'Temperature = ', temperature
!     ! print *, 'time = ', time
!     ! print *, 'Temp = ', temp_heat
!     ! print *, 'Duration = ', duration
!     ! print *, 'Nheating = ', ndur
!     print *, 'RDmodel = ', RDmodel
!     print *, 'D0 = ', D0
!     print *, 'Ea = ', Ea
!   else
!     print *, 'Apatite'
!     print *, 'Size = ', grainsize
!     print *, 'Uppm = ', Uppm
!     print *, 'Thppm = ', Thppm
!     print *, 'rmr0 = ', rmr0
!     ! print *, 'Temperature = ', temperature
!     ! print *, 'time = ', time
!     ! print *, 'Temp = ', temp_heat
!     ! print *, 'Duration = ', duration
!     ! print *, 'Nheating = ', ndur
!     print *, 'RDmodel = ', RDmodel
!     print *, 'D0 = ', D0
!     print *, 'Ea = ', Ea
!   endif
 
 !#########################################################################
 !######################### Compute alpha ejection ########################
 !#########################################################################
 ! Alpha stopping distances (apatite, µm)
 s238 = 0
 s235 = 0
 s232 = 0
 if (mineral.eq.0) then 
     Temp_crit = 200 ! critical temperature for He apatite

     if (Alpha_Ejec_flag.eq.0) then
     ! No computation
          s238 = 0
          s235 = 0
          s232 = 0
      elseif (Alpha_Ejec_flag.eq.1) then
        ! Farley 1996 (µm)
          s238 = 19.68
          s235 = 22.83
          s232 = 22.46
      elseif (Alpha_Ejec_flag.eq.2) then
        ! Ketcham (2011)
          s238 = 18.81
          s235 = 21.80
          s232 = 22.25
      else
          stop 'option not supported in Hediff.f90 - Alpha stopping distances'
      endif
  elseif (mineral.eq.1) then ! Zircon
      Temp_crit = 250 ! critical temperature for He zircon
      if (Alpha_Ejec_flag.eq.0) then
      ! No computation
           s238 = 0
           s235 = 0
           s232 = 0
       elseif (Alpha_Ejec_flag.eq.1) then
         ! Farley 1996 (µm)
           s238 = 16.65
           s235 = 19.64
           s232 = 19.32
       elseif (Alpha_Ejec_flag.eq.2) then
         ! Ketcham (2011)
           s238 = 15.55
           s235 = 18.05
           s232 = 18.43
       else
           stop 'option not supported in Hediff.f90 - Alpha stopping distances'
       endif
  endif


  !#### Correction for alpha ejection #####
  allocate (aE238_rad(nrad),aE232_rad(nrad),aE235_rad(nrad))
  allocate (Eprod_U238(nrad),Eprod_U235(nrad),Eprod_Th232(nrad))
  
  Eprod_U238(:) = U238_array
  Eprod_U235(:) = U235_array
  Eprod_Th232(:) = Th232_array
  
  ! Alpha ejection for each isotope along radius (no zonation)
  do i=1,nrad
    X0 = (i-0.5)*dr;
    if (X0.ge.grainsize-s238) then
        Xs = (X0**2 + grainsize**2 - s238**2) / (2*X0)
        Fr = 0.5 + (Xs-X0) / (2*s238) !Fraction retained
        D = Eprod_U238(i) / U238_array(i)
        aE238_rad(i) = Fr * D * U238_array(i)
    else
        aE238_rad(i) = U238_array(i)
    endif
  
    if (X0.ge.grainsize-s235) then
        Xs = (X0**2 + grainsize**2 - s235**2) / (2*X0)
        Fr = 0.5 + (Xs-X0) / (2*s235) !Fraction retained
        D = Eprod_U235(i) / U235_array(i)
        aE235_rad(i) = Fr * D * U235_array(i)
    else
        aE235_rad(i) = U235_array(i)
    endif
    
    if (X0.ge.grainsize-s232) then
        Xs = (X0**2 + grainsize**2 - s232**2) / (2*X0)
        Fr = 0.5 + (Xs-X0) / (2*s232) !Fraction retained
        D = Eprod_Th232(i) / Th232_array(i)
        aE232_rad(i) = Fr * D * Th232_array(i)
    else
        aE232_rad(i) = Th232_array(i)
    endif
  enddo 
  
  ! Calculate the total amount of isotope in each X (no zonation)
  innerVol = 0.
  rad_pos = 0.
  total238 = 0.
  total235 = 0.
  total232 = 0.
  
  do i=1,nrad
    rad_pos = rad_pos + dr
    outerVol = rad_pos**3
    total238 = total238 + U238_array(i) * (outerVol-innerVol)
    total235 = total235 + U235_array(i) * (outerVol-innerVol)
    total232 = total232 + Th232_array(i) * (outerVol-innerVol)
    innerVol = outerVol
  enddo
    
  ! Scale by He production and volume with a base of /(1.333*pi)
  total238 = 8*total238
  total235 = 7*total235
  total232 = 6*total232
 
  ! Alpha ejection correction factors - Farley et al. 1996
  innerVol = 0.
  rad_pos = 0.
  d232 = 0.
  d235 = 0.
  d238 = 0.
  nd232 = 0.
  nd235 = 0.
  nd238 = 0.
  do i=1,nrad
    rad_pos = rad_pos + dr
    outerVol = rad_pos**3
    vol = outerVol - innerVol
    d232 = d232 + vol * aE232_rad(i)
    d235 = d235 + vol * aE235_rad(i)
    d238 = d238 + vol * aE238_rad(i)
    nd232 = nd232 + vol * Th232_array(i)
    nd235 = nd235 + vol * U235_array(i)
    nd238 = nd238 + vol * U238_array(i)
    innerVol = outerVol
  enddo
  

    if (Th232_array(1).eq.0) then
      ft232 = 0
    else 
      ft232 = d232/nd232
    endif
    ft235 = d235/nd235
    ft238 = d238/nd238

  
  !#########################################################################
  !##########  4He/3He and Compute  4He profile - ejection only  ###########
  !#########################################################################
  ! Make duration - heat schedule a minute time step
  if (He_flag.eq.1) then
      allocate (time_step_array(5000), temp_step_array(5000))
      ! add one heating step with high temperature to ensure all degasing
      duration2(1:ndur) = duration
      temp_heat2(1:ndur) = temp_heat
      duration2(ndur+1) = 2 ! the step is one hour
      temp_heat2(ndur+1) = 1000 ! at 1000°C 
      dt0 = 60.d0 ! a minute step heating (s)
      stepCounter = 0
      time_step_array = 0.d0
      temp_step_array = 0.d0
      time_heat = 0.
      ! Get time and temperature arrays for a minute step heating
      do i=1,ndur+1
          nstep=max(1,int((duration2(i)*3600.d0+tiny(dt0))/dt0)) !duration is in hours
          dti=duration2(i)/nstep*3600.d0 !get in second with a minute step
          do istep=1,nstep
              stepCounter = stepCounter+1
              fstep=float(istep)/(nstep)
              time_step_array(stepCounter) = dti
              temp_step_array(stepCounter) = temp_heat2(i)
          enddo
      enddo

      ! Make cumulated time
      do i=2,stepCounter
          time_step_array(i) = time_step_array(i-1)+time_step_array(i)
      enddo
      ! reverse time
      time_step_array = maxval(time_step_array(1:stepCounter)) - time_step_array(1:stepCounter)
      temp_step_array = temp_step_array(1:stepCounter)
      nstep_He43= size(time_step_array)
  endif

  ! Compute for ejection only
  ! only when computing 4He/3He
  if (He_flag.eq.1) then
  
      allocate (heProfile(nrad))  
      heProfile(:) = 0.
      
      do tstep=1,ntime-1
            t1 = time(tstep)*1e6   !older time (yr)
            t2 = time(tstep+1)*1e6 !younger time  
            dt0 = abs(t1-t2)*365.25*24*60*60 !seconds
    
            ! Production for the first node (atoms/g)
            He_prod = 8*aE238_rad(1)*(exp(lambda238*t1)-exp(lambda238*t2)) &
                    + 7*aE235_rad(1)*(exp(lambda235*t1)-exp(lambda235*t2)) &
                    + 6*aE232_rad(1)*(exp(lambda232*t1)-exp(lambda232*t2)) 
            heProfile(1) = heProfile(1) + He_prod
    
            ! Dirichlet boundary condition
            He_prod = 8*aE238_rad(nrad)*(exp(lambda238*t1)-exp(lambda238*t2)) &
                    + 7*aE235_rad(nrad)*(exp(lambda235*t1)-exp(lambda235*t2)) &
                    + 6*aE232_rad(nrad)*(exp(lambda232*t1)-exp(lambda232*t2)) 
            heProfile(nrad) = heProfile(nrad) + He_prod
    
            ! Compute for intermediate nodes
            do i=2,nrad-1
                He_prod = 8*aE238_rad(i)*(exp(lambda238*t1)-exp(lambda238*t2)) &
                        + 7*aE235_rad(i)*(exp(lambda235*t1)-exp(lambda235*t2)) &
                        + 6*aE232_rad(i)*(exp(lambda232*t1)-exp(lambda232*t2)) 
                heProfile(i) = heProfile(i) + He_prod
            enddo
      enddo
      ! Store He profile of the geological model
      HeSphereEjec_stored(:) = heProfile(:)
      
      deallocate (heProfile)
  endif

  !########################################################################################
  !############################## Compute diffusivities ###################################
  !########################################################################################
  
  ! Radiation Damage models parameters
  if ((RDmodel.eq.2).or.(RDmodel.eq.3).or.(RDmodel.eq.4).or.(RDmodel.eq.6)) then
      
      ts = 1e6*365.25*24*3600 !time in s
      
      ntime_highres = MAX_NUM_TIME_STEPS
      allocate(time_out(MAX_NUM_TIME_STEPS), temp_out(MAX_NUM_TIME_STEPS))
      status = 1

      time_out = 0.d0
      temp_out = 0.d0
      ! Input time and temperature has to be in ascending order (from 0Ma to X Ma)
      ! a minimum timstep of 0.2 Ma is set
      time_reversed = time(ntime:1:-1)
      temperature_reversed = temperature(ntime:1:-1)
      call InterpolateTTPathKet(ntime, time_reversed, temperature_reversed, MAX_NUM_TIME_STEPS,&
                 time_out, temp_out, ntime_highres,Temp_crit, status)

      ! Maximum length of time-temperature vector is set by MAX_NUM_TIME_STEPS
      ! if less than that, reduce time_out to ntime_highres
      allocate(time_highres(ntime_highres), temp_highres(ntime_highres))

      ! reverse time-temperature back in descending order (from XMa to 0 Ma)
      time_highres = time_out(1:ntime_highres)
      temp_highres = temp_out(1:ntime_highres)

      allocate(rho_r(ntime_highres,ntime_highres),ageTime_highres(ntime_highres),dtf(ntime_highres))
      rho_r = 0.

      ! print *, 'time : ', time_highres
      ! print *, 'temperature: ', temp_highres
      ! print *, 'ntime: ', ntime

      ! get time step for each time interval
      dtf(2:ntime_highres) = abs(time_highres(2:ntime_highres)-time_highres(1:ntime_highres-1))
     
      call RD09(time_highres,temp_highres,ntime_highres,rho_r,dtf*ts,rmr0)
      
      deallocate(dtf)
      
  elseif (RDmodel.gt.6) then
      stop 'option not supported in Hediff.f90 - Radiation Damage model'
  else
      allocate(time_highres(ntime), temp_highres(ntime), ageTime_highres(ntime))
      time_highres = time
      temp_highres = temperature
      ntime_highres = ntime
  endif
  
  ! Geological model - Track Generation
  allocate(diffusivities(ntime_highres), rhov(ntime_highres))
  rhov = 0.
  do i=2,ntime_highres !time increases
    t2 = time_highres(i)*1e6
    t1 = time_highres(i-1)*1e6
    !convert to atoms/volume with apatite density of 3.19 g/cm3
    !and then make it proportional to each isotope
    rhov(i) = ((U238_array(1)*3.19)*(exp(lambda238*t1)-exp(lambda238*t2))) + &
        (7./8.*(U235_array(1)*3.19)*(exp(lambda235*t1)-exp(lambda235*t2))) + &
        (6./8.*(Th232_array(1)*3.19)*(exp(lambda232*t1)-exp(lambda232*t2)))
  enddo
  ! Get diffusivity
  diffusivities = 0.
  call RD_model(diffusivities,rhov,RDmodel,ntime_highres,D0,Ea,rho_r,temp_highres,time_highres,0,&
                  Uppm_array,Thppm_array,nrad,data_ADAM_array,grainsize,eU)
  ! open (134,file='diffusivities.csv',status='unknown')
  ! do kk=1,ntime_highres
  !    write (134,'(g12.6,3(",",g12.6))') time_highres(kk),temp_highres(kk),diffusivities(kk)/1e8,rhov(kk)
  ! enddo
  ! close(134)
  if ((RDmodel.eq.2).or.(RDmodel.eq.3).or.(RDmodel.eq.4).or.(RDmodel.eq.6)) then             
    deallocate (rho_r)
  endif

  !########################################################################
  !######################## Diffusion model ###############################
  !########################################################################
  allocate (diag(nrad), upper(nrad), lower(nrad), f(nrad), HeSphere(nrad))
  allocate (heProfile(nrad), integrand(nrad))
  diag = 0.
  f = 0.
  fejec = 0.
  HeSphere_prev = 100.
  counter = 0
  ! Loop through time
  do tstep=1,ntime_highres-1
    !---------------------- Geological model --------------------------------------
        t1 = time_highres(tstep)*1e6   !older time (yr)
        t2 = time_highres(tstep+1)*1e6 !younger time  
        dt0 = abs(t1-t2)*365.25*24*60*60 !seconds
        DT = diffusivities(tstep)
        beta = (2*dr**2)/(DT*dt0)
        upper = 1.
        lower = 1.
        do i=1,nrad
            diag(i) = -2-beta
        enddo
        HeSphere = 0.

        ! Production for the first node (atoms/g)
        He_prod = 8*aE238_rad(1)*(exp(lambda238*t1)-exp(lambda238*t2)) &
                + 7*aE235_rad(1)*(exp(lambda235*t1)-exp(lambda235*t2)) &
                + 6*aE232_rad(1)*(exp(lambda232*t1)-exp(lambda232*t2)) 
    
        ! Newman boundary condition
        diag(1) = -3-beta
        f(1) = (3-beta)*HeSphere_prev(1)-HeSphere_prev(2)-He_prod*(1-0.5)*dr*beta
    
        ! Dirichlet boundary condition
        He_prod = 8*aE238_rad(nrad)*(exp(lambda238*t1)-exp(lambda238*t2)) &
                + 7*aE235_rad(nrad)*(exp(lambda235*t1)-exp(lambda235*t2)) &
                + 6*aE232_rad(nrad)*(exp(lambda232*t1)-exp(lambda232*t2)) 
        f(nrad) = -HeSphere_prev(nrad-1)+(2-beta)*HeSphere_prev(nrad)-He_prod*(nrad-0.5)*dr*beta
    
        ! Compute for intermediate nodes
        do i=2,nrad-1
            He_prod = 8*aE238_rad(i)*(exp(lambda238*t1)-exp(lambda238*t2)) &
                    + 7*aE235_rad(i)*(exp(lambda235*t1)-exp(lambda235*t2)) &
                    + 6*aE232_rad(i)*(exp(lambda232*t1)-exp(lambda232*t2)) 
            f(i) = (2-beta)*HeSphere_prev(i)-HeSphere_prev(i-1)-HeSphere_prev(i+1)-He_prod*(i-0.5)*dr*beta
        enddo
    
        ! Tridiagonal matrix algorithm
        call tridag (lower,diag,upper,f,HeSphere,nrad)

        HeSphere_prev(:) = HeSphere(:)
    
        ! Convert sphere u array to He
        minHe = HeSphere(1)/(0.5*dr)
        heProfile = 0.
        integrand = 0.
        do i=1,nrad
            rad_pos = (i-0.5)*dr
            heProfile(i) = HeSphere(i)/rad_pos
            integrand(i) = heProfile(i)*4*pi*rad_pos**2
        enddo
    
        ! Store He profile of the geological model
        if (tstep.eq.ntime_highres-1) then
            HeSphere_stored(:) = heProfile(:)
        endif
        call Romberg_integration(integrand,grainsize,res,nrad)
    
        ! Fill in crystal center 
        totalHe = res !atoms/g
        totalHe = totalHe + 0.5*dr*(integrand(1)*55/24-integrand(2)*59/24 &
                + integrand(3)*37/24-integrand(4)*9/24)
            
        ! Convert to same basis as the other isotope totals
        totalHe = totalHe/(4*pi/3)
        if (minHe.lt.-heProfile(1)*0.05) then
            totalHe = 0.
        endif
        ! Compute He age
        ! iterative age calculation
        leftSum = total238*ft238 + total235*ft235 + total232*ft232 + totalHe
        loAge = 0.
        hiAge = time_highres(1)*1e6
        ageConv = 100
        do while ((hiAge-loAge).gt.ageConv)
            midAge = (hiAge+loAge)/2
            midVal = total238*ft238*exp(lambda238*midAge)+ total235*ft235*exp(lambda235*midAge) &
                    + total232*ft232*exp(lambda232*midAge)
            if (midVal.lt.leftSum) then
                loAge = midAge
            else
                hiAge = midAge
            endif
        enddo
        heModelAge = (hiAge+loAge)/2
        heModelAge = heModelAge/1e6
        ! Store ages in array
        ! if (Pecube.eq.0) then
        ageTime_highres(tstep+1) = heModelAge
    ! else
    Apatite_age = heModelAge
    ! endif
  enddo
  if (He_flag.eq.0.and.Pecube.eq.0) then
      print '(/,"     Age (Ma) = ",f5.2)', heModelAge
  endif
  
  ! interpolate high-res age vector to size of input time vector
  if ((RDmodel.eq.2).or.(RDmodel.eq.3).or.(RDmodel.eq.4).or.(RDmodel.eq.6)) then
    ageTime_dp = 0.d0
    call interpolate2_1D(time_highres,time,ntime_highres,ageTime_highres,ageTime_dp,ntime)
    ageTime = real(ageTime_dp, 4)
  else
    ageTime = ageTime_highres
  endif
  deallocate (diffusivities,ageTime_highres)
  
  
  !-------------  Check for 4He/3He prediction ---------------------------
if (He_flag.eq.1) then
  ! Get diffusivity
  allocate (diffusivities(nstep_He43),rhov_stored(nstep_He43),rho_r(nstep_He43,nstep_He43))
  allocate(dtf(nstep_He43))
  
  diffusivities = 0.
  rho_r = 0.
  rhov_stored = rhov(ntime_highres) ! Assumes no He generation during step heating
  dtf(2:nstep_He43) = abs(time_step_array(2:nstep_He43)-time_step_array(1:nstep_He43-1)) ! in seconds
  
  ! Compute damages annealing
  call RD09(time_step_array,temp_step_array,nstep_He43,rho_r,dtf,rmr0)
  
  ! Compute diffusivities
  call RD_model(diffusivities,rhov_stored,RDmodel,nstep_He43,D0,Ea,rho_r,temp_step_array,time_step_array/(3600*24*365.25),1,&
                  Uppm_array,Thppm_array,nrad,data_ADAM_array,grainsize,eU)
  
  ! Start with He profile from the geological model
  HeSphere_prev = HeSphere_stored(:)
  HeSphereEjec_prev = HeSphereEjec_stored(:)
  do i=1,nrad !convert to u
    rad_pos = (i-0.5)*dr
    HeSphere_prev(i) = HeSphere_prev(i)*rad_pos
    He3Sphere_prev(i) = 1e10*rad_pos !Uniform 3He profile
        HeSphereEjec_prev(i) = HeSphereEjec_prev(i)*rad_pos
  enddo
  ! Initial Total amount of 3He
  heProfile = 0.
  integrand = 0.
  do i=1,nrad
    rad_pos = (i-0.5)*dr
    heProfile(i) = He3Sphere_prev(i)/rad_pos
    integrand(i) = heProfile(i)*4*pi*rad_pos**2
  enddo
  res = 0.
  call Romberg_integration(integrand,grainsize,res,nrad)
  He3Sphere_init = res + 0.5*dr*(integrand(1)*55/24-integrand(2)*59/24 &
    + integrand(3)*37/24-integrand(4)*9/24)
  
  ! Number of steps to records (for 4He/3He)
  counterStep = 2 ! start at 2 to let first entry 0
  He4_frac_cum = 0.
  He3_frac_cum = 0.
  He3_released = 0.
  He3_frac_step = He3Sphere_init
  He4Ejec_frac_cum = 0.
  TotHe3_released_temp = 0.
  He_ratio_temp = 0.
  
  ! Loop through time
  do tstep=1,nstep_He43-1
    t2 = time_step_array(tstep)   !older time (s)
    t1 = time_step_array(tstep+1) !younger time
    dt0 = abs(t1-t2)! time step seconds

    ! 4He release fraction at the begining of the time step 
    He4_frac = 0.

    ! Convert sphere u array to He
    heProfile = 0.
    integrand = 0.
    do i=1,nrad
        rad_pos = (i-0.5)*dr
        heProfile(i) = HeSphere_prev(i)/rad_pos
        integrand(i) = heProfile(i)*4*pi*rad_pos**2
    enddo
    res = 0.
    call Romberg_integration(integrand,grainsize,res,nrad)
    totalHe = res
    He4_frac = totalHe + 0.5*dr*(integrand(1)*55/24-integrand(2)*59/24 &
        + integrand(3)*37/24-integrand(4)*9/24)
          
    ! 4He release fraction at the begining of the time step - Ejec only
    He4_frac_Ejec = 0.

    ! Convert sphere u array to He
    heProfile = 0.
    integrand = 0.
    do i=1,nrad
        rad_pos = (i-0.5)*dr
        heProfile(i) = HeSphereEjec_prev(i)/rad_pos
        integrand(i) = heProfile(i)*4*pi*rad_pos**2
    enddo
    res = 0.
    call Romberg_integration(integrand,grainsize,res,nrad)
    totalHe = res
    He4_frac_Ejec = totalHe + 0.5*dr*(integrand(1)*55/24-integrand(2)*59/24 &
        + integrand(3)*37/24-integrand(4)*9/24)
          
    ! 3He release fraction at the begining of the time step
    He3_frac = 0.

    ! Convert sphere u array to He
    he3Profile = 0.
    integrand3 = 0.
    do i=1,nrad
        rad_pos = (i-0.5)*dr
        he3Profile(i) = He3Sphere_prev(i)/rad_pos
        integrand3(i) = he3Profile(i)*4*pi*rad_pos**2
    enddo
    res = 0.
    call Romberg_integration(integrand3,grainsize,res,nrad)
    totalHe3 = res
    He3_frac = totalHe3 + 0.5*dr*(integrand3(1)*55/24-integrand3(2)*59/24 &
        + integrand3(3)*37/24-integrand3(4)*9/24)
    
    ! Get diffusivity
    DT = diffusivities(tstep)
    beta = (2*dr**2)/(DT*dt0)
    upper = 1.
    lower = 1.
    do i=1,nrad
        diag(i) = -2-beta
    enddo
    HeSphere = 0.
    He3Sphere = 0.
    HeSphereEjec = 0.
    
    ! Newman boundary condition
    diag(1) = -3-beta
    f(1) = (3-beta)*HeSphere_prev(1)-HeSphere_prev(2)
    if (He_flag.eq.1) then !Repeat for 3He and ejec only
        f3(1) = (3-beta)*He3Sphere_prev(1)-He3Sphere_prev(2)
        fejec(1) = (3-beta)*HeSphereEjec_prev(1)-HeSphereEjec_prev(2)
    endif
    
    ! Dirichlet boundary condition
    f(nrad) = -HeSphere_prev(nrad-1)+(2-beta)*HeSphere_prev(nrad)
    if (He_flag.eq.1) then !Repeat for 3He
        f3(nrad) = -He3Sphere_prev(nrad-1)+(2-beta)*He3Sphere_prev(nrad)
        fejec(nrad) = -HeSphereEjec_prev(nrad-1)+(2-beta)*HeSphereEjec_prev(nrad)
    endif   
    
    ! Compute for intermediate nodes
    do i=2,nrad-1
        f(i) = (2-beta)*HeSphere_prev(i)-HeSphere_prev(i-1)-HeSphere_prev(i+1)
        if (He_flag.eq.1) then !Repeat for 3He
            f3(i) = (2-beta)*He3Sphere_prev(i)-He3Sphere_prev(i-1)-He3Sphere_prev(i+1)
            fejec(i) = (2-beta)*HeSphereEjec_prev(i)-HeSphereEjec_prev(i-1)-HeSphereEjec_prev(i+1)
        endif
    enddo
    
    ! Tridiagonal matrix algorithm
    call tridag (lower,diag,upper,f,HeSphere,nrad)
    
    HeSphere_prev(:) = HeSphere(:)
    
    ! Convert sphere u array to He
    minHe = HeSphere(1)/(0.5*dr)
    heProfile = 0.
    integrand = 0.
    do i=1,nrad
        rad_pos = (i-0.5)*dr
        heProfile(i) = HeSphere(i)/rad_pos
        integrand(i) = heProfile(i)*4*pi*rad_pos**2
    enddo
    
    call Romberg_integration(integrand,grainsize,res,nrad)
    
    ! Fill in crystal center 
    totalHe = res !atoms/g
    totalHe = totalHe + 0.5*dr*(integrand(1)*55/24-integrand(2)*59/24 &
            + integrand(3)*37/24-integrand(4)*9/24)
    
    !---------------------------------
    ! Integrated fraction of 4He at the end of time step
    ! 4He fraction at the begining of the next time step
    He4_frac_next = totalHe
    
       ! Tridiagonal matrix algorithm - Ejec only
    call tridag (lower,diag,upper,fejec,HeSphereEjec,nrad)
    
    HeSphereEjec_prev(:) = HeSphereEjec(:)
    
    ! Convert sphere u array to He
    minHe = HeSphereEjec(1)/(0.5*dr)
    heProfile = 0.
    integrand = 0.
    do i=1,nrad
        rad_pos = (i-0.5)*dr
        heProfile(i) = HeSphereEjec(i)/rad_pos
        integrand(i) = heProfile(i)*4*pi*rad_pos**2
    enddo
    
    call Romberg_integration(integrand,grainsize,res,nrad)
    
    ! Fill in crystal center 
    totalHe = res !atoms/g
    totalHe = totalHe + 0.5*dr*(integrand(1)*55/24-integrand(2)*59/24 &
            + integrand(3)*37/24-integrand(4)*9/24)
    
    !---------------------------------
    ! Integrated fraction of 4He at the end of time step
    ! 4He fraction at the begining of the next time step
    He4_fracEjec_next = totalHe
    
    ! Tridiagonal matrix algorithm - He3
    call tridag (lower,diag,upper,f3,He3Sphere,nrad)
    
    He3Sphere_prev(:) = He3Sphere(:)
    
    ! Convert sphere u array to He
    he3Profile = 0.
    integrand3 = 0.
    do i=1,nrad
        rad_pos = (i-0.5)*dr
        he3Profile(i) = He3Sphere(i)/rad_pos
        integrand3(i) = he3Profile(i)*4*pi*rad_pos**2
    enddo
    res3 = 0.
    call Romberg_integration(integrand3,grainsize,res3,nrad)
    
    ! Fill in crystal center 
    totalHe3 = res3 !atoms/g
    He3_frac_next = totalHe3 + 0.5*dr*(integrand3(1)*55/24-integrand3(2)*59/24 &
            + integrand3(3)*37/24-integrand3(4)*9/24)
    

    if (He4_frac.lt.0) He4_frac = 1e-10
    if (He4_frac_next.lt.0) then
         He4_frac_next = 1e-10
         HeSphere_prev(:) = 1e-10
    endif
    if (He4_frac_Ejec.lt.0) He4_frac_Ejec = 1e-10
    if (He4_fracEjec_next.lt.0) He4_fracEjec_next = 1e-10
    if (He3_frac.lt.0) He3_frac = 1e-10
    if (He3_frac_next.lt.0) then
        He3_frac_next = 1e-10
        He3Sphere_prev(:) = 1e-10
    endif
    
    ! Calculate cumulative release fraction
    He4_frac_cum = He4_frac_cum + He4_frac - He4_frac_next
    He4Ejec_frac_cum = He4Ejec_frac_cum + He4_frac_Ejec - He4_fracEjec_next
    He3_frac_cum = He3_frac_cum + He3_frac-He3_frac_next
    
    TotHe3_released_temp(counterStep) = (He3Sphere_init - He3_frac_next) / He3Sphere_init
    He3_frac_step = He3_frac_next
    He_ratio_temp(counterStep) = (He4_frac_cum+1e-10) / (He3_frac_cum+1e-10)
    He3_released(counterStep) = He3_frac_cum
    He_ratioEjec_temp(counterStep) = (He4Ejec_frac_cum+1e-10) / (He3_frac_cum+1e-10)
    He4_frac_cum = 0.
    He3_frac_cum = 0.
    He4Ejec_frac_cum = 0.
        
    counterStep = counterStep+1

    enddo

    TotHe3_released_temp(1) = TotHe3_released_temp(2) 
    He_ratio_temp(1) = He_ratio_temp(2) 
    He3_released(1) = He3_released(2)
    He_ratioEjec_temp(1) = He_ratioEjec_temp(2) 
    
    deallocate (diffusivities,dtf)
  endif ! end time loop
    

  if (He_flag.eq.1) then
    bulk = sum(He_ratio_temp(1:counterStep)*He3_released(1:counterStep))/sum(He3_released(1:counterStep))
      ! print *, 'bulk = ', bulk,  sum(He_ratio_temp(:)*He3_released(:))
      ! print *, 'He_ratio_temp = ', He_ratio_temp
      ! print *, 'He3_released = ', He3_released
    He_ratio_temp(1:counterStep) = He_ratio_temp(1:counterStep)/bulk
    bulk = sum(He_ratioEjec_temp(1:counterStep)*He3_released(1:counterStep))/sum(He3_released(1:counterStep))
      He_ratioEjec_temp(1:counterStep) = He_ratioEjec_temp(1:counterStep)/bulk

      ! Interpolate back to the same SF3He released observed (if observations)
      if (Pecube.eq.1) then
      ! check first entry equal zero
      if (SF3Heobs(1).eq.0d0) then
        SF3Heobs(1) = 1e-6
      endif
      !  print *, 'step 1', ndur, nstep_He43
      !  print *, '4He/3he step pred: ', He_ratio_temp(1:counterStep)
      !  print *, 'sumF3He pred: ', TotHe3_released_temp(1:counterStep)
      !  print *, 'He_ratio pred: ', He_ratio
      !  print *, 'sumF3He obs: ', SF3Heobs
      !  print *, 'step 2', ndur, nstep_He43, TotHe3_released
      !  print *, 'step 3', ndur, nstep_He43
        call interp_1D(TotHe3_released_temp,SF3Heobs,nstep_He43,He_ratio_temp,He_ratio,ndur)
        call interp_1D(TotHe3_released_temp,SF3Heobs,nstep_He43,TotHe3_released_temp,TotHe3_released,ndur)
        call interp_1D(TotHe3_released_temp,SF3Heobs,nstep_He43,He_ratioEjec_temp,He_ratioEjec,ndur)
      
        ! write 4He/3He profiles predicted
        ! open (134,file='He43_predictions.csv',status='unknown',position='append',action='write')
        ! do kk=1,ndur
        !    write (134,'(g12.6,1(",",g12.6))') TotHe3_released(kk),He_ratio(kk)
        ! enddo
        ! close(134)

      !!!! Debug !!!!
      ! open (134,file='He43spectrum.csv',status='unknown')
      ! do kk=1,nstep_He43
      !    write (134,'(g12.6,3(",",g12.6))') time_step_array(kk)/(3600),temp_step_array(kk),TotHe3_released_temp(kk),He_ratio_temp(kk)
      ! enddo
    else
        He_ratio = He_ratio_temp(:)
        TotHe3_released = TotHe3_released_temp(:)
    endif   
        
      deallocate (time_step_array,temp_step_array,rho_r,rhov_stored)
  endif
  
 deallocate (radius, U238_array,U235_array,Th232_array)
 deallocate (Uppm_array,Thppm_array)
 deallocate (aE238_rad, aE232_rad,aE235_rad)
 deallocate (Eprod_U238, Eprod_U235, Eprod_Th232)
 deallocate (diag, upper, lower, f)
 deallocate (heProfile,integrand)
 deallocate (rhov,time_highres,temp_highres)
 if ((RDmodel.eq.2).or.(RDmodel.eq.3).or.(RDmodel.eq.4).or.(RDmodel.eq.6)) then
   deallocate (time_out,temp_out)
 endif

 
end subroutine Hediff




!----------------------------------------------------------------------------------------------------
subroutine Romberg_integration(integrand,grainsize,res,nrad)

    implicit none
    integer decdigs, maxSt, minSt, st,j,i,nrad
    real, dimension(:), allocatable:: romall
    real, dimension(:,:), allocatable::rom
    real h,sumRom,integrand(nrad), res
    real *4 grainsize
    decdigs = 10

    allocate (romall(nrad))
    allocate (rom(2,decdigs))
    
    !Romberg integration algorithms
    rom = 0.
    romall(:) = integrand(:)
    h = grainsize - 0
    rom(1,1) = h*(romall(1)+romall(nrad))/2
    do i=2,decdigs
        st = 2**(decdigs-i+1)
        maxSt = 2**(decdigs-1)  
        minSt = st/2 + 1
        sumRom = 0.
        do j=minSt,maxSt,st
            sumRom = sumRom + romall(j)
        enddo
        rom(2,1) = (rom(1,1)+h*sumRom)/2
        do j=1,i-1
            rom(2,j+1) = ((4**j)*rom(2,j)-rom(1,j))/((4**j)-1)
        enddo
        rom(1,1:i) = rom(2,1:i)
        h=h/2
    enddo
    res = rom(1,decdigs)
    
    deallocate (romall, rom)
    
end subroutine 




