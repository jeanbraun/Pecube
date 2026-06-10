!----------------------------------------------------------------------------------------------------
subroutine FissionTrackAge(time, ztemp,ntime, rho_red,ageAFT, Pecube,rhoST,rmr0,fdist)

 !Compute AFT ages of each grain
 !Rho_r is normalized fission-track density
 
 real rho_std, ageAFT,ts,rho_red(ntime)
 integer ntime, rr,Pecube
 double precision rho_r(ntime,ntime), time(ntime),dtf(ntime),time_interp(ntime),temp_interp(ntime)
 double precision ztemp(ntime), fdist(200),temp_interp2(ntime)
 double precision rmr0
 double precision rhoST
 
 !Estimated fission track density reduction for durango apatite
 !spontaneous mean track length = 14.47 µm
 !mean induced track length = 16.21 µm

 rho_red = 0.
 !Time step
 dtf(2:ntime) = abs(time(2:ntime)-time(1:ntime-1))
 ts = 1e6*365.25*24*3600 !time in s
 time_interp =0.
 temp_interp = 0.
 time_interp(1:ntime-1) = time(1:ntime-1)-dtf(2:ntime)/2

 !interpolate temperature onto grid
 call interpolate(ztemp,temp_interp,time,time_interp,ntime)
 temp_interp2(3:ntime) = temp_interp(2:ntime-1)
 temp_interp = temp_interp2
 temp_interp(2) = (ztemp(1)+ztemp(2))/2
 
 call RD09(time_interp,temp_interp,ntime,rho_r,dtf*ts,rmr0)

 do rr=2,ntime !See Ketcham (2005)
    rho_red(rr) = rho_red(rr-1) + sum(rho_r(2:rr,rr))/(rr-2+1)*dtf(rr)
 enddo

 !AFT age 
 ageAFT =  rho_red(ntime)/rhoST

 if (Pecube.eq.0) then
     write (*,*) 'AFT age (Ma) = ', ageAFT
 endif
 
 
end subroutine FissionTrackAge

!----------------------------------------------------------------------------------------------------
subroutine Compute_rmr0(rmr0, kinID)
 ! calculate rmr0 from Dpar, Cl, or OH values
 double precision rmr0
 integer kinID
 
 if (kinID.eq.0) then ! from Dpar
      if (rmr0.le.1.75) then
        rmr0 = 0.84
      elseif (rmr0.ge.4.58) then
        rmr0 = 0.0
      else 
        rmr0 = 1.0-exp(0.647*(rmr0-1.75)-1.834)
      endif
 elseif (kinID.eq.1) then ! from Cl apfu
    if (abs(rmr0-1.0).le.0.130) then
        rmr0 = 0.0
    else
        rmr0 = 1.0-exp(2.107*(1.0-abs(rmr0-1.0))-1.834);
    endif
 elseif (kinID.eq.2) then ! from OH apfu
    rmr0 = 0.84*(1.0-(1.0-abs(rmr0-1.0))**4.5);
 elseif (kinID.eq.3) then ! Cl (wt %)
    rmr0 = rmr0 *0.2978
    if (abs(rmr0-1.0).le.0.130) then
        rmr0 = 0.0
    else
        rmr0 = 1.0-exp(2.107*(1.0-abs(rmr0-1.0))-1.834);
    endif
elseif (kinID.eq.4) then ! from rmr0
    rmr0 = rmr0
endif

 end subroutine Compute_rmr0

!----------------------------------------------------------------------------------------------------
subroutine RD09(time, temperature, ntime, rho_r, dtf, rmr0)

 ! Damage annealing model based on Ketcham et al. (2007) 
 ! based on Greg Falco's script -- Berkeley Geochronology Center
 ! Adapted by Maxime Bernard -- 07 june 2021
 ! dtf in seconds

 implicit none
 integer ntime, i, j, nstep, k
 double precision time(ntime), temperature(ntime)
 double precision rho_r(ntime,ntime),dtf(ntime)
 real c0, c1, c2, c3, alpha, kappa, teq, temp, rp,last_rs
 double precision rmr0, dt, dt_Myr
 real,dimension(:,:),allocatable::rcb2,rc
 
 !Set coefficients
 c0 = 0.39528
 c1 = 0.01073
 c2 = -65.12969
 c3 = -7.91715
 alpha = 0.04672
 kappa = 1.04 - rmr0

 
 !Create damage array to hold tracks
 allocate(rcb2(ntime,ntime),rc(ntime,ntime))
 
 rcb2 = 0.
 rc = 0.
 
 ! annealing at the step of initial track production (main diagonals)
 do i=2,ntime
    temp = (temperature(i)+temperature(i-1))/2  + 273.15 ! take average temp within time step
    rcb2(i,i) = 1/( (c0 + c1* (log(dtf(i))-c2) / ((log(1./temp))-c3) )**(1./alpha) + 1 )
 enddo

 ! Row are tracks generations and colums time/temperature history
 do i=3,ntime !Move forward in time
    temp = (temperature(i)+temperature(i-1))/2  + 273.15 ! take average temp within time step
    do j=2,i-1 ! move through track population
        last_rs = rcb2(j,i-1)
        teq = exp(c2 + ((log(1/temp)-c3)/c1)*((((1/last_rs)-1)**alpha)-c0))
        rcb2(j,i) = 1/( (c0 + c1* (log(dtf(i)+teq)-c2) / ((log(1./temp))-c3) )**(1./alpha) + 1 )
        ! endif
    enddo
 enddo

 rc = ((rcb2 - rmr0)/(1-rmr0))**kappa

 !Suppress situations where rcB2 < rmr0
 !Presumably, this means total annealing, so rc = 0
 do i=1,ntime
    do j=1,ntime
        if (rcb2(i,j).lt.rmr0) then
            rc(i,j) = 0.
        endif
    enddo
  enddo

 !Convert to reduced spontaneous density
 rp = 0.5274435 !Max roots of polynomial 9.205x²-9.157x+2.269
 do i=1,ntime
    do j=1,ntime
        if (rc(i,j).ge.0.765) then
            rho_r(i,j) = 1.6*rc(i,j)-0.6
        else
            rho_r(i,j) = 9.205*rc(i,j)**2-9.157*rc(i,j)+2.269
        endif
        if (rc(i,j).lt.rp) then
            rho_r(i,j) = 0.
        endif
    enddo
 enddo
 deallocate(rc,rcb2)
 
end subroutine RD09