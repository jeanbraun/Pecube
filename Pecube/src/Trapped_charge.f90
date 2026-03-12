subroutine GOK_FAD_model(tim, temp, nstep, ddot, D0, a, Et, s, b,&
 rhop,imax, nnf_array,age_array, nobs, iobs)

! function to compute the relative trap content from a given time-temperature
! history (tim-temp) of length nstep, given a set of kinetic parameters 
! General order kinetic model with fading (GOK_FAD)
! see Guralnik et al. (2015)


implicit none

integer nstep, nrp, nobs, iobs, ntime
double precision tim(nstep), temp(nstep)

double precision, dimension(:), allocatable :: nNf, rprime, pr, inv_tauath, inv_tauth, xkd, xk, T
double precision, dimension(:), allocatable :: nN_temp,Df
double precision, dimension(:,:), allocatable :: nN
double precision Ma, ddot, D0, a, Et,s, b, rhop, agemax
double precision ddot_out, s_out, rhop_out
double precision kb, Hs, magic_ratio, dt, npr, imax
double precision Temperature, rate
double precision nnf_array(nstep,nobs), age_array(nstep,nobs)
double precision T1, T2, time1, time2, delta_time, delta_temp, dtemp_max, max_dt
integer i, j, nsub

Ma = 3600*24*365*1.d6
! print *, 'Doser = ',ddot
! print *, 'D0 = ',D0
! print *, 's = ',s
! print *, 'Et = ',Et
! print *, 'a = ',a
! print *, 'b = ', b
! print *, 'imax = ',imax
! print *, 'rhop = ',rhop

! Extract parameters
ddot_out = ddot/(1.d3*365.d0*24.d0*3600.d0)!  % Gy/ka to Gy/seconds
s_out = 10.d0**s
rhop_out = 10.d0**rhop
dtemp_max = 1.0d0 ! 1°C maximum temperature change per time step
max_dt = 0.01d0 ! minimum time step for gok
nsub = 1

! Define constants
kb = 8.617343d-5 ! Boltzman constant
Hs = 3.d15 !s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot_out/D0
nrp = 100

! Define rprime range and tau athermic
allocate (rprime(nrp), pr(nrp), inv_tauath(nrp), inv_tauth(nrp))

do i=1,nrp
  rprime(i)= 2.5d0*float(i-1)/(nrp-1)
enddo
pr = 3.d0*rprime**2.d0*exp(-rprime**3.d0)
npr = sum(pr)

inv_tauath = (Hs*exp(-(rhop_out**(-1.d0/3.d0))*rprime))

allocate (T(nstep), nNf(nstep), Df(nstep))
T = temp+273.15d0

! computes nN for a Tt path
allocate (nN(nrp,nstep), nN_temp(nrp), xk(nrp), xkd(nrp))
nN=0.d0
nNf=0.d0
nN_temp =0.d0
dt = 0.d0
do i = 2,nstep
  T1 = T(i-1)
  T2 = T(i)
  time1 = tim(i-1)
  time2 = tim(i)
  delta_time = time1 - time2
  rate = (T2 - T1) / (delta_time * Ma)

  if (temp(i-1).lt.300) then ! We decrease time step for temperature lower than 200°C

    ! first max time step
    if (delta_time.gt.max_dt) then
        ntime = max(1, ceiling(delta_time/ max_dt))
        dt = (delta_time / dble(ntime)) * Ma ! time step in seconds
        
        ! then check for maximum temperature change within one time step
        delta_temp = abs(T2 - T1) / dble(ntime)
        nsub = max(1, ceiling(delta_temp/ dtemp_max))
        ntime = ntime * nsub
        dt = dt / nsub ! time step in seconds
    endif

  else
        dt = delta_time * Ma
        ntime = 1
  endif

  nN_temp = nN(:,i-1)
  do j = 1,ntime
    Temperature = T1 + rate * (dt*(j-0.5)) !Mid-point
    ! Calculate inverse thermal lifetime
    inv_tauth = s_out*exp(-Et/(kb*Temperature))
    ! Derivatives for integration
    xkd = -a*magic_ratio*(1.d0-nN_temp)**(a-1)-b*(nN_temp)**(b-1)*inv_tauth-inv_tauath
    xk = magic_ratio*(1.d0-nN_temp)**a-(nN_temp)**b*inv_tauth-inv_tauath*nN_temp
    nN_temp = nN_temp+xk*dt/(1.-dt*xkd)
  enddo
  nN(:,i) = nN_temp
  nNf(i) = sum(nN(:,i)*pr)
end do
nNf = nNf/npr
   
  ! Calculate equivalent dose
  Df = 0.d0
  do i = 2,nstep
    if (nNf(i).ge.imax) then
        Df(i) = abs( D0 * log(1 - 0.9999))
    else
        Df(i) = abs( D0 * log(1 - nNf(i)/imax))
    endif
  enddo

  ! Save results to output arrays
  nnf_array(:,iobs) = nNf
  agemax = (2*D0)/ddot/1000 ! Max age in Ma

  do i = 2,nstep
    age_array(i,iobs) = (Df(i)/ddot_out)/Ma
    if (age_array(i,iobs).ge.agemax) then
        age_array(i,iobs) = agemax
    endif
  enddo

! Cleanup
deallocate (rprime, pr, inv_tauath, inv_tauth, T, nN, xk, xkd, nNf,nN_temp, Df)

return

end subroutine GOK_FAD_model


!#########################################################################################
subroutine GOK_Gauss_FAD_model(tim, temp, nstep, age_array,ddot, D0, s, Et, sigmaEt,rhop,&
                            a,imax,nnf_array, nobs, iobs)

! function to compute the relative trap content from a given time-temperature
! history (tim-temp) of length nstep, given a set of kinetic parameters 
! trapping GOK Gauss model with fading (athermal detrapping) 
implicit none

integer nstep, nrp, nobs, iobs, ntime
double precision tim(nstep), temp(nstep), inv_tauth

double precision, dimension(:), allocatable :: nNf, T,Df
double precision, dimension(:), allocatable :: Ea, pEa,rprime, pr, inv_tauath,xkd, xk
double precision, dimension(:), allocatable :: nN_temp2
double precision, dimension(:,:), allocatable :: nN2, nN_temp
double precision, dimension(:,:,:), allocatable :: nN
double precision Temperature, rate, nN_temp3

double precision Ma, ddot, D0, Et, s, sigmaEt, rhop, agemax
double precision ddot_out, s_out, rhop_out, a, Hs, imax
double precision kb, magic_ratio, dt, pi, npEa, npr
double precision nnf_array(nstep,nobs), age_array(nstep,nobs)
integer i,j,k, nEa,istep, nsub
double precision T1, T2, time1, time2, delta_time, delta_temp, dtemp_max,max_dt

Ma = 3600*24*365*1.d6
! Extract parameters
ddot_out = ddot*1000/Ma ! Gy/s
!ddot = params(1)
s_out = 10.d0**s
rhop_out = 10.d0**rhop
if (rhop_out.lt.0) rhop_out = 0.d0
dtemp_max = 1.0d0 ! 1°C maximum temperature change per time step
max_dt = 0.01 ! maximum time step is 10 000 years
nsub = 1
! print *, 'Doser = ',ddot, ddot_out 
! print *, 'D0 = ',D0
! print *, 's = ',s, s_out
! print *, 'Et = ',Et
! print *, 'sigmaEt = ',sigmaEt
! print *, 'a = ',a
! print *, 'rhop = ',rhop
! print *, 'imax = ', imax
! Define constants
kb = 8.617343d-5 
Hs = 3.d15 !s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot_out/D0
nrp = 100

! Define rprime range and tau athermic
allocate (rprime(nrp), pr(nrp), inv_tauath(nrp))
  do i=1,nrp
  rprime(i)= 2.5d0*float(i-1)/(nrp-1)
  enddo
pr = 3.d0*rprime**2.d0*exp(-rprime**3.d0)
npr = sum(pr)

inv_tauath = (Hs*exp(-(rhop_out**(-1.d0/3.d0))*rprime))

! Create variables for GAUSS model Lambert et al (Submitted)
nEa = 100
allocate (Ea(nEa), pEa(nEa))
Ea =  (/((5.d0*i)/nEa, i=1,nEa)/)
pi = 4.d0*atan(1.d0)
pEa = exp(-0.5d0*((Ea-Et)/sigmaEt)**2)/(sigmaEt*sqrt(2.d0*pi))
npEa = sum(pEa)

allocate (T(nstep), nNf(nstep), Df(nstep))
T = temp+273.15d0


! computes nN for the random Tt path
allocate (nN(nrp,nEa,nstep), nN2(nrp,nstep), xk(nrp), xkd(nrp))
allocate (nN_temp(nrp,nEa), nN_temp2(nrp))
nN=0.d0
nNf=0.d0
nN_temp=0.d0
nN_temp3 = 0.d0
dt = 0.d0
do i = 2,nstep
  T1 = T(i-1)
  T2 = T(i)
  time1 = tim(i-1)
  time2 = tim(i)
  delta_time = time1 - time2
  rate = (T2 - T1) / (delta_time * Ma)

  if (temp(i-1).lt.250) then ! We decrease time step for temperature lower than 200°C

    ! first max time step
    if (delta_time.gt.max_dt) then
        ntime = max(1, ceiling(delta_time/ max_dt))
        dt = (delta_time / dble(ntime)) * Ma ! time step in seconds
        
        ! then check for maximum temperature change within one time step
        delta_temp = abs(T2 - T1) / dble(ntime)
        nsub = max(1, ceiling(delta_temp/ dtemp_max))
        ntime = ntime * nsub
        dt = dt / nsub ! time step in seconds
    endif

  else
        dt = delta_time * Ma
        ntime = 1
  endif

    nN_temp = nN(:,:,i-1)
    nN_temp2 = 0.d0
    do istep = 1,ntime
      Temperature = T(i-1) + rate * (dt*(istep-0.5)) !Mid-point
      do j=1,nEa
        inv_tauth = s_out*exp(-Ea(j)/(kb*Temperature))
        xkd = -a*magic_ratio*(1.d0-nN_temp(:,j))**(a-1)-inv_tauth-inv_tauath
        xk = magic_ratio*(1.d0-nN_temp(:,j))**a-(nN_temp(:,j))*inv_tauth-inv_tauath*nN_temp(:,j)
        nN_temp(:,j) = (nN_temp(:,j)+xk*dt/(1.-dt*xkd))
    enddo

    do k=1,nrp
        nN_temp2(k) = sum(nN_temp(k,:)*pEa)/npEa
    enddo

    nN_temp3 =  sum(nN_temp2(:)*pr)/npr

    enddo
    nN(:,:,i) = nN_temp
    nNf(i) = nN_temp3
end do

Df = 0.d0
do i = 2,nstep
    if (nNf(i).ge.imax) then
        Df(i) = abs( D0 * log(1 - 0.9999))
    else
        Df(i) = abs( D0 * log(1 - nNf(i)/imax))
    endif
enddo


nnf_array(:,iobs) = nNf

! The age is modelled from the calculated n/N value and by calculating an equivalent dose df, this is to remain as close as possible to the definition of observed ages
 agemax = (2*D0)/ddot/1000 ! in Ma
 do i = 2,nstep
    age_array(i,iobs) = (Df(i)/ddot_out)/Ma
    if (age_array(i,iobs).ge.agemax) then
        age_array(i,iobs) = agemax
    endif
enddo


deallocate (rprime, pr, inv_tauath, Ea, pEa, T, nN, nNf,xkd, xk, Df)
deallocate (nN_temp,nN_temp2)

return
end subroutine GOK_Gauss_FAD_model




!#########################################################################################
subroutine SSE_Gauss_FAD_model(tim, temp, nstep, age_array,ddot, D0, s, Et, sigmaEt,rhop,&
                            imax,nnf_array, nobs, iobs)

! function to compute the relative trap content from a given time-temperature
! history (tim-temp) of length nstep, given a set of kinetic parameters 
! The rate equation uses the single saturating exponential function and 
! the gaussian model for detrapping (Lambert, 2018), this also include
! fading term as described in Kars et al., 2008.
                    

implicit none

integer nstep, nrp, nobs, iobs, ntime
double precision tim(nstep), temp(nstep), inv_tauth

double precision, dimension(:), allocatable :: nNf, T,Df,alpha1
double precision, dimension(:), allocatable :: Ea, pEa,rprime, pr, inv_tauath,xkd, xk
double precision, dimension(:), allocatable :: nN_temp2
double precision, dimension(:,:), allocatable :: nN2, nN_temp
double precision, dimension(:,:,:), allocatable :: nN
double precision Temperature, rate, nN_temp3

double precision Ma, ddot, D0, Et, s, sigmaEt, rhop, agemax
double precision ddot_out, s_out, rhop_out, Hs, imax
double precision kb, magic_ratio, dt, pi, npEa, npr
double precision nnf_array(nstep,nobs), age_array(nstep,nobs)
integer i,j,k, nEa,istep,nsub
double precision T1, T2, time1, time2, delta_time, delta_temp, dtemp_max, max_dt

Ma = 3600*24*365*1.d6
! Extract parameters
ddot_out = ddot*1000/Ma ! Gy/s
!ddot = params(1)
s_out = 10.d0**s
rhop_out = 10.d0**rhop
if (rhop_out.lt.0) rhop_out = 0.d0
dtemp_max = 1.0d0 ! 1°C maximum temperature change per time step
max_dt = 0.01 ! maximum time step is 10 000 years
nsub = 1
! Define constants
kb = 8.617343d-5 
Hs = 3.d15 !s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot_out/D0
nrp = 100

! Define rprime range and tau athermic
allocate (rprime(nrp), pr(nrp), inv_tauath(nrp))
  do i=1,nrp
  rprime(i)= 2.5d0*float(i-1)/(nrp-1)
  enddo
pr = 3.d0*rprime**2.d0*exp(-rprime**3.d0)
npr = sum(pr)

inv_tauath = (Hs*exp(-(rhop_out**(-1.d0/3.d0))*rprime))

! Create variables for GAUSS model Lambert et al (Submitted)
nEa = 100
allocate (Ea(nEa), pEa(nEa))
Ea =  (/((5.d0*i)/nEa, i=1,nEa)/)
pi = 4.d0*atan(1.d0)
pEa = exp(-0.5d0*((Ea-Et)/sigmaEt)**2)/(sigmaEt*sqrt(2.d0*pi))
npEa = sum(pEa)

allocate (T(nstep), nNf(nstep), Df(nstep))
T = temp+273.15d0


! computes nN for the random Tt path
allocate (nN(nrp,nEa,nstep), nN2(nrp,nstep), xk(nrp), xkd(nrp))
allocate (nN_temp(nrp,nEa), nN_temp2(nrp),alpha1(nrp))
nN=0.d0
nNf=0.d0
nN_temp=0.d0
nN_temp3 = 0.d0
alpha1 = 0.d0
dt = 0.d0
do i = 2,nstep
  T1 = T(i-1)
  T2 = T(i)
  time1 = tim(i-1)
  time2 = tim(i)
  delta_time = time1 - time2
  rate = (T2 - T1) / (delta_time * Ma)

  if (temp(i-1).lt.150) then ! We decrease time step for temperature lower than 200°C

    ! first max time step
    if (delta_time.gt.max_dt) then
        ntime = max(1, ceiling(delta_time/ max_dt))
        dt = (delta_time / dble(ntime)) * Ma ! time step in seconds
        
        ! then check for maximum temperature change within one time step
        delta_temp = abs(T2 - T1) / dble(ntime)
        nsub = max(1, ceiling(delta_temp/ dtemp_max))
        ntime = ntime * nsub
        dt = dt / nsub ! time step in seconds
    endif

  else
        dt = delta_time * Ma
        ntime = 1
  endif
    
  nN_temp = nN(:,:,i-1)
  nN_temp2 = 0.d0
  do istep = 1,ntime
    Temperature = T(i-1) + rate * (dt*(istep-0.5)) !Mid-point
    do j=1,nEa
        inv_tauth = s_out*exp(-Ea(j)/(kb*Temperature))
        alpha1=magic_ratio+inv_tauth+inv_tauath
        nN_temp(:,j) = (nN_temp(:,j)+magic_ratio*dt)*(1/(1+dt*alpha1))
    enddo

    do k=1,nrp
        nN_temp2(k) = sum(nN_temp(k,:)*pEa)/npEa
    enddo

    nN_temp3 =  sum(nN_temp2(:)*pr)/npr

  enddo
  nN(:,:,i) = nN_temp
  nNf(i) = nN_temp3
end do

Df = 0.d0
do i = 2,nstep
    if (nNf(i).ge.imax) then
        Df(i) = abs( D0 * log(1 - 0.9999))
    else
        Df(i) = abs( D0 * log(1 - nNf(i)/imax))
    endif
enddo


nnf_array(:,iobs) = nNf

! The age is modelled from the calculated n/N value and by calculating an equivalent dose df, this is to remain as close as possible to the definition of observed ages
 agemax = (2*D0)/ddot/1000 ! in Ma
 do i = 2,nstep
    age_array(i,iobs) = (Df(i)/ddot_out)/Ma
    if (age_array(i,iobs).ge.agemax) then
        age_array(i,iobs) = agemax
    endif
enddo


deallocate (rprime, pr, inv_tauath, Ea, pEa, T, nN, nNf,xkd, xk, Df)
deallocate (nN_temp,nN_temp2)

return
end subroutine SSE_Gauss_FAD_model



!#########################################################################################
subroutine SSE_BTS_FAD_model(tim, temp, nstep,age_array, ddot, D0, Et, s,&
 rhop, Eu, nnf_array,imax, nobs, iobs)

! function to compute the relative trap content from a given time-temperature
! history (tim-temp) of length nstep, given a set of kinetic parameters 
! the rate equation uses the single saturating exponential growth and
! the band-tail state model for decay (Li and Li (2013))

implicit none

integer nstep, nrp,neb, nobs, iobs
double precision tim(nstep), temp(nstep)

double precision, dimension(:), allocatable :: nNf, rprime, pr, T, eb,peb,Df
double precision, dimension(:), allocatable :: inv_tauath,inv_tauth, alpha1
double precision, dimension(:), allocatable :: nN_temp2
double precision, dimension(:,:), allocatable :: nN_temp
double precision, dimension(:,:,:), allocatable :: nN
double precision Temperature, rate

double precision Ma, ddot, D0, Et,s, rhop, Eu
double precision ddot_out, s_out, rhop_out, imax, agemax
double precision kb, Hs, magic_ratio, dt, npr,npeb,inv_tauth1
double precision nnf_array(nstep,nobs), age_array(nstep,nobs)
integer i,j,istep,ntime,nsub
double precision T1, T2, time1, time2, delta_time, delta_temp, dtemp_max, max_dt

Ma = 3600*24*365*1.d6
! Extract parameters
ddot_out = ddot*1000/Ma ! Gy/ka to Gy/Ma (in seconds)
s_out = 10.d0**s
rhop_out = 10.d0**rhop
if (rhop_out.lt.0) rhop_out = 0.d0
dtemp_max = 1.0d0 ! 1°C maximum temperature change per time step
max_dt = 0.01 ! maximum time step is 10 000 years
nsub = 1
! print *, 'Doser = ',ddot, ddot_out 
! print *, 'D0 = ',D0
! print *, 's = ',s, s_out
! print *, 'Et = ',Et
! print *, 'Eu = ',Eu
! print *, 'rhop = ',rhop, rhop_out
! print *, 'imax = ',imax

! Define constants
kb = 8.617343d-5 ! Boltzman constant
Hs = 3.d15 !s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot_out/D0
nrp = 100
neb = 100

! Define rprime range and tau athermic
allocate (rprime(nrp), pr(nrp), inv_tauath(nrp), inv_tauth(nrp),eb(neb))
allocate (peb(neb))
do i=1,nrp
  rprime(i)= 2.5d0*float(i-1)/(nrp-1)
  inv_tauath(i) = Hs*exp(-(rhop_out**-(1./3))*rprime(i))
enddo
pr = 3.d0*rprime**2.d0*exp(-rprime**3.d0)
npr = sum(pr)

! Distribution of band-tail states
do j=1,neb
 eb(j)=Et*(j-1)/float(neb-1)
 peb(j)=exp(-eb(j)/Eu)
enddo
npeb=sum(peb)

allocate (T(nstep), nNf(nstep),Df(nstep))
T = temp+273.15d0

! computes nN for a Tt path
allocate (nN(neb,nrp,nstep),nN_temp2(neb),nN_temp(neb,nrp),alpha1(nrp))
nN=0.d0
nNf=0.d0
nN_temp = 0.d0
nN_temp2 = 0.d0
dt = 0.d0
do i = 2,nstep

  T1 = T(i-1)
  T2 = T(i)
  time1 = tim(i-1)
  time2 = tim(i)
  delta_time = time1 - time2
  rate = (T2 - T1) / (delta_time * Ma)

  if (temp(i-1).lt.150) then ! We decrease time step for temperature lower than 200°C

    ! first max time step
    if (delta_time.gt.max_dt) then
        ntime = max(1, ceiling(delta_time/ max_dt))
        dt = (delta_time / dble(ntime)) * Ma ! time step in seconds
        
        ! then check for maximum temperature change within one time step
        delta_temp = abs(T2 - T1) / dble(ntime)
        nsub = max(1, ceiling(delta_temp/ dtemp_max))
        ntime = ntime * nsub
        dt = dt / nsub ! time step in seconds
    endif

  else
        dt = delta_time * Ma
        ntime = 1
  endif

    nN_temp = nN(:,:,i-1)
    do istep = 1, ntime
        Temperature = T(i-1) + rate * (dt*(istep-0.5)) ! Mid-point
        do j = 1,neb
            inv_tauth1 = s_out*exp(-(Et-eb(j))/(kb*Temperature));
            alpha1=magic_ratio+inv_tauth1+inv_tauath;
            nN_temp(j,:) = (nN_temp(j,:)+magic_ratio*dt)*(1/(1+dt*alpha1));
            nN_temp2(j) = sum(nN_temp(j,:) * pr);
        enddo
    enddo
    nN_temp2 = peb * nN_temp2;
    nN(:,:,i) = nN_temp;
    nNf(i) = sum(nN_temp2)
enddo

nNf = nNf/(npeb*npr)
Df = 0.d0
do i = 2,nstep
    if (nNf(i).ge.imax) then
        Df(i) = abs( D0 * log(1 - 0.9999))
    else
        Df(i) = abs( D0 * log(1 - nNf(i)/imax))
    endif
enddo

nnf_array(:,iobs) = nNf

! The age is modelled from the calculated n/N value and by calculating an equivalent dose df, this is to remain as close as possible to the definition of observed ages
 agemax = (2*D0)/ddot/1000 ! in Ma
 do i = 2,nstep
    age_array(i,iobs) = (Df(i)/ddot_out)/Ma
    if (age_array(i,iobs).ge.agemax) then
        age_array(i,iobs) = agemax
    endif
enddo


deallocate (rprime, pr, inv_tauath, inv_tauth, T, nN, peb,eb,nNf,Df)
deallocate (nN_temp2,nN_temp)

return
end subroutine SSE_BTS_FAD_model




!#########################################################################################
subroutine SSE_Gauss_model(tim, temp, nstep, ddot, D0, s, Et, sigmaEt,&
imax,nnf_array, ageESR_array, nobs, iobs)

! function to compute the relative trap content from a given time-temperature
! history (tim-temp) of length nstep, given a set of parameters provided by TL
! experiment (code loosely translated from Rabiul's matlab code)
! trapping_Gauss_Fad (but no FAD for ESR) for Single Saturation Exponential (SSE) model

implicit none

integer nstep, nobs, iobs
double precision tim(nstep), temp(nstep)

double precision, dimension(:), allocatable :: nNf, T,Df
double precision, dimension(:), allocatable :: Ea, pEa, inv_tauth, alpha, alpha2
double precision, dimension(:), allocatable :: nN_temp
double precision, dimension(:,:), allocatable :: nN, D_array
double precision Temperature, rate

double precision Ma, ddot, D0, Et, s, sigmaEt, agemax
double precision ddot_out, s_out, imax
double precision kb, magic_ratio, dEa, dt, pi, npEa
double precision nnf_array(nstep,nobs), ageESR_array(nstep,nobs)
double precision T1, T2, time1, time2, delta_time, delta_temp, dtemp_max, max_dt
integer i, nEa,istep,ntime,nsub

Ma = 3600*24*365*1.d6
! Extract parameters
ddot_out = ddot*1000/Ma ! Gy/s
!ddot = params(1)
s_out = 10.d0**s
dtemp_max = 1.0d0 ! 1°C maximum temperature change per time step
max_dt = 0.01 ! maximum time step is 10 000 years
nsub = 1
! print *, 'Doser = ',ddot, ddot_out 
! print *, 'D0 = ',D0
! print *, 's = ',s, s_out
! print *, 'Et = ',Et
! print *, 'sigmaEt = ',sigmaEt
! print *, 'imax = ',imax

! Define constants
kb = 8.617343d-5 
magic_ratio = ddot_out/D0
dEa = 0.01d0 

! Create variables for GAUSS model Lambert et al (Submitted)
nEa = 100
allocate (Ea(nEa), pEa(nEa))
Ea =  (/((5.d0*i)/nEa, i=1,nEa)/)
nEa = size(Ea)
pi = 4.d0*atan(1.d0)
pEa = exp(-0.5d0*((Ea-Et)/sigmaEt)**2)/(sigmaEt*sqrt(2.d0*pi))
npEa = sum(pEa)

allocate (T(nstep), nNf(nstep), Df(nstep))
T = temp+273.15d0

! computes nN for the random Tt path
allocate (nN(nea,nstep),inv_tauth(nEa),alpha(nEa),alpha2(nEa),D_array(nea,nstep))
allocate (nN_temp(nEa))
nN=0.d0
nNf=0.d0
nN_temp = 0.d0
dt = 0.d0
do i = 2,nstep
  T1 = T(i-1)
  T2 = T(i)
  time1 = tim(i-1)
  time2 = tim(i)
  delta_time = time1 - time2
 rate = (T2 - T1) / (delta_time * Ma)

  if (temp(i-1).lt.250) then ! We decrease time step for temperature lower than 200°C

    ! first max time step
    if (delta_time.gt.max_dt) then
        ntime = max(1, ceiling(delta_time/ max_dt))
        dt = (delta_time / dble(ntime)) * Ma ! time step in seconds
        
        ! then check for maximum temperature change within one time step
        delta_temp = abs(T2 - T1) / dble(ntime)
        nsub = max(1, ceiling(delta_temp/ dtemp_max))
        ntime = ntime * nsub
        dt = dt / nsub ! time step in seconds
    endif

  else
        dt = delta_time * Ma
        ntime = 1
  endif

    nN_temp = nN(:,i-1)
    do istep = 1, ntime
        Temperature = T1 + rate * (dt*(istep-0.5)) !Mid-point
        inv_tauth = s_out*exp(-Ea/(kb*Temperature))
        alpha = magic_ratio + inv_tauth
        ! alpha2 = ddot_out + D0*inv_tauth
        nN_temp = (nN_temp+magic_ratio*dt)*(1.d0/(1.d0+dt*alpha))
        ! D_array(:,i) = (D_array(:,i-1)+dt*ddot_out)*((dt*ddot_out)/((dt*ddot_out)+dt*alpha))
        ! D_array(:,i) = (D_array(:,i-1)+dt*ddot_out)*(1.d0/(1.d0+dt*alpha))
        ! Df(i) = dot_product(pEa,D_array(:,i))
    enddo
    nN(:,i) = nN_temp
    nNf(i) = dot_product(pEa, nN_temp)
end do

nNf = nNf/npEa
! Calculate equivalent dose
Df = 0.d0
do i = 2,nstep
    if (nNf(i).ge.imax) then
        Df(i) = abs( D0 * log(1 - 0.9999))
    else
        Df(i) = abs( D0 * log(1 - nNf(i)/imax))
    endif
enddo
 ! print *, 'nnf = ', nNf(nstep), imax, D0, Df(nstep), ddot_out, npEa
! print *, 'nnf = ', nNf
! ESRModel = nNf(nstep)
 nnf_array(:,iobs) = nNf
 agemax = (2*D0)/ddot/1000 ! in Ma
 do i = 2,nstep
    ageESR_array(i,iobs) = (Df(i)/ddot_out)/Ma
    if (ageESR_array(i,iobs).ge.agemax) then
        ageESR_array(i,iobs) = agemax
    endif
enddo
! ageESR_array(:,iobs) = -1/(ddot/D0) * log(1 - nNf) / 1000 ! Ma
deallocate (Ea, pEa, T, nN, nNf, inv_tauth, alpha, alpha2,Df,D_array)

return
end subroutine SSE_Gauss_model


!#######################################################################################################################
subroutine GOK_model(tim, temp, nstep, ddot, D0, a, Et, s, b,&
 imax, nnf_array, ageESR_array, nobs, iobs)

! function to compute the relative trap content from a given time-temperature
! history (tim-temp) of length nstep, given a set of kinetic parameters 
! The rate equation uses the general order kinetic model (Guralnik et al., 2015)
! without fading.

implicit none

! Inputs
integer nstep, nobs, iobs
double precision tim(nstep), temp(nstep)
double precision ddot, D0, a, Et, s, b, imax

! Outputs
double precision nnf_array(nstep,nobs), ageESR_array(nstep,nobs)

! Local parameters
double precision, dimension(:), allocatable :: nNf, T
double precision, dimension(:), allocatable :: Df, nN
double precision Temperature, rate, nN_temp, inv_tauth, xkd, xk

double precision Ma, Ma_inv, agemax
double precision ddot_out, s_out
double precision kb, magic_ratio, dt

integer i, istep, ntime, nsub
double precision T1, T2, time1, time2, delta_time, delta_temp, dtemp_max,max_dt

! Constants
kb = 8.617343d-5 ! Boltzman constant eV/K
Ma = 3600.d0*24.d0*365.d0*1.d6 ! Million year to second
Ma_inv = 1.d0 / Ma
dtemp_max = 1.0d0 ! 1°C maximum temperature change per time step
max_dt = 0.01d0 ! max time step for GOK
nsub = 1
! print *, 'Doser = ',ddot, ddot_out 
! print *, 'D0 = ',D0
! print *, 's = ',s, s_out
! print *, 'Et = ',Et
! print *, 'a = ',a
! print *, 'b = ', b
! print *, 'imax = ',imax

! Extract parameters
ddot_out = ddot/(1.d3*365.d0*24.d0*3600.d0) ! Gy/ka to Gy/seconds
s_out = 10.d0**s
magic_ratio = ddot_out/D0

! Allocate arrays
allocate (T(nstep), nNf(nstep), Df(nstep))
allocate (nN(nstep))

! Convert temperature to Kelvin
T = temp+273.15d0
! Initialize arrays
nN=0.d0
nNf=0.d0
nN_temp = 0.d0
dt = 0.d0
! Time integration
do i = 2,nstep
  T1 = T(i-1)
  T2 = T(i)
  time1 = tim(i-1)
  time2 = tim(i)
  delta_time = time1 - time2
  delta_temp = abs(T2 - T1)
  rate = (T2 - T1) / (delta_time * Ma)

  if (temp(i-1).lt.300) then ! We decrease time step for temperature lower than 200°C

    ! first max time step
    if (delta_time.gt.max_dt) then
        ntime = max(1, ceiling(delta_time/ max_dt))
        dt = (delta_time / dble(ntime)) * Ma ! time step in seconds
        
        ! then check for maximum temperature change within one time step
        delta_temp = abs(T2 - T1) / dble(ntime)
        nsub = max(1, ceiling(delta_temp/ dtemp_max))
        ntime = ntime * nsub
        dt = dt / nsub ! time step in seconds
    endif

  else
        dt = delta_time * Ma
        ntime = 1
  endif

    ! Start with previous step nN values
    nN_temp = nN(i-1)

    ! Sub-time stepping loop
    do istep = 1, ntime
        Temperature = T1 + rate * (dt*(istep-0.5)) !Mid-point
        ! Calculate inverse thermal lifetime
        inv_tauth = s_out * exp(-Et / (kb*Temperature))
        ! Derivatives for integration
        xkd = -a*magic_ratio*(1.d0-nN_temp)**(a-1)-b*(nN_temp)**(b-1)*inv_tauth
        xk = magic_ratio*(1.d0-nN_temp)**a - (nN_temp)**b*inv_tauth
        ! Update trap occupancy using implicit Euler method
        nN_temp = nN_temp + xk*dt / (1.-dt*xkd)
    enddo
    ! Save current trap occupancy
    nN(i) = nN_temp
    nNf(i) = nN(i)
    ! print *, 'nN = ',nN(i), Temperature

  end do

! Calculate equivalent dose

Df = 0.d0
do i = 2,nstep
    if (nNf(i).ge.imax) then
        Df(i) = abs( D0 * log(1 - 0.9999))
    else
        Df(i) = abs( D0 * log(1 - nNf(i)/imax))
    endif
enddo

 ! Save results to output arrays
 nnf_array(:,iobs) = nNf
 agemax = (2*D0)/ddot/1000 ! Max age in Ma

 do i = 2,nstep
    ageESR_array(i,iobs) = (Df(i)/ddot_out)/Ma
    if (ageESR_array(i,iobs).ge.agemax) then
        ageESR_array(i,iobs) = agemax
    endif
enddo

! Cleanup
deallocate ( T, nN, nNf, Df)

return

end subroutine GOK_model