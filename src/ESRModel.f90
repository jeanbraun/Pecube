double precision function ESRModel (tim, temp, nstep, params)

! function to compute the relative trap content from a given time-temperature
! history (tim-temp) of length nstep, given a set of parameters provided by TL
! experiment (code loosely translated from Rabiul's matlab code)

implicit none

integer nstep, nrp,neb
double precision tim(nstep), temp(nstep), params(5)

double precision, dimension(:), allocatable :: nNf, T, eb, peb
double precision, dimension(:), allocatable :: Ea, pEa, inv_tauth, alpha
double precision, dimension(:,:), allocatable :: nN

double precision Ma, ddot, D0, Et, s, sigmaEt,res
double precision kb, magic_ratio, dEa, dt, pi, npEa
integer i,j,k, nEa

Ma = 3600*24*365*1.d6
! Extract parameters
ddot = params(1)*1000/Ma
!ddot = params(1)
D0 = params(2)
s = 10.d0**params(3)
Et = params(4)
sigmaEt = params(5)

! Define constants
kb = 8.617343d-5 
magic_ratio = ddot/D0
dEa = 0.01d0 

! Create variables for GAUSS model Lambert et al (Submitted)
nEa = 500
allocate (Ea(nEa), pEa(nEa))
Ea =  (/((5.d0*i)/nEa, i=1,nEa)/)
nEa = size(Ea)
pi = 4.d0*atan(1.d0)
pEa = exp(-0.5d0*((Ea-Et)/sigmaEt)**2)/(sigmaEt*sqrt(2.d0*pi))
npEa = sum(pEa)

allocate (T(nstep), nNf(nstep))
T = temp+273.15d0

! computes nN for the random Tt path
allocate (nN(nea,nstep),inv_tauth(nEa),alpha(nEa))
nN=0.d0
nNf=0.d0

do i = 2,nstep
dt = (tim(i-1)-tim(i))*Ma
inv_tauth = s*exp(-Ea/(kb*T(i-1)))
alpha = magic_ratio + inv_tauth
nN(:,i) = (nN(:,i-1)+magic_ratio*dt)*(1.d0/(1.d0+dt*alpha))
nNf(i) = dot_product(pEa,nN(:,i))
end do

nNf = nNf/npEa

ESRModel = nNf(nstep)

deallocate (Ea, pEa, T, nN, nNf, inv_tauth, alpha)

return
end function ESRModel
