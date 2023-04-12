double precision function OSLModel (tim, temp, nstep, params)

! function to compute the relative trap content from a given time-temperature
! history (tim-temp) of length nstep, given a set of parameters provided by TL
! experiment (code loosely translated from Rabiul's matlab code)

implicit none

integer nstep, nrp,neb
double precision tim(nstep), temp(nstep), params(8)

double precision, dimension(:), allocatable :: nNf, rprime, pr, T, eb,peb
double precision, dimension(:), allocatable :: inv_tauath,inv_tauth
double precision, dimension(:,:,:), allocatable :: nN

double precision Ma, ddot, D0, a, Et,s, b, rhop, Eu
double precision kb, Hs, magic_ratio, dt, npr,npeb,inv_tauth1,alpha1
integer i,j,k

Ma = 3600*24*365*1.d6
! Extract parameters
ddot = params(1)*1000/Ma
D0 = params(2)
Et = params(3)
Eu = params(4)
s = 10.d0**params(5)
rhop = 10.d0**params(6)

! Define constants
kb = 8.617343d-5
Hs = 3.d15 !s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0
nrp = 100
neb = 50

! Define rprime range and tau athermic
allocate (rprime(nrp), pr(nrp), inv_tauath(nrp), inv_tauth(nrp),eb(neb))
allocate (peb(neb))
do i=1,nrp
  rprime(i)= 2.5d0*float(i-1)/(nrp-1)
  inv_tauath(i) =Hs*exp(-(rhop**-(1./3))*rprime(i))
enddo
pr = 3.d0*rprime**2.d0*exp(-rprime**3.d0)
npr = sum(pr)

do j=1,neb
eb(j)=Et*(j-1)/float(neb-1)
peb(j)=exp(-eb(j)/Eu)
enddo
npeb=sum(peb)

allocate (T(nstep), nNf(nstep))
T = temp+273.15d0

! computes nN for a Tt path
allocate (nN(neb,nrp,nstep))
nN=0.d0
nNf=0.d0

do i = 2,nstep
    dt = (tim(i-1)-tim(i))*Ma
    do k = 1,nrp
        do j = 1,neb
            inv_tauth1 = s*exp(-(Et-eb(j))/(kb*T(i-1)));
            alpha1=magic_ratio+inv_tauth1+inv_tauath(k);
            nN(j,k,i) =(nN(j,k,i-1)+magic_ratio*dt)*(1/(1+dt*alpha1));
            nNf(i)=nNf(i)+(peb(j)*(nN(j,k,i))*pr(k));
        enddo
    enddo
enddo

nNf = nNf/(npeb*npr)

OSLModel = nNf(nstep)

deallocate (rprime, pr, inv_tauath, inv_tauth, T, nN, peb,eb,nNf)

return
end
