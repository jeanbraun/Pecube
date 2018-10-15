double precision function TLModel (tim, temp, nstep, params)

! function to compute the relative trap content from a given time-temperature
! history (tim-temp) of length nstep, given a set of parameters provided by TL
! experiment (code loosely translated from Rabiul's matlab code)

implicit none

integer nstep, nrp
double precision tim(nstep), temp(nstep), params(8)

double precision, dimension(:), allocatable :: nNf, rprime, pr, inv_tauath, inv_tauth, xkd, xk, T
double precision, dimension(:,:), allocatable :: nN

double precision Ma, ddot, D0, a, Et,s, b, rhop
double precision kb, Hs, magic_ratio, dt, npr
integer i

Ma = 3600*24*365*1.d6
! Extract parameters
ddot = params(1)/(1.d3*365.d0*24.d0*3600.d0)!                  % Gy/ka to Gy/seconds
D0 = params(2)
a  = params(3)
Et = params(4)
s = 10.d0**params(6)
b = params(7)
rhop = 10.d0**params(8)

! Define constants
kb = 8.617343d-5
Hs = 3.d15 !s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0
nrp = 100
dt = ((maxval(tim)-minval(tim))/(nstep-1))*Ma

! Define rprime range and tau athermic
allocate (rprime(nrp), pr(nrp), inv_tauath(nrp), inv_tauth(nrp))
  do i=1,nrp
  rprime(i)= 2.5d0*float(i-1)/(nrp-1)
  enddo
pr = 3.d0*rprime**2.d0*exp(-rprime**3.d0)
npr = sum(pr)

inv_tauath = (Hs*exp(-(rhop**(-1.d0/3.d0))*rprime))

allocate (T(nstep), nNf(nstep))
T = temp+273.15d0

! computes nN for a Tt path
allocate (nN(nrp,nstep), xk(nrp), xkd(nrp))
nN=0.d0
nNf=0.d0

  do i = 2,nstep
  inv_tauth = s*exp(-Et/(kb*T(i-1) ))
  xkd = -a*magic_ratio*(1.d0-nN(:,i-1))**(a-1)-b*(nN(:,i-1))**(b-1)*inv_tauth-inv_tauath
  xk = magic_ratio*(1.d0-nN(:,i-1))**a-(nN(:,i-1))**b*inv_tauth-inv_tauath*nN(:,i-1)
  nN(:,i) = nN(:,i-1)+xk*dt/(1.-dt*xkd)
  nNf(i) = sum(nN(:,i)*pr)
  end do
nNf = nNf/npr

TLModel = nNf(nstep)

deallocate (rprime, pr, inv_tauath, inv_tauth, T, nN, xk, xkd, nNf)

return

end function TLModel

