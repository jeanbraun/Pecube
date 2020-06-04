!----------------------------

subroutine read_in_fault_parameters (fault,nfault,xlon1,xlat1,xlon2,xlat2,xl,yl,zl,timeend, &
                                     nd,range,param,run)

! this subroutines reads the fault information from the fault_parameters.txt file
! and generates the necessary complementary information to be used for fault
! kinematics and movement

use Pecube

implicit none

integer nfault,i,j,k,icol,jcol,kcol,jd
type(faulttype) fault(nfault)
type (parameters) p
character*5 run

double precision xn,yn,xyn,x1,y1,x2,y2,timeend
double precision xlon1,xlat1,xlon2,xlat2,xl,yl,zl
character line*1024
logical logvelo

real*4 range(2,*),param(*)
integer nd

nd = 0
call read_input_file (run//'/input/Pecube.in', 0, p, nd, range, param)

logvelo = .false.
if (p%logarithmic_velocity.ne.0) logvelo = .true.

if (nfault.eq.0) return

x1 = p%x1
y1 = p%y1
x2 = p%x2
y2 = p%y2

x1=(x1-xlon1)/(xlon2-xlon1)*xl
y1=(y1-xlat1)/(xlat2-xlat1)*yl
x2=(x2-xlon1)/(xlon2-xlon1)*xl
y2=(y2-xlat1)/(xlat2-xlat1)*yl

  do i=1,nfault

  fault(i)%n = p%npoint(i)

  if (fault(i)%n.lt.0) then
  allocate (fault(i)%x(4),fault(i)%y(4))
  fault(i)%x(1) = p%bottom_left
  fault(i)%x(2) = p%bottom_right
  fault(i)%x(4) = p%upper_right
  fault(i)%x(3) = p%upper_left
  else
  allocate (fault(i)%x(fault(i)%n),fault(i)%y(fault(i)%n))
    do k=1,fault(i)%n
    fault(i)%x(k) = p%r(k,i)
    fault(i)%y(k) = p%s(k,i)
    fault(i)%y(k)=fault(i)%y(k)+zl
    enddo
  fault(i)%x1=x1;fault(i)%y1=y1;fault(i)%x2=x2;fault(i)%y2=y2
  endif

  fault(i)%nstep = p%nstep(i)
  fault(i)%static = p%static(i)

  allocate (fault(i)%timestart(fault(i)%nstep),fault(i)%timeend(fault(i)%nstep), &
            fault(i)%velo(fault(i)%nstep))
    do k=1,fault(i)%nstep
    fault(i)%timestart(k) = p%time_start(k,i)
    fault(i)%timeend(k) = p%time_end(k,i)
    fault(i)%velo(k) = p%velo(k,i)
    enddo

    do k=1,fault(i)%nstep
    fault(i)%timestart(k)=timeend-fault(i)%timestart(k)
    fault(i)%timeend(k)=timeend-fault(i)%timeend(k)
    enddo

  enddo

  if (logvelo) then
    do i=1,nfault
      do k=1,fault(i)%nstep
      fault(i)%velo(k)=10.d0**fault(i)%velo(k)
      enddo
    enddo
  endif
 
  do i=1,nfault
  if (fault(i)%n.gt.0) &
  allocate (fault(i)%xs(fault(i)%n-1),fault(i)%ys(fault(i)%n-1))
  enddo

call calculate_fault_parameters (fault,nfault)
  
  do i=1,nfault
  if (fault(i)%n.gt.0) then
  if (fault(i)%x(fault(i)%n).gt.fault(i)%x(1) .and. &
      fault(i)%y(fault(i)%n).lt.fault(i)%y(1)) &
          fault(i)%velo(1:fault(i)%nstep)=-fault(i)%velo(1:fault(i)%nstep)
  if (fault(i)%x(fault(i)%n).lt.fault(i)%x(1) .and. &
      fault(i)%y(fault(i)%n).gt.fault(i)%y(1)) &
          fault(i)%velo(1:fault(i)%nstep)=-fault(i)%velo(1:fault(i)%nstep)
  endif
  enddo
  
  do i=1,nfault
  xn=fault(i)%y2-fault(i)%y1
  yn=fault(i)%x1-fault(i)%x2
  xyn=sqrt(xn**2+yn**2)
  fault(i)%xn=xn/xyn
  fault(i)%yn=yn/xyn
  enddo

return

end subroutine read_in_fault_parameters
