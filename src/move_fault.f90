!----------------

subroutine move_fault (fault,faultp,nfault,dt,time,xmin,xmax,ymin,ymax)

use Pecube

implicit none

type (faulttype) fault(nfault),faultp(nfault)
double precision dt,distance,eps,vx,vy,vz,time
double precision xmin,xmax,ymin,ymax
integer nfault,i,j,k
logical move

faultp=fault

eps=1.d-10

do i=1,nfault

  do j=1,nfault

  if (fault(i)%n.lt.0) then
  move=.false.
  else
  move=.false.
    do k=1,fault(i)%n
    call find_distance (faultp(i)%x(k),0.d0,faultp(i)%y(k),distance,faultp(j),1)
    if (distance.gt.eps) move=.true.
    enddo
  endif

  if (move) then
    do k=1,faultp(i)%n
    call find_velocity (faultp(i)%x(k),0.d0,faultp(i)%y(k),vx,vy,vz,faultp(j),1,time,0,xmin,xmax,ymin,ymax)
    call find_velocity (faultp(i)%x(k)+vx*dt/2.d0,0.d0,faultp(i)%y(k)+vz*dt/2.d0,vx,vy,vz,faultp(j),1,time+dt/2.d0, &
                        0,xmin,xmax,ymin,ymax)
    fault(i)%x(k)=fault(i)%x(k)+dt*vx
    fault(i)%y(k)=fault(i)%y(k)+dt*vz
    enddo
  endif

  enddo

enddo

call calculate_fault_parameters (fault,nfault)

return

end subroutine move_fault
