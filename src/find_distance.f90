!---------------------

subroutine find_distance (x,y,z,distance,fault,nfault)
    
! Beware this routine does not rotate the coordinate system
! it works in 'absolute' fault coordinates
    
use Pecube
    
implicit none
    
type(faulttype) fault(nfault)
double precision x,y,z,eps,distance,dist,xp,yp
integer nfault,i,k
    
! Note (x-x0)*xn+(y-y0)*yn is positive when x,y is on top of the plane
! passing through (x0,y0) and having (xn,yn) as a normal
    
eps=0.d0
    
distance=-1.d0
    
  do i=1,nfault 
  xp=x
  yp=z
! first check segment by segment
    do k=1,fault(i)%n-1
    if ((xp-fault(i)%x(k))*fault(i)%xs(k)+(yp-fault(i)%y(k))*fault(i)%ys(k).lt.eps) goto 111
    if ((xp-fault(i)%x(k))*fault(i)%ys(k)-(yp-fault(i)%y(k))*fault(i)%xs(k).lt.eps) goto 111
    if ((xp-fault(i)%x(k+1))*fault(i)%xs(k)+(yp-fault(i)%y(k+1))*fault(i)%ys(k).gt.eps) goto 111
    dist=(xp-fault(i)%x(k))*fault(i)%ys(k)-(yp-fault(i)%y(k))*fault(i)%xs(k)
    if (distance.lt.0.d0) distance=dist
    distance=min(distance,dist)
111 continue
    enddo
! then point by point
    do k=1,fault(i)%n
    if (k.ne.1) then
    if ((xp-fault(i)%x(k))*fault(i)%ys(k-1)-(yp-fault(i)%y(k))*fault(i)%xs(k-1).lt.eps) goto 222
    if ((xp-fault(i)%x(k))*fault(i)%xs(k-1)+(yp-fault(i)%y(k))*fault(i)%ys(k-1).lt.eps) goto 222
    endif
    if (k.ne.fault(i)%n) then
    if ((xp-fault(i)%x(k))*fault(i)%xs(k)+(yp-fault(i)%y(k))*fault(i)%ys(k).gt.eps) goto 222
    if ((xp-fault(i)%x(k))*fault(i)%ys(k)-(yp-fault(i)%y(k))*fault(i)%xs(k).lt.eps) goto 222
    endif
    if (k.eq.1) then
    dist=(xp-fault(i)%x(k))*fault(i)%ys(k)-(yp-fault(i)%y(k))*fault(i)%xs(k)
    elseif (k.eq.fault(i)%n) then
    dist=(xp-fault(i)%x(k))*fault(i)%ys(k-1)-(yp-fault(i)%y(k))*fault(i)%xs(k-1)
    else
    dist=(xp-fault(i)%x(k))**2+(y-fault(i)%y(k))**2
    if (dist.gt.0.d0) dist=sqrt(dist)
    endif
    if (distance.lt.0.d0) distance=dist
    distance=min(distance,dist)
222 continue
    enddo
  enddo

return

end subroutine find_distance
