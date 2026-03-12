!---------------

subroutine calculate_fault_parameters (fault,nfault)

! given the fault geometry recalculates fault parameters needed by other routines

use Pecube

implicit none

type(faulttype) fault(nfault)
integer nfault,i,k
double precision xn,yn,xyn

  do i=1,nfault
    do k=1,fault(i)%n-1
    xn=fault(i)%x(k+1)-fault(i)%x(k)
    yn=fault(i)%y(k+1)-fault(i)%y(k)
    xyn=sqrt(xn**2+yn**2)
    xn=xn/xyn
    yn=yn/xyn
    fault(i)%xs(k)=xn
    fault(i)%ys(k)=yn
    enddo
  enddo

return

end subroutine calculate_fault_parameters
