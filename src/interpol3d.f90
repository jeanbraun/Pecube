      subroutine interpol3d (x,y,z,f,xi,yi,zi,fi,ni,in)

      implicit none

      double precision xi(ni),yi(ni),zi(ni),fi(ni)
      double precision,dimension(:),allocatable :: ri,si,ti
      logical in
      integer ni,i,niter
      double precision x,y,z,f
      double precision rp,sp,tp,sum1,sum2,sum,r,s,t,drst

      allocate (ri(ni),si(ni),ti(ni))

      if (ni.eq.8) then
      ri=(/-1.,1.,1.,-1.,-1.,1.,1.,-1./)
      si=(/-1.,-1.,1.,1.,-1.,-1.,1.,1./)
      ti=(/-1.,-1.,-1.,-1.,1.,1.,1.,1./)
      r=0.
      s=0.
      t=0.

      niter=0
    1 continue
      rp=r
      sp=s
      tp=t
      sum1=0.
      sum2=0.
        do i=1,ni
        sum=(1.+s*si(i))*(1.+t*ti(i))*xi(i)
        sum1=sum1+sum
        sum2=sum2+sum*ri(i)
        enddo
      r=(8.*x-sum1)/sum2
      sum1=0.
      sum2=0.
        do i=1,ni
        sum=(1.+t*ti(i))*(1.+r*ri(i))*yi(i)
        sum1=sum1+sum
        sum2=sum2+sum*si(i)
        enddo
      s=(8.*y-sum1)/sum2
      sum1=0.
      sum2=0.
        do i=1,ni
        sum=(1.+r*ri(i))*(1.+s*si(i))*zi(i)
        sum1=sum1+sum
        sum2=sum2+sum*ti(i)
        enddo
      t=(8.*z-sum1)/sum2
      drst=(r-rp)**2+(s-sp)**2+(t-tp)**2
      niter=niter+1
      if (drst.gt.1.d-10.and.niter.lt.10) goto 1

      f=0.
        do i=1,ni
        f=f+fi(i)*(1.+ri(i)*r)*(1.+si(i)*s)*(1.+ti(i)*t)/8.
        enddo

      else

      r=((x-xi(1))*(yi(3)-yi(1))-(y-yi(1))*(xi(3)-xi(1)))/ &
        ((xi(2)-xi(1))*(yi(3)-yi(1))-(yi(2)-yi(1))*(xi(3)-xi(1)))
      s=((x-xi(1))*(yi(2)-yi(1))-(y-yi(1))*(xi(2)-xi(1)))/ &
        ((xi(3)-xi(1))*(yi(2)-yi(1))-(yi(3)-yi(1))*(xi(2)-xi(1)))
      t=(-2.*z+(1.-r-s)*(zi(1)+zi(4))+r*(zi(2)+zi(5))+s*(zi(3)+zi(6)))/ &
        ((1.-r-s)*(zi(1)-zi(4))+r*(zi(2)-zi(5))+s*(zi(3)-zi(6)))
      f=(1.-r-s)*(1.-t)/2.*fi(1)+r*(1.-t)/2.*fi(2)+s*(1.-t)/2.*fi(3) &
       +(1.-r-s)*(1.+t)/2.*fi(4)+r*(1.+t)/2.*fi(5)+s*(1.+t)/2.*fi(6)

      endif

      in=.false.
      if ((t+1.)*(t-1.).le.0.) in=.true.

      deallocate (ri,si,ti)

      return
      end
