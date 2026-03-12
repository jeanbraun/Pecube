      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n
      REAL a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      REAL bet
	real,dimension(:),allocatable::gam
	allocate (gam(n))
      if(b(1).eq.0.) stop 'in tridag'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
        print*,'tridag failed'
        stop
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      deallocate (gam)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.
