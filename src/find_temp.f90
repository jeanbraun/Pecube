      subroutine find_temp (t,z,n,depth,tnow)

      implicit none

! interpolation routine

      double precision t(n),z(n)
      double precision depth,tnow,x0,x1,r
      integer n,i

      tnow=t(n)
      if (depth.gt.z(n)) return

      tnow=t(1)
      if (depth.lt.z(1)) return

        do i=1,n-1
        x0=z(i)
        x1=z(i+1)
          if ((depth-x0)*(depth-x1).le.0.) then
          r=(depth-x0)/(x1-x0)
          tnow=t(i)+(t(i+1)-t(i))*r
          return
          endif
        enddo

      return
      end

