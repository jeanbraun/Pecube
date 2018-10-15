      subroutine find_velo (x,y,z,time,dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,forward)

! this extra layer between Pecube and the velocity function geometry ensures
! that the advection of the "lagrangian particles" for which we track the thermal
! history (to calculate the age) is second order accurate; it uses a mid-point
! algorithm.

! Note that this utility has been turned off in this new version of Pecube
! this routine is useless but has been kept for compatibility with ancient versions

      use Pecube

      implicit none

      double precision x,y,z,time,dt,vx,vy,vz,xmin,xmax,ymin,ymax
      integer nfault
      logical forward

      type (faulttype) fault(nfault)

      if (dt.gt.0.d0) then
      call find_velocity (x,y,z,vx,vy,vz,fault,nfault,time-dt*1.d-8,1,xmin,xmax,ymin,ymax)
!      call find_velocity (x+vx*dt/2.d0,y+vy*dt/2.d0,z+vz*dt/2.d0,vx,vy,vz,fault,nfault,time+dt/2.,1,xmin,xmax,ymin,ymax)
      else
      call find_velocity (x,y,z,vx,vy,vz,fault,nfault,time-dt*1.d-8,1,xmin,xmax,ymin,ymax)
!      call find_velocity (x+vx*dt/2.d0,y+vy*dt/2.d0,z+vz*dt/2.d0,vx,vy,vz,fault,nfault,time-dt/2.,1,xmin,xmax,ymin,ymax)
      endif

      return

      end
