      subroutine find_dt (zl,diffusivity,nsurf,zsurf,zsurfp, &
                          nz,timesurf,timesurfp,istep,eps,ilog,nd, &
                          dt,ntime,istatic,fault,nfault,usurf)

      use Pecube

      implicit none

!      implicit real*8 (a-h,o-z)

      type (faulttype) fault(nfault)

      double precision zsurf(nsurf),zsurfp(nsurf),usurf(nsurf)
      double precision zl,diffusivity,timesurf,timesurfp,eps,dt
      double precision dt1,dt2,dt3,Pecletz,Pesurf,dzmax
      integer nsurf,nz,istep,ilog,ntime,istatic,nfault,i,k,nd

! first constrain on time step from conduction

      dt1=zl**2/abs(diffusivity)/100.

! second constrain on time step from advection

      dt2=dt1
      Pecletz=0.d0
        do i=1,nfault
          do k=1,fault(i)%nstep
          if ((timesurf-fault(i)%timestart(k))*(timesurf-fault(i)%timeend(k)).le.0.d0) &
                  Pecletz=max(Pecletz,abs(fault(i)%velo(k)))
          if ((timesurfp-fault(i)%timestart(k))*(timesurfp-fault(i)%timeend(k)).le.0.d0) &
                  Pecletz=max(Pecletz,abs(fault(i)%velo(k)))
          if ((fault(i)%timestart(k)-timesurf)*(fault(i)%timestart(k)-timesurfp).le.0.d0) &
                  Pecletz=max(Pecletz,abs(fault(i)%velo(k)))
          enddo
        enddo
      Pecletz=max(Pecletz,maxval(abs(usurf)))
      if (abs(Pecletz).gt.eps) dt2=zl/abs(Pecletz)/100.

! third constrain on time step from surface lowering

      dt3=dt1
        if (istep.ne.0) then
        dzmax=0.
          do i=1,nsurf
          dzmax=max(dzmax,zsurfp(i)-zsurf(i))
          enddo
        Pesurf=dzmax/(timesurf-timesurfp)
        if (abs(Pesurf).gt.eps) dt3=zl/Pesurf/5.
        endif

! find optimum time step and number of steps

      dt=min(dt1,dt2)
      dt=min(dt,dt3)
      if (istep.ne.0) dt=min(dt,timesurf-timesurfp)

      ntime=int((timesurf-timesurfp)/dt)
      ntime=ntime+1
      dt=(timesurf-timesurfp)/ntime
      istatic=0

        if (istep.eq.0) then
        ntime=1
        dt=0.
        istatic=1
        endif

      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'ntime/dt/dt1/dt2/dt3= ',ntime,dt,dt1,dt2,dt3

      return
      end
