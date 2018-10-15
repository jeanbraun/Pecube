      subroutine calculate_TL (nrec,nsurf,nz,istep,nN,TL_doser,TL_D0,TL_a,TL_b,TL_Et,TL_logs,TL_logrho)
!__________________________________________________________________________
! modified version of calculate_ages.f90 for TL
!__________________________________________________________________________

      implicit none

      integer nsurf,nrec,nz,istep,i,ij,irec,k
      real*4 nN(nsurf)
      real*4,dimension(:),allocatable::depth,ztime,ztemp,ztime_3,ztemp_3    
      real*4 zjunk,timenow,ftld(17),ftldmean,ftldsd,grainsize(nsurf)
      double precision fdist(200),alo,oldest_age,final_age,fmean
      double precision ftdist(20,nsurf),tprevious
      double precision TL_doser,TL_D0,TL_a,TL_b,TL_Et,TL_logs,TL_logrho,params(8),res
      integer aftmodel
      integer iproc,nd

      external TLModel
      double precision TLModel

      allocate (ztime(nrec),ztemp(nrec),depth(nrec),ztime_3(nrec),ztemp_3(nrec))

      nN=0.

        do i=1,nsurf
        if (iproc.eq.0.and.nd.eq.0) call screen_counter (i,nsurf,0)
        ij=(i-1)*nz+1
            do irec=1,nrec
            read (100+istep,rec=irec) ztime(irec),(zjunk,k=1,i-1),ztemp(irec), &
                                    (zjunk,k=i+1,nsurf), &
                                    (zjunk,k=1,i-1),depth(irec)
            enddo
        timenow=ztime(nrec)
        ztime=ztime-timenow

        params(1) = TL_doser
        params(2) = TL_D0
        params(3) = TL_a
        params(4) = TL_Et
        params(5) = 0.
        params(6) = TL_logs
        params(7) = TL_b
        params(8) = TL_logrho
        res = TLModel (ztime, ztemp, nrec, params)

        nN(i) = res

        enddo

      if (iproc.eq.0.and.nd.eq.0) print*,''

      deallocate (ztime,ztime_3,ztemp,ztemp_3,depth)

      return
      end
