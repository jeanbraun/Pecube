      subroutine calculate_OSL (nrec,nsurf,nz,istep,nN,OSL_doser,OSL_D0,OSL_Et,OSL_Eu,OSL_logs,OSL_logrho)
!__________________________________________________________________________
! modified version of calculate_ages.f90 for OSL
!__________________________________________________________________________

      implicit none

      integer nsurf,nrec,nz,istep,i,ij,irec,k
      real*4 nN(nsurf)
      real*4,dimension(:),allocatable::depth,ztime,ztemp,ztime_3,ztemp_3    
      real*4 zjunk,timenow,ftld(17),ftldmean,ftldsd,grainsize(nsurf)
      double precision fdist(200),alo,oldest_age,final_age,fmean
      double precision ftdist(20,nsurf),tprevious
      double precision OSL_doser,OSL_D0,OSL_Et,OSL_Eu,OSL_logs,OSL_logrho,params(6),res
      integer aftmodel
      integer iproc,nd

      external OSLModel
      double precision OSLModel

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

        params(1) = OSL_doser
        params(2) = OSL_D0
        params(3) = OSL_Et
        params(4) = OSL_Eu
        params(5) = OSL_logs
        params(6) = OSL_logrho
        res = OSLModel (dble(ztime), dble(ztemp), nrec, params)

        nN(i) = res

        enddo

      if (iproc.eq.0.and.nd.eq.0) print*,''

      deallocate (ztime,ztime_3,ztemp,ztemp_3,depth)

      return
      end
