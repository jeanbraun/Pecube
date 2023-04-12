      subroutine calculate_ESR (nrec,nsurf,nz,istep,nN,ESR_doser,ESR_D0,ESR_s,ESR_Et,ESR_sigmaEt)
!__________________________________________________________________________
! modified version of calculate_ages.f90 for ESR
!__________________________________________________________________________

      implicit none

      integer nsurf,nrec,nz,istep,i,ij,irec,k
      real*4 nN(nsurf)
      real*4,dimension(:),allocatable::depth,ztime,ztemp,ztime_3,ztemp_3    
      real*4 zjunk,timenow
      double precision ESR_doser,ESR_D0,ESR_s,ESR_Et,ESR_sigmaEt,params(6),res
      integer aftmodel
      integer iproc,nd

      external ESRModel
      double precision ESRModel

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

        params(1) = ESR_doser
        params(2) = ESR_D0
        params(3) = ESR_s
        params(4) = ESR_Et
        params(5) = ESR_sigmaEt
        res = ESRModel (dble(ztime), dble(ztemp), nrec, params)

        nN(i) = res

        enddo

      if (iproc.eq.0.and.nd.eq.0) print*,''

      deallocate (ztime,ztime_3,ztemp,ztemp_3,depth)

      return
      end
