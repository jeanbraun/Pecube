      subroutine calculate_OSL (nrec,nsurf,nz,istep,nN,ageOSL,OSL_doser,&
                    OSL_D0,OSL_Et,OSL_Eu,&
                    OSL_logs,OSL_logrho,OSL_Lmax,OSL_a,OSL_Model)
!__________________________________________________________________________
! modified version of calculate_ages.f90 for OSL
! Calculate OSL n/N or ages signals (test version).
! Trapped-charge models are all stored in Trapped-charge.f90 file.
! 
! Author: Maxime Bernard
!_________________________________________________________________________

      implicit none

      integer nsurf,nrec,nz,istep,i,ij,irec,k
      real*4 nN(nsurf), ageOSL(nsurf)
      real*4,dimension(:),allocatable::depth,ztime,ztemp,ztime_3,ztemp_3    
      real*4 zjunk,timenow
      double precision OSL_doser,OSL_D0,OSL_Et,OSL_Eu,OSL_logs,OSL_logrho,params(8),res
      double precision OSL_a,OSL_Lmax, OSL_b
      integer iproc,nd, OSL_Model

      external Luminescence
      double precision nnf_array(nrec,1), ageOSL_array(nrec,1)
      integer no, io

      allocate (ztime(nrec),ztemp(nrec),depth(nrec),ztime_3(nrec),ztemp_3(nrec))

      nN=0.
      ageOSL = 0.

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

        no = 1
        io = 1
        nnf_array(:,1) = 0.d0
        ageOSL_array(:,1) = 0.d0
        params(1) = OSL_doser
        params(2) = OSL_D0
        params(3) = OSL_Et
        params(4) = OSL_Eu ! Eu or SigmaEt according to the model used
        params(5) = OSL_logs
        params(6) = OSL_logrho
        params(7) = OSL_a
        params(8) = OSL_Lmax
        if (OSL_Model.eq.0) then
            call SSE_BTS_FAD_model (dble(ztime), dble(ztemp), nrec,ageOSL_array, OSL_doser, OSL_D0, OSL_Et, OSL_logs,&
                OSL_logrho, OSL_Eu,nnf_array,OSL_Lmax, no, io)
                
        elseif (OSL_Model.eq.1) then
            call GOK_Gauss_FAD_model (dble(ztime), dble(ztemp), nrec,ageOSL_array, OSL_doser, OSL_D0, OSL_logs, OSL_Et, OSL_Eu, &
                                        OSL_logrho,OSL_a,OSL_Lmax,nnf_array, no, io)

        elseif (OSL_Model.eq.2) then ! GOK Fad model
            call GOK_FAD_model (dble(ztime), dble(ztemp), nrec, OSL_doser, OSL_D0, OSL_a, OSL_Et, OSL_logs, OSL_Eu, OSL_logrho,&
                         OSL_Lmax,nnf_array,ageOSL_array, no, io)

        elseif (OSL_Model.eq.3) then ! Gauss FAD model (to test )
            call SSE_Gauss_FAD_model (dble(ztime), dble(ztemp), nrec, ageOSL_array, OSL_doser, OSL_D0, OSL_logs, OSL_Et, OSL_Eu, &
                                        OSL_logrho,OSL_Lmax,nnf_array,no, io)
        endif

        nN(i) = nnf_array(nrec,1)
        ageOSL(i) = ageOSL_array(nrec,1)
        ! print *, 'nNf = ', nnf_array(nrec,1)

        enddo

      if (iproc.eq.0.and.nd.eq.0) print*,''

      deallocate (ztime,ztime_3,ztemp,ztemp_3,depth)

      return
      end
