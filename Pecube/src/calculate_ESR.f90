      subroutine calculate_ESR (nrec,nsurf,nz,istep,nN,ageESR,&
                ESR_doser,ESR_D0,ESR_s,ESR_Et,ESR_sigmaEt,ESR_Lmax,&
                ESR_a, ESR_b, ESR_Model)
!__________________________________________________________________________
! modified version of calculate_ages.f90 for ESR.
! Calculate OSL n/N or ages signals (test version).
! Trapped-charge models are all stored in Trapped-charge.f90 file.
! 
! Author: Maxime Bernard
!__________________________________________________________________________

      implicit none

      integer nsurf,nrec,nz,istep,i,ij,irec,k
      real*4 nN(nsurf),ageESR(nsurf)
      real*4,dimension(:),allocatable::depth,ztime,ztemp,ztime_3,ztemp_3    
      real*4 zjunk,timenow
      double precision ESR_doser,ESR_D0,ESR_s,ESR_Et,ESR_sigmaEt,params(6),res
      double precision ESR_a, ESR_b, ESR_logrho, ESR_Lmax
      integer aftmodel
      integer iproc,nd

      external Luminescence
      double precision nnf_array(nrec,1), ageESR_array(nrec,1)
      integer no, io, ESR_Model, ESRage, Agemax

      allocate (ztime(nrec),ztemp(nrec),depth(nrec),ztime_3(nrec),ztemp_3(nrec))

      nN=0.
      ageESR = 0.

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
        ageESR_array(:,1) = 0.d0
        params(1) = ESR_doser
        params(2) = ESR_D0
        params(3) = ESR_s
        params(4) = ESR_Et
        params(5) = ESR_sigmaEt
        ! print *, 'ESR_model = ', ESR_Model
        if (ESR_Model.eq.0) then
            call SSE_Gauss_model (dble(ztime), dble(ztemp), nrec, ESR_doser, ESR_D0, ESR_s, ESR_Et,&
                                ESR_sigmaEt,ESR_Lmax,nnf_array,ageESR_array, no, io)
            ! res = ESRModel (dble(ztime), dble(ztemp), nrec, params)
        elseif (ESR_Model.eq.1) then
            call GOK_model (dble(ztime), dble(ztemp), nrec, ESR_doser, ESR_D0, ESR_a, ESR_Et,&
                            ESR_s, ESR_b,ESR_Lmax,nnf_array,ageESR_array, no, io)
        endif
        nN(i) = nnf_array(nrec,1)
        ageESR(i) = ageESR_array(nrec,1)

        enddo
        ! print *, 'Age ESR = ', ageESR

      if (iproc.eq.0.and.nd.eq.0) print*,''

      deallocate (ztime,ztime_3,ztemp,ztemp_3,depth)

      return
      end
