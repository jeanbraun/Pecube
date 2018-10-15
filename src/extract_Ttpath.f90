      subroutine extract_ttpath (nrec,nsurf,nz,istep,thist,temphist, &
                                 iproc,nd)
!__________________________________________________________________________
!This file has been changed to use ketcham program. 
!original version is found in calculate_ages_old.f90 in this same folder
!by VKP
!__________________________________________________________________________

      implicit none

      integer nsurf,nrec,nz,istep,i,ij,irec,k
      double precision thist(nrec,nsurf),temphist(nrec,nsurf)
      real*4 zjunk,timenow,t1,t2
      integer iproc,nd
      double precision rat,depth

        do i=1,nsurf
        ij=(i-1)*nz+1
            do irec=1,nrec
            read (100+istep,rec=irec) t1,(zjunk,k=1,i-1),t2, &
                                    (zjunk,k=i+1,nsurf), &
                                    (zjunk,k=1,i-1),zjunk
            thist(irec,i)=t1
            temphist(irec,i)=t2
            enddo
        timenow=thist(nrec,i)
        thist=thist-timenow
        enddo

      if (iproc.eq.0.and.nd.eq.0) print*,''

      return
      end
