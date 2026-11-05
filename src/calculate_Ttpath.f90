      subroutine calculate_ttpath (nrec,nsurf,nz,istep,nhist,thist,temphist, &
                                   iproc,nd)
!--------------------------------------------------------------------------
! subroutine to return the temperature for a given set of time values
!--------------------------------------------------------------------------

      implicit none

      integer nsurf,nrec,nz,istep,i,ij,irec,k
      integer nhist(nsurf)
      double precision thist(1024,nsurf),temphist(1024,nsurf)
      real*4,dimension(:),allocatable::depth,ztime,ztemp 
      real*4 zjunk,timenow
      integer iproc,nd
      double precision rat
      double precision,dimension(:,:),allocatable::mem

      allocate (ztime(nrec),ztemp(nrec),depth(nrec))

        do i=1,nsurf
        ij=(i-1)*nz+1
            do irec=1,nrec
            read (100+istep,rec=irec) ztime(irec),(zjunk,k=1,i-1),ztemp(irec), &
                                    (zjunk,k=i+1,nsurf), &
                                    (zjunk,k=1,i-1),depth(irec)
            enddo
        timenow=ztime(nrec)
        ztime=ztime-timenow

          do k=1,nhist(i)
          temphist(k,i)=ztemp(nrec)
            if (thist(k,i).gt.ztime(1))then
            temphist(k,i)=ztemp(1)
            temphist(k,i)=-9999.d0
            else
              do irec=1,nrec-1
                if ((thist(k,i)-ztime(irec))*(thist(k,i)-ztime(irec+1)).le.0.d0) then
                rat=(thist(k,i)-ztime(irec+1))/(ztime(irec)-ztime(irec+1))
                temphist(k,i)=ztemp(irec+1)+(ztemp(irec)-ztemp(irec+1))*rat
                endif
              enddo
            endif
          enddo

        enddo

      if (iproc.eq.0.and.nd.eq.0) print*,''

      deallocate (ztime,ztemp,depth)

      return
      end
