      subroutine find_neighbours (icon,neighbour,npe,nelem,nnode)

! This subroutine computes the neighbour array which contains the "name" of
! each neighbouring element at the surface of the mesh (2D)

      implicit none

      integer icon(npe,nelem),neighbour(npe,nelem)

      integer,dimension(:),allocatable :: nc,kc,lc,ncs
      integer npe,nelem,nnode
      integer ie,k,ic,nctot,i,kp,icp,nstore,ke,km,jc,je

      allocate (nc(nnode),ncs(nnode))

! nc contains the number of elements that each node is connected to

      nc=0
        do ie=1,nelem
          do k=1,npe
          ic=icon(k,ie)
          nc(ic)=nc(ic)+1
          enddo
        enddo

! nctot is the total number of node-element connections
! ncs is the index array that stores where the list of elements connected
! to a given node starts in kc and lc

      nctot=sum(nc)
      allocate (kc(nctot),lc(nctot))
      ncs(1)=0
        do i=1,nnode-1
        ncs(i+1)=ncs(i)+nc(i)
        enddo

! kc and lc are the lists of elements to which each node is connected to and
! the number ("name") of the node in the corresponding element

      nc=0
        do ie=1,nelem
          do k=1,npe
          ic=icon(k,ie)
          nc(ic)=nc(ic)+1
          nstore=ncs(ic)+nc(ic)
          kc(nstore)=ie
          lc(nstore)=k
          enddo
        enddo

! neighbour is the neighbour list

      neighbour=0
        do ie=1,nelem
          do k=1,npe
          ic=icon(k,ie)
          kp=1+mod(k,npe)
          icp=icon(kp,ie)
            do je=1,nc(icp)
            nstore=ncs(icp)+je
            ke=kc(nstore)
              km=1+mod(lc(nstore),npe)
              jc=icon(km,ke)
                if (jc.eq.ic) then
                neighbour(k,ie)=ke
                endif
            enddo
          enddo
        enddo

      deallocate (nc,ncs,kc,lc)

      return
      end
