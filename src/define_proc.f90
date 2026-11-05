      subroutine define_proc (nproc,nsurf,nz,proc)

      implicit none

      integer nproc,nsurf,nz,iprc,node1,node2,k,i
      integer proc(nsurf*nz)

        do iprc=0,nproc-1
        node1=(iprc*nsurf)/nproc+1
        node2=((iprc+1)*nsurf)/nproc
        if (iprc.eq.nproc) node2=nsurf
          do i=node1,node2
            do k=1,nz
            proc((i-1)*nz+k)=iprc
            enddo
          enddo
        enddo

      return
      end
