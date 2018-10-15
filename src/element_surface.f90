      double precision function  element_surface (npe,icon,xg,yg,nnode)

      use Pecube

      implicit none

      double precision xg(nnode),yg(nnode)
      double precision jcb(2,2)
      integer icon(npe)
      double precision,dimension(:),allocatable :: x,y,dhdr,dhds

      integer npe,k,nnode
      double precision eps,r,s,w

      allocate (x(npe),y(npe))
      allocate (dhdr(npe),dhds(npe))

      eps=tiny(eps)

        do k=1,npe
        x(k)=xg(icon(k))
        y(k)=yg(icon(k))
        enddo

        r=0.d0
        s=0.d0
        w=4.d0

          if (npe.eq.4) then
          dhdr(1)=-(1.d0-s)/4.d0
          dhdr(2)=(1.d0-s)/4.d0
          dhdr(3)=(1.d0+s)/4.d0
          dhdr(4)=-(1.d0+s)/4.d0
          dhds(1)=-(1.d0-r)/4.d0
          dhds(2)=-(1.d0+r)/4.d0
          dhds(3)=(1.d0+r)/4.d0
          dhds(4)=(1.d0-r)/4.d0
          else
          dhdr(1)=-1.d0/2.d0
          dhdr(2)=1.d0/2.d0
          dhdr(3)=0.d0
          dhds(1)=-1.d0/2.d0
          dhds(2)=0.d0
          dhds(3)=1.d0/2.d0
          endif

        jcb=0.d0
          do k=1,npe
          jcb(1,1)=jcb(1,1)+dhdr(k)*x(k)
          jcb(1,2)=jcb(1,2)+dhdr(k)*y(k)
          jcb(2,1)=jcb(2,1)+dhds(k)*x(k)
          jcb(2,2)=jcb(2,2)+dhds(k)*y(k)
          enddo

        element_surface=(jcb(1,1)*jcb(2,2)-jcb(1,2)*jcb(2,1))*w

          if (element_surface.le.eps) then
          print*,'volume= ',element_surface
          stop 'Element bow-tied or collapsed'
          endif

      return
      end
