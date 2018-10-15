      subroutine find_element (xdepthi,ydepthi,zdepth,tnow,x,y,z,t, &
                               xsurf,ysurf,zsurf,ielsurf,neighbour, &
                               iconsurf,icon,nelemsurf,nelem,nsurf,nz,nnode,npe,mpe, &
                               interpolate)

      implicit none

      double precision x(nnode),y(nnode),z(nnode),t(nnode)
      double precision xsurf(nsurf),ysurf(nsurf),zsurf(nsurf)
      integer iconsurf(npe,nelemsurf),icon(mpe,nelem)
      integer neighbour(npe,nelemsurf)
      logical interpolate,in
      double precision,dimension(:),allocatable :: xi,yi,zi,fi

      double precision xdepth,ydepth,zdepth,tnow
      double precision xdepthi,ydepthi
      double precision eps,side
      integer ielsurf,nelemsurf,nelem,nsurf,nz,nnode,npe,mpe
      integer i,k,i1,i2,iel3D,ic

      allocate (xi(mpe),yi(mpe),zi(mpe),fi(mpe))

      xdepth=max(xdepthi,minval(xsurf))
      xdepth=min(xdepthi,maxval(xsurf))
      ydepth=max(ydepthi,minval(ysurf))
      ydepth=min(ydepthi,maxval(ysurf))

      interpolate=.TRUE.

      eps=tiny(eps)
      if (ielsurf.eq.0) ielsurf=1

1     continue
        do k=1,npe
        i1=iconsurf(k,ielsurf)
        i2=iconsurf(1+mod(k,npe),ielsurf)
        side=(xdepth-xsurf(i1))*(ysurf(i2)-ysurf(i1))-(ydepth-ysurf(i1))*(xsurf(i2)-xsurf(i1))
          if (side.gt.eps) then
            if (neighbour(k,ielsurf).eq.0) then
            if (k.ne.npe) goto 2
            interpolate=.FALSE.
            return
            endif
          ielsurf=neighbour(k,ielsurf)
          goto 1
          endif
2       continue
        enddo

! we found the element

        do k=1,nz-1
        iel3D=(nz-1)*(ielsurf-1)+k
          do i=1,mpe
          ic=icon(i,iel3D)
          xi(i)=x(ic)
          yi(i)=y(ic)
          zi(i)=z(ic)
          fi(i)=t(ic)
          enddo
        call interpol3d (xdepth,ydepth,zdepth,tnow,xi,yi,zi,fi,mpe,in)
        if (in) return
        enddo

      interpolate=.FALSE.
      return

      end
