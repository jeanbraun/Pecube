      subroutine make_matrix (mpe,ael,bel,icon,xg,yg,zg,ug,xgp,ygp,zgp,ugp, & !VTK
                              kfix,diffusivity,heatproduction,&
                              x1f,y1f,x2f,y2f,def,dif, &
                              alpha0,dt0,time,tempp,nnode,istatic, &
                              fault,nfault,rebound, &
                              xming,xmaxg,yming,ymaxg,zming,zmaxg,friction, &
                              lsf,lsfx,lsfy,nx,ny,nz,tfinal)

      use Pecube

      implicit none

! this is where the FE matrices and RHS vectors are built

      type (faulttype) fault(nfault)

      double precision xg(nnode),yg(nnode),zg(nnode),ug(nnode) !VKP
      double precision xgp(nnode),ygp(nnode),zgp(nnode),ugp(nnode) !VKP
      double precision tempp(nnode),rebound(nnode)
      double precision ael(mpe,mpe),bel(mpe)
      double precision jcb(3,3),jcbi(3,3),jcbp(3,3),jcbip(3,3)
      integer icon(mpe),kfix(nnode),nx,ny,nz
      double precision,dimension(:),allocatable :: rr,ss,tt,ww
      double precision,dimension(:,:),allocatable :: aelp,b,bp,vv
      double precision,dimension(:),allocatable :: x,y,z,u,xp,yp,zp,up,h,ht,dhdr,dhds,dhdt,rb,temp !VKP

      integer mpe,nint,k,i,iint,j,nnode,istatic,nfault
      double precision diffusivity,heatproduction,x1f,y1f,x2f,y2f,def,dif,alpha0,dt0,time
      double precision eps,rhoc,cond,heat,dynamic,dt,alpha,timep,xint,yint,zint,rint,uint,tint !VKP
      double precision veloz,volume,volumep,velox,veloy,uvnorm,xmin,xmax,ymin,ymax,velou !VKP
      double precision xming,xmaxg,yming,ymaxg,zming,zmaxg,friction,specific_heat
      double precision zmin,zmax,dx,dy,dz,tau,xcond,xmass,xcondp,xmassp,fixt
      double precision r,s,t,w,vx,vy,vz,epsij(3,3),e2d
      double precision lsf(nx*ny,nz),lsfx(nx*ny,nz),lsfy(nx*ny,nz),xinit,yinit,zinit,tfinal
      double precision heat_production,thermal_diffusivity,hinterpolate
      external heat_production,thermal_diffusivity,hinterpolate

      allocate (aelp(mpe,mpe),b(mpe,3),bp(mpe,3),rb(mpe))
      allocate (x(mpe),y(mpe),z(mpe),u(mpe),xp(mpe),yp(mpe),zp(mpe),up(mpe),temp(mpe)) !VKP
      allocate (vv(mpe,3))
      allocate (h(mpe),ht(mpe),dhdr(mpe),dhds(mpe),dhdt(mpe))

      eps=tiny(eps)

        if (mpe.eq.8) then
        nint=8
        else
        nint=6
        endif

      rhoc=1.
      specific_heat=1.d3
      cond=diffusivity
      heat=heatproduction
      dynamic=1.
      dt=dt0
      alpha=alpha0
        if (istatic.eq.1) then
        dynamic=0.
        dt=1.
        alpha=1.
        endif

      allocate (rr(nint),ss(nint),tt(nint),ww(nint))

      if (nint.eq.6) then
      rr(1)=0.16667
      ss(1)=0.16667
      tt(1)=-.57735
      ww(1)=1.
      rr(2)=0.66667
      ss(2)=0.16667
      tt(2)=-.57735
      ww(2)=1.
      rr(3)=0.16667
      ss(3)=0.66667
      tt(3)=-.57735
      ww(3)=1.
      rr(4)=0.16667
      ss(4)=0.16667
      tt(4)=.57735
      ww(4)=1.
      rr(5)=0.66667
      ss(5)=0.16667
      tt(5)=.57735
      ww(5)=1.
      rr(6)=0.16667
      ss(6)=0.66667
      tt(6)=.57735
      ww(6)=1.
      elseif (nint.eq.8) then
      rr(1)=-.57735
      ss(1)=-.57735
      tt(1)=-.57735
      ww(1)=1.
      rr(2)=.57735
      ss(2)=-.57735
      tt(2)=-.57735
      ww(2)=1.
      rr(3)=.57735
      ss(3)=.57735
      tt(3)=-.57735
      ww(3)=1.
      rr(4)=-.57735
      ss(4)=.57735
      tt(4)=-.57735
      ww(4)=1.
      rr(5)=-.57735
      ss(5)=-.57735
      tt(5)=.57735
      ww(5)=1.
      rr(6)=.57735
      ss(6)=-.57735
      tt(6)=.57735
      ww(6)=1.
      rr(7)=.57735
      ss(7)=.57735
      tt(7)=.57735
      ww(7)=1.
      rr(8)=-.57735
      ss(8)=.57735
      tt(8)=.57735
      ww(8)=1.
      endif

        do k=1,mpe
        x(k)=xg(icon(k))
        y(k)=yg(icon(k))
        z(k)=zg(icon(k))
        u(k)=ug(icon(k)) !VKP
        temp(k)=tempp(icon(k))
        rb(k)=rebound(icon(k))
        xp(k)=xgp(icon(k))
        yp(k)=ygp(icon(k))
        zp(k)=zgp(icon(k))
        up(k)=ugp(icon(k)) !VKP
        enddo

      ael=0.
      aelp=0.
      bel=0.

      timep=time-dt

        do i=1,mpe
        call find_velo (x(i),y(i),z(i),time,dt,fault,nfault,vv(i,1),vv(i,2),vv(i,3),xming,xmaxg,yming,ymaxg,.true.)
        enddo

        do iint=1,nint

        r=rr(iint)
        s=ss(iint)
        t=tt(iint)
        w=ww(iint)

          if (mpe.eq.8) then
          h(1)=(1.-r)*(1.-s)*(1.-t)/8.
          h(2)=(1.+r)*(1.-s)*(1.-t)/8.
          h(3)=(1.+r)*(1.+s)*(1.-t)/8.
          h(4)=(1.-r)*(1.+s)*(1.-t)/8.
          h(5)=(1.-r)*(1.-s)*(1.+t)/8.
          h(6)=(1.+r)*(1.-s)*(1.+t)/8.
          h(7)=(1.+r)*(1.+s)*(1.+t)/8.
          h(8)=(1.-r)*(1.+s)*(1.+t)/8.
          else
          h(1)=(1.-r-s)*(1.-t)/2.
          h(2)=r*(1.-t)/2.
          h(3)=s*(1.-t)/2.
          h(4)=(1.-r-s)*(1.+t)/2.
          h(5)=r*(1.+t)/2.
          h(6)=s*(1.+t)/2.
          endif

        xint=0.
        yint=0.
        zint=0.
        rint=0.
        uint=0.d0 !VKP
        tint=0.d0
          do k=1,mpe
          xint=xint+h(k)*x(k)
          yint=yint+h(k)*y(k)
          zint=zint+h(k)*z(k)
          rint=rint+h(k)*rb(k)
          uint=uint+h(k)*u(k) !VKP
          tint=tint+h(k)*temp(k)
          enddo

        veloz=rint/dt
        velou=uint !VKP

          if (mpe.eq.8) then
          dhdr(1)=-(1.-s)*(1.-t)/8.
          dhdr(2)=(1.-s)*(1.-t)/8.
          dhdr(3)=(1.+s)*(1.-t)/8.
          dhdr(4)=-(1.+s)*(1.-t)/8.
          dhdr(5)=-(1.-s)*(1.+t)/8.
          dhdr(6)=(1.-s)*(1.+t)/8.
          dhdr(7)=(1.+s)*(1.+t)/8.
          dhdr(8)=-(1.+s)*(1.+t)/8.
          dhds(1)=-(1.-r)*(1.-t)/8.
          dhds(2)=-(1.+r)*(1.-t)/8.
          dhds(3)=(1.+r)*(1.-t)/8.
          dhds(4)=(1.-r)*(1.-t)/8.
          dhds(5)=-(1.-r)*(1.+t)/8.
          dhds(6)=-(1.+r)*(1.+t)/8.
          dhds(7)=(1.+r)*(1.+t)/8.
          dhds(8)=(1.-r)*(1.+t)/8.
          dhdt(1)=-(1.-r)*(1.-s)/8.
          dhdt(2)=-(1.+r)*(1.-s)/8.
          dhdt(3)=-(1.+r)*(1.+s)/8.
          dhdt(4)=-(1.-r)*(1.+s)/8.
          dhdt(5)=(1.-r)*(1.-s)/8.
          dhdt(6)=(1.+r)*(1.-s)/8.
          dhdt(7)=(1.+r)*(1.+s)/8.
          dhdt(8)=(1.-r)*(1.+s)/8.
          else
          dhdr(1)=-(1.-t)/2.
          dhdr(2)=(1.-t)/2.
          dhdr(3)=0.
          dhdr(4)=-(1.+t)/2.
          dhdr(5)=(1.+t)/2.
          dhdr(6)=0.
          dhds(1)=-(1.-t)/2.
          dhds(2)=0.
          dhds(3)=(1.-t)/2.
          dhds(4)=-(1.+t)/2.
          dhds(5)=0.
          dhds(6)=(1.+t)/2.
          dhdt(1)=-(1.-r-s)/2.
          dhdt(2)=-r/2.
          dhdt(3)=-s/2.
          dhdt(4)=(1.-r-s)/2.
          dhdt(5)=r/2.
          dhdt(6)=s/2.
          endif

        jcb=0.
        jcbp=0.
          do k=1,mpe
          jcb(1,1)=jcb(1,1)+dhdr(k)*x(k)
          jcb(1,2)=jcb(1,2)+dhdr(k)*y(k)
          jcb(1,3)=jcb(1,3)+dhdr(k)*z(k)
          jcb(2,1)=jcb(2,1)+dhds(k)*x(k)
          jcb(2,2)=jcb(2,2)+dhds(k)*y(k)
          jcb(2,3)=jcb(2,3)+dhds(k)*z(k)
          jcb(3,1)=jcb(3,1)+dhdt(k)*x(k)
          jcb(3,2)=jcb(3,2)+dhdt(k)*y(k)
          jcb(3,3)=jcb(3,3)+dhdt(k)*z(k)
          jcbp(1,1)=jcbp(1,1)+dhdr(k)*xp(k)
          jcbp(1,2)=jcbp(1,2)+dhdr(k)*yp(k)
          jcbp(1,3)=jcbp(1,3)+dhdr(k)*zp(k)
          jcbp(2,1)=jcbp(2,1)+dhds(k)*xp(k)
          jcbp(2,2)=jcbp(2,2)+dhds(k)*yp(k)
          jcbp(2,3)=jcbp(2,3)+dhds(k)*zp(k)
          jcbp(3,1)=jcbp(3,1)+dhdt(k)*xp(k)
          jcbp(3,2)=jcbp(3,2)+dhdt(k)*yp(k)
          jcbp(3,3)=jcbp(3,3)+dhdt(k)*zp(k)
          enddo

        volume=jcb(1,1)*jcb(2,2)*jcb(3,3)+jcb(1,2)*jcb(2,3)*jcb(3,1) &
              +jcb(2,1)*jcb(3,2)*jcb(1,3) &
              -jcb(1,3)*jcb(2,2)*jcb(3,1)-jcb(1,2)*jcb(2,1)*jcb(3,3) &
              -jcb(2,3)*jcb(3,2)*jcb(1,1)
        volumep=jcbp(1,1)*jcbp(2,2)*jcbp(3,3)+jcbp(1,2)*jcbp(2,3)*jcbp(3,1) &
               +jcbp(2,1)*jcbp(3,2)*jcbp(1,3) &
               -jcbp(1,3)*jcbp(2,2)*jcbp(3,1)-jcbp(1,2)*jcbp(2,1)*jcbp(3,3) &
               -jcbp(2,3)*jcbp(3,2)*jcbp(1,1)

          if (volume.le.eps) then
          print*,'volume= ',volume
          stop 'Element bow-tied or collapsed'
          endif

          if (volumep.le.eps) then
          print*,'volumep= ',volumep
          stop 'Element bow-tied or collapsed'
          endif

        jcbi(1,1)=(jcb(2,2)*jcb(3,3)-jcb(2,3)*jcb(3,2))/volume
        jcbi(2,1)=(jcb(2,3)*jcb(3,1)-jcb(2,1)*jcb(3,3))/volume
        jcbi(3,1)=(jcb(2,1)*jcb(3,2)-jcb(2,2)*jcb(3,1))/volume
        jcbi(1,2)=(jcb(1,3)*jcb(3,2)-jcb(1,2)*jcb(3,3))/volume
        jcbi(2,2)=(jcb(1,1)*jcb(3,3)-jcb(1,3)*jcb(3,1))/volume
        jcbi(3,2)=(jcb(1,2)*jcb(3,1)-jcb(1,1)*jcb(3,2))/volume
        jcbi(1,3)=(jcb(1,2)*jcb(2,3)-jcb(1,3)*jcb(2,2))/volume
        jcbi(2,3)=(jcb(1,3)*jcb(2,1)-jcb(1,1)*jcb(2,3))/volume
        jcbi(3,3)=(jcb(1,1)*jcb(2,2)-jcb(1,2)*jcb(2,1))/volume
        jcbip(1,1)=(jcbp(2,2)*jcbp(3,3)-jcbp(2,3)*jcbp(3,2))/volumep
        jcbip(2,1)=(jcbp(2,3)*jcbp(3,1)-jcbp(2,1)*jcbp(3,3))/volumep
        jcbip(3,1)=(jcbp(2,1)*jcbp(3,2)-jcbp(2,2)*jcbp(3,1))/volumep
        jcbip(1,2)=(jcbp(1,3)*jcbp(3,2)-jcbp(1,2)*jcbp(3,3))/volumep
        jcbip(2,2)=(jcbp(1,1)*jcbp(3,3)-jcbp(1,3)*jcbp(3,1))/volumep
        jcbip(3,2)=(jcbp(1,2)*jcbp(3,1)-jcbp(1,1)*jcbp(3,2))/volumep
        jcbip(1,3)=(jcbp(1,2)*jcbp(2,3)-jcbp(1,3)*jcbp(2,2))/volumep
        jcbip(2,3)=(jcbp(1,3)*jcbp(2,1)-jcbp(1,1)*jcbp(2,3))/volumep
        jcbip(3,3)=(jcbp(1,1)*jcbp(2,2)-jcbp(1,2)*jcbp(2,1))/volumep

          do k=1,mpe
            do i=1,3
            b(k,i)=jcbi(i,1)*dhdr(k)+jcbi(i,2)*dhds(k)+jcbi(i,3)*dhdt(k)
            bp(k,i)=jcbip(i,1)*dhdr(k)+jcbip(i,2)*dhds(k)+jcbip(i,3)*dhdt(k)
            enddo
          enddo

          do j=1,3
            do i=1,3
            epsij(i,j)=0.
              do k=1,mpe
              epsij(i,j)=epsij(i,j)+b(k,i)*vv(k,j)
              enddo
      	    enddo
          enddo
        e2d=epsij(1,1)**2+epsij(2,2)**2+epsij(3,3)**2 &
           +(epsij(1,2)+epsij(2,1))**2/4.d0+(epsij(1,3)+epsij(3,1))**2/4.d0+(epsij(2,3)+epsij(3,2))**2/4.d0
        if (e2d.gt.0.d0) e2d=sqrt(e2d/2.d0)
        heat=heat+abs(e2d)*1.d8*friction/(specific_heat/3.d3)

        call find_velo (xint,yint,zint,time,dt,fault,nfault,vx,vy,vz,xming,xmaxg,yming,ymaxg,.true.)

        velox=vx
        veloy=vy
        veloz=veloz+vz+velou !VKP

        uvnorm=velox**2+veloy**2+veloz**2
        xmin=minval(x)
        xmax=maxval(x)
        dx=xmax-xmin
        ymin=minval(y)
        ymax=maxval(y)
        dy=ymax-ymin
        zmin=minval(z)
        zmax=maxval(z)
        dz=zmax-zmin
	  if (uvnorm.le.tiny(uvnorm)) then
	  tau=0.
	  else
    	  tau=(abs(velox)*dx+abs(veloy)*dy+abs(veloz)*dz)/sqrt(15.)/uvnorm
    	  endif

	  do k=1,mpe
	  ht(k)=h(k)+tau*(velox*b(k,1)+veloy*b(k,2)+veloz*b(k,3))
	  enddo

          if (heatproduction.lt.0.d0) then
          xinit=hinterpolate (lsfx,xint,yint,zint,xming,xmaxg,yming,ymaxg,zming,zmaxg,nx,ny,nz)
          yinit=hinterpolate (lsfy,xint,yint,zint,xming,xmaxg,yming,ymaxg,zming,zmaxg,nx,ny,nz)
          zinit=hinterpolate (lsf,xint,yint,zint,xming,xmaxg,yming,ymaxg,zming,zmaxg,nx,ny,nz)
          heat=heat_production(xinit,yinit,zinit,tint,tfinal-time)
          endif

          if (diffusivity.lt.0.d0) then
          xinit=hinterpolate (lsfx,xint,yint,zint,xming,xmaxg,yming,ymaxg,zming,zmaxg,nx,ny,nz)
          yinit=hinterpolate (lsfy,xint,yint,zint,xming,xmaxg,yming,ymaxg,zming,zmaxg,nx,ny,nz)
          zinit=hinterpolate (lsf,xint,yint,zint,xming,xmaxg,yming,ymaxg,zming,zmaxg,nx,ny,nz)
          cond=thermal_diffusivity (xinit,yinit,zinit,tint,tfinal-time)
          endif

          do j=1,mpe
          bel(j)=bel(j)+h(j)*dt*((1.-alpha)*heat*volumep+alpha*heat*volume)*w
            do i=1,mpe
            xcond=cond*(b(i,1)*b(j,1)+b(i,2)*b(j,2)+b(i,3)*b(j,3))*volume &
                 +rhoc*ht(i)*(velox*b(j,1)+veloy*b(j,2)+veloz*b(j,3))*volume
            xmass=h(i)*rhoc*h(j)*volume*dynamic
            xcondp=cond*(bp(i,1)*bp(j,1)+bp(i,2)*bp(j,2)+bp(i,3)*bp(j,3))*volumep &
                 +rhoc*ht(i)*(velox*bp(j,1)+veloy*bp(j,2)+veloz*bp(j,3))*volumep
            xmassp=h(i)*rhoc*h(j)*volumep*dynamic
            ael(i,j)=ael(i,j)+(xmass+dt*alpha*xcond)*w
            aelp(i,j)=aelp(i,j)+(xmass-dt*(1.-alpha)*xcondp)*w
            enddo
          enddo

        enddo

        do i=1,mpe
          do j=1,mpe
          bel(i)=bel(i)+aelp(i,j)*tempp(icon(j))
          enddo
        enddo

        do i=1,mpe
          if (kfix(icon(i)).eq.1) then
          fixt=tempp(icon(i))
            do j=1,mpe
            bel(j)=bel(j)-ael(j,i)*fixt
            ael(i,j)=0.
            ael(j,i)=0.
            enddo
          ael(i,i)=1.
          bel(i)=fixt
          endif
        enddo

      deallocate (rr,ss,tt,ww)
      deallocate (aelp,b,bp,vv)
      deallocate (x,y,z,xp,yp,zp,rb)
      deallocate (h,ht,dhdr,dhds,dhdt)

      return
      end
