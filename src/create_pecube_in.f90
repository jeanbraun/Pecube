      subroutine create_pecube_in (fault,nfault,nd,range,param,run,iproc)

      use Pecube
      use DEM

      implicit none

! Subroutine to create the input file Pecube.in to be used by Pecube

! The user can modify this file at will but the output it produces (Pecube.in)
! must obey rules that are described in the user guide and in Pecube.f90

      integer i1,j1,nfault,ieobs,ilog,icol,jd,jcol,iproc

      type (parameters) p
      type (faulttype) fault(nfault)

      double precision,dimension(:,:),allocatable::zNZ
      double precision,dimension(:),allocatable::z,xz,yz,zsmooth
      double precision,dimension(:),allocatable::timek,topomag,topooffset
      integer,dimension(:),allocatable::iout
      integer,dimension(:,:),allocatable::iconz,neighz
      integer i,j,k,nfnme,nx0,ny0,nskip,nstep,istep,isoflag,nxiso,nyiso,nz,ij,it
      integer nx,ny,nsurf,nelem,nne,nobsfile,iterative,interpol,icon1,icon2,icon3,icon4,nobs,nobs1,nobs2,nobs3,nobs4,nobs5,nobs6
      integer ftlflag,mftflag,fltflag,ageflag(9),nhist,nheating
      double precision dx,dy,ddx,ddy,xlon,xlat,tau,rhoc,rhom,young,poisson,thickness,tprevious
      double precision crustal_thickness,diffusivity,tmax,tmsl,tlapse,heatproduction,friction
      double precision xl,yl,zl,xlon1,xlat1,xlon2,xlat2,x,y,r,s
      double precision xlonobs,xlatobs,heightobs,ageheobs,dageheobs,ageftobs,dageftobs
      double precision wobs1,wobs2,wobs3,wobs4,a1,b1,c1,a2,b2,c2,a3,b3,c3,surf2
      double precision ageheZobs,dageheZobs,ageftZobs,dageftZobs,ageKarobs,dageKarobs
      double precision ageBarobs,dageBarobs,ageMarobs,dageMarobs,ageHarobs,dageHarobs
      double precision ftdist(20),grainsize
      double precision thist(1024),temphist(1024),errortemphist(1024)
      double precision theating(1024),released(1024),dreleased(1024),agereleased(1024),dagereleased(1024),duration(1024)
      double precision size43He,age43He,dage43He
      double precision doser,d0,radius,et,logs,b,logrho,nn,eu,sigmaet

      character run*5,fnme*300,obsfile*300,line*1024,c5*5,PecubeFnme*10,cproc*4

      real*4 range(2,*),param(*),xf
      integer nd,nd0,nc5
      logical vivi,xyz !VKP

      nd = 0
      call read_input_file (run//'/input/Pecube.in', 0, p, nd, range, param)

      p%run_name = run

      fnme = p%topo_file_name

        do i=1,300
        if (fnme(i:i).ne.' ') nfnme=i
        enddo
      vivi=.FALSE. !VKP
      if (fnme(nfnme:nfnme).eq.'/') vivi=.TRUE. !VKP
      if (vivi) nfnme=nfnme-1 !VKP
      xyz=.false.
      if (fnme(nfnme:nfnme).eq.'$') xyz=.TRUE.
      if (xyz) nfnme=nfnme-1

! nx0 and ny0 are the number of points in the input topographic file
! dx and dy are the data spacings in the x- and y-direction in the input topographic
! file (in degrees)
! nskip is the number of points that are skipped when sampling the initial
! topo file
! xlat and xlon are the latitude/longitude of the bottom left corner of the
! data set (in deg.)
      nx0 = p%nx
      ny0 = p%ny
      dx = p%dlon
      dy = p%dlat
      nskip = p%nskip
      xlon = p%lon0
      xlat = p%lat0

! nstep is the number of time step
! tau is the erosion time scale (assuming an exponential decvrease in topography) (in Myr)
      nstep = p%ntime
      tau = p%erosional_time_scale

! for each time step + 1
! timek is the time (in Myr)
! topomag is the topographic amplification factor at this step
! topooffset is the topographic offset
      allocate (timek(nstep+1),topomag(nstep+1),topooffset(nstep+1),iout(nstep+1))

        do istep=1,nstep+1
        timek(istep) = p%time_topo(istep)
        topomag(istep) = p%amplification(istep)
        topooffset(istep) = p%offset(istep)
        iout(istep) = p%output(istep)
        enddo

! converts geological time into model time
        do istep=nstep+1,1,-1
        timek(istep)=timek(1)-timek(istep)
        enddo

! isostasy flag (0 no isostasy, 1 isostasy on)
! rhoc and rhom are the densities for the crust and mantle, respectively (in kg/m3)
! these values are used in the isostatic calculations
! young is the elastic plate young modulus (in Pa)
! poisson is poisson's ratio (dimensionless)
! thickness is the elastic thickness of the plate (in km)
! nxiso and nyiso are the resolutions in the x- and y-directions of the grid on
! which the isostatic (flexural) calculations are performed (including the FFT)
! note that these numbers must be powers of two.
      isoflag = p%isostasy
      rhoc = p%rho_crust
      rhom = p%rho_asthenosphere
      young = p%young_modulus
      poisson = p%poisson_ratio
      thickness = p%EET
      nxiso = p%nx_isostasy
      nyiso = p%ny_isostasy

! crustal thickness is the averaged crustal thickness (i.e. the depth at which the
! temperature is assumed to be constant) (in km)
! nz is the number of points in the z-direction
! diffusivity is the heat diffusivity (in km2/Myr)
! tmax is the basal temperature (in C)
! tmsl is the temperature at the top of the model (at z=0)
! tlapse is the lapse rate (or change of temperature with height in the atmosphere)
! (in C/km)
! heatproduction is the rate of heat production (in C/Myr)
      crustal_thickness = p%thickness
      nz = p%nz
      diffusivity = p%thermal_diffusivity
      tmax = p%basal_temperature
      tmsl = p%sea_level_temperature
      tlapse = p%lapse_rate
      heatproduction = p%heat_production

      crustal_thickness=crustal_thickness*1.d3

! obsfile is the name of the observation file
        do i=1,300
        obsfile(i:i)=' '
        enddo
      obsfile = p%data_folder
        do i=1,300
        if (obsfile(i:i).ne.' ') nobsfile=i
        enddo

! new parameters added in Pecube2

      tprevious=timek(nstep+1)
      ftlflag=0
      mftflag=0
      fltflag=0
      friction=0.d0
      ageflag=1

      tprevious = p%default_age
      if (tprevious.lt.0.) tprevious=timek(nstep+1)

      ftlflag = p%FT_code_flag

      mftflag = p%misfit_slope

      fltflag = p%fault_advect_flag

      friction = p%shear_heating

      ageflag(1) = p%age_AHe_flag
      ageflag(2) = p%age_ZHe_flag
      ageflag(3) = p%age_AFT_flag
      ageflag(4) = p%age_ZFT_flag
      ageflag(5) = p%age_KAr_flag
      ageflag(6) = p%age_BAr_flag
      ageflag(7) = p%age_MAr_flag
      ageflag(8) = p%age_HAr_flag
      ageflag(9) = p%age_FTL_flag

! reads in topography

if (.not.vivi) then !VKP

      if (nx0.gt.0) then

      allocate (zNZ(nx0,ny0))
        if (fnme(1:nfnme).eq.'Nil') then
        zNZ=0.d0
        elseif (fnme(1:nfnme).eq.'Topo30') then
        dx = 360.d0/43200
        dy = dx
        write (cproc,'(i4)') iproc
        if (iproc.lt.10) cproc(1:3)='000'
        if (iproc.lt.100) cproc(1:2)='00'
        if (iproc.lt.1000) cproc(1:1)='0'
        PecubeFnme = trim(fnme(1:nfnme)//cproc)
        call ExtractDEM (xlon, xlat, nx0, ny0, PecubeFnme = run//'/data/'//trim(PecubeFnme))
        open (8,file=run//'/data/'//trim(PecubeFnme)//'.dat',status='old')
        read (8,*) zNZ
        close (8)
        call system ('rm '//run//'/data/'//trim(PecubeFnme)//'.dat')
        else
        open (8,file=run//'/data/'//fnme(1:nfnme),status='old')
          if (xyz) then
            do j=ny0,1,-1
              do i=1,nx0
              read (8,*) x,y,zNZ(i,j)
              enddo
            enddo
          else
          read (8,*) zNZ
          endif
        close (8)
        endif

      nx=(nx0-1)/nskip+1
      ny=(ny0-1)/nskip+1
      allocate (z(nx*ny))
      ij=0
        do j=1,ny0,nskip
          do i=1,nx0,nskip
          ij=ij+1
          z(ij)=zNZ(i,j)
          enddo
        enddo

      deallocate (zNZ)

      xlon1=xlon
      xlat1=xlat
      xlon2=xlon+(nx-1)*dx*nskip
      xlat2=xlat+(ny-1)*dy*nskip

      else

      allocate (z(-nx0),xz(-nx0),yz(-nx0),iconz(3,-ny0))
        open (8,file=run//'/data/'//fnme(1:nfnme),status='old')
          do i=1,-nx0
          read (8,*) z(i)
          enddo
        close (8)
        open (8,file=run//'/data/'//fnme(1:nfnme)//'.geometry',status='old')
          do i=1,-nx0
          read (8,*) xz(i),yz(i)
          enddo
          do i=1,-ny0
          read (8,*) (iconz(k,i),k=1,3)
          enddo
        close (8)

      xlon1=minval(xz)
      xlat1=minval(yz)
      xlon2=maxval(xz)
      xlat2=maxval(yz)
      xlon=xlon1
      xlat=xlat1

      endif

else !VKP

      if (nx0.gt.0) then

      nx=(nx0-1)/nskip+1 !VKP
      ny=(ny0-1)/nskip+1 !VKP
      allocate (z(nx*ny)) !VKP
      topomag=1.d0 !VKP
      topooffset=0.d0 !VKP

      xlon1=xlon
      xlat1=xlat
      xlon2=xlon+(nx-1)*dx*nskip
      xlat2=xlat+(ny-1)*dy*nskip

      else

      allocate (z(-nx0),xz(-nx0),yz(-nx0),iconz(3,-ny0))
        open (8,file=run//'/data/'//fnme(1:nfnme)//'/geometry',status='old')
          do i=1,-nx0
          read (8,*) xz(i),yz(i)
          enddo
          do i=1,-ny0
          read (8,*) (iconz(k,i),k=1,3)
          enddo
        close (8)
      topomag=1.d0 !VKP
      topooffset=0.d0 !VKP
      xlon1=minval(xz)
      xlat1=minval(yz)
      xlon2=maxval(xz)
      xlat2=maxval(yz)
      xlon=xlon1
      xlat=xlat1

      endif

endif !VKP

      if (nx0.gt.0) then

      xl=dx*(nx-1)*nskip*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
      yl=dy*(ny-1)*nskip*111.11
      zl=crustal_thickness/1.e3
      z=z/crustal_thickness*zl

      else

      zl=crustal_thickness/1.e3
      z=z/crustal_thickness*zl
      xl=(xlon2-xlon1)*111.11*cos((xlat+(xlat2-xlat1)/2.)*3.141592654/180.)
      yl=(xlat2-xlat1)*111.11

      endif

! reads in fault definitions

      call read_in_fault_parameters (fault,nfault,xlon1,xlat1,xlon2,xlat2,xl,yl,zl,timek(nstep+1), &
                                     nd,range,param,run)

      if (nd0.gt.0) nd=nd0

      open (7,status='scratch')
!      open (7,file='JUNK.txt')

      ilog=p%debug
      iterative=1
      interpol=1
      nsurf=nx*ny
      nelem=(nx-1)*(ny-1)
      nne=4
      ddx=xl/(nx-1)*1.d3
      ddy=yl/(ny-1)*1.d3
        if (nx0.lt.0) then
        nsurf=-nx0
        nelem=-ny0
        nne=3
        ddx=xl/sqrt(dble(nsurf))*1.d3
        ddy=yl/sqrt(dble(nsurf))*1.d3
	endif

      write (7,'(a)') run
        if (vivi) then !VKP
        write (7,*) nne,-nsurf,nz,nelem,zl,diffusivity,heatproduction,friction !VKP
        else !VKP
        write (7,*) nne,nsurf,nz,nelem,zl,diffusivity,heatproduction,friction
        endif !VKP
      write (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol,tprevious,ftlflag,mftflag,fltflag
      write (7,*) isoflag,tau,rhoc,rhom
      write (7,*) nx,ny,nxiso,nyiso
      write (7,*) ddx,ddy,young,poisson,thickness*1.d3
      write (7,*) xlon1,xlon2,xlat1,xlat2

      if (nx0.gt.0) then

        do j=1,ny
          do i=1,nx
          x=xl*float(i-1)/float(nx-1)
          y=yl*float(j-1)/float(ny-1)
          write (7,*) x,y
          enddo
        enddo

        do j=1,ny-1
          do i=1,nx-1
          icon1=(j-1)*nx+i
          icon2=icon1+1
          icon3=icon1+nx+1
          icon4=icon1+nx
          write (7,*) icon1,icon2,icon3,icon4
          enddo
        enddo

      else

        do i=1,nsurf
        x=(xz(i)-xlon1)*111.11*cos((xlat+(xlat2-xlat1)/2.)*3.141592654/180.)
        y=(yz(i)-xlat1)*111.11
        write (7,*) x,y
        enddo

        do i=1,nelem
        write (7,*) (iconz(k,i),k=1,3)
        enddo

      endif

        do istep=0,nstep
        write (7,*) timek(istep+1),iout(istep+1)
          if (vivi) then !VKP
            if (istep.lt.10) then !VKP
            write (c5(1:1),'(i1)') istep !VKP
            nc5=1 !VKP
            elseif (istep.lt.100) then !VKP
            write (c5(1:2),'(i2)') istep !VKP
            nc5=2 !VKP
            elseif (istep.lt.1000) then !VKP
            write (c5(1:3),'(i3)') istep !VKP
            nc5=3 !VKP
            elseif (istep.lt.10000) then !VKP
            write (c5(1:4),'(i4)') istep !VKP
            nc5=4 !VKP
            else !VKP
            write (c5(1:5),'(i5)') istep !VKP
            nc5=5 !VKP
            endif !VKP
            if (nx0.gt.0) then
            allocate (zNZ(nx0,ny0))
            open (67,file=run//'/data/'//fnme(1:nfnme)//'/topo'//c5(1:nc5),status='old') !VKP
            read (67,*) zNZ !VKP
            ij=0 !VKP
              do j=1,ny0,nskip  !VKP
                do i=1,nx0,nskip  !VKP
                ij=ij+1 !VKP
                z(ij)=zNZ(i,j) !VKP
                enddo !VKP
              enddo !VKP
            z=z/1.e3 !VKP
            close (67) !VKP
            else
            open (67,file=run//'/data/'//fnme(1:nfnme)//'/topo'//c5(1:nc5),status='old') !VKP
            read (67,*) z
            z=z/1.e3
            close (67)
            endif
          endif !VKP
          if (topomag(istep+1).lt.0.) then
            allocate(zsmooth(nsurf))
            call smooth_topo(z, zsmooth, nx, ny, int(-topomag(istep+1)), int(topooffset(istep+1)))
            write (7,*) (zsmooth(k), k=1,nsurf)
            deallocate(zsmooth)
          else
            write (7,*) (z(k)*topomag(istep+1)+topooffset(istep+1),k=1,nsurf)
          endif
          if (vivi) then !VKP
            if (nx0.gt.0) then
            open (67,file=run//'/data/'//fnme(1:nfnme)//'/uplift'//c5(1:nc5),status='old') !VKP
            read (67,*) zNZ !VKP
            ij=0 !VKP
              do j=1,ny0,nskip !VKP
                do i=1,nx0,nskip !VKP
                ij=ij+1 !VKP
                z(ij)=zNZ(i,j) !VKP
                enddo !VKP
              enddo !VKP
            close (67) !VKP
            write (7,*) (z(k),k=1,nsurf) !VKP
            open (67,file=run//'/data/'//fnme(1:nfnme)//'/temp'//c5(1:nc5),status='old') !VKP
            read (67,*) zNZ !VKP
            ij=0 !VKP
              do j=1,ny0,nskip !VKP
                do i=1,nx0,nskip !VKP
                ij=ij+1 !VKP
                z(ij)=zNZ(i,j) !VKP
                enddo !VKP
              enddo !VKP
            close (67) !VKP
            write (7,*) (z(k),k=1,nsurf) !VKP
            deallocate (zNZ)
            else
            open (67,file=run//'/data/'//fnme(1:nfnme)//'/uplift'//c5(1:nc5),status='old') !VKP
            read (67,*) z
            close (67)
            write (7,*) (z(k),k=1,nsurf) !VKP
            open (67,file=run//'/data/'//fnme(1:nfnme)//'/temp'//c5(1:nc5),status='old') !VKP
            read (67,*) z
            close (67)
            write (7,*) (z(k),k=1,nsurf)
            endif
          endif !VKP
        enddo

! observations

      if (obsfile(1:nobsfile).eq.'Nil') then
      nobs=0
      write (7,*) nobs,0,0,0
      else
      call read_data_folder (run//'/data/'//obsfile(1:nobsfile), .true., xlon1, xlon2, xlat1, xlat2, iproc, nd, &
                             nobs1, nobs2, nobs3, nobs4, nobs5, nobs6)
        if (nx0.lt.0) then
        allocate (neighz(3,nelem))
        call neighbours (iconz,neighz,nelem)
        it=1
        endif
      write (cproc,'(i4)') iproc
      if (iproc.lt.10) cproc(1:3)='000'
      if (iproc.lt.100) cproc(1:2)='00'
      if (iproc.lt.1000) cproc(1:1)='0'
      open (8,file=run//'/data/'//obsfile(1:nobsfile)//cproc//'.txt',status='old')
      if (p%save_PTT_paths.ne.0) nobs1 = -nobs1
      write (7,*) nobs1,nobs2,nobs3,nobs4,nobs5,nobs6
      if (nobs1.lt.0) nobs1=-nobs1
        do i=1,nobs1
        read (8,*) xlonobs,xlatobs,heightobs,ageheobs,dageheobs,ageftobs,dageftobs, &
                   ageheZobs,dageheZobs,ageftZobs,dageftZobs,ageKarobs,dageKarobs, &
                   ageBarobs,dageBarobs,ageMarobs,dageMarobs,ageHarobs,dageHarobs, &
                   ftdist,grainsize
        if (nx0.gt.0) then
        i1=int((xlonobs-xlon)/(dx*nskip))+1 !VKP
        if (i1.eq.nx) i1=nx-1
        j1=int((xlatobs-xlat)/(dy*nskip))+1 !VKP
        if (j1.eq.ny) j1=ny-1
        ieobs=i1+(j1-1)*(nx-1)
        r=(xlonobs-(i1-1)*dx*nskip-xlon)/(dx*nskip) !VKP
        r=-1.+2.*r
        s=(xlatobs-(j1-1)*dy*nskip-xlat)/(dy*nskip) !VKP
        s=-1.+2.*s
        wobs1=(1.-r)*(1.-s)/4.
        wobs2=(1.+r)*(1.-s)/4.
        wobs3=(1.+r)*(1.+s)/4.
        wobs4=(1.-r)*(1.+s)/4.
        else
        call find_triangle (xlonobs,xlatobs,xz,yz,iconz,neighz,nelem,nsurf,it)
        ieobs=it
        surf2=xz(iconz(1,it))*yz(iconz(2,it))+xz(iconz(2,it))*yz(iconz(3,it))+xz(iconz(3,it))*yz(iconz(1,it)) &
             -yz(iconz(1,it))*xz(iconz(2,it))-yz(iconz(2,it))*xz(iconz(3,it))-yz(iconz(3,it))*xz(iconz(1,it))
        a1=xz(iconz(2,it))*yz(iconz(3,it))-xz(iconz(3,it))*yz(iconz(2,it))
        b1=yz(iconz(2,it))-yz(iconz(3,it))
        c1=xz(iconz(3,it))-xz(iconz(2,it))
        a2=xz(iconz(3,it))*yz(iconz(1,it))-xz(iconz(1,it))*yz(iconz(3,it))
        b2=yz(iconz(3,it))-yz(iconz(1,it))
        c2=xz(iconz(1,it))-xz(iconz(3,it))
        a3=xz(iconz(1,it))*yz(iconz(2,it))-xz(iconz(2,it))*yz(iconz(1,it))
        b3=yz(iconz(1,it))-yz(iconz(2,it))
        c3=xz(iconz(2,it))-xz(iconz(1,it))
        wobs1=(a1+b1*xlonobs+c1*xlatobs)/surf2
        wobs2=(a2+b2*xlonobs+c2*xlatobs)/surf2
        wobs3=(a3+b3*xlonobs+c3*xlatobs)/surf2
        wobs4=0.
        endif
        write (7,*) xlonobs,xlatobs,ageheobs,dageheobs,ageftobs,dageftobs, &
                    ageheZobs,dageheZobs,ageftZobs,dageftZobs,ageKarobs,dageKarobs, &
                    ageBarobs,dageBarobs,ageMarobs,dageMarobs,ageHarobs,dageHarobs, &
                    ftdist,heightobs,grainsize,ieobs,wobs1,wobs2,wobs3,wobs4
        enddo
! cooling curves
        do i=1,nobs2
        read (8,*) xlonobs,xlatobs,heightobs,nhist, &
                   (thist(j),temphist(j),errortemphist(j),j=1,nhist)
        if (nx0.gt.0) then
        i1=int((xlonobs-xlon)/(dx*nskip))+1 !VKP
        if (i1.eq.nx) i1=nx-1
        j1=int((xlatobs-xlat)/(dy*nskip))+1 !VKP
        if (j1.eq.ny) j1=ny-1
        ieobs=i1+(j1-1)*(nx-1)
        r=(xlonobs-(i1-1)*dx*nskip-xlon)/(dx*nskip) !VKP
        r=-1.+2.*r
        s=(xlatobs-(j1-1)*dy*nskip-xlat)/(dy*nskip) !VKP
        s=-1.+2.*s
        wobs1=(1.-r)*(1.-s)/4.
        wobs2=(1.+r)*(1.-s)/4.
        wobs3=(1.+r)*(1.+s)/4.
        wobs4=(1.-r)*(1.+s)/4.
        else
        call find_triangle (xlonobs,xlatobs,xz,yz,iconz,neighz,nelem,nsurf,it)
        ieobs=it
        surf2=xz(iconz(1,it))*yz(iconz(2,it))+xz(iconz(2,it))*yz(iconz(3,it))+xz(iconz(3,it))*yz(iconz(1,it)) &
             -yz(iconz(1,it))*xz(iconz(2,it))-yz(iconz(2,it))*xz(iconz(3,it))-yz(iconz(3,it))*xz(iconz(1,it))
        a1=xz(iconz(2,it))*yz(iconz(3,it))-xz(iconz(3,it))*yz(iconz(2,it))
        b1=yz(iconz(2,it))-yz(iconz(3,it))
        c1=xz(iconz(3,it))-xz(iconz(2,it))
        a2=xz(iconz(3,it))*yz(iconz(1,it))-xz(iconz(1,it))*yz(iconz(3,it))
        b2=yz(iconz(3,it))-yz(iconz(1,it))
        c2=xz(iconz(1,it))-xz(iconz(3,it))
        a3=xz(iconz(1,it))*yz(iconz(2,it))-xz(iconz(2,it))*yz(iconz(1,it))
        b3=yz(iconz(1,it))-yz(iconz(2,it))
        c3=xz(iconz(2,it))-xz(iconz(1,it))
        wobs1=(a1+b1*xlonobs+c1*xlatobs)/surf2
        wobs2=(a2+b2*xlonobs+c2*xlatobs)/surf2
        wobs3=(a3+b3*xlonobs+c3*xlatobs)/surf2
        wobs4=0.
        endif
        write (7,*) xlonobs,xlatobs,nhist,(thist(j),temphist(j),errortemphist(j),j=1,nhist), &
                    heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
        enddo
! 4-3He data
        do i=1,nobs3
        read (8,*) xlonobs,xlatobs,heightobs,size43He,age43He,dage43He,nheating, &
                   (theating(j),duration(j),released(j),dreleased(j), &
                    agereleased(j),dagereleased(j),j=1,nheating)
        if (nx0.gt.0) then
        i1=int((xlonobs-xlon)/(dx*nskip))+1 !VKP
        if (i1.eq.nx) i1=nx-1
        j1=int((xlatobs-xlat)/(dy*nskip))+1 !VKP
        if (j1.eq.ny) j1=ny-1
        ieobs=i1+(j1-1)*(nx-1)
        r=(xlonobs-(i1-1)*dx*nskip-xlon)/(dx*nskip) !VKP
        r=-1.+2.*r
        s=(xlatobs-(j1-1)*dy*nskip-xlat)/(dy*nskip) !VKP
        s=-1.+2.*s
        wobs1=(1.-r)*(1.-s)/4.
        wobs2=(1.+r)*(1.-s)/4.
        wobs3=(1.+r)*(1.+s)/4.
        wobs4=(1.-r)*(1.+s)/4.
        else
        call find_triangle (xlonobs,xlatobs,xz,yz,iconz,neighz,nelem,nsurf,it)
        ieobs=it
        surf2=xz(iconz(1,it))*yz(iconz(2,it))+xz(iconz(2,it))*yz(iconz(3,it))+xz(iconz(3,it))*yz(iconz(1,it)) &
             -yz(iconz(1,it))*xz(iconz(2,it))-yz(iconz(2,it))*xz(iconz(3,it))-yz(iconz(3,it))*xz(iconz(1,it))
        a1=xz(iconz(2,it))*yz(iconz(3,it))-xz(iconz(3,it))*yz(iconz(2,it))
        b1=yz(iconz(2,it))-yz(iconz(3,it))
        c1=xz(iconz(3,it))-xz(iconz(2,it))
        a2=xz(iconz(3,it))*yz(iconz(1,it))-xz(iconz(1,it))*yz(iconz(3,it))
        b2=yz(iconz(3,it))-yz(iconz(1,it))
        c2=xz(iconz(1,it))-xz(iconz(3,it))
        a3=xz(iconz(1,it))*yz(iconz(2,it))-xz(iconz(2,it))*yz(iconz(1,it))
        b3=yz(iconz(1,it))-yz(iconz(2,it))
        c3=xz(iconz(2,it))-xz(iconz(1,it))
        wobs1=(a1+b1*xlonobs+c1*xlatobs)/surf2
        wobs2=(a2+b2*xlonobs+c2*xlatobs)/surf2
        wobs3=(a3+b3*xlonobs+c3*xlatobs)/surf2
        wobs4=0.
        endif
        dreleased=max(0.01d0,dreleased)
        dagereleased=max(0.1d0,dagereleased)
        write (7,*) xlonobs,xlatobs,size43He,age43He,dage43He,nheating, &
                   (theating(j),duration(j),released(j),dreleased(j), &
                    agereleased(j),dagereleased(j),j=1,nheating), &
                    heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
        enddo
! TSL data
        do i=1,nobs4
        read (8,*) xlonobs,xlatobs,heightobs,doser,d0,radius,et,logs,b,logrho,nn
        if (nx0.gt.0) then
        i1=int((xlonobs-xlon)/(dx*nskip))+1 !VKP
        if (i1.eq.nx) i1=nx-1
        j1=int((xlatobs-xlat)/(dy*nskip))+1 !VKP
        if (j1.eq.ny) j1=ny-1
        ieobs=i1+(j1-1)*(nx-1)
        r=(xlonobs-(i1-1)*dx*nskip-xlon)/(dx*nskip) !VKP
        r=-1.+2.*r
        s=(xlatobs-(j1-1)*dy*nskip-xlat)/(dy*nskip) !VKP
        s=-1.+2.*s
        wobs1=(1.-r)*(1.-s)/4.
        wobs2=(1.+r)*(1.-s)/4.
        wobs3=(1.+r)*(1.+s)/4.
        wobs4=(1.-r)*(1.+s)/4.
        else
        call find_triangle (xlonobs,xlatobs,xz,yz,iconz,neighz,nelem,nsurf,it)
        ieobs=it
        surf2=xz(iconz(1,it))*yz(iconz(2,it))+xz(iconz(2,it))*yz(iconz(3,it))+xz(iconz(3,it))*yz(iconz(1,it)) &
          -yz(iconz(1,it))*xz(iconz(2,it))-yz(iconz(2,it))*xz(iconz(3,it))-yz(iconz(3,it))*xz(iconz(1,it))
        a1=xz(iconz(2,it))*yz(iconz(3,it))-xz(iconz(3,it))*yz(iconz(2,it))
        b1=yz(iconz(2,it))-yz(iconz(3,it))
        c1=xz(iconz(3,it))-xz(iconz(2,it))
        a2=xz(iconz(3,it))*yz(iconz(1,it))-xz(iconz(1,it))*yz(iconz(3,it))
        b2=yz(iconz(3,it))-yz(iconz(1,it))
        c2=xz(iconz(1,it))-xz(iconz(3,it))
        a3=xz(iconz(1,it))*yz(iconz(2,it))-xz(iconz(2,it))*yz(iconz(1,it))
        b3=yz(iconz(1,it))-yz(iconz(2,it))
        c3=xz(iconz(2,it))-xz(iconz(1,it))
        wobs1=(a1+b1*xlonobs+c1*xlatobs)/surf2
        wobs2=(a2+b2*xlonobs+c2*xlatobs)/surf2
        wobs3=(a3+b3*xlonobs+c3*xlatobs)/surf2
        wobs4=0.
        endif
        dreleased=max(0.01d0,dreleased)
        dagereleased=max(0.1d0,dagereleased)
        write (7,*) xlonobs,xlatobs, &
          doser,d0,radius,et,logs,b,logrho,nn, &
          heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
        enddo

! OSL data
        do i=1,nobs5
          read (8,*) xlonobs,xlatobs,heightobs,doser,d0,et,logs,logrho,eu,nn
          if (nx0.gt.0) then
          i1=int((xlonobs-xlon)/(dx*nskip))+1 !VKP
          if (i1.eq.nx) i1=nx-1
          j1=int((xlatobs-xlat)/(dy*nskip))+1 !VKP
          if (j1.eq.ny) j1=ny-1
          ieobs=i1+(j1-1)*(nx-1)
          r=(xlonobs-(i1-1)*dx*nskip-xlon)/(dx*nskip) !VKP
          r=-1.+2.*r
          s=(xlatobs-(j1-1)*dy*nskip-xlat)/(dy*nskip) !VKP
          s=-1.+2.*s
          wobs1=(1.-r)*(1.-s)/4.
          wobs2=(1.+r)*(1.-s)/4.
          wobs3=(1.+r)*(1.+s)/4.
          wobs4=(1.-r)*(1.+s)/4.
          else
          call find_triangle (xlonobs,xlatobs,xz,yz,iconz,neighz,nelem,nsurf,it)
          ieobs=it
          surf2=xz(iconz(1,it))*yz(iconz(2,it))+xz(iconz(2,it))*yz(iconz(3,it))+xz(iconz(3,it))*yz(iconz(1,it)) &
            -yz(iconz(1,it))*xz(iconz(2,it))-yz(iconz(2,it))*xz(iconz(3,it))-yz(iconz(3,it))*xz(iconz(1,it))
          a1=xz(iconz(2,it))*yz(iconz(3,it))-xz(iconz(3,it))*yz(iconz(2,it))
          b1=yz(iconz(2,it))-yz(iconz(3,it))
          c1=xz(iconz(3,it))-xz(iconz(2,it))
          a2=xz(iconz(3,it))*yz(iconz(1,it))-xz(iconz(1,it))*yz(iconz(3,it))
          b2=yz(iconz(3,it))-yz(iconz(1,it))
          c2=xz(iconz(1,it))-xz(iconz(3,it))
          a3=xz(iconz(1,it))*yz(iconz(2,it))-xz(iconz(2,it))*yz(iconz(1,it))
          b3=yz(iconz(1,it))-yz(iconz(2,it))
          c3=xz(iconz(2,it))-xz(iconz(1,it))
          wobs1=(a1+b1*xlonobs+c1*xlatobs)/surf2
          wobs2=(a2+b2*xlonobs+c2*xlatobs)/surf2
          wobs3=(a3+b3*xlonobs+c3*xlatobs)/surf2
          wobs4=0.
          endif
          dreleased=max(0.01d0,dreleased)
          dagereleased=max(0.1d0,dagereleased)
          write (7,*) xlonobs,xlatobs, &
            doser,d0,et,logs,logrho,eu,nn, &
            heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
          enddo
  
! ESR data
          do i=1,nobs6
            read (8,*) xlonobs,xlatobs,heightobs,doser,d0,logs,et,sigmaet,nn
            if (nx0.gt.0) then
            i1=int((xlonobs-xlon)/(dx*nskip))+1 !VKP
            if (i1.eq.nx) i1=nx-1
            j1=int((xlatobs-xlat)/(dy*nskip))+1 !VKP
            if (j1.eq.ny) j1=ny-1
            ieobs=i1+(j1-1)*(nx-1)
            r=(xlonobs-(i1-1)*dx*nskip-xlon)/(dx*nskip) !VKP
            r=-1.+2.*r
            s=(xlatobs-(j1-1)*dy*nskip-xlat)/(dy*nskip) !VKP
            s=-1.+2.*s
            wobs1=(1.-r)*(1.-s)/4.
            wobs2=(1.+r)*(1.-s)/4.
            wobs3=(1.+r)*(1.+s)/4.
            wobs4=(1.-r)*(1.+s)/4.
            else
            call find_triangle (xlonobs,xlatobs,xz,yz,iconz,neighz,nelem,nsurf,it)
            ieobs=it
            surf2=xz(iconz(1,it))*yz(iconz(2,it))+xz(iconz(2,it))*yz(iconz(3,it))+xz(iconz(3,it))*yz(iconz(1,it)) &
              -yz(iconz(1,it))*xz(iconz(2,it))-yz(iconz(2,it))*xz(iconz(3,it))-yz(iconz(3,it))*xz(iconz(1,it))
            a1=xz(iconz(2,it))*yz(iconz(3,it))-xz(iconz(3,it))*yz(iconz(2,it))
            b1=yz(iconz(2,it))-yz(iconz(3,it))
            c1=xz(iconz(3,it))-xz(iconz(2,it))
            a2=xz(iconz(3,it))*yz(iconz(1,it))-xz(iconz(1,it))*yz(iconz(3,it))
            b2=yz(iconz(3,it))-yz(iconz(1,it))
            c2=xz(iconz(1,it))-xz(iconz(3,it))
            a3=xz(iconz(1,it))*yz(iconz(2,it))-xz(iconz(2,it))*yz(iconz(1,it))
            b3=yz(iconz(1,it))-yz(iconz(2,it))
            c3=xz(iconz(2,it))-xz(iconz(1,it))
            wobs1=(a1+b1*xlonobs+c1*xlatobs)/surf2
            wobs2=(a2+b2*xlonobs+c2*xlatobs)/surf2
            wobs3=(a3+b3*xlonobs+c3*xlatobs)/surf2
            wobs4=0.
            endif
            dreleased=max(0.01d0,dreleased)
            dagereleased=max(0.1d0,dagereleased)
            write (7,*) xlonobs,xlatobs, &
              doser,d0,logs,et,sigmaet,nn, &
              heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
            enddo
    
  
      close (8)
      call system ('rm '//run//'/data/'//obsfile(1:nobsfile)//cproc//'.txt')
      endif
      write (7,*) ageflag

      deallocate (z)
      if (nx0.lt.0) deallocate (xz,yz,iconz)
      deallocate (timek,topomag,topooffset)

      return
      end
