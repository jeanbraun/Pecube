      subroutine create_pecube_in (fault,nfault,nd,range,param,run,iproc)

      use Pecube
      use DEM

      implicit none

! Subroutine to create the input file Pecube.in to be used by Pecube

! The user can modify this file at will but the output it produces (Pecube.in)
! must obey rules that are described in the user guide and in Pecube.f90
! Also modify Parameter structure in module_Pecube.f90

! Parameters for sample-specific predictions added by Maxime Bernard: 02/06/2021

      integer i1,j1,nfault,ieobs,ilog,icol,jd,jcol,iproc

      type (parameters) p
      type (faulttype) fault(nfault)

      double precision,dimension(:,:),allocatable::zNZ
      double precision,dimension(:),allocatable::z,ztemp,xz,yz,zsmooth, Zincision, amplification, NewTopo
      double precision,dimension(:),allocatable::timek_temp2,topomag_temp2,topooffset_temp2,x_array,phase_array_temp2
      double precision,dimension(:),allocatable::timek_temp,topomag_temp,topooffset_temp,phase_array_temp
      double precision,dimension(:),allocatable::timek,topomag,topooffset,phase_array
      integer,dimension(:),allocatable::iout,index_array,iout_temp,index_array2,iout_temp2
      integer,dimension(:,:),allocatable::iconz,neighz
      integer i,j,k,nfnme,nx0,ny0,nskip,nstep,istep,isoflag,nxiso,nyiso,nz,ij,it
      integer nx,ny,nsurf,nelem,nne,nobsfile,iterative,interpol,icon1,icon2,icon3,icon4,nobs,nobs1,nobs2,nobs3,nobs4,nobs5,nobs6
      integer ftlflag,mftflag,fltflag,ageflag(9),nhist,nheating
      integer agecomput, Headward,res ! Maxime
    integer heatprod_flag
      double precision topoRef !Maxime
      double precision dx,dy,ddx,ddy,xlon,xlat,tau,rhoc,rhom,young,poisson,thickness,tprevious
      double precision crustal_thickness,diffusivity,tmax,tmsl,tlapse,heatproduction,friction
      double precision D_heat
      double precision xl,yl,zl,xlon1,xlat1,xlon2,xlat2,x,y,r,s
      double precision xlonobs,xlatobs,heightobs,ageheobs,dageheobs,ageftobs,dageftobs
      double precision wobs1,wobs2,wobs3,wobs4,a1,b1,c1,a2,b2,c2,a3,b3,c3,surf2
      double precision ageheZobs,dageheZobs,ageftZobs,dageftZobs,ageKarobs,dageKarobs,fmeano,dfmeano
      double precision ageBarobs,dageBarobs,ageMarobs,dageMarobs,ageHarobs,dageHarobs
      double precision ftdist(20)
      double precision time_incision_start, time_incision_stop, minZ, maxZ, tauH, Di, minW
      double precision Final_offset, Final_amp, Final_phase
      real *4 grainsize
      double precision thist(1024),temphist(1024),errortemphist(1024)
      ! Helium thermochronometry
      real *4 AUppm, AThppm,ASize ! Maxime
      real *4 D0_spec, Ea_spec
      double precision kinFTL_AHe, KinFTL_AFT !Maxime
      real *4 ZUppm, ZThppm,ZSize ! Maxime
      real *4 D0z, Eaz ! Maxime
      real *4 D0k, Eak ! Maxime
      real *4 D0b, Eab ! Maxime
      real *4 D0m, Eam ! Maxime
      real *4 D0h, Eah ! Maxime
      integer Alpha_Ejec_Flag, RDmodel,RDmodelz !Maxime
      ! 4He/3He theermochronometry 
      double precision theating(1024),released(1024),dreleased(1024),agereleased(1024),dagereleased(1024),duration(1024)
      double precision size43He,age43He,dage43He
      real Uppm43, Thppm43, rmr043, rhoST !Maxime 
      ! TL/OSL thermochronometry
      integer  OSL_Model, TL_Model, ESR_Model
      double precision doser,d0,radius,et,logs,b,logrho,nn,dnn,eu,sigmaet,imax,a_coef, GOK_a, GOK_b
      character run*5,fnme*300,obsfile*300,line*1024,c5*7,PecubeFnme*10,cproc*4
      ! General parameters
      character*100 :: sample
      character*50 :: sampleID

      real*4 range(2,*),param(*),xf
      integer nd,nd0,nc5
      logical vivi,xyz !VKP
      logical saveCr !Maxime
      logical is_unix !Maxime
      
      ! Check Os system
      call GetOperatingSystem(is_unix)

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
      nstep = p%ntime+sum(p%nstep)*2 ! Maxime
      tau = p%erosional_time_scale
      Headward = p%do_Headward
      tauH = p%tauH
      Di = p%depth_incision
      time_incision_start = p%time_incision_start
      time_incision_stop = p%time_incision_stop
      Final_offset = p%topo_offset ! vertical offset for synthetic topography (sinusoidal)
      Final_amp = p%topo_amp ! Relief amplification for final time step for synthetic topography (sinusoidal)
      Final_phase = p%topo_phase ! Phase shift at the end of simulation


! In the following, we ensure that all times are recorded in time topo. This is to prevent
! miscalculation of initial depth of samples and to ensure that samples end
! en up with their surface temperature. Some miscalculation arises because
! a tectonic time (i.e., time start or time end), can be within a time step 'dt'.
! It results that the teconic uplift rate applies for longer time that wanted, and thus
! the a greater initial depth for samples. To prevent such issue, we keep all times input
! by the user and interpolate the amplification and offset value for the time not set
! by time_topo. (Maxime)

! for each time step + 1
! timek is the time (in Myr)
! topomag is the topographic amplification factor at this step
! topooffset is the topographic offset
    allocate (timek_temp(p%ntime+1),topomag_temp(p%ntime+1),topooffset_temp(p%ntime+1),iout_temp(p%ntime+1),&
                phase_array_temp(p%ntime+1))
        
        ! Record input time topo        
        do istep=1,p%ntime+1
            timek_temp(istep) = p%time_topo(istep)
            topomag_temp(istep) = p%amplification(istep)
            topooffset_temp(istep) = p%offset(istep)
            iout_temp(istep) = p%output(istep)
            phase_array_temp(istep) = p%phase(istep)
        enddo

      allocate (timek_temp2(nstep+1),topomag_temp2(nstep+1),topooffset_temp2(nstep+1),iout_temp2(nstep+1),&
                    phase_array_temp2(nstep+1))

       !!!! 1. get all time topo !!!!
       ij = 0
       timek_temp2 (:) = 0.d0
       topomag_temp2 (:) = 0.d0
       topooffset_temp2 (:) = 0.d0
       iout_temp2 (:) = 0
       do istep=1,p%ntime+1
          ij = ij +1
          timek_temp2(istep) = p%time_topo(istep)
          iout_temp2(istep) = p%output(istep)
       enddo
       
       !!!!!! 2. get all times from tectonic scenario  !!!!!
       do k=1, p%nfault
            do istep=1,p%nstep(k)
                ! check if time_start is already in time topo
                res = 0
                call isinarray(timek_temp2,p%time_start(istep,k), nstep,res)
                if (res.eq.0) then 
                    ij = ij+1
                    timek_temp2(ij) = p%time_start(istep,k)
                    iout_temp2(ij) = 0
                endif
                ! check if time_end is already in time topo
                res = 0
                call isinarray(timek_temp2,p%time_end(istep,k), nstep,res)
                if (res.eq.0) then 
                    ij = ij+1
                    timek_temp2(ij) = p%time_end(istep,k)
                    iout_temp2(ij) = 0
                endif
            enddo
        enddo
        
        !!! 3. reduce array size to number of times recorded !!!!
        allocate (timek(ij),topomag(ij),topooffset(ij),iout(ij),&
                    phase_array(ij),index_array(ij),index_array2(p%ntime+1))
                    
        timek(:) = 0.d0
        timek(:) = timek_temp2(1:ij)
        topomag(:) = 0.d0
        topooffset(:) = 0.d0
        iout(:) = 0
        iout(:) = iout_temp2(1:ij)
        !!!!!  4. sort the array  !!!!!!
        index_array(:) = 0
        call sort_index(timek,index_array,ij,1) !ascendent order
        call sort_index(timek_temp,index_array2,p%ntime+1,1) !ascendent order
        iout(:) = iout(index_array)

        !!!! 5. interpolate topo values for times !!!!
        ! We assume a linear interpolation between time topo step
        call interpolate_1D(timek_temp,timek,p%ntime+1,topomag_temp(index_array2),topomag,ij)
        call interpolate_1D(timek_temp,timek,p%ntime+1,topooffset_temp(index_array2),topooffset,ij)
        call interpolate_1D(timek_temp,timek,p%ntime+1,phase_array_temp(index_array2),phase_array,ij)
        
        !!!! 6. sort back to descendant order
        call sort_index(timek,index_array,ij,0) !descendent order
        iout(:) = iout(index_array)
        topomag = topomag(index_array)
    topooffset = topooffset(index_array)
        ! do k=1, ij ! Make sure amplification is not negative
        !   if (topomag(k).lt.0.d0) then
        !       topomag(k) = 0.d0
        !   endif
      !   enddo
        ! print *, 'timek = ',timek
        ! print *, 'topomag = ',topomag
        ! print *, 'topooffset = ',topooffset
        ! print *, 'iout = ',iout

        
        nstep = ij-1
       deallocate (timek_temp,topomag_temp,topooffset_temp,iout_temp,phase_array_temp)
       deallocate (timek_temp2,topomag_temp2,topooffset_temp2,iout_temp2,phase_array_temp2)
       
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
      heatprod_flag = p%heatprod_flag
      D_heat = p%D_heat

      crustal_thickness=crustal_thickness*1.d3

! Parameters for AHe/AFT age computation - Maxime
      AUppm = p%AUPPM
      AThppm = p%AThPPM
      RDmodel = p%RDmodel
      Alpha_Ejec_Flag = p%Alpha_Ejec_Flag
      D0_spec = p%D0_AHe
      Ea_spec = p%Ea_AHe
      kinFTL_AFT = p%Kinetic_FTL_Parameter_value_AFT
    kinFTL_AHe = p%Kinetic_FTL_Parameter_value_AHe
      rhoST = p%rhoST
      D0z = p%D0_ZHe
      Eaz = p%Ea_ZHe
      ZUppm = p%ZUppm
      ZThppm = p%ZThppm
      ZSize = p%ZSize
      RDmodelz = p%RDmodelz
      D0k = p%D0_KAr
      Eak = p%Ea_KAr
      D0b = p%D0_BAr
      Eab = p%Ea_BAr
      D0m = p%D0_MAr
      Eam = p%Ea_MAr
      D0h = p%D0_HAr
      Eah = p%Ea_HAr
      OSL_Model = p%OSL_Model
      TL_Model = p%TL_Model
      ESR_Model = p%ESR_Model
      
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

! ! In output parameters - Sample-specific flag - Maxime B
      agecomput = 0
      saveCr = .FALSE.
      agecomput = p%age_computation
      saveCr = p%save_cooling_rates

! reads in topography

if (.not.vivi) then !VKP

      if (nx0.gt.0) then

        allocate (zNZ(nx0,ny0))
        if (fnme(1:nfnme).eq.'Nil' .or. p%topo_wavelength.gt.0.d0) then
          zNZ=0.d0
          
        elseif (fnme(1:nfnme).eq.'Topo30') then
           dx = 360.d0/43200
           dy = dx
           write (cproc,'(i4)') iproc
           if (iproc.lt.10) cproc(1:3)='000'
           if (iproc.lt.100) cproc(1:2)='00'
           if (iproc.lt.1000) cproc(1:1)='0'
           PecubeFnme = trim(fnme(1:nfnme))
           call ExtractDEM (xlon, xlat, nx0, ny0, PecubeFnme = run//'/data/'//trim(PecubeFnme))
           open (8,file=run//'/data/'//trim(PecubeFnme)//'.dat',status='old')
           read (8,*) zNZ
           close (8)
           if (is_unix) then
                call system ('rm '//run//'/data/'//trim(PecubeFnme)//'.dat')
            else !Assume it is Windows
                call system ('del /Q '//run//'/data/'//trim(PecubeFnme)//'.dat')
            endif
           
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
            !Ensure lat0, lon0, dx, dy are positive - Maxime
            dx = abs(p%dlon)
            dy = abs(p%dlat)
            xl = xlon+(nx0-1)*dx 
            yl = xlat+(ny0-1)*dy
            if (p%lon0.lt.0) xlon = p%lon0 + 360
          endif
        close (8)
        endif

      nx=(nx0-1)/nskip+1
      ny=(ny0-1)/nskip+1
      allocate (z(nx*ny),ztemp(nx*ny),Zincision(nx*ny),NewTopo(nx*ny))
      
      ! Degraded topography (nskip)
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

        allocate (z(-nx0),xz(-nx0),yz(-nx0),iconz(3,-ny0),Zincision(-nx0),NewTopo(-nx0))
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
      allocate (z(nx*ny),Zincision(nx*ny),NewTopo(nx*ny)) !VKP
      topomag=1.d0 !VKP
      topooffset=0.d0 !VKP

      xlon1=xlon
      xlat1=xlat
      xlon2=xlon+(nx-1)*dx*nskip
      xlat2=xlat+(ny-1)*dy*nskip

      else

      allocate (z(-nx0),xz(-nx0),yz(-nx0),iconz(3,-ny0), Zincision(-nx0), NewTopo(-nx0))
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

        xl=dx*(nx0-1)*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
        yl=dy*(ny0-1)*111.11
        zl=crustal_thickness/1.e3
        z=z/crustal_thickness*zl
        ztemp = z
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

      ilog=p%debug
      iterative=1
      interpol=1
      nsurf=nx*ny
      nelem=(nx-1)*(ny-1)
      nne=4
      ddx=xl/(nx-1)*1.d3
      ddy=yl/(ny-1)*1.d3
      allocate (x_array(nsurf)) ! array for x values (used in synthetic topography)
      if (nx0.lt.0) then
        nsurf=-nx0
        nelem=-ny0
        nne=3
        ddx=xl/sqrt(dble(nsurf))*1.d3
        ddy=yl/sqrt(dble(nsurf))*1.d3
      endif
    
      !Start to write in the file
      write (7,'(a)') run
      if (vivi) then !VKP
          write (7,*) nne,-nsurf,nz,nelem,zl,diffusivity,heatproduction,friction,heatprod_flag !VKP
      else !VKP
          write (7,*) nne,nsurf,nz,nelem,zl,diffusivity,heatproduction,friction,heatprod_flag
      endif !VKP
      write (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol,tprevious,ftlflag,mftflag,fltflag,D_heat
      write (7,*) isoflag,tau,rhoc,rhom
      write (7,*) nx,ny,nxiso,nyiso
      write (7,*) ddx,ddy,young,poisson,thickness*1.d3
      write (7,*) xlon1,xlon2,xlat1,xlat2

      if (nx0.gt.0) then
        ij = 0
        do j=1,ny
          do i=1,nx
          ij = ij +1
          x=xl*float(i-1)/float(nx-1)
          y=yl*float(j-1)/float(ny-1)
          x_array(ij) = x !Maxime
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
      
      
      ! !Write parameters for sample-specific prediction. Added by Maxime Bernard
      write (7,*) agecomput,saveCr!,diffmodel,alphaejec,alphadist,ngrains
      write (7,*) D0_spec,Ea_spec,RDmodel,Alpha_Ejec_Flag,D0z,Eaz,RDmodelz
    
     ! Initiate amplication matrix (Maxime)
      allocate (amplification(nx*ny))
      amplification(:) =  topomag(1)
      
    ! Topography evolution - loop time steps
      do istep=0,nstep
        write (7,*) timek(istep+1),iout(istep+1)
          if (vivi) then !VKP
            if (istep.lt.10) then !VKP
              write (c5(1:3),'(a2,i1)') '00',istep !VKP
              nc5=3 !VKP
            elseif (istep.lt.100) then !VKP
              write (c5(1:3),'(a1,i2)') '0',istep !VKP
              nc5=3 !VKP
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
         
        !Maxime - Check the reference elevation from which to compute topographic evolution
          if (istep.eq.0) then
            if (p%topo_ref.eq.0) then !from sea level
              topoRef = 0
            elseif (p%topo_ref.eq.1) then !from min elevation
              topoRef = minval(z(:))
            elseif (p%topo_ref.eq.2) then !from max elevation
              topoRef = maxval(z(:))
            elseif (p%topo_ref.eq.3) then !from mean elevation
              topoRef = sum(z(:))/(nx*ny)
            elseif (p%topo_ref.eq.4) then !from custom elevation
              topoRef = p%topo_ref_custom/1e3
            elseif (p%topo_ref.eq.5) then !from mean elevation but varying amplitude
              topoRef = sum(z(:))/(nx*ny)
            endif !end Maxime
            
            ! Get topography before incision (in case of headward propagation scenario)
            call get_initial_topo(topomag, topooffset,  maxval(timek)-timek, nstep, topoRef, z, nsurf,&
                                    time_incision_start,time_incision_stop, Zincision)
            minZ = minval(Zincision(:))
            maxZ = maxval(Zincision(:))
            NewTopo(:) = topoRef-(topomag(1)*(topoRef-z))+topooffset(1) ! Topo for headward propagation scenario
          endif
        
          if (topomag(istep+1).lt.0.) then ! Smoothed topography
          allocate(zsmooth(nsurf))
          call smooth_topo(z, zsmooth, nx, ny, int(-topomag(istep+1)), int(topooffset(istep+1)))
              write (7,*) (zsmooth(k), k=1,nsurf)
          ! print *, zsmooth
          deallocate(zsmooth)
      else
           if (p%topo_wavelength.gt.0.d0) then ! build synusoidal topography
                call Synthetic_topography(topomag,p%topo_ref, z, topooffset, istep+1, nstep, nsurf,&
                    x_array,p%topo_wavelength,Final_offset,phase_array,Final_amp,Final_phase,p%topo_ref_custom/1e3,nx,ny)
                    write (7,*) (z(k),k=1,nsurf)

             else ! Others topography
                if (p%topo_ref.lt.5.and.Headward.eq.0) then
                    write (7,*) (topoRef-(topomag(istep+1)*(topoRef-z(k)))+topooffset(istep+1),k=1,nsurf)
                    
                elseif (p%topo_ref.lt.5.and.Headward.eq.1) then
                    if (istep.gt.0) then 
                        call Headward_propagation(topomag,topoRef, NewTopo, topooffset, istep, nstep, amplification,&
                                                  nsurf,Di, time_incision_start,time_incision_stop, Zincision, z,&
                                                  maxval(timek(:))-timek(istep+1), maxval(timek(:))-timek(istep),&
                                                  minZ, maxZ, tauH)
                    endif
                                                
                    write (7,*) (NewTopo(k),k=1,nsurf)
                else
                    ! Compute new topo
                    do k=1,nsurf
                       if (z(k).lt.topoRef.and.topomag(istep+1).ne.0) then
                          ztemp(k) = topoRef-(1/topomag(istep+1)*(topoRef-z(k)))+topooffset(istep+1)
    !                       if (ztemp(k).gt.maxval(z(:))) then
    !                          ztemp(k) = maxval(z(:))
    !                       endif
                       else
                         ztemp(k) = topoRef-(topomag(istep+1)*(topoRef-z(k)))+topooffset(istep+1)
                       endif
                    enddo
                    write (7,*) (ztemp(k),k=1,nsurf)
                endif
             endif
     endif  
          
        
        !write (7,*) (z(k)*topomag(istep+1)+topooffset(istep+1),k=1,nsurf)
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
        write (7,*) nobs,0,0,0,0,0
      else
        ! Read the observed files
        call read_data_folder (run//'/data/'//obsfile(1:nobsfile), p%echo_input_file, p%lon0, p%lon0+(nx-1)*dx*nskip,&
        p%lat0, p%lat0+(ny-1)*dy*nskip, iproc, nd, nobs1, nobs2, nobs3, nobs4, nobs5, nobs6) !Maxime - changed xlon1,xlon2,xlat1,xlat2 into p%**
    
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
      open (112,file=run//'/data/'//obsfile(1:nobsfile)//cproc//'samples.txt',status='old') ! Maxime

      if (p%save_PTT_paths.ne.0) nobs1 = -nobs1
        write (7,*) nobs1,nobs2,nobs3,nobs4,nobs5,nobs6
      if (nobs1.lt.0) nobs1=-nobs1
      
      ! ! Print number of conventional ages observed
      ! if (nd.eq.0.and.iproc.eq.0) print*,'Number of Observations:',nobs1+nobs2+nobs3+nobs4+nobs5+nobs6

      do i=1,nobs1
        read (8,*) xlonobs,xlatobs,heightobs,ageheobs,dageheobs,ageftobs,dageftobs, &
                   ageheZobs,dageheZobs,ageftZobs,dageftZobs,ageKarobs,dageKarobs, &
                   ageBarobs,dageBarobs,ageMarobs,dageMarobs,ageHarobs,dageHarobs, &
                   ftdist,grainsize,ASize,AUppm,AThppm,kinFTL_AFT,kinFTL_AHe,ZUppm,ZThppm,ZSize,fmeano,&
                            dfmeano!Maxime (Uppm,Thppm,rmr0)

        if (nx0.gt.0) then
            !Ensure xlonobs is correct relative to xlon - Maxime
            xl = p%lon0+(p%nx-1)*dx 
            yl = p%lat0+(p%ny-1)*dy
            xlonobs = xlon + abs(xl - p%lon0) * ((xlonobs - p%lon0) / (xl - p%lon0))
            xlatobs = xlat + abs(yl - p%lat0) * ((xlatobs - p%lat0) / (yl - p%lat0))
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
        ! Read sample ID
        read (112,*) sampleID
        write (7,*) sampleID,xlonobs,xlatobs,ageheobs,dageheobs,ageftobs,dageftobs, &
                    ageheZobs,dageheZobs,ageftZobs,dageftZobs,ageKarobs,dageKarobs, &
                    ageBarobs,dageBarobs,ageMarobs,dageMarobs,ageHarobs,dageHarobs, &
                    ftdist,heightobs,fmeano,dfmeano,grainsize,ieobs,wobs1,wobs2,wobs3,wobs4,ASize,AUppm,AThppm,& 
                    kinFTL_AFT,kinFTL_AHe,ZUppm,ZThppm,ZSize !Maxime
        enddo
        close (112)
! cooling curves
        do i=1,nobs2
        read (8,*) xlonobs,xlatobs,heightobs,nhist, &
                   (thist(j),temphist(j),errortemphist(j),j=1,nhist)
        if (nx0.gt.0) then
        !Ensure xlonobs is correct relative to xlon - Maxime
        if (fnme(1:nfnme).ne.'Nil'.and.fnme(1:nfnme).ne.'Topo30') then
            xl = p%lon0+(p%nx-1)*dx 
            yl = p%lat0+(p%ny-1)*dy
            xlonobs = xlon + abs(xl - p%lon0) * ((xlonobs - p%lon0) / (xl - p%lon0))
            xlatobs = xlat + abs(yl - p%lat0) * ((xlatobs - p%lat0) / (yl - p%lat0))
        endif !Maxime
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
        read (8,*) sample,xlonobs,xlatobs,heightobs,size43He,age43He,dage43He,nheating, &
                   Uppm43,Thppm43,rmr043,(theating(j),duration(j),released(j),dreleased(j), &
                    agereleased(j),dagereleased(j),j=1,nheating) !Maxime added Uppm43 and Thppm43
        if (nx0.gt.0) then
        !Ensure xlonobs is correct relative to xlon - Maxime
        if (fnme(1:nfnme).ne.'Nil'.and.fnme(1:nfnme).ne.'Topo30') then
            xl = p%lon0+(p%nx-1)*dx 
            yl = p%lat0+(p%ny-1)*dy
            xlonobs = xlon + abs(xl - p%lon0) * ((xlonobs - p%lon0) / (xl - p%lon0))
            xlatobs = xlat + abs(yl - p%lat0) * ((xlatobs - p%lat0) / (yl - p%lat0))
        endif !Maxime
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
        write (7,*) sample,xlonobs,xlatobs,size43He,age43He,dage43He,nheating, &
                   Uppm43,Thppm43,rmr043,(theating(j),duration(j),released(j),dreleased(j), &
                    agereleased(j),dagereleased(j),j=1,nheating), &
                    heightobs,ieobs,wobs1,wobs2,wobs3,wobs4 !Maxime added Uppm43, Thppm43
        enddo
        
! TSL data
        do i=1,nobs4
        read (8,*) sampleID,xlonobs,xlatobs,heightobs,doser,d0,radius,et,logs,b,logrho,nn,dnn
        if (nx0.gt.0) then
        !Ensure xlonobs is correct relative to xlon - Maxime
        if (fnme(1:nfnme).ne.'Nil'.and.fnme(1:nfnme).ne.'Topo30') then
            xl = p%lon0+(p%nx-1)*dx 
            yl = p%lat0+(p%ny-1)*dy
            xlonobs = xlon + abs(xl - p%lon0) * ((xlonobs - p%lon0) / (xl - p%lon0))
            xlatobs = xlat + abs(yl - p%lat0) * ((xlatobs - p%lat0) / (yl - p%lat0))
        endif !Maxime
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
        write (7,*) sampleID,xlonobs,xlatobs, &
          doser,d0,radius,et,logs,b,logrho,nn,dnn, &
          heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
        enddo
        
! OSL data
        do i=1,nobs5
          read (8,*) sampleID,xlonobs,xlatobs,heightobs,doser,d0,et,logs,logrho,eu,nn,dnn,a_coef,imax
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
          write (7,*) sampleID,xlonobs,xlatobs, &
            doser,d0,et,logs,logrho,eu,nn,dnn,a_coef,imax, &
            heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
          enddo
        
! ESR data
          do i=1,nobs6
            read (8,*) sampleID,xlonobs,xlatobs,heightobs,doser,d0,logs,et,sigmaet,nn,dnn,GOK_a,GOK_b,imax
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
            write (7,*) sampleID,xlonobs,xlatobs, &
              doser,d0,logs,et,sigmaet,nn,dnn,GOK_a,GOK_b,imax, &
              heightobs,ieobs,wobs1,wobs2,wobs3,wobs4
            enddo
      close (8)

      if (is_unix) then
          call system ('rm '//run//'/data/'//obsfile(1:nobsfile)//cproc//'.txt')
          call system ('rm '//run//'/data/'//obsfile(1:nobsfile)//cproc//'samples.txt')
      else !Assume it is Windows
          call system ('del /Q '//run//'\data\'//obsfile(1:nobsfile)//cproc//'.txt')
          call system ('del /Q '//run//'\data\'//obsfile(1:nobsfile)//cproc//'samples.txt')
      endif
      endif
      write (7,*) ageflag
      
      deallocate (z,ztemp, Zincision, amplification, NewTopo)
      deallocate (x_array)
      
      if (nx0.lt.0) deallocate (xz,yz,iconz)
      deallocate (timek,topomag,topooffset)
      
      return
      end
