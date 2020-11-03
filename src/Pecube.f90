!--------------------------------------------------------------------------------------------

program Main

use Pecube
type (version) :: vers
type (parameters) :: p, p0
integer :: nd
real*4, dimension(1024) :: par
real*4, dimension(2,1024) :: range
real*4 :: misfit
character*5 run

character*9 arg
integer :: numarg

integer k
character*9 cparam(1024)


numarg=command_argument_count()
if (numarg.eq.0) then
write (*,*) 'You need to specify a run directory (i.e. RUN00 for example)'
stop 'End of run'
else
   call getarg (1,arg)
   if (arg.eq.'--version') then
      write (*,*) vers%str
      stop
   else
      run = arg
   endif
endif

! checks that the run directory exists, otherwise exits

open (7,file=run//'/input/Pecube.in',status='unknown',err=999)
goto 998
999 write (*,*) 'Cant open/create '//run//'/input/Pecube.in'
write (*,*) 'Check that the directory '//run//'/input exists'
stop
998 continue
close (7)

call system ('mkdir "'//run//'/output"')
call system ('mkdir "'//run//'/data"')
call system ('mkdir "'//run//'/NA"')
call system ('mkdir "'//run//'/LOG"')
call system ('rm '//run//'/output/*')

nd = 0
call read_input_file (run//'/input/Pecube.in', 0, p0, nd, range, par)

if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'
if (nd.eq.0) write (*,*) '-------------------------------Pecube--------------------------------------------'
if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'
if (nd.eq.0) write (*,*) 'version:',vers%str
if (p0%echo_input_file.gt.0) then
if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'
if (nd.eq.0) write (*,*) '-----------------------Echoing input parameters----------------------------------'
if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'
endif

nd = 0
par = 0.
range = 0.
if (p0%echo_input_file.eq.4) write (*,*) 'Complete list of model parameters'
call read_input_file (run//'/input/Pecube.in', p0%echo_input_file, p, nd, range, par)
if (p0%echo_input_file.eq.4) stop

if (p0%echo_input_file.gt.0) then
if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'
endif

  if (nd.ne.0) then
    do k=1,1024
    write (cparam(k),'(a5,i4)') 'Param',k
    if (k.lt.10) cparam(k)(6:8)='000'
    if (k.lt.100) cparam(k)(6:7)='00'
    if (k.lt.1000) cparam(k)(6:6)='0'
    enddo
  open (71,file=run//'/NA/NA_int_res.csv',status='unknown')
  write (71,'(1025a)') 'Misfit',(','//cparam(k),k=1,nd)
  close (71)
  open  (71,file=run//'/NA/na.in',status='unknown')
  write (71,'(a)') '#'
  write (71,'(a)') '#       Neighbourhood Algorithm input options file'
  write (71,'(a)') '#'
  write (71,'(a)') '0         : Algorithm type (NA or Uniform MC: 1=MC,0=NA)'
  write (71,'(i5,a)') p%maximum_number_of_iterations,'         : Maximum number of iterations'
  write (71,'(i5,a)') p%sample_size_for_first_iteration,'        : Sample size for first iteration'
  write (71,'(i5,a)') p%sample_size_for_all_other_iterations,'         : Sample size for all other iterations'
  write (71,'(i5,a)') p%number_of_cells_to_resample,'         : Number of cells to re-sample'
  write (71,'(a)') 'n,210728  : Use Quasi random number generator ? (y/n);random seed'
  write (71,'(a)') '0         : Type of initial sample (0=random;1=read in a NAD file)'
  write (71,'(a)') '1         : Output information level (0=silent,1=summary info,2=1+models)'
  write (71,'(a)') 'y         : Turn timing mode on ? (y/n)'
  write (71,'(a)') 'n         : Turn debug mode on ? (y/n)'
  close (71)
    if (p%save_ages_inversion.eq.1) then
      if (p%age_AHe_flag.ne.0.) then
      open (71,file=run//'/NA/AgeApatiteHelium.csv',status='unknown',position='rewind')
      write (71,'("Predicted Apatite Helium Ages")')
      close (71)
      endif
      if (p%age_AFT_flag.ne.0.) then
      open (71,file=run//'/NA/AgeApatiteFT.csv',status='unknown',position='rewind')
      write (71,'("Predicted Apatite FT Ages")')
      close (71)
      endif
      if (p%age_ZHe_flag.ne.0.) then
      open (71,file=run//'/NA/AgeZirconHelium.csv',status='unknown',position='rewind')
      write (71,'("Predicted Zircon Helium Ages")')
      close (71)
      endif
      if (p%age_ZFT_flag.ne.0.) then
      open (71,file=run//'/NA/AgeZirconFT.csv',status='unknown',position='rewind')
      write (71,'("Predicted Zircon FT Ages")')
      close (71)
      endif
      if (p%age_KAr_flag.ne.0.) then
      open (71,file=run//'/NA/AgeKSparArgon.csv',status='unknown',position='rewind')
      write (71,'("Predicted K-spar Argon Ages")')
      close (71)
      endif
      if (p%age_BAr_flag.ne.0.) then
      open (71,file=run//'/NA/AgeBiotiteArgon.csv',status='unknown',position='rewind')
      write (71,'("Predicted Biotite Argon Ages")')
      close (71)
      endif
      if (p%age_MAr_flag.ne.0.) then
      open (71,file=run//'/NA/AgeMuscoviteArgon.csv',status='unknown',position='rewind')
      write (71,'("Predicted Muscovite Argon Ages")')
      close (71)
      endif
      if (p%age_HAr_flag.ne.0.) then
      open (71,file=run//'/NA/AgeHornblendeArgon.csv',status='unknown',position='rewind')
      write (71,'("Predicted Hornblende Argon Ages")')
      close (71)
      endif
    endif
  endif
  if (nd.eq.0) then
  call forward (nd, par, misfit, run, 0)
  else
  call NA (nd, range, run)
  endif

end program Main

!--------------------------------------------------------------------------------------------

subroutine user_init (nd,range,scales)

implicit none

real*4 range(2,*)
real*4 scales(*)
integer nd

scales(1:nd)=-1.

end subroutine user_init

!--------------------------------------------------------------------------------------------

subroutine writemodels (nd,ntot,models,misfit,ns1,ns2,itmax,nh_max,nh,header,run)

implicit none

integer :: nd,ntot,ns1,ns2,itmax,nh_max,nh
real*4 :: models (nd,*)
real*4 :: misfit (ntot)
character*(*) :: header
character*5 run
character*9 cparam(1024)

integer :: i,k

  do k=1,1024
  write (cparam(k),'(a5,i4)') 'param',k
  if (k.lt.10) cparam(k)(6:8)='000'
  if (k.lt.100) cparam(k)(6:7)='00'
  if (k.lt.1000) cparam(k)(6:6)='0'
  enddo

open (87,file=run//'/NA/NA_results.csv',status='unknown')
write (87,'(1025a)') 'Misfit',(','//cparam(k),k=1,nd)
  do i=1,ntot
  write (87,'(g15.9,1024(a,g15.9))') misfit(i),(',',models(k,i),k=1,nd)
  enddo
close (87)

end subroutine writemodels

!--------------------------------------------------------------------------------------------

!subroutine forwardp (nd, par, misfit)

!use Pecube

!implicit none

!type (parameters) :: p
!integer nd
!real*4 par(nd), misfit
!real*4 range(2,1024)

!nd=0
!call read_input_file ('Pecube.in', 0, p, nd, range, par)

! here I should be able to run Pecube using a single line such as:
! call Pecube (p, nd, misfit)
! that would take the input parameters stored in p and return a misfit

!return
!end

!--------------------------------------------------------------------------------------------

      subroutine forward (nd,param,misfit,run, iproc)

      use Pecube

!      implicit real*8 (a-h,o-z)
      implicit none

      type (version) :: vers
      type (parameters) p
      type (faulttype),dimension(:),allocatable::fault,faultp

      double precision,dimension(:),allocatable :: x,y,z,xp,yp,zp
      double precision,dimension(:),allocatable :: u,up !by VKP
      double precision,dimension(:),allocatable :: t,tp,f
      integer,dimension(:,:),allocatable :: icon
      double precision,dimension(:,:,:),allocatable :: ael
      double precision,dimension(:,:),allocatable :: bel
      integer, dimension(:),allocatable :: kfix,ielsurf,ieo
      double precision,dimension(:),allocatable :: xsurf,ysurf,zsurf,zsurfp
      double precision,dimension(:),allocatable :: usurf,usurfp,tempsurf,tempsurfp !by VKP
      double precision,dimension(:),allocatable :: topoa,topob,rsurf,rebound
      integer,dimension(:,:),allocatable :: iconsurf,neighbour
      double precision,dimension(:,:),allocatable :: xdepth,ydepth,zdepth
      double precision,dimension(:,:),allocatable::ftdist,ftdisto
      double precision,dimension(:),allocatable :: xdeptho,ydeptho,zdeptho,zsurfo,zsurfpo
      double precision,dimension(:,:),allocatable :: tprev
      double precision,dimension(:),allocatable :: tprevo
      double precision, dimension(:,:),allocatable:: xexhumation,yexhumation,zexhumation
      double precision, dimension(:),allocatable:: xexhumationo,yexhumationo,zexhumationo
      double precision,dimension(:),allocatable::we1,we2,we3,we4,het,diff,zinit,xinit,yinit
      integer,dimension(:),allocatable :: proc,ice
      real*4,dimension(:),allocatable :: t2,t2p,t3,t4,t5
      real*4,dimension(:),allocatable::vxx,vyy,vzz,exhumation
      integer, dimension(:), allocatable :: nhist,nheating
      double precision, dimension(:,:), allocatable :: thist,temphist,errortemphist,temphistp
      double precision, dimension(:), allocatable :: size43He,age43He,dage43He
      double precision, dimension(:,:), allocatable :: theating,released,dreleased,agereleased,dagereleased,duration
      double precision, dimension(:,:), allocatable :: time43He,temperature43He
      double precision, dimension(:), allocatable :: releasedp,agereleasedp
      double precision, dimension(:,:), allocatable :: doser,d0,radius,et,logs,b,logrho,nn
      double precision, dimension(:,:), allocatable :: timeTSL,temperatureTSL
      double precision age43Hep,nnf,paramsTSL(8)

      double precision, dimension(:,:),allocatable::lsf,lsfn,lsfx,lsfy,lsfxn,lsfyn

      character run*5,cstep*3,iq*1

      real*4 times(1000)
      real*4,dimension(:),allocatable :: age1,age2,age3,age4,age5,age6,age7,age8,grainsize,grainobs
      real*4 t1

      integer nproc,iproc,nfault,npe,mpe,nsurf,nz,nelemsurf,nstep,ilog,iterative,interpol
      integer isoflag,nx,ny,nxiso,nyiso,niteradvec,istep,kstep
      integer i1,i2,i3,i4,iout,istatic,ntime,nnode,nelem,ie,iesurf,i,j,k,kk,ido,nelemloc,je
      integer irec,jrec,in,itime,ij,iu,ic,ieobs,inp,iobs,krec,ktime,niter,nobs,nobs1,nobs2,nobs3,nobs4 !VKP
      integer ftlflag,mftflag,fltflag
      double precision x1,x2,x3,x4,x5,dt,time,alpha,zsurfmax,fact,tsurf,zh,ftime,ftimep,tprevious,tfake
      double precision eps,zl,diffusivity,heatproduction,zsurf1,zsurf2
      double precision tempsurf1,tempsurf2 !by VKP
      double precision tmax,tmsl,tlapse,tau,rhoc,rhom,xstep,ystep,young,poisson,thickness
      double precision xlonmin,xlonmax,xlatmin,xlatmax,xxx,yyy,zzz,vx,vy,vz,dxxx,dyyy,dzzz
      double precision xmin,xmax,ymin,ymax,timesurf,timesurfp,tfinal,dtimesurf
      double precision xlonobs,xlatobs,tnow
      double precision ageheobs,dageheobs
      double precision ageftobs,dageftobs
      double precision ageheZobs,dageheZobs
      double precision ageftZobs,dageftZobs
      double precision agearKobs,dagearKobs
      double precision agearBobs,dagearBobs
      double precision agearMobs,dagearMobs
      double precision agearHobs,dagearHobs
      double precision tlobs(20),length(20)
      real*4,dimension(:),allocatable::ageft,agehe,ageheZ,ageftZ,agearK,agearB,agearM,agearH,hei,nNs
      real*4,dimension(:),allocatable::lonobs,latobs
      double precision wobs1,wobs2,wobs3,wobs4,x1f,x2f,y1f,y2f
      double precision def,dif,heightobs,friction
      double precision,dimension(:),allocatable::aheoreg,hheoreg,ahepreg
      double precision,dimension(:),allocatable::aftoreg,hftoreg,aftpreg
      double precision,dimension(:),allocatable::aheZoreg,hheZoreg,aheZpreg
      double precision,dimension(:),allocatable::aftZoreg,hftZoreg,aftZpreg
      double precision,dimension(:),allocatable::aarKoreg,harKoreg,aarKpreg
      double precision,dimension(:),allocatable::aarBoreg,harBoreg,aarBpreg
      double precision,dimension(:),allocatable::aarMoreg,harMoreg,aarMpreg
      double precision,dimension(:),allocatable::aarHoreg,harHoreg,aarHpreg
      double precision,dimension(:),allocatable::hpreg
      integer nhe,nft,nheZ,nftZ,narK,narB,narM,narH,ageflag(9),ageflagp(9)
      double precision aheom,hheom,heoi,heos
      double precision aftom,hftom,ftoi,ftos
      double precision aheZom,hheZom,heZoi,heZos
      double precision aftZom,hftZom,ftZoi,ftZos
      double precision aarKom,harKom,arKoi,arKos
      double precision aarBom,harBom,arBoi,arBos
      double precision aarMom,harMom,arMoi,arMos
      double precision aarHom,harHom,arHoi,arHos
      double precision ahepm,hhepm,hepi,heps
      double precision aftpm,hftpm,ftpi,ftps
      double precision aheZpm,hheZpm,heZpi,heZps
      double precision aftZpm,hftZpm,ftZpi,ftZps
      double precision aarKpm,harKpm,arKpi,arKps
      double precision aarBpm,harBpm,arBpi,arBps
      double precision aarMpm,harMpm,arMpi,arMps
      double precision aarHpm,harHpm,arHpi,arHps
      double precision volume_eroded(10000),surf_element,time_eroded(10000)
      integer nvolume
      double precision element_surface,volume_total
      external element_surface

      logical interpol_success,vivi,saveTt

      integer nd,nmisfit,nmisfitp,nmisfit1,nmisfit2,nmisfit3,nmisfit4,nmisfit1a,nmisfit3a
      real*4 param(nd)
      real*4 misfit,misfitp
      real*4 range(2,1024)
      real*4 misfit1, misfit2, misfit3, misfit4,misfit1a,misfit3a

      external hinterpolate,thermal_diffusivity,heat_production,TLModel
      double precision hinterpolate,zmax,zlsf,thermal_diffusivity,heat_production,TLModel,dummy

      call cpu_time (times(1))
      nproc=1

      eps=tiny(eps)

        do i=1,17
        length(i)=i
        enddo

! Pecube is a Finite Element solver of the 3D, transient heat transfer
! equation that allows for conduction, vertical advection and production of
! heat. Pecube also allows for the surface geometry to vary with time.

! This version is an improvement (we hope) of the original version
! available by download from the author's website in that it reads an
! input (topography) file as well as a set of cooling ages (He and FT in
! apatite of known location and calculates synthetic ages at the same
! locations for comparison.

! A flexural isostatic model has also been included to
! calculate the effect of isostasy in amplifying the exhumation caused
! by erosional unloading.

! Two dimensional faults have also been included; they can be orientated
! in any direction; more than one can be active at any given time and
! they advect each other (this only works for simple scenario)

! To understand how this version works, read the information in the file
! README

! comment this line out if you want Pecube to read your own input file
! (named Pecube.in)

      nd=0
      p%run_name = run
      call read_input_file (run//'/input/Pecube.in', 0, p, nd, range, param)
      if (nd.ne.0 .and. iproc.eq.0) then
      write (*,*) '---------------------------------------------------------------------------------'
      write (*,*) '------------------------Pecube inversion mode------------------------------------'
      write (*,*) '---------------------------------------------------------------------------------'
      write (*,*) 'Number of model parameters:',nd
      endif
      nfault = p%nfault
      if (nfault.eq.0) nfault=1

      if (nfault.gt.0) allocate (fault(nfault),faultp(nfault))
      call create_pecube_in (fault,nfault,nd,range,param,run,iproc)
      call cpu_time (times(2))

! opens input files

      rewind (7)

! read in general information

! first line:
! run: 5 character string that will determine the name of the folder where the input file
!      will be copied to and where the output files (Pecube.ou and Pecube.ptt) will be stored.

      read (7,'(a)') run

! second line
! npe: number of nodes per surface (2D) elements (3 = triangular elements - 4 = rectangular elements)
! nsurf: number of nodes defining the surface geometry
! nz: number of nodes in the vertical direction
! nelemsurf: number of surface 2D elements
! zl: thickness of crustal layer (= depth at which t is fixed) (in km)
! diffusivity: thermal diffusivity (in km^2/Myr)
! heatproduction: heat production (in degC/Myr)

      read (7,*) npe,nsurf,nz,nelemsurf,zl,diffusivity,heatproduction,friction
      vivi=.FALSE.
        if (nsurf.lt.0) then
        nsurf=-nsurf
        vivi=.TRUE.
        endif
      mpe=2*npe

! second line:
! tmax: basal temperature (at z=-zl) (in degC)
! tmsl: temperature at mean sea level (z=0) (in degC)
! tlapse: lapse rate (indegC/km)
! nstep: number of stages (the surface topo and exhumatoin rate can be specified at each time stage)
! ilog: redundant (dont use)
! iterative: iterative solver method (1 = Gauss Siedel with overrelaxation - 2 = Conjugate gradient)
! interpol: interpolation used (1 = linear - 2 = quadratic, not recommended)

      read (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol,tprevious,ftlflag,mftflag,fltflag
      read (7,*) isoflag,tau,rhoc,rhom
      read (7,*) nx,ny,nxiso,nyiso
      read (7,*) xstep,ystep,young,poisson,thickness

      allocate (iconsurf(npe,nelemsurf),ielsurf(nsurf))
      allocate (neighbour(npe,nelemsurf))
      allocate (xsurf(nsurf),ysurf(nsurf),zsurf(nsurf),zsurfp(nsurf))
      allocate (usurf(nsurf),usurfp(nsurf),tempsurf(nsurf),tempsurfp(nsurf)) !VKP
      allocate (topoa(nsurf),topob(nsurf),rsurf(nsurf))
      allocate (tprev(nsurf,nstep))
      allocate (xdepth(nsurf,nstep),ydepth(nsurf,nstep),zdepth(nsurf,nstep))
      allocate (xexhumation(nsurf,nstep),yexhumation(nsurf,nstep),zexhumation(nsurf,nstep))
      allocate (lsf(nsurf,nz),lsfx(nsurf,nz),lsfy(nsurf,nz))

      tprev=0.d0
      usurf=0.d0 !VKP
      usurfp=0.d0 !VKP
      tempsurf=0.d0 !VKP
      tempsurfp=0.d0 !VKP

      if (nd.eq.0.and.ilog.eq.1) open (9,file=run//'/LOG/Pecube.log',status='unknown')
      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'This is Pecube version ',vers%str
      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Reading Geometry'

! read in nodal geometry

      read (7,*) xlonmin,xlonmax,xlatmin,xlatmax

! xsurf, ysurf: x- and y-locations of the surface nodes
! iconsurf: connectivity matrix describing the 2D surface elements

      read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
      read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)

      zmax=0.d0
          do kstep=0,nstep
          read (7,*) timesurfp,iout
          read (7,*) (zsurfp(i),i=1,nsurf)
          zmax=max(maxval(zsurfp),zmax)
          if (vivi) read (7,*) (usurfp(i),i=1,nsurf)  !VKP
          if (vivi) read (7,*) (tempsurfp(i),i=1,nsurf)  !VKP
          enddo

      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Advecting rocks'

      xmin=minval(xsurf)
      xmax=maxval(xsurf)
      ymin=minval(ysurf)
      ymax=maxval(ysurf)

        do k=1,nsurf
        lsfx(k,:)=xsurf(k)
        lsfy(k,:)=ysurf(k)
        enddo

        do k=1,nz
        lsf(:,k)=zl-(zl+zmax)*float(k-1)/(nz-1)
        enddo

      nvolume=0

! go through the input file a first time to determine the depth of the
! points that will end up at the surface at the end of the model run

if (nd.eq.0) write (*,*) '-------------------Backward advecting surface nodes------------------------------'
if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'

      niteradvec=3
        do istep=1,nstep
        xexhumation(1:nsurf,istep)=xsurf
        yexhumation(1:nsurf,istep)=ysurf
        zexhumation(1:nsurf,istep)=zl
        enddo
        do istep=nstep,1,-1
        if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Time interval ',istep
        timesurfp=0.
        rewind (7)
        read (7,*)
        read (7,*) i1,i2,i3,i4,x1,x2,x3,x4
        read (7,*) x1,x2,x3,i1,i2,i3,i4,x1,i2,i3,i4
        read (7,*) i1,x1,x2,x3
        read (7,*) i1,i2,i3,i4
        read (7,*) x1,x2,x3,x4,x5
        read (7,*) x1,x2,x3,x4
        read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
        read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)
          do kstep=0,istep-1
          read (7,*) timesurfp,iout
          read (7,*) (zsurfp(i),i=1,nsurf)
          if (vivi) read (7,*) (usurfp(i),i=1,nsurf)  !VKP
          if (vivi) read (7,*) (tempsurfp(i),i=1,nsurf)  !VKP
          enddo
        read (7,*) timesurf,iout
        read (7,*) (zsurf(i),i=1,nsurf)
        if (vivi) read (7,*) (usurf(i),i=1,nsurf) !VKP
        if (vivi) read (7,*) (tempsurf(i),i=1,nsurf) !VKP
          if (istep.eq.nstep) then
          read (7,*) nobs1,nobs2,nobs3,nobs4
          saveTt=.FALSE.
            if (nobs1.lt.0.and.nd.eq.0) then
            saveTt=.TRUE.
            endif
          if (nobs1.lt.0) nobs1=-nobs1
          if (nobs2.lt.0) nobs2=-nobs2
          if (nobs3.lt.0) nobs3=-nobs3
          if (nobs4.lt.0) nobs4=-nobs4
          nobs=nobs1+nobs2+nobs3+nobs4
          allocate (xdeptho(nobs),ydeptho(nobs),zdeptho(nobs))
          allocate (zsurfo(nobs),zsurfpo(nobs))
          allocate (xexhumationo(nobs),yexhumationo(nobs),zexhumationo(nobs),grainobs(nobs))
          allocate (tprevo(nobs),ieo(nobs),we1(nobs),we2(nobs),we3(nobs),we4(nobs))
          allocate (nhist(nobs),thist(1024,nobs),temphist(1024,nobs),errortemphist(1024,nobs))
          allocate (nheating(nobs),theating(1024,nobs),released(1024,nobs),dreleased(1024,nobs),duration(1024,nobs))
          allocate (agereleased(1024,nobs),dagereleased(1024,nobs))
          allocate (size43He(nobs),age43He(nobs),dage43He(nobs))
          allocate (doser(1024,nobs),d0(1024,nobs),radius(1024,nobs),et(1024,nobs),logs(1024,nobs),b(1024,nobs))
          allocate (logrho(1024,nobs),nn(1024,nobs))
            do i=1,nobs1
            read (7,*) xlonobs,xlatobs,ageheobs,dageheobs,ageftobs,dageftobs, & ! By Xav
                       ageheZobs,dageheZobs,ageftZobs,dageftZobs,agearKobs,dagearKobs, &
                       agearBobs,dagearBobs,agearMobs,dagearMobs,agearHobs,dagearHobs, &
                       tlobs,heightobs,grainobs(i),ieo(i),we1(i),we2(i),we3(i),we4(i)   ! By Xav
            xexhumationo(i)=we1(i)*xsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*xsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*xsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*xsurf(iconsurf(npe,ieo(i)))
            yexhumationo(i)=we1(i)*ysurf(iconsurf(1,ieo(i))) &
                           +we2(i)*ysurf(iconsurf(2,ieo(i))) &
                           +we3(i)*ysurf(iconsurf(3,ieo(i))) &
                           +we4(i)*ysurf(iconsurf(npe,ieo(i)))
            zexhumationo(i)=we1(i)*zsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*zsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*zsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*zsurf(iconsurf(npe,ieo(i)))+zl+min(0.d0,heightobs/1.e3)
            enddo
            do i=nobs1+1,nobs1+nobs2
            read (7,*) xlonobs,xlatobs,nhist(i),(thist(j,i),temphist(j,i),errortemphist(j,i),j=1,nhist(i)), &
                       heightobs,ieo(i),we1(i),we2(i),we3(i),we4(i)   ! By Xav
            xexhumationo(i)=we1(i)*xsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*xsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*xsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*xsurf(iconsurf(npe,ieo(i)))
            yexhumationo(i)=we1(i)*ysurf(iconsurf(1,ieo(i))) &
                           +we2(i)*ysurf(iconsurf(2,ieo(i))) &
                           +we3(i)*ysurf(iconsurf(3,ieo(i))) &
                           +we4(i)*ysurf(iconsurf(npe,ieo(i)))
            zexhumationo(i)=we1(i)*zsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*zsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*zsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*zsurf(iconsurf(npe,ieo(i)))+zl+min(0.d0,heightobs/1.e3)
            enddo
            do i=nobs1+nobs2+1,nobs1+nobs2+nobs3
            read (7,*) xlonobs,xlatobs,size43He(i),age43He(i),dage43He(i),nheating(i), &
                      (theating(j,i),duration(j,i),released(j,i),dreleased(j,i), &
                       agereleased(j,i),dagereleased(j,i),j=1,nheating(i)), &
                       heightobs,ieo(i),we1(i),we2(i),we3(i),we4(i)   ! By Xav
            xexhumationo(i)=we1(i)*xsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*xsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*xsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*xsurf(iconsurf(npe,ieo(i)))
            yexhumationo(i)=we1(i)*ysurf(iconsurf(1,ieo(i))) &
                           +we2(i)*ysurf(iconsurf(2,ieo(i))) &
                           +we3(i)*ysurf(iconsurf(3,ieo(i))) &
                           +we4(i)*ysurf(iconsurf(npe,ieo(i)))
            zexhumationo(i)=we1(i)*zsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*zsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*zsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*zsurf(iconsurf(npe,ieo(i)))+zl+min(0.d0,heightobs/1.e3)
            enddo
            do i=nobs1+nobs2+nobs3+1,nobs1+nobs2+nobs3+nobs4
            read (7,*) xlonobs,xlatobs,nheating(i), &
                      (doser(j,i),d0(j,i),radius(j,i),et(j,i), &
                      logs(j,i),b(j,i),logrho(j,i),nn(j,i),j=1,nheating(i)), &
                      heightobs,ieo(i),we1(i),we2(i),we3(i),we4(i)   ! By Xav
            xexhumationo(i)=we1(i)*xsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*xsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*xsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*xsurf(iconsurf(npe,ieo(i)))
            yexhumationo(i)=we1(i)*ysurf(iconsurf(1,ieo(i))) &
                           +we2(i)*ysurf(iconsurf(2,ieo(i))) &
                           +we3(i)*ysurf(iconsurf(3,ieo(i))) &
                           +we4(i)*ysurf(iconsurf(npe,ieo(i)))
            zexhumationo(i)=we1(i)*zsurf(iconsurf(1,ieo(i))) &
                           +we2(i)*zsurf(iconsurf(2,ieo(i))) &
                           +we3(i)*zsurf(iconsurf(3,ieo(i))) &
                           +we4(i)*zsurf(iconsurf(npe,ieo(i)))+zl+min(0.d0,heightobs/1.e3)
            enddo
          read (7,*) ageflag
          endif
! isostatic rebound
        topoa=zsurfp
        topob=zsurf
          if (isoflag.eq.1) then
          call isostatic_rebound (topoa,topob,xsurf,ysurf,nsurf,rsurf, &
                                  rhoc,rhom,nx,ny,nxiso,nyiso, &
                                  xstep,ystep,young,poisson,thickness)
          else
          rsurf=0.
          endif
        zexhumation(1:nsurf,istep)=zexhumation(1:nsurf,istep)+zsurf
        if (istep.eq.nstep) tfinal=timesurf
        call find_dt (zl,diffusivity,nsurf,zsurf,zsurfp, &
                      nz,timesurf,timesurfp,istep,eps,ilog,nd, &
                      dt,ntime,istatic,fault,nfault,usurf)
        if (nd.eq.0) write (*,*) ntime,'time steps required'
          do ktime=1,ntime
          time=timesurf-dt*ktime
!          time=max(time,1.d-6)
            do kstep=istep,nstep
              do i=1,nsurf
              xxx=xexhumation(i,kstep)
              yyy=yexhumation(i,kstep)
              zzz=zexhumation(i,kstep)
                do k=1,niteradvec
                call find_velo (xxx,yyy,zzz,time,-dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,.true.)
                dxxx=xexhumation(i,kstep)-dt*vx/2.-xxx
                dyyy=yexhumation(i,kstep)-dt*vy/2.-yyy
                dzzz=zexhumation(i,kstep)-dt*vz/2.-zzz
                xxx=xxx+dxxx
                yyy=yyy+dyyy
                zzz=zzz+dzzz
                enddo
              xexhumation(i,kstep)=xexhumation(i,kstep)-dt*vx
              yexhumation(i,kstep)=yexhumation(i,kstep)-dt*vy
              zexhumation(i,kstep)=zexhumation(i,kstep)-dt*vz-rsurf(i)/ntime-dt*usurf(i) !VKP
              enddo
            enddo

! lines to compute volume eroded
        nvolume=nvolume+1
        if (nvolume.gt.10000) stop 'nvolume > 10000'
        volume_eroded(nvolume)=0.d0
          do i=1,nsurf
          surf_element=element_surface(npe,iconsurf,xsurf,ysurf,nsurf)
          dtimesurf=timesurf-timesurfp
          if (dtimesurf.gt.0.d0) volume_eroded=volume_eroded-(zsurf(i)-zsurfp(i))/dtimesurf*surf_element
          call find_velo (xsurf(i),ysurf(i),zl+zsurf(i),time,dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,.true.)
          volume_eroded(nvolume)=volume_eroded(nvolume)+(vz+rsurf(i)/dtimesurf+usurf(i))*surf_element
          enddo
        time_eroded(nvolume)=tfinal-time
! end lines to compute volume eroded

              do i=1,nobs
              xxx=xexhumationo(i)
              yyy=yexhumationo(i)
              zzz=zexhumationo(i)
                do k=1,niteradvec
                call find_velo (xxx,yyy,zzz,time,-dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,.true.)
                dxxx=xexhumationo(i)-dt*vx/2.-xxx
                dyyy=yexhumationo(i)-dt*vy/2.-yyy
                dzzz=zexhumationo(i)-dt*vz/2.-zzz
                xxx=xxx+dxxx
                yyy=yyy+dyyy
                zzz=zzz+dzzz
                enddo
              xexhumationo(i)=xexhumationo(i)-dt*vx
              yexhumationo(i)=yexhumationo(i)-dt*vy
              zexhumationo(i)=zexhumationo(i)-dt*vz-we1(i)*(rsurf(iconsurf(1,ieo(i)))/ntime+dt*usurf(iconsurf(1,ieo(i)))) &
                                                   -we2(i)*(rsurf(iconsurf(2,ieo(i)))/ntime+dt*usurf(iconsurf(2,ieo(i)))) &
                                                   -we3(i)*(rsurf(iconsurf(3,ieo(i)))/ntime+dt*usurf(iconsurf(3,ieo(i)))) &
                                                   -we4(i)*(rsurf(iconsurf(npe,ieo(i)))/ntime+dt*usurf(iconsurf(npe,ieo(i))))
              enddo
          if (fltflag.eq.1) call move_fault (fault,faultp,nfault,-dt,time,xmin,xmax,ymin,ymax)
          enddo
        zsurfp=zsurf
        timesurfp=timesurf
        enddo

        if (p%save_eroded_volume.eq.1) then
        open (43,file=run//'/output/VolumeEroded.csv',status='unknown')
        write (43,'(a)') 'Interval_start,Interval_End,Flux,Accumulated_Volume'
        volume_total=0.d0
          do i=nvolume,2,-1
          volume_total=volume_total+volume_eroded(i)*(time_eroded(i)-time_eroded(i-1))
          write (43,'(g12.6,",",g12.6,",",g12.6,",",g12.6)') time_eroded(i),time_eroded(i-1),volume_eroded(i),volume_total
          enddo
        volume_total=volume_total+volume_eroded(1)*time_eroded(1)
        write (43,'(g12.6,",",g12.6,",",g12.6,",",g12.6)') time_eroded(1),0.d0,volume_eroded(1),volume_total
        close (43)
        endif

! depth is the initial depth of the rocks that end up at the location of the
! surface nodes at the end of the model experiment

      xdepth=xexhumation
      ydepth=yexhumation
      zdepth=zexhumation
      xdeptho=xexhumationo
      ydeptho=yexhumationo
      zdeptho=zexhumationo

if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'

! reset the input file to its proper position

      if (nd.eq.0) write (*,*) ''

      rewind (7)
      read (7,'(a)') run
      read (7,*) npe,nsurf,nz,nelemsurf,zl,diffusivity,heatproduction,friction
      nsurf=iabs(nsurf)
      read (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol,tprevious,ftlflag,mftflag,fltflag
      nstep=abs(nstep)
      read (7,*) isoflag,tau,rhoc,rhom
      read (7,*) nx,ny,nxiso,nyiso
      read (7,*) xstep,ystep,young,poisson,thickness
      read (7,*) xlonmin,xlonmax,xlatmin,xlatmax
      read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
      read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)

      nnode=nsurf*nz
      nelem=nelemsurf*(nz-1)
      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'nnode/nelem= ',nnode,nelem

! opens output files

! Pecube.out contains the temperature field at the end of each stage

      if (iproc.eq.0.and.nd.eq.0) open (8,file=run//'/output/Pecube.out',status='unknown',access='direct', &
            recl=4*(4+12*nnode+mpe*nelem))

! Ages.out contains the ages at the end of each stage

      if (iproc.eq.0.and.nd.eq.0) open (11,file=run//'/output/Ages.out',status='unknown',access='direct', &
            recl=4*(4+14*nsurf+npe*nelemsurf))

! Pecube.ptt contains the depth-temperture-paths of all surface nodes

    if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Start of time stepping'

      do istep=1,nstep
      open (100+istep,status='scratch',access='direct', &
            recl=4*(1+nsurf*4))
      enddo

    open (500,status='scratch',access='direct', &
          recl=4*(1+nobs*4))

      allocate (x(nnode),y(nnode),z(nnode),t(nnode))
      allocate (u(nnode)) !VKP
      allocate (xp(nnode),yp(nnode),zp(nnode),tp(nnode))
      allocate (up(nnode)) !VKP
      allocate (icon(mpe,nelem))
      allocate (kfix(nnode))
      allocate (f(nnode),rebound(nnode))

      u=0.d0
      up=0.d0

! build 3D element connectivity

      ie=0
        do iesurf=1,nelemsurf
          do k=1,(nz-1)
          ie=ie+1
            do kk=1,npe
            icon(kk,ie)=(iconsurf(kk,iesurf)-1)*nz+k
            icon(kk+npe,ie)=icon(kk,ie)+1
            enddo
          enddo
        enddo

      if (ie.ne.nelem) then
      stop 'nelem mismatch'
      endif

! finds processor topology

      allocate (proc(nnode))
      proc=0
      nelemloc = nelem

      allocate (ael(mpe,mpe,nelemloc),bel(mpe,nelemloc),ice(nelemloc))

      je=0
        do ie=1,nelem
        je=je+1
        ice(je)=ie
        enddo

      if (nd.eq.0.and.ilog.eq.1) then
      write (9,*) 'icon'
        do ie=1,nelem
        write (9,*) (icon(k,ie),k=1,mpe)
        enddo
      endif

! finds neighbour conectivity matrix

      call find_neighbours (iconsurf,neighbour,npe,nelemsurf,nsurf)
      ielsurf=1

! initialize global parameters
! alpha is the time integration parameter

      alpha=0.5
      time=0.
      timesurf=0.
      irec=0
      jrec=0

if (nd.eq.0) write (*,*) '-----------------------------Time stepping---------------------------------------'
if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'

! begining of surface stepping (stageing)

      do istep=0,nstep

      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'istep= ',istep

      if (nd.eq.0) write (6,*) 'Interval : ',istep,' of ',nstep

      if (istep.ne.0) zsurfp=zsurf
      if (istep.ne.0) usurfp=usurf !VKP
      if (istep.ne.0) tempsurfp=tempsurf !VKP
      timesurfp=timesurf

! read the step information

      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Reading step information'

! timesurf is the time since the begining of the experiment
! Peclet is the exhumation velocity since the end of the last stage
! iout indicates if the T information is to be saved at the end of the stage

      read (7,*) timesurf,iout
        if (istep.eq.0) then
          if (timesurf.gt.eps.and.iproc.eq.0) then
          write (6,*) 'timesurf= ',timesurf
          write (6,*) 'first topography record must be at time zero ...'
          stop
          endif
        endif
      read (7,*) (zsurf(i),i=1,nsurf)
      if (vivi) read (7,*) (usurf(i),i=1,nsurf) !VKP
      if (vivi) read (7,*) (tempsurf(i),i=1,nsurf) !VKP

! initial (or zeroth) step

        if (istep.eq.0) then
        zsurfp=zsurf
        usurfp=usurf !VKP
        tempsurfp=tempsurf !VKP

        zsurfmax=maxval(zsurf)
        in=0
          do i=1,nsurf
            do k=1,nz
            in=in+1
            fact=float(k-1)/float(nz-1)
            if (interpol.eq.2.and.fact.gt.eps) fact=sqrt(fact)
!
            zh=zsurf(i)+zl
            zp(in)=zh*fact
! Kendra
!!!            zh=zl
!!!            zp(in)=zh*fact
!!!              if (k.eq.nz) zp(in)=zl+zsurf(i)
!
            xp(in)=xsurf(i)
            x(in)=xp(in)
            yp(in)=ysurf(i)
            y(in)=yp(in)
            up(in)=usurf(i)
            u(in)=up(in)
            kfix(in)=0
            if (k.eq.1 .or. k.eq.nz) kfix(in)=1
            enddo
          enddo

! calculates initial temperature

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Calculatin ginitial temperature'

          do i=1,nsurf
          tsurf=-zsurf(i)*tlapse+tmsl
          if (vivi) tsurf = tempsurf(i)
            do k=1,nz
            in=(i-1)*nz+k
            zh=zsurf(i)+zl
            tp(in)=tsurf+(tmax-tsurf)*(zh-zp(in))/zh
            enddo
          enddo

        t=tp

          if (nd.eq.0.and.ilog.eq.1) then
          write (9,*) 'Nodal geometry'
            do i=1,nnode
            write (9,*) i,xp(i),yp(i),zp(i),tp(i),kfix(i)
            enddo
          endif

        endif

      call find_dt (zl,diffusivity,nsurf,zsurf,zsurfp, &
                    nz,timesurf,timesurfp,istep,eps,ilog,nd, &
                    dt,ntime,istatic,fault,nfault,usurf)

! beginning of time stepping

      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Begining of time stepping'

        do itime=1,ntime
        call cpu_time (times(3))
        time=time+dt
        ftime=float(itime)/float(ntime)
        ftimep=float(itime-1)/float(ntime)
! Mod made by Jean on 15/3/2010
! now one should use tau=0 to indicate a linear change in topography
! a positive value for a change late in the time step
! a negative value for a change early in the time step (most meaningful geomorphologically)
          if (tau.ne.0.d0) then
          ftime=(1.-exp(-ftime*tau/(timesurf-timesurfp)))/(1.-exp(-tau/(timesurf-timesurfp)))
          ftimep=(1.-exp(-ftimep*tau/(timesurf-timesurfp)))/(1.-exp(-tau/(timesurf-timesurfp)))
      endif

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'itime,time= ',itime,time

! build new node geometry

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Build new node geometry'

        z=zp

          do i=1,nsurf
          in=i*nz
          inp=i*nz-1
          zsurf1=zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime+zl
          zsurf2=zsurfp(i)+(zsurf(i)-zsurfp(i))*ftimep+zl
          topoa(i)=zsurf2
          topob(i)=zsurf1
          z(in)=z(in)+zsurf1-zsurf2
          z(inp)=z(inp)+(zsurf1-zsurf2)/2.
! Mod made by Jean on 23/3/2010
! This insures that the surface temperature follows the lapse rate at every time step
! not just at t=0
          tp(in)=-(zsurfp(i)+(zsurf(i)-zsurfp(i))*ftimep)*tlapse+tmsl
          t(in)=-(zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)*tlapse+tmsl
            if (vivi) then !VKP
            tempsurf1 = tempsurfp(i) + (tempsurf(i)-tempsurfp(i))*ftime !VKP
            tempsurf2 = tempsurfp(i) + (tempsurf(i)-tempsurfp(i))*ftimep !VKP
!            t(in) = t(in) + tempsurf1 - tempsurf2 !VKP
            tp(in)=tempsurf2
            t(in)=tempsurf1
            endif !VKP
          enddo

       in=0
          do i=1,nsurf !VKP
            do k=1,nz !VKP
                in=in+1 !VKP
                u(in)=usurf(i) !VKP
            enddo !VKP
          enddo !VKP

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Isostatic rebound'

          if (isoflag.eq.1) then
          call isostatic_rebound (topoa,topob,xsurf,ysurf,nsurf,rsurf, &
                                  rhoc,rhom,nx,ny,nxiso,nyiso, &
                                  xstep,ystep,young,poisson,thickness)
          else
          rsurf=0.
          endif

          do i=1,nsurf
            do k=1,nz
            rebound(k+(i-1)*nz)=rsurf(i)
            enddo
          enddo

        if (in.ne.nnode) then
        stop 'nnode mismatch'
        endif

! build local FE matrices

      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Building matrix'

      call cpu_time (times(7))
      je=0
          do je=1,nelemloc
          ie=ice(je)
          call make_matrix (mpe,ael(1,1,je),bel(1,je),icon(1,ie), &
                            x,y,z,u,xp,yp,zp,up, & !VKP
                            kfix,diffusivity,heatproduction, &
                            x1f,y1f,x2f,y2f,def,dif, &
                            alpha,dt,time,tp,nnode,istatic, &
                            fault,nfault,rebound, &
                            xmin,xmax,ymin,ymax,0.d0,zl+zmax,friction, &
                            lsf,lsfx,lsfy,nx,ny,nz,tfinal)
          enddo
      call cpu_time (times(8))

! build global RHS vector

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Building RHS vector'

        f=0.d0
          do je=1,nelemloc
          ie=ice(je)
            do k=1,mpe
            ic=icon(k,ie)
            f(ic)=f(ic)+bel(k,je)
            enddo
          enddo

! solve global FE equations


      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Solving matrix'

      call cpu_time (times(5))
        if (iterative.eq.1) then
        call solve_iterative (mpe,1,ael,f,t,kfix,icon, &
                              nnode,nelemloc,niter,proc,ice)
        if (nd.eq.0.and.ilog.eq.1) write (9,*) niter,' iterations'
        elseif (iterative.eq.2) then
        stop 'solution strategy not implemented'
        else
        stop 'solution strategy not implemented'
        endif
      call cpu_time (times(6))

      if (iproc.eq.0.and.nd.eq.0) call screen_counter (itime,ntime,niter)

! stretch old grid

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Stretch old grid'

          do i=1,nsurf
            do k=1,nz
            in=(i-1)*nz+k
            fact=float(k-1)/float(nz-1)
            if (interpol.eq.2.and.fact.gt.eps) fact=sqrt(fact)
            zsurf1=zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime
!
            zh=zsurf1+zl
            zp(in)=zh*fact
! Kendra
!!!            zh=zl
!!!            zp(in)=zh*fact
!!!              if (k.eq.nz) zp(in)=zl+zsurf1
!

            enddo
          enddo

! update fault geometry

      if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Update fault geometry'

      if (fltflag.eq.1) call move_fault (fault,faultp,nfault,dt,time,xmin,xmax,ymin,ymax)

! interpolate result onto undeformed mesh

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Interpolate onto undeformed mesh'


          do i=1,nsurf
          ij=(i-1)*nz+1
          call interpolate (t(ij),tp(ij),z(ij),zp(ij),nz)
          enddo

! update ptt

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Update PTT'


        zsurfo=we1*zsurf(iconsurf(1,ieo))+we2*zsurf(iconsurf(2,ieo))+ &
               we3*zsurf(iconsurf(4,ieo))+we4*zsurf(iconsurf(4,ieo))
        zsurfpo=we1*zsurfp(iconsurf(1,ieo))+we2*zsurfp(iconsurf(2,ieo))+ &
                we3*zsurfp(iconsurf(4,ieo))+we4*zsurfp(iconsurf(4,ieo))

        do kstep=max(1,istep),nstep
          do i=1,nsurf
! Jean unchanged the position of the following lines after the computations for the Temperature
! as it appeared that the depth (and thus ages) where offset by one time step
          xxx=xdepth(i,kstep)
          yyy=ydepth(i,kstep)
          zzz=zdepth(i,kstep)
          do k=1,niteradvec
          call find_velo (xxx,yyy,zzz,time,dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,.true.)
          dxxx=xdepth(i,kstep)+dt*vx/2.-xxx
          dyyy=ydepth(i,kstep)+dt*vy/2.-yyy
          dzzz=zdepth(i,kstep)+dt*vz/2.-zzz
          xxx=xxx+dxxx
          yyy=yyy+dyyy
          zzz=zzz+dzzz
          enddo
          xdepth(i,kstep)=xdepth(i,kstep)+vx*dt
          ydepth(i,kstep)=ydepth(i,kstep)+vy*dt
          zdepth(i,kstep)=zdepth(i,kstep)+vz*dt+rsurf(i)+dt*usurf(i) ! VKP
          call find_element (xdepth(i,kstep),ydepth(i,kstep),zdepth(i,kstep),tnow,x,y,z,t, &
                             xsurf,ysurf,zsurf,ielsurf(i),neighbour, &
                             iconsurf,icon,nelemsurf,nelem, &
                             nsurf,nz,nnode,npe,mpe, &
                             interpol_success)
            if (interpol_success) then
            tprev(i,kstep)=tnow
            else
              if (zdepth(i,kstep).gt.zl) then
              tprev(i,kstep)=0.d0
              elseif (zdepth(i,kstep).lt.0.d0) then
              tprev(i,kstep)=tmax
              else
!              tprev(i,kstep)=tmax-tmax*zdepth(i,kstep)/zl
              tprev(i,kstep)=tmax-tmax*zdepth(i,kstep)/(zl+zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)
              endif
            endif
          enddo
        enddo

          do i=1,nobs
! Jean unchanged the position of the following lines after the computations for the Temperature
! as it appeared that the depth (and thus ages) where offset by one time step
          xxx=xdeptho(i)
          yyy=ydeptho(i)
          zzz=zdeptho(i)
          do k=1,niteradvec
          call find_velo (xxx,yyy,zzz,time,dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,.true.)
          dxxx=xdeptho(i)+dt*vx/2.-xxx
          dyyy=ydeptho(i)+dt*vy/2.-yyy
          dzzz=zdeptho(i)+dt*vz/2.-zzz
          xxx=xxx+dxxx
          yyy=yyy+dyyy
          zzz=zzz+dzzz
          enddo
          xdeptho(i)=xdeptho(i)+vx*dt
          ydeptho(i)=ydeptho(i)+vy*dt
          zdeptho(i)=zdeptho(i)+vz*dt+we1(i)*(rsurf(iconsurf(1,ieo(i)))+dt*usurf(iconsurf(1,ieo(i)))) &
                                      +we2(i)*(rsurf(iconsurf(2,ieo(i)))+dt*usurf(iconsurf(2,ieo(i)))) &
                                      +we3(i)*(rsurf(iconsurf(3,ieo(i)))+dt*usurf(iconsurf(3,ieo(i)))) &
                                      +we4(i)*(rsurf(iconsurf(npe,ieo(i)))+dt*usurf(iconsurf(npe,ieo(i))))
          call find_element (xdeptho(i),ydeptho(i),zdeptho(i),tnow,x,y,z,t, &
                             xsurf,ysurf,zsurf,ielsurf(i),neighbour, &
                             iconsurf,icon,nelemsurf,nelem, &
                             nsurf,nz,nnode,npe,mpe, &
                             interpol_success)
            if (interpol_success) then
            tprevo(i)=tnow
            else
              if (zdeptho(i).gt.zl) then
              tprevo(i)=0.d0
              elseif (zdeptho(i).lt.0.d0) then
              tprevo(i)=tmax
              else
!              tprevo(i)=tmax-tmax*zdeptho(i)/zl
              tprevo(i)=tmax-tmax*zdeptho(i)/(zl+zsurfpo(i)+(zsurfo(i)-zsurfpo(i))*ftime)
              endif
            endif
          enddo

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Update Level Set Functions'

        allocate (lsfn(nsurf,nz),lsfxn(nsurf,nz),lsfyn(nsurf,nz))
          do k=1,nz
            do i=1,nsurf
            xxx=xsurf(i)
            yyy=ysurf(i)
            zzz=(zl+zmax)*float(k-1)/(nz-1)
            call find_velo (xxx,yyy,zzz,time,dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,.false.)
            xxx=xxx-vx*dt
            yyy=yyy-vy*dt
            zzz=zzz-vz*dt-rsurf(i)-dt*usurf(i)
            lsfxn(i,k)=hinterpolate(lsfx,xxx,yyy,zzz,xmin,xmax,ymin,ymax,0.d0,zl+zmax,nx,ny,nz)
            lsfyn(i,k)=hinterpolate(lsfy,xxx,yyy,zzz,xmin,xmax,ymin,ymax,0.d0,zl+zmax,nx,ny,nz)
            lsfn(i,k)=hinterpolate(lsf,xxx,yyy,zzz,xmin,xmax,ymin,ymax,0.d0,zl+zmax,nx,ny,nz)
            enddo
          enddo
        lsf=lsfn
        lsfx=lsfxn
        lsfy=lsfyn
        deallocate (lsfn,lsfxn,lsfyn)

! saves tt and pt

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Save Time-Temperature paths'

        jrec=jrec+1
        tfake=tfinal-time
        if (jrec.eq.1) tfake=tprevious-time
          do kstep=max(1,istep),nstep
          write (100+kstep,rec=jrec) sngl(tfake),sngl(tprev(1:nsurf,kstep)), &
                sngl(zl+zsurfp(1:nsurf)+(zsurf(1:nsurf)-zsurfp(1:nsurf))*ftime-zdepth(1:nsurf,kstep)), &
                sngl(xdepth(1:nsurf,kstep)),sngl(ydepth(1:nsurf,kstep))
          enddo
        write (500,rec=jrec) sngl(tfake),sngl(tprevo(1:nobs)), &
              sngl(zl+zsurfpo(1:nobs)+(zsurfo(1:nobs)-zsurfpo(1:nobs))*ftime-zdeptho(1:nobs)), &
              sngl(xdeptho(1:nobs)),sngl(ydeptho(1:nobs))

! end of time stepping

        call cpu_time (times(4))

        enddo

! clean up

      if (istep.ne.0) then

      allocate(t2(nsurf),t2p(nsurf),t3(nsurf),t4(nsurf),t5(nsurf))
      t2p=0.
        do krec=jrec,1,-1
        read (100+istep,rec=krec) t1,t2,t3,t4,t5
          do i=1,nsurf
          if (t2(i).lt.0.) t2(i)=t2p(i)
          t2p(i)=t2(i)
          enddo
        write (100+istep,rec=krec) t1,t2,t3,t4,t5
        enddo
      deallocate (t2,t2p,t3,t4,t5)

!      jrec=jrec+1
      write (100+istep,rec=jrec+1) sngl(-1.d0),sngl(tprev(1:nsurf,istep)), &
           sngl(zl+zsurfp+(zsurf-zsurfp)*ftime-zdepth(1:nsurf,istep)), &
           sngl(xdepth(1:nsurf,istep)),sngl(ydepth(1:nsurf,istep))

      endif

! calculates ages

        if (nd.eq.0.and.ilog.eq.1) write (9,*) 'Calculate ages'

        allocate (age1(nsurf),age2(nsurf),age3(nsurf),age4(nsurf),age5(nsurf),nNs(nsurf))
        allocate (age6(nsurf),age7(nsurf),age8(nsurf),ftdist(20,nsurf),grainsize(nsurf))

        if (istep.eq.0.or.iout.eq.0) then

        age1=tprevious;age2=tprevious;age3=tprevious;age4=tprevious;age5=tprevious
        age6=tprevious;age7=tprevious;age8=tprevious;ftdist=0.;nNs=0.

        else

        if (nd.eq.0) then
        write (6,*) ''
        write (6,*) 'Calculating ages'
        endif

        grainsize = 100.
        call calculate_ages (jrec,nsurf,nz,istep,age1,age2,age3,age4,age5,age6,age7,age8, &
                             ftdist,tprevious,ftlflag,ageflag,iproc,nd,.FALSE.,grainsize)

        if (nd.eq.0.and.p%age_TL_flag.ne.0) then
        write (6,*) ''
        write (6,*) 'Calculating TL/OSL'
        call calculate_TL (jrec,nsurf,nz,istep,nNs,p%TL_doser,p%TL_D0,p%TL_a,p%TL_b,p%TL_Et,p%TL_logs,p%TL_logrho)
        endif

        endif

! write output

        if (iout.eq.1) then

        if (iproc.eq.0.and.nd.eq.0) then
        write (6,*) 'Saving at step ',istep
        endif
        irec=irec+1

        allocate (vxx(nnode),vyy(nnode),vzz(nnode),diff(nnode),het(nnode),zinit(nnode),xinit(nnode),yinit(nnode))
        ij=0
        iu=0 !VKP
            do i=1,nsurf
            iu=iu+1 !VKP
              do k=1,nz
              ij=ij+1
              call find_velo (x(ij),y(ij),z(ij),time,dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,.true.)
              vxx(ij)=vx;vyy(ij)=vy;vzz(ij)=vz+usurf(iu)
              if (dt.gt.tiny(dt)) vzz(ij)=vzz(ij)+rebound(ij)/dt !VKP
              xinit(ij)=hinterpolate (lsfx,x(ij),y(ij),z(ij),xmin,xmax,ymin,ymax,0.d0,zl+zmax,nx,ny,nz)
              yinit(ij)=hinterpolate (lsfy,x(ij),y(ij),z(ij),xmin,xmax,ymin,ymax,0.d0,zl+zmax,nx,ny,nz)
              zinit(ij)=hinterpolate (lsf,x(ij),y(ij),z(ij),xmin,xmax,ymin,ymax,0.d0,zl+zmax,nx,ny,nz)
                if (heatproduction.lt.0.d0) then
                het(ij)=heat_production(xinit(ij),yinit(ij),zinit(ij),t(ij),tfinal-time)
                else
                het(ij)=heatproduction
                endif
                if (diffusivity.lt.0.d0) then
                diff(ij)=thermal_diffusivity(xinit(ij),yinit(ij),zinit(ij),t(ij),tfinal-time)
                else
                diff(ij)=diffusivity
                endif
              enddo
            enddo
        if (nd.eq.0) write (8,rec=irec) nnode,nelem,mpe,sngl(tfinal-time),sngl(x),sngl(y),sngl(z),sngl(t), &
                           vxx,vyy,vzz,sngl(xinit),sngl(yinit),sngl(zinit),sngl(diff),sngl(het),((icon(j,i),j=1,mpe),i=1,nelem)
        deallocate (vxx,vyy,vzz,diff,het,xinit,yinit,zinit)

        allocate (exhumation(nsurf))
        exhumation=0.d0
          do i=1,nsurf
          dtimesurf=timesurf-timesurfp
          if (dtimesurf.gt.0.d0) exhumation(i)=-(zsurf(i)-zsurfp(i))/dtimesurf
          call find_velo (xsurf(i),ysurf(i),zl+zsurf(i),time,dt,fault,nfault,vx,vy,vz,xmin,xmax,ymin,ymax,.true.)
          exhumation(i)=exhumation(i)+vz+usurf(i) !VKP
          enddo
        if (nd.eq.0) write (11,rec=irec) nsurf,nelemsurf,npe,sngl(tfinal-time),sngl(xsurf),sngl(ysurf),sngl(zl+zsurf), &
                            exhumation,age1,age2,age3,age4,age5,age6,age7,age8, &
                            (sngl(sum(length*ftdist(:,j))),j=1,nsurf),nNs,iconsurf
        deallocate (exhumation)

        write (cstep,'(i3)') istep
        if (istep.lt.10) cstep(1:2)='00'
        if (istep.lt.100) cstep(1:1)='0'
        if (nd.eq.0) then
        open (12,file=run//'/output/Ages'//cstep//'.csv',status='unknown')
        write (12,'(a)') 'Longitude,Latitude,Height,HeApatite,HeZircon,FTApatite,' &
                             //'FTZircon,ArKFeldspar,ArBiotite,ArMuscovite,ArHornblend,' &
                             //'FTApMeanTL,nN'
          do i=1,nsurf
          xxx=xlonmin+(xlonmax-xlonmin)*(xsurf(i)-xmin)/(xmax-xmin)
          yyy=xlatmin+(xlatmax-xlatmin)*(ysurf(i)-ymin)/(ymax-ymin)
          write (12,'(g12.6,12(",",g12.6))') xxx,yyy,zsurf(i)*1000,age1(i),age2(i),age3(i),age4(i),age5(i), &
                                       age6(i),age7(i),age8(i),sngl(sum(length*ftdist(:,i))),nNs(i)
          enddo
        close (12)
  !      call vtkp (lsf,nx,ny,nz,xmax-xmin,ymax-ymin,zl+zmax,istep)
        endif

        endif

      if (istep.ne.nstep) deallocate (age1,age2,age3,age4,age5,age6,age7,age8,ftdist,grainsize,nNs)

! end of surface stepping

      if (istep.ne.0) close (100+istep)

      enddo

if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'
if (nd.eq.0) write (*,*) '------------------Computing misfit and writing output----------------------------'
if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'

    zsurfo=we1*zsurf(iconsurf(1,ieo))+we2*zsurf(iconsurf(2,ieo))+ &
           we3*zsurf(iconsurf(4,ieo))+we4*zsurf(iconsurf(4,ieo))
    zsurfpo=we1*zsurfp(iconsurf(1,ieo))+we2*zsurfp(iconsurf(2,ieo))+ &
            we3*zsurfp(iconsurf(4,ieo))+we4*zsurfp(iconsurf(4,ieo))
    write (500,rec=jrec+1) sngl(-1.d0),sngl(tprevo(1:nobs)), &
           sngl(zl+zsurfpo(1:nobs)+(zsurfo(1:nobs)-zsurfpo(1:nobs))*ftime-zdeptho(1:nobs)), &
           sngl(xdeptho(1:nobs)),sngl(ydeptho(1:nobs))

      if (iproc.eq.0.and.nd.eq.0) then
      close (8)
      open (8,file=run//'/output/Pecube.out',status='unknown',access='direct', &
            recl=4)
      write (8,rec=irec*(4+12*nnode+mpe*nelem)+1) -1
      close (8)
      endif

      if (iproc.eq.0.and.nd.eq.0) then
      close (11)
      open (11,file=run//'/output/Ages.out',status='unknown',access='direct', &
            recl=4)
      write (11,rec=irec*(4+14*nsurf+npe*nelemsurf)+1) -1
      close (11)
      endif

      if (nd.eq.0.and.ilog.eq.1) close (9)

! calculate misfit

!      if (iproc.eq.0) then
! added by Jean to change misfit to slope rather than exact ages
! (4/6/2008)
      if (nd.eq.0.and.nobs1.gt.0) open (13,file=run//'/output/CompareAGE.csv',status='unknown')
      misfit=0.
      nmisfit1=0
      misfit1 = 0.d0
      nhe=0
      nft=0
      nheZ=0
      nftZ=0
      narK=0
      narB=0
      narM=0
      narH=0
      read (7,*) nobs1,nobs2,nobs3,nobs4
      if (nobs1.lt.0) nobs1=-nobs1
      if (nobs2.lt.0) nobs2=-nobs2
      if (nobs3.lt.0) nobs3=-nobs3
      if (nobs4.lt.0) nobs4=-nobs4
      nobs=nobs1+nobs2+nobs3+nobs4
      allocate (lonobs(nobs),latobs(nobs))
      allocate (aheoreg(nobs),ahepreg(nobs),hheoreg(nobs))
      allocate (aftoreg(nobs),aftpreg(nobs),hftoreg(nobs))
      allocate (aheZoreg(nobs),aheZpreg(nobs),hheZoreg(nobs))
      allocate (aftZoreg(nobs),aftZpreg(nobs),hftZoreg(nobs))
      allocate (aarKoreg(nobs),aarKpreg(nobs),harKoreg(nobs))
      allocate (aarBoreg(nobs),aarBpreg(nobs),harBoreg(nobs))
      allocate (aarMoreg(nobs),aarMpreg(nobs),harMoreg(nobs))
      allocate (aarHoreg(nobs),aarHpreg(nobs),harHoreg(nobs))
      allocate (hpreg(nobs))
      allocate (hei(nobs),agehe(nobs),ageheZ(nobs),ageft(nobs),ageftZ(nobs), &
                agearB(nobs),agearK(nobs),agearM(nobs),agearH(nobs),ftdisto(20,nobs), &
                temphistp(1024,nobs))
      if (nd.eq.0.and.nobs1.gt.0) write (13,'(128(a,","))') 'LON','LAT','HEIGHTOBS','HEIGHTPRED', &
                  'AHEOBS','AHEPRED','AFTOBS','AFTPRED','ZHEOBS','ZHEPRED','ZFTOBS','ZFTPRED', &
                  'KAROBS','KARPRED','BAROBS','BARPRED','MAROBS','MARPRED','HAROBS','HARPRED', &
                  'FT01OBS','FT01PRED','FT02OBS','FT02PRED','FT03OBS','FT03PRED','FT04OBS','FT04PRED', &
                  'FT05OBS','FT05PRED','FT06OBS','FT06PRED','FT07OBS','FT07PRED','FT08OBS','FT08PRED', &
                  'FT09OBS','FT09PRED','FT10OBS','FT10PRED','FT11OBS','FT11PRED','FT12OBS','FT12PRED', &
                  'FT13OBS','FT13PRED','FT14OBS','FT14PRED','FT15OBS','FT15PRED','FT16OBS','FT16PRED', &
                  'FT17OBS','FT17PRED','FT18OBS','FT18PRED','FT19OBS','FT19PRED','FT20OBS','FT20PRED'
      if (saveTt.and.nd.eq.0) open (82,file=run//'/output/TimeTemperaturePaths.csv',status='unknown')
      ageflagp=1
        do iobs=1,nobs1
        read (7,*) dummy,dummy,dummy,dummy,dummy,dummy, &
        dummy,dummy,dummy,dummy,dummy,dummy, &
        dummy,dummy,dummy,dummy,dummy,dummy, &
        (dummy,i=1,20),dummy,grainsize(iobs),ieobs,wobs1,wobs2,wobs3,wobs4
        enddo
        do iobs=1,nobs1
        backspace (7)
        enddo
      call calculate_ages (jrec,nobs1,nz,400,agehe,ageheZ,ageft,ageftZ,agearK,agearB,agearM,agearH, &
                           ftdisto,tprevious,ftlflag,ageflag,iproc,nd,saveTt,grainsize)
      if (saveTt) close (82)

! compute misfit for ages

      nmisfit1 = 0
      misfit1 = 0.
        do iobs=1,nobs1
        read (7,*) xlonobs,xlatobs,ageheobs,dageheobs,ageftobs,dageftobs, & ! By Xav
                   ageheZobs,dageheZobs,ageftZobs,dageftZobs,agearKobs,dagearKobs, &
                   agearBobs,dagearBobs,agearMobs,dagearMobs,agearHobs,dagearHobs, &
                   tlobs,heightobs,dummy,ieobs,wobs1,wobs2,wobs3,wobs4   ! By Xav
        hei(iobs)=wobs1*zsurf(iconsurf(1,ieobs)) &
                 +wobs2*zsurf(iconsurf(2,ieobs)) &
                 +wobs3*zsurf(iconsurf(3,ieobs)) &
                 +wobs4*zsurf(iconsurf(npe,ieobs))
        hei(iobs)=hei(iobs)*1000.+min(0.d0,heightobs)
        lonobs(iobs)=xlonobs
        latobs(iobs)=xlatobs
        if (ageheobs.gt.0.) then
        nhe=nhe+1
        aheoreg(nhe)=ageheobs
        hheoreg(nhe)=hei(iobs)
        endif
        if (ageftobs.gt.0.) then
        nft=nft+1
        aftoreg(nft)=ageftobs
        hftoreg(nft)=hei(iobs)
        endif
        if (ageheZobs.gt.0.) then
        nheZ=nheZ+1
        aheZoreg(nheZ)=ageheZobs
        hheZoreg(nheZ)=hei(iobs)
        endif
        if (ageftZobs.gt.0.) then
        nftZ=nftZ+1
        aftZoreg(nftZ)=ageftZobs
        hftZoreg(nftZ)=hei(iobs)
        endif
        if (agearKobs.gt.0.) then
        narK=narK+1
        aarKoreg(narK)=agearKobs
        harKoreg(narK)=hei(iobs)
        endif
        if (agearBobs.gt.0.) then
        narB=narB+1
        aarBoreg(narB)=agearBobs
        harBoreg(narB)=hei(iobs)
        endif
        if (agearMobs.gt.0.) then
        narM=narM+1
        aarMoreg(narM)=agearMobs
        harMoreg(narM)=hei(iobs)
        endif
        if (agearHobs.gt.0.) then
        narH=narH+1
        aarHoreg(narH)=agearHobs
        harHoreg(narH)=hei(iobs)
        endif
        ahepreg(iobs)=agehe(iobs)
        aftpreg(iobs)=ageft(iobs)
        aheZpreg(iobs)=ageheZ(iobs)
        aftZpreg(iobs)=ageftZ(iobs)
        aarKpreg(iobs)=agearK(iobs)
        aarBpreg(iobs)=agearB(iobs)
        aarMpreg(iobs)=agearM(iobs)
        aarHpreg(iobs)=agearH(iobs)
        hpreg(iobs)=hei(iobs)
        if (nd.eq.0.and.nobs1.gt.0) write (13,'(128(g12.6,","))') xlonobs,xlatobs,heightobs,hei(iobs), &
                          ageheobs,agehe(iobs),ageftobs,ageft(iobs),ageheZobs,ageheZ(iobs), &
                          ageftZobs,ageftZ(iobs), &
                          agearKobs,agearK(iobs),agearBobs,agearB(iobs),agearMobs,agearM(iobs),agearHobs,agearH(iobs), &
                          (tlobs(i),ftdisto(i,iobs),i=1,20)
        if (ageheobs.gt.0.) then
        misfit1=misfit1+(ageheobs-agehe(iobs))**2/dageheobs**2
        nmisfit1=nmisfit1+1
        endif
        if (ageftobs.gt.0.) then
        misfit1=misfit1+(ageftobs-ageft(iobs))**2/dageftobs**2
        nmisfit1=nmisfit1+1
        endif
        if (ageheZobs.gt.0.) then
        misfit1=misfit1+(ageheZobs-ageheZ(iobs))**2/dageheZobs**2
        nmisfit1=nmisfit1+1
        endif
        if (ageftZobs.gt.0.) then
        misfit1=misfit1+(ageftZobs-ageftZ(iobs))**2/dageftZobs**2
        nmisfit1=nmisfit1+1
        endif
        if (agearKobs.gt.0.) then
        misfit1=misfit1+((agearKobs-agearK(iobs))/dagearKobs)**2
        nmisfit1=nmisfit1+1
        endif
        if (agearBobs.gt.0.) then
        misfit1=misfit1+(agearBobs-agearB(iobs))**2/dagearBobs**2
        nmisfit1=nmisfit1+1
        endif
        if (agearMobs.gt.0.) then
        misfit1=misfit1+(agearMobs-agearM(iobs))**2/dagearMobs**2
        nmisfit1=nmisfit1+1
        endif
        if (agearHobs.gt.0.) then
        misfit1=misfit1+(agearHobs-agearH(iobs))**2/dagearHobs**2
        nmisfit1=nmisfit1+1
        endif
        nmisfit1a=0
        misfit1a=0.
          if (tlobs(1).gt.-1.d-3) then
            do j=1,20
            misfit1a=misfit1a+(tlobs(j)-ftdisto(j,iobs))**2
            nmisfit1a=nmisfit1a+1
            enddo
          endif
        enddo

! if the flag mftflag is set, misfit1 is replaced by a misfit on the slope of the age-elevation relationship for each chronometer

  if (mftflag.eq.1) then
  nmisfit1=0
  misfit1=0.
  if (nhe.gt.1) call regression (aheoreg,hheoreg,nhe,aheom,hheom,heoi,heos)
  if (nft.gt.1) call regression (aftoreg,hftoreg,nft,aftom,hftom,ftoi,ftos)
  if (nheZ.gt.1) call regression (aheZoreg,hheZoreg,nheZ,aheZom,hheZom,heZoi,heZos)
  if (nftZ.gt.1) call regression (aftZoreg,hftZoreg,nftZ,aftZom,hftZom,ftZoi,ftZos)
  if (narK.gt.1) call regression (aarKoreg,harKoreg,narK,aarKom,harKom,arKoi,arKos)
  if (narB.gt.1) call regression (aarBoreg,harBoreg,narB,aarBom,harBom,arBoi,arBos)
  if (narM.gt.1) call regression (aarMoreg,harMoreg,narM,aarMom,harMom,arMoi,arMos)
  if (narH.gt.1) call regression (aarHoreg,harHoreg,narH,aarHom,harHom,arHoi,arHos)
  if (nobs.gt.1) call regression (ahepreg,hpreg,nobs,ahepm,hhepm,hepi,heps)
  if (nobs.gt.1) call regression (aftpreg,hpreg,nobs,aftpm,hftpm,ftpi,ftps)
  if (nobs.gt.1) call regression (aheZpreg,hpreg,nobs,aheZpm,hheZpm,heZpi,heZps)
  if (nobs.gt.1) call regression (aftZpreg,hpreg,nobs,aftZpm,hftZpm,ftZpi,ftZps)
  if (nobs.gt.1) call regression (aarKpreg,hpreg,nobs,aarKpm,harKpm,arKpi,arKps)
  if (nobs.gt.1) call regression (aarBpreg,hpreg,nobs,aarBpm,harBpm,arBpi,arBps)
  if (nobs.gt.1) call regression (aarMpreg,hpreg,nobs,aarMpm,harMpm,arMpi,arMps)
  if (nobs.gt.1) call regression (aarHpreg,hpreg,nobs,aarHpm,harHpm,arHpi,arHps)
  misfit=0.d0
    if (nhe.gt.1) then
    misfit1=misfit1+(abs(heos-heps)**2/abs(heos)**2+abs(aheom-ahepm)**2/abs(aheom)**2)
    nmisfit1=nmisfit1+2
    endif
    if (nft.gt.1) then
    misfit1=misfit1+(abs(ftos-ftps)**2/abs(ftos)**2+abs(aftom-aftpm)**2/abs(aftom)**2)
    nmisfit1=nmisfit1+2
    endif
    if (nheZ.gt.1) then
    misfit1=misfit1+(abs(heZos-heZps)**2/abs(heZos)**2+abs(aheZom-aheZpm)**2/abs(aheZom)**2)
    nmisfit1=nmisfit1+2
    endif
    if (nftZ.gt.1) then
    misfit1=misfit1+(abs(ftZos-ftZps)**2/abs(ftZos)**2+abs(aftZom-aftZpm)**2/abs(aftZom)**2)
    nmisfit1=nmisfit1+2
    endif
    if (narK.gt.1) then
    misfit1=misfit1+(abs(arKos-arKps)**2/abs(arKos)**2+abs(aarKom-aarKpm)**2/abs(aarKom)**2)
    nmisfit1=nmisfit1+2
    endif
    if (narB.gt.1) then
    misfit1=misfit1+(abs(arBos-arBps)**2/abs(arBos)**2+abs(aarBom-aarBpm)**2/abs(aarBom)**2)
    nmisfit1=nmisfit1+2
    endif
    if (narM.gt.1) then
    misfit1=misfit1+(abs(arMos-arMps)**2/abs(arMos)**2+abs(aarMom-aarMpm)**2/abs(aarMom)**2)
    nmisfit1=nmisfit1+2
    endif
    if (narH.gt.1) then
    misfit1=misfit1+(abs(arHos-arHps)**2/abs(arHos)**2+abs(aarHom-aarHpm)**2/abs(aarHom)**2)
    nmisfit1=nmisfit1+2
    endif
  endif

      if (nd.eq.0.and.nobs1.gt.0) close (13)

! compute misfit for thermal histories

      if (nd.eq.0.and.nobs2.gt.0) then
      open (13,file=run//'/output/CompareTT.csv',status='unknown')
      write (13,'(128(a,","))') 'LON','LAT','TIME','TEMP','TEMPPRED'
      endif

      nhist=0
        do iobs=nobs1+1,nobs1+nobs2
        read (7,*) xlonobs,xlatobs,nhist(iobs),(thist(j,iobs),temphist(j,iobs),errortemphist(j,iobs),j=1,nhist(iobs)), &
                   heightobs,ieobs,wobs1,wobs2,wobs3,wobs4   ! By Xav
        hei(iobs)=wobs1*zsurf(iconsurf(1,ieobs)) &
                 +wobs2*zsurf(iconsurf(2,ieobs)) &
                 +wobs3*zsurf(iconsurf(3,ieobs)) &
                 +wobs4*zsurf(iconsurf(npe,ieobs))
        hei(iobs)=hei(iobs)*1000.+min(0.d0,heightobs)
        lonobs(iobs)=xlonobs
        latobs(iobs)=xlatobs
        enddo

      call calculate_ttpath (jrec,nobs,nz,400,nhist,thist,temphistp,iproc,nd)

      nmisfit2 = 0
      misfit2 = 0.
        do iobs=nobs1+1,nobs1+nobs2
        if (nd.eq.0.and.nobs2.gt.0) write (13,'(g15.9,",",g15.9)') lonobs(iobs),latobs(iobs)
          do j=1,nhist(iobs)
            if (temphistp(j,iobs).gt.-9998.d0) then
            misfit2=misfit2+(temphist(j,iobs)-temphistp(j,iobs))**2/errortemphist(j,iobs)**2
            nmisfit2=nmisfit2+1
            endif
          if (nd.eq.0.and.nobs2.gt.0) write (13,'(",",3(",",g15.9))') thist(j,iobs),temphist(j,iobs),temphistp(j,iobs)
          enddo
        enddo

      if (nd.eq.0.and.nobs2.gt.0) close (13)

! compute misfit for 43He

      if (nd.eq.0.and.nobs3.gt.0) then
      open (13,file=run//'/output/Compare43HE.csv',status='unknown')
      write (13,'(128(a,","))') 'LON','LAT','AGEHE','AGEHEPRED','%RELEASED','%RELEASEDPRED','AGERELEASED','AGERELEASEDPRED'
      endif

      nheating=0
        do iobs=nobs1+nobs2+1,nobs1+nobs2+nobs3
        read (7,*) xlonobs,xlatobs,size43He(iobs),age43He(iobs),dage43He(iobs),nheating(iobs), &
                  (theating(j,iobs),duration(j,iobs),released(j,iobs),dreleased(j,iobs), &
                   agereleased(j,iobs),dagereleased(j,iobs),j=1,nheating(iobs)), &
                   heightobs,ieobs,wobs1,wobs2,wobs3,wobs4   ! By Xav
        hei(iobs)=wobs1*zsurf(iconsurf(1,ieobs)) &
                 +wobs2*zsurf(iconsurf(2,ieobs)) &
                 +wobs3*zsurf(iconsurf(3,ieobs)) &
                 +wobs4*zsurf(iconsurf(npe,ieobs))
        hei(iobs)=hei(iobs)*1000.+min(0.d0,heightobs)
        lonobs(iobs)=xlonobs
        latobs(iobs)=xlatobs
        enddo

      allocate (time43He(jrec,nobs),temperature43He(jrec,nobs))
      call extract_ttpath (jrec,nobs,nz,400,time43He,temperature43He,iproc,nd)

      nmisfit3 = 0
      misfit3 = 0.
      nmisfit3a = 0
      misfit3a = 0.
        do iobs=nobs1+nobs2+1,nobs1+nobs2+nobs3
        allocate (releasedp(nheating(iobs)),agereleasedp(nheating(iobs)))
        call Diffusion1D (time43He(1,iobs),temperature43He(1,iobs),jrec,theating(1,iobs),duration(1,iobs),nheating(iobs), &
                          age43Hep,releasedp,agereleasedp,size43He(iobs))
        if (nd.eq.0.and.nobs2.gt.0) write (13,'(g15.9,",",g15.9",",g15.9",",g15.9)') &
              lonobs(iobs),latobs(iobs),age43He(iobs),age43Hep
          do i=1,nheating(iobs)
          misfit3a=misfit3a+((releasedp(i)-released(i,iobs))**2/dreleased(i,iobs)**2+ &
                           (agereleasedp(i)-agereleased(i,iobs))**2/dagereleased(i,iobs)**2)
          nmisfit3a=nmisfit3a+1
          if (nd.eq.0.and.nobs3.gt.0) write (13,'(",,,",4(",",g15.9))') released(i,iobs),releasedp(i), &
                        agereleased(i,iobs),agereleasedp(i)
          enddo
        misfit3=misfit3+(age43Hep-age43He(iobs))**2/dage43He(iobs)**2
        nmisfit3=nmisfit3+1
        deallocate (releasedp,agereleasedp)
        enddo

      deallocate (time43He,temperature43He)

      if (nd.eq.0.and.nobs3.gt.0) close (13)

! compute misfit for TL

      if (nd.eq.0.and.nobs4.gt.0) then
      open (13,file=run//'/output/CompareTL.csv',status='unknown')
      if (nd.eq.0.and.nobs4.gt.0) write (13,'(128(a,","))') 'LON','LAT','DURATION','TEMPERATURE','NN','NNPRED'
      endif

      nheating=0
        do iobs=nobs1+nobs2+nobs3+1,nobs1+nobs2+nobs3+nobs4
        read (7,*) xlonobs,xlatobs,nheating(iobs), &
                  (doser(j,iobs),d0(j,iobs),radius(j,iobs),et(j,iobs), &
                  logs(j,iobs),b(j,iobs),logrho(j,iobs),nn(j,i),j=1,nheating(iobs)), &
                  heightobs,ieobs,wobs1,wobs2,wobs3,wobs4   ! By Xav
        hei(iobs)=wobs1*zsurf(iconsurf(1,ieobs)) &
                  +wobs2*zsurf(iconsurf(2,ieobs)) &
                  +wobs3*zsurf(iconsurf(3,ieobs)) &
                  +wobs4*zsurf(iconsurf(npe,ieobs))
        hei(iobs)=hei(iobs)*1000.+min(0.d0,heightobs)
        lonobs(iobs)=xlonobs
        latobs(iobs)=xlatobs
        enddo

      allocate (timeTSL(jrec,nobs),temperatureTSL(jrec,nobs))
      call extract_ttpath (jrec,nobs,nz,400,timeTSL,temperatureTSL,iproc,nd)

      misfit4 = 0.
      nmisfit4 = 0
        do iobs=nobs1+nobs2+nobs3+1,nobs1+nobs2+nobs3+nobs4
        write (13,'(g15.9,",",g15.9)') lonobs(iobs),latobs(iobs)
          do i=1,nheating(iobs)
          paramsTSL(1)=doser(i,iobs)
          paramsTSL(2)=d0(i,iobs)
          paramsTSL(3)=radius(i,iobs)
          paramsTSL(4)=et(i,iobs)
          paramsTSL(5)=0.
          paramsTSL(6)=logs(i,iobs)
          paramsTSL(7)=b(i,iobs)
          paramsTSL(8)=logrho(i,iobs)
          nnf = TLModel (timeTSL(1,iobs), temperatureTSL(1,iobs), jrec, paramsTSL)
          misfit4=misfit4+(nn(i,iobs) - nnf)**2/(nn(i,iobs)*0.1d0)**2
          nmisfit4=nmisfit4+1
          if (nd.eq.0.and.nobs4.gt.0) write (13,'(",",4(",",g15.9))') timeTSL(i,iobs),temperatureTSL(i,iobs),nn(i,iobs),nnf
          enddo
        enddo

      deallocate (timeTSL,temperatureTSL)

      if (nd.eq.0.and.nobs4.gt.0) close (13)

      deallocate (lonobs,latobs)

      if (nmisfit1.ne.0) misfit1 = sqrt(misfit1/nmisfit1)
      if (nmisfit1a.ne.0) misfit1a = sqrt(misfit1a/nmisfit1a)
      if (nmisfit2.ne.0) misfit2 = sqrt(misfit2/nmisfit2)
      if (nmisfit3.ne.0) misfit3 = sqrt(misfit3/nmisfit3)
      if (nmisfit3a.ne.0) misfit3a = sqrt(misfit3a/nmisfit3a)
      if (nmisfit4.ne.0) misfit4 = sqrt(misfit4/nmisfit4)
      misfit = misfit1*p%misfit_weight_AGE + misfit1a * p%misfit_weight_FTLD + misfit2*p%misfit_weight_TH + &
            (misfit3 + misfit3a)*p%misfit_weight_43HE + misfit4*p%misfit_weight_TL
        if (p%misfit_corrected.eq.1) then
          if (nmisfit1 + nmisfit1a + nmisfit2 + nmisfit3 + nmisfit3a + nmisfit4 - nd -1 .gt. 0) then
          misfit = misfit/(nmisfit1 + nmisfit1a + nmisfit2 + nmisfit3 + nmisfit3a + nmisfit4 - nd - 1)
          else
            if (iproc.eq.0) then
            write (*,*) 'Warning - Pecube will not correct the misfit (misfit_corrected = 1) when the total number'
            write (*,*) 'of data points is smaller than the number of parameters being inverted  + 1'
            write (*,*) 'Currently it is:',nmisfit1 + nmisfit1a + nmisfit2 + nmisfit3 + nmisfit3a + nmisfit4
            write (*,*) 'and the number of parameters is ',nd
            endif
          endif
        endif

        if (nd.ne.0) then
        open (71,file=run//'/NA/NA_int_res.csv',status='unknown',position='append')
        write (71,'(g15.9,99(a1,g15.9))') misfit,(",",param(i),i=1,nd)
        close (71)
          if (p%save_ages_inversion.eq.1) then
            if (p%age_AHe_flag.ne.0.) then
            open (71,file=run//'/NA/AgeApatiteHelium.csv',status='unknown',position='append')
            write (71,'(g15.9,99(a1,g15.9))') misfit,(",",agehe(i),i=1,nhe)
            close (71)
            endif
            if (p%age_AFT_flag.ne.0.) then
            open (71,file=run//'/NA/AgeApatiteFT.csv',status='unknown',position='append')
            write (71,'(g15.9,99(a1,g15.9))') misfit,(",",ageft(i),i=1,nft)
            close (71)
            endif
            if (p%age_ZHe_flag.ne.0.) then
            open (71,file=run//'/NA/AgeZirconHelium.csv',status='unknown',position='append')
            write (71,'(g15.9,99(a1,g15.9))') misfit,(",",ageheZ(i),i=1,nheZ)
            close (71)
            endif
            if (p%age_ZFT_flag.ne.0.) then
            open (71,file=run//'/NA/AgeZirconFT.csv',status='unknown',position='append')
            write (71,'(g15.9,99(a1,g15.9))') misfit,(",",ageftZ(i),i=1,nftZ)
            close (71)
            endif
            if (p%age_KAr_flag.ne.0.) then
            open (71,file=run//'/NA/AgeKSparArgon.csv',status='unknown',position='append')
            write (71,'(g15.9,99(a1,g15.9))') misfit,(",",agearK(i),i=1,narK)
            close (71)
            endif
            if (p%age_BAr_flag.ne.0.) then
            open (71,file=run//'/NA/AgeBiotiteArgon.csv',status='unknown',position='append')
            write (71,'(g15.9,99(a1,g15.9))') misfit,(",",agearB(i),i=1,narB)
            close (71)
            endif
            if (p%age_MAr_flag.ne.0.) then
            open (71,file=run//'/NA/AgeMuscoviteArgon.csv',status='unknown',position='append')
            write (71,'(g15.9,99(a1,g15.9))') misfit,(",",agearM(i),i=1,narM)
            close (71)
            endif
            if (p%age_HAr_flag.ne.0.) then
            open (71,file=run//'/NA/AgeHornblendeArgon.csv',status='unknown',position='append')
            write (71,'(g15.9,99(a1,g15.9))') misfit,(",",agearH(i),i=1,narH)
            close (71)
            endif
          endif
        endif
      if (nd.eq.0) write (6,*) 'Misfit : ',misfit
      if (nd.ne.0) write (*,*) 'Misfit:',misfit,misfit1,misfit2,misfit3,misfit4,' - Params',param(1:nd)
      close (13)
! added by Jean to change misfit to slope rather than exact ages
! (4/6/2008)
      deallocate (aheoreg,ahepreg,hheoreg)
      deallocate (aftoreg,aftpreg,hftoreg)
      deallocate (aheZoreg,aheZpreg,hheZoreg)
      deallocate (aftZoreg,aftZpreg,hftZoreg)
      deallocate (aarKoreg,aarKpreg,harKoreg)
      deallocate (aarBoreg,aarBpreg,harBoreg)
      deallocate (aarMoreg,aarMpreg,harMoreg)
      deallocate (aarHoreg,aarHpreg,harHoreg)
      deallocate (hei,agehe,ageheZ,ageft,ageftZ,agearB,agearK,agearM,agearH,ftdisto)
      deallocate (hpreg)
!      endif

      deallocate (age1,age2,age3,age4,age5,age6,age7,age8,grainsize)

      deallocate (iconsurf,ielsurf,neighbour)
      deallocate (xsurf,ysurf,zsurf,zsurfp)
      deallocate (topoa,topob,rsurf,rebound)
      deallocate (xdepth,ydepth,zdepth,tprev)
      deallocate (xexhumation,yexhumation,zexhumation)

      deallocate (x,y,z,t)
      deallocate (xp,yp,zp,tp)
      deallocate (icon,ael,bel)
      deallocate (kfix)
      deallocate (f)
      deallocate (proc,ice)

      call cpu_time (times(9))

if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'
      if (iproc.eq.0.and.nd.eq.0) write (6,*) 'Total compute time : ',times(9)-times(1),' sec'

      close (7)

if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'
if (nd.eq.0) write (*,*) '-------------------------End of Pecube run---------------------------------------'
if (nd.eq.0) write (*,*) '---------------------------------------------------------------------------------'

      return

      end

!--------------------------------------------------------------------------------------------

      subroutine regression (x,y,n,xmean,ymean,intercept,slope)

! this routine is used to calculate the means, intercept and slope
! of the regression between two datasets of length n stored in x and y

      implicit none

      integer n
      double precision x(n),y(n),xmean,ymean,intercept,slope

      xmean=sum(x)/n
      ymean=sum(y)/n
      slope=sum((x-xmean)*(y-ymean))/sum((x-xmean)**2)
      intercept=ymean-xmean*slope

      return
      end
