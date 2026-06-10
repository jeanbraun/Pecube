subroutine testTopo (p,nx,ny,nz,xl,yl,zl,timesteps,dosteps,nstep,xlon1,xlat1,xlon2,xlat2,fltflag)

use Pecube
use DEM

implicit none

type (parameters) p
character line*1024,run*5,fnme*300,obsfile*300,cs*3,dum_char*50
integer nx0,ny0,nskip,nstep,istep,isoflag,nxiso,nyiso,nz,i,j,k,ii,jj,nx,ny,nobs,ij,nobs1,nobs2,nobs3,nobs4,nobs5,nobs6
integer nfnme,nobsfile,iunit,icon1,icon2,icon3,icon4,itime,fltflag,mftflag,ftlflag
double precision dx,dy,xlon,xlat,tau,rhoc,rhom,young,poisson,thickness,dum,tprevious
double precision crustal_thickness,diffusivity,tmax,tmsl,tlapse,heatproduction
double precision zl,xl,yl,x,y,h
double precision topoRef !Maxime
double precision,dimension(:),allocatable::timek,topomag,topooffset,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
double precision,dimension(:),allocatable::timek_temp2,topomag_temp2,topooffset_temp2,x_array,phase_array_temp2
double precision,dimension(:),allocatable::timek_temp,topomag_temp,topooffset_temp,phase_array_temp
double precision,dimension(:),allocatable::phase_array
integer,dimension(:),allocatable::iout,index_array,iout_temp,index_array2,iout_temp2
double precision,dimension(:),allocatable::Zincision, amplification, NewTopo
double precision time_incision_start, time_incision_stop, minZ, maxZ, tauH, Di !Headward propagation
integer, dimension(:),allocatable::io
double precision,dimension(:,:),allocatable::zNZ
double precision,dimension(:),allocatable::z,xz,yz,zsmooth,ztemp
double precision timesteps(1000)
integer dosteps(1000), Headward, res
double precision xlon1,xlat1,xlon2,xlat2,xl0,yl0
integer,dimension(:,:),allocatable::iconz

logical vivi,xyz,topoDEM
character c5*5
integer nc5

run = p%run_name

  do i=1,300
  fnme(i:i)=' '
  enddo
fnme = p%topo_file_name
  do i=1,300
  if (fnme(i:i).ne.' ') nfnme=i
  enddo
vivi=.FALSE.
if (fnme(nfnme:nfnme).eq.'/') vivi=.TRUE.
if (vivi) nfnme=nfnme-1
xyz=.FALSE.
if (fnme(nfnme:nfnme).eq.'$') xyz=.TRUE.
if (xyz) nfnme=nfnme-1
topoDEM = .FALSE.
nx0 = p%nx
ny0 = p%ny
dx = p%dlon
dy = p%dlat
nskip = p%nskip
xlon = p%lon0
xlat = p%lat0
nz = p%nz

nstep = p%ntime+sum(p%nstep)*2 ! Maxime
tau = p%erosional_time_scale

allocate (io(nstep+1))

! ! For tests, output topography and fault for all times provided (topo and tectonic)
! do istep=1,nstep+1
!   timek(istep) = p%time_topo(istep)
!   topomag(istep) = p%amplification(istep)
!   topooffset(istep) = p%offset(istep)
!   dosteps(istep) = 1 !p%output(istep) - Maxime
!   io(istep) = p%output(istep)
! enddo

      
! ! converts geological time into model time
! do istep=nstep+1,1,-1
!   timek(istep)=timek(1)-timek(istep)
!   timesteps(istep)=timek(istep)
! enddo

! From create_input
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
    iout_temp(istep) = 1
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
    iout_temp2(istep) = 1
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
        iout_temp2(ij) = 1
      endif
      ! check if time_end is already in time topo
      res = 0
      call isinarray(timek_temp2,p%time_end(istep,k), nstep,res)
      if (res.eq.0) then 
          ij = ij+1
        timek_temp2(ij) = p%time_end(istep,k)
        iout_temp2(ij) = 1
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
  iout(:) = 1
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
  io(:) = iout(index_array)
  dosteps(1:nstep+1) = io(:)

  topomag = topomag(index_array)
  topooffset = topooffset(index_array)
  ! do k=1, ij ! Make sure amplification is not negative
  !     if (topomag(k).lt.0.d0) then
  !         topomag(k) = 0.d0
  !     endif
  !   enddo
  ! print *, 'timek = ',timek
  ! print *, 'topomag = ',topomag
  ! print *, 'topooffset = ',topooffset
  ! print *, 'iout = ',io

  nstep = ij-1

  deallocate (timek_temp,topomag_temp,topooffset_temp,iout_temp,phase_array_temp)
  deallocate (timek_temp2,topomag_temp2,topooffset_temp2,iout_temp2,phase_array_temp2)
    
! converts geological time into model time
  do istep=nstep+1,1,-1
    timek(istep)=timek(1)-timek(istep)
    timesteps(istep)=timek(istep)
  enddo


! Headward erosion ?
Headward = p%do_Headward
tauH = p%tauH
Di = p%depth_incision
time_incision_start = p%time_incision_start
time_incision_stop = p%time_incision_stop


crustal_thickness = p%thickness
crustal_thickness=crustal_thickness*1.d3

  do i=1,300
  obsfile(i:i)=' '
  enddo
  obsfile = p%data_folder
  do i=1,300
  if (obsfile(i:i).ne.' ') nobsfile=i
  enddo

ftlflag = p%FT_code_flag
mftflag = p%misfit_slope
fltflag = p%fault_advect_flag
tprevious = p%default_age
if (tprevious.eq.0.) tprevious=timek(nstep+1)

999 close (77)

if (.not.vivi) then

  if (nx0.gt.0) then

  allocate (zNZ(nx0,ny0))
    if (fnme(1:nfnme).eq.'Nil') then
    zNZ=0.d0
    elseif (fnme(1:nfnme).eq.'Topo30') then
    dx = 360.d0/43200
    dy = dx
    call ExtractDEM (xlon, xlat, nx0, ny0, PecubeFnme = run//'/data/ExtractedDEM')
    open (8,file=run//'/data/ExtractedDEM.dat',status='old')
    read (8,*) zNZ
    close (8)
    else
    open (8,file = run//'/data/'//fnme(1:nfnme),status='old')
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
      topoDEM= .TRUE.
      endif
    close (8)
    endif

  nx=(nx0-1)/nskip+1
  ny=(ny0-1)/nskip+1
  allocate (z(nx*ny),ztemp(nx*ny),Zincision(nx*ny),NewTopo(nx*ny))
  allocate (zsmooth(nx*ny))                    
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
  xlon2=xlon+(nx0-1)*dx
  xlat2=xlat+(ny0-1)*dy
  
  else

  allocate (z(-nx0),xz(-nx0),yz(-nx0),iconz(3,-ny0),Zincision(-nx0),NewTopo(-nx0))
    open (8,file=run//'data/'//fnme(1:nfnme),status='old')
      do i=1,-nx0
      read (8,*) z(i)
      enddo
    close (8)
    open (8,file=run//'data/'//fnme(1:nfnme)//'.geometry',status='old')
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

else

  if (nx0.gt.0) then

  nx=(nx0-1)/nskip+1 !VKP
  ny=(ny0-1)/nskip+1 !VKP
  allocate (zNZ(nx0,ny0)) !VKP
  allocate (z(nx*ny),ztemp(nx*ny),Zincision(nx*ny),NewTopo(nx*ny)) !VKP
  allocate (zsmooth(nx*ny))                    
  topomag=1.d0 !VKP
  topooffset=0.d0 !VKP

  xlon1=xlon
  xlat1=xlat
  xlon2=xlon+(nx0-1)*dx
  xlat2=xlat+(ny0-1)*dy

  else

  allocate (z(-nx0),xz(-nx0),yz(-nx0),iconz(3,-ny0),Zincision(-nx0),NewTopo(-nx0))
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

endif

  if (nx0.gt.0) then

  xl=dx*(nx-1)*nskip*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
  yl=dy*(ny-1)*nskip*111.11
  zl=crustal_thickness/1.e3
  z=z/crustal_thickness*zl
  ztemp = z

  else

  zl=crustal_thickness/1.e3
  z=z/crustal_thickness*zl
  ztemp = z
  xl=(xlon2-xlon1)*111.11*cos((xlat+(xlat2-xlat1)/2.)*3.141592654/180.)
  yl=(xlat2-xlat1)*111.11

  endif


if (istep.eq.0) then
        ! Initiate amplication matrix (Maxime)
    allocate (amplification(nx*ny))
    do i=1,nx*ny
        amplification(i) =  topomag(1)
    enddo
    
    !Maxime - Check the reference elevation from which to compute topographic evolution
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
    elseif (p%topo_ref.eq.5) then !from mean elevation
      topoRef = sum(z(:))/(nx*ny)
    endif !end Maxime
    ! Get topography before incision (in case of headward propagation scenario)
    if (Headward.eq.1) then ! Headward propagation of erosion
        call get_initial_topo(topomag, topooffset, maxval(timek)-timek, nstep, topoRef, z, nx*ny,&
                                time_incision_start,time_incision_stop, Zincision)
        minZ = minval(Zincision(:))
        maxZ = maxval(Zincision(:))
        NewTopo(:) = topoRef-(topomag(1)*(topoRef-z))+topooffset(1) ! Topo for headward propagation scenario
    endif   
 endif
  
  
iunit=30
do istep=0,nstep
if (io(istep+1).eq.1) then
write(cs,'(i3)') istep
if (istep.lt.10) cs(1:2)='00'
if (istep.lt.100) cs(1:1)='0'
itime=int(timesteps(istep+1))

  if (vivi) then !VKP
      if (istep.lt.10) then !VKP
      write (c5(1:3),'(a2,i1)') '00',istep !VKP
      nc5=3 !VKP
      elseif (istep.lt.100) then !VKP
      write (c5(1:3),'(a1,i2)') '00',istep !VKP
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
    open (67,file=run//'/data/'//fnme(1:nfnme)//'/topo'//c5(1:nc5),status='old') !VKP
    read (67,*) zNZ !VKP
    ij=0
      do j=1,ny0,nskip  !VKP
        do i=1,nx0,nskip  !VKP
        ij=ij+1
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

print*,'creating Topo file'
if (topomag(istep+1).lt.0) call smooth_topo(z, zsmooth, nx, ny, int(-topomag(istep+1)),&
                                            int(topooffset(istep+1)))                                                                                                                
open(unit=iunit,file=run//'/VTK/Topo'//cs//'.vtk',status='unknown')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'surface'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'

if (nx0.gt.0) then
   write(iunit,'(a7,i10,a6)')'POINTS ',nx*ny,' float'
    
    if (p%topo_ref.lt.5.and.Headward.eq.1) then
        if (istep.gt.0) then 
            call Headward_propagation(topomag,topoRef, NewTopo, topooffset, istep, nstep, amplification,&
                                    nx*ny, Di, time_incision_start,time_incision_stop, Zincision,&
                                    z, maxval(timek(:))-timek(istep+1), maxval(timek(:))-timek(istep), minZ, maxZ, tauH)
        endif
    endif
    
    ij=0
    do j=1,ny
        do i=1,nx
            x=xl*float(i-1)/float(nx-1)
            y=yl*float(j-1)/float(ny-1)
            ij=ij+1
            if (topomag(istep+1).lt.0) then
              write(iunit,'(3f16.11)') x,y,zsmooth(ij)+zl
            else                             
               if (p%topo_ref.lt.5.and.Headward.eq.0) then
                   write(iunit,'(3f16.11)') x,y,topoRef-(topomag(istep+1)*(topoRef-z(ij)))+topooffset(istep+1)+zl !Maxime
               elseif (p%topo_ref.lt.5.and.Headward.eq.1) then
                write(iunit,'(f18.13)') x,y,NewTopo(ij)
               else
                   ! Compute new topo
                   if (z(ij).gt.topoRef.and.topomag(istep+1).ne.0) then
                      ztemp(ij) = topoRef-(1/topomag(istep+1)*(topoRef-z(ij)))+topooffset(istep+1) +zl
        !               if (ztemp(ij).gt.maxval(z(:))) then
        !                  ztemp(ij) = maxval(z(:))
        !               endif
                   else
                        ztemp(ij) = topoRef-(topomag(istep+1)*(topoRef-z(ij)))+topooffset(istep+1) +zl
                   endif
                   write(iunit,'(3f16.11)') x,y,ztemp(ij)   
               endif
            endif
        enddo
  enddo
  
  write(iunit,'(A6, 2I10)') 'CELLS ',(nx-1)*(ny-1),(1+4)*(nx-1)*(ny-1)
  do j=1,ny-1
    do i=1,nx-1
        icon1=(j-1)*nx+i-1
        icon2=icon1+1
        icon3=icon1+nx+1
        icon4=icon1+nx
        write (iunit,'(9I10)') 4,icon1,icon2,icon3,icon4
    enddo
  enddo
  write(iunit,'(A11, I10)') 'CELL_TYPES ',(nx-1)*(ny-1)
  do k=1,(nx-1)*(ny-1)
    write(iunit,'(I2)')9 ! rectangles
  enddo
  write(iunit,'(a11,i10)')'POINT_DATA ',nx*ny
  write(iunit,'(a)')'SCALARS Topo float 1'
  write(iunit,'(a)')'LOOKUP_TABLE default'
    
  ij=0
  do j=1,ny
    do i=1,nx
    ij=ij+1
    if (topomag(istep+1).lt.0.) then
      write(iunit,'(f18.13)') zsmooth(ij)
    else
        if (p%topo_ref.lt.5.and.Headward.eq.0) then
            write(iunit,'(f18.13)') topoRef-(topomag(istep+1)*(topoRef-z(ij)))+topooffset(istep+1)!Maxime
        elseif (p%topo_ref.lt.5.and.Headward.eq.1) then
            write(iunit,'(f18.13)') NewTopo(ij)
        else
            ! Compute new topo
            if (z(ij).gt.topoRef.and.topomag(istep+1).ne.0) then
               ztemp(ij) = topoRef-(1/topomag(istep+1)*(topoRef-z(ij)))+topooffset(istep+1) 
!                if (ztemp(ij).gt.maxval(z(:))) then
!                   ztemp(ij) = maxval(z(:))
!                endif
            else
                 ztemp(ij) = topoRef-(topomag(istep+1)*(topoRef-z(ij)))+topooffset(istep+1) 
            endif
            write(iunit,'(f18.13)') ztemp(ij)
       endif                                
    endif 
    enddo
  enddo
else
  write(iunit,'(a7,i10,a6)')'POINTS ',-nx0,' float'
  do i=1,-nx0
  if (topomag(istep+1).lt.0.) then
    write(iunit,'(3f16.11)') (xz(i)-xlon1)/(xlon2-xlon1)*xl, &
                            (yz(i)-xlat1)/(xlat2-xlat1)*yl, &
                            zsmooth(i)+zl
  else                          
    write(iunit,'(3f16.11)') (xz(i)-xlon1)/(xlon2-xlon1)*xl, &
                           (yz(i)-xlat1)/(xlat2-xlat1)*yl, &
                           topoRef-(topomag(istep+1)*(topoRef-z(ij)))+topooffset(istep+1)+zl
  endif
  enddo
write(iunit,'(A6, 2I10)') 'CELLS ',-ny0,(1+3)*(-ny0)
    do i=1,-ny0
    write (iunit,'(9I10)') 3,iconz(:,i)-1
    enddo
write(iunit,'(A11, I10)') 'CELL_TYPES ',-ny0
  do k=1,-ny0
  write(iunit,'(I2)')5 ! rectangles
  enddo
write(iunit,'(a11,i10)')'POINT_DATA ',-nx0
write(iunit,'(a)')'SCALARS Topo float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  if (topomag(istep+1).lt.0.) then                            
    do i=1,-nx0
    write(iunit,'(f18.13)') zsmooth(i)
    enddo
  else
    do i=1,-nx0                                    
    write(iunit,'(f18.13)') topoRef-(topomag(istep+1)*(topoRef-z(ij)))+topooffset(istep+1)
    enddo
  endif  
endif
close(iunit)
print*,'Step',istep,'done'
endif
enddo

if (obsfile(1:nobsfile).eq.'Nil') return

!call read_data_folder ( run//'/data/'//obsfile(1:nobsfile), .FALSE., xlon1, xlon2, xlat1, xlat2, 0, 0)
call read_data_folder (run//'/data/'//obsfile(1:nobsfile), .FALSE., p%lon0, p%lon0+(nx-1)*dx*nskip,&
      p%lat0, p%lat0+(ny-1)*dy*nskip, 0, 0, nobs1, nobs2, nobs3, nobs4, nobs5, nobs6) !Maxime - changed xlon1,xlon2,xlat1,xlat2 into p%**

open (88, file = run//'/data/'//obsfile(1:nobsfile)//'0000.txt',status='old')

!open (88,file=run//'/data/'//obsfile(1:nobsfile)//'.txt',status='old')
!read (88,*) nobs1
!nobs1=iabs(nobs1)
!  do i=1,nobs1
!  read (88,*)
!  enddo
!nobs2=0
!read (88,*,end=1111) nobs2
!1111 nobs2=iabs(nobs2)
nobs=nobs1+nobs2+nobs3+nobs4+nobs5+nobs6
!rewind (88)
!read (88,*)
allocate (a1(nobs),a2(nobs),a3(nobs),a4(nobs),a5(nobs),a6(nobs),a7(nobs),a8(nobs),a9(nobs),a10(nobs),a11(nobs))

a1=-1.d0
a2=-1.d0
a3=-1.d0
a4=-1.d0
a5=-1.d0
a6=-1.d0
a7=-1.d0
a8=-1.d0
a9=-1.d0

open(unit=iunit,file=run//'/VTK/Data.vtk')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'AgeData'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
write(iunit,'(a7,i10,a6)')'POINTS ',nobs,' float'
  do i=1,nobs1
  read (88,*) x,y,h,a1(i),dum,a2(i),dum,a3(i),dum,a4(i),dum,a5(i),dum,a6(i),dum,a7(i),dum,a8(i),dum
  !print *, 'h = ', h
  !Ensure xlonobs is correct relative to xlon - Maxime
  if (topoDEM) then
      xl0 = xl
      yl0 = yl
      xl = p%lon0+(p%nx-1)*dx 
      yl = p%lat0+(p%ny-1)*dy
      x = xlon1 + abs(xl - p%lon0) * ((x - p%lon0) / (xl - p%lon0))
      y = xlat1 + abs(yl - p%lat0) * ((y - p%lat0) / (yl - p%lat0))
      xl = xl0
      yl = yl0
  endif
  x=(x-xlon1)/(xlon2-xlon1)*xl
  y=(y-xlat1)/(xlat2-xlat1)*yl
  write(iunit,'(3f16.11)') x,y,abs(h)/1.e3+zl ! elevation h is absolute value when there is non-zero topography, for borehole take true value
  enddo
a1=max(a1,-1.d0)
a2=max(a2,-1.d0)
a3=max(a3,-1.d0)
a4=max(a4,-1.d0)
a5=max(a5,-1.d0)
a6=max(a6,-1.d0)
a7=max(a7,-1.d0)
a8=max(a8,-1.d0)
!read (88,*,end=1112)
!1112 continue

  do i=nobs1+1,nobs1+nobs2+nobs3
  read (88,*) x,y,h
  !Ensure xlonobs is correct relative to xlon - Maxime
  if (topoDEM) then
      xl0 = xl
      yl0 = yl
      xl = p%lon0+(p%nx-1)*dx 
      yl = p%lat0+(p%ny-1)*dy
      x = xlon1 + abs(xl - p%lon0) * ((x - p%lon0) / (xl - p%lon0))
      y = xlat1 + abs(yl - p%lat0) * ((y - p%lat0) / (yl - p%lat0))
      xl = xl0
      yl = yl0
  endif
  x=(x-xlon1)/(xlon2-xlon1)*xl
  y=(y-xlat1)/(xlat2-xlat1)*yl
  write(iunit,'(3f16.11)') x,y,abs(h)/1.e3+zl
  enddo

  do i=nobs1+nobs2+nobs3+1,nobs1+nobs2+nobs3+nobs4
    read (88,*) dum_char,x,y,h,dum,dum,dum,dum,dum,dum,dum,a9(i),dum
    x=(x-xlon1)/(xlon2-xlon1)*xl
    y=(y-xlat1)/(xlat2-xlat1)*yl
    write(iunit,'(3f16.11)') x,y,abs(h)/1.e3+zl
    enddo
a9=max(a9,-1.d0)
do i=nobs1+nobs2+nobs3+nobs4+1,nobs1+nobs2+nobs3+nobs4+nobs5
  read (88,*) dum_char,x,y,h,dum,dum,dum,dum,dum,dum,dum,a10(i),dum,dum
  x=(x-xlon1)/(xlon2-xlon1)*xl
  y=(y-xlat1)/(xlat2-xlat1)*yl
  write(iunit,'(3f16.11)') x,y,abs(h)/1.e3+zl
  enddo
a10=max(a10,-1.d0)

do i=nobs1+nobs2+nobs3+nobs4+nobs5+1,nobs1+nobs2+nobs3+nobs4+nobs5+nobs6
  read (88,*) dum_char,x,y,h,dum,dum,dum,dum,dum,a11(i),dum,dum
  x=(x-xlon1)/(xlon2-xlon1)*xl
  y=(y-xlat1)/(xlat2-xlat1)*yl
  write(iunit,'(3f16.11)') x,y,abs(h)/1.e3+zl
  enddo
a11=max(a11,-1.d0)        
write(iunit,'(A6, 2I10)') 'CELLS ',1,nobs+1
write (iunit,'(256I10)') nobs,(i-1,i=1,nobs)
write(iunit,'(A11, I10)') 'CELL_TYPES ',1
write (iunit,'(I10)') 2
write(iunit,'(a11,i10)')'POINT_DATA ',nobs
write(iunit,'(a)')'SCALARS ApatiteHeAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a1(i)
  enddo
write(iunit,'(a)')'SCALARS ApatiteFTAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a2(i)
  enddo
write(iunit,'(a)')'SCALARS ZirconHeAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a3(i)
  enddo
write(iunit,'(a)')'SCALARS ZirconFTAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a4(i)
  enddo
write(iunit,'(a)')'SCALARS KSparArAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a5(i)
  enddo
write(iunit,'(a)')'SCALARS BiotiteArAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a6(i)
  enddo
write(iunit,'(a)')'SCALARS MuscoviteArAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a7(i)
  enddo
write(iunit,'(a)')'SCALARS HornblendeArAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a8(i)
  enddo
write(iunit,'(a)')'SCALARS NnTL float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a9(i)
  enddo
write(iunit,'(a)')'SCALARS NnOSL float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nobs
  write (iunit,'(f18.5)') a10(i)
  enddo
  write(iunit,'(a)')'SCALARS NnESR float 1'
  write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nobs
    write (iunit,'(f18.5)') a11(i)
    enddo
close(iunit)
print *,'To the end'
deallocate (ztemp,Zincision,NewTopo,amplification)

end

