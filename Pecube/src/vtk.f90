program vtk

implicit none

real,dimension(:),allocatable::x,y,z,t,vx,vy,vz,xinit,yinit,zinit,diffusivity,heat
real,dimension(:),allocatable::xs,ys,zs,ex,a1,a2,a3,a4,a5,a6,a7,a8,He43,ftlm,nNtl,nNosl,nNesr
integer,dimension(:,:),allocatable::icon,iconsurf
real tim
integer nnode,nelem,nsurf,nelemsurf,mpe,npe,irec,nrec,itime,i,j,k,iunit,i1,i2,istep
integer nnx
character cs*3,ct*3,run*5,iq*1
integer numarg
logical is_unix

numarg=command_argument_count()
if (numarg.eq.0) then
print*,'You need to specify a run directory (i.e. RUN00 for example)'
stop 'End of run'
else
call getarg (1,run)
endif

! Check Os system
call GetOperatingSystem(is_unix)
if (is_unix) then
    call system ('mkdir -p tmp')
else
    call RemoveDirectory(run//'\VTK')
    call system ('mkdir '//run//'\VTK')
endif

write (*,*) '---------------------------------------------------------------------------------'
write (*,*) '---------------------------- Make VTK files -------------------------------------'
write (*,*) '---------------------------------------------------------------------------------'

print*,'Doing temperature and velocity'

open (77,file=run//'/output/Pecube.out',status='old',access='direct',recl=16)
read (77,rec=1) nnode,nelem,mpe,tim
close (77)

open (77,file=run//'/output/Pecube.out',status='old',access='direct',recl=4)
nrec=1
1 read (77,rec=(nrec-1)*(4+12*nnode+mpe*nelem)+1) nnx
if (nnx.lt.0) goto 999
print*,'Record',nrec
nrec=nrec+1
goto 1
999 nrec=nrec-1
close (77)

open (77,file=run//'/output/Pecube.out',status='old',access='direct',recl=4*(4+12*nnode+mpe*nelem))

allocate (x(nnode),y(nnode),z(nnode),t(nnode))
allocate (vx(nnode),vy(nnode),vz(nnode))
allocate (heat(nnode),diffusivity(nnode),xinit(nnode),yinit(nnode),zinit(nnode))
allocate (icon(mpe,nelem))

do irec=0,nrec-1

print*,'doing record ',irec,'out of',nrec

read (77,rec=irec+1) nnode,nelem,mpe,tim,(x(i),i=1,nnode),(y(i),i=1,nnode),(z(i),i=1,nnode),(t(i),i=1,nnode), &
                     (vx(i),i=1,nnode),(vy(i),i=1,nnode),(vz(i),i=1,nnode),(xinit(i),i=1,nnode),(yinit(i),i=1,nnode), &
                     (zinit(i),i=1,nnode),(diffusivity(i),i=1,nnode),(heat(i),i=1,nnode),(icon(1:mpe,i),i=1,nelem)

write(cs,'(i3)') irec
if (irec.lt.10) cs(1:2)='00'
if (irec.lt.100) cs(1:1)='0'
itime=int(tim)
write(ct,'(i3)') itime
if (itime.lt.10) ct(1:2)='00'
if (itime.lt.100) ct(1:1)='0'

iunit=30

open(unit=iunit,file=run//'/VTK/Pecube'//cs//'.vtk')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'velocities'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
write(iunit,'(a7,i10,a6)')'POINTS ',nnode,' float'

  do i=1,nnode
  write(iunit,'(3f16.11)') x(i),y(i),z(i)
  enddo

print*,mpe

  if (mpe.eq.8) then

  write(iunit,'(A6, 2I10)') 'CELLS ',nelem,9*nelem
    do i=1,nelem
    write(iunit,'(9I10)')8,icon(1,i)-1,icon(2,i)-1,icon(3,i)-1,icon(4,i)-1, &
                           icon(5,i)-1,icon(6,i)-1,icon(7,i)-1,icon(8,i)-1
    enddo
  write(iunit,'(A11, I10)') 'CELL_TYPES ',nelem
    do k=1,nelem
    write(iunit,'(I2)')12 ! octree  (8 nodes)
    enddo

  else

  write(iunit,'(A6, 2I10)') 'CELLS ',nelem,7*nelem
    do i=1,nelem
    write(iunit,'(9I10)')6,icon(1,i)-1,icon(3,i)-1,icon(2,i)-1, &
                           icon(4,i)-1,icon(6,i)-1,icon(5,i)-1
    enddo
  write(iunit,'(A11, I10)') 'CELL_TYPES ',nelem
    do k=1,nelem
    write(iunit,'(I2)')13 ! octree  (6 nodes)
    enddo

 endif

write(iunit,'(a11,i10)')'POINT_DATA ',nnode

write(iunit,'(a)')'SCALARS Temperature float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nnode
    write(iunit,'(f18.13)') t(i)
    enddo

write(iunit,'(a)')'SCALARS Xinit float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nnode
    write(iunit,'(f18.13)') xinit(i)
    enddo

write(iunit,'(a)')'SCALARS Yinit float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nnode
    write(iunit,'(f18.13)') yinit(i)
    enddo

write(iunit,'(a)')'SCALARS Zinit float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nnode
    write(iunit,'(f18.13)') zinit(i)
    enddo

write(iunit,'(a)')'SCALARS Diffusivity float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nnode
    write(iunit,'(f18.13)') diffusivity(i)
    enddo

write(iunit,'(a)')'SCALARS Heat float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nnode
    write(iunit,'(f18.13)') heat(i)
    enddo

write(iunit,'(a)')'VECTORS velo float'
    do i=1,nnode
    write(iunit,'(3f18.13)') vx(i),vy(i),vz(i)
    enddo

close (iunit)

enddo

deallocate (x,y,z,t,vx,vy,vz,xinit,yinit,zinit,diffusivity,heat)

close (77)

print*,'Doing Ages and topo'

open (77,file=run//'/output/Ages.out',status='old',access='direct',recl=16)
read (77,rec=1) nsurf,nelemsurf,npe,tim
close (77)

open (77,file=run//'/output/Ages.out',status='old',access='direct',recl=4)
nrec=1
2 read (77,rec=(nrec-1)*(4+17*nsurf+npe*nelemsurf)+1) nnx
if (nnx.lt.0) goto 998
print*,'Record',nrec
nrec=nrec+1
goto 2
998 nrec=nrec-1
close (77)

open (77,file=run//'/output/Ages.out',status='old',access='direct',recl=4*(4+17*nsurf+npe*nelemsurf))

allocate (xs(nsurf),ys(nsurf),zs(nsurf),ex(nsurf))
allocate (a1(nsurf),a2(nsurf),a3(nsurf),a4(nsurf),a5(nsurf))
allocate (a6(nsurf),a7(nsurf),a8(nsurf),nNtl(nsurf),nNosl(nsurf),nNesr(nsurf))
allocate (ftlm(nsurf),He43(nsurf))
allocate (iconsurf(npe,nelemsurf))

do irec=0,nrec-1

print*,'doing record ',irec,'out of',nrec

read (77,rec=irec+1) nsurf,nelemsurf,npe,tim,xs,ys,zs,ex,a1,a2,a3,a4,a5,a6,a7,a8,He43,ftlm,nNtl,nNosl,nNesr,iconsurf

write(cs,'(i3)') irec
if (irec.lt.10) cs(1:2)='00'
if (irec.lt.100) cs(1:1)='0'
itime=int(tim)
write(ct,'(i3)') itime
if (itime.lt.10) ct(1:2)='00'
if (itime.lt.100) ct(1:1)='0'

iunit=30

open(unit=iunit,file=run//'/VTK/Ages'//cs//'.vtk')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'velocities'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
write(iunit,'(a7,i10,a6)')'POINTS ',nsurf,' float'

    do i=1,nsurf
    write(iunit,'(3f16.10)') xs(i),ys(i),zs(i)
    enddo

  if (npe.eq.4) then

  write(iunit,'(A6, 2I10)') 'CELLS ',nelemsurf,5*nelemsurf
    do i=1,nelemsurf
    write(iunit,'(9I10)')4,iconsurf(1:npe,i)-1
    enddo

  write(iunit,'(A11, I10)') 'CELL_TYPES ',nelemsurf
    do k=1,nelemsurf
    write(iunit,'(I2)')9 ! octree  (8 nodes)
    enddo

  else

  write(iunit,'(A6, 2I10)') 'CELLS ',nelemsurf,4*nelemsurf
    do i=1,nelemsurf
    write(iunit,'(9I10)')3,iconsurf(1:npe,i)-1
    enddo

  write(iunit,'(A11, I10)') 'CELL_TYPES ',nelemsurf
    do k=1,nelemsurf
    write(iunit,'(I2)')5 ! octree  (8 nodes)
    enddo

  endif

write(iunit,'(a11,i10)')'POINT_DATA ',nsurf

write(iunit,'(a)')'SCALARS ExhumationRate float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
    write(iunit,'(f18.13)') ex(i)
    enddo

write(iunit,'(a)')'SCALARS ApatiteHeAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (a1(i).lt.1e-4.or.a1(i).gt.1e3) then
        a1(i) = 0.d0
      endif
    write(iunit,'(f18.13)') a1(i)
    enddo

write(iunit,'(a)')'SCALARS ZirconHeAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (a2(i).lt.1e-4.or.a2(i).gt.1e3) then
        a2(i) = 0.d0
      endif
    write(iunit,'(f18.13)') a2(i)
    enddo

write(iunit,'(a)')'SCALARS ApatiteFTAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (a3(i).lt.1e-4.or.a3(i).gt.1e3) then
        a3(i) = 0.d0
      endif
    write(iunit,'(f18.13)') a3(i)
    enddo

write(iunit,'(a)')'SCALARS ZirconFTAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (a4(i).lt.1e-4.or.a4(i).gt.1e3) then
        a4(i) = 0.d0
      endif
    write(iunit,'(f18.13)') a4(i)
    enddo

write(iunit,'(a)')'SCALARS KsparArgonAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (a5(i).lt.1e-4.or.a5(i).gt.1e3) then
        a5(i) = 0.d0
      endif
    write(iunit,'(f18.13)') a5(i)
    enddo

write(iunit,'(a)')'SCALARS BiotiteArgonAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (a6(i).lt.1e-4.or.a6(i).gt.1e3) then
        a6(i) = 0.d0
      endif
    write(iunit,'(f18.13)') a6(i)
    enddo

write(iunit,'(a)')'SCALARS MuscoviteArgonAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (a7(i).lt.1e-4.or.a7(i).gt.1e3) then
        a7(i) = 0.d0
      endif
    write(iunit,'(f18.13)') a7(i)
    enddo

write(iunit,'(a)')'SCALARS HornblendeArgonAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (a8(i).lt.1e-4.or.a8(i).gt.1e3) then
        a8(i) = 0.d0
      endif
    write(iunit,'(f18.13)') a8(i)
    enddo
    
write(iunit,'(a)')'SCALARS EdgeAge float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (He43(i).lt.1e-4.or.He43(i).gt.1e3) then
        He43(i) = 0.d0
      endif
    write(iunit,'(f18.13)') He43(i)
    enddo

write(iunit,'(a)')'SCALARS ApatiteMeanFTLength float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nsurf
    if (ftlm(i).lt.1e-4.or.ftlm(i).gt.1e3) then
        ftlm(i) = 0.d0
      endif
  write(iunit,'(f18.13)') ftlm(i)
  enddo

write(iunit,'(a)')'SCALARS TLnN float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do i=1,nsurf
    if (nNtl(i).lt.1e-4.or.nNtl(i).gt.1e3) then
        nNtl(i) = 0.d0
      endif
  write(iunit,'(f18.13)') nNtl(i)
  enddo

  write(iunit,'(a)')'SCALARS OSLnN float 1'
  write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (nNosl(i).lt.1e-4.or.nNosl(i).gt.1e3) then
        nNosl(i) = 0.d0
      endif
    write(iunit,'(f18.13)') real(nNosl(i))
    enddo
    
  write(iunit,'(a)')'SCALARS ESRnN float 1'
  write(iunit,'(a)')'LOOKUP_TABLE default'
    do i=1,nsurf
      if (nNesr(i).lt.1e-4.or.nNesr(i).gt.1e3) then
        nNesr(i) = 0.d0
      endif
    write(iunit,'(f18.13)') nNesr(i)
    enddo           
   
close (iunit)

enddo

close (77)

end program vtk
