module definitions

  type faulttype

! A fault vertical cross-sectional (2D) geometry is defined by a set of points (x,y)
! Its 3D location is defined by two points in the horizontal plane (x1,y1),(x2,y2)
! Those points limit its lateral extent
! In the (x,y) coordinate system (in which the faul vertical geometry is defined)
! x lies to the right of the (x1,y1),(x2,y2) line and z is vertical upwards
! and zero at the surface

! n is the number of points used to define the fault
! x(n) and y(n) are the coordinates of the points
! (x1,y1) and (x2,y2) are the coordinates of the trace of the fault at the surface
! nstep is the number of time intevals used to define the motion story of the fault
! per interval: timestart,timeend a,d velo
! timestart and timeend are the start and end time of the interval
! velo is the velocity across the fault (along the plane of the fault) for the interval

! By convention as one goes along the fault from point 1 to point n
! it is the block to the right that moves

! The sign of velo determines whether the sign of the x-component of velocity

! from that information we compute xs(n-1),ys(n-1), the direction
! of each of the (n-1) segments connecting the n points

! as well as xn,yn the normal to the (x1,y1)-(x2,y2) line in the z=0 plane

  integer n,nstep
  double precision,dimension(:),pointer::x,y
  double precision x1,y1,x2,y2
  double precision,dimension(:),pointer::timestart,timeend,velo
  double precision,dimension(:),pointer::xs,ys
  double precision xn,yn

  end type faulttype

!=====[EDGE]====================================================================
!this type is to store edges in a trianglulation
! it is used to update (in a generalized Delaunay sense)
! the triangulation of the 3D points on the surfaces
! for each edge:
! n1, n2 are the node numbers defining the edge
! t1, t2 are the triangle numbers on either side of the edge
! going from n1 to n2, t1 is to the left and t2 is to the right
! m1, m2 are the node numbers of the two other nodes making t1 and t2

      type edge
      integer n1,n2,m1,m2,t1,t2
      end type edge

end module definitions

!---------------------------

program showNA

      use definitions

      type (faulttype) fault(1024)
      real*4 param(1024),range(2,1024)
      integer nd,nfault,nd0
      character iq*1,ci*2,cj*2,ck*2

      double precision junk
      integer nr,i,j,k,ii
      double precision,dimension(:),allocatable::misfit
      double precision,dimension(:,:),allocatable::val
      iproc=0
      nproc=1

      if (iproc.eq.0) then
      nd=0
      if (iproc.eq.0) then
      open (77,file='input/fault_parameters.txt',status='old')
11    read (77,'(a1)') iq
      if (iq.eq.'$' .or. iq.eq.' ') goto 11
      backspace (77)
      read (77,*) nfault
      if (nfault.eq.0) nfault=1
      close (77)
      endif
      call create_pecube_in (fault,nfault,nd,range,param)
      close (7)
      endif

      if (nd.eq.0) then
      print*,'No NA results to show...'
      else
      print*,nd,'parameters'
      print*,'range:'
        do i=1,nd
        print*,range(1,i),range(2,i)
        enddo
      endif

open (7,file='NA/NA_results.txt',status='old')
nr=0
1 read (7,*,end=2) junk
nr=nr+1
goto 1
2 print*,nr,'model runs',nd,'model parameters'
allocate (misfit(nr),val(max(3,nd),nr))
rewind (7)
  do i=1,nr
  read (7,*) misfit(i),(val(j,i),j=1,nd)
  enddo
close (7)
  do i=1,nr
    do j=1,nd
    val(j,i)=(val(j,i)-range(1,j))/(range(2,j)-range(1,j))
    enddo
  enddo

nd0=nd
4 continue
if (nd.lt.3) then
nd=nd+1
val(nd,:)=0.5
goto 4
endif

do i=1,nd
do j=i+1,nd
do k=j+1,nd

write (ci,'(i2)') i
if (i.lt.10) ci(1:1)='0'
write (cj,'(i2)') j
if (j.lt.10) cj(1:1)='0'
write (ck,'(i2)') k
if (k.lt.10) ck(1:1)='0'

iunit=77

if (nd0.eq.1) then
open(unit=iunit,file='VTK/NA_'//ci//'.vtk')
elseif (nd0.eq.1) then
open(unit=iunit,file='VTK/NA_'//ci//'_'//cj//'.vtk')
else
open(unit=iunit,file='VTK/NA_'//ci//'_'//cj//'_'//ck//'.vtk')
endif
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'NA Results parameters '//ci//cj//ck
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
write(iunit,'(a7,i10,a6)')'POINTS ',nr,' float'
  do ii=1,nr
  write(iunit,'(3f16.11)') val(i,ii),val(j,ii),val(k,ii)
  enddo
write(iunit,'(A6, 2I10)') 'CELLS ',1,nr+1
write (iunit,'(256I10)') nr,(ii-1,ii=1,nr)
write(iunit,'(A11, I10)') 'CELL_TYPES ',1
write (iunit,'(I10)') 2
write(iunit,'(a11,i10)')'POINT_DATA ',nr
write(iunit,'(a)')'SCALARS Misfit float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
  do ii=1,nr
  write (iunit,'(f18.5)') misfit(ii)
  enddo
close(iunit)

enddo
enddo
enddo

end program showNA
