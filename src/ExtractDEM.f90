!-----------------------------------------------------------------

subroutine extractDEM (lonmin0, latmin0, nx, ny, demout, PecubeFnme, DEMFnme)

! subroutine to extract a DEM from the SRTM30-PLUS dataset
! given a range of longitudes and latitudes
! all of the output options are optional:
! if demout is given, a double precision array of the required size is created
! and returned; it contains the DEM
! if PecubeFnme is given, a file is created that is named PecubeFnme.dat
! and contains the DEM in Pecube format
! if DEMFnme is given, two files are created that are named DEMFnme.DEM and DEMFnme.HDR
! that contain the dem and a header to locate it in lat long

implicit none

double precision, intent(in) :: lonmin0, latmin0
integer, intent(in) :: nx,ny
double precision, dimension(:,:), pointer, optional, intent(out) :: demout
character*(*), optional, intent(out) :: PecubeFnme, DEMFnme

integer*2, dimension(:), allocatable :: row
double precision, dimension(:,:), allocatable :: dem
integer :: ncol = 43200
integer :: nrow = 21600
integer :: i, j, imin, imax, jmin, jmax
double precision :: lonmin,latmin

lonmin = lonmin0
latmin = latmin0

if (lonmin < 0.) lonmin = lonmin + 360.d0

imin = int(lonmin*120.d0) + 1
imin = min (imin, ncol)

jmin = int((90.d0 - latmin)*120.d0) + 1
jmin = min (jmin, nrow)

!nx = imax - imin + 1
!ny = jmin - jmax + 1
imax = imin + nx - 1
jmax = jmin - ny + 1

allocate (row(ncol))

allocate (dem(nx, ny))

open (7, file='topo/topo30', status='old', access='direct', recl=ncol*2, &
          convert='big_endian')
  do j = jmin, jmax, -1
  read (7, rec = j) row
  dem(1:nx, j-jmax+1) = row(imin:imax)
  end do
close (7)

  if (present(demout)) then

  allocate (demout(nx,ny))
  demout = dem

  endif

  if (present(PecubeFnme)) then

  open (10, file = PecubeFnme//'.dat', status='unknown')

    do j = ny, 1, -1
      do i= 1, nx
      write (10, '(f12.6)') dem(i,j)
      end do
    end do

  close (10)

  endif

  if (present(DEMFnme)) then

  open (8 , file = DEMFnme//'.DEM', status='unknown', access = 'direct', &
            recl = 2*nx, convert='big_endian')
    do j = 1, ny
    write (8, rec = j) (int2(dem( i, j)), i =1, nx)
    end do
  close (8)

  open (9, file = DEMFnme//'.HDR', status='unknown')
  write (9, '(a)') 'BYTEORDER      M'
  write (9, '(a)') 'LAYOUT       BIL'
  write (9, '(a, i5)') 'NROWS         ', ny
  write (9, '(a, i5)') 'NCOLS         ', nx
  write (9, '(a)') 'NBANDS        1'
  write (9, '(a)') 'NBITS         16'
  write (9, '(a, i6)') 'BANDROWBYTES         ', 2*ny
  write (9, '(a, i6)') 'TOTALROWBYTES        ', 2*nx
  write (9, '(a)') 'BANDGAPBYTES         0'
  write (9, '(a)') 'NODATA        -9999'
  write (9, '(a, f16.10)') 'ULXMAP        ', lonmin
  write (9, '(a, f16.10)') 'ULYMAP        ', latmin + 0.00833333333333 * (ny-1)
  write (9, '(a)') 'XDIM          0.00833333333333'
  write (9, '(a)') 'YDIM          0.00833333333333'
  close (9)

  endif

deallocate (row, dem)

end subroutine ExtractDEM
