!----------------------------

module read_data_module

!----------------------------

interface read_data_files

subroutine read_data_files (fnme,t,nx,ny,f,s)

character*(*) :: fnme
character*100, dimension(:,:), allocatable :: t
character*100, dimension(:), allocatable :: f,s
integer :: nx,ny

end subroutine read_data_files

end interface

!----------------------------

interface ucase

function ucase(in) result (out)

character (*), intent(in)  :: in
character(:), allocatable  :: out

end function ucase

end interface

!----------------------------

end module read_data_module

!----------------------------

module read_string_module

!----------------------------

interface read_string

!----------------------------

function read_string (unit, istring, jstring) result (out)

integer, intent(in)  :: unit, istring, jstring
character(:), allocatable  :: out

end function read_string

end interface

!----------------------------

end module read_string_module

!----------------------------

subroutine read_data_folder ( folder, echo, xlon1, xlon2, xlat1, xlat2, iproc, nd, &
                              nobs1, nobs2, nobs3, nobs4, nobs5, nobs6)

! subroutines goes throuh the files in forlder "folder" and attempts to read the data in each of them

use read_data_module

implicit none

character(len=*), intent(in) :: folder
logical, intent(in) :: echo
double precision, intent(in) :: xlon1, xlon2, xlat1, xlat2
integer, intent(in) :: iproc, nd

character(len=1024) :: fnme
logical french
character*100, dimension(:,:), allocatable :: t
character*100, dimension(:), allocatable :: f,s,sample
integer :: nx,ny,nsample,nfield
integer :: i,j
double precision, dimension(:,:), allocatable :: obs
character*6, dimension(:), allocatable :: field_name
character*4 :: cproc

integer nobs1, nobs2, nobs3, nobs4,nobs5,nobs6

if (nd.eq.0) print*,'-----------------------------Reading data----------------------------------------'
if (nd.eq.0) print*,'---------------------------------------------------------------------------------'

nfield = 40
allocate (field_name(nfield))

field_name = (/'LON   ','LAT   ','HEIGHT','AHE   ','DAHE  ','AFT   ','DAFT  ','ZHE   ','DZHE  ', &
'ZFT   ','DZFT  ','KAR   ','DKAR  ','BAR   ','DBAR  ','MAR   ','DMAR  ','HAR   ','DHAR  ', &
'TL1   ','TL2   ','TL3   ','TL4   ','TL5   ','TL6   ','TL7   ','TL8   ','TL9   ','TL10  ', &
'TL11  ','TL12  ','TL13  ','TL14  ','TL15  ','TL16  ','TL17  ','TL18  ','TL19  ','TL20  ', &
'SIZE  '/)

french=.false.

write (cproc,'(i4)') iproc
if (iproc.lt.10) cproc(1:3)='000'
if (iproc.lt.100) cproc(1:2)='00'
if (iproc.lt.1000) cproc(1:1)='0'

call system ('mkdir -p tmp')

call system ('ls '//folder//' > tmp/data_lst'//cproc//'.txt')

  if (nd.eq.0) then
  print*,'Reading data from folder: ',folder
  endif

! first find number of samples

nsample = 0
allocate (sample(1000))
open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=999)
  do
  read (8,'(a)', end=998) fnme
  call read_data_files (folder//'/'//fnme, t, nx, ny, f, s)
    do i = 1, ny
      do j = 1, nsample
      if (trim(s(i)) == sample(j)) goto 997
      enddo
    nsample = nsample + 1
    sample(nsample) = trim(s(i))
997 continue
    enddo
  deallocate (t,f,s)
  enddo
998 continue
close (8)
999 continue

! second fill the obs array with data from the files

allocate (obs(nfield,nsample))
obs = -9999.

open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=996)
  do
  read (8,'(a)', end=995) fnme
  call read_data_files (folder//'/'//fnme, t, nx, ny, f, s)
  call add_to_obs (t, nx, ny, f, s, obs, sample, nsample, fnme, field_name, nfield)
  deallocate (t,f,s)
  enddo
995 continue
close (8)
996 continue

! third check that data is complete

if (iproc.eq.0) open (123, file = folder//'_Messages.txt', status = 'unknown')
call check_obs (nsample, sample, obs, field_name, nfield, xlon1, xlon2, xlat1, xlat2, iproc)
if (iproc.eq.0) close (123)

! HERE WE WRITE DATA TO PECUBE INPUT FILE

open (111, file = folder//cproc//'.txt')
!write (111,*) nsample
  do i = 1, nsample
  write (111,*) (obs(j,i), j = 1, nfield)
  enddo
close (111)
nobs1 = nsample

  if (nd.eq.0.and.iproc.eq.0) print*,'Number of Age Samples:',nsample

! now scans folder for thermal histories

call read_data_files_for_thermal_histories (folder,cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs2)

! now scans folder for thermal histories

call read_data_files_for_43He (folder,cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs3)

! now scans for TSL data

call read_data_files_for_TSL (folder,cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs4)

! now scans for OSL data
call read_data_files_for_OSL (folder,cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs5)

! now scans for ESR data
call read_data_files_for_ESR (folder,cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs6)

deallocate (sample, field_name)

if (nd.eq.0) then
print*,'---------------------------------------------------------------------------------'
endif

return

end subroutine read_data_folder

!----------------------------

subroutine add_to_obs (t, nx, ny, f, s, obs, sample, nsample, fnme, field_name, nfield)

use read_data_module

integer :: nx, ny, nsample, nfield
character*100, dimension(nx,ny) :: t
character*100, dimension(nx) :: f
character*100, dimension(ny) :: s
character*100, dimension(nsample) :: sample
character(len=1024) :: fnme
double precision, dimension(nfield, nsample) :: obs
integer :: i,j,jsample,it
character*6, dimension(nfield) :: field_name
double precision val
integer, dimension(:,:), allocatable :: nobs

allocate (nobs(nfield, nsample))
nobs = 0

! gets data from file table into observation array

  do j = 1, ny
    do jsample = 1, nsample
      if (trim(s(j))==trim(sample(jsample))) then
        do i=1, nx
          do it = 1, nfield
            if (trim(ucase(f(i))).eq.trim(field_name(it))) then
              read (t(i,j),*,end=999,err=998) val
              if (obs(it,jsample).lt.-9998.) then
              obs(it,jsample) = val
              nobs(it,jsample) = 1
              else
              obs(it,jsample) = obs(it,jsample) + val
              nobs(it,jsample) = nobs(it,jsample) + 1
              endif
  999       continue
            endif
          enddo
        enddo
      endif
    enddo
  enddo

where (nobs.ne.0) obs = obs/nobs

deallocate (nobs)

!  do jsample = 1, nsample
!  if (obs(1,jsample).lt.0.d0) obs(1,jsample) = obs(1,jsample) + 360.d0
!  enddo

return

998 continue
print*,'Error in reading '//field_name(it)//' for sample '//sample(jsample)
stop 'Wrong format'

end subroutine add_to_obs

!----------------------------

! Here we run a few tests on the data to make sure that they can be localized, that each data has its uncertainty and that
! uncertainty has its data

subroutine check_obs (nsample, sample, obs, field_name, nfield, xlon1, xlon2, xlat1, xlat2, iproc)

integer :: nsample, nfield
double precision, dimension(nfield,nsample) :: obs
double precision, intent(in) :: xlon1, xlon2, xlat1, xlat2

integer :: jsample, i, it
character*6, dimension(nfield) :: field_name
character*100, dimension(nsample) :: sample

! remove samples that are not in the Pecube box defined by xlon1, xlon2, xlat1, xlat2
1111 continue
  do jsample = 1, nsample
    if ((obs(1,jsample)-xlon1)*(obs(1,jsample)-xlon2).gt.0.d0 .or. &
        (obs(2,jsample)-xlat1)*(obs(2,jsample)-xlat2).gt.0.d0) then
    if (iproc.eq.0) write (123,'(a)') 'Warning: sample '//trim(sample(jsample))//' Outside lat-lon box for this Pecube run'
    if (iproc.eq.0) write (123,'(a)') 'Lat-Lon:'
    if (iproc.eq.0) write (123,*) obs(2,jsample), obs(1,jsample)
    if (iproc.eq.0) write (123,'(a)') 'Bounds:'
    if (iproc.eq.0) write (123,*) xlat1, xlat2, xlon1, xlon2
      do isample = jsample, nsample - 1
      obs(:,isample) = obs(:,isample + 1)
      sample(isample) = sample(isample + 1)
      enddo
    nsample = nsample - 1
    goto 1111
    endif
  enddo

! sends warning and error messages

  do jsample = 1, nsample
  if (obs(1,jsample).lt.-9998..and.iproc.eq.0) write (123,'(a)') 'Error: NO LONGITUDE, cannot use sample '&
                                //trim(sample(jsample))//' that is not located'
  if (obs(2,jsample).lt.-9998..and.iproc.eq.0) write (123,'(a)') 'Error: NO LATITUDE, cannot use sample '&
                                //trim(sample(jsample))//' that is not located'
  if (obs(3,jsample).lt.-9998..and.iproc.eq.0) write (123,'(a)') 'Warning: No Height provided, set to zero for sample '&
                                //trim(sample(jsample))
    do i = 4, 19, 2
    if (obs(i,jsample).gt.-9998. .and. obs(i+1,jsample).lt.-9998..and.iproc.eq.0) write (123,'(a)') &
    'Warning: Adding uncertainty for ' &
    //field_name(i)//' in sample '//trim(sample(jsample))
    if (obs(i,jsample).lt.-9998. .and. obs(i+1,jsample).gt.-9998..and.iproc.eq.0) write (123,'(a)') &
    'Error: Uncertainty is specified' &
    //'but not data for '//field_name(i)//' in sample '//trim(sample(jsample))
    enddo
  enddo

! normalizes track distributions
  do jsample = 1, nsample
  if (maxval(obs(20:39,jsample)).gt.0.d0) then
  where (obs(20:39,jsample).lt.-9998.) obs(20:39,jsample)=0.d0
  obs(20:39,jsample)=obs(20:39,jsample)/sum(obs(20:39,jsample))
  endif
  enddo

! set height to 0 and size to 100 if not given
  do jsample = 1, nsample
  if (obs(3,jsample).lt.-9998.) obs(3,jsample) = 0.d0
  if (obs(40,jsample).lt.-9998.) obs(40,jsample) = 100.d0
  enddo

! stop execution if errors detected

  do jsample = 1, nsample
    if (obs(1,jsample).lt.-9998.) then
    print*,'Error: NO LONGITUDE, cannot use sample '//trim(sample(jsample))//' that is not located'
    stop
    endif
    if (obs(2,jsample).lt.-9998.) then
    print*,'Error: NO LATITUDE, cannot use sample '//trim(sample(jsample))//' that is not located'
    stop
    endif
    do i = 4, 19, 2
    if (obs(i,jsample).gt.-9998. .and. obs(i+1,jsample).lt.-9998.) obs(i+1,jsample) = obs(i,jsample)/10.d0
      if (obs(i,jsample).lt.-9998. .and. obs(i+1,jsample).gt.-9998.) then
      print*, 'Error: Uncertainty is specified but data not for '//field_name(i)//' in sample '//trim(sample(jsample))
      stop
      endif
    enddo
  enddo

end subroutine check_obs

!----------------------------

subroutine read_data_files (fnme, t, nx, ny, f, s)

! this routine reads in information from a data file named fnme
! are returns a table of characters, t, of size nx by ny, as well
! as a character arrays of station name (s(ny)) and field names (f(nx))

use read_string_module

implicit none

character*(*) :: fnme
character*100, dimension(:,:), allocatable :: t
character*100, dimension(:), allocatable :: f,s
integer :: nx,ny

character*100 :: name, field, field_name(1000), sample_name(1000)
integer :: locfield(1000), locsample(1000)
integer :: isample, jsample, nsample, nfield, i, j, iname
integer :: iftl, jftl, jftlp, ftld(20)
double precision ftl
character c1*1, c2*2
character*100, dimension(:,:), allocatable :: tp
character*100, dimension(:), allocatable :: fp

open (77, file=trim(fnme), status='old')

call find_string (77, 'SAMPLE', isample, jsample)

nsample = 0
  do j = 1, 1000
  name = read_string (77, isample, jsample + j)
  if (name.eq.'NotFound') goto 111
    if (trim(name).ne.'') then
    nsample = nsample + 1
    sample_name(nsample) = trim(name)
    locsample(nsample) = j
    endif
  enddo
111 continue

nfield = 0
  do i = 1, 1000
  field = read_string (77, isample+i, jsample)
  if (field.eq.'NotFound') goto 110
    if (field.ne.'') then
    nfield = nfield + 1
    field_name(nfield) = trim(field)
    locfield(nfield) = i
    endif
  enddo
110 continue

allocate (fp(nfield), s(nsample))
nx = nfield
ny = nsample
fp(1:nx) = field_name(1:nfield)
s(1:ny) = sample_name(1:nsample)

allocate (tp(nx,ny))
do i = 1, nfield
do j = 1, nsample
tp(i,j) = read_string(77, isample + locfield(i), jsample + locsample(j))
enddo
enddo

! change for fission track length measurements

call find_string (77, 'FTL', iftl, jftl)
  if (iftl.ne.0 .and. jftl.ne.0) then
    do i = 1, 20
    nfield = nfield + 1
      if (i.lt.10) then
      write (c1,'(i1)') i
      field_name(nfield) = 'TL'//c1
      else
      write (c2,'(i2)') i
      field_name(nfield) = 'TL'//c2
      endif
    enddo
  allocate (f(nfield),t(nfield,ny))
  f(1:nfield) = field_name(1:nfield)
  t(1:nx,1:ny) = tp(1:nx,1:ny)
  name=''
  jftlp = jftl
  ftld = 0
    do while (name.ne.'NotFound')
    jftlp = jftlp + 1
    name = read_string (77, isample, jftlp)
      if (trim(name).ne.''.and.trim(name).ne.'NotFound') then
! find the sample number to which the FTL distribution will be added
      iname = 1
        do while (trim(sample_name(iname)).ne.trim(name))
        iname = iname + 1
        enddo
      endif
    field = read_string(77, iftl, jftlp)
      if ((trim(field).eq.'NotFound'.or.trim(field).eq.'').and.sum(ftld).ne.0) then
      write (t(nx+1:nfield,iname),'(i10)') ftld
      ftld = 0
      endif
      if (trim(field).ne.'NotFound'.and.trim(field).ne.'') then
      read (field,*,err=998) ftl
      if (ftl*(ftl-20.d0).gt.0.d0) goto 997
      ftld(1+int(ftl)) = ftld(1+int(ftl)) + 1
      endif
    enddo
  nx = nfield
  else
  allocate (f(nx),t(nx,ny))
  f = fp
  t = tp
  endif

close (77)

return

998 continue
print*,'Error in reading fission track length in file '//trim(fnme)//' at line ',jftlp
stop 'Wrong format'

997 continue
print*,'Fission track length of ',ftl,' in '//trim(fnme)//' at line ',jftlp
print*,'ouside permissible bounds (0-20 micons)'
stop 'Change value'

end subroutine read_data_files

!----------------------------

subroutine read_data_files_for_thermal_histories (folder, cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs)

! this subroutine finds and reads the thermal histories from all files contained in folder
! and writes it to unit 111 (the data input file for Pecube)

use read_string_module

implicit none

double precision, intent(in) :: xlon1, xlon2, xlat1, xlat2
character*(*) :: folder
integer, intent(in) :: nd
character(len=1024) :: fnme
integer :: itemph, jtemph, isample, jsample, itimeh, jtimeh, idtemph, jdtemph
integer :: ilat, jlat, ilon, jlon, iheight, jheight
integer :: nobs, nhist
character*100 :: sample, field
double precision :: lat, lon, height, timeh(10000), temph(10000), dtemph(10000)
integer :: i
character*4 :: cproc

! first count the number of thermal histories in all files

nobs = 0

open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=996)

  do

  read (8,'(a)', end=995) fnme
  open (9, file = folder//'/'//fnme, status = 'old')
  call find_string (9, 'TIMEH', itimeh, jtimeh)
    if (itimeh.ne.0) then
    call find_string (9, 'SAMPLE', isample, jsample)

      if (isample.ne.0) then
      sample=''
        do while (trim(sample).ne.'NotFound')
        jsample = jsample + 1
        sample = read_string(9, isample, jsample)
        if (sample.ne.''.and.sample.ne.'NotFound') then
        call find_string (9, 'LAT', ilat, jlat)
          if (ilat.eq.0) then
          print*,'Need LAT to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
          stop
          endif
        field = read_string(9, ilat, jsample)
        read(field,*) lat
        call find_string (9, 'LON', ilon, jlon)
          if (ilon.eq.0) then
          print*,'Need LON to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
          stop
          endif
        field = read_string(9, ilon, jsample)
        read(field,*) lon
!        if (lon.lt.0.d0) lon = lon + 360.d0
        if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) nobs = nobs + 1
        endif
        enddo
      endif
    endif

  close (9)
  enddo

995 close (8)
996 continue

if (nd.eq.0) print*,'Number of Thermal Histories:',nobs

if (nobs.eq.0) return

open (111, file = folder//cproc//'.txt', access = 'append')

! then reads the lat, lon, height and thermal histories of all samples in all files

open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=994)

  do

  read (8,'(a)', end=993) fnme

  open (9, file = folder//'/'//fnme, status = 'old')
  call find_string (9, 'SAMPLE', isample, jsample)
  call find_string (9, 'TIMEH', itimeh, jtimeh)
  call find_string (9, 'TEMPH', itemph, jtemph)
  call find_string (9, 'DTEMPH', idtemph, jdtemph)
  call find_string (9, 'LAT', ilat, jlat)
    if (ilat.eq.0) then
    print*,'Need LAT to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
    stop
    endif
  field = read_string(9, ilat, jsample+1)
  read(field,*) lat
  call find_string (9, 'LON', ilon, jlon)
    if (ilon.eq.0) then
    print*,'Need LON to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
    stop
    endif
  field = read_string(9, ilon, jsample+1)
  read(field,*) lon
!  if (lon.lt.0.d0) lon = lon + 360.d0

    if (itimeh.ne.0) then
    call find_string (9, 'SAMPLE', isample, jsample)

      if (itemph.eq.0) then
      print*,'Found TIMEH field but not TEMPH field in file ',trim(folder//'/'//fnme)
      stop
      endif

      if (isample.ne.0) then
      sample=''
      nhist=0

        do while (trim(sample).ne.'NotFound')
        jsample = jsample + 1
        sample = read_string(9, isample, jsample)

          if (sample.ne.'') then
            if (nhist.ne.0) then
            if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) &
              write (111,*) lon,lat,height,nhist,(timeh(i),temph(i),dtemph(i),i=1,nhist)
            endif
            if (sample.ne.'NotFound') then
            call find_string (9, 'LAT', ilat, jlat)
              if (ilat.eq.0) then
              print*,'Need LAT to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
              stop
              endif
            field = read_string(9, ilat, jsample)
            read(field,*) lat
            call find_string (9, 'LON', ilon, jlon)
              if (ilon.eq.0) then
              print*,'Need LON to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
              stop
              endif
            field = read_string(9, ilon, jsample)
            read(field,*) lon
!            if (lon.lt.0.d0) lon = lon + 360.d0
            call find_string (9, 'HEIGHT', iheight, jheight)
            height = 0.
              if (iheight.ne.0) then
              field = read_string(9, iheight, jsample)
              read(field,*) height
              endif
            nhist = 0
            endif
          endif

        field = read_string(9, itimeh, jsample)
          if (field.ne.''.and.field.ne.'NotFound') then
          nhist = nhist + 1
          field = read_string(9, itimeh, jsample)
          read (field,*,err=999) timeh(nhist)
          field = read_string(9, itemph, jsample)
          read (field,*,err=998) temph(nhist)
            if (idtemph.ne.0) then
            field = read_string(9, idtemph, jsample)
            read (field,*,err=997) dtemph(nhist)
            else
            dtemph(nhist) = temph(nhist)*0.1d0
            endif
          endif

        enddo

      endif

    endif

  close (9)
  enddo

993 close (8)
994 continue

close (111)

return

999 continue
print*,'Error in reading time value in temperature history in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

998 continue
print*,'Error in reading temperature value in temperature history in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

997 continue
print*,'Error in reading temperature uncertainty in temperature history in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

end subroutine read_data_files_for_thermal_histories

!----------------------------

subroutine read_data_files_for_43He (folder, cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs)

! this subroutine finds and reads the thermal histories from all files contained in folder
! and writes it to unit 111 (the data input file for Pecube)

use read_string_module

implicit none

double precision, intent(in) :: xlon1, xlon2, xlat1, xlat2
character*(*) :: folder
integer, intent(in) :: nd
character(len=1024) :: fnme
integer :: idurh, jdurh, itemph, jtemph, irel, jrel, idrel, jdrel, iarel, jarel, idarel, jdarel
integer :: ilat, jlat, ilon, jlon, iheight, jheight, iage, jage, idage, jdage, isize, jsize, isample, jsample
integer :: nobs, nhist
character*100 :: sample, field
double precision :: lat, lon, height,size,age,dage,temph(10000),durh(10000),rel(10000),drel(10000),arel(10000),darel(10000)
integer :: i
character*4 :: cproc

! first count the number of thermal histories in all files

nobs = 0

open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=996)

  do

  read (8,'(a)', end=995) fnme
  open (9, file = folder//'/'//fnme, status = 'old')
  call find_string (9, 'DUR43', idurh, jdurh)
    if (idurh.ne.0) then
    call find_string (9, 'SAMPLE', isample, jsample)

      if (isample.ne.0) then
      sample=''
        do while (trim(sample).ne.'NotFound')
        jsample = jsample + 1
        sample = read_string(9, isample, jsample)
        if (sample.ne.''.and.sample.ne.'NotFound') then
        call find_string (9, 'LAT', ilat, jlat)
          if (ilat.eq.0) then
          print*,'Need LAT to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
          stop
          endif
        field = read_string(9, ilat, jsample+1)
        read(field,*) lat
        call find_string (9, 'LON', ilon, jlon)
          if (ilon.eq.0) then
          print*,'Need LON to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
          stop
          endif
        field = read_string(9, ilon, jsample+1)
        read(field,*) lon
!        if (lon.lt.0.d0) lon = lon + 360.d0
        if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) nobs = nobs + 1
        endif
        enddo
      endif
    endif

  close (9)
  enddo

995 close (8)
996 continue

if (nd.eq.0) print*,'Number of 43He Age Data:',nobs

if (nobs.eq.0) return
  
open (111, file = folder//cproc//'.txt', access = 'append')

! then reads the lat, lon, height and thermal histories of all samples in all files

open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=994)

  do

  read (8,'(a)', end=993) fnme

  open (9, file = folder//'/'//fnme, status = 'old')
  call find_string (9, 'SAMPLE', isample, jsample)
  call find_string (9, 'DUR43', idurh, jdurh)
  call find_string (9, 'TEMP43', itemph, jtemph)
  call find_string (9, '%rel', irel, jrel)
  call find_string (9, 'd%rel', idrel, jdrel)
  call find_string (9, 'arel', iarel, jarel)
  call find_string (9, 'darel', idarel, jdarel)

    if (idurh.ne.0) then
    call find_string (9, 'SAMPLE', isample, jsample)

      if (itemph.eq.0) then
      print*,'Found DUR43 field but not TEMP43 field in file ',trim(folder//'/'//fnme)
      stop
      endif

      if (isample.ne.0) then
      sample=''
      nhist=0

        do while (trim(sample).ne.'NotFound')
        jsample = jsample + 1
        sample = read_string(9, isample, jsample)

          if (sample.ne.'') then
            if (nhist.ne.0) then
            if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) &
              write (111,*) lon,lat,height,size,age,dage,nhist,(temph(i),durh(i),rel(i),drel(i),arel(i),darel(i),i=1,nhist)
            endif
            if (sample.ne.'NotFound') then
            call find_string (9, 'LAT', ilat, jlat)
              if (ilat.eq.0) then
              print*,'Need LAT to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
              stop
              endif
            call find_string (9, 'LON', ilon, jlon)
              if (ilat.eq.0) then
              print*,'Need LON to locate sample ',sample,' in file ',trim(folder//'/'//fnme)
              stop
              endif
            call find_string (9, 'SIZE', isize, jsize)
              if (isize.eq.0) then
              print*,'Need GRAIN SIZE for 43He data in sample ',sample,' in file ',trim(folder//'/'//fnme)
              stop
              endif
            call find_string (9, 'AGE43', iage, jage)
              if (iage.eq.0) then
              print*,'Need AGE43 for 43He data in sample ',sample,' in file ',trim(folder//'/'//fnme)
              stop
              endif
            call find_string (9, 'HEIGHT', iheight, jheight)
            call find_string (9, 'DAGE43', idage, jdage)
            field = read_string(9, ilat, jsample)
            read(field,*) lat
            field = read_string(9, ilon, jsample)
            read(field,*) lon
!            if (lon.lt.0.d0) lon = lon + 360.d0
            field = read_string(9, isize, jsample)
            read(field,*,err=899) size
            field = read_string(9, iage, jsample)
            read(field,*,err=898) age
            height = 0.
              if (iheight.ne.0) then
              field = read_string(9, iheight, jsample)
              read(field,*) height
              endif
            dage = age*0.1d0
              if (idage.ne.0) then
              field = read_string(9, idage, jsample)
              read(field,*,err=897) dage
              endif
            nhist = 0
            endif
          endif

        field = read_string(9, idurh, jsample)
          if (field.ne.''.and.field.ne.'NotFound') then
          nhist = nhist + 1
          field = read_string(9, idurh, jsample)
          read (field,*,err=896) durh(nhist)
          field = read_string(9, itemph, jsample)
          read (field,*,err=895) temph(nhist)
          field = read_string(9, irel, jsample)
          read (field,*,err=894) rel(nhist)
            if (idrel.ne.0) then
            field = read_string(9, idrel, jsample)
            read (field,*,err=893) drel(nhist)
            else
            drel(nhist) = rel(nhist)*0.1d0
            endif
          field = read_string(9, iarel, jsample)
          read (field,*,err=892) arel(nhist)
            if (idarel.ne.0) then
            field = read_string(9, idarel, jsample)
            read (field,*,err=891) darel(nhist)
            else
            darel(nhist) = arel(nhist)*0.1d0
            endif
          endif

        enddo

      endif

    endif

  close (9)
  enddo

993 close (8)
994 continue

close (111)

return

899 continue
print*,'Error in reading grain size in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

898 continue
print*,'Error in reading 43He age in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

897 continue
print*,'Error in reading uncertainty in 43He age in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

896 continue
print*,'Error in reading 43He heating step duration in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

895 continue
print*,'Error in reading 43He heating step temperature in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

894 continue
print*,'Error in reading percent argon released in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

893 continue
print*,'Error in reading uncertainty in percent argon released in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

892 continue
print*,'Error in reading age argon released in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

891 continue
print*,'Error in reading uncertainty in age argon released in 43He data in file '//trim(folder//'/'//fnme)//' at line ',jsample
stop 'Wrong format'

end subroutine read_data_files_for_43He

!----------------------------

subroutine read_data_files_for_TSL (folder, cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs)

! this subroutine finds and reads the TSL data from all files contained in folder
! and writes it to unit 111 (the data input file for Pecube)

use read_string_module

implicit none

double precision, intent(in) :: xlon1, xlon2, xlat1, xlat2
character*(*) :: folder
integer, intent(in) :: nd
character(len=1024) :: fnme
integer :: isample
integer :: nobs
double precision :: lat, lon, height, doser, d0, a, et, logs, b, logrho, nn, dnn
integer :: i,j
character*4 cproc
character*128 tag0(12)
character*128 sample_names(10000)
double precision sample_values(12, 10000)
integer nsamples

! first count the number of thermal histories in all files

nobs = 0

tag0(1)='LON'
tag0(2)='LAT'
tag0(3)='HEIGHT'
tag0(4)='DOSER'
tag0(5)='D0'
tag0(6)='ATL'
tag0(7)='ET'
tag0(8)='LOGS'
tag0(9)='BTL'
tag0(10)='LOGRHO'
tag0(11)='N/N'
tag0(12)='DN/N'

open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=996)

do
  read (8,'(a)', end=995) fnme
  call read_table(folder//'/'//fnme, tag0, 12, sample_names, sample_values, nsamples, 10000)
  do isample = 1, nsamples
    if (abs(sample_values(6,isample)+9999).gt.tiny(0.d0)) then
      do j = 1,12
        if (abs(sample_values(j,isample)+9999).lt.tiny(0.d0) .and. j.ne.3) then
          if (nd.eq.0) print*,'TLSample ',trim(sample_names(isample)),' needs the ',trim(tag0(j)),' field'
          stop
        endif
      enddo
      lon = sample_values(1,isample)
      lat = sample_values(2,isample)
      if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) &
        nobs = nobs + 1
    endif
  enddo
enddo

995 close (8)
996 continue

if (nd.eq.0) print*,'Number of TL Ages Data:', nobs

if (nobs.eq.0) return

open (111, file = folder//cproc//'.txt', access = 'append')

! then reads the lat, lon, height and TSL data of all samples in all files

open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=994)

do
  read (8,'(a)', end=993) fnme
  call read_table(folder//'/'//fnme, tag0, 12, sample_names, sample_values, nsamples, 10000)
  do isample = 1, nsamples
    if (abs(sample_values(6,isample)+9999).gt.tiny(0.d0)) then
      lon = sample_values(1,isample)
      lat = sample_values(2,isample)
      height = sample_values(3,isample)
      if (abs(height+9999.).lt.tiny(height)) height = 0
      doser = sample_values(4,isample)
      d0 = sample_values(5,isample)
      a = sample_values(6,isample)
      et = sample_values(7,isample)
      logs = sample_values(8,isample)
      b = sample_values(9,isample)
      logrho = sample_values(10,isample)
      nn = sample_values(11,isample)
      dnn = sample_values(12,isample)
      if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) &
      write (111,*) lon,lat,height,doser,d0,a,et,logs,b,logrho,nn,dnn
    endif
  enddo
enddo

993 close (8)
994 continue

close (111)

return

end subroutine read_data_files_for_TSL

!----------------------------

subroutine read_data_files_for_OSL (folder, cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs)

  ! this subroutine finds and reads the TSL data from all files contained in folder
  ! and writes it to unit 111 (the data input file for Pecube)
  
  use read_string_module
  
  implicit none
  
  double precision, intent(in) :: xlon1, xlon2, xlat1, xlat2
  character*(*) :: folder
  integer, intent(in) :: nd
  character(len=1024) :: fnme
  integer :: isample
  integer :: nobs
  double precision :: lat, lon, height, doser, d0, et, logs, logrho, eu, nn, dnn
  integer :: i,j
  character*4 cproc
  character*128 tag0(12)
  character*128 sample_names(10000)
  double precision sample_values(12, 10000)
  integer nsamples
  
  ! first count the number of thermal histories in all files
  
  nobs = 0
  
  tag0(1)='LON'
  tag0(2)='LAT'
  tag0(3)='HEIGHT'
  tag0(4)='DOSER'
  tag0(5)='D0'
  tag0(6)='ET'
  tag0(7)='LOGS'
  tag0(8)='LOGRHO'
  tag0(9)='EU'
  tag0(10)='N/N'
  tag0(11)='DN/N'
  
  open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=996)
  
  do
    read (8,'(a)', end=995) fnme
    call read_table(folder//'/'//fnme, tag0, 11, sample_names, sample_values, nsamples, 10000)
    do isample = 1, nsamples
      if (abs(sample_values(9,isample)+9999).gt.tiny(0.d0)) then
        do j = 1,11
          if (abs(sample_values(j,isample)+9999).lt.tiny(0.d0) .and. j.ne.3) then
            if (nd.eq.0) print*,'OSLSample ',trim(sample_names(isample)),' needs the ',trim(tag0(j)),' field'
            stop
          endif
        enddo
        lon = sample_values(1,isample)
        lat = sample_values(2,isample)
        if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) &
          nobs = nobs + 1
      endif
    enddo
  enddo
  
  995 close (8)
  996 continue
  
  if (nd.eq.0) print*,'Number of OSL Ages Data:', nobs
  
  if (nobs.eq.0) return
  
  open (111, file = folder//cproc//'.txt', access = 'append')
  
  ! then reads the lat, lon, height and TSL data of all samples in all files
  
  open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=994)
  
  do
    read (8,'(a)', end=993) fnme
    call read_table(folder//'/'//fnme, tag0, 11, sample_names, sample_values, nsamples, 10000)
    do isample = 1, nsamples
      if (abs(sample_values(9,isample)+9999).gt.tiny(0.d0)) then
        lon = sample_values(1,isample)
        lat = sample_values(2,isample)
        height = sample_values(3,isample)
        if (abs(height+9999.).lt.tiny(height)) height = 0
        doser = sample_values(4,isample)
        d0 = sample_values(5,isample)
        et = sample_values(6,isample)
        logs = sample_values(7,isample)
        logrho = sample_values(8,isample)
        eu = sample_values(9,isample)
        nn = sample_values(10,isample)
        dnn = sample_values(11,isample)
        if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) &
          write (111,*) lon,lat,height,doser,d0,et,logs,logrho,eu,nn,dnn
      endif
    enddo
  enddo
  
  993 close (8)
  994 continue
  
  close (111)
  
  return
  
  end subroutine read_data_files_for_OSL
  
  !----------------------------

  subroutine read_data_files_for_ESR (folder, cproc, xlon1, xlon2, xlat1, xlat2, nd, nobs)

    ! this subroutine finds and reads the TSL data from all files contained in folder
    ! and writes it to unit 111 (the data input file for Pecube)
    
    use read_string_module
    
    implicit none
    
    double precision, intent(in) :: xlon1, xlon2, xlat1, xlat2
    character*(*) :: folder
    integer, intent(in) :: nd
    character(len=1024) :: fnme
    integer :: isample
    integer :: nobs
    double precision :: lat, lon, height, doser, d0, logs, et, sigmaet, nn, dnn
    integer :: i,j
    character*4 cproc
    character*128 tag0(12)
    character*128 sample_names(10000)
    double precision sample_values(12, 10000)
    integer nsamples
    
    ! first count the number of thermal histories in all files
    
    nobs = 0
    
    tag0(1)='LON'
    tag0(2)='LAT'
    tag0(3)='HEIGHT'
    tag0(4)='DOSER'
    tag0(5)='D0'
    tag0(6)='LOGS'
    tag0(7)='ET'
    tag0(8)='SIGMAET'
    tag0(9)='N/N'
    tag0(10)='DN/N'
    
    open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=996)
    
    do
      read (8,'(a)', end=995) fnme
      call read_table(folder//'/'//fnme, tag0, 10, sample_names, sample_values, nsamples, 10000)
      do isample = 1, nsamples
        if (abs(sample_values(8,isample)+9999).gt.tiny(0.d0)) then
          do j = 1,10
            if (abs(sample_values(j,isample)+9999).lt.tiny(0.d0) .and. j.ne.3) then
              if (nd.eq.0) print*,'ESR Sample ',trim(sample_names(isample)),' needs the ',trim(tag0(j)),' field'
              stop
            endif
          enddo
          lon = sample_values(1,isample)
          lat = sample_values(2,isample)
          if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) &
            nobs = nobs + 1
        endif
      enddo
    enddo
    
    995 close (8)
    996 continue
    
    if (nd.eq.0) print*,'Number of ESR Ages Data:', nobs
    
    if (nobs.eq.0) return
    
    open (111, file = folder//cproc//'.txt', access = 'append')
    
    ! then reads the lat, lon, height and TSL data of all samples in all files
    
    open (8, file='tmp/data_lst'//cproc//'.txt', status='old', err=994)
    
    do
      read (8,'(a)', end=993) fnme
      call read_table(folder//'/'//fnme, tag0, 10, sample_names, sample_values, nsamples, 10000)
      do isample = 1, nsamples
        if (abs(sample_values(8,isample)+9999).gt.tiny(0.d0)) then
          lon = sample_values(1,isample)
          lat = sample_values(2,isample)
          height = sample_values(3,isample)
          if (abs(height+9999.).lt.tiny(height)) height = 0
          doser = sample_values(4,isample)
          d0 = sample_values(5,isample)
          logs = sample_values(6,isample)
          et = sample_values(7,isample)
          sigmaet = sample_values(8,isample)
          nn = sample_values(9,isample)
          dnn = sample_values(10,isample)
          if ((lon - xlon1)*(lon - xlon2).le.0.d0 .and. (lat - xlat1)*(lat - xlat2).le.0.d0) &
            write (111,*) lon,lat,height,doser,d0,logs,et,sigmaet,nn,dnn
        endif
      enddo
    enddo
    
    993 close (8)
    994 continue
    
    close (111)
    
    return
    
    end subroutine read_data_files_for_ESR
    
    !----------------------------
    
    function read_string (unit, istring, jstring) result (out)

! Returns the string located at location istring,jstring in a csv file
! istring is column number and jstring is line or row number

implicit none

integer, intent(in)  :: unit, istring, jstring
character(:), allocatable  :: out

character*10000 line
character*1 del
integer :: i, eof, start, end, comma

out='NotFound'

rewind (unit)
  do  i = 1, jstring
  read (unit, '(a)', iostat = eof) line
  if (eof.ne.0) return
  enddo

del = ','
start = 1

  do i = 1, istring-1
  comma = index(line(start:), del)
  if (comma.eq.0) return
  start = start + comma
  enddo
comma = index(line(start:), del)
out = ''
if (comma.eq.1) return
end  = start + comma - 2
if (end.lt.start) end = len(trim(line))

out = line(start:end)

end function read_string

!-----------------------------

subroutine find_string (unit, word, iword, jword)

! Returns the position (istring, jstring) of a string from a csv file (unit=unit)
! jstring is line (row) number and istring is column number

use read_data_module

implicit none

character*(*) :: word
integer :: iword, jword, unit

character*1 del
character*10000 line
integer eof, pos, comma, start, end, ends

rewind (unit)

del = ','

eof = 0

jword=0

pos = 0

  do while (eof.ne.-1.and.pos.eq.0)
  read (unit, '(a)', iostat = eof) line
  jword = jword + 1
  pos = index (ucase(line), ucase(word))
    if (pos.ne.0) then
    comma = -1
    iword = 0
    end = pos - 1
    start = 1
      do while (comma.ne.0)
      comma = index(line(start:end), del)
      iword = iword + 1
      start = start + comma
      enddo
    ends = index(line(pos:), del)
      if (ends.eq.0) then
      ends = len_trim(line)
      else
      ends = pos + ends - 2
      endif
    if (ucase(line(start:ends)).eq.ucase(word)) return
    endif
  enddo

iword = 0
jword = 0

return

end subroutine find_string

!-----------------------------

function ucase(in) result (out)

! returns the uppercase version of the string in

implicit none

character (*), intent(in)  :: in
character(:), allocatable  :: out

integer                    :: i, j

out = trim(in)
do i = 1, len_trim(out)
j = iachar(out(i:i))
if ((j-97)*(j-122).le.0) out(i:i) = achar(j-32)
end do

end function ucase

!-----------------------------

!program test

!use read_data_module

!character*100, dimension(:,:), allocatable :: t
!character*100, dimension(:), allocatable :: f,s
!integer nx,ny
!integer i,j,k

!call read_data_files ('Test.csv',t,nx,ny,f,s)
!print*, 'Fields: ',(trim(f(k)),' ', k = 1, nx)
!print*, 'samples: ',(trim(s(k)),' ', k = 1, ny)
!  do j = 1, ny
!  print*, (trim(t(i,j)),' ', i = 1, nx)
!  enddo

!end program test

!program test

!call read_data_folder ('tmp', .true.)

!end program test

