!---------------------------------------

character*1024 function cget (table, label, sample)

use Pecube

implicit none

type (csv), intent(in) :: table
character*(*), intent(in) :: label, sample

integer i, j, jrow

if (table%ncol == 0) stop 'Error in function cget: csv table not initialized'

jrow = 0
do j=1,table%nrow
if (trim(adjustl(sample))==trim(adjustl(table%sample(j)))) jrow = j
end do

if (jrow==0) stop 'Error in function cget: Sample not found'

do i=1, table%ncol
if (table%label(i) == label) then
cget = ''
if (table%field((jrow-1)*table%ncol+i)/='!^') cget = table%field((jrow-1)*table%ncol+i)
return
endif
end do

stop 'Error in function cget: Label not found'

end function cget

!---------------------------------------

double precision function dget (table, label, sample)

use Pecube

implicit none

type (csv) table
character*(*) label, sample
integer i, j, jrow

if (table%ncol == 0) stop 'Error in function dget: csv table not initialized'

jrow = 0
do j=1,table%nrow
if (trim(adjustl(sample))==trim(adjustl(table%sample(j)))) jrow = j
end do

if (jrow==0) stop 'Error in function dget: Sample not found'

do i=1, table%ncol
if (table%label(i) == label) then
dget = 0.d0
if (table%field((jrow-1)*table%ncol+i)/='!^') read (table%field((jrow-1)*table%ncol+i),*) dget
return
endif
end do

stop 'Error in function dget: Label not found'

end function dget

!---------------------------------------

integer function iget (table, label, sample)

use Pecube

implicit none

type (csv) table
character*(*) label, sample
integer i, j, jrow

if (table%ncol == 0) stop 'Error in function iget: csv table not initialized'

jrow = 0
do j=1,table%nrow
if (trim(adjustl(sample))==trim(adjustl(table%sample(j)))) jrow = j
end do

if (jrow==0) stop 'Error in function iget: Sample not found'

do i=1, table%ncol
if (table%label(i) == label) then
iget = 0
if (table%field((jrow-1)*table%ncol+i)/='!^') read (table%field((jrow-1)*table%ncol+i),*) iget
return
endif
end do

stop 'Error in function iget: Label not found'

end function iget

!----------------------------------------------------------

subroutine remove_commas (line)

character*(*) line
integer ico

ico = -1
do while (ico.ne.0)
ico = index (line, ',')
if (ico.ne.0) line(ico:ico) = '.'
end do

end subroutine remove_commas

!----------------------------------------------------------

logical function empty (st, del)

character*(*) st
character*1 del

empty = .true.
if (verify(trim(adjustl(st)), del) == 0) return
empty = .false.

end function empty

!----------------------------------------------------------

subroutine ReadCSV (fnme, table, french)

! subroutine to extract table from CSV file

use Pecube
use MCSV

implicit none

character*(*), intent(in) :: fnme
logical, intent(in) :: french
type (csv), intent(inout) :: table

integer io, nlines, ico, istart, iend, ncol, icol, count, colstation
character*4096 line
logical label_found
character*1 del

del = ','

open (7,file=trim(adjustl(fnme)),status='old')

io = 0

nlines = 0

label_found = .false.

  do while (io == 0)
  read (7, '(a)', iostat = io) line
    if (.not.empty(line, del)) then
    if (french) call remove_commas (line)
!    print*,trim(adjustl(line))
    line=trim(adjustl(line))//del
    ico = -1
    istart = 0
    ncol = 0
      do while (ico.ne.0)
		  ico = index (line(istart+1:), del)
    	iend = istart+ico
      ncol = ncol + 1
      count = count +1
	     	if (istart+1.le.iend-1) then
          if (label_found) then
          table%field(count) = line(istart+1:iend-1)
          if (ncol == colstation) table%sample(nlines) = line(istart+1:iend-1)
!          print*,count,trim(adjustl(table%field(count)))
          else
					table%label(ncol)=line(istart+1:iend-1)
          table%ncol=ncol
          endif
        endif
      istart = iend
      end do
    count = count -1
      if (.not.label_found) then
        do icol =1,ncol
          if (ucase(table%label(icol)).eq.'SAMPLE') then
          nlines = 0
          count = 0
          label_found = .true.
          colstation = icol
!          print*,'Sample',colstation
          end if
        end do
      end if
    nlines = nlines + 1
!    print*,nlines
    endif
  enddo

nlines = nlines - 1

table%nrow = nlines

close (7)

end subroutine ReadCSV

!----------------------------------------------------------

!program Test

!use Pecube
!use MCSV

!implicit none

!type (csv) table
!character*1024 field1
!integer i, j

!call ReadCSV ('Test.csv', table, .false.)

!print*,'--------------'
!print*,table%nrow,table%ncol
!do j=1,table%nrow
!print*,'Sample: '//trim(adjustl(table%sample(j)))
!print*,trim(cget (table, 'SAMPLE', trim(adjustl(table%sample(j)))))
!do i=1,table%ncol
!print*,'Colum:'//trim(adjustl(table%label(i)))
!if (trim(adjustl(table%label(i))).ne.'SAMPLE') print*,dget(table, trim(adjustl(table%label(i))), trim(table%sample(j)))
!end do
!end do

!end program Test

