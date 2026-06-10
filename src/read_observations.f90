!----------------------------

module read_data_module2

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

end module read_data_module2

!----------------------------

module read_string_module2

!----------------------------

interface read_string

!----------------------------

function read_string (unit, istring, jstring) result (out)

integer, intent(in)  :: unit, istring, jstring
character(:), allocatable  :: out

end function read_string

end interface

!----------------------------

end module read_string_module2

!--------------------------------------------------------------------------

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



!--------------------------------------------------------------------------

subroutine find_string (unit, word, iword, jword)

! Returns the position (istring, jstring) of a string from a csv file (unit=unit)
! jstring is line (row) number and istring is column number

use read_data_module2

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