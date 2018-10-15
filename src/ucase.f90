function ucase(in) result (out)

! function that forces a string to be in upper case

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

!program testucase

!use Pecube

!print*,ucase('Jean Braun')

!end program testucase
