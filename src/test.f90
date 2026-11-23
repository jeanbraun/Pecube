program test

use Pecube

implicit none

type (parameters) p
integer nd
real*4 range(2,1024),par(1024)

nd = 0
par = 0.
range = 0.
call read_input_file ('Pecube.in', 1, p, nd, range, par)

print*,p%time_start

end program test
