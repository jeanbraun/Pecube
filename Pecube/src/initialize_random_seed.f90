subroutine initialize_random_seed(proc)

implicit none

integer proc
integer :: i, n, tim
integer, dimension(:), allocatable :: seed

call random_seed (size = n)
allocate (seed(n))
call system_clock (count=tim)
seed = tim + proc + 37 * (/ (i - 1, i = 1, n) /)
call random_seed (put = seed)
deallocate (seed)

return

end subroutine initialize_random_seed
