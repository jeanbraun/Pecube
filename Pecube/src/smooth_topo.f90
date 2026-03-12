subroutine smooth_topo (zin, z, nx,ny, nrad, direction)

implicit none

double precision zin(nx,ny), z(nx,ny)
integer direction, nx, ny, nrad, nradp, nnrad
integer i, j, ii, jj, iii, jjj, k
double precision dzmax, weight, weightn, dz
double precision, dimension(:,:), allocatable :: zn

allocate (zn(nx,ny))

z = zin

nradp = 10
nnrad = nrad

do k=1, nnrad

do j=1,ny
    do i=1,nx
        dzmax=-100
        do jj=-nradp,nradp
            jjj = jj+j
            if (jjj.ge.1.and.jjj.le.ny) then
                do ii=-nradp,nradp
                    iii = ii+i
                    if (iii.ge.1.and.iii.le.nx) then
                        dzmax = max(dzmax,abs(z(iii,jjj)-z(i,j)))
                    endif
                enddo
            endif
        enddo
        zn(i,j) = 0
        weightn = 0.d0
        do jj=-nradp,nradp
            jjj = jj+j
            if (jjj.ge.1.and.jjj.le.ny) then
                do ii=-nradp,nradp
                    iii = ii+i
                    if (iii.ge.1.and.iii.le.nx) then
                        weight = (i-iii)**2 + (j-jjj)**2
                        dz = z(iii,jjj) - z(i,j)
                        weight = 10**(dz/dzmax*direction)*exp(-weight**2/(nrad)**2)
                        zn(i,j)= zn(i,j) + weight*z(iii,jjj)
                        weightn = weightn + weight
                    endif
                enddo
            endif
        enddo
        if (weightn.eq.0) print*,weightn,i,j
        zn(i,j) = zn(i,j)/weightn
    enddo
enddo

z = zn

enddo

deallocate (zn)

end subroutine smooth_topo