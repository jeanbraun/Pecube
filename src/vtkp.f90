subroutine vtkp (h,nx,ny,nz,xl,yl,zl,istep)

implicit none

integer nx,ny,nz,istep
double precision h(nx,ny,nz),xl,yl,zl
character*4 cstep
integer i,j,k,iunit

iunit=71

write (cstep,'(i4)') istep
if (istep.lt.10) cstep(1:3)='000'
if (istep.lt.100) cstep(1:2)='00'
if (istep.lt.1000) cstep(1:1)='0'

open(unit=iunit,file='AdvectionLSFRK'//cstep//'.vtk')
write(iunit,'(a)')'# vtk DataFile Version 3.0'
write(iunit,'(a)')'LSF'
write(iunit,'(a)')'ASCII'
write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
write(iunit,'(a7,i10,a6)')'POINTS ',nx*ny*nz,' float'

do k=1,nz
  do j=1,ny
    do i=1,nx
    write(iunit,'(3f16.11)') xl*float(i-1)/(nx-1),yl*float(j-1)/(ny-1),zl*float(k-1)/(nz-1)   
    enddo
  enddo
enddo

  write(iunit,'(A6, 2I10)') 'CELLS ',(nx-1)*(ny-1)*(nz-1),9*(nx-1)*(ny-1)*(nz-1)
  do k=1,nz-1
    do j=1,ny-1
      do i=1,nx-1
      write(iunit,'(9I10)') 8,(k-1)*nx*ny+(j-1)*nx+(i-1),(k-1)*nx*ny+(j-1)*nx+i, &
                              (k-1)*nx*ny+j*nx+i,(k-1)*nx*ny+j*nx+(i-1), & 
                              k*nx*ny+(j-1)*nx+(i-1),k*nx*ny+(j-1)*nx+i, &
                              k*nx*ny+j*nx+i,k*nx*ny+j*nx+(i-1)
      enddo
    enddo
  enddo
  write(iunit,'(A11, I10)') 'CELL_TYPES ',(nx-1)*(ny-1)*(nz-1)
    do k=1,(nx-1)*(ny-1)*(nz-1)
    write(iunit,'(I2)') 12
    enddo

write(iunit,'(a11,i10)')'POINT_DATA ',nx*ny*nz

write(iunit,'(a)')'SCALARS H float 1'
write(iunit,'(a)')'LOOKUP_TABLE default'
do k=1,nz
  do j=1,ny
    do i=1,nx
    write(iunit,'(f18.13)') h(i,j,k)
    enddo
  enddo
enddo

close (iunit)

return
end subroutine vtkp
