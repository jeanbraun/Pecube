double precision function hinterpolate (h,x,y,z,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz)

! this is an internal routine to interpolate a 3D field h known
! on a regular grid of dimension nx,ny,nz and extent (xmin-xmax,
! ymin-ymax and zmin-zmax) at point x,y,z by bilinear local interpolation

implicit none

integer nx,ny,nz
double precision h(nx*ny,nz),x,y,z,xmin,xmax,ymin,ymax,zmin,zmax

integer i,j,k
double precision r,s,t

i=max(1,min(nx-1,1+int((x-xmin)*(nx-1)/(xmax-xmin))))
j=max(1,min(ny-1,1+int((y-ymin)*(ny-1)/(ymax-ymin))))
k=max(1,min(nz-1,1+int((z-zmin)*(nz-1)/(zmax-zmin))))

r=max(-1.d0,min(1.d0,((x-xmin)*(nx-1)/(xmax-xmin)-(i-1))*2.d0-1.d0))
s=max(-1.d0,min(1.d0,((y-ymin)*(ny-1)/(ymax-ymin)-(j-1))*2.d0-1.d0))
t=max(-1.d0,min(1.d0,((z-zmin)*(nz-1)/(zmax-zmin)-(k-1))*2.d0-1.d0))

hinterpolate=((1.d0-r)*(1.d0-s)*(1.d0-t)*h(i+(j-1)*nx,k) &
             +(1.d0+r)*(1.d0-s)*(1.d0-t)*h(i+1+(j-1)*nx,k) &
             +(1.d0+r)*(1.d0+s)*(1.d0-t)*h(i+1+j*nx,k) &
             +(1.d0-r)*(1.d0+s)*(1.d0-t)*h(i+j*nx,k) &
             +(1.d0-r)*(1.d0-s)*(1.d0+t)*h(i+(j-1)*nx,k+1) &
             +(1.d0+r)*(1.d0-s)*(1.d0+t)*h(i+1+(j-1)*nx,k+1) &
             +(1.d0+r)*(1.d0+s)*(1.d0+t)*h(i+1+j*nx,k+1) &
             +(1.d0-r)*(1.d0+s)*(1.d0+t)*h(i+j*nx,k+1))/8.d0


end function hinterpolate
