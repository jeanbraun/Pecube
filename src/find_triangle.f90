subroutine neighbours (icon,neigh,nt)

use Pecube

implicit none

integer nt,np,nedge,nnmax,iedge,jedge,i,j,k,ie,inode,kp,kpp,inodep,inodepp
integer t1,t2
integer icon(3,nt),neigh(3,nt)
type (edge),dimension(:),allocatable::ed
integer,dimension(:,:),allocatable::nodenodenumber,nodeedgenumber
integer,dimension(:),allocatable::nedgepernode

np=maxval(icon)
nedge=nt+np-1
nnmax=24

! first build edge information

allocate (ed(nedge))
allocate (nedgepernode(np))
allocate (nodenodenumber(nnmax,np),nodeedgenumber(nnmax,np))

ed%t1=0
ed%t2=0

nedgepernode=0

iedge=0
do ie=1,nt
   do k=1,3
   inode=icon(k,ie)
   kp=mod(k,3)+1
   kpp=mod(k+1,3)+1
   inodep=icon(kp,ie)
   inodepp=icon(kpp,ie)
   do j=1,nedgepernode(inodep)
      if (nodenodenumber(j,inodep).eq.inode) then
         jedge=nodeedgenumber(j,inodep)
         ed(jedge)%m2=inodepp
         ed(jedge)%t2=ie
         goto 111
      endif
   enddo
   iedge=iedge+1
   if (iedge.gt.nedge) stop 'too many edges'
   nedgepernode(inode)=nedgepernode(inode)+1
   if (nedgepernode(inode).gt.nnmax) then
      print*,'node',inode,'has',nedgepernode(inode),'neighbours'
      stop 'too many neighbours'
   endif
   nodenodenumber(nedgepernode(inode),inode)=inodep
   nodeedgenumber(nedgepernode(inode),inode)=iedge
   ed(iedge)%n1=inode
   ed(iedge)%n2=inodep
   ed(iedge)%m1=inodepp
   ed(iedge)%t1=ie
   111 continue
   enddo
enddo

! we know what the nmber of edges should be and we therefore check it

if (iedge.lt.nedge) stop 'too few edges'
if (iedge.gt.nedge) stop 'too many edges'

! now compute the neighbours information

do iedge=1,nedge
t1=ed(iedge)%t1
t2=ed(iedge)%t2
  if (t1.gt.0) then
    do j=1,3
    if (icon(j,t1).eq.ed(iedge)%n1) neigh(j,t1)=t2
    enddo
  endif
  if (t2.gt.0) then
    do j=1,3
    if (icon(j,t2).eq.ed(iedge)%n2) neigh(j,t2)=t1
    enddo
  endif
enddo

deallocate (ed,nedgepernode,nodenodenumber,nodeedgenumber)

return
end subroutine neighbours

!----------------------

subroutine find_triangle(xp,yp,x,y,icon,neigh,nt,np,it)

implicit none

integer np,nt,it,i,j
double precision xp,yp,x1,x2,y1,y2
double precision x(np),y(np)
integer icon(3,nt),neigh(3,nt)

2 continue

do i=1,3
x1=x(icon(i,it))
y1=y(icon(i,it))
j=mod(i,3)+1
x2=x(icon(j,it))
y2=y(icon(j,it))
if ((x2-x1)*(yp-y1)-(y2-y1)*(xp-x1).lt.0.d0) goto 1
enddo
return

1 continue
it=neigh(i,it)
if (it.eq.0) stop 'error in find_triangle'
goto 2

end subroutine find_triangle
