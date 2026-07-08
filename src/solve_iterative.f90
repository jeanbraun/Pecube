      subroutine solve_iterative (n,ndof,ael,f,x,kfix,icon, &
                                  nnode,nelem,niter,proc,ice)

      implicit none

! Gauss-Seidel iterative solver

! note that this routine has only been checked for ndof=1

      integer n,ndof,nnode,nelem,niter
      double precision ael(n*ndof,n*ndof,nelem),f(ndof*nnode),x(ndof*nnode)
      integer icon(n,nelem),kfix(nnode*ndof),ice(nelem),proc(nnode)

      integer,dimension(:),allocatable::jcmax,mask0,request,requestp
      integer,dimension(:,:),allocatable::jc,jcloc,kc,length
      integer,dimension(:,:,:),allocatable::mask
      double precision,dimension(:),allocatable::diag,x0
      double precision,dimension(:,:),allocatable::dx,work
      logical,dimension(:),allocatable::in

      real*4 time1,time2,time3,time4,time5,time6
      integer nproc,iproc
      integer jelem,ielem,ic,melem,inode,je,k,ie,nnodemax,iprc,i,j,jprc,length0
      integer lproc,mproc,jproc
      integer itmax,it,ielemloc,knode,jnode
      double precision tol,beta,timec,sup,sumrefall,sum1,sum0,sumref

      call cpu_time (time5)
      nproc=1
      iproc=0

! prepare inverse connectivity arrays

      allocate (jcmax(nnode),diag(nnode*ndof),x0(nnode*ndof))

      jcmax=0
        do jelem=1,nelem
        ielem=ice(jelem)
          do inode=1,n
          ic=icon(inode,ielem)
          jcmax(ic)=jcmax(ic)+1
          enddo
        enddo

      melem=maxval(jcmax)

      allocate (jc(melem,nnode),jcloc(melem,nnode),kc(melem,nnode))

      jcmax=0
        do jelem=1,nelem
        ielem=ice(jelem)
          do inode=1,n
          ic=icon(inode,ielem)
          jcmax(ic)=jcmax(ic)+1
          jc(jcmax(ic),ic)=ielem
          jcloc(jcmax(ic),ic)=jelem
          kc(jcmax(ic),ic)=inode
          enddo
        enddo

! prepare diagonal

      diag=0.
        do je=1,nelem
        ie=ice(je)
          do k=1,n
          ic=icon(k,ie)
          diag(ic)=diag(ic)+ael(k,k,je)
          enddo
        enddo

! prepare communication masks

      nnodemax=nnode/nproc
      allocate (length(0:nproc-1,0:nproc-1),mask(nnodemax,0:nproc-1,0:nproc-1))
      allocate (mask0(nnode),request(0:nproc-1),in(0:nproc-1),requestp(0:nproc-1))

! proc(i),i=1,nnode is the processor in which node i is calculated

! mask(i,jproc,iproc) is the mask containing which node stored in processor
! jproc is necessary for the computation performed in another processor iproc
! length(jproc,iproc) is the length of the mask

      length=0
        do iprc=0,nproc-1
          do i=1,nnode
            if (proc(i).eq.iprc) then
              do je=1,jcmax(i)
              ie=jc(je,i)
                do j=1,n
                ic=icon(j,ie)
                jprc=proc(ic)
                  if (jprc.ne.iprc) then
                    do k=1,length(jprc,iprc)
                    if (ic.eq.mask(k,jprc,iprc)) goto 4
                    enddo
                  length(jprc,iprc)=length(jprc,iprc)+1
                  if (length(jprc,iprc).gt.nnodemax) then
                  print*,'nnodemax needs to be increased...'
                  stop
                  endif
                  mask(length(jprc,iprc),jprc,iprc)=ic
4                 continue
                  endif
                enddo
              enddo
            endif
          enddo
        enddo

! mask0 is the mask containing the nodes stored in the current processor
! (iproc) but not needed by another processor
! lengtho is the length of mask0

      length0=0
        do i=1,nnode
        if (proc(i).eq.iproc) then
          do jproc=0,nproc-1
            do k=1,length(iproc,jproc)
            if (mask(k,iproc,jproc).eq.i) goto 3
            enddo
          enddo
        length0=length0+1
        mask0(length0)=i
3       continue
        endif
        enddo

! calculates the memory required for message passing and allocates arrays

      lproc=0
      mproc=0
        do jproc=0,nproc-1
        if (length(iproc,jproc).ne.0) lproc=lproc+1
        if (length(jproc,iproc).ne.0) mproc=mproc+1
        enddo

      allocate (dx(nnode*ndof,lproc),work(nnode*ndof,mproc))

! start of iterations

      itmax=100000
      tol=1.e-6
      beta=1.3
      x0=x
      timec=0.d0

      call cpu_time (time3)

      it=0
1     it=it+1
      if (it.gt.itmax) then
      print*,'too many iterations'
      stop
      endif

      call cpu_time (time1)
      call cpu_time (time2)
      timec=timec+time2-time1

! first performs the matrix multiplication for nodes that need to be shared
! and sends the information to the other processors

      lproc=0
        do jproc=0,nproc-1
        if (length(iproc,jproc).ne.0) lproc=lproc+1
        do k=1,length(iproc,jproc)
        inode=mask(k,iproc,jproc)
        sup=f(inode)
          do jelem=1,jcmax(inode)
          ielem=jc(jelem,inode)
          ielemloc=jcloc(jelem,inode)
          knode=kc(jelem,inode)
            do jnode=1,n
            ic=icon(jnode,ielem)
            sup=sup-ael(knode,jnode,ielemloc)*x(ic)
            enddo
          enddo
        if (kfix(inode).ne.1) x(inode)=x0(inode)+beta*sup/diag(inode)
        dx(k,lproc)=x(inode)-x0(inode)
        enddo
          if (length(iproc,jproc).ne.0) then
          call cpu_time (time1)
          call cpu_time (time2)
          timec=timec+time2-time1
          endif
        enddo

! then the remainder

        do k=1,length0
        inode=mask0(k)
        sup=f(inode)
          do jelem=1,jcmax(inode)
          ielem=jc(jelem,inode)
          ielemloc=jcloc(jelem,inode)
          knode=kc(jelem,inode)
            do jnode=1,n
            ic=icon(jnode,ielem)
            sup=sup-ael(knode,jnode,ielemloc)*x(ic)
            enddo
          enddo
        if (kfix(inode).ne.1) x(inode)=x0(inode)+beta*sup/diag(inode)
        enddo

      call cpu_time (time1)
      call cpu_time (time2)
      timec=timec+time2-time1

! check for convergence

      sum1=0.
      sum0=0.
        do inode=1,nnode
        sum1=max(sum1,abs(x(inode)))
        sum0=max(sum0,abs((x(inode)-x0(inode))))
        enddo

      sumref=sum0/sum1
      x0=x

      call cpu_time (time1)
      sumrefall=sumref
      call cpu_time (time2)
      timec=timec+time2-time1

      if(sumrefall.gt.tol) goto 1

      deallocate (jcmax,jc,jcloc,kc,diag,x0,dx,work,mask,length,mask0,request,requestp,in)

      niter=it

      return
      end
