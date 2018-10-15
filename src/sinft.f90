! SINFT

! taken from numerical recipes (see book for further information)

! subroutines called:
! NONE

      SUBROUTINE SINFT(Y,N)
      common /vocal/ ivocal
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      real*8 Y(N)
      THETA=3.14159265358979D0/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      Y(1)=0.0
      M=N/2
      DO 11 J=1,M
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
        Y1=WI*(Y(J+1)+Y(N-J+1))
        Y2=0.5*(Y(J+1)-Y(N-J+1))
        Y(J+1)=Y1+Y2
        Y(N-J+1)=Y1-Y2
11    CONTINUE
      CALL REALFT(Y,M,+1)
      SUM=0.0
      Y(1)=0.5*Y(1)
      Y(2)=0.0
      DO 12 J=1,N-1,2
        SUM=SUM+Y(J)
        Y(J)=Y(J+1)
        Y(J+1)=SUM
12    CONTINUE
      RETURN
      END

