c
c Vanka (1986), JCP, 1986, Vol. 65, p. 138
c
      include 'incl.h'
      DIMENSION UC(NX,NY),VC(NX,NY)
! define # of grid points
      NI=20
      NJ=20
      RE=100.
! define convergence criterion
      SORMAX=0.0001
! input total number of iterations and Reynolds number
      write(6,*)'RE,N,MAXIT=?'
      MAXIT = 10000
      read( 5,*) RE,NI,MAXIT
      NJ=NI
      write(6,*)'URFU?'
      read( 5,*) URFU
      URFV=URFU
! mesh generation for a 1x1 square cavity
      CALL GRID
! specify initial condition
      CALL INIT
! Vanka's SCGS algorithm
      CALL CPU_TIME(t_start)
      CALL SCGS
      CALL CPU_TIME(t_end) 
      t_elapsed = t_end - t_start 
      WRITE(6,*) 'Wall-clock time (s): ', t_elapsed
      WRITE(6,*) 'RE=',RE, ', N=', NI, ', URFU=', URFU
! post-processing
! calculate U,V at the centroid of P-CV for plotting purpose
! 
      DO I=1,NI
      DO J=1,NJ
      UC(I,J)=0.5*(U(I,J)+U(I-1,J))
      VC(I,J)=0.5*(V(I,J)+V(I,J-1))
      END DO
      END DO
c
c tecplot output
c
      open(50,file='result.dat')
      WRITE(50,*) 'RE=',RE, ', MAXIT=', MAXIT, ', URFU=', URFU
      WRITE(50,*) 'Wall-clock time (s): ', t_elapsed
      WRITE(50,*)'VARIABLES = "X", "Y", "U", "V", "P"'
c
      JMAX=NJ
      IMAX=NI
C
      WRITE(50,*)'ZONE F=POINT, I=',IMAX, ', J=',JMAX
      DO J=1,NJ
      DO I=1,NI
      XC=0.5*DX+(I-1)*DX
      YC=0.5*DY+(J-1)*DY
      WRITE(50,*)XC,YC,UC(I,J),VC(I,J),P(I,J)
      END DO
      END DO
      close(50)
c
      STOP
      END
c
      SUBROUTINE GRID
      include 'incl.h'
      DX=1./FLOAT(NI)  !uniform grid in x and y
      DY=1./FLOAT(NJ)

c
c including halo data region
c
! specify number of "halo-data" layers
      Ihalo=2
      Jhalo=2
      DO I=-Ihalo,NI+Ihalo
      DO J=-Jhalo,NJ+Jhalo
      X(I,J)=DX*I
      Y(I,J)=DY*J
      END DO
      END DO
      RETURN
      END
c
      SUBROUTINE INIT
      include 'incl.h'
c
      DO I=-Ihalo,NI+Ihalo
      DO J=-Jhalo,NJ+Jhalo
      U(I,J)=0.0
      V(I,J)=0.0
      P(I,J)=0.0
      END DO
      END DO
      RETURN
      END
c
      SUBROUTINE SCGS
      include 'incl.h'
      REAL FEU(0:1), FWU(0:1), FNU(0:1), FSU(0:1)
      REAL FEV(0:1), FWV(0:1), FNV(0:1), FSV(0:1)

c
      DO IT=1,MAXIT
      RESORU=0.0
      RESORV=0.0
      RESORM=0.0
!block SOR
      DO J=1,NJ
      DO I=1,NI
c
c for u(i  ,j)
c
c Face Flux for U-CV at (i,j):
      ! FACE FLUX FOR U-CV AT (I,J):
      FEU(1) = DY*(0.5*U(I + 1, J) + 0.5*U(I, J))
      FWU(1) = DY*(0.5*U(I - 1, J) + 0.5*U(I, J))
      FNU(1) = DX*(0.5*V(I + 1, J) + 0.5*V(I, J))
      FSU(1) = DX*(0.5*V(I + 1, J - 1) + 0.5*V(I, J - 1))

      ! A COEFFICIENTS FOR U-CV AT (I,J):
      AEU(1) = (1/RE)*DY/DX + 0.75*MAX(0.0, -FEU(1)) - 0.375*MAX(0.0, FEU(1)) + 0.125*MAX(0.0, -FWU(1))
      AEEU(1) = -0.125*MAX(0.0, -FEU(1))
      AWU(1) = (1/RE)*DY/DX + 0.125*MAX(0.0, FEU(1)) - 0.375*MAX(0.0, -FWU(1)) + 0.75*MAX(0.0, FWU(1))
      AWWU(1) = -0.125*MAX(0.0, FWU(1))
      ANU(1) = (1/RE)*DX/DY + 0.75*MAX(0.0, -FNU(1)) - 0.375*MAX(0.0, FNU(1)) + 0.125*MAX(0.0, -FSU(1))
      ANNU(1) = -0.125*MAX(0.0, -FNU(1))
      ASU(1) = (1/RE)*DX/DY + 0.125*MAX(0.0, FNU(1)) - 0.375*MAX(0.0, -FSU(1)) + 0.75*MAX(0.0, FSU(1))
      ASSU(1) = -0.125*MAX(0.0, FSU(1))
      APU(1) = 2*(1/RE)*DX/DY + 2*(1/RE)*DY/DX - 0.375*MAX(0.0, -FEU(1))
     &    + 0.75*MAX(0.0, FEU(1)) - 0.375*MAX(0.0, -FNU(1)) + 0.75*MAX(0.0, FNU(1))
     &    + 0.75*MAX(0.0, -FSU(1)) - 0.375*MAX(0.0, FSU(1)) + 0.75*MAX(0.0, -FWU(1)) - 0.375*MAX(0.0, FWU(1))

      ! B RESIDUAL FOR U-CV AT (I,J):
      BU(1)=DY*(-P(I + 1, J) + P(I, J))
     &    - APU(1)*U(I, J)
     &    + AEU(1)*U(I + 1, J)
     &    + AEEU(1)*U(I + 2, J)
     &    + AWU(1)*U(I - 1, J)
     &    + AWWU(1)*U(I - 2, J)
     &    + ANU(1)*U(I, J + 1)
     &    + ANNU(1)*U(I, J + 2)
     &    + ASU(1)*U(I, J - 1)
     &    + ASSU(1)*U(I, J - 2)
c
c for u(i-1,j)
c
      ! FACE FLUX FOR U-CV AT (I - 1,J):
      FEU(0) = DY*(0.5*U(I - 1, J) + 0.5*U(I, J))
      FWU(0) = DY*(0.5*U(I - 1, J) + 0.5*U(I - 2, J))
      FNU(0) = DX*(0.5*V(I - 1, J) + 0.5*V(I, J))
      FSU(0) = DX*(0.5*V(I - 1, J - 1) + 0.5*V(I, J - 1))

      ! A COEFFICIENTS FOR U-CV AT (I - 1,J):
      AEU(0) = (1/RE)*DY/DX + 0.75*MAX(0.0, -FEU(0)) - 0.375*MAX(0.0, FEU(0)) + 0.125*MAX(0.0, -FWU(0))
      AEEU(0) = -0.125*MAX(0.0, -FEU(0))
      AWU(0) = (1/RE)*DY/DX + 0.125*MAX(0.0, FEU(0)) - 0.375*MAX(0.0, -FWU(0)) + 0.75*MAX(0.0, FWU(0))
      AWWU(0) = -0.125*MAX(0.0, FWU(0))
      ANU(0) = (1/RE)*DX/DY + 0.75*MAX(0.0, -FNU(0)) - 0.375*MAX(0.0, FNU(0)) + 0.125*MAX(0.0, -FSU(0))
      ANNU(0) = -0.125*MAX(0.0, -FNU(0))
      ASU(0) = (1/RE)*DX/DY + 0.125*MAX(0.0, FNU(0)) - 0.375*MAX(0.0, -FSU(0)) + 0.75*MAX(0.0, FSU(0))
      ASSU(0) = -0.125*MAX(0.0, FSU(0))
      APU(0) = 2*(1/RE)*DX/DY + 2*(1/RE)*DY/DX - 0.375*MAX(0.0, -FEU(0))
     &    + 0.75*MAX(0.0, FEU(0)) - 0.375*MAX(0.0, -FNU(0)) + 0.75*MAX(0.0, FNU(0))
     &    + 0.75*MAX(0.0, -FSU(0)) - 0.375*MAX(0.0, FSU(0)) + 0.75*MAX(0.0, -FWU(0)) - 0.375*MAX(0.0, FWU(0))

      ! B RESIDUAL FOR U-CV AT (I - 1,J):
      BU(0)=DY*(P(I - 1, J) - P(I, J))
     &    - APU(0)*U(I - 1, J)
     &    + AEU(0)*U(I, J)
     &    + AEEU(0)*U(I + 1, J)
     &    + AWU(0)*U(I - 2, J)
     &    + AWWU(0)*U(I - 3, J)
     &    + ANU(0)*U(I - 1, J + 1)
     &    + ANNU(0)*U(I - 1, J + 2)
     &    + ASU(0)*U(I - 1, J - 1)
     &    + ASSU(0)*U(I - 1, J - 2)
c
c for v(i,j  )
c
      ! FACE FLUX FOR V-CV AT (I,J):
      FEV(1) = DY*(0.5*U(I, J + 1) + 0.5*U(I, J))
      FWV(1) = DY*(0.5*U(I - 1, J + 1) + 0.5*U(I - 1, J))
      FNV(1) = DX*(0.5*V(I, J + 1) + 0.5*V(I, J))
      FSV(1) = DX*(0.5*V(I, J - 1) + 0.5*V(I, J))

      ! A COEFFICIENTS FOR V-CV AT (I,J):
      AEV(1) = (1/RE)*DY/DX + 0.75*MAX(0.0, -FEV(1)) - 0.375*MAX(0.0, FEV(1)) + 0.125*MAX(0.0, -FWV(1))
      AEEV(1) = -0.125*MAX(0.0, -FEV(1))
      AWV(1) = (1/RE)*DY/DX + 0.125*MAX(0.0, FEV(1)) - 0.375*MAX(0.0, -FWV(1)) + 0.75*MAX(0.0, FWV(1))
      AWWV(1) = -0.125*MAX(0.0, FWV(1))
      ANV(1) = (1/RE)*DX/DY + 0.75*MAX(0.0, -FNV(1)) - 0.375*MAX(0.0, FNV(1)) + 0.125*MAX(0.0, -FSV(1))
      ANNV(1) = -0.125*MAX(0.0, -FNV(1))
      ASV(1) = (1/RE)*DX/DY + 0.125*MAX(0.0, FNV(1)) - 0.375*MAX(0.0, -FSV(1)) + 0.75*MAX(0.0, FSV(1))
      ASSV(1) = -0.125*MAX(0.0, FSV(1))
      APV(1) = 2*(1/RE)*DX/DY + 2*(1/RE)*DY/DX - 0.375*MAX(0.0, -FEV(1))
     &    + 0.75*MAX(0.0, FEV(1)) - 0.375*MAX(0.0, -FNV(1)) + 0.75*MAX(0.0, FNV(1))
     &    + 0.75*MAX(0.0, -FSV(1)) - 0.375*MAX(0.0, FSV(1)) + 0.75*MAX(0.0, -FWV(1)) - 0.375*MAX(0.0, FWV(1))

      ! B RESIDUAL FOR V-CV AT (I,J):
      BV(1)=DX*(-P(I, J + 1) + P(I, J))
     &    - APV(1)*V(I, J)
     &    + AEV(1)*V(I + 1, J)
     &    + AEEV(1)*V(I + 2, J)
     &    + AWV(1)*V(I - 1, J)
     &    + AWWV(1)*V(I - 2, J)
     &    + ANV(1)*V(I, J + 1)
     &    + ANNV(1)*V(I, J + 2)
     &    + ASV(1)*V(I, J - 1)
     &    + ASSV(1)*V(I, J - 2)
c
c for v(i,j-1)
c
       ! FACE FLUX FOR V-CV AT (I,J - 1):
      FEV(0) = DY*(0.5*U(I, J - 1) + 0.5*U(I, J))
      FWV(0) = DY*(0.5*U(I - 1, J - 1) + 0.5*U(I - 1, J))
      FNV(0) = DX*(0.5*V(I, J - 1) + 0.5*V(I, J))
      FSV(0) = DX*(0.5*V(I, J - 1) + 0.5*V(I, J - 2))

      ! A COEFFICIENTS FOR V-CV AT (I,J - 1):
      AEV(0) = (1/RE)*DY/DX + 0.75*MAX(0.0, -FEV(0)) - 0.375*MAX(0.0, FEV(0)) + 0.125*MAX(0.0, -FWV(0))
      AEEV(0) = -0.125*MAX(0.0, -FEV(0))
      AWV(0) = (1/RE)*DY/DX + 0.125*MAX(0.0, FEV(0)) - 0.375*MAX(0.0, -FWV(0)) + 0.75*MAX(0.0, FWV(0))
      AWWV(0) = -0.125*MAX(0.0, FWV(0))
      ANV(0) = (1/RE)*DX/DY + 0.75*MAX(0.0, -FNV(0)) - 0.375*MAX(0.0, FNV(0)) + 0.125*MAX(0.0, -FSV(0))
      ANNV(0) = -0.125*MAX(0.0, -FNV(0))
      ASV(0) = (1/RE)*DX/DY + 0.125*MAX(0.0, FNV(0)) - 0.375*MAX(0.0, -FSV(0)) + 0.75*MAX(0.0, FSV(0))
      ASSV(0) = -0.125*MAX(0.0, FSV(0))
      APV(0) = 2*(1/RE)*DX/DY + 2*(1/RE)*DY/DX - 0.375*MAX(0.0, -FEV(0))
     &    + 0.75*MAX(0.0, FEV(0)) - 0.375*MAX(0.0, -FNV(0)) + 0.75*MAX(0.0, FNV(0))
     &    + 0.75*MAX(0.0, -FSV(0)) - 0.375*MAX(0.0, FSV(0)) + 0.75*MAX(0.0, -FWV(0)) - 0.375*MAX(0.0, FWV(0))

      ! B RESIDUAL FOR V-CV AT (I,J - 1):
      BV(0)=DX*(P(I, J - 1) - P(I, J))
     &    - APV(0)*V(I, J - 1)
     &    + AEV(0)*V(I + 1, J - 1)
     &    + AEEV(0)*V(I + 2, J - 1)
     &    + AWV(0)*V(I - 1, J - 1)
     &    + AWWV(0)*V(I - 2, J - 1)
     &    + ANV(0)*V(I, J)
     &    + ANNV(0)*V(I, J + 1)
     &    + ASV(0)*V(I, J - 2)
     &    + ASSV(0)*V(I, J - 3)
c
      a11=APU(0)
      a15=DY
      a22=APU(1)
      a25=-DY
      a33=APV(0)
      a35=DX
      a44=APV(1)
      a45=-DX
! under-relaxation
      a11=a11/URFU
      a22=a22/URFU
      a33=a33/URFV
      a44=a44/URFV
c
c for p(i,j)
c
      a51=-DY
      a52=DY
      a53=-DX
      a54=DX
c
      b1=BU(0)
      b2=BU(1)
      b3=BV(0)
      b4=BV(1)
      ! continuity residual at (I,J): (u_E - u_W)*DY + (v_N - v_S)*DX
      b5 = -(U(I,J)-U(I-1,J))*DY-(V(I,J)-V(I,J-1))*DX

c
c BCs for u', v'
c
!!!
      IF(I.EQ.1) THEN ! West boundary
      a51=0. ! continuity
      a11=1.
      a15=0.
      b1=0.
      END IF
      IF(I.EQ.NI) THEN ! East boundary
      a52=0. ! continuity
      a22=1.
      a25=0.
      b2=0.
      END IF
      IF(J.EQ.1) THEN ! South boundary 
      a53=0. ! continuity
      a33=1.
      a35=0.
      b3=0.
      END IF
      IF(J.EQ.NJ) THEN ! North boundary
      a54=0. ! continuity
      a44=1.
      a45=0.
      b4=0.
      END IF
!!!
      CALL VANKA(
     &a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,
     &b1,b2,b3,b4,b5,
     &x1,x2,x3,x4,x5)
! correcting U,V and P
      U(I-1,J  )=U(I-1,J  )+x1
      U(I  ,J  )=U(I  ,J  )+x2
      V(I  ,J-1)=V(I  ,J-1)+x3
      V(I  ,J  )=V(I  ,J  )+x4
      P(I  ,J  )=P(I  ,J  )+x5
! calculate residuals at centroid of P-CV
      RESORU=RESORU+(abs(b1)+abs(b2))/2.
      RESORV=RESORV+(abs(b3)+abs(b4))/2.
      RESORM=RESORM+abs(b5)
      END DO
      END DO
      ! WRITE(6,*) 'U at lid:', U(1,NJ), U(NI/2,NJ), U(NI,NJ)

! specify BCs for U and V
      CALL BCMOD
      WRITE(6,51) IT, RESORU, RESORV, RESORM

      IF(MOD(IT,100).EQ.0) THEN
            WRITE(6,51)IT,RESORU,RESORV,RESORM
      END IF
c
      IF(MAX(RESORU,RESORV,RESORM).LE.SORMAX.AND.IT.GT.1) GO TO 1000
      END DO
 1000 CONTINUE
   51 FORMAT(' IT=',I6,'  RESU=',1PE12.4,'  RESV=',1PE12.4,'  RESM=',1PE12.4)




      RETURN
      END

c
      SUBROUTINE VANKA(
     &a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,
     &b1,b2,b3,b4,b5,
     &x1,x2,x3,x4,x5)

      r1=a51/a11
      r2=a52/a22
      r3=a53/a33
      r4=a54/a44
      DEN=r1*a15+r2*a25+r3*a35+r4*a45
      x5=(r1*b1+r2*b2+r3*b3+r4*b4-b5)/DEN
      x1=(b1-a15*x5)/a11
      x2=(b2-a25*x5)/a22
      x3=(b3-a35*x5)/a33
      x4=(b4-a45*x5)/a44
      RETURN
      END
c
      SUBROUTINE BCMOD
      include 'incl.h'
! specify BCs for U and V
!!!
      REAL ULID
      ULID = 1.0

      DO J=1, NJ
            ! West boundary
            V(0,J)=-V(1,J) ! no-slip
            V(-1,J)=-V(2,J)
            ! East boundary
            V(NI+1,J)=-V(NI,J)
            V(NI+2,J) = -V(NI-1,J) 
      END DO

      DO I=1, NI
            ! South boundary 
            U(I,0)=-U(I,1)
            U(I,-1)  = -U(I,2) 
            ! North boundary
            U(I,NJ+1)=2*ULID-U(I,NJ)
            U(I,NJ+2) = 2*ULID - U(I,NJ-1)
      END DO

      RETURN
      END