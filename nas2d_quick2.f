c
c Ref: NAST2D (fractional-step method)
c
      include 'incl.h'
      NI=40
      NJ=40
      RE=100.
      URFP=0.3
      NSWPP=10000
      SORMAX=0.0001
!MAXIT=number of maximum time steps
      write(6,*)'MAXIT,DT,RE=?'
      read( 5,*) MAXIT,DT,RE
c
      CALL GRID !mesh generation for a 1x1 square cavity
      CALL INIT
      CALL CPU_TIME(t_start)
      CALL NAST2D
      CALL CPU_TIME(t_end) 
      t_elapsed = t_end - t_start 
      WRITE(6,*) 'Wall-clock time (s): ', t_elapsed
      WRITE(6,*) 'RE=',RE, ', N=', NI, ', DT=', DT
C
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
      rewind 50
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
      Ihalo=2
      Jhalo=2
      DO I=-Ihalo,NI+Ihalo
      DO J=-Jhalo,NJ+Jhalo
      X(I,J)=DX*I
      Y(I,J)=DY*J
      END DO
      END DO
c
      RETURN
      END
c
      SUBROUTINE INIT
      include 'incl.h'
c
      DO I=-Ihalo,NI+Ihalo
      DO J=-Jhalo,NJ+Jhalo
      U(I,J)=0.
      V(I,J)=0.
      P(I,J)=0.
      F(I,J)=0.
      G(I,J)=0.
      END DO
      END DO
      RETURN
      END
c
      SUBROUTINE NAST2D
      include 'incl.h'
c
      DO IT=1,MAXIT
c
      CALL CALCU
      CALL CALCV
      CALL CALCP
c
      IF(MOD(IT,25).EQ.0) THEN
      WRITE(*,*)IT,N,RESORM,RESORU,RESORV
      END IF
c
c check steady state solutions
c
      IF(MAX(RESORU,RESORV).LE.0.0005) GO TO 1000
      END DO
 1000 CONTINUE
      RETURN
      END
c
      SUBROUTINE BCMOD
      include 'incl.h'
c
      ENTRY MODU
      REAL ULID
      ULID = 1.0
      DO I=1, NI
            ! South boundary 
            U(I,0)=-U(I,1) ! no-slip
            U(I,-1)  = -U(I,2) 
            ! North boundary
            U(I,NJ+1)=2*ULID-U(I,NJ)
            U(I,NJ+2) = 2*ULID - U(I,NJ-1)
      END DO 
      DO J=1, NJ
            ! West boundary
            U(0,J)=0 ! no-slip
            U(-1,J)=0
            ! East boundary
            U(NI,J)=0
            U(NI+1,J)=0
      END DO
      RETURN
c
      ENTRY MODV
      DO I=1, NI
            ! South boundary 
            V(I,0)=0
            V(I,-1)=0
            ! North boundary
            V(I,NJ)=0
            V(I,NJ+1)=0
      END DO 
      DO J=1, NJ
            ! West boundary
            V(0,J)=-V(1,J) ! no-slip
            V(-1,J)=-V(2,J)
            ! East boundary
            V(NI+1,J)=-V(NI,J)
            V(NI+2,J) = -V(NI-1,J)
      END DO
      RETURN
c
      ENTRY MODP
c
      DO J=1,NJ
            AEP(NI,J)=0.
            AWP(1 ,J)=0.
      END DO
c
      DO I=1,NI
            ANP(I,NJ)=0.
            ASP(I,1 )=0.
      END DO
c
      RETURN
      END
      SUBROUTINE CALCU
      include 'incl.h'
c
      DO J=1,NJ
      DO I=1,NI-1
c
      ! FACE FLUX FOR U-CV AT (I,J):
      FEU = DY*(0.5*U(I + 1, J) + 0.5*U(I, J))
      FWU = DY*(0.5*U(I - 1, J) + 0.5*U(I, J))
      FNU = DX*(0.5*V(I + 1, J) + 0.5*V(I, J))
      FSU = DX*(0.5*V(I + 1, J - 1) + 0.5*V(I, J - 1))

      ! A COEFFICIENTS FOR U-CV AT (I,J):
      AEU(I,J) = (1/RE)*DY/DX + 0.75*MAX(0.0, -FEU) - 0.375*MAX(0.0, FEU) + 0.125*MAX(0.0, -FWU)
      AEEU(I,J) = -0.125*MAX(0.0, -FEU)
      AWU(I,J) = (1/RE)*DY/DX + 0.125*MAX(0.0, FEU) - 0.375*MAX(0.0, -FWU) + 0.75*MAX(0.0, FWU)
      AWWU(I,J) = -0.125*MAX(0.0, FWU)
      ANU(I,J) = (1/RE)*DX/DY + 0.75*MAX(0.0, -FNU) - 0.375*MAX(0.0, FNU) + 0.125*MAX(0.0, -FSU)
      ANNU(I,J) = -0.125*MAX(0.0, -FNU)
      ASU(I,J) = (1/RE)*DX/DY + 0.125*MAX(0.0, FNU) - 0.375*MAX(0.0, -FSU) + 0.75*MAX(0.0, FSU)
      ASSU(I,J) = -0.125*MAX(0.0, FSU)
c
      END DO
      END DO
c
      CALL MODU
c
      DO J=1,NJ
      DO I=1,NI-1
      APU(I,J)=AEU(I,J)+AEEU(I,J)+AWU(I,J)+AWWU(I,J)
     &        +ANU(I,J)+ANNU(I,J)+ASU(I,J)+ASSU(I,J)
      F(I,J)= DT/(DX*DY)*(
     &    (-APU(I,J)+(DX*DY)/DT)*U(I,J)
     &    + AEU(I,J)*U(I + 1, J)
     &    + AEEU(I,J)*U(I + 2, J)
     &    + AWU(I,J)*U(I - 1, J)
     &    + AWWU(I,J)*U(I - 2, J)
     &    + ANU(I,J)*U(I, J + 1)
     &    + ANNU(I,J)*U(I, J + 2)
     &    + ASU(I,J)*U(I, J - 1)
     &    + ASSU(I,J)*U(I, J - 2)
     &  )     
      END DO
      END DO
c
      RETURN
      END
      SUBROUTINE CALCV
      include 'incl.h'
c
      DO J=1,NJ-1
      DO I=1,NI
      ! FACE FLUX FOR V-CV AT (I,J):
      FEV = DY*(0.5*U(I, J + 1) + 0.5*U(I, J))
      FWV = DY*(0.5*U(I - 1, J + 1) + 0.5*U(I - 1, J))
      FNV = DX*(0.5*V(I, J + 1) + 0.5*V(I, J))
      FSV = DX*(0.5*V(I, J - 1) + 0.5*V(I, J))
      AEV(I,J) = (1/RE)*DY/DX + 0.75*MAX(0.0, -FEV) - 0.375*MAX(0.0, FEV) + 0.125*MAX(0.0, -FWV)
      AEEV(I,J) = -0.125*MAX(0.0, -FEV)
      AWV(I,J) = (1/RE)*DY/DX + 0.125*MAX(0.0, FEV) - 0.375*MAX(0.0, -FWV) + 0.75*MAX(0.0, FWV)
      AWWV(I,J) = -0.125*MAX(0.0, FWV)
      ANV(I,J) = (1/RE)*DX/DY + 0.75*MAX(0.0, -FNV) - 0.375*MAX(0.0, FNV) + 0.125*MAX(0.0, -FSV)
      ANNV(I,J) = -0.125*MAX(0.0, -FNV)
      ASV(I,J) = (1/RE)*DX/DY + 0.125*MAX(0.0, FNV) - 0.375*MAX(0.0, -FSV) + 0.75*MAX(0.0, FSV)
      ASSV(I,J) = -0.125*MAX(0.0, FSV)
c
      END DO
      END DO
c
      CALL MODV
c
      DO J=1,NJ-1
      DO I=1,NI
      APV(I,J)=AEV(I,J)+AEEV(I,J)+AWV(I,J)+AWWV(I,J)
     &        +ANV(I,J)+ANNV(I,J)+ASV(I,J)+ASSV(I,J)
        G(I,J)= DT/(DX*DY)*(
     &  (-APV(I,J)+(DX*DY)/DT)*V(I,J)
     &    + AEV(I,J)*V(I + 1, J)
     &    + AEEV(I,J)*V(I + 2, J)
     &    + AWV(I,J)*V(I - 1, J)
     &    + AWWV(I,J)*V(I - 2, J)
     &    + ANV(I,J)*V(I, J + 1)
     &    + ANNV(I,J)*V(I, J + 2)
     &    + ASV(I,J)*V(I, J - 1)
     &    + ASSV(I,J)*V(I, J - 2)
     &  )  
      END DO
      END DO
c
      RETURN
      END
      SUBROUTINE CALCP
      include 'incl.h'
c
      DO J=1,NJ
      DO I=1,NI
      AEP(I,J) = (DY/DX)*DT
      AWP(I,J) = (DY/DX)*DT
      ANP(I,J) = (DX/DY)*DT
      ASP(I,J) = (DX/DY)*DT
      END DO
      END DO
c
      CALL MODP
c
      DO J=1,NJ

      DO I=1,NI
      APP(I,J)=AEP(I,J)+AWP(I,J)+ANP(I,J)+ASP(I,J)
       BP(I,J)=
     &-(F(I,J)-F(I-1,J))*DY
     &-(G(I,J)-G(I,J-1))*DX
c
      END DO
      END DO
c
c SOR (successive over relaxation)
c
      DO N=1,NSWPP
      RESORM=0.
      DO J=1,NJ
      DO I=1,NI
      P(I,J)=
     &(AEP(I,J)*P(I+1,J  )+AWP(I,J)*P(I-1,J  )
     &+ANP(I,J)*P(I  ,J+1)+ASP(I,J)*P(I  ,J-1)
     &+BP(I,J))*URFP/APP(I,J)+(1.-URFP)*P(I,J)
c
      RESORM=RESORM+ABS(
     & APP(I,J)*P(I  ,J  )
     &-AEP(I,J)*P(I+1,J  )-AWP(I,J)*P(I-1,J  )
     &-ANP(I,J)*P(I  ,J+1)-ASP(I,J)*P(I  ,J-1)
     &-BP(I,J) )
      END DO
      END DO
      IF(RESORM.LE.SORMAX) GO TO 100
      END DO
 100  CONTINUE
c
      RESORU=0.
      DO J=1,NJ
      DO I=1,NI-1
      RESORU=RESORU+ABS(
     &U(I,J)-(F(I,J)+DT/DX*(P(I,J)-P(I+1,J)))
     &)
      U(I,J)=F(I,J)+DT/DX*(P(I,J)-P(I+1,J))
      END DO
      END DO
c
      RESORV=0.
      DO J=1,NJ-1
      DO I=1,NI
      RESORV=RESORV+ABS(
     &V(I,J)-(G(I,J)+DT/DY*(P(I,J)-P(I,J+1)))
     &)
      V(I,J)=G(I,J)+DT/DY*(P(I,J)-P(I,J+1))
      END DO
      END DO
c
      RETURN
      END