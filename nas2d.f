
      SUBROUTINE NAST2D
      USE global_data
c
      DO IT=1,MAXIT
c
      CALL CALCU
      CALL CALCV
      CALL CALCP
c --- transient snapshot every NTRANS iterations
      NTRANS = 100
      IF (MOD(IT,NTRANS).EQ.0) THEN
         CALL WRITE_SNAPSHOT(IT)
      END IF
c
      IF(MOD(IT,25).EQ.0) THEN
         WRITE(*,*) IT,N,RESORM,RESORU,RESORV
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
      SUBROUTINE CALCU
      USE global_data
c
      DO J=1,NJ
      DO I=1,NI-1
c
      FEU = DY*(0.5*U(I + 1, J) + 0.5*U(I, J))
      FWU = DY*(0.5*U(I - 1, J) + 0.5*U(I, J))
      FNU = DX*(0.5*V(I + 1, J) + 0.5*V(I, J))
      FSU = DX*(0.5*V(I + 1, J - 1) + 0.5*V(I, J - 1))
      AEU(I,J)=(1/RE)*DY/DX + MAX(0.0, -FEU)
      AWU(I,J)=(1/RE)*DY/DX + MAX(0.0, FWU)
      ANU(I,J)=(1/RE)*DX/DY + MAX(0.0, -FNU)
      ASU(I,J)=(1/RE)*DX/DY + MAX(0.0, FSU)
c
      END DO
      END DO
c
      CALL MODU
c
      DO J=1,NJ
      DO I=1,NI-1
      APU(I,J)=AEU(I,J)+AWU(I,J)+ANU(I,J)+ASU(I,J)
        F(I,J)= DT/(DX*DY)*(
     &    (-APU(I,J)+(DX*DY)/DT)*U(I,J)
     &    + AEU(I,J)*U(I + 1, J)
     &    + AWU(I,J)*U(I - 1, J)
     &    + ANU(I,J)*U(I, J + 1)
     &    + ASU(I,J)*U(I, J - 1) 
     &  )     
      END DO
      END DO
c
      RETURN
      END
      SUBROUTINE CALCV
      USE global_data
c
      DO J=1,NJ-1
      DO I=1,NI
      FEV = DY*(0.5*U(I, J + 1) + 0.5*U(I, J))
      FWV = DY*(0.5*U(I - 1, J + 1) + 0.5*U(I - 1, J))
      FNV = DX*(0.5*V(I, J + 1) + 0.5*V(I, J))
      FSV = DX*(0.5*V(I, J - 1) + 0.5*V(I, J))
      AEV(I,J)= (1/RE)*DY/DX + MAX(0.0, -FEV)
      AWV(I,J)= (1/RE)*DY/DX + MAX(0.0, FWV)
      ANV(I,J)= (1/RE)*DX/DY + MAX(0.0, -FNV)
      ASV(I,J)= (1/RE)*DX/DY + MAX(0.0, FSV)
c
      END DO
      END DO
c
      CALL MODV
c
      DO J=1,NJ-1
      DO I=1,NI
      APV(I,J)=AEV(I,J)+AWV(I,J)+ANV(I,J)+ASV(I,J)
        G(I,J)= DT/(DX*DY)*(
     &  (-APV(I,J)+(DX*DY)/DT)*V(I,J)
     &    + AEV(I,J)*V(I + 1, J)
     &    + AWV(I,J)*V(I - 1, J)
     &    + ANV(I,J)*V(I, J + 1)
     &    + ASV(I,J)*V(I, J - 1)
     &  )  
      END DO
      END DO
c
      RETURN
      END
      SUBROUTINE CALCP
      USE global_data
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

