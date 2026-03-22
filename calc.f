      SUBROUTINE CALCU
      USE global_data
c
      DO J=1,NJ
      DO I=1,NI-1
      FEU = DY*(0.5*U(I + 1, J) + 0.5*U(I, J))
      FWU = DY*(0.5*U(I - 1, J) + 0.5*U(I, J))
      FNU = DX*(0.5*V(I + 1, J) + 0.5*V(I, J))
      FSU = DX*(0.5*V(I + 1, J - 1) + 0.5*V(I, J - 1))
      SELECT CASE (SCHEME_ID)
      CASE (1)  ! UDS
         AEU(I,J) = (1.0/RE)*DY/DX + MAX(0.0, -FEU)
         AEEU(I,J) = 0.0
         AWU(I,J) = (1.0/RE)*DY/DX + MAX(0.0, FWU)
         AWWU(I,J) = 0.0
         ANU(I,J) = (1.0/RE)*DX/DY + MAX(0.0, -FNU)
         ANNU(I,J) = 0.0
         ASU(I,J) = (1.0/RE)*DX/DY + MAX(0.0, FSU)
         ASSU(I,J) = 0.0
      CASE (2)  ! QUICK
         AEU(I,J)  = (1.0/RE)*DY/DX + 0.75*MAX(0.0, -FEU)
     &             - 0.375*MAX(0.0, FEU) + 0.125*MAX(0.0, -FWU)
         AEEU(I,J) = -0.125*MAX(0.0, -FEU)
         AWU(I,J)  = (1.0/RE)*DY/DX + 0.125*MAX(0.0, FEU)
     &             - 0.375*MAX(0.0, -FWU) + 0.75*MAX(0.0, FWU)
         AWWU(I,J) = -0.125*MAX(0.0, FWU)
         ANU(I,J)  = (1.0/RE)*DX/DY + 0.75*MAX(0.0, -FNU)
     &             - 0.375*MAX(0.0, FNU) + 0.125*MAX(0.0, -FSU)
         ANNU(I,J) = -0.125*MAX(0.0, -FNU)
         ASU(I,J)  = (1.0/RE)*DX/DY + 0.125*MAX(0.0, FNU)
     &             - 0.375*MAX(0.0, -FSU) + 0.75*MAX(0.0, FSU)
         ASSU(I,J) = -0.125*MAX(0.0, FSU)
      END SELECT
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
     &    + AEU(I,J)*U(I + 1, J) + AEEU(I,J)*U(I + 2, J)
     &    + AWU(I,J)*U(I - 1, J) + AWWU(I,J)*U(I - 2, J)
     &    + ANU(I,J)*U(I, J + 1) + ANNU(I,J)*U(I, J + 2)
     &    + ASU(I,J)*U(I, J - 1) + ASSU(I,J)*U(I, J - 2)
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
      SELECT CASE (SCHEME_ID)
      CASE (1)  ! UDS
            AEV(I,J)= (1/RE)*DY/DX + MAX(0.0, -FEV)
            AEEV(I,J) = 0.0
            AWV(I,J)= (1/RE)*DY/DX + MAX(0.0, FWV)
            AWWV(I,J) = 0.0
            ANV(I,J)= (1/RE)*DX/DY + MAX(0.0, -FNV)
            ANNV(I,J) = 0.0
            ASV(I,J)= (1/RE)*DX/DY + MAX(0.0, FSV)
            ASSV(I,J) = 0.0
      CASE (2) ! QUICK
         AEV(I,J) = (1.0/RE)*DY/DX + 0.75*MAX(0.0, -FEV)
     &           - 0.375*MAX(0.0, FEV) + 0.125*MAX(0.0, -FWV)
         AEEV(I,J) = -0.125*MAX(0.0, -FEV)
         AWV(I,J) = (1.0/RE)*DY/DX + 0.125*MAX(0.0, FEV)
     &          - 0.375*MAX(0.0, -FWV) + 0.75*MAX(0.0, FWV)
         AWWV(I,J) = -0.125*MAX(0.0, FWV)
         ANV(I,J) = (1.0/RE)*DX/DY + 0.75*MAX(0.0, -FNV)
     &          - 0.375*MAX(0.0, FNV) + 0.125*MAX(0.0, -FSV)
         ANNV(I,J) = -0.125*MAX(0.0, -FNV)
         ASV(I,J) = (1.0/RE)*DX/DY + 0.125*MAX(0.0, FNV)
     &          - 0.375*MAX(0.0, -FSV) + 0.75*MAX(0.0, FSV)
         ASSV(I,J) = -0.125*MAX(0.0, FSV)
      END SELECT
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
c
      SUBROUTINE CALCP(IT)
      USE global_data

      DO J=1,NJ
         DO I=1,NI
            AEP(I,J) = (DY/DX)*DT
            AWP(I,J) = (DY/DX)*DT
            ANP(I,J) = (DX/DY)*DT
            ASP(I,J) = (DX/DY)*DT
         END DO
      END DO

      CALL MODP

      DO J=1,NJ
         DO I=1,NI
            APP(I,J)=AEP(I,J)+AWP(I,J)+ANP(I,J)+ASP(I,J)
            BP(I,J)=-(F(I,J)-F(I-1,J))*DY
     &             -(G(I,J)-G(I,J-1))*DX
         END DO
      END DO

      CALL SOR(IT)
      
      RETURN
      END

