      SUBROUTINE SIMPLE_CALCU
      USE global_data
c
      DO J=1,NJ
      DO I=1,NI-1
!define coefficients for u-equation
      FEU = DY*(0.5*U(I + 1, J) + 0.5*U(I, J))
      FWU = DY*(0.5*U(I - 1, J) + 0.5*U(I, J))
      FNU = DX*(0.5*V(I + 1, J) + 0.5*V(I, J))
      FSU = DX*(0.5*V(I + 1, J - 1) + 0.5*V(I, J - 1))
      SELECT CASE (SCHEME_NAME)
      CASE ('uds')  ! UDS
         AEU(I,J) = (1.0/RE)*DY/DX + MAX(0.0, -FEU)
         AEEU(I,J) = 0.0
         AWU(I,J) = (1.0/RE)*DY/DX + MAX(0.0, FWU)
         AWWU(I,J) = 0.0
         ANU(I,J) = (1.0/RE)*DX/DY + MAX(0.0, -FNU)
         ANNU(I,J) = 0.0
         ASU(I,J) = (1.0/RE)*DX/DY + MAX(0.0, FSU)
         ASSU(I,J) = 0.0
      CASE ('quick')  ! QUICK
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
      END DO
      END DO
c
!B.C.
      CALL MODU
c
!calculate residuals
      RESORU=0.
      DO J=1,NJ
      DO I=1,NI-1
      APU(I,J)=AEU(I,J)+AEEU(I,J)+AWU(I,J)+AWWU(I,J)
     &        +ANU(I,J)+ANNU(I,J)+ASU(I,J)+ASSU(I,J)
       BU(I,J)=
     &    + AEU(I,J)*U(I + 1, J) + AEEU(I,J)*U(I + 2, J)
     &    + AWU(I,J)*U(I - 1, J) + AWWU(I,J)*U(I - 2, J)
     &    + ANU(I,J)*U(I, J + 1) + ANNU(I,J)*U(I, J + 2)
     &    + ASU(I,J)*U(I, J - 1) + ASSU(I,J)*U(I, J - 2)
     &    - APU(I,J)*U(I  , J)
     &    + (P(I,J)-P(I+1,J))*DY
      RESORU=RESORU+ABS(BU(I,J))
      END DO
      END DO
c
c SOR
c
!: number of "inner iterations" for u-equation
      DO N=1,NSWPU
      DO J=1,NJ
      DO I=1,NI-1
         U_new =(
     &    AEU(I,J)*U(I + 1, J) + AEEU(I,J)*U(I + 2, J)
     &    + AWU(I,J)*U(I - 1, J) + AWWU(I,J)*U(I - 2, J)
     &    + ANU(I,J)*U(I, J + 1) + ANNU(I,J)*U(I, J + 2)
     &    + ASU(I,J)*U(I, J - 1) + ASSU(I,J)*U(I, J - 2)
     &     + (P(I,J)-P(I+1,J))*DY
     &   ) / APU(I,J)
         U(I,J) = URFU * U_new + (1.0-URFU)*U(I,J)
      END DO
      END DO
      END DO
c
      RETURN
      END
      SUBROUTINE SIMPLE_CALCV
      USE global_data
c
      DO J=1,NJ-1
      DO I=1,NI
      FEV = DY*(0.5*U(I, J + 1) + 0.5*U(I, J))
      FWV = DY*(0.5*U(I - 1, J + 1) + 0.5*U(I - 1, J))
      FNV = DX*(0.5*V(I, J + 1) + 0.5*V(I, J))
      FSV = DX*(0.5*V(I, J - 1) + 0.5*V(I, J))
      SELECT CASE (SCHEME_NAME)
      CASE ('uds')  ! UDS
            AEV(I,J)= (1/RE)*DY/DX + MAX(0.0, -FEV)
            AEEV(I,J) = 0.0
            AWV(I,J)= (1/RE)*DY/DX + MAX(0.0, FWV)
            AWWV(I,J) = 0.0
            ANV(I,J)= (1/RE)*DX/DY + MAX(0.0, -FNV)
            ANNV(I,J) = 0.0
            ASV(I,J)= (1/RE)*DX/DY + MAX(0.0, FSV)
            ASSV(I,J) = 0.0
      CASE ('quick') ! QUICK
         AEV(I,J)  = (1.0/RE)*DY/DX + 0.75*MAX(0.0, -FEV)
     &             - 0.375*MAX(0.0, FEV) + 0.125*MAX(0.0, -FWV)
         AEEV(I,J) = -0.125*MAX(0.0, -FEV)
         AWV(I,J)  = (1.0/RE)*DY/DX + 0.125*MAX(0.0, FEV)
     &             - 0.375*MAX(0.0, -FWV) + 0.75*MAX(0.0, FWV)
         AWWV(I,J) = -0.125*MAX(0.0, FWV)
         ANV(I,J)  = (1.0/RE)*DX/DY + 0.75*MAX(0.0, -FNV)
     &             - 0.375*MAX(0.0, FNV) + 0.125*MAX(0.0, -FSV)
         ANNV(I,J) = -0.125*MAX(0.0, -FNV)
         ASV(I,J)  = (1.0/RE)*DX/DY + 0.125*MAX(0.0, FNV)
     &             - 0.375*MAX(0.0, -FSV) + 0.75*MAX(0.0, FSV)
         ASSV(I,J) = -0.125*MAX(0.0, FSV)
      END SELECT

      END DO
      END DO
c
!B.C.
      CALL MODV
c
!calculate residuals
      RESORV=0.
      DO J=1,NJ-1
      DO I=1,NI
      APV(I,J)= AEV(I,J)+AEEV(I,J)+AWV(I,J)+AWWV(I,J)
     &        + ANV(I,J)+ANNV(I,J)+ASV(I,J)+ASSV(I,J)
       BV(I,J)=
     &    + AEV(I,J)*V(I + 1, J) + AEEV(I,J)*V(I + 2, J)
     &    + AWV(I,J)*V(I - 1, J) + AWWV(I,J)*V(I - 2, J)
     &    + ANV(I,J)*V(I, J + 1) + ANNV(I,J)*V(I, J + 2)
     &    + ASV(I,J)*V(I, J - 1) + ASSV(I,J)*V(I, J - 2)
     &    - APV(I,J)*V(I  ,J  )
     &    + (P(I,J)-P(I,J+1))*DX
      RESORV=RESORV+ABS(BV(I,J))
      END DO
      END DO
c
c SOR
c
!NSWPV: number of "inner iterations" for v-equation
      DO N=1,NSWPV
      DO J=1,NJ-1
      DO I=1,NI
         V_new =(
     &      AEV(I,J)*V(I + 1, J) + AEEV(I,J)*V(I + 2, J)
     &    + AWV(I,J)*V(I - 1, J) + AWWV(I,J)*V(I - 2, J)
     &    + ANV(I,J)*V(I, J + 1) + ANNV(I,J)*V(I, J + 2)
     &    + ASV(I,J)*V(I, J - 1) + ASSV(I,J)*V(I, J - 2)
     &    + (P(I,J)-P(I,J+1))*DX
     &    ) / APV(I,J)
         V(I,J) = URFV * V_new + (1.0-URFV)*V(I,J)
      END DO
      END DO
      END DO
c
      RETURN
      END
      SUBROUTINE SIMPLE_CALCP
      USE global_data
c
      DO J=1,NJ
      DO I=1,NI
        AEP(I,J)= DY**2/APU(I,J)
        AWP(I,J)= DY**2/APU(I-1,J)    
        ANP(I,J)= DX**2/ANV(I,J)
        ASP(I,J)= DX**2/ANV(I,J-1)
!define coefficients for p'-equation
      END DO
      END DO
c
!B.C.
      CALL MODP
c
!calculate residuals
      RESORM=0.
      DO J=1,NJ
      DO I=1,NI
      PP(I,J)=0.
      APP(I,J)=AEP(I,J)+AWP(I,J)+ANP(I,J)+ASP(I,J)
      BP(I,J)=
     &-(U(I,J)-U(I-1,J))*DY
     &-(V(I,J)-V(I,J-1))*DX
      RESORM=RESORM+ABS(BP(I,J))
      END DO
      END DO
c
c SOR
c
!NSWPP: number of "inner iterations" for p'-equation
      DO N=1,NSWPP
      DO J=1,NJ
      DO I=1,NI
        PP(I,J)=(
     &  AEP(I,J)*PP(I+1,J)+AWP(I,J)*PP(I-1,J)
     &  +ANP(I,J)*PP(I,J+1)+ASP(I,J)*PP(I,J-1)
     &  +BP(I,J))/APP(I,J)
      END DO
      END DO
      END DO
c
!correct u-velocity
      DO J=1,NJ
      DO I=1,NI-1
    !   DU = DY / APU(I,J)
      U(I,J) = U(I,J) + DY/APU(I,J) * (PP(I,J) - PP(I+1,J))
      END DO
      END DO
c
!correct v-velocity
      DO J=1,NJ-1
      DO I=1,NI

      V(I,J) = V(I,J) + DX/APV(I,J) * (PP(I,J) - PP(I,J+1))
      END DO
      END DO
c
!correct p (the pressure)
      DO J=1,NJ
      DO I=1,NI
      P(I,J) = P(I,J) + URFP * PP(I,J)
      END DO
      END DO
c
      RETURN
      END