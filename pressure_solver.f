      SUBROUTINE SOR(IT)
      USE global_data

      INTEGER I, J, ITER

      DO ITER=1,NSWPP

         RESORM = 0.0

         DO J=1,NJ
            DO I=1,NI

               P(I,J) =
     &         ( AEP(I,J)*P(I+1,J) + AWP(I,J)*P(I-1,J)
     &         + ANP(I,J)*P(I,J+1) + ASP(I,J)*P(I,J-1)
     &         + BP(I,J) ) * URFP / APP(I,J)
     &         + (1.0-URFP)*P(I,J)

               RESORM = RESORM + ABS(
     &           APP(I,J)*P(I,J)
     &         - AEP(I,J)*P(I+1,J)
     &         - AWP(I,J)*P(I-1,J)
     &         - ANP(I,J)*P(I,J+1)
     &         - ASP(I,J)*P(I,J-1)
     &         - BP(I,J) )

            END DO
         END DO

      IF(RESORM.LE.SORMAX) GO TO 100
      END DO
 100  CONTINUE
c
      RESORU=0.
      DO J=1,NJ
      DO I=1,NI-1
         IF (TIME_SCHEME.EQ.2 .AND. IT.GT.1) THEN ! AB2
            F_eff(I,J)=(3.0/2.0)*F(I,J) - (1.0/2.0)*F_old(I,J)
         ELSE
            F_eff(I,J)=F(I,J)
         ENDIF
      RESORU=RESORU+ABS(
     &U(I,J)-(F_eff(I,J)+DT/DX*(P(I,J)-P(I+1,J)))
     &)
      U(I,J)=F_eff(I,J)+DT/DX*(P(I,J)-P(I+1,J))
      END DO
      END DO
c
      RESORV=0.
      DO J=1,NJ-1
      DO I=1,NI
         IF (TIME_SCHEME.EQ.2 .AND. IT.GT.1) THEN
            G_eff(I,J)=(3.0/2.0)*G(I,J) - (1.0/2.0)*G_old(I,J)
         ELSE
            G_eff(I,J)=G(I,J)
         ENDIF
      RESORV=RESORV+ABS(
     &V(I,J)-(G_eff(I,J)+DT/DY*(P(I,J)-P(I,J+1)))
     &)
      V(I,J)=G_eff(I,J)+DT/DY*(P(I,J)-P(I,J+1))
      END DO
      END DO
c
      RETURN
      END