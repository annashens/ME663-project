      SUBROUTINE SOR
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

         IF (RESORM .LE. SORMAX) EXIT

      END DO

      RETURN
      END