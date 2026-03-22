      SUBROUTINE WRITE_TECPLOT(FNAME, UNIT)
      USE global_data

      CHARACTER*(*) FNAME
      INTEGER UNIT
      INTEGER I, J
      REAL XC, YC

      OPEN(UNIT, FILE=FNAME)

      WRITE(UNIT,*) 'VARIABLES="X","Y","U","V","P"'
      WRITE(UNIT,*) 'ZONE F=POINT, I=',NI,', J=',NJ

      DO J=1,NJ
         DO I=1,NI

            UC(I,J) = 0.5*(U(I,J) + U(I-1,J))
            VC(I,J) = 0.5*(V(I,J) + V(I,J-1))

            XC = 0.5*DX + (I-1)*DX
            YC = 0.5*DY + (J-1)*DY

            WRITE(UNIT,*) XC, YC, UC(I,J), VC(I,J), P(I,J)

         END DO
      END DO

      CLOSE(UNIT)

      RETURN
      END
c
      SUBROUTINE WRITE_SNAPSHOT(IT)
      USE global_data

      CHARACTER*30 FNAME

      WRITE(FNAME,'(A,"_snap_",I6.6,".dat")') TRIM(SCHEME_NAME), IT
      CALL WRITE_TECPLOT(FNAME, 99)

      RETURN
      END
c
      SUBROUTINE WRITE_RESULT
      USE global_data
      CHARACTER*30 FNAME
      WRITE(FNAME,'(A,"_result.dat")') TRIM(SCHEME_NAME)
      CALL WRITE_TECPLOT(FNAME, 50)
            
      RETURN
      END