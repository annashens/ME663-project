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
      IF (TRIM(SOLVER_NAME) == 'fs') THEN
            WRITE(FNAME,'(A,A,"_",A,"_result.dat")') TRIM(SOLVER_NAME), TIME_SCHEME, TRIM(SCHEME_NAME)
      ELSE IF (TRIM(SOLVER_NAME) == 'simple') THEN
            WRITE(FNAME,'(A,"_",A,"_result.dat")') TRIM(SOLVER_NAME), TRIM(SCHEME_NAME)
      ELSE
            WRITE(FNAME,'(A)') 'result.dat'
      END IF
      CALL WRITE_TECPLOT(FNAME, 50)
            
      RETURN
      END
c
      SUBROUTINE WRITE_DIAGNOSTICS
      USE global_data
      WRITE(*,*) 'Solver: ', SOLVER_NAME, ', Convection scheme: ', SCHEME_NAME
      SELECT CASE (SOLVER_NAME)
        CASE ('fs')
            WRITE(*,*) 'Selected Time scheme: ', TIME_SCHEME
            WRITE(*,*) 'DT=', DT
            TIME = IT*DT
            WRITE(*,*) 'Reached steady state at Time = ', TIME
        CASE ('simple')
            WRITE(*,*) 'URFU=',URFUs
      END SELECT
      WRITE(*,*) 'NSWPP=', NSWPP
      WRITE(*,*) 'URFP=',URFP
      WRITE(*,*) 'Wall-clock time (s): ', t_elapsed
      WRITE(*,*) 'RE=', RE
      WRITE(*,*) 'N=', NI 
      END