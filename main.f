      PROGRAM CAVITY
      USE global_data

      REAL t_start, t_end, t_elapsed

c --- Input parameters ---
      NI = 40
      NJ = 40
      RE = 100.0
      URFP = 0.3
      NSWPP = 10000
      SORMAX = 0.0001

      WRITE(6,*) 'MAXIT, DT, RE=?'
      READ(5,*) MAXIT, DT, RE

c --- Setup ---
      CALL GRID
      CALL INIT

c --- Solve ---
      CALL CPU_TIME(t_start)
      CALL NAST2D
      CALL CPU_TIME(t_end)

      t_elapsed = t_end - t_start

c --- Diagnostics ---
      WRITE(6,*) 'Wall-clock time (s): ', t_elapsed
      WRITE(6,*) 'RE=', RE, ', N=', NI, ', DT=', DT
      TIME = IT*DT
      WRITE(*,*) 'Reached steady state at IT =',IT,'  TIME =',TIME

c --- Post-processing ---
      DO I=1,NI
         DO J=1,NJ
            UC(I,J) = 0.5*(U(I,J) + U(I-1,J))
            VC(I,J) = 0.5*(V(I,J) + V(I,J-1))
         END DO
      END DO

c --- Output ---
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