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

      WRITE(*,*) 'Enter scheme (uds / quick):'
      READ(*,*) SCHEME_NAME
      IF (SCHEME_NAME == 'uds') THEN
        SCHEME_ID = 1
      ELSE IF (SCHEME_NAME == 'quick') THEN
        SCHEME_ID = 2
      ELSE
        WRITE(*,*) 'Unknown scheme'
      STOP
      
      END IF
      WRITE(*,*) 'Selected scheme: ', SCHEME_NAME, SCHEME_ID
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

c
      CALL WRITE_RESULT

      STOP
      END