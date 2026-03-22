      PROGRAM CAVITY
      USE global_data
      REAL t_start, t_end, t_elapsed

c --- Input ---
      NI = 40
      NJ = 40
      RE = 100.0
      URFP = 0.3 ! pressure under-relaxation factor
      NSWPP = 10000 ! max pressure iterations per time step
      SORMAX = 0.0001 ! poisson convergence criteria

      WRITE(6,*) 'MAXIT, DT, RE=?'
      READ(5,*) MAXIT, DT, RE

      WRITE(*,*) 'Enter scheme (uds / quick)'
      READ(*,*) SCHEME_NAME, TIME_SCHEME

      IF (TRIM(SCHEME_NAME) == 'uds') THEN
        SCHEME_ID = 1
      ELSE IF (TRIM(SCHEME_NAME) == 'quick') THEN
        SCHEME_ID = 2
      ELSE
        WRITE(*,*) 'Unknown scheme'
      STOP
      END IF

      SOLVER_NAME = 'simple'
      IF (TRIM(SOLVER_NAME) == 'fs') THEN
        SCHEME_ID = 1
      ELSE IF (TRIM(SOLVER_NAME) == 'simple') THEN
        SCHEME_ID = 2
      ELSE
        WRITE(*,*) 'Unknown scheme'
      STOP
      END IF

      SOLVER_ID = 2
      SELECT CASE (SOLVER_ID)
      CASE (1)
        WRITE(*,*) 'Enter time scheme (1 = Euler, 2 = Adams-Bashforth):'
        READ(*,*) TIME_SCHEME
      CASE (2)  
        NI=80
        NJ=80
        RE=100.
        URFU=0.5
        URFV=0.5
        URFP=0.3
        NSWPU=4
        NSWPV=4
        NSWPP=8
        TIME_SCHEME=2
       write(6,*)'URFU?'
        read( 5,*) URFU
        URFV=URFU
      END SELECT
c --- Setup ---
      CALL GRID
      CALL INIT
c --- Solve ---
      CALL CPU_TIME(t_start)
      SELECT CASE (SOLVER_ID)
        CASE (1)
            CALL FS
        CASE (2)
            CALL SIMPLE
      END SELECT
      CALL CPU_TIME(t_end)
      t_elapsed = t_end - t_start
c --- Diagnostics ---
      WRITE(6,*) 'Wall-clock time (s): ', t_elapsed
      WRITE(6,*) 'N=', NI, 'URFP=',URFP,'NSWPP=', NSWPP, 'SORMAX=', SORMAX
      WRITE(*,*) 'Selected scheme: ', SCHEME_NAME, SCHEME_ID, 'Selected time scheme: ', TIME_SCHEME
      WRITE(6,*) 'RE=', RE, ', N=', NI, ', DT=', DT
      TIME = IT*DT
      WRITE(*,*) 'Reached steady state at IT =',IT,'  TIME =',TIME
c --- Post-process ---
      CALL WRITE_RESULT

      STOP
      END