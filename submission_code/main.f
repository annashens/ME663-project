      PROGRAM CAVITY
      USE global_data

c --- Input ---
      NSWPP = 10000 ! max pressure iterations per time step
      SORMAX = 0.0001 ! poisson convergence criteria
      URFP = 0.3 ! pressure under-relaxation factor

      WRITE(*,*) '(uds/quick), (fs/simple)=?'
      READ(*,*) SCHEME_NAME, SOLVER_NAME

      WRITE(6,*) 'MAXIT, RE, N=?'
      READ(5,*) MAXIT,  RE, NI

      SOLVER_NAME = TRIM(SOLVER_NAME)
      SCHEME_NAME = TRIM(SCHEME_NAME)

      SELECT CASE (SOLVER_NAME)
        CASE ("fs")
            WRITE(*,*) 'Enter time scheme (1 = Euler, 2 = Adams-Bashforth), DT:'
            READ(*,*) TIME_SCHEME, DT
            NTRANS=10000
            ! NI=128
        CASE ("simple")  
            ! URFU=0.5
            NSWPU=4
            NSWPV=4
            NSWPP=8
            write(6,*)'URFU=?'
            read( 5,*) URFU
            URFV=URFU
      END SELECT
      NJ = NI
c --- Setup ---
      CALL GRID
      CALL INIT
c --- Solve ---
      CALL CPU_TIME(t_start)
      SELECT CASE (SOLVER_NAME)
        CASE ('fs')
            CALL FS
        CASE ('simple')
            CALL SIMPLE
      END SELECT
      CALL CPU_TIME(t_end)
      t_elapsed = t_end - t_start
c --- Diagnostics ---
      CALL WRITE_DIAGNOSTICS
      
c --- Post-process ---
      CALL WRITE_RESULT

      STOP
      END