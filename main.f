      PROGRAM CAVITY
      USE global_data

c --- Input ---
      URFP = 0.3 ! pressure under-relaxation factor
      NSWPP = 10000 ! max pressure iterations per time step
      SORMAX = 0.0001 ! poisson convergence criteria

      ! WRITE(6,*) 'MAXIT, RE, N=?'
      ! READ(5,*) MAXIT,  RE, NI
      
      ! WRITE(*,*) '(uds/quick), (fs/simple)=?'
      ! READ(*,*) SCHEME_NAME, SOLVER_NAME
      ! SOLVER_NAME = TRIM(SOLVER_NAME)
      ! SCHEME_NAME = TRIM(SCHEME_NAME)
      SOLVER_NAME='fs'
      SCHEME_NAME='quick'
      RE=1000
      MAXIT=999999
      SELECT CASE (SOLVER_NAME)
        CASE ("fs")
            ! WRITE(*,*) 'Enter time scheme (1 = Euler, 2 = Adams-Bashforth), DT:'
            ! READ(*,*) TIME_SCHEME, DT
            
            TIME_SCHEME=2
            WRITE(*,*) 'NTRANS, DT=?'
            READ(*,*)   NTRANS,  DT
            NI=64
        CASE ("simple")  
            URFU=0.5
            URFV=URFU
            NSWPU=4
            NSWPV=4
            NSWPP=8
            write(6,*)'URFU?'
            read( 5,*) URFU
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