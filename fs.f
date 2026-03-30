
      SUBROUTINE FS
      USE global_data
c
      LOGICAL :: SNAPSHOT_DONE = .FALSE.
      DO IT=1,MAXIT
c
      CALL CALCU(IT)
      CALL CALCV(IT)
      CALL CALCP(IT)
c --- transient snapshot every NTRANS iterations
      IF (MOD(IT,NTRANS).EQ.0 ) THEN
         CALL WRITE_SNAPSHOT(IT)
      !    WRITE(*,*) SNAPSHOT_DONE, DT*IT
      END IF
      IF(.NOT. SNAPSHOT_DONE .AND. ABS(DT*IT - 2.50) < 0.005) THEN
         CALL WRITE_SNAPSHOT(IT)
         SNAPSHOT_DONE = .TRUE.
      !    WRITE(*,*) 'Snapshot done'
      !    WRITE(*,*) DT*IT
      END IF
c`
      IF(MOD(IT,20).EQ.0) THEN
         WRITE(*,*) IT, N, RESORM,RESORU,RESORV
      END IF
      ! --- store old RHS ---
      DO J=1,NJ
      DO I=1,NI-1
            F_old(I,J) = F(I,J)
      END DO
      END DO

      DO J=1,NJ-1
      DO I=1,NI
            G_old(I,J) = G(I,J)
      END DO
      END DO
c
c check steady state solutions
c
      IF(MAX(RESORU,RESORV).LE.0.0005) THEN 
            WRITE(*,*) IT, N, RESORM,RESORU,RESORV
            CALL WRITE_SNAPSHOT(IT)
            GO TO 1000
      END IF
      END DO
      IF(RESORU > 10.0 .OR. RESORV > 10.0) THEN 
            diverged = .TRUE.
            WRITE(*,*) 'Diverged, DT=', DT
            GO TO 1000
      END IF
 1000 CONTINUE
      RETURN
      END
c
