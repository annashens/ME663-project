
      SUBROUTINE NAST2D
      USE global_data
c
      DO IT=1,MAXIT
c
      CALL CALCU
      CALL CALCV
      CALL CALCP
c --- transient snapshot every NTRANS iterations
      NTRANS = 100
      IF (MOD(IT,NTRANS).EQ.0) THEN
         CALL WRITE_SNAPSHOT(IT)
      END IF
c
      IF(MOD(IT,25).EQ.0) THEN
         WRITE(*,*) IT,N,RESORM,RESORU,RESORV
      END IF
c
c check steady state solutions
c
      IF(MAX(RESORU,RESORV).LE.0.0005) GO TO 1000
      END DO
 1000 CONTINUE
      RETURN
      END
