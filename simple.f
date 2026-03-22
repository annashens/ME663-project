      SUBROUTINE SIMPLE
      USE global_data
c
      DO IT=1,MAXIT
c
!solve u-equation
      CALL SIMPLE_CALCU
!solve v-equation
      CALL SIMPLE_CALCV
!solve p'-equation
      CALL SIMPLE_CALCP
c
      IF(MOD(IT,50).EQ.0) THEN
      WRITE(6,51)IT,RESORU,RESORV,RESORM
      END IF
c
!define convergence criterion
      IF(MAX(RESORU,RESORV,RESORM).LE.SORMAX.AND.IT.GT.1) GO TO 1000
c
      END DO
 1000 CONTINUE
 50   FORMAT(8F7.3)
 51   FORMAT(I5,3F14.5)
      RETURN
      END