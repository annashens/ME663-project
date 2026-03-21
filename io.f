      SUBROUTINE WRITE_SNAPSHOT(IT)
      USE global_data
      CHARACTER*30 FNAME
      INTEGER I,J

      WRITE(FNAME,'("snap_",I6.6,".dat")') IT
      OPEN(99,FILE=FNAME)
      WRITE(99,*) 'VARIABLES="X","Y","U","V","P"'
      WRITE(99,*) 'ZONE F=POINT, I=',NI,', J=',NJ

      DO J=1,NJ
      DO I=1,NI
         UC(I,J)=0.5*(U(I,J)+U(I-1,J))
         VC(I,J)=0.5*(V(I,J)+V(I,J-1))
         XC = 0.5*DX + (I-1)*DX
         YC = 0.5*DY + (J-1)*DY
         WRITE(99,*) XC,YC,UC(I,J),VC(I,J),P(I,J)
      END DO 
      END DO

      CLOSE(99)
      RETURN
      END