c     
      SUBROUTINE GRID
      USE global_data
      DX=1./FLOAT(NI)  !uniform grid in x and y
      DY=1./FLOAT(NJ)
c
c including halo data region
c
      Ihalo=2
      Jhalo=2
      DO I=-Ihalo,NI+Ihalo
        DO J=-Jhalo,NJ+Jhalo
            X(I,J)=DX*I
            Y(I,J)=DY*J
      END DO
      END DO
c
      RETURN
      END
c
      SUBROUTINE INIT
      USE global_data
c
      DO I=-Ihalo,NI+Ihalo
      DO J=-Jhalo,NJ+Jhalo
      U(I,J)=0.
      V(I,J)=0.
      P(I,J)=0.
      F(I,J)=0.
      G(I,J)=0.
      END DO
      END DO
      RETURN
      END