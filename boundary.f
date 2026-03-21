      SUBROUTINE BCMOD
      USE global_data
c
      ENTRY MODU
      REAL ULID
      ULID = 1.0
      DO I=1, NI
            ! South boundary 
            U(I,0)=-U(I,1) ! no-slip
            ! North boundary
            U(I,NJ+1)=2*ULID-U(I,NJ)
      END DO 
      DO J=1, NJ
            ! West boundary
            U(0,J)=0 ! no-slip
            ! East boundary
            U(NI,J)=0
      END DO
      RETURN
c
      ENTRY MODV
      DO I=1, NI
            ! South boundary 
            V(I,0)=0
            ! North boundary
            V(I,NJ)=0
      END DO 
      DO J=1, NJ
            ! West boundary
            V(0,J)=-V(1,J) ! no-slip
            ! East boundary
            V(NI+1,J)=-V(NI,J)
      END DO
      RETURN
c
      ENTRY MODP
c
      DO J=1,NJ
            AEP(NI,J)=0.
            AWP(1 ,J)=0.
      END DO
c
      DO I=1,NI
            ANP(I,NJ)=0.
            ASP(I,1 )=0.
      END DO
c
      RETURN
      END