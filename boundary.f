      SUBROUTINE MODU
      USE global_data
      REAL    :: ULID
      ULID = 1.0
      DO I = 1, NI
         ! South boundary 
         U(I,0) = -U(I,1)   ! no-slip
         ! North boundary
         U(I,NJ+1) = 2.0 * ULID - U(I,NJ)
         IF (SCHEME_ID .EQ. 2) THEN ! quick
            U(I,-1)   = -U(I,2)
            U(I,NJ+2) = 2.0 * ULID - U(I,NJ-1)
         END IF
      END DO 
      DO J = 1, NJ
         ! West boundary
         U(0,J) = 0.0   ! no-slip
         ! East boundary
         U(NI,J) = 0.0
         IF (SCHEME_ID .EQ. 2) THEN ! quick
            U(-1,J)   = 0.0
            U(NI+1,J) = 0.0
         END IF
      END DO

      RETURN
      END
c
      SUBROUTINE MODV
      USE global_data
      DO I = 1, NI
         ! South boundary 
         V(I,0) = 0.0
         ! North boundary
         V(I,NJ) = 0.0

         IF (SCHEME_ID .EQ. 2) THEN 
            V(I,-1)   = 0.0
            V(I,NJ+1) = 0.0
         END IF
      END DO 

      DO J = 1, NJ
         ! West boundary
         V(0,J) = -V(1,J) ! no-slip
         ! East boundary
         V(NI+1,J) = -V(NI,J)

         IF (SCHEME_ID .EQ. 2) THEN 
            V(-1,J)   = -V(2,J)
            V(NI+2,J) = -V(NI-1,J)
         END IF
      END DO

      RETURN
      END
c
      SUBROUTINE MODP
      USE global_data
      DO J = 1, NJ
         AEP(NI,J) = 0.0
         AWP(1,J)  = 0.0
      END DO

      DO I = 1, NI
         ANP(I,NJ) = 0.0
         ASP(I,1)  = 0.0
      END DO

      RETURN
      END