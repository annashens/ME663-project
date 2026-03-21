      MODULE global_data
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: NX=100, NY=100

      REAL U(-2:NX,-2:NY), V(-2:NX,-2:NY), P(-2:NX,-2:NY)
      REAL F(-2:NX,-2:NY), G(-2:NX,-2:NY)

      REAL X(-2:NX,-2:NY), Y(-2:NX,-2:NY)

      REAL UC(NX,NY), VC(NX,NY)

      REAL APU(-2:NX,-2:NY), APV(-2:NX,-2:NY), APP(-2:NX,-2:NY)
      REAL AEU(-2:NX,-2:NY), AEV(-2:NX,-2:NY)
      REAL AEEU(-2:NX,-2:NY), AEEV(-2:NX,-2:NY), AEP(-2:NX,-2:NY)
      REAL AWU(-2:NX,-2:NY), AWV(-2:NX,-2:NY)
      REAL AWWU(-2:NX,-2:NY), AWWV(-2:NX,-2:NY), AWP(-2:NX,-2:NY)
      REAL ANU(-2:NX,-2:NY), ANV(-2:NX,-2:NY)
      REAL ANNU(-2:NX,-2:NY), ANNV(-2:NX,-2:NY), ANP(-2:NX,-2:NY)
      REAL ASU(-2:NX,-2:NY), ASV(-2:NX,-2:NY)
      REAL ASSU(-2:NX,-2:NY), ASSV(-2:NX,-2:NY), ASP(-2:NX,-2:NY)

      REAL BU(-2:NX,-2:NY), BV(-2:NX,-2:NY), BP(-2:NX,-2:NY)

      INTEGER NI, NJ, Ihalo, Jhalo, MAXIT, N
      REAL DX, DY, DT
      REAL RE, URFU, URFV, URFP
      INTEGER NSWPU, NSWPV, NSWPP
      REAL RESORU, RESORV, RESORM, SORMAX, TIME

      END MODULE global_data