      PARAMETER(NX=150,NY=150)
      COMMON/VARS/
     &U(-2:NX,-2:NY),V(-2:NX,-2:NY),P(-2:NX,-2:NY)
      COMMON/GEOM/
     &X(-2:NX,-2:NY),Y(-2:NX,-2:NY)
      COMMON/COEF/
     &APU(0:1),APV(0:1),  !0: left, 1: right
     &AEU(0:1),AEV(0:1),AEEU(0:1),AEEV(0:1),
     &AWU(0:1),AWV(0:1),AWWU(0:1),AWWV(0:1),
     &ANU(0:1),ANV(0:1),ANNU(0:1),ANNV(0:1),
     &ASU(0:1),ASV(0:1),ASSU(0:1),ASSV(0:1),
     & BU(0:1), BV(0:1)
      COMMON/PARA_I/
     &NI,NJ,Ihalo,Jhalo,MAXIT
      COMMON/PARA_R/
     &DX,DY,
     &RE,URFU,URFV,URFP,
     &RESORU,RESORV,RESORM,SORMAX
      COMMON/DEBG/
     &RESU(-2:NX,-2:NY),RESV(-2:NX,-2:NY),RESM(-2:NX,-2:NY)
