**********************************************************************
*     -------------------------
      FUNCTION NEMODSEL(IPAR,E)
*     -------------------------
*     
*     ( purpose )
*        SELECT INTERACTION MODE
*
*     ( input ) 
*       IPAR    : PARICLE CODE
*       E       : NEUTRINO ENERGY ( GEV )
*
*     ( output ) 
*       NEMODESEL
*
*     ( creation date and author )
*       1983.11.26 ; First version by M.Nakahata
*       1985.01.12 ; M.Nakahata   Fermi mothion and pauli principle are
*                                 considered.      
*       1987.08.13 ; N.Sato       include tau
*       1988.10.08 ; T.Kajita     NC-EL. INT & Nuclear effects in 1 pi prod.
*       1988.10.10 ; T.Kajita     IFREE ofr elastic & 1 pai
*                                 ( INO16 is on at apivct for multi pi)      
*       1995.02.25 ; K.Kaneyuki for S.K.
*                    MODSEL -> NEMODSEL
*       1996.01.10 ; K.Kaneyuki
*                    use card
*       1996.02.05 ; Add COHERENT PI-PRODUCTIONS
*                    ( 1994.09.06 T.KAJITA & K.Okumura )
*       1997.12.01 ; Add single eta production by J.Kameda 
*       1998.02.25 ; Add single  K  production by J.Kameda 
*       1998.03.02 ; Add deep inelastic scattering using JET SET by M.Shiozawa
*       1998.03.10 ; Bug fix by J.Kameda
*       2006.06.12 ; Use IBOUND ; Not to use INO16 or IFREE ; Y.Hayato
*       2007.08.23 ; G.Mitsuka add delta->gamma decay
*       2007.11.05 ; G.Mitsuka support target nucleus besides 16O
*       2007.11.10 ; T.Tanaka add upmu mode
*       2013.02.25 ; A.Redij  added npnh mode
*       2014.12.25 ; Y.Hayato added strangeness violated K prod. mode
*                      ( interaction part has not been implemented, though )
*     
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C    MODSEL   : INTERACTION MODE
C
C
C        ######### NEUTRINO MODE #########
C
C            ***** CHARGED CURRENT *****
C
C               -- ELASTIC --
C           1 : NEU,N --> LEPTON-,P
C
C           2 : NEU,N+X --> LEPTON-,P+X  (X=(N or P))
C               -- SINGLE PI FROM DELTA RESONANCE --
C          11 : NEU,P --> LEPTON-,P,PI+
C          12 : NEU,N --> LEPTON-,P,PI0
C          13 : NEU,N --> LEPTON-,N,PI+
C
C          15 : NEU,P --> LEPTON-,P,PI+  ( diffractive )
C
C          16 : NEU,O(16) --> LEPTON-,O(16),PI+
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C          17 : NEU,N --> LEPTON-,P,GAMMA
C
C               -- SINGLE K : Strangeness violated mode --
C                                     (added 14/12/25 Y. Hayato)
C          18 : NEU,N --> LEPTON-,N,K+
C          19 : NEU,N --> LEPTON-,P,K0
C          20 : NEU,N --> LEPTON-,P,K+
C
C               -- MULTI PI (1.3 < W < 2.0 GeV) --
C          21 : NEU,(N OR P) --> LEPTON-,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C          22 : NEU,N --> LEPTON-,P,ETA0
C
C               -- SINGLE K FROM DELTA RESONANCE --
C                                     (added 98/02/25 J.Kameda)
C          23 : NEU,N --> LEPTON-,LAMBDA,K+
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C          26 : NEU,(N OR P) --> LEPTON-,(N OR P),MESONS
C
C            ***** NEUTAL CURRENT *****
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C          31 : NEU,N --> NEU,N,PI0
C          32 : NEU,P --> NEU,P,PI0
C          33 : NEU,N --> NEU,P,PI-
C          34 : NEU,P --> NEU,N,PI+
C
C          15 : NEU,P --> LEPTON-,P,PI0 ( diffractive )
C
C          36 : NEU,O(16) --> NEU,O(16),PI0
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C          38 : NEU,N --> NEU,N,GAMMA
C          39 : NEU,P --> NEU,P,GAMMA
C
C               -- MULTI PI (1.3 GeV < W < 2.0 GeV) --
C          41 : NEU,(N OR P) --> NEU,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C          42 : NEU,N --> NEU,N,ETA0
C          43 : NEU,P --> NEU,P,ETA0
C
C               -- SINGLE  K  FROM DELTA RESONANCE --
C                                     (added 98/02/20 J.Kameda)
C          44 : NEU,N --> NEU,LAMBDA,K0
C          45 : NEU,P --> NEU,LAMBDA,K+
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C          46 : NEU,(N OR P) --> NEU,(N OR P),MESONS
C
C               -- ELASTIC --
C          51 : NEU,P --> NEU,P
C          52 : NEU,N --> NEU,N
C
C
C        ######### ANTI NEUTRINO MODE #########
C
C            ***** CHARGED CURRENT *****
C
C               -- ELASTIC --
C          -1 : NEUBAR,P --> LEPTON+,N
C          -2 : NEU,P+X --> LEPTON-,N+X  (X=(N or P))
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C         -11 : NEUBAR,N --> LEPTON+,N,PI-
C         -12 : NEUBAR,P --> LEPTON+,N,PI0
C         -13 : NEUBAR,P --> LEPTON+,P,PI-
C
C         -15 : NEUBAR,P --> LEPTON+,P,PI-  ( diffractive )
C
C         -16 : NEUBAR,O(16) --> LEPTON+,O(16),PI-
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C         -17 : NEUBAR,P --> LEPTON+,N,GAMMA
C
C               -- MULTI PI (W > 1.4 GEV) --
C         -21 : NEUBAR,(N OR P) --> LEPTON+,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C         -22 : NEUBAR,P --> LEPTON+,N,ETA0
C
C               -- SINGLE  K FROM DELTA RESONANCE --
C                                     (added 98/02/25 J.Kameda)
C         -23 : NEUBAR,P --> LEPTON+,LAMBDA,K0
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C         -26 : NEUBAR,(N OR P) --> LEPTON+,(N OR P),MESONS
C
C               ** NEUTAL CURRENT **
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C         -31 : NEUBAR,N --> NEUBAR,N,PI0
C         -32 : NEUBAR,P --> NEUBAR,P,PI0
C         -33 : NEUBAR,N --> NEUBAR,P,PI-
C         -34 : NEUBAR,P --> NEUBAR,N,PI+
C
C         -35 : NEUBAR,P --> LEPTON+,P,PI0  ( diffractive )
C
C         -36 : NEUBAR,O(16) --> NEUBAR,O(16),PI0
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C         -38 : NEUBAR,N --> NEUBAR,N,GAMMA
C         -39 : NEUBAR,P --> NEUBAR,P,GAMMA
C
C               -- MULTI PI (W > 1.4 GEV) --
C         -41 : NEUBAR,(N OR P) --> NEUBAR,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C         -42 : NEUBAR,N --> NEUBAR,N,ETA0
C         -43 : NEUBAR,P --> NEUBAR,P,ETA0
C
C               -- SINGLE  K  FROM DELTA RESONANCE --
C                                     (added 98/02/20 J.Kameda)
C          -44 : NEUBAR,N --> NEUBAR,LAMBDA,K0
C          -45 : NEUBAR,P --> NEUBAR,LAMBDA,K+
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C         -46 : NEUBAR,(N OR P) --> NEUBAR,(N OR P),MESONS
C
C               -- ELASTIC --
C         -51 : NEUBAR,P --> NEUBAR,P
C         -52 : NEUBAR,N --> NEUBAR,N
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C    IPAR : NEUTRINO TYPE
C        = 12  NUE
C        =-12  NUEBAR
C        = 14  NUMU
C        =-14  NUMUBAR
C    E    : ENERGY OF NEUTRINO (UNIT : GEV)
C
C    MEANING OF DATA
C    RATN   : NEUTRON RATIO IN TARGET
C    RATP   : PROTON RATIO IN TARGET
C    RATFRE : FREE PROTON RATIO IN TARGET
C    IFREE  : =0 BOUND NUCLEON, =1 FREE PROTON
C
C    IBOUND : =0 Free proton  , =1 Bound nucleon
C     
      IMPLICIT NONE

      integer*4 nemodsel
      integer*4 ipar
      real*4    e

#include "necard.h"
#include "vcwork.h"
#include "posinnuc.h"
#include "neutcrs.h"

      real*4   fnelscrs,fnelscrn,fnelspau_select,fnelspan,cohcrsne,fnapicrs,fnmec_select
      external fnelscrs,fnelscrn,fnelspau_select,fnelspan,cohcrsne,fnapicrs,fnmec_select
      real*4   rsspicrs,rssetcrs,rsskcrs,rssgcrs
      external rsspicrs,rssetcrs,rsskcrs,rssgcrs

      real*4   difcrsnerein
      external difcrsnerein

      real*4   rlu
      external rlu

      integer*4 imode
      PARAMETER (IMODE=31)
      REAL*4    BUF(IMODE)

      integer*4 i
      real*4    ran,dum

C      COMMON/FREPRO/IFREE

      INTEGER MODNEU(IMODE-1),MODNEUB(IMODE-1)
      DATA MODNEU /  1, 11, 12, 13, 21, 31, 32, 33, 34, 41,
     &              51, 51, 52, 16, 36, 22, 42, 43, 23, 44,
     &              45,  0, 26, 46, 17, 38, 39, 2,  15, 35/
      DATA MODNEUB/ -1,-11,-12,-13,-21,-31,-32,-33,-34,-41,
     &              -1,-51,-51,-52,-16,-36,-22,-42,-43,-23,
     &              -44,-45,-26,-46,-17,-38,-39,-2,-15,-35/

      Real APICRSCC, APICRSNC, R1320CC, R1320NC

      REAL*4 RATN, RATP, RATFRE, H2O
      RATN  =float(NUMBNDN)/float(NUMBNDN+NUMBNDP+NUMFREP)
      RATP  =float(NUMBNDP+NUMFREP)/float(NUMBNDN+NUMBNDP+NUMFREP)
      RATFRE=float(NUMFREP)/float(NUMBNDN+NUMBNDP+NUMFREP)
      H2O   =float(NUMBNDN+NUMBNDP+NUMFREP)
#ifdef UPMU_RCK
      RATN  =0.5
      RATP  =0.5
      RATFRE =0.0
      H2O = 18.
#endif
C     
C      IFREE=0
      IBOUND=1

      NEMODSEL=0
C     
C --  FOR NEUTRINO
C     
      IF (IPAR.GT.0) THEN
#ifdef UPMU_RCK
         APICRSCC=FNAPICRS(E,IPAR,1,RATP,R1320CC)
         APICRSNC=FNAPICRS(E,IPAR,0,RATN,R1320NC)
#else
         APICRSCC=FNAPICRS(E,IPAR,1,-1.,R1320CC)
         APICRSNC=FNAPICRS(E,IPAR,0,-1.,R1320NC)
#endif
         BUF(1)=FNELSPAU_SELECT(E,IPAR)*RATN*CRSNEUT(1)
         BUF(2)=BUF(1)+RSSPICRS(E,IPAR,1)*RATP*CRSNEUT(2)
         BUF(3)=BUF(2)+RSSPICRS(E,IPAR,2)*RATN*CRSNEUT(3)
         BUF(4)=BUF(3)+RSSPICRS(E,IPAR,3)*RATN*CRSNEUT(4)
C         BUF(5)=BUF(4)+FNAPICRS(E,IPAR,1,-1.,R1320)*CRSNEUT(5) ! old
         BUF(5)=BUF(4)+APICRSCC*R1320CC*CRSNEUT(5)
         BUF(6)=BUF(5)+RSSPICRS(E,IPAR,11)*RATN*CRSNEUT(6)
         BUF(7)=BUF(6)+RSSPICRS(E,IPAR,12)*RATP*CRSNEUT(7)
         BUF(8)=BUF(7)+RSSPICRS(E,IPAR,13)*RATN*CRSNEUT(8)
         BUF(9)=BUF(8)+RSSPICRS(E,IPAR,14)*RATP*CRSNEUT(9)
C         BUF(10)=BUF(9)+FNAPICRS(E,IPAR,0,-1.,R1320)*CRSNEUT(10) ! old
         BUF(10)=BUF(9)+APICRSNC*R1320NC*CRSNEUT(10)
         BUF(11)=BUF(10)+FNELSCRN(E,IPAR,2212)*RATFRE*CRSNEUT(11)
         BUF(12)=BUF(11)+FNELSPAN(E,IPAR,2212)*(RATP-RATFRE)
     $                                               *CRSNEUT(12)
         BUF(13)=BUF(12)+FNELSPAN(E,IPAR,2112)*RATN*CRSNEUT(13)
C                         (=0:NC,   NE.0:CC)   T.KAJITA 06-SEP-94

         BUF(14)=BUF(13)+ COHCRSNE(E,IPAR,1) /H2O *CRSNEUT(14)
         BUF(15)=BUF(14)+ COHCRSNE(E,IPAR,0) /H2O *CRSNEUT(15)

         BUF(16)=BUF(15)+RSSETCRS(E,IPAR,1)*RATN*CRSNEUT(16)
         BUF(17)=BUF(16)+RSSETCRS(E,IPAR,11)*RATN*CRSNEUT(17)
         BUF(18)=BUF(17)+RSSETCRS(E,IPAR,12)*RATP*CRSNEUT(18)

         BUF(19)=BUF(18)+RSSKCRS(E,IPAR,1)*RATN*CRSNEUT(19)
         BUF(20)=BUF(19)+RSSKCRS(E,IPAR,11)*RATN*CRSNEUT(20)
         BUF(21)=BUF(20)+RSSKCRS(E,IPAR,12)*RATP*CRSNEUT(21)

         BUF(22)=BUF(21)+0.*CRSNEUT(22)
         BUF(23)=BUF(22)+APICRSCC*(1.-R1320CC)*CRSNEUT(23)
         BUF(24)=BUF(23)+APICRSNC*(1.-R1320NC)*CRSNEUT(24)

         BUF(25)=BUF(24)+RSSGCRS(E,IPAR, 1)*RATN*CRSNEUT(25)
         BUF(26)=BUF(25)+RSSGCRS(E,IPAR,11)*RATN*CRSNEUT(26)
         BUF(27)=BUF(26)+RSSGCRS(E,IPAR,12)*RATP*CRSNEUT(27)
         BUF(28)=BUF(27)+FNMEC_SELECT(E,IPAR)*CRSNEUT(28)
         BUF(29)=BUF(28)+DIFCRSNEREIN(E,IPAR,1)*RATFRE*CRSNEUT(29)
         BUF(30)=BUF(29)+DIFCRSNEREIN(E,IPAR,0)*RATFRE*CRSNEUT(30)

         TOTCRSNE = BUF(30)

 10      RAN=RLU(DUM)*BUF(30)
         DO 20 I=1,30
            IF (BUF(I).GE.RAN) THEN
               IF (NEMODFLG.NE.0 .AND. NEMODFLG.NE.-1 .AND.
     &             NEMODFLG.NE.ABS(MODNEU(I))) GOTO 10
               NEMODSEL=MODNEU(I)
C              IF (I.EQ.11 .OR. I.EQ.14 .OR. I.EQ.15) IFREE=1
              IF (I.EQ.11 .OR. I.EQ.14 .OR. I.EQ.15 .OR. I.EQ.29
     $              .OR. I.EQ.30 ) IBOUND=0
              IF (I.EQ.2 .OR. I.EQ.7 .OR. I.EQ.9 .OR. I.EQ.18 .OR.
     &         I.EQ.21 .OR. I.EQ.27) THEN
C                  IF (RLU(DUM).LT.(RATFRE/RATP)) IFREE = 1
                  IF (RLU(DUM).LT.(RATFRE/RATP)) IBOUND = 0
               ENDIF
C               WRITE(*,*) "MODE",NEMODSEL,E,IPAR
               RETURN
            ENDIF
 20      CONTINUE
         GO TO 10
C     
C -- FOR ANTI NEUTRINO
C     
      ELSE IF (IPAR.LT.0) THEN
#ifdef UPMU_RCK
         APICRSCC=FNAPICRS(E,IPAR,1,RATP,R1320CC)
         APICRSNC=FNAPICRS(E,IPAR,0,RATN,R1320NC)
#else
         APICRSCC=FNAPICRS(E,IPAR,1,-1.,R1320CC)
         APICRSNC=FNAPICRS(E,IPAR,0,-1.,R1320NC)
#endif
         BUF(1)=FNELSCRS(E,IPAR)*RATFRE*CRSNEUTB(1)
         BUF(2)=BUF(1)+RSSPICRS(E,IPAR,4)*RATN*CRSNEUTB(2)
         BUF(3)=BUF(2)+RSSPICRS(E,IPAR,5)*RATP*CRSNEUTB(3)
         BUF(4)=BUF(3)+RSSPICRS(E,IPAR,6)*RATP*CRSNEUTB(4)
C         BUF(5)=BUF(4)+FNAPICRS(E,IPAR,1,-1.,R1320)*CRSNEUTB(5) ! old
         BUF(5)=BUF(4)+APICRSCC*R1320CC*CRSNEUTB(5)
         BUF(6)=BUF(5)+RSSPICRS(E,IPAR,15)*RATN*CRSNEUTB(6)
         BUF(7)=BUF(6)+RSSPICRS(E,IPAR,16)*RATP*CRSNEUTB(7)
         BUF(8)=BUF(7)+RSSPICRS(E,IPAR,17)*RATN*CRSNEUTB(8)
         BUF(9)=BUF(8)+RSSPICRS(E,IPAR,18)*RATP*CRSNEUTB(9)
C         BUF(10)=BUF(9)+FNAPICRS(E,IPAR,0,-1.,R1320)*CRSNEUTB(10) !old
         BUF(10)=BUF(9)+APICRSNC*R1320NC*CRSNEUTB(10)
         BUF(11)=BUF(10)+FNELSPAU_SELECT(E,IPAR)*(RATP-RATFRE)*CRSNEUTB(11)
         BUF(12)=BUF(11)+FNELSCRN(E,IPAR,2212)*RATFRE*CRSNEUTB(12)
         BUF(13)=BUF(12)+FNELSPAN(E,IPAR,2212)*(RATP-RATFRE)
     $                                                 *CRSNEUTB(13)
         BUF(14)=BUF(13)+FNELSPAN(E,IPAR,2112)*RATN*CRSNEUTB(14)
C                         (=0:NC,   NE.0:CC)   T.KAJITA 06-SEP-94

         BUF(15)=BUF(14)+ COHCRSNE(E,IPAR,1) /H2O *CRSNEUTB(15)
         BUF(16)=BUF(15)+ COHCRSNE(E,IPAR,0) /H2O *CRSNEUTB(16)

         BUF(17)=BUF(16)+RSSETCRS(E,IPAR,2)*RATP*CRSNEUTB(17)
         BUF(18)=BUF(17)+RSSETCRS(E,IPAR,13)*RATN*CRSNEUTB(18)
         BUF(19)=BUF(18)+RSSETCRS(E,IPAR,14)*RATP*CRSNEUTB(19)

         BUF(20)=BUF(19)+RSSKCRS(E,IPAR,2)*RATP*CRSNEUTB(20)
         BUF(21)=BUF(20)+RSSKCRS(E,IPAR,13)*RATN*CRSNEUTB(21)
         BUF(22)=BUF(21)+RSSKCRS(E,IPAR,14)*RATP*CRSNEUTB(22)

         BUF(23)=BUF(22)+APICRSCC*(1.-R1320CC)*CRSNEUTB(23)
         BUF(24)=BUF(23)+APICRSNC*(1.-R1320NC)*CRSNEUTB(24)

         BUF(25)=BUF(24)+RSSGCRS(E,IPAR, 2)*RATP*CRSNEUTB(25)
         BUF(26)=BUF(25)+RSSGCRS(E,IPAR,13)*RATN*CRSNEUTB(26)
         BUF(27)=BUF(26)+RSSGCRS(E,IPAR,14)*RATP*CRSNEUTB(27)
         BUF(28)=BUF(27)+FNMEC_SELECT(E,IPAR)*CRSNEUTB(28)  
         BUF(29)=BUF(28)+DIFCRSNEREIN(E,IPAR,1)*RATFRE*CRSNEUTB(29)
         BUF(30)=BUF(29)+DIFCRSNEREIN(E,IPAR,0)*RATFRE*CRSNEUTB(30)

         TOTCRSNE = BUF(30)

 110     RAN=RLU(DUM)*BUF(30)
         DO 120 I=1,30
            IF (BUF(I).GE.RAN) THEN
               IF (NEMODFLG.NE.0 .AND. NEMODFLG.NE.-1 .AND.
     &             NEMODFLG.NE.ABS(MODNEUB(I))) GOTO 110
               NEMODSEL=MODNEUB(I)
C               IF (I.EQ.1 .OR. I.EQ.12. OR. I.EQ.15 .OR. I.EQ.16)IFREE=1
               IF (I.EQ.1 .OR. I.EQ.12. OR. I.EQ.15 .OR. I.EQ.16
     $              .OR. I.eq.29 .OR. I.eq.30) THEN
                  IBOUND=0
               endif
               IF (I.EQ.3 .OR. I.EQ.4 .OR. I.EQ.7 .OR. I.EQ.9 .OR. 
     &           I.EQ.17 .OR. I.EQ.19 .OR. I.EQ.20 .OR. I.EQ.22
     &              .OR. I.EQ.25 .OR. I.EQ.27) THEN
C                  IF (RLU(DUM).LT.(RATFRE/RATP)) IFREE = 1
                  IF (RLU(DUM).LT.(RATFRE/RATP)) IBOUND = 0
               ENDIF
C               WRITE(*,*) "MODE",NEMODSEL,E,IPAR
               RETURN
            ENDIF
 120     CONTINUE
         GO TO 110
      ENDIF
      
      RETURN
      END
