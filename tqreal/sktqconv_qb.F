************************************************************************
*     ------------------------------------------
CKO      SUBROUTINE SKTQCONV_qb(ICAB, IT0, IT, IQRAW, T, Q)
      SUBROUTINE SKTQCONV_qb(ICAB, IT0, IT, IQRAW, Qshift, IFLAG, T, Q)
*     ------------------------------------------
*
*     (Purpose)  convert TDC/ADC count into real T/Q for QB
*
* 1. apply QB constant
* 2. apply TQMAP
* 3. apply pc2pe factor
* 4. apply PGAIN 
*
*     (Input)
*
*       ICAB  : serial cable number (1-MAXPM4)
*       IT0   : T0 from software trigger
*       IT    : TDC count (10bit)
*       IQRAW : raw QDC count (10bit) + range flag
*
*       Qshift : correction for Q count due to QTC multi-ch hit effect
*       IFLAG  : == -1  ->  do not apply pc2pe/PGAIN (Q in pC)
*                != -1  ->  apply pc2pe/PGAIN (Q in p.e.)
*
*     (Output)
*       
*       T    : T (nsec)
*       Q    : Q (p.e.)
*
*     input from common
*
*       NRUNSK
*
*
*     07-JUL-2007 ported from sktqconv.F by Y.Takeuchi
*       o modified for QB 
*     11-JUL-2007 modified by Y.Takeuchi
*       o skip after TQMAP, if MA_TQMAP.ge.99
*     13-JUL-2007 modified by Y.Takeuchi
*       o new T/Q calc from Nakayama-san
*       o use /skofl/const/qbconst_XX.YY
*     18-JUL-2007 modified by Y.Takeuchi
*       o comment out pedq=0
*     20-JUL-2007 modified by Y.Takeuchi
*       o move T offset (11700) after TQREAL to before TQREAL
*     17-SEP-2008 modified by Y.Takeuchi
*       o just skip TQMAP and/or PGAIN when the major version=99
*     26-DEC-2008 modified by S.Nakayama/K.Okumura
*       o added Qshift and IFLAG (for QTC multi-ch hit correction)
* 
************************************************************************
      IMPLICIT NONE

#include "skparm.h"
#include "skcalib.h"
#include "skodcalib.h"
#include "skhead.h"
#include "skpdst.h"

      INTEGER ICAB, IT, IQRAW, IT0
      REAL T, Q, X
      INTEGER IFLAG
      REAL Qshift, RQ

      real c, dped
      integer IQ, irange, nev_error
      data nev_error/-1/
      save nev_error

      REAL PC2PE_SK1, PC2PE_SK2, PC2PE_SK3, PC2PE_SK4, PC2PE_SK5, PC2PE
      PARAMETER (PC2PE_SK1 = 2.055, PC2PE_SK2 = 2.297, PC2PE_SK3 = 2.243,
     &     PC2PE_SK4 = 1., PC2PE_SK5 = 1.)
      DATA PC2PE/0./
      SAVE PC2PE

      INTEGER I, ISTA, IMIN, NRUNSK_LAST 

      INTEGER MA_QBCONST1, MA_QBCONST_LAST, MA_PGAIN1, MA_PGAIN_LAST
      INTEGER MA_TQMAP1, MA_TQMAP_LAST
      INTEGER MI_QBCONST1, MI_QBCONST_LAST, MI_PGAIN1, MI_PGAIN_LAST
      INTEGER MI_TQMAP1, MI_TQMAP_LAST
      DATA MA_QBCONST_LAST /-1/
      DATA MA_TQMAP_LAST /-1/
      DATA MA_PGAIN_LAST /-1/
      DATA MI_QBCONST_LAST /-1/
      DATA MI_TQMAP_LAST /-1/
      DATA MI_PGAIN_LAST /-1/
      save ma_qbconst_last, ma_tqmap_last, ma_pgain_last
      save mi_qbconst_last, mi_tqmap_last, mi_pgain_last
      SAVE NRUNSK_LAST


      integer ma_od1,ma_od_last
      integer mi_od1,mi_od_last
      data ma_od_last /-1/
      data mi_od_last /-1/
      save ma_od_last
      save mi_od_last

c*** init.
      T = -9999.
      Q = -9999.

c*** check cable #
      if (icab.lt.0.or.icab.gt.MAXPM4) then
         print *, 'sktqconv_qb: cable error! icab=', icab
         return
      endif

c*** to avoid icab=0, 2011/04/15 added by y.takeuchi
c*** (I have not investigated why icab=0...) 
      if (icab.eq.0) return

c*** run chekck
      IF(NRUNSK.lt.40000) then
         print *, 'sktqconv_qb: invalid run number NRUNSK = ', nrunsk
         return
      endif

      IF (NRUNSK .NE. NRUNSK_LAST) THEN
         NRUNSK_LAST = NRUNSK

C*** set PC2PE
         if(NRUNSK.ge.40000.and.NRUNSK.lt.60000) then
            PC2PE=PC2PE_SK3
         else if(NRUNSK.ge.60000.and.NRUNSK.lt.78100) then
            PC2PE=PC2PE_SK4
         else if(NRUNSK.ge.78100) then
            PC2PE=PC2PE_SK5
         else
            write(6,*) 'Not supported pc2pe run number = ', NRUNSK
            stop
         endif

c*** set calibration table versions from run number
         CALL SKCALTABVER(NRUNSK, MA_QBCONST1, MI_QBCONST1, 
     $        MA_TQMAP1, MI_TQMAP1, MA_PGAIN1, MI_PGAIN1)

         WRITE(6,10) NRUNSK, MA_QBCONST1, MI_QBCONST1, 
     $        MA_TQMAP1, MI_TQMAP1, MA_PGAIN1, MI_PGAIN1
 10      FORMAT(1X, 'sktqconv_qb: nrunsk=', I6, 
     $        '  table version = ', 3(I2.2,'.',I2.2,' '))
C
C*** set QBCONST table
 20      FORMAT(1X, 'sktqconv_qb: setting ',A5,' version ',
     $        I2.2, '.', I2.2,' ...')
 21      FORMAT(1X, 'sktqconv_qb: setting ',A7,' version ',
     $        I2.2, '.', I2.2,' ...')


C
C*** set QB const
         IF (MA_QBCONST1 .NE. MA_QBCONST_LAST .OR. 
     $       MI_QBCONST1 .NE. MI_QBCONST_LAST) THEN
            WRITE(6,21) 'QBCONST', MA_QBCONST1, 0
            CALL SKCALTABSET_QBCONST(MA_QBCONST1, 0, ISTA)
            WRITE(6,*) 'finish!'
            IF (MI_QBCONST1 .GT. 0) THEN
               WRITE(6,21) 'QBCONST', MA_QBCONST1, MI_QBCONST1
               CALL SKCALTABSET_QBCONST(MA_QBCONST1, MI_QBCONST1, ISTA)
               WRITE(6,*) 'finish!'
            ENDIF
         ENDIF
         MA_QBCONST_LAST = MA_QBCONST1
         MI_QBCONST_LAST = MI_QBCONST1

C
C*** set TQ-map
         IF (MA_TQMAP1 .NE. MA_TQMAP_LAST .OR. 
     $       MI_TQMAP1 .NE. MI_TQMAP_LAST) THEN
            WRITE(6,20) 'TQMAP', MA_TQMAP1, 0
c*** 11-JUL-2007 MA_TQMAP.ge.99 check is added by y.t.
            if (MA_TQMAP1.ge.99) then
               MA_TQMAP = MA_TQMAP1 
               print *, 'skip TQMAP MA_TQMAP = ',MA_TQMAP
            else
               CALL SKCALTABSET_TQMAP(MA_TQMAP1, 0, ISTA)
               WRITE(6,*) 'finish!'
               IF (MI_TQMAP1 .GT. 0) THEN
                  WRITE(6,20) 'TQMAP', MA_TQMAP1, MI_TQMAP1
                  CALL SKCALTABSET_TQMAP(MA_TQMAP1, MI_TQMAP1, ISTA)
                  WRITE(6,*) 'finish!'
               ENDIF
            endif
         ENDIF
         MA_TQMAP_LAST = MA_TQMAP1
         MI_TQMAP_LAST = MI_TQMAP1
C
C*** set PMT gain table
         IF (MA_PGAIN1 .NE. MA_PGAIN_LAST .OR. 
     $       MI_PGAIN1 .NE. MI_PGAIN_LAST) THEN
            WRITE(6,20) 'PGAIN', MA_PGAIN1, 0
c*** 26-SEP-2007 MA_PGAIN.ge.99 check is added by y.t.
            if (MA_PGAIN1.ge.99) then
               MA_PGAIN = MA_PGAIN1 
               print *, 'skip PGAIN MA_PGAIN = ',MA_PGAIN
            else
               CALL SKCALTABSET_PGAIN(MA_PGAIN1, 0, ISTA)
               WRITE(6,*) 'finish!'
               IF (MI_PGAIN1 .GT. 0) THEN
                  WRITE(6,20) 'PGAIN', MA_PGAIN1, MI_PGAIN1
                  CALL SKCALTABSET_PGAIN(MA_PGAIN1, MI_PGAIN1, ISTA)
                  WRITE(6,*) 'finish!'
               ENDIF
            endif
         ENDIF
         MA_PGAIN_LAST = MA_PGAIN1
         MI_PGAIN_LAST = MI_PGAIN1


c*** set calibration table versions from run number for OD
         CALL SKCALTABVER_OD(NRUNSK, MA_OD1, MI_OD1)

         WRITE(6,10) NRUNSK, MA_OD1, MI_OD1
 11      FORMAT(1X, 'sktqconv_qb: nrunsk=', I6, 
     $        'od table version = ', I2.2,'.',I2.2,' ')

C
C*** set OD tables
         IF (MA_OD1 .NE. MA_OD_LAST .OR. 
     $       MI_OD1 .NE. MI_OD_LAST) THEN
            WRITE(6,20) 'OD', MA_OD1, 0
            CALL SKCALTABSET_OD(MA_OD1, 0, ISTA)
            WRITE(6,*) 'finish!'
            IF (MI_OD1 .GT. 0) THEN
               WRITE(6,20) 'OD', MA_OD1, MI_OD1
               CALL SKCALTABSET_OD(MA_OD1, MI_OD1, ISTA)
               WRITE(6,*) 'finish!'
            ENDIF
         ENDIF
         MA_OD_LAST = MA_OD1
         MI_OD_LAST = MI_OD1

      ENDIF ! End of new run if



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C*** QBCONST (all channels)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c*** check pedestal (range 0) iqbcalmode is checked in caltabset_qb()
      if (IQBPDSTQ(icab,0,iqbcalmode).eq.0) then
         if (nevsk.ne.nev_error) then
ccyt            print *,'sktqconv_qb: QB pedq=0 at icab,nevsk=', icab,nevsk
            nev_error = nevsk
         endif
         return
c         stop
      endif

      if (IQBPDSTT(icab,0,iqbcalmode).eq.0) then
         if (nevsk.ne.nev_error) then
            print *,'sktqconv_qb: QB pedt=0 at icab,nevsk=', icab,nevsk
            nev_error = nevsk
         endif
         return
c         stop
      endif

c*** not used now
c      if (IQBCALQ(icab,0,iqbcalmode).eq.0) then
c         print *, 'sktqconv_qb: QBCAL=0 at icab=', icab
c         stop
c      endif

c*** Q

c*** separate range bit and ADC count
      call SKQDIV(IQRAW, IQ, IRANGE)
      RQ = float(IQ)+Qshift

c*** apply temperature correction
CKO      c = (qb_cal0mean(icab,irange,iqbcalmode) 
CKO     $     - qb_ped0mean(icab,irange,iqbcalmode))
CKO     $     * (float(IQ) - float(IQBPDSTQ(icab,irange,iqbcalmode))/10.0) 
CKO     $     / (float(IQBCALQ(icab,irange,iqbcalmode) 
CKO     $     - IQBPDSTQ(icab,irange,iqbcalmode))/10.0)
      c = (qb_cal0mean(icab,irange,iqbcalmode) 
     $     - qb_ped0mean(icab,irange,iqbcalmode))
     $     * (RQ - float(IQBPDSTQ(icab,irange,iqbcalmode))/10.0) 
     $     / (float(IQBCALQ(icab,irange,iqbcalmode) 
     $     - IQBPDSTQ(icab,irange,iqbcalmode))/10.0)

c*** tentative 13-JUL-2007 (y.t.)
c      dped = float(IQBPDSTQ(icab,irange,iqbcalmode))/10.0 
c     $     - qb_ped0mean(icab,irange,iqbcalmode)
c      c = (1.0 - 0.01208*dped + 0.0003375*(dped**2))
c     $     * (float(IQ) - float(IQBPDSTQ(icab,irange,iqbcalmode))/10.0) 

c*** check range 
      if (irange.eq.0) then
cc         if (c.ge.qb_cntmin(icab,-1) .and. c.le.qb_cntmax(icab,-1)) then
         if (c.le.qb_cntmax(icab,-1)) then
            irange = -1
cc         else if(c.ge.qb_cntmin(icab,0).and.c.le.qb_cntmax(icab,0)) then
         else
            irange = 0
cc         else
cc            print *, 'sktqconv_qb: range error! count=', c
         endif
      endif

c*** apply function
      Q = ((((qb_a5(icab,irange)*c + qb_a4(icab,irange))*c 
     $     + qb_a3(icab,irange))*c + qb_a2(icab,irange))*c
     $     + qb_a1(icab,irange))*c + qb_a0(icab,irange)

c      print *,nevsk,icab,irange,iqbcalmode,c,Q,dped,
c     $     IQBPDSTQ(icab,irange,iqbcalmode),
c     $     qb_ped0mean(icab,irange,iqbcalmode)

c*** T
c      T = (IT - IT0) / count_per_nsec   ! 1/60e6/32(sec) = 1/1.92(nsec) 

c*** updated 13-JUL-2007 (y.t.)
      T = (float(IT-IT0)-float(IQBPDSTT(icab,irange,iqbcalmode))/10.000) 
     $     / count_per_nsec

c*** add offset for monitor/veto channel
      if(icab.gt.maxpm) then
         if (icab.lt.20000) then
            T = T + float(IQBPDSTT(15001,0,1))/10.000/count_per_nsec + 1000
         else
            T = T + 1538.
         endif
         goto 100
      endif

**** For PMT Timing

c*** add offset
c         T = T - 200.0    ! to make T0 arround 700   
C         T = T +11700.0    ! to make T0 arround 700  for QB test data y.t. 
      if (nrunsk.lt.50791) then
         T = T +3900.           ! to make T0 arround 1000  for QB test data koshio 
      else if (nrunsk.eq.52026) then
         T = T +5705.           ! to make T0 arround 1000  for QB test data koshio 
      else if (nrunsk.lt.60001) then
         T = T +5205.           ! to make T0 arround 1000  for QB test data shoei
      else
c*** Merge offset 1 and 2 from Obayashi-san's request
c         T = T +1538.           ! to make T0 arround 1000 (T-pedestal changed) shoei

      endif


      if (IFLAG.eq.-1) return   ! do not apply pc2pe/PGAIN


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C*** TQMAP (only ID PMTs)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c*** skip TQMAP if the major version = 99 (to make TQMAP)
      if (MA_TQMAP.ne.99) then

c*** make Q-bin
         if (icab.le.MAXPM) then
ccc x-axis = Q(pc)
            IF (Q.GE.10.000) THEN
               X = LOG10(Q) * 50.000 
            ELSE
               X = Q * 5.0000
            ENDIF

c*** limit x-axis from 0-175   1500p.e. log(3000)*50= 173.8 (tentative)
            if (x.lt.0.0)   x = 0.0
cc            if (x.gt.175.0) x = 175.0
c*** limit x-axis for SK-III TQMAP 140bin
            if (x.gt.170.0) x = 170.0

c*** apply TQMAP function
          
	    if(MA_TQMAP.eq.50) then ! P7P7P7 function
              if (x.lt.10.0) then
                  
               T = ((((((tqmap_qb(icab,8)*X + tqmap_qb(icab,7))*X 
     $              + tqmap_qb(icab,6))*X + tqmap_qb(icab,5))*X
     $              + tqmap_qb(icab,4))*X + tqmap_qb(icab,3))*X
     $              + tqmap_qb(icab,2))*X + tqmap_qb(icab,1) + T 
               
              else if(x.lt. 50.0) then
               
               T = ((((((tqmap_qb(icab,8)*10.0 + tqmap_qb(icab,7))*10.0 
     $              + tqmap_qb(icab,6))*10.0 + tqmap_qb(icab,5))*10.0
     $              + tqmap_qb(icab,4))*10.0 + tqmap_qb(icab,3))*10.0
     $              + tqmap_qb(icab,2))*10.0 + tqmap_qb(icab,1) 
     $              + (((((((tqmap_qb(icab,15)*(x-10.0) + tqmap_qb(icab,14))*(x-10.0))
     $              + tqmap_qb(icab,13))*(x-10.0) + tqmap_qb(icab,12))*(x-10.0)
     $              + tqmap_qb(icab,11))*(x-10.0) + tqmap_qb(icab,10))*(x-10.0)
     $              + tqmap_qb(icab, 9))*(x-10.0) + T 
               
              else
               
               T = ((((((tqmap_qb(icab,8)*10.0 + tqmap_qb(icab,7))*10.0 
     $              + tqmap_qb(icab,6))*10.0 + tqmap_qb(icab,5))*10.0
     $              + tqmap_qb(icab,4))*10.0 + tqmap_qb(icab,3))*10.0
     $              + tqmap_qb(icab,2))*10.0 + tqmap_qb(icab,1) 
     $              + (((((((tqmap_qb(icab,15)*(50.0-10.0) + tqmap_qb(icab,14))*(50.0-10.0))
     $              + tqmap_qb(icab,13))*(50.0-10.0) + tqmap_qb(icab,12))*(50.0-10.0)
     $              + tqmap_qb(icab,11))*(50.0-10.0) + tqmap_qb(icab,10))*(50.0-10.0)
     $              + tqmap_qb(icab, 9))*(50.0-10.0)                      
     $              + (((((((tqmap_qb(icab,22)*(x-50.0) + tqmap_qb(icab,21))*(x-50.0))
     $              + tqmap_qb(icab,20))*(x-50.0) + tqmap_qb(icab,19))*(x-50.0)
     $              + tqmap_qb(icab,18))*(x-50.0) + tqmap_qb(icab,17))*(x-50.0)
     $              + tqmap_qb(icab,16))*(x-50.0) + T
              endif

	    else ! P3P5P7 function for MA_TQMAP.ge.51 
              if (x.lt.10.0) then
               T = ((tqmap_qb(icab,4)*X + tqmap_qb(icab,3))*X
     $              +tqmap_qb(icab,2))*X + tqmap_qb(icab,1) + T


              else if(x.lt. 50.0 ) then
               T = ((tqmap_qb(icab,4)*10.0 + tqmap_qb(icab,3))*10.0
     $              + tqmap_qb(icab,2))*10.0 + tqmap_qb(icab,1)
     $              + ((((tqmap_qb(icab, 8)*(x-10.0) + tqmap_qb(icab,7))*(x-10.0)
     $              + tqmap_qb(icab, 6))*(x-10.0) + tqmap_qb(icab, 5))*(x-10.0)
     $              + 3*tqmap_qb(icab,4)*10.*10.+2*tqmap_qb(icab,3)*10.+tqmap_qb(icab,2))*(x-10.0) + T

              else
               T =  ((tqmap_qb(icab,4)*10.0 + tqmap_qb(icab,3))*10.0
     $              + tqmap_qb(icab,2))*10.0 + tqmap_qb(icab,1)
     $              + ((((tqmap_qb(icab, 8)*(50.0-10.0) + tqmap_qb(icab,7))*(50.0-10.0)
     $              + tqmap_qb(icab, 6))*(50.0-10.0) + tqmap_qb(icab, 5))*(50.0-10.0)
     $              + 3*tqmap_qb(icab,4)*10.*10.+2*tqmap_qb(icab,3)*10.+tqmap_qb(icab,2))*(50.0-10.0)
     $              + (((((((tqmap_qb(icab,15)*(x-50.0) + tqmap_qb(icab,14))*(x-50.0))
     $              + tqmap_qb(icab,13))*(x-50.0) + tqmap_qb(icab,12))*(x-50.0)
     $              + tqmap_qb(icab,11))*(x-50.0) + tqmap_qb(icab,10))*(x-50.0)
     $              + tqmap_qb(icab,9))*(x-50.0) + T
              endif
            endif


c*** Merge offset 1 and 2 from Obayashi-san's request
c*** add offset 2
c            T = T - 1255  ! to make the peak at T=1000 for Run 52130 with tqmap_31.00
c            T = T + 283.  ! removed on 2008/10/28 and taken into account in TQmap. Y.Obayashi
            if (nrunsk.eq.52026) T = T + 290 ! fine tune for Ni run

         endif

      endif

 100  continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C*** pc2pe & PMT GAIN (only ID PMTs)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c*** skip PGAIN if the major version = 99 (to check TQMAP)
      if (MA_PGAIN.ne.99) then

         if (icab.le.MAXPM) then
            if(PC2PE.eq.0. .or. PGAIN_qb(ICAB).eq.0.) then
               write(6,*) 'sktqconv_qb: Error ICAB/PC2PE/PGAIN=', 
     $              icab, pc2pe, pgain_qb(icab)
               stop
            else
               Q = Q / PC2PE/ PGAIN_qb(ICAB)
            endif
         endif

      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C*** OD PMTs
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (icab.gt.20000.and.icab.le.20000+maxpma) then

            Q = Q/ODGAIN(ICAB-20000)
            T = T + ODTOFF(ICAB-20000)

      endif


      RETURN
      END
