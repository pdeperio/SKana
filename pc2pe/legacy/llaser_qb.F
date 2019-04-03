***   ==============                                                 ***
      PROGRAM LLASER
***   ==============                                                 ***
*     
*     (Purpose)
*        Read large Q laser data and calc Qmean for pc2pe analysis for QBee
*
*     (Input)
*        LUNI ; laser run data
*
*     (Output) 
*        ******.dat (****** = run #)
*
*     (Creation Date and Author)
*       2006.08.11  by T Iida
*
*       Copied from: /home/iida/calib/pc2pe/pc2pe_08oct/llaser_qb.F
*
****************************************************************

      IMPLICIT NONE

#include "skparm.h"
#include "skhead.h"
#include "sktq.h"
#include "skpdst.h"
#include "geopmt.h"      
#include "skveto.h"
#include "geotnk.h"
#include "pawc.h"
#include "dtdeadtime.h"

      common/pawc/pawc(200000)

*     LUN
      INTEGER LUNI,LUNA,LUNCORRA,LUNNOA,LUNALLA

      PARAMETER (LUNI      = 10)
      PARAMETER (LUNA      = 21)
      PARAMETER (LUNCORRA  = 25)
      PARAMETER (LUNNOA    = 27)
      PARAMETER (LUNALLA   = 29)

      CHARACTER*100 FILENAME,FILENAME2,FILENAME0,FILEPLACE,
     &     dataname,name,name2,fname_in
      
*     some parameters
      real Qmean_ton(11146),Qmean_toff(11146),Qmean(11146)
      INTEGER I, J, IERR ,select(10), iflag,pmtflg
      integer maxqeflag(11146),ifile, narg, npmt

      real ball_pos(3),dist,ttof,c,expd,tiskzmont
      real ton_max,toff_max,ton_min,toff_min,tiskz2

c      PARAMETER (ton_min=950.0, ton_max=1300.0)
c      PARAMETER (toff_min=600.0, toff_max=950.0)

C***  parameter for SK-IV
      PARAMETER (ton_min=-1100.0, ton_max=1450.0)
c      PARAMETER (ton_min=-1100.0, ton_max=1300.0)
c      PARAMETER (ton_min=600.0, ton_max=850.0)
      PARAMETER (toff_min=0.0, toff_max=300.0)
      PARAMETER (c = 21.6438)

      real qisk_ton(maxpm),qisk_toff(maxpm),qismsk_ton,qismsk_toff,r
      integer nev_ton(maxpm),nev_toff(maxpm),nlaser,itshift
      integer hid,hid2, nev

      integer cab,address1,hut,tko,atm_board,ntotal
      integer atm_ch,hv_crate,hv_board,hv_ch,a,pre,sk
      integer sk1(maxpm),hut1(maxpm),pre1(maxpm),ID
      real hv_hama,hv_sk2
      real skb,hv_pre
      character*10 ser,address2
      real*8 beta8

*     PAW
      real pawc

      CALL KZINIT
      call hlimit(-500000)
c      SK_GEOMETRY=4
      call skoptn('31,30,19,18')
      call skcontabrd

*     ball position
      
      ball_pos(1)=35.3
      ball_pos(2)=-70.7
      ball_pos(3)=0.
      
*     get PMT information
      i=0
      open(13,file='sk3_hv.table4')
 30   read (13,*,end=300) cab,ser,address1,address2
     $     ,hut,tko,atm_board,atm_ch,hv_crate
     $     ,hv_board,hv_ch,skb,hv_hama,hv_sk2,hv_pre
     $     ,a,pre,sk,beta8
      sk1(cab) = sk
      hut1(cab) = hut
      pre1(cab) = pre
      i=i+1
      if (i.eq.11146) goto 300
      goto 30
 300  continue
      close (13)
      
*     --------------------------
*     file name
*     --------------------------
      call getarg(1,dataname)
      data name/'.dat'/
      data name2/'.hbk'/
      FILENAME  = dataname(1:6)//name
      FILENAME2  = dataname(1:6)//name2
      
*     ---------------------------
*     save file
*     ---------------------------
      open(LUNA,FILE=FILENAME)

*     ------
*     HBOOK
*     ------
      call hbook1(10,' T-TOF ',3000,0.,3000., 0.)
      call hbook1(11,' T-TOF ',3000,0.,3000., 0.)
      call hbook1(16,' Tiskz ',3000,0.,3000., 0.)
      call hbook1(17,' Tiskz ',3000,-10000.,20000., 0.)
      call hbook1(19,' Tisk ',3000,0.,3000., 0.)
      call hbook1(101,' Tisk ',3000,0.,3000., 0.)
      call hbook1(125,' Tisk ',3000,0.,3000., 0.)
      call hbook2(200,' Qiskz vs Tiskz ',3000,0.,3000., 
     $     100, 0., 100. ,0.)

      call hbook1(40,' QISKZ ',300,-5.,10., 0.)
      call hbook1(50,' total Q ',500,0.,1000000., 0.)
      call hbook1(51,' total Q ',500,300000.,1300000., 0.)
      call hbook1(52,' total Q ',200,300000.,500000., 0.)


*     ----------
*     Initialize
*     ----------

      nev=0
      DO J=1,11146

         Qmean_ton(J)  = 0.0
         Qmean_toff(J)  = 0.0
         maxqeflag(J) = 0
      ENDDO
      
      narg=iargc()
      do ifile=2,narg

      call getarg(ifile,fname_in)

*     ----------
*     Open Files
*     ----------
      call set_rflist(LUNI,fname_in,
     $     'DISK',' ','RED',' ',' ',
     $     'recl=5670 status=old',' ',' ')

      CALL SKOPENF(LUNI,0,'Z',IERR)
c      call skoptn('31,30,29,19,18')

 1000 continue
      CALL SKrawREAD(10,*1002,*1001,*1000,*1000)
      CALL SKREAD(-10,*1002,*1001,*1000,*1000)
      if(nevsk.eq.0) goto 1000
c*** skip pedestal(is not ready) events
      if (iand(idtgsk, 2**30).ne.0) goto 1000

c*** skip pedestal on
      if (iand(ifevsk, 2**9).ne.0) goto 1000

c*** skip incomplete
      if (iand(ifevsk, 2**20).ne.0) goto 1000

c*** skip inner off
      if (iand(ifevsk, 2**28).ne.0) goto 1000

c*** skip slow data
      if (iand(ifevsk, 2**14).ne.0) goto 1000

c*** skip run info
      if (iand(ifevsk, 2**15).ne.0) goto 1000

C*****check if laser monitor hit exist 
      pmtflg=0
      do i = 1,NQISKZ           ! NQISKZ is defined in sktq.h
         npmt = ICABIZ(i)       ! ICABIZ(I) : hit cable number
         itshift = itiskz(i)
         tiskzmont = tiskz(i)
         if ( npmt .eq. 15012 .and. itshift.eq.0) then
            write(*,*)'!!!!! select laser !!!!!',itiskz(i),tiskz(i)
            goto 555
         endif
      enddo
      goto 1000
 555  write(*,*) nevsk,qismsk

C***  Set T0 for LSI laser
c      call set_timing_gate(it0sk + itshift)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     check IFEVSK flag
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     check run information or slow data

      select(1)=select(1)+1

      If (iand(ifevsk,dtruninformation).ne.0
     &     .or.iand(ifevsk,dtslowdata).ne.0) then
         select(2)=select(2)+1
         goto 1000
      endif
      
C     check unreasonabel flag
      If (iand(ifevsk,dtunreasonable).ne.0) then
         select(3)=select(3)+1
         goto 1000
      endif
      
C     now pedestal taking period?
      If (iand(ifevsk,dtpedestal_taking).ne.0) then
         select(4)=select(4)+1
         goto 1000
      endif
      
C     check inner ws on/off flag
      If (iand(ifevsk,dtinner_ws_off).ne.0
     &     .or.iand(ifevsk,dtuncompleted).ne.0) then
         select(5)= select(5)+1
         goto 1000
      endif

**********************************************

      do J=1,NQISK
         call hf1(19,tisk(J),1.0)
      enddo
      do J=1,NQISKZ

         if ( ICABIZ(J).ge.1. and. ICABIZ(J).le.11146) then
            tiskz2=tiskz(J)
            call hf1(16,tiskz2,1.0)
            call hf1(17,tiskz2,1.0)
         endif
      enddo

      ntotal = ntotal + 1
      if(mod(ntotal,100).eq.0) then
         write(*,*) ifile,nrunsk,nevsk,qismsk
      endif

      goto 200
      
c     ----------------
c     make TvsQ sample 
c     ----------------
      
      if(mod(nevsk,100).eq.0) then 
         write(dataname,*)nevsk
         dataname = dataname(4:10)
         FILEPLACE = '/net/sukcald1/work1/iida/pc2pe_newtq/sample/'
         
         if (nevsk.lt.1000) then
            FILENAME0 = FILEPLACE(1:44)//dataname(1:3)//name
         else if (nevsk.lt.5000) then
            FILENAME0 = FILEPLACE(1:44)//dataname(1:4)//name
         else 
            goto 200
         endif
         
         write(*,*)FILENAME0
         open(999,FILE=FILENAME0)

      endif
      
 200  continue
      
      nlaser = nlaser + 1
      do i = 1,NQISKZ           ! NQISKZ is defined in sktq.h

         cab=ICABIZ(i)
         if (cab.lt.1 .or. cab.gt.11146) goto 400

         dist = sqrt((ball_pos(1)-xyzpm(1,cab))**2+(ball_pos(2)
     $        -xyzpm(2,cab))**2+(ball_pos(3)-xyzpm(3,cab))**2)
         ttof = real(tiskz(i)) - dist/c
         r = sqrt(xyzpm(1,cab)**2+xyzpm(2,cab)**2)
         
         call hf1(40,qiskz(i),1.0)
c         write(80,*) qiskz(i)
         if (qiskz(i).ne.0.) call hf1(10,ttof,1.0)
         if (qiskz(i).eq.0.) call hf1(11,ttof,1.0)
         call hf2(200,ttof,qiskz(i),1.0)
         if (cab.eq.1.) call hf1(101,ttof,1.0)
         if (cab.eq.25) call hf1(125,ttof,1.0)
         if (ttof.gt.1700. .and. ttof.lt.1750.) 
     $        write(57,*) nevsk,cab,qiskz(i)
         
         if(qiskz(i).ge.0.0)then
            if (ttof.gt.ton_min.and.ttof.lt.ton_max) then
               nev_ton(cab) = nev_ton(cab) + 1
               qisk_ton(cab) = qisk_ton(cab)+ qiskz(i)
c               call hf1(40,qiskz(i),1.0)
            else if (ttof.gt.toff_min.and.ttof.lt.toff_max) then
               nev_toff(cab) = nev_toff(cab) + 1
               qisk_toff(cab) = qisk_toff(cab) + qiskz(i)
            endif
         endif
         
 400     continue
      ENDDO
      call hf1(50,qismsk,1.0)
      call hf1(51,qismsk,1.0)
      call hf1(52,qismsk,1.0)
      
      CALL FLUSH(6)
      CALL FLUSH(0)
      goto 1000
      
 1002 CONTINUE
      WRITE(0,*) 'llaser: Read error'
      
 1001 CONTINUE
      
      WRITE(0,*) '**** end of read ****'
      CALL SKCLOSEF(LUNI)
      enddo
      
      write(*,*)'# of event = ',select(1)
      write(*,*)'rejected event = ',select(2),select(3)
     $     ,select(4),select(5)
      
*     ---------------
*     Calc. Mean of Q
*     ---------------
      DO J=1,11146
         
         Qmean_ton(J) = qisk_ton(J)/real(nlaser)
         Qmean_toff(J) = qisk_toff(J)/real(nlaser)
c         Qmean(J) = Qmean_ton(J) - Qmean_toff(J)

***   old caluculaton method   ***
c         if (nev_ton(J).ge.1) then
c            Qmean_ton(J) = qisk_ton(J)/real(nev_ton(J))
c         else if (nev_ton(J).le.0) then
c            Qmean_ton(J) = 0.0
c         endif
c         if (nev_toff(J).ge.1) then
c            Qmean_toff(J) = qisk_toff(J)/real(nev_toff(J))
c         else if (nev_toff(J).le.0) then
c            Qmean_toff(J) = 0.0
c         endif

      ENDDO

*     -------------
*     Write Results
*     -------------
      
      DO I=1,11146
         WRITE(LUNA,'(I6,I6,F10.5)')
c     &        I,nev_ton(I),Qmean(I)
     &        I,nev_ton(I),Qmean_ton(I)

c         WRITE(99,'(I6,I6,F10.5)')
c     &        I,nev_toff(I),Qmean_toff(I)
         call hf1(110,real(I),Qmean_ton(I))
         call hf1(120,real(I),Qmean_toff(I))
      END DO
      
*     -----------
*     Close Files
*     -----------      
c      CALL SKCLOSEF(LUNI)
c      call hrput(0,'ltisk_tail.hbk','nt')      

      CLOSE(LUNA)
      call hrput(0,FILENAME2,'N')

      stop
      END

      real function scntunifrm(radang)
      implicit none
      real x, radang
      x = radang
      scntunifrm = 1.
      return
      end
      


