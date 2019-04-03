***   ==============                                                 ***
      PROGRAM SLASER_QB
***   ==============                                                 ***
*     
*     (Purpose)
*        Read small Q laser data and calc hit rate for pc2pe analysis
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
*       Copied from: /home/iida/calib/pc2pe/pc2pe_08oct/slaser_qb.F
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

      CHARACTER*100 FILENAME,FILENAME0,FILEPLACE,
     &     dataname,name,name2, fname_in,FILENAME2,
     $     FILENAME3,name3

*     some parameters
      integer nhit_ton(11146),nhit_toff(11146),narg,ifile
      INTEGER I, J, IERR ,select(10), iflag
      integer maxqeflag(11146),cab,nev_laser,nhit(maxpm)
      real rhit_ton(maxpm),rhit_toff(maxpm),rhit(maxpm)
      real rhit_ton_occu(maxpm),rhit_toff_occu(maxpm)
     $     ,rhit_occu(maxpm),nhit_occu

      integer pmtflg, npmt,itshift
      real ball_pos(3),dist,ttof,c,tiskzmont
      real ton_max,toff_max,ton_min,toff_min
c      PARAMETER (ton_min=950.0, ton_max=1300.0)
c      PARAMETER (toff_min=600.0, toff_max=950.0)

C     *** parameter for SK-IV
c      PARAMETER (ton_min=1150.0, ton_max=1250.0)
c      PARAMETER (toff_min=500.0, toff_max=600.0)
      PARAMETER (ton_min=1100.0, ton_max=1450.0)
      PARAMETER (toff_min=750.0, toff_max=1100.0)
c      PARAMETER (toff_min=1450.0, toff_max=1800.0)
      PARAMETER (c = 21.6438)
      
      real qisk_ton(maxpm),qisk_toff(maxpm),qismsk_ton,qismsk_toff
      integer nev_ton(maxpm),nev_toff(maxpm)

*     PAW
      real pawc

      CALL KZINIT
      call hlimit(-200000)
      call skoptn('31,30,19,18')
      call skcontabrd
c      SK_GEOMETRY=3


*     ball position

      ball_pos(1)=35.3
      ball_pos(2)=-70.7
      ball_pos(3)=0.

*     --------------------------
*     file name
*     --------------------------
*     --------------------------
*     file name
*     --------------------------
      call getarg(1,dataname)
      data name/'.dat'/
      data name2/'.hbk'/
      data name3/'.dat2'/
      FILENAME  = dataname(1:6)//name
      FILENAME2  = dataname(1:6)//name2
      FILENAME3  = dataname(1:6)//name3

*     ---------------------------
*     save file
*     ---------------------------

      open(LUNA,FILE=FILENAME)
      open(77,FILE=FILENAME3)
      
*     ------
*     HBOOK
*     ------

      call hbook1(201,' Cab#1 ',2000,0.,2000., 0.)
      call hbook1(225,' Cab#25 ',2000,0.,2000., 0.)
      call hbook1(252,' Cab#8452 ',2000,0.,2000., 0.)
      call hbook1(10,' all ',2000,0.,2000., 0.)
      call hbook1(11,' Q>0 ',500,0.,2000., 0.)
      call hbook1(12,' Q=0 ',500,0.,2000., 0.)
      call hbook1(17,' Tiskz ',3000,-10000.,20000., 0.)
      call hbook1(16,' Tiskz ',2000,0.,2000., 0.)
c      call hbook1(11,' barrel ',2000,0.,2000., 0.)
c      call hbook1(12,' top ',2000,0.,2000., 0.)
c      call hbook1(13,' bottom ',2000,0.,2000., 0.)
      call hbook1(20,' cab# 6157 ',2000,0.,2000., 0.)
      call hbook1(50,' number of hit ',100,0.,1000., 0.)
      call hbook1(100,' nhit ontime ',11146,0.5,11146.5, 0.)
      call hbook1(110,' nhit offtime ',11146,0.5,11146.5, 0.)

*     ----------
*     Initialize
*     ----------

         DO J=1,11146
            nhit_ton(J)  = 0.0
            nhit_toff(J)  = 0.0
            maxqeflag(J) = 0
         ENDDO
         
         narg=iargc()
         do ifile=2,narg

            call getarg(ifile,fname_in)
*     ----------
*     Open Files
*     ----------
            call set_rflist(LUNI,fname_in,
     $           'DISK',' ','RED',' ',' ',
     $           'recl=5670 status=old',' ',' ')
            
            CALL SKOPENF(LUNI, 0, 'Z', IERR)

 1000       continue
            CALL SKrawREAD(10,*1002,*1001,*1000,*1000)
            CALL SKREAD(-10,*1002,*1001,*1000,*1000)
            if(nevsk.eq.0) goto 1000
          
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!'
c*** skip pedestal(is not ready) events
            if (iand(idtgsk, 2**30).ne.0) goto 1000
            
c***  skip pedestal on
            if (iand(ifevsk, 2**9).ne.0) goto 1000
            
c***  skip incomplete
            if (iand(ifevsk, 2**20).ne.0) goto 1000
            
c***  skip inner off
            if (iand(ifevsk, 2**28).ne.0) goto 1000
            
c***  skip slow data
            if (iand(ifevsk, 2**14).ne.0) goto 1000

c***  skip run info
            if (iand(ifevsk, 2**15).ne.0) goto 1000

C***  check if laser monitor hit exist 
      pmtflg=0
      do i = 1,NQISKZ           ! NQISKZ is defined in sktq.h
         npmt = ICABIZ(i)       ! ICABIZ(I) : hit cable number
         itshift = itiskz(i)
         tiskzmont = tiskz(i)
         if ( npmt .eq. 15012 .and. itshift.eq.0) then
            write(*,*)'!!!!! select laser !!!!!',tiskz(i)
            goto 555
         endif
      enddo
      goto 1000
 555  write(*,*) nevsk,qismsk
            
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

c     check unreasonabel flag
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
*     --------------
*     Calc. hit rate
*     --------------
      nev_laser = nev_laser + 1

      if(mod(nev_laser,100).eq.0) print*,ifile,nevsk,nqisk

c      write(*,*)nqisk
      call hf1(50,real(nqisk),1.0)
      DO J=1,NQISKZ

         cab = ICABIZ(J)
         if (cab.lt.1 .or. cab.gt.11146) goto 400

         dist = sqrt((ball_pos(1)-xyzpm(1,cab))**2+(ball_pos(2)
     $        -xyzpm(2,cab))**2+(ball_pos(3)-xyzpm(3,cab))**2)
         ttof = tiskz(J) - dist/c

         if (ttof.gt.ton_min.and.ttof.lt.ton_max) then
            nhit_ton(cab) = nhit_ton(cab) + 1
         else if (ttof.gt.toff_min.and.ttof.lt.toff_max) then
            nhit_toff(cab) = nhit_toff(cab) + 1
         endif

         call hf1(10,ttof,1.0)
         call hf1(201,ttof,1.0)
         call hf1(225,ttof,1.0)
         call hf1(252,ttof,1.0)

         if(qisk(J).gt.0.) call hf1(11,ttof,1.0)
         if(qisk(J).eq.0.) call hf1(12,ttof,1.0)
         call hf1(17,real(tiskz(J)),1.0)
         call hf1(16,real(tiskz(J)),1.0)
         
 400     continue
         
      ENDDO
      
      CALL FLUSH(6)
      CALL FLUSH(0)
      
      goto 1000
      
 1002 CONTINUE
      WRITE(0,*) 'xen_z: Read error'
      
 1001 CONTINUE

      CALL SKCLOSEF(LUNI)
      enddo

      write(*,*)'# of event = ',select(1)
      write(*,*)'rejected event = ',select(2),select(3)
     $     ,select(4),select(5)
      
*     ---------------
*     Calc. Hit Rate
*     ---------------
      DO J=1,11146
         rhit_ton(J) = nhit_ton(J)/real(nev_laser)
         rhit_toff(J) = nhit_toff(J)/real(nev_laser)
         nhit(J) = nhit_ton(J)-nhit_toff(J)
         rhit(J) = real(nhit(J)) / real(nev_laser)
      ENDDO
*     -------------------
*     ocupancy correction
*     -------------------
      DO J=1,11146
         rhit_occu(J)=(real(nev_laser) - real(nhit(J)))
     $        / real(nev_laser)
         rhit_ton_occu(J)=(real(nev_laser) - real(nhit_ton(J)))
     $        / real(nev_laser)
         rhit_toff_occu(J)=(real(nev_laser) - real(nhit_toff(J)))
     $        / real(nev_laser)
c         write(99,*)J,nhit(J),rhit(J)
         rhit_occu(J)=-log(rhit_occu(J))
         rhit_ton_occu(J)=-log(rhit_ton_occu(J))
         rhit_toff_occu(J)=-log(rhit_toff_occu(J))
      ENDDO

*     -------------
*     Write Results
*     -------------
      
      DO I=1,11146
c         WRITE(88,'(I6,I6,I6)')
c     &        I,nhit_ton(I),nhit_toff(I)
         WRITE(77,'(I6,4f9.5)')
     &        I,rhit_ton(I),rhit_toff(I)
     $        ,rhit_ton_occu(I),rhit_toff_occu(I)

         nhit_occu = rhit_occu(I)*real(nev_laser)   
         WRITE(LUNA,'(I6,F12.3,F10.5)')
     &        I,nhit_occu,rhit_occu(I)
c         WRITE(LUNA,'(I6,I6,F10.5)')
c     &        I,nhit(I),rhit(I)

         call hf1(100,real(I),real(nhit_ton(I)))
         call hf1(110,real(I),real(nhit_toff(I)))

      END DO
      
*     -----------
*     Close Files
*     -----------      
c      CALL SKCLOSEF(LUNI)
c      call hrput(0,'stisk.hbk','nt')      
      call hrput(0,FILENAME2,'nt')      

      CLOSE(LUNA)
      CLOSE(77)
      CLOSE(88)
      CLOSE(99)

c      call hrput(0,'hist.hbk','N')

      stop
      END

      real function scntunifrm(radang)
      implicit none
      real x, radang
      x = radang
      scntunifrm = 1.
      return
      end
      
