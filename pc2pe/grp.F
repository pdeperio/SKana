***   ==============                                                 ***
      PROGRAM XEMEAN
***   ==============                                                 ***
C
C     Cut from /home/iida/calib/pc2pe/pc2pe_08oct/grp/grp.F
C
      IMPLICIT NONE

#include "skparm.h"
#include "skhead.h"
#include "sktq.h"
#include "skpdst.h"
#include "geopmt.h"      
#include "skveto.h"
#include "geotnk.h"

*     LUN
      INTEGER LUNGRP
      character*10 ser,address2
      
      PARAMETER (LUNGRP   = 31)

      REAL BALLX,BALLY,BALLZ,R
      integer I, J, K, grp, ipmt, KAI, REM
    
      SK_GEOMETRY=4
      call geoset

      open(LUNGRP,FILE='pmt_groups.txt')


************************ grouping   ****************************************
      DO K = 1,7
         
         write(*,*)'K=',K
 
         REM = 0
         KAI = 0
         
         DO ipmt = 1,MAXPM

            IF(KIJPM(1,ipmt).eq.1)THEN ! barrel
               REM=mod(ipmt,51)
               IF(REM.eq.0)THEN
                  REM=51
               END IF 
               grp=int(real(REM-1)/3.)+1

               write(LUNGRP,*) ipmt,grp

            ELSE IF(KIJPM(1,ipmt).eq.2)THEN ! top
               r=sqrt((xyzpm(1,ipmt)-ballx*100.)**2+(xyzpm(2,ipmt)-bally*100.)**2)               

               if(r.gt.RINTK) then ! oosoto        
                  
                  write(LUNGRP,*) ipmt,25
                   
                  GOTO 20
                  
               endif
               
               DO I=1,8
                  J=I+17
                  
                  if(r.le.I*RINTK/8) then ! group I
                    
                     write(LUNGRP,*) ipmt,J

                     GOTO 20
                  end if
               END DO   
               
            ELSE IF(KIJPM(1,ipmt).eq.3)THEN ! bottom   
               r=sqrt((xyzpm(1,ipmt)-ballx*100.)**2+(xyzpm(2,ipmt)-bally*100.)**2)
               
               if(r.gt.RINTK) then ! oosoto
                  
                  write(LUNGRP,*) ipmt,33
                     
                  GOTO 20                  
                  
               endif

               DO I=1,8
                  J=I+25
                  if(r.le.I*RINTK/8) then ! group I

                     write(LUNGRP,*) ipmt,J
                        
                     GOTO 20
                     
                  end if
               END DO               
               
            END IF
 20         CONTINUE  
         END DO
      END DO
   
      CLOSE(LUNGRP)

      stop
      END

      real function scntunifrm(radang)
      implicit none
      real x, radang
      x = radang
      scntunifrm = 1.
      return
      end
    



