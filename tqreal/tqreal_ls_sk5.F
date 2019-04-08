
      program sample

      implicit none
      
      integer i, ierr, icab, nevent, nevent2

#include "skhead.h"
#include "skparm.h"
#include "sktq.h"

      integer lun, iret, narg, iargc
      parameter (lun=10)
      character*200 fname_in, fname_out

      character*200 temp
      integer it,nintoskip,nouttoskip
      parameter (nintoskip=6)
      parameter (nouttoskip=4)
      character*200 in_branches_to_skip(nintoskip)
      character*200 out_branches_to_skip(nouttoskip)

      character SK_GEOMETRY_TMP

      data in_branches_to_skip /'SPACERS','QBEESTATUS','DBSTATUS',
     $     'MISMATCHEDHITS','ATMPD','UPMU'/
      data out_branches_to_skip /'SOFTWARETRG','EVENTTRAILER',
     $     'HWTRGLIST','PEDESTALS'/

      integer num/0/, i15001,it15001,i15012
      real t15001

c*** check arguments
      narg = iargc()
      if (narg .lt. 3) then
         print *, 'Usage: tqreal output_root input_root SK_GEOMETRY'
         call exit(1)
      endif

c*** set SK_GEOMETRY
      call getarg(3,SK_GEOMETRY_TMP)
      read(SK_GEOMETRY_TMP, *) SK_GEOMETRY
      print *, 'SK_GEOMETRY=', SK_GEOMETRY

c***  set table version
      call skcaltabfile('./table_version.99.00')

c*** input file (SKROOT)
c      fname_in='rfm_run062267.000001.root'

c*** get output file name
      call getarg(1, fname_out)
      call skroot_open(lun,fname_out(1:lnblnk(fname_out)))

c*** prepare input files
      call getarg(2, fname_in) 
      call skroot_set_input_file(lun,fname_in(1:lnblnk(fname_in)))

c*** skip any unneeded branches; this must be called before skroot_init
      do it=1,nintoskip
         temp = in_branches_to_skip(it)
         call skroot_zero_branch(lun,0,temp(1:lnblnk(temp)))
      enddo

      do it=1,nouttoskip
         temp = out_branches_to_skip(it)
         call skroot_zero_branch(lun,1,temp(1:lnblnk(temp)))
      enddo

c*** initilalize skroot files
      call skroot_init(lun)

c*** set options
      call skoptn('31,30')

 1000 continue
      CALL SKRAWREAD(lun,*1002,*1001,*1000,*1000)

      if (SK_GEOMETRY.eq.5) then

         CALL SKREAD(-lun,*1002,*1001,*1000,*1000)

         if(idtgsk.ne.16) goto 1000

         i15001=0
         i15012=0
         do i=1,nqiskz
c            print*,i,icabiz(i),tiskz(i),tbuf_raw(i)
            if(icabiz(i).eq.15001) then
c               print*,icabiz(i),itiskz(i),tiskz(i)
               i15001 = i15001+1
               it15001 = itiskz(i)
               t15001 = tiskz(i)-1000.
            endif
            if(icabiz(i).eq.15012) i15012 = i15012+1
         enddo

         if(i15001.ne.2 .or. i15012.ne.1) goto 1000
c         if(i15001.ne.2) goto 1000

         num=num+1
c         call SET_TIMING_GATE(it15001)
         call SET_TIMING_GATE_NSEC(t15001)
c         print*,'NEVSK: ',num, nevsk, nqisk, ntimsk(2), ntimsk(3)
         do i=1,nqiskz
            if(icabiz(i).le.MAXPM) tbuf_raw(i) = tiskz(i)
         enddo

      endif

c***  store header & TQ info.
      call skroot_set_tree(lun)
         
c***  output root file 
      call skroot_fill_tree(lun) ! write to disk

      goto 1000
      
 1002 print *,'Read Error'
 1001 continue

c*** close file
      call skclosef(lun)
      call skroot_end

      stop
      end

