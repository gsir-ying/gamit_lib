      Subroutine readc1 (lunit,ntext,text)

      implicit none

      include '../includes/dimpar.h'

      integer*4          nversn,iflag,ntext
      character*4        buf4
      character*16       cfname
      character*80       text(maxtxt),prog_name
      character*256      message
      integer*4          lunit,ioerr,inqerr,len,rcpar,i

c     Read the first part of a C-file which is record #, version #

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)

c     first read only the first two integers, to make sure we have the right version
                     
      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag,nversn                    
c rwk 060905: check for binary compatibility
      if( iflag.ne.1 ) 
     .   call report_stat('FATAL',prog_name,'lib/readc1',' ',
     .'Error reading first record of c-file; binary incompatibility?',0)
     
c rwk 050131: no longer check for ancient (<1995) versions, but rather for 9.80 vs 10.2.
* MOD TAH 200126: Check to see if 1071 version.
      if( nversn.ne.1071 ) then    
        write(message,'(a,i4,a)') 'Old version of C-file (',nversn,')'
        call report_stat('FATAL',prog_name,'lib/readc1',' ',message,0)
      endif
      rewind ( lunit )

c     now read the first three integers, checking that ntext isn't too big for dimensions

      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag,nversn,ntext
      if( ntext.gt.maxtxt ) then
c     Get the calling module name for report_stat
        len =rcpar(0,prog_name)
        write(message,'(a,i4,a,i4)')
     .    'READC1: ntext (',ntext,') > maxtxt (',maxtxt
        call report_stat('FATAL',prog_name,'lib/readc1',' ',message,0)
      endif
      rewind ( lunit )

c     now read the whole record for real
      read (unit    =  lunit
     .,     iostat  =  ioerr
     .,     end     =  1000
     .,     err     =  1000)
     .      iflag,nversn,ntext,(text(i),i=1,ntext)

 1000 if (ioerr .ne. 0) then
           inquire ( unit=lunit,name=cfname,iostat=inqerr )
           call report_stat('FATAL',prog_name,'lib/readc1',cfname
     .                     ,'Error reading first record of C-file--is it
     . empty from a previous step?',ioerr)
      endif

      if (iflag .ne. 1) then
         write(buf4,'(i4)') iflag
         call report_stat('FATAL',prog_name,'lib/readc1',buf4
     .                   , 'Wrong iflag: ',0)
      endif
c      print *,'READC1: ntext = ',ntext,' text follows: '
c      do 2000 i=1,ntext
c         print *,text(i)
c 2000 continue
      return
      end
