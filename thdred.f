
      Subroutine thdred( iut,iscrn,iprnt,nsat,gnss,itsat,satnam
     .                 , jdb,tb,jdf,tf,sdelt,nepcht,jde,te
     .                 , nics,tsatic,nintrs,icsnam
     .                 , precmod,nutmod,gravmod,frame,srpmod,eradmod
     .                 , antradmod )
C
C     Read the headers on the input ephemeris (T-) file

c     Input:

C       IUT    = unit number for T-file
C       ISCRN  = unit number for screen
C       IPRNT  = unit number for print
c       frame = frame of tfile ( 'J2000', 'B1950', 'EFIXD', 'INERT' or 'UNKWN')
c     Output:

C       NSAT   = number of satellites on T-file
c       gnss   = 1-character GNSS designation (G R E C J I)
C       ITSAT  = PRN numbers of satellites on T-file
c       satnam = 16-character satellite name 
C       jdb    = start (PEP) Julian Day of ephemeris
C       tb     = start seconds of day of ephemeris
C       jdf    = stop JD
C       tf     = stop seconds of day
C       SDELT  = tabular interval of ephemeris in seconds
c       NEPCHT = number of epochs on ephemeris
c       JDE    = (PEP) JD of ephemeris epoch
c       TE     = seconds of day of ephemeris epoch
c       NICS   = number of initial conditions plus non-grav parameters
c       TSATIC = initial conditions for satellites
c       NINTRS = number of elements on T-file (3=pos only, 6=pos&vel, 30,48=partials)
c       ICSNAM  = the names of the tfile parameters (IC's)
c       precmod = precession model used on tfile
c       nutmod = nutation model
c       gravmod = gravity model
c       eradmod = Earth radiation model
c       antradmod = antenna radation model

      implicit none

      include '../includes/dimpar.h'

      character*1  earthf,inert,yawbias,buf1,gnss
      character*2  buf2,prn,navstr,ns
      character*4  icsnam(maxorb),icsdef(19)
      character*5 precmod,nutmod,gravmod,frame,srpmod,eradmod,antradmod
      character*16 satnam(maxsat),buf16
      character*80 head,comt(3),prog_name
      character*256 message

      integer*4 julday,jde,jdb,jdf,ite(3),itb(3),itstp(3)
     .        , nepcht,itsat(maxsat),nics,isat,nsat,insn,nintrs
     .        , iscrn,iprnt,iut,id,iday,im,iyr,ihr,imin,iblock,i,j,k
     .        , nclen 

      real*8 tsatic(maxorb,maxsat),tee(3),tbb(3),tstp(3),tb,te,tf
* MOD TAH 200126: Removed svantdx(3) because not used.
     .     , sdelt,sbmass,yawrate,utcoff,taiutc
      real*8 antpwr   ! Antenna transmit power (W)

      integer*4 len,rcpar

c     function
      integer*4 idoy

* MOD TAH 190701: Added svnstart,svnstop which were missing from 
*     svnav_read call.  (Presumably not needed here)
      integer*4 svnstart(5),svnstop(5)  ! Start and stop times from 
                        ! satellite meta data returned from svnav_read.
                        ! Year, DOY, H, M, S (second as integer)

      data icsdef/'X   ','Y   ','Z   ','XDOT','YDOT','ZDOT'
     .           ,'RAD1','RAD2','RAD3','DCOS','DSIN','YCOS'
     .           ,'YSIN','BCOS','BSIN','D2CS','D2SN','D4CS','D4SN'/

      data prn,navstr,ns/'PR','NA','NS'/
      data earthf,inert/'E','I'/

      call uppers(prn)
      call uppers(navstr)
      call uppers(ns)
      call uppers(earthf)
      call uppers(inert)

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)
        
      do 5 I=1,maxsat
5     itsat(I)=0
C
C
C        Read the first header record of the T-file
C
      rewind iut
      read(iut) head,ite,tee,itb,tbb,itstp,tstp,sdelt,nepcht
      call check_y2k(ite(3))
      call check_y2k(itb(3))
      call check_y2k(itstp(3)) 


      write (buf2,'(a2)') head(79:80)
      read  (buf2,'(i2)') nics
      if( nics.gt.maxorb ) then
         write(message,'(a,i2,a,i2)') 'Parameter number on T-file ('
     .         ,nics ,') exceeds MAXORB (',maxorb,')'
         call report_stat('FATAL',prog_name,'lib/thdred',' '
     .                      ,message,0)
       endif
C     for old T files without nics in the header
      if (nics.eq.0) nics=6
c **  temporary to handle NOAA t-files
      if (nics.eq.8) nics = 9

C        Convert the start, stop times to JD, seconds of day
           
      iyr= ite(3)
      im=ite(1)
      id=ite(2)
      ihr = tee(1)
      imin = tee(2)    
cd      print *,'THDRED iyr im id ihr ',iyr,im,id,ihr
cd      print *,'   itb ',itb 
cd      print *,'   itstp ',itstp 
      jde= julday(im,id,iyr)
      jdb= julday(itb(1),itb(2),itb(3))
      jdf= julday(itstp(1),itstp(2),itstp(3))    
      te= 3600.d0*tee(1) + 60.d0*tee(2) + tee(3)
      tb= 3600.d0*tbb(1) + 60.d0*tbb(2) + tbb(3)
      tf= 3600.d0*tstp(1) + 60.d0*tstp(2) + tstp(3)

C        Read the second header record of the T-file

      read(iut) comt,nsat,nintrs
     .        , (satnam(isat),(tsatic(i,isat),i=1,nics),isat=1,nsat)
                       
      buf1= comt(1)(21:21)
      call uppers(buf1)
      comt(1)(21:21) = buf1

c        Extract the orbit IC parameters from the comment line if they exist.
                     
      nclen = nics
      if( nics.gt.15) nclen = 15
      k = 1
      do j = 1,nclen
        if ( comt(2)(k:k+3) .eq. '    ') then
          icsnam(j) = icsdef(j)
          comt(2)(k:k+3) = icsdef(j)
        else
          icsnam(j) = comt(2)(k:k+3)
        endif
        k = k+5
      enddo    
      if( nics.gt.15) then 
        k=1
        do j=16,nics
          if( comt(3)(k:k+3).eq.'    ') then
            icsnam(j) = icsdef(j)
            comt(3)(k:k+3) = icsdef(j)
          else
            icsnam(j) = comt(3)(k:k+3)
          endif
          k= k + 5 
        enddo 
      endif 
      do isat = 1, nsat 
         call uppers(satnam(isat))
      enddo

c        Convert the times to GPST if necessary

      if( comt(1)(11:14).ne.'GPST' ) then
         utcoff = taiutc(jde) - 19.d0
         call timinc(jde,te,utcoff)
         call timinc(jdb,tb,utcoff)
         call timinc(jdf,tf,utcoff)
      endif

c  PT970321: rearrange the checks here to use the "frame" variable
c            to determine/check the ref. frame of tfile  
c  first check inertial vs efixed
      if((frame.eq.'J2000'.or.frame.eq.'B1950'.or.frame.eq.'INERT').and.
     .    comt(1)(21:21).eq.earthf.or
     .              .frame.eq.'EFIXD'.and.comt(1)(21:21).eq.inert)then
          write(message,10) comt(1)(21:35),frame
10        format('Reference System for T-file is ',a15,
     1      ' but expected system is ',a5)
           call report_stat('FATAL',prog_name,'lib/thdred',' '
     .                   ,message,0)     
      endif 

c now check cases and assign more specific values for "frame" where appropriate
      if(frame.eq.'INERT')then
        if( index(comt(1)(30:40),'5').ne.0)  then
          frame='B1950'
        elseif( index(comt(1)(30:40),'2').ne.0)  then
          frame='J2000'
        elseif(comt(1)(36:40).eq.'     ' ) then
          frame='B1950'
        else
          write(message,'(a,1x,a5)')
     .             'THDRED: Unknown inertial frame in T-file header: '
     .             ,comt(1)(36:40)
           call report_stat('FATAL',prog_name,'lib/thdred',' '
     .                   , message,0)
        endif
      elseif(frame.eq.'J2000')then     
        if( index(comt(1)(30:40),'2').eq.0)  then
          write(message,'(a,1x,a5,a,a5)')
     .             'THDRED: T-file frame is: ',comt(1)(36:40)
     .            ,' but expected frame is ',frame
           call report_stat('FATAL',prog_name,'lib/thdred',' '
     .                   , message,0)
        endif
      elseif(frame.eq.'B1950')then
        if( index(comt(1)(30:40),'5').eq.0)  then
          write(message,'(a,1x,a5,a,a5)')
     .             'THDRED: T-file frame is: ',comt(1)(36:40)
     .            ,' but expected frame is ',frame
           call report_stat('FATAL',prog_name,'lib/thdred',' '
     .                   , message,0)
        endif
c PT95072: if there is no frame passed into thdred then assign B1950
c           as the default frame     
c PT970321: there will now always be a value passed in to thdred. If it is
c           EFIXED then write message that no inertial frame will be assigned.
c       if(frame.eq.'     ')then
c rwk 970529: Skip this message since it's not really an error
c      elseif(frame.eq.'EFIXD')then
c         write(message,'(2a)')'Earth-fixed T-file - no inertial frame '
c     .,         'assigned'
c         call report_stat('WARNING',prog_name,'lib/thdred',' '
c     .                   , message,0)


c assign a precession model to EFIXD tfile
c PT970321: Why do we do this? Is it needed? Is there a need to define an
c           inertial frame as well?
        if(comt(1)(42:46).eq.'IAU68')  then
          precmod = 'IAU68'
        elseif(comt(1)(42:46).eq.'IAU76')  then
          precmod = 'IAU76'
        endif  

c  assign a value to "frame" if 'UNKWN' was passed in
      elseif(frame.eq.'UNKWN')then
        if(comt(1)(21:21).eq.earthf)then
          frame = 'EFIXD'
        else
          if( index(comt(1)(30:40),'5').ne.0)  then
            frame='B1950'
          elseif( index(comt(1)(30:40),'2').ne.0)  then
            frame='J2000' 
          elseif(frame.eq.'          ') then
            frame = 'B1950'
          endif
        endif
      endif

c PT 950415: read precession, nutation and gravity models used
      precmod = comt(1)(42:46)
      if(precmod.eq.'     '.and.frame.eq.'B1950') precmod = 'IAU68'
      if(precmod.eq.'     '.and.frame.eq.'J2000') precmod = 'IAU76'
      nutmod = comt(1)(48:52)
      gravmod = comt(1)(54:58)
c      write(*,4276)precmod,nutmod,gravmod
c4276  format(1x,' prec model ',a5,' nut model ',a5,' grav model ',a5)


c        Set model for radiation pressure and non-gravitational forces

      srpmod = comt(1)(60:64)
      if( srpmod.eq.'     ') srpmod = 'SPHRC'
      call uppers(srpmod) 
      eradmod = comt(1)(66:70)  
      if( eradmod.eq.'     ') eradmod = 'NONE '
      call uppers(eradmod)
      antradmod = comt(1)(72:76)              
      if( antradmod.eq.'     ')  antradmod = 'NONE '
      call uppers(antradmod)

c       Determine PRN numbers from the T-file labels

      do 50 isat=1,nsat
         write(buf16,'(a16)') satnam(isat)   
         if( satnam(isat)(1:1) .eq. 'G' .or.
     .       satnam(isat)(1:1) .eq. 'R' .or.
     .       satnam(isat)(1:1) .eq. 'E' .or.
     .       satnam(isat)(1:1) .eq. 'C' .or.
     .       satnam(isat)(1:1) .eq. 'J' .or.
     .       satnam(isat)(1:1) .eq. 'I' ) then
c            label is 1-char GNSS system + PN
             read(buf16,'(a1,i2)') gnss,itsat(isat)
         elseif( satnam(isat)(1:2) .eq. prn ) then
C           Label is 'PRN nn'
            gnss = 'G'
            read(buf16,'(4x,i2)') itsat(isat)
         else    
            gnss = 'G'
            if( satnam(isat)(1:2) .eq. navstr ) then
C              Label is 'NAVSTAR nn'
               read(buf16,'(8x,i2)') insn
            elseif ( satnam(isat)(1:2) .eq. ns ) then
C              Label is 'NSnn'  
               read(buf16,'(2x,i2)') insn
c              need to allow for left-justified names
               if( satnam(isat)(4:4).eq.' ') then
                  if( insn.gt.20 ) insn= insn/10
               endif
            else
              write(message,96) satnam(isat)
96            format('THDRED: satnam=',a16,' not allowed')
              call report_stat('FATAL',prog_name,'lib/thdred',' '
     .                        ,message,0)
            endif 
            iday = idoy(iyr,im,id) 
* MOD TAH 190701: Added to two missing arguments from the call svnstart,svnstop
* MOD TAH 190702: Added antpwr to snav_read call

cwanghutongji20230304
c            call svnav_read( 1,iyr,iday,ihr,imin
c     .                     , insn, itsat(isat)
c     .                     , iblock
c     .                     , sbmass, yawbias, yawrate, antpwr,
c     .                       svnstart, svnstop )

            call svnav_read_old( 1,iyr,iday,ihr,imin, gnss
     .                     , insn, itsat(isat)
c    .                     , iblock
     .                     , sbmass, yawbias, yawrate, antpwr,
     .                       svnstart, svnstop )
cwanghutongji20230304

cccccwanghutongji20170425
c            call svnav_read( 1,iyr,iday,ihr,imin
c     .                     , insn, itsat(isat)
c     .                     , iblock
c     .                     , sbmass,yawbias,yawrate )
cccccwanghutongji20170425
c            call svnav_read(1,iyr,iday,ihr,imin,gnss
c     .	                    ,  insn, itsat(isat)
c     .
c     .                     ,  sbmass,yawbias,yawrate) 
         endif
50    continue
C        Print the header information

      if( iprnt.gt.0 ) then
* MOD TAH 190704: Add comt(3)  'D2CS D2SN D4CS D4SN' string comt(2)
*       so that the Orbital parameter types are correctly report
        write( iprnt,70) head,ite,tee,itb,tbb,itstp,tstp,sdelt,
     .                nepcht,comt,nsat,nintrs,
     .                trim(comt(2)) // ' ' // trim(comt(3))
   70   format(//,
     .        '-------------------------------------------------------',
     .        /,
     .        ' ** T-FILE HEADER INFORMATION **',//,1x,a80,//,
     .        ' Epoch of initial conditions:',3i5,2f4.0,f7.3,/,
     .        ' Ephemeris start            :',3i5,2f4.0,f7.3,/,
     .        ' Ephemeris end              :',3i5,2f4.0,f7.3,/,
     .        ' Tabular interval (sec)     :',f10.3,/,
     .        ' T-file epochs              :',i5,/,
     .        ' Comments                   :',/,3(/,1x,a80),/,
     .        ' Number of satellites       :',1x,i2,/,
     .        ' NINTRS                     :',i3,/,
     .        ' Orbital parameter types    :',1x,a80,//,
     .        ' Satellite name and ics     :')
        do isat = 1,nsat
          write( iprnt,80) satnam(isat),(tsatic(i,isat),i=1,nics)
        enddo
* MOD TAH 190704: Updated format to allow for 13 radition parameters
*       in ECOMC (used to be max 9).  
   80     format(1x,a16,1x,3f14.6,3f8.4,13f9.4)
      endif
C-------------------------------------------------------------
      return
      end

cwanghutongji20230304
c      Subroutine svnav_read( idir,iyr,iday,ihr,imin, gnss,
c     .             isat ,jsat,
c     .             frqchn,antbody,
c     .             sbmass,yawbias,yawrate, antpwr,
c     .             start,stop ) 

      Subroutine svnav_read_old( idir,iyr,iday,ihr,imin, gnss,
     .             isat ,jsat,
c    .             frqchn,antbody,
     .             sbmass,yawbias,yawrate, antpwr,
     .             start,stop ) 
cwanghutongji20230304

c
c PURPOSE: Determine the GNSS SVN from PRN or vica versa, and the antenna/body
c          -type, frequency (for GLONASS), satellite mass, maximum yaw rate, 
c           and sun sensor bias status.
c
c PARAMETERS:                        
c         IN: idir    : idir = 1 SVN to PRN, idir = -1 PRN to SVN  I*4
c             iyr     : 4-digit year                               I*4
c             iday    : 3-digit day of year                        I*4
c             ihr     : 2-digit hr                                 I*4
c             imin    : 2-digit minute                             I*4
c             gnss    : 1-character GNSS code (G R E C J I)        C*1
c             isat    : input SVN or PRN number                    I*4
c
c        OUT: jsat    : output SVN or PRN number                   I*4
c             frqchn  : frquency channel (GLONASS only)            I*4
c             antbody : antenna/body-type (rcvant_tab/ANTEX std)   C*20
c             sbmass  : S/C mass in grams (if old-style will be 0) R*8
c             yawbias : Yaw bias                                   C*1
c             yawrate : maximum yaw rate of the input prn or svn   R*8 
c             antpwr  : Transmission power from igs_metadata       R*8
c
c
c There are three formats possible:
c   Old-style  table with header begining '  NSN/PRN', no longer supported
c   New-style table with header 'svnav.dat  Version  2.0'
c   IGS Satellite Metadata SINEX with header beginning '%=SNX' 
c Since there are no new-style tables beyond version 2.0, we'll use 
c the version variable for subsequent release of the SINEX file 

c The new-style format is:

c svnav.dat  Version  2.0         
c SYS SVN  PRN CHAN ANT/BODY               MASS(G) YAW BIAS  YAW RATE  START           STOP             COMMENTS
c  G     1   4   0  BLOCK I                 453800.     U      0.1999  1978  53  0  0  2100   1  0  0                                           
c  R   701   6  99  GLONASS-M              1450000.     U      0.2500  2003 344  0  0  2010 118  0  0   #channel 99 indicates unknown  
c  R   731  22  -3  GLONASS-M              1450000.     U      0.2500  2010  60  0  0  2100   1  0  0
c  D     1  30   0  BEIDOU-2M              1000000.     Y      0.3000  2007 103  0  0  2100   1  0  0  #mass unknown
c  E   101  11   0  GALILEO-1               700000.     Y      0.4000  2011 294  0  0  2100   1  0  0
 
c Non-blank first column is a header or comment.   
                                                                                 
c NOTE 1: Since ARC integrations using cross a day boundary to provide interpolation,
c         there is a logical problem if the SVN/PRN assignment changes on the day
c         boundary. Since we'll want to keep the assignment that corresponds to 
c         the observation day (IC epoch), if a subsequent entry is found, we'll ignore 
c         the change.

* MOD TAH 190702: Added return of transmission power from igs_metadata.snx.  If
*     svnav.dat linked, 0.0 is returned since value is not known.

      implicit none

                  
      character*1  gnss,gnss_tmp,yawbias         
      character*10 svid 
      character*20 antbody 
      character*80 prog_name   
      character*128 record
      character*256 message

      integer*4 lun,idir,isat,jsat,frqchn
     .         ,iyr,iday,ihr,imin,time(5),start(5),stop(5),power
     .         ,iprn_tmp,isvn_tmp,len,rcpar,itemp,ioerr,i
c
      real*8 sbmass,yawrate,version
      real*8 antpwr    ! Tranmission power (W).
                    
      logical eof,found   

* MOD TAH 190627: Cleaned up code indenation and made use of svnav.dat
*     or igs_metadata_WWWW.snx more explicit.
      logical svnav_file  ! Set true is svnav.dat linked to svnav.dat
                          ! as opposed to igs_metadata_WWWW.snx file.  

c     function
      integer*8 itimdif

c Get calling program name for report_stat
      len = rcpar(0,prog_name)
                  
cd      print *,'SVNAV_READ idir iyr iday ihr imin isat '
cd     .      ,  idir,gnss,iyr,iday,ihr,imin,isat  

c Put the requested time into an array for comparing with file entires
      time(1) = iyr
      time(2) = iday
      time(3) = ihr
      time(4) = imin
      time(5) = 0                            
c     times read from file have no seconds
      start(5) = 0
      stop(5) = 0 

c Open the svnav.dat file
      lun = 69 
      open(unit=lun,file='svnav.dat',form='formatted',status='old'
     .    ,iostat=ioerr ) 
      If( ioerr.ne.0 )  call report_stat('FATAL',prog_name
     .   ,'lib/svnav_read','svnav.dat','File not found ',ioerr)
                                       
c Get the version number from the first head line
      read(lun,'(a)') record               
cd      print *,'1st record ',record 
      if( record(3:5).eq.'NSN') then
        call report_stat( 'FATAL',prog_name,'lib/svnav_read'
     .    ,'svnav.dat','Old version of file not supported',0)
      elseif( record(12:18).eq.'Version') then
        read( record(20:23),'(f4.0)') version
c       this must be 2.0 (Make explicit test)
* MOD TAH 190627: Make explicit this is svnav.dat
        if( version.ne.2.0 ) then
           call report_stat('FATAL',prog_name,'lib/snvav_read',
     .          'svnav.dat','Wrong version of snvav.dat ',
     .          int(version*10))
        else
          svnav_file = .true.
        endif 
      elseif( record(3:5).eq.'SNX') then
        read( record(7:10),'(f4.0)')  version 
* MOD TAH 190627: Make explicit this is not svnav.dat
        svnav_file = .false.
      else
        call report_stat( 'FATAL',prog_name, 'lib/svnav_read',
     .               'svnav.dat','Unknown version of file ',0)
      endif
          
c See if a GAMIT-style svnav.dat file
* MOD TAH 190627: Explicitly test for svnav.dat (old 
C     if( version.eq.2.0 ) then 
      if( svnav_file ) then 

c        Skip the column head line
         read(lun,'(a)') record
         if (ioerr .ne. 0) 
     .   call report_stat('FATAL',prog_name,'lib/snvav_read',
     .             'svnav.dat','Error reading svnav.dat file',ioerr)
                       

c        Loop over svnav.dat file to find the requested entry
         eof = .false.    
         found = .false.
         do while (.not.found .and. .not.eof )
           read(lun,'(a)',iostat=ioerr) record   
           if( ioerr.eq.-1 ) then
             eof = .true.
           elseif( ioerr.ne.0 ) then
              call report_stat('FATAL',prog_name,'lib/snvav_read'
     .         ,'svnav.dat','Error reading svnav.dat file',ioerr)
           elseif (record(1:1).ne.' ') then
c            comment 
             continue 
           else
             read(record,'(1x,a1,2x,i4,2x,i2,2x,i2,2x,a20,f11.0,5x
     .              , a1,f12.4,2(i6,i4,2i3))')  
     .         gnss_tmp,isvn_tmp,iprn_tmp,frqchn,antbody,sbmass
     .       , yawbias, yawrate
     .       , (start(i),i=1,4),(stop(i),i=1,4)

c            check if prn or svn is being requested
             if ( idir.lt.0 ) then  
c              SVN from PRN   
               if ( gnss_tmp.eq.gnss .and. iprn_tmp.eq.isat ) then
                 if( itimdif(start,time).le.0 .and. 
     .            itimdif(stop,time).ge.0 ) then           
                   found = .true.
cd                   print *,'found ',found 
                   jsat = isvn_tmp
                 endif
               endif
             elseif (idir.gt.0 ) then 
c              PRN from SVN
c              Beidou does not have official SVNS but svnav.dat has 
c              MGEX-assigned values
c              if( gnss.eq.'C' )  call report_stat('FATAL',prog_name
c    .   ,'lib/svnav_read','svnav.dat','Beidou has no SVN defined',0)
               if ( gnss_tmp.eq.gnss .and. isvn_tmp.eq.isat ) then
                 if( itimdif(start,time).ge.0 .and. 
     .            itimdif(stop,time).le.0 ) then
                   found = .true.
                   jsat = iprn_tmp
                 endif
               endif
             else
               call report_stat('FATAL',prog_name,'lib/svnav_read',' '
     .                      ,'idir=0, something wrong',0)
             endif
           endif
         enddo
                                   
c Stop if SV not found in svnav.dat file
         if( .not.found ) then 
           if( idir.lt.0 ) then 
             write(message,'(a,a1,i2,a)') 'PRN ',gnss,isat
     .           ,' not found on svnav.dat file'
           elseif( idir.gt.0 ) then
             write(message,'(a,a1,i3,a)') 'SVN  ',gnss,isat
     .           ,' not found on svnav.dat file'
           endif             
cd           print *,'prog_name(1:7)',prog_name(1:7)
           if( prog_name(1:7).eq.'dcbtab2' .or. 
     .         prog_name(1:3).eq.'arc') then
               call report_stat('WARNING',prog_name,'lib/svnav_read',
     .                     ' ', message,0)
           else
               call report_stat('FATAL',prog_name,'lib/svnav_read',
     .                    ' ', message,0)
           endif
         endif  

* MOD TAH 1900702: Set antpwr to 0 since not in svnav.dat
         antpwr = 0.d0

c If not a GAMIT-style file, read a SINEX files
  
      ELSE 
        call read_svsinex( lun,idir,iyr,iday,ihr,imin,gnss,isat
     .                   , jsat,frqchn,antbody,svid,sbmass
     .                   , yawbias,yawrate,power ) 
*       Convert power I4 to R8 antpwr
        antpwr = power
      endif
      close(lun)    

c     DEBUG
c     stop

      end










