      Subroutine check_y2k(iyear)

c     Check for 2 or 4-digit years.  If 2-digits, convert to 4 by guessing
c     the date.

      implicit none

      integer*4 iyear,len,rcpar

       character*80 prog_name

      if( iyear.le.1900) then 
                  
c        temporary to catch problems
c        get the program name to report errors
         len = rcpar(0,prog_name)
         call report_stat('WARNING',prog_name,'lib/check_y2k','  '
     .                   , 'Unexpected 2-digit year',0)
         if( iyear.gt.60 ) then
            iyear = iyear + 1900
         else
            iyear = iyear + 2000
         endif  

      endif

      return
      end
