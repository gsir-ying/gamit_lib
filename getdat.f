      subroutine getdat(iyear,imonth,iday)

c     return the run time date

      integer date(5)
      integer iyear,imonth,iday
      real*8 secs(2)
      external idate

       call systime (date,secs)
       iyear =  date(1)
       call fix_y2k(iyear)
       imonth = date(2)
       iday = date(3)

      return
      end
