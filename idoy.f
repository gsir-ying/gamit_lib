      integer function idoy (iyear,imonth,iday)
c
c     determines the day of year from the calendar month and day

      integer julday, imonth, iyear,iday

      idoy = julday(imonth,iday,iyear) - julday(1,1,iyear) + 1

      return
      end

