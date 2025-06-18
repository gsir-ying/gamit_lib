      subroutine gettim(ihr,imn,isec,ihnsec)
      integer i3(3)
      integer ihr,imn,isec,ihnsec

c     return operating system time in hour,min,sec

      call itime(i3)
      ihr = i3(1)
      imn = i3(2)
      isec = i3(3)
      ihnsec = 0

      return
      end
