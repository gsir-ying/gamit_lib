      subroutine timinc(jd,sec,delt)
c
      implicit none

      integer*4 jd,inc,nss
      real*8 sec,delt              
cd     real*8 test,test2
      data nss/86400/
                          
cd      print *,'TIMINC sec delt ',sec,delt
      sec = sec + delt                   
cd      print *,'  sec ',sec 
cd      test = sec/nss + 1.d-15
      inc = int(sec/nss+1.d-15)     
cd      inc = dint(test)  
cd      print *,'test ',test
cd      print *,'inc ',inc
      sec = sec - inc*nss  
cd      test2 = inc*nss
cd      print *,'test2 sec ',test2,sec
      jd = jd + inc
      if (sec .ge. 0) return
      jd = jd - 1
      sec = nss + sec
      return
      end

