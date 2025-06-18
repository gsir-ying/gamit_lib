
      SUBROUTINE DAYJUL( JD,iyear,IDOY )
C     WRITTEN BY R.KING 1 AUG 83
C     GIVEN JD, RETURNS DAY OF YEAR AND YEAR    

      implicit none
      integer*4 jd,jd1900,id,idoy,iyear
      data jd1900/2415020/
C
      ID= (JD - jd1900)*100
      iyear= int(ID/36525) 
      IDOY= ID/100 - iyear*365 - (iyear-1)/4   
      if( iyear.gt.60 ) then
         iyear = iyear + 1900
      else
         iyear = iyear + 2000
      endif  
      RETURN
      END
