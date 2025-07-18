
SUBROUTINE time_TT (mjd , mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC)

! ----------------------------------------------------------------------
! Subroutine:	time_TT.f90
! ----------------------------------------------------------------------
! Purpose:
!  Time Systems transformation
! ----------------------------------------------------------------------
! Input arguments:
! - mjd:			Modified Julian Day number in TT (including fraction of the day)
!
! Output arguments:
! - mjd_TT:			MJD number in TT (including the fraction of the day)
! - mjd_GPS:		MJD number in GPS (including the fraction of the day)
! - mjd_TAI:		MJD number in TAI (including the fraction of the day)
! - mjd_UTC:		MJD number in UTC (including the fraction of the day)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
      IMPLICIT NONE
  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = 8), INTENT(IN) :: mjd
! OUT
      REAL (KIND = 8), INTENT(OUT) :: mjd_TT, mjd_GPS, mjd_TAI, mjd_UTC
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = 8) :: jd0
      INTEGER IY, IM, ID, J_flag
      DOUBLE PRECISION FD  
      DOUBLE PRECISION TAI_UTC, TT_TAI, TAI_GPS
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      jd0 = 2400000.5D0
      CALL iau_JD2CAL ( jd0, mjd, IY, IM, ID, FD, J_flag )
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Leap Seconds: TAI-UTC difference
! ----------------------------------------------------------------------
!SUBROUTINE iau_DAT ( IY, IM, ID, FD, DELTAT, J )
!*     DELTAT   d     TAI minus UTC, seconds
      CALL iau_DAT ( IY, IM, ID, FD, TAI_UTC, J_flag )
! ----------------------------------------------------------------------
      !PRINT *,"TAI_UTC", TAI_UTC

! ----------------------------------------------------------------------
! TT
      mjd_TT = mjd
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! TAI
      TT_TAI = 32.184D0
      mjd_TAI = mjd_TT - TT_TAI / 86400D0
! ----------------------------------------------------------------------
  
! ----------------------------------------------------------------------
! UTC
      mjd_UTC = mjd_TAI - TAI_UTC / 86400D0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! GPS Time scale
      TAI_GPS = 19D0
      mjd_GPS = mjd_TAI - TAI_GPS / 86400D0
! ----------------------------------------------------------------------


END
