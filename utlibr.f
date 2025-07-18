      SUBROUTINE UTLIBR (RMJD, DUT1, DLOD)
*+
*  - - - - - - - - 
*   U T L I B R
*  - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
**
** Downloaded from ftp:/tai.bipm.org/iers/conventions2010/chapter5
** and added to GAMIT by R. King 1 July 2013.  The only change is
** to append the subroutine FUNDARG.F from convetions/2010/chapter8
** to the bottom of this file.  This is nearly duplicative with 
** gamit/lib/funarg.f and kf/gen_util/fund_arg.f, which have 
** different calling arguments, but to stay consistent with the IERS 
** and to avoid some extra checking (needed eventually), we'll not 
** now unify these routines.
** 
*  This subroutine evaluates the model of subdiurnal libration
*  in the axial component of rotation, expressed by UT1 and LOD.
*  This effect is due to the influence of tidal gravitation on the
*  departures of the Earth's mass distribution from the rotational
*  symmetry, expressed by the non-zonal components of geopotential.
*  The amplitudes have been computed for an elastic Earth with liquid
*  core. The adopted truncation level is 0.033 microseconds in UT1
*  corresponding to the angular displacement of 0.5 microarcseconds
*  or to 0.015 mm at the planet surface. With this truncation level
*  the model contains 11 semidiurnal terms. The coefficients of
*  the model are given in Table 5.1b of the IERS Conventions (2010).
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status:  Class 3 model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as
*     a Class 1, 2, or 3 model.
*
*  Given:
*     rmjd        d      Time expressed as modified Julian date
*
*  Returned:
*     dUT1        d      Incremental UT1 in microseconds
*     dLOD        d      Incremental LOD in microseconds per day
*
*  Notes:
*  1) The procedure FUNDARG.F is the same as used by the program PMSDNUT2.F
*     which implements the corresponding model of the lunisolar libration in
*     polar motion.
*
*  Called:
*     FUNDARG             Compute the angular fundamental arguments
*
*  Test cases:
*     given input:  rmjd_a = 44239.1 ( January 1, 1980 2:24.00 )
*                   rmjd_b = 55227.4 ( January 31, 2010 9:35.59 )
*
*     expected output: dUT1_a =   2.441143834386761746D0 mus;
*                      dLOD_a = -14.78971247349449492D0 mus / day
*                      dUT1_b = - 2.655705844335680244D0 mus;
*                      dLOD_b =  27.39445826599846967D0 mus / day
*
*  References:
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2010 May       A.Brzezinski   Original code
*  2010 June  1   B.E.Stetzler   Initial changes to code
*  2010 June  2   B.E.Stetzler   Provided test case
*  2010 June  2   B.E.Stetzler   Capitalized all variables for FORTRAN
*                                77 compatibility
*  2010 June  2   B.E.Stetzler   Replaced call to PMARGS to FUNDARG
*                                for universal fundamental argument
*                                subroutine
*  2010 June  2   B.E.Stetzler   Validated test case using internally
*                                computed GMST and call to FUNDARG
*                                matched previous external call to
*                                PMARGS for all six parameters
*  2010 June  23  B.E.Stetzler   Modified coefficients of semi-diurnal
*                                variations in UT1 and LOD due to
*                                libration for a non-rigid Earth to
*                                coincide with Table 5.1b
*-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION RMJD, DUT1, DLOD

*         ----------------------------
*           D E F I N I T I O N S
*         ----------------------------
*  iarg   - array defining for each of the 11 trigonometric terms a set
*           of 6 integer multipliers of the fundamental angular arguments
*  arg    - vector of the following 6 fundamental arguments used to
*           compute the angular argument of the trigonometric functions
*           arg(1:6) = [ GMST+pi, el, elp, f, d, om ]; this vector is
*           evaluated by the subroutine FUNDARG which is called as an 
*           external subroutine.  Originally evaluated by the subroutine
*           PMARGS. 
*  period - array of periods of the trigonometric terms of expansion, in
*           mean solar days; only for a check - not used in computations
*  dUT1s, dUT1c - sine and cosine coefficients of dUT1, in microseconds
*  dLODs, dLODc - sine and cosine coefficients of dLOD, in microseconds
*                 per day
*  angle  - angular argument of the trigonometric functions
*           angle = Sum(i=1:6) iarg(i,j)*arg(i), for j=1,11

      INTEGER I, J
      INTEGER IARG(6,11)
      DOUBLE PRECISION T, GMST, L, LP, F, D, OM
      DOUBLE PRECISION ARG(6)
      DOUBLE PRECISION PER(11), DUT1S(11), DUT1C(11), DLODS(11),
     .                 DLODC(11) 
      DOUBLE PRECISION ANGLE

* Set constants

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  rmjd0   - modified Julian date of J2000
*  twopi   - 2*pi

      DOUBLE PRECISION RMJD0, PI, TWOPI
      PARAMETER ( RMJD0   = 51544.5D0                )
      PARAMETER ( PI      = 3.141592653589793238462643D0 )
      PARAMETER ( TWOPI   = 6.283185307179586476925287D0 )

*  Radians to seconds
      DOUBLE PRECISION RAD2SEC
      PARAMETER ( RAD2SEC = 86400D0/TWOPI            )

* Coefficients of the quasi semidiurnal terms in dUT1, dLOD 
* Source: IERS Conventions (2010), Table 5.1b

      DATA 
     .((IARG(I,J),I=1,6), PER(J), DUT1S(J),DUT1C(J), DLODS(J), DLODC(J),
     .                                                           J=1,11)
     ./2, -2,  0, -2,  0, -2, 0.5377239,  0.05, -0.03,  -0.3,  -0.6,
     . 2,  0,  0, -2, -2, -2, 0.5363232,  0.06, -0.03,  -0.4,  -0.7,
     . 2, -1,  0, -2,  0, -2, 0.5274312,  0.35, -0.20,  -2.4,  -4.1,
     . 2,  1,  0, -2, -2, -2, 0.5260835,  0.07, -0.04,  -0.5,  -0.8,
     . 2,  0,  0, -2,  0, -1, 0.5175645, -0.07,  0.04,   0.5,   0.8,
     . 2,  0,  0, -2,  0, -2, 0.5175251,  1.75, -1.01, -12.2, -21.3,
     . 2,  1,  0, -2,  0, -2, 0.5079842, -0.05,  0.03,   0.3,   0.6,
     . 2,  0, -1, -2,  2, -2, 0.5006854,  0.04, -0.03,  -0.3,  -0.6,
     . 2,  0,  0, -2,  2, -2, 0.5000000,  0.76, -0.44,  -5.5,  -9.6,
     . 2,  0,  0,  0,  0,  0, 0.4986348,  0.21, -0.12,  -1.5,  -2.6,
     . 2,  0,  0,  0,  0, -1, 0.4985982,  0.06, -0.04,  -0.4,  -0.8/

* Compute the harmonic model of dUT1 and dLOD 
* dUT1 and dLOD are set to zero first 
      DUT1 = 0D0
      DLOD = 0D0

* Evaluate the vector of the fundamental arguments
* arg(1:6) = [ GMST+pi, el, elp, f, d, om ] at t = rmjd

*  Convert the input epoch to Julian centuries of TDB since J2000
      T = (RMJD-RMJD0)/36525D0

*  Compute GMST + pi
      GMST = MOD (   67310.54841D0 +
     .               T*( (8640184.812866D0 + 3155760000D0) +
     .               T*( 0.093104D0 +
     .               T*( -0.0000062 ))), 86400D0 )

      CALL FUNDARG ( T, L, LP, F, D, OM )

      ARG(1) = GMST / RAD2SEC + PI
      ARG(1) = DMOD( ARG(1), TWOPI )
      ARG(2) = L
      ARG(3) = LP
      ARG(4) = F
      ARG(5) = D
      ARG(6) = OM 

      DO 20 J=1,11

* For the j-th term of the trigonometric expansion, compute the angular
* argument angle of sine and cosine functions as a linear integer
* combination of the 6 fundamental arguments
        ANGLE = 0D0
        DO 10 I=1,6
          ANGLE = ANGLE + IARG(I,J) * ARG(I)
   10   CONTINUE
        ANGLE = DMOD( ANGLE, TWOPI )

* Compute contribution from the j-th term of expansion to dUT1 and dLOD 
        DUT1 = DUT1 + DUT1S(J)*DSIN(ANGLE) + DUT1C(J)*DCOS(ANGLE)
        DLOD = DLOD + DLODS(J)*DSIN(ANGLE) + DLODC(J)*DCOS(ANGLE)
   20 CONTINUE
      RETURN

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*
*-----------------------------------------------------------------------
      END
      
*************************************************************************

      SUBROUTINE FUNDARG ( T, L, LP, F, D, OM )
*+
*  - - - - - - - - - - -
*   F U N D A R G 
*  - - - - - - - - - - -
*
*  This routine is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This subroutine computes the lunisolar fundamental arguments.
*  The model used is from Simon et al. (1994) as recommended by the IERS
*  Conventions (2010).  Refer to IERS Conventions (2010) Chapter 5 
*  Sections 5.7.1 - 5.7.2 (pp. 57 - 59).
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
* 
*  Status: Canonical model
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as a
*     Class 1, 2, or 3 model.
*
*  Given:
*     T           d      TT, Julian centuries since J2000 (Note 1)
*
*  Returned:
*     L           d      Mean anomaly of the Moon (Note 2)
*     LP          d      Mean anomaly of the Sun (Note 2)
*     F           d      L - OM (Notes 2 and 3)
*     D           d      Mean elongation of the Moon from the Sun
*                                                         (Note 2)
*     OM          d      Mean longitude of the ascending node of
*                                                the Moon (Note 2)
*
*  Notes:
*
*  1) Though T is strictly TDB, it is usually more convenient to use
*     TT, which makes no significant difference.  Julian centuries since
*     J2000 is (JD - 2451545.0)/36525.
*
*  2) The expression used is as adopted in IERS Conventions (2010) and
*     is from Simon et al. (1994).  Arguments are in radians.
*
*  3) L in this instance is the Mean Longitude of the Moon. OM is the 
*     Mean longitude of the ascending node of the Moon.
*
*  Test case:
*     given input: T = 0.07995893223819302 Julian centuries since J2000
*                  (MJD = 54465)
*     expected output:  L = 2.291187512612069099 radians
*                       LP = 6.212931111003726414 radians
*                       F = 3.658025792050572989 radians
*                       D = 4.554139562402433228 radians
*                       OM = -0.5167379217231804489 radians
*
*  References:
*
*     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*     Francou, G., Laskar, J., 1994, Astron.Astrophys. 282, 663-683
*
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2008 January 18 B.E.Stetzler  Initial changes to header
*               and used 2PI instead of PI as parameter
*  2008 January 25 B.E. Stetzler Additional changes to header
*               and defined fundamental arguments
*  2008 January 28 B.E. Stetzler Additional changes to header
*  2008 March   12 B.E. Stetzler Applied changes to wording of notes.
*  2008 April   03 B.E. Stetzler Provided example test case.
*  2009 February 11 B.E. Stetzler Corrected term in OM from 6962890.2665
*                                 to 6962890.5431 and updated test case
*  2009 May     07 B.E. Stetzler Code formatting changes based on 
*                                client recommendations
*  2009 May     07 B.E. Stetzler Updated test case due to above changes
*  2010 February 25 B.E. Stetzler Recalculation of fundamental arguments
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION T, L, LP, F, D, OM

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Arcseconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Compute the fundamental argument L.
      L = MOD (       485868.249036D0 +
     .                  T*( 1717915923.2178D0 +
     .                  T*(         31.8792D0 +
     .                  T*(          0.051635D0 +
     .                  T*(        - 0.00024470D0 )))), TURNAS ) * DAS2R

*  Compute the fundamental argument LP.
      LP = MOD (       1287104.79305D0 +
     .            T*( 129596581.0481D0 +
     .            T*(       - 0.5532D0 +
     .            T*(         0.000136D0 +
     .            T*(       - 0.00001149D0 )))), TURNAS ) * DAS2R

*  Compute the fundamental argument F.
      F  = MOD (       335779.526232D0 +
     .                  T*( 1739527262.8478D0 +
     .                  T*(       - 12.7512D0 +
     .                  T*(       -  0.001037D0 +
     .                  T*(          0.00000417D0 )))), TURNAS ) * DAS2R

*  Compute the fundamental argument D.
      D = MOD (        1072260.70369D0 +
     .          T*( 1602961601.2090D0 +
     .          T*(        - 6.3706D0 +
     .          T*(          0.006593D0 +
     .          T*(        - 0.00003169D0 )))), TURNAS ) * DAS2R

*  Compute the fundamental argument OM.
      OM = MOD (       450160.398036D0 +
     .             T*( - 6962890.5431D0 +
     .             T*(         7.4722D0 +
     .             T*(         0.007702D0 +
     .             T*(       - 0.00005939D0 )))), TURNAS ) * DAS2R

*  Finished.

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
* 
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution. 
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*
*-----------------------------------------------------------------------
      END



