      subroutine gpt (dmjd,dlat,dlon,dhgt,pres,temp,undu)

C     This subroutine determines Global Pressure and Temperature
C     based on Spherical Harmonics up to degree and order 9
C
C     input data
C     ----------
C     dmjd: modified julian date
C     dlat: latitude in radians
C     dlon: longitude in radians
C     dhgt: ellipsoidal height in m
C
C     output data
C     -----------
C     pres: pressure in hPa
C     temp: temperature in Celsius
C     undu: Geoid undulation in m (from a 9x9 EGM based model)
C
C     geoid undulation is accounted for
C
C     Reference: A spherical harmonic expansion 
C     of annual pressure and temperature variations for geodetic applications
c
c
c      implicit double precision (a-h,o-z)
c
c      dimension dfac(20),P(10,10),aP(55),bP(55),
c     .          ap_mean(55),bp_mean(55),ap_amp(55),bp_amp(55),
c     .          at_mean(55),bt_mean(55),at_amp(55),bt_amp(55)

       implicit none
       real*8  dfac(20),P(10,10),aP(55),bP(55),
     .          ap_mean(55),bp_mean(55),ap_amp(55),bp_amp(55)
     .          ,at_mean(55),bt_mean(55),at_amp(55),bt_amp(55)
     .          ,a_geoid(55),b_geoid(55)
     .          ,pi,doy,dmjd,pres,temp,dlat,dlon,dhgt,t,apm,apa
     .          ,pres0,atm,ata,temp0,sum,hort,undu

       integer i,j,k,m,n,ir



      pi = 3.14159265359d0

C     reference day is 28 January
C     this is taken from Niell (1996) to be consistent
      doy = dmjd  - 44239.d0 + 1 - 28

      data (a_geoid(i),i=1,55)/
     .-5.6195d-001,-6.0794d-002,-2.0125d-001,-6.4180d-002,-3.6997d-002,
     .+1.0098d+001,+1.6436d+001,+1.4065d+001,+1.9881d+000,+6.4414d-001,
     .-4.7482d+000,-3.2290d+000,+5.0652d-001,+3.8279d-001,-2.6646d-002,
     .+1.7224d+000,-2.7970d-001,+6.8177d-001,-9.6658d-002,-1.5113d-002,
     .+2.9206d-003,-3.4621d+000,-3.8198d-001,+3.2306d-002,+6.9915d-003,
     .-2.3068d-003,-1.3548d-003,+4.7324d-006,+2.3527d+000,+1.2985d+000,
     .+2.1232d-001,+2.2571d-002,-3.7855d-003,+2.9449d-005,-1.6265d-004,
     .+1.1711d-007,+1.6732d+000,+1.9858d-001,+2.3975d-002,-9.0013d-004,
     .-2.2475d-003,-3.3095d-005,-1.2040d-005,+2.2010d-006,-1.0083d-006,
     .+8.6297d-001,+5.8231d-001,+2.0545d-002,-7.8110d-003,-1.4085d-004,
     .-8.8459d-006,+5.7256d-006,-1.5068d-006,+4.0095d-007,-2.4185d-008/

      data (b_geoid(i),i=1,55)/
     .+0.0000d+000,+0.0000d+000,-6.5993d-002,+0.0000d+000,+6.5364d-002,
     .-5.8320d+000,+0.0000d+000,+1.6961d+000,-1.3557d+000,+1.2694d+000,
     .+0.0000d+000,-2.9310d+000,+9.4805d-001,-7.6243d-002,+4.1076d-002,
     .+0.0000d+000,-5.1808d-001,-3.4583d-001,-4.3632d-002,+2.2101d-003,
     .-1.0663d-002,+0.0000d+000,+1.0927d-001,-2.9463d-001,+1.4371d-003,
     .-1.1452d-002,-2.8156d-003,-3.5330d-004,+0.0000d+000,+4.4049d-001,
     .+5.5653d-002,-2.0396d-002,-1.7312d-003,+3.5805d-005,+7.2682d-005,
     .+2.2535d-006,+0.0000d+000,+1.9502d-002,+2.7919d-002,-8.1812d-003,
     .+4.4540d-004,+8.8663d-005,+5.5596d-005,+2.4826d-006,+1.0279d-006,
     .+0.0000d+000,+6.0529d-002,-3.5824d-002,-5.1367d-003,+3.0119d-005,
     .-2.9911d-005,+1.9844d-005,-1.2349d-006,-7.6756d-009,+5.0100d-008/

      data (ap_mean(i),i=1,55)/ 
     .+1.0108d+003,+8.4886d+000,+1.4799d+000,-1.3897d+001,+3.7516d-003,
     .-1.4936d-001,+1.2232d+001,-7.6615d-001,-6.7699d-002,+8.1002d-003,
     .-1.5874d+001,+3.6614d-001,-6.7807d-002,-3.6309d-003,+5.9966d-004,
     .+4.8163d+000,-3.7363d-001,-7.2071d-002,+1.9998d-003,-6.2385d-004,
     .-3.7916d-004,+4.7609d+000,-3.9534d-001,+8.6667d-003,+1.1569d-002,
     .+1.1441d-003,-1.4193d-004,-8.5723d-005,+6.5008d-001,-5.0889d-001,
     .-1.5754d-002,-2.8305d-003,+5.7458d-004,+3.2577d-005,-9.6052d-006,
     .-2.7974d-006,+1.3530d+000,-2.7271d-001,-3.0276d-004,+3.6286d-003,
     .-2.0398d-004,+1.5846d-005,-7.7787d-006,+1.1210d-006,+9.9020d-008,
     .+5.5046d-001,-2.7312d-001,+3.2532d-003,-2.4277d-003,+1.1596d-004,
     .+2.6421d-007,-1.3263d-006,+2.7322d-007,+1.4058d-007,+4.9414d-009/ 

      data (bp_mean(i),i=1,55)/ 
     .+0.0000d+000,+0.0000d+000,-1.2878d+000,+0.0000d+000,+7.0444d-001,
     .+3.3222d-001,+0.0000d+000,-2.9636d-001,+7.2248d-003,+7.9655d-003,
     .+0.0000d+000,+1.0854d+000,+1.1145d-002,-3.6513d-002,+3.1527d-003,
     .+0.0000d+000,-4.8434d-001,+5.2023d-002,-1.3091d-002,+1.8515d-003,
     .+1.5422d-004,+0.0000d+000,+6.8298d-001,+2.5261d-003,-9.9703d-004,
     .-1.0829d-003,+1.7688d-004,-3.1418d-005,+0.0000d+000,-3.7018d-001,
     .+4.3234d-002,+7.2559d-003,+3.1516d-004,+2.0024d-005,-8.0581d-006,
     .-2.3653d-006,+0.0000d+000,+1.0298d-001,-1.5086d-002,+5.6186d-003,
     .+3.2613d-005,+4.0567d-005,-1.3925d-006,-3.6219d-007,-2.0176d-008,
     .+0.0000d+000,-1.8364d-001,+1.8508d-002,+7.5016d-004,-9.6139d-005,
     .-3.1995d-006,+1.3868d-007,-1.9486d-007,+3.0165d-010,-6.4376d-010/ 

      data (ap_amp(i),i=1,55)/ 
     .-1.0444d-001,+1.6618d-001,-6.3974d-002,+1.0922d+000,+5.7472d-001,
     .-3.0277d-001,-3.5087d+000,+7.1264d-003,-1.4030d-001,+3.7050d-002,
     .+4.0208d-001,-3.0431d-001,-1.3292d-001,+4.6746d-003,-1.5902d-004,
     .+2.8624d+000,-3.9315d-001,-6.4371d-002,+1.6444d-002,-2.3403d-003,
     .+4.2127d-005,+1.9945d+000,-6.0907d-001,-3.5386d-002,-1.0910d-003,
     .-1.2799d-004,+4.0970d-005,+2.2131d-005,-5.3292d-001,-2.9765d-001,
     .-3.2877d-002,+1.7691d-003,+5.9692d-005,+3.1725d-005,+2.0741d-005,
     .-3.7622d-007,+2.6372d+000,-3.1165d-001,+1.6439d-002,+2.1633d-004,
     .+1.7485d-004,+2.1587d-005,+6.1064d-006,-1.3755d-008,-7.8748d-008,
     .-5.9152d-001,-1.7676d-001,+8.1807d-003,+1.0445d-003,+2.3432d-004,
     .+9.3421d-006,+2.8104d-006,-1.5788d-007,-3.0648d-008,+2.6421d-010/ 

      data (bp_amp(i),i=1,55)/ 
     .+0.0000d+000,+0.0000d+000,+9.3340d-001,+0.0000d+000,+8.2346d-001,
     .+2.2082d-001,+0.0000d+000,+9.6177d-001,-1.5650d-002,+1.2708d-003,
     .+0.0000d+000,-3.9913d-001,+2.8020d-002,+2.8334d-002,+8.5980d-004,
     .+0.0000d+000,+3.0545d-001,-2.1691d-002,+6.4067d-004,-3.6528d-005,
     .-1.1166d-004,+0.0000d+000,-7.6974d-002,-1.8986d-002,+5.6896d-003,
     .-2.4159d-004,-2.3033d-004,-9.6783d-006,+0.0000d+000,-1.0218d-001,
     .-1.3916d-002,-4.1025d-003,-5.1340d-005,-7.0114d-005,-3.3152d-007,
     .+1.6901d-006,+0.0000d+000,-1.2422d-002,+2.5072d-003,+1.1205d-003,
     .-1.3034d-004,-2.3971d-005,-2.6622d-006,+5.7852d-007,+4.5847d-008,
     .+0.0000d+000,+4.4777d-002,-3.0421d-003,+2.6062d-005,-7.2421d-005,
     .+1.9119d-006,+3.9236d-007,+2.2390d-007,+2.9765d-009,-4.6452d-009/ 

      data (at_mean(i),i=1,55)/ 
     .+1.6257e+001,+2.1224e+000,+9.2569e-001,-2.5974e+001,+1.4510e+000,
     .+9.2468e-002,-5.3192e-001,+2.1094e-001,-6.9210e-002,-3.4060e-002,
     .-4.6569e+000,+2.6385e-001,-3.6093e-002,+1.0198e-002,-1.8783e-003,
     .+7.4983e-001,+1.1741e-001,+3.9940e-002,+5.1348e-003,+5.9111e-003,
     .+8.6133e-006,+6.3057e-001,+1.5203e-001,+3.9702e-002,+4.6334e-003,
     .+2.4406e-004,+1.5189e-004,+1.9581e-007,+5.4414e-001,+3.5722e-001,
     .+5.2763e-002,+4.1147e-003,-2.7239e-004,-5.9957e-005,+1.6394e-006,
     .-7.3045e-007,-2.9394e+000,+5.5579e-002,+1.8852e-002,+3.4272e-003,
     .-2.3193e-005,-2.9349e-005,+3.6397e-007,+2.0490e-006,-6.4719e-008,
     .-5.2225e-001,+2.0799e-001,+1.3477e-003,+3.1613e-004,-2.2285e-004,
     .-1.8137e-005,-1.5177e-007,+6.1343e-007,+7.8566e-008,+1.0749e-009/ 

      data (bt_mean(i),i=1,55)/ 
     .+0.0000e+000,+0.0000e+000,+1.0210e+000,+0.0000e+000,+6.0194e-001,
     .+1.2292e-001,+0.0000e+000,-4.2184e-001,+1.8230e-001,+4.2329e-002,
     .+0.0000e+000,+9.3312e-002,+9.5346e-002,-1.9724e-003,+5.8776e-003,
     .+0.0000e+000,-2.0940e-001,+3.4199e-002,-5.7672e-003,-2.1590e-003,
     .+5.6815e-004,+0.0000e+000,+2.2858e-001,+1.2283e-002,-9.3679e-003,
     .-1.4233e-003,-1.5962e-004,+4.0160e-005,+0.0000e+000,+3.6353e-002,
     .-9.4263e-004,-3.6762e-003,+5.8608e-005,-2.6391e-005,+3.2095e-006,
     .-1.1605e-006,+0.0000e+000,+1.6306e-001,+1.3293e-002,-1.1395e-003,
     .+5.1097e-005,+3.3977e-005,+7.6449e-006,-1.7602e-007,-7.6558e-008,
     .+0.0000e+000,-4.5415e-002,-1.8027e-002,+3.6561e-004,-1.1274e-004,
     .+1.3047e-005,+2.0001e-006,-1.5152e-007,-2.7807e-008,+7.7491e-009/ 

      data (at_amp(i),i=1,55)/ 
     .-1.8654e+000,-9.0041e+000,-1.2974e-001,-3.6053e+000,+2.0284e-002,
     .+2.1872e-001,-1.3015e+000,+4.0355e-001,+2.2216e-001,-4.0605e-003,
     .+1.9623e+000,+4.2887e-001,+2.1437e-001,-1.0061e-002,-1.1368e-003,
     .-6.9235e-002,+5.6758e-001,+1.1917e-001,-7.0765e-003,+3.0017e-004,
     .+3.0601e-004,+1.6559e+000,+2.0722e-001,+6.0013e-002,+1.7023e-004,
     .-9.2424e-004,+1.1269e-005,-6.9911e-006,-2.0886e+000,-6.7879e-002,
     .-8.5922e-004,-1.6087e-003,-4.5549e-005,+3.3178e-005,-6.1715e-006,
     .-1.4446e-006,-3.7210e-001,+1.5775e-001,-1.7827e-003,-4.4396e-004,
     .+2.2844e-004,-1.1215e-005,-2.1120e-006,-9.6421e-007,-1.4170e-008,
     .+7.8720e-001,-4.4238e-002,-1.5120e-003,-9.4119e-004,+4.0645e-006,
     .-4.9253e-006,-1.8656e-006,-4.0736e-007,-4.9594e-008,+1.6134e-009/ 

      data (bt_amp(i),i=1,55)/ 
     .+0.0000e+000,+0.0000e+000,-8.9895e-001,+0.0000e+000,-1.0790e+000,
     .-1.2699e-001,+0.0000e+000,-5.9033e-001,+3.4865e-002,-3.2614e-002,
     .+0.0000e+000,-2.4310e-002,+1.5607e-002,-2.9833e-002,-5.9048e-003,
     .+0.0000e+000,+2.8383e-001,+4.0509e-002,-1.8834e-002,-1.2654e-003,
     .-1.3794e-004,+0.0000e+000,+1.3306e-001,+3.4960e-002,-3.6799e-003,
     .-3.5626e-004,+1.4814e-004,+3.7932e-006,+0.0000e+000,+2.0801e-001,
     .+6.5640e-003,-3.4893e-003,-2.7395e-004,+7.4296e-005,-7.9927e-006,
     .-1.0277e-006,+0.0000e+000,+3.6515e-002,-7.4319e-003,-6.2873e-004,
     .-8.2461e-005,+3.1095e-005,-5.3860e-007,-1.2055e-007,-1.1517e-007,
     .+0.0000e+000,+3.1404e-002,+1.5580e-002,-1.1428e-003,+3.3529e-005,
     .+1.0387e-005,-1.9378e-006,-2.7327e-007,+7.5833e-009,-9.2323e-009/ 

C     parameter t
      t = dsin(dlat)

C     degree n and order m
      n = 9
      m = 9

C determine n!  (faktorielle)  moved by 1
      dfac(1) = 1
      do i = 1,(2*n + 1)
        dfac(i+1) = dfac(i)*i
      end do

C     determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
      do i = 0,n
        do j = 0,min(i,m)
          ir = int((i - j)/2)
          sum = 0
          do k = 0,ir
            sum = sum + (-1)**k*dfac(2*i - 2*k + 1)/dfac(k + 1)/
     .       dfac(i - k + 1)/dfac(i - j - 2*k + 1)*t**(i - j - 2*k)
          end do
C         Legendre functions moved by 1
          P(i + 1,j + 1) = 1.d0/2**i*dsqrt((1 - t**2)**(j))*sum
        end do
      end do

C     spherical harmonics
      i = 0
      do n = 0,9
        do m = 0,n
          i = i + 1
          aP(i) = P(n+1,m+1)*dcos(m*dlon)
          bP(i) = P(n+1,m+1)*dsin(m*dlon)
        end do
      end do

C     Geoidal height
      undu = 0.d0
      do i = 1,55
        undu = undu + (a_geoid(i)*aP(i) + b_geoid(i)*bP(i))
      end do

C     orthometric height
      hort = dhgt - undu

C     Surface pressure
      apm = 0.d0
      apa = 0.d0
      do i = 1,55
        apm = apm + (ap_mean(i)*aP(i) + bp_mean(i)*bP(i))
        apa = apa + (ap_amp(i) *aP(i) + bp_amp(i) *bP(i))
      end do
      pres0  = apm + apa*dcos(doy/365.25d0*2.d0*pi)

C     height correction for pressure
      pres = pres0*(1.d0-0.0000226d0*hort)**5.225d0

C     Surface temperature
      atm = 0.d0
      ata = 0.d0
      do i = 1,55
        atm = atm + (at_mean(i)*aP(i) + bt_mean(i)*bP(i))
        ata = ata + (at_amp(i) *aP(i) + bt_amp(i) *bP(i))
      end do
      temp0 =  atm + ata*dcos(doy/365.25d0*2*pi)

C     height correction for temperature
      temp = temp0 - 0.0065d0*hort
          
      return
      end
