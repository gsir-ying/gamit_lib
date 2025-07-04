      subroutine geodms_xyz (dlat,mlat,seclat,latflag,dlon,mlon,seclon
     .                      ,lonflag,radius,xstn,ystn,zstn)
c
c     Subroutine to compute space rectangular coordinates from
c     either geocentric coordinates or geodetic coordinates. The
c     latter are assumed to be WGS84.
c     Coordinate system is assumed to be right handed, with
c     positive longitude to the East
c

      integer*4 dlat,mlat,dlon,mlon
      real*8 seclat,seclon,radius,xstn,ystn,zstn
      real*8 conv, a, flattening, flat, datan
      real*8 latitude, longitude, n, nplush, esquare, pi

      character latflag*1, lonflag*1,south,west

c     case conversion function
      character lowerc

      data a          /6378137./,
     .     flattening /298.25722201/,
     .     south      /'S'/,
     .     west       /'W'/

      pi= 4.0d0*datan(1.0d0)

      conv =180.d0/pi
      latitude = (dble(dlat)+dble(mlat)/60.d0+seclat/3600.d0)/conv
      longitude= (dble(dlon)+dble(mlon)/60.d0+seclon/3600.d0)/conv

      if(lowerc(latflag).eq.lowerc(south)) latitude=-latitude
      if(lowerc(lonflag).eq.lowerc(west))  longitude=-longitude
      flat = 1.d0/flattening
      esquare = 2.d0*flat-flat*flat
      if(radius.gt.2000000.d0) esquare = 0.d0
      n = a/dsqrt(1.d0-esquare*dsin(latitude)**2)
      if ( radius .gt. 2000000.d0 ) n = 0.d0
      nplush = n+radius
      xstn = nplush * dcos(latitude) * dcos(longitude)
      ystn = nplush * dcos(latitude) * dsin(longitude)
      zstn = (n * (1.d0-esquare) + radius) * dsin(latitude)
      return
      end
