PROGRAM ZONER
    IMPLICIT NONE
!*--ZONER4
!*** Start of declarations inserted by SPAG
    REAL  :: along , BCOlat1 , BCOlat2 , BLOng1 , BLOng2 , colat , cp ,   &
		   & ct , deg , delta , phi , sp , st , step , theta, x(3), p(3, 6), &
		   & coord(2), q(2, 2), r, long, normR
	!REAL :: thetaComp(2), phiComp(2), normR
    INTEGER :: i , IBOX , iflag , j , jmax , kount , nu, lkount
!*** End of declarations inserted by SPAG
!  Starting afresh Aug 2021,  C. Roberts with C. Constable's version
!
!
!	
!  Generates approximately equally spaced unit vectors on a sphere with
!  user-specified separation DELTA.  Start with north pole, then increment colat by delta and go around
!  latitude bands - , end with SP.
!
!
!  Divide vectors into two regions specified by  a colat, long bounding box
!  e.g. 0,90, -180, 180 gives N. hemisphere, flag 1, then SH will have flag -1, boundary will be 0
!  0, 180, -90, 90 gives "front" of Earth
!
!  Write only the vectors within the bounding box to the
!
!
        DATA deg/0.01745329/
		DATA r/3486.0E3/
! specify north and south poles and 4 equidistant equatorial spots  as first 6 points
        DATA ((p(i,j),i=1,3),j=1,6)/0 , 0 , 1 , 0 , 0 , -1 , -1 , -1 , 0 ,&
			& 1 , -1 , 0 , 1 , 1 , 0 , -1 , 1 , 0/
		DATA ((q(i, j), i=1, 2), j=1, 2)/90.0, 0.0, -90.0, 0.0/
        COMMON /BOXY  / BCOlat1 , BCOlat2 , BLOng1 , BLOng2
! Enter bounding box, colat1,2 ,long1,2
        PRINT * , 'Enter colat1,2 and long1,2 for bounding box'
        READ * , BCOlat1 , BCOlat2 , BLOng1 , BLOng2
        PRINT * , BCOlat1 , BCOlat2 , BLOng1 , BLOng2
! Open files for writing
		OPEN (UNIT=3,FILE='tessel')
		OPEN (UNIT=30, FILE='llTessel')
		! OPEN (UNIT=300, FILE='thetaCompare')
		! OPEN (UNIT=3000, FILE='phiCompare')
		! OPEN (UNIT=30000, FILE='unity')
! is  north pole in box? If not, value still needs to be fixed or it will take on random default value
        IF ( BCOlat1.EQ.0.0 ) THEN
			iflag = 1
			lkount = 1
		ELSE
			iflag = -1
!			lkount = 0
		ENDIF
        kount = 0
        DO j = 1 , 1
            WRITE (3,*) (p(i,j),i=1,3) , iflag
			WRITE (30, *) (q(i, 1), i=1, 2), r , iflag
			! WRITE (300, *) 0, 0
			! WRITE (30000, *) 1.00
        ENDDO
        kount = 1
        PRINT * , ' Enter DELTA in degrees'
        READ * , delta
        PRINT * ,' assigning vertices for a tesselation with ', &
		& 'angular distance', delta
!      nu=90.0/delta -0.5
!      delta=180.0/(2.0*nu + 1.0)
!      print '(a,f9.3)',' DELTA modified to be ',delta
        nu = 179/delta
!	These are all unit vectors, so we are effectively taking r=1
!	print *, 'nu= ',nu
        DO i = 1 , nu
		    colat = i*delta
			theta = deg*colat
			ct = COS(theta)
			st = SIN(theta)
			step = 360.0/INT(st*360.0/delta)
			jmax = 360/step
!	print *, 'jmax= ',jmax
			DO j = 1 , jmax
			  along = (j-1)*step
			  long = (along+0.4771*colat)
			  phi = deg*long
!	uncomment  next line to align all trainsgle on zero longitude
!	phi = deg*along
			  cp = COS(phi)
			  sp = SIN(phi)
			  x(1) = st*cp
			  x(2) = st*sp
			  x(3) = ct
			  iflag = IBOX(colat,along)
			  normR = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
			  IF (long.gt.180) long = long - 360
			  coord(1) = 90 - acos(x(3)/normR)/deg
			  coord(2) = atan(x(2)/x(1))/deg
			!   thetaComp(1) = colat
			!   thetaComp(2) = acos(x(3)/normR)/deg
			!   phiComp(1) = long
			!   phiComp(2) = atan(x(2)/x(1))/deg
			  WRITE (30, *) coord, r, iflag
			!   WRITE (300, *) thetaComp
			!   WRITE (3000, *) phiComp
			!   WRITE (30000, *) normR
			  WRITE (3, *) x , iflag
			  kount = 1 + kount
		    ENDDO
		ENDDO
! is south pole in box?
		IF ( BCOlat2.EQ.180.0 ) THEN
			iflag = 1
!			lkount = lkount + 1
		ELSE
			iflag=-1
		ENDIF
		DO j = 2 , 2
		   WRITE (3,*) (p(i,j),i=1,3) , iflag
		   WRITE (30, *) (q(i, 2), i=1, 2), r, iflag
		!    WRITE (300, *) 180, 180
		!    WRITE (30000, *) 1.00
		ENDDO
		kount = kount + 1
		PRINT * , kount , 'vectors written to tessel'
		PRINT * , kount , 'vectors written to llTessel'
END PROGRAM 
!*==IBOX.spg  processed by SPAG 6.72Dc at 00:33 on 24 Aug 2021
!_________________________________________________
!
		INTEGER FUNCTION IBOX(Colat,Along)
		IMPLICIT NONE
!*--IBOX85
!*** Start of declarations inserted by SPAG
		REAL Along , BCOlat1 , BCOlat2 , BLOng1 , BLOng2 , Colat
		INTEGER inlat , inlong
!*** End of declarations inserted by SPAG
! Check if colat,along is inside, outside or on boundary of box specified in boxy
! Box in degrees  (0,180,0,360)
		COMMON /BOXY  / BCOlat1 , BCOlat2 , BLOng1 , BLOng2
!
		IF ( Colat.GE.BCOlat1 .AND. Colat.LE.BCOlat2 ) THEN
		   inlat = 1
		ELSE
		   inlat = -1
		ENDIF
		IF ( Along.GE.BLOng1 .AND. Along.LE.BLOng2 ) THEN
		   inlong = 1
		ELSE
		   inlong = -1
		ENDIF
		IF ( inlat.EQ.1 .AND. inlong.EQ.1 ) THEN
		   IBOX = 1
		ELSE
		   IBOX = -1
		ENDIF
		END
!
!
! __________________________________________________________________________