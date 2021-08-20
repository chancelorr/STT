	program zoner
c  Starting afresh April 2020,  C. Constable
c
c  Generates approximately equally spaced unit vectors on a sphere with
c  user-specified separation DELTA.  Start with north pole, then increment colat by delta and go around 
c  latitude bands - , end with SP.
c  
c  
c  Divide vectors into two regions specified by  a colat, long bounding box
c  e.g. 0,90, -180, 180 gives N. hemisphere, flag 1, then SH will have flag -1, boundary will be 0
c  0, 180, -90, 90 gives "front" of Earth
c
      dimension x(3),p(3,6)
      data deg/0.01745329/
c specify north and south poles and 4 equidistant equatorial spots  as first 6 points
	data ((p(i,j),i=1,3),j=1,6)/0,0,1, 0,0,-1, -1,-1,0, 1,-1,0, 1,1,0, 
     $ -1,1,0 /
	common /boxy/ bcolat1,bcolat2,blong1,blong2
c
c Enter bounding box, colat1,2 ,long1,2
	print *, 'Enter colat1,2 and long1,2 for bounding box'
	read *, bcolat1, bcolat2, blong1, blong2 
	print *, bcolat1, bcolat2, blong1, blong2 
c is  north pole in box?
	if(bcolat1.eq.0.0)iflag =1
      kount=0
      open (unit=3, file='tessel')
	do 9 j=1,1
	write(3,*)(p(i,j),i=1,3),iflag
9	continue
	kount=2
10      print *,' Enter DELTA in degrees'
      read *,delta
	print *,' assigning vertices for a tesselation 
     $  with angular distance', delta 
c      nu=90.0/delta -0.5
c      delta=180.0/(2.0*nu + 1.0)
c      print '(a,f9.3)',' DELTA modified to be ',delta
	nu = 179/delta
c	print *, 'nu= ',nu
      do 1500 i=1,nu
	colat=I*delta
        theta=deg*colat
        ct=cos(theta)
        st=sin(theta)
        step=360.0/int(st*360.0/delta)
	jmax = 360/step 
c	print *, 'jmax= ',jmax
        do 1400 j= 1, jmax
	  along= (j-1)*step
          phi=deg*(along + 0.4771*colat)
c uncomment  next line to align all trainsgle on zero longitude
c	  phi = deg*along
          cp=cos(phi)
          sp=sin(phi)
          x(1)=cp*st
          x(2)=sp*st
          x(3)=ct
	  iflag = ibox(colat,along)
          write(3,*) x, iflag
          kount=1 + kount
 1400   continue
 1500 continue
c is pole in box?
	if(bcolat2.eq.180.0)iflag =1
	do 19 j=2,2
	write(3,*)(p(i,j),i=1,3),iflag
19	continue
      print *,kount,'vectors written to tessel'
      stop
      end
c_________________________________________________
c
	integer function ibox(colat,along)
c Check if colat,along is inside, outside or on boundary of box specified in boxy 
c Box in degrees  (0,180,0,360)
	common /boxy/ bcolat1,bcolat2,blong1,blong2
c
	if(colat.ge.bcolat1.and.colat.le.bcolat2) then
		inlat=1
	else 
		inlat =-1
	endif
	if(along.ge.blong1.and.along.le.blong2) then
		inlong=1
	else 
		inlong =-1
	endif
	if (inlat.eq.1.and.inlong.eq.1)then
	  ibox=1
	else
	  ibox=-1
	endif
	return
	end
