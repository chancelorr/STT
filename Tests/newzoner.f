c  Generates approximately equally spaced unit vectors on a sphere with
c  user-specified separation DELTA.  Simply goes round in constant
c  latitude bands - start in NH with NP, end with SP.
c  
c  modified June 2009 C. Constable
c  Divide vectors into two regions specified by  a colat, long bounding box
c  e.g. 0,90, -180, 180 gives N. hemisphere, flag 1, then SH will have flag -1, boundary will be 0
c  0, 180, -90, 90 gives "front" of Earth
c
      dimension x(3),p(3,6)
      data deg/0.01745329/
	data ((p(i,j),i=1,3),j=1,2)/0,0,1,0,0,-1/
	common /boxy/ bcolat1,bcolat2,blong1,blong2
c
c Enter bounding box, colat1,2 ,long1,2
	print *, 'Enter colat1,1 and long1,2 for bounding box'
	read *, bcolat1, bcolat2, blong1, blong2 
c is pole in box?
	if(bcolat1.eq.0.0)iflag =1
      kount=0
      open (unit=3, file='tessel')
	do 9 j=1,1
	write(3,*)(p(i,j),i=1,3),iflag
9				continue
	kount=2
10      print *,' Enter DELTA in degrees'
      read *,delta
      nu=90.0/delta -0.5
      delta=180.0/(2.0*nu + 1.0)
      print '(a,f9.3)',' DELTA modified to be ',delta
      do 1500 colat=delta, 179.0, delta
        theta=deg*colat
        ct=cos(theta)
        st=sin(theta)
        step=360.0/int(st*360.0/delta)
        do 1400 along=0.0, 359.0, step
          phi=deg*(along + 0.4771*colat)
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
