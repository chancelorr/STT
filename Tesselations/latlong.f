c  reads vectors defining tesseltaion on core, outputs lat, long ,pt
c  no. Just acts as filter from standard in to standard out
	dimension x(3)
	data drad/57.29578/
	i=0
1	read(5,*,end=10,err=10)x(1),x(2),x(3)
	i=i+1
	rlat=90. - drad*acos(x(3))
	rlong =drad*atan2(x(2),x(1))
	write(6,*)rlong,rlat,i
	goto 1
10	stop
	end
