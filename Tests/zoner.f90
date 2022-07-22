! This version of zoner is rewritten in f90. 
! by Chance adapted from C. Constable
! most up to date version of zoner as of 7/15/2022
! Some modifications and additions are made.
! The additions: exporting a tessel file written in terms of latitude, longitude, and core radius.
! This export is titled llTessel
program zoner
    implicit none

    real :: long, bColat1 , bColat2 , bLong1 , bLong2 , colat , cp , &
    & ct , toRad , phi , sp , st , theta, x(3), p(3, 6), &
    & coord(2), q(2, 2), r, dPhi, delta, deltaTwiddle, dPhiTwiddle
    integer :: i, j, inBox , iflag , kount , mu, nu, n

!   Generates approximately equally spaced unit vectors on a sphere with
!   user-specified separation DELTA.  Start with north pole, then increment colat by delta and go around
!   latitude bands - , end with SP.
! specify north and south poles and 4 equidistant equatorial spots  as first 6 points
    data ((p(i,j),i=1,3),j=1,6)/0 , 0 , 1 , 0 , 0 , -1 , -1 , -1 , 0 , &
    & 1 , -1 , 0 , 1 , 1 , 0 , -1 , 1 , 0/
    data ((q(i, j), i=1, 2), j=1, 2)/90.0, 0.0, -90.0, 0.0/

! Set some values
    data toRad/0.01745329/
    data r/3486.0E3/
    data n/6/

! Inputs to exist as common block boxy
    common /boxy  /bColat1 , bColat2 , bLong1 , bLong2

! Enter bounding box, colat1,2 ,long1,2
    print *, 'Enter colatitudes 1&2 (North to South), ', & 
    & 'longitudes 1&2 (West (-) to East (+)) for bounding box'
    print * , '(e.g., for "front" of earth, 0 180 -90 90)'
    read * , bColat1 , bColat2 , bLong1 , bLong2
    if (bLong2.gt.180) then
        print *, 'Upper longitudinal bound must be < 180 degrees'
        stop
    elseif (blong1.lt.-180) then
        print *, 'Lower longitudinal bound must be > -180 degrees'
        stop
    endif
    print * , bColat1 , bColat2 , bLong1 , bLong2

! Request interval Delta in degrees
    print *, 'Enter Delta in degrees'
    print *, 'note: number of lat/long divisions will be rounded down'
    read *, delta

! Open files for writing
    open (unit=3,file='tessel')
    open (unit=30, file='llTessel')	

! Need to manually write the poles. Starting with north
! iflag must be fix either way, otherwise it returns a crazy number
    if (BCOlat1.eq.0.0) then
        iflag = 1
    else
        iflag = -1
    endif
    kount = 0
do j=1, 1
    write (3, *) (p(i, j), i=1, 3), iflag
    write (30, *) (q(i, j), i=1, 2), r, iflag
enddo
kount = 1

! Loop over CMB via latitude and longitude using specified delta  
! Define a colatitude, latitude, and cos/sin, which require radians
! Spherical convention: 
! Keep in mind 
mu = nint(180.0/delta)
deltaTwiddle = 180.0/mu
if (delta.ne.deltaTwiddle) then
    print *, "To ensure equal spacing between tesselation vertices, "
    print *, "Delta has been assigned value", deltaTwiddle
else
    print *, 'Assigning vertices for a tesselation with '
    print *, 'angular distance', deltaTwiddle
endif

! print *, 'mu is ', mu
! nu must be looped from 1, otherwise 0 and 360 will be same point
! Every other latitudinal line should be shifted by delta/2 to make triangles
! Shift every even j-indexed long by delta/2
do i=1, (mu-1)
    colat = i*deltaTwiddle
    theta = toRad * colat
    ct = cos(theta)
    st = sin(theta)
    ! Want to ensure equal area divisions
    dPhi = deltaTwiddle/st
    ! nu is an integer
    nu = nint(360.0/dPhi)
    dPhiTwiddle = 360.0/nu
!    print *, 'nu times dPhiTwiddle is ', nu*dPhiTwiddle
!    print *, 'colat is ', colat
    do j=0, (nu-1)
        long = j*dPhiTwiddle
    ! shifting odd i is useful for illustration, but even i minimizes the downsized triangles
    !    if (mod(i, 2).eq.0) long = long - dPhiTwiddle/2.0
        phi = toRad * long
        cp = cos(phi)
        sp = sin(phi)
        x(1) = st*cp
        x(2) = st*sp
        x(3) = ct
!        print *, iflag
        if (long.gt.180) long = long - 360
        coord(1) = 90 - colat
        coord(2) = long
        iflag = inBox(colat, long)
        write (3, *) x, iflag
        write (30, *) coord, r, iflag
        kount = kount + 1
    enddo
enddo

! South pole too
if (BCOlat2.eq.180.0) then
    iflag = 1
else
    iflag = -1
endif
do j=2, 2
write (3, *) (p(i, j), i=1, 3), iflag
write (30, *) (q(i, j), i=1, 2), r, iflag
enddo
kount = kount + 1

! wrap it up
PRINT * , kount , 'vectors written to tessel and lltessel'

! dont need tessel file anymore
! CALL execute_command_line('rm tessel' )

end program zoner

!
! _________________________________________________

integer function inBox(colat, long)
implicit none

! long was the longitude input. 
! It is different from long since there are no formatting constraints
real, intent(in) :: colat, long
real :: bColat1 , bColat2 , bLong1 , bLong2
integer :: inLat, inLong

common /boxy /bColat1 , bColat2 , bLong1 , bLong2

if (colat.ge.bColat1 .and. colat.le.bColat2) then
    inLat = 1
else
    inLat = -1
endif
 

if (long.ge.bLong1 .and. long.le.bLong2) then
    inLong = 1
else
    inLong = -1
endif

if (inLat.eq.1 .and. inLong.eq.1) then
    inBox=1
else
    inBox=-1
endif
!print *, 'for ', '(', colat, ', ', long, '): '
!print *, 'inLat is ', inLat
!print *, 'inLong is ', inLong
!print *, 'inBox is ', inBox
return
end function inBox

!
! _________________________________________________

! Truncate number to n decimals
real function truncate(num, n)
implicit none

integer :: n, numBig
real :: num

! numBig is an integer, this is the truncation step
numBig = num * (10.0 ** n)
num = numBig/(10.0 ** n)

return

end function truncate