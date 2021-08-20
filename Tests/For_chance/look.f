	program look
c
c  Looks at triangulations in various projections with aid of plotxy
c  Triangles are read as groups of three vectors, output from synops
c  will do.
c  Centroids of triangles are written to a separate file.
c
      character*64 name
      dimension x(3),y(3),z(3),u(2,500),w(3)
      common /sides/ ns,s(3,500)
c
      print *,' Enter file name of vectors'
      read '(a64)', name
      print *,' Enter 0 for orthographic, 1 for equal-area projection,
     $ 2 for radian'
      read *,jproj
c
      open(unit=9, file=name)
      do 1500 i=1, 10000
        read (9,*, end=1600) x,y,z
        call norm(x)
        call norm(y)
        call norm(z)
        ns=0
        call inter(x, y)
        call inter(y, z)
        call inter(z, x)
        call inter(x, x)
        if (jproj .eq. 0) then
c  Northern hemisphere with Greenwhich on -y axis
          if (x(3).ge.0 .and. y(3).ge.0 .and. z(3).ge.0)
     $    write(1,'(2f8.4)')  (s(2,k),-s(1,k), k=1,ns)
c  Southern hemisphere with Greenwhich on +y axis 
          if (x(3).le.0 .and. y(3).le.0 .and. z(3).le.0)
     $    write(2,'(2f8.4)') (s(2,k),s(1,k), k=1,ns)
        else if (jproj.eq.1) then
          call lamb( 1, ns, s, u)
          write(1,'(2f8.4)') (u(1,k),u(2,k), k=1,ns)
          call averag(x, y, z, w)
          call lamb( 1, 2, w, u(1,1))
          if (u(1,1)**2+u(2,1)**2 .le. 1.1)
     $     write(3,'(2f8.4,i4)') u(1,1),u(2,1),i
c
          call lamb(-1, ns, s, u)
          write(2,'(2f8.4)') (u(1,k),u(2,k), k=1,ns)
          call lamb(-1, 2, w, u(1,1))
          if (u(1,1)**2+u(2,1)**2 .le. 1.1)
     $     write(4,'(2f8.4,i4)') u(1,1),u(2,1),i
	else
	  call radian(ns,s,u)  
          write(1,'(2f8.4)') (u(1,k),u(2,k), k=1,ns)
c          call averag(x, y, z, w)
c          call radian(2, w, u(1,1))
          write(2,'(2f8.4,i4)') (u(1,k),u(2,k),k,k=1,ns)
        endif
 1500 continue
 1600 stop
      end
c____________________________________________________
      subroutine lamb(ihemi, ns, s, u)
c$$$$ calls no other routines
c  Performs mapping into Lambert equal; area projection with North or
c  South pole as map origin.  
c  The circle u(1)**2+u(2)**2=1 is the image of the equator.
      dimension s(3,ns),u(2,ns)
c
c  Northern hemisphere
      if (ihemi .gt. 0) then
      do 1100 i=1, ns-1
        R=sqrt(1.0/((1.0 + s(3,i))+1.0e-6))
        u(1,i)= s(2,i)*R
        u(2,i)=-s(1,i)*R
 1100 continue
      else
c  Southern hemisphere
      do 1200 i=1, ns-1
        R=sqrt(1.0/((1.0 - s(3,i))+1.0e-6))
        u(1,i)= s(2,i)*R
        u(2,i)= s(1,i)*R
 1200 continue
      endif
c  Supress large jumps that arise from sides near the antipode
      do 1500 i=2, ns-1
        if (abs(u(1,i)-u(1,i-1))+abs(u(2,i)-u(2,i-1)) .gt. 1.0 .and.
     $  u(1,i-1) .lt. 2.0) u(1,i)=9.0
 1500 continue
      u(1,ns)=9.0
      u(2,ns)=9.0
      return
      end
c_____________________________________________________________________
      subroutine copy(x, y)
c$$$$ calls no other routines
c  Copies 3-vector x into y
      dimension x(3),y(3)
      y(1)=x(1)
      y(2)=x(2)
      y(3)=x(3)
      return
      end
c____________________________________________________
      subroutine norm(x)
c$$$$ calls no other routines
c  Normalizes 3-vector x to unit length
      dimension x(3)
      xx=1.0/sqrt(x(1)**2 + x(2)**2 + x(3)**2)
      x(1)=x(1)*xx
      x(2)=x(2)*xx
      x(3)=x(3)*xx
      return
      end
c____________________________________________________
      subroutine inter(x, y)
c$$$$ calls copy, norm
c  Generates a sequence of unit vectors  s(*,1),s(*,2) ... s(*,ns) in
c  common /sides/ lying in the c  plane of the vectors x, y.
c  s(*,1)=x,  s(*,2), s(*,3) are convex interpolants.  Note that
c  unless x=y, s(*,ns) is not in fact y.
c
      dimension x(3),y(3),z(3),q(3)
      common /sides/ ns,s(3,500)
      data q/9.0,9.0,0.0/
c
      dot=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
      if (dot .gt. 0.99999) then
         call copy(x, s(1,ns+1))
         call copy(q, s(1,ns+2))
         ns=ns+2
      else
         m=(1.0 - dot)*30.0 + 1.0
         d=1.0/m
         do 1200 i=1, m
          a=(i-1)*d
          z(1)=x(1)*(1.0 - a) + y(1)*a
          z(2)=x(2)*(1.0 - a) + y(2)*a
          z(3)=x(3)*(1.0 - a) + y(3)*a
          call norm(z)
          call copy(z, s(1,ns+i))
 1200   continue
        ns=ns + m
      endif
      return
      end
c____________________________________________________
      subroutine averag(x, y, z, w)
c$$$$ calls norm
      dimension x(*),y(*),z(*),w(*)
c
      w(1)=x(1) + y(1) + z(1)
      w(2)=x(2) + y(2) + z(2)
      w(3)=x(3) + y(3) + z(3)
      call norm(w)
      return
      end
c______________________________________________________
	subroutine radian(ns,s,u)
c $$$calls no other routines 
c  Maps onto radian projection u(1,i)=longitude, u(2,i)=lat
c
c  Unit vectors assumed
	dimension s(3,ns),u(2,ns)
	data pi/3.1415927/
	do 1100 i=1,ns-1
	  u(1,i)= atan2(s(2,i),s(1,i))
	  u(2,i)= asin(s(3,i))
1100	continue
c  inhibit wrapping in longitude
	do 1200 i=2,ns-1
	if (abs(u(2,i)).le.1.4.and.abs(u(2,i-1)).le.1.4)then
  	   if((u(1,i)-u(1,i-1)).gt.2.5)u(1,i)=u(1,i)-2*pi
	   if((u(1,i)-u(1,i-1)).lt.-2.5)u(1,i)=u(1,i)+2*pi
	endif
1200	continue
	u(2,ns)=9.0
	u(1,ns)=9.0
	return
	end
