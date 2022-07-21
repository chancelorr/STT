c_______________________________________________________________________
      subroutine synops(new, n, p, ntria, iv)
	implicit double precision (a-h,o-z)
c$$$$ calls dotty
c  prints number of triangles in tesselation and outputs list of
c  triangles to fort.1 if new>=2
c
      dimension p(4,*),iv(3,ntria)
      	common /io/ inp,iout,iprint
 	common /bound/dmin,dmax,amin,amax
	sarea=0.
	if (new.eq.0)then
	amin=12.
	amax=0.
     	do 60 i=1,ntria
	i1=iv(1,i)
	i2=iv(2,i)
	i3=iv(3,i)
	x=area(p(1,i1),p(1,i2),p(1,i3))
	sarea=sarea + x
	if(x.lt.amin)amin=x
	if(x.gt.amax)amax=x
	p12= dotty(p(1,i1),p(1,i2))
	dmax=max(p12,dmax)
	dmin=min(p12,dmin)
	p13=dotty(p(1,i1),p(1,i3))
	dmax=max(p13,dmax)
	dmin=min(dmin,p13)
	p23=dotty(p(1,i2),p(1,i3))
	dmax=max(p23,dmax)
	dmin=min(dmin,p23)
 60     continue
	else if(new.eq.2) then
	open(unit=4,file='fort.4')
	open(unit=5,file='fort.5')
     	do 70 i=1,ntria
	i1=iv(1,i)
	i2=iv(2,i)
	i3=iv(3,i)
	x=area(p(1,i1),p(1,i2),p(1,i3))
	write(4,*)x
	sarea=sarea + x
	if(x.lt.amin)amin=x
	if(x.gt.amax)amax=x
	p12= dotty(p(1,i1),p(1,i2))
	dmax=max(p12,dmax)
	dmin=min(p12,dmin)
	p13=dotty(p(1,i1),p(1,i3))
	dmax=max(p13,dmax)
	dmin=min(dmin,p13)
	p23=dotty(p(1,i2),p(1,i3))
	dmax=max(p23,dmax)
	dmin=min(dmin,p23)
	write(5,*)p12,p13,p23
70     continue
	endif
c
c  Describe triangle model
      write(iout,'(//a,i6)')' Number of triangles in body ',ntria
c
      write(iout, '(/a,2f8.2,a)')' Triangle sides range between ',
     $ 57.296*acos(dmax),57.296*acos(dmin),' degrees'
      write(iout, '(/a,f8.4,a,f8.4)')' Triangle areas range between ',
     $ amax,' and ',amin
c
      if (new .le. 1) return
c
c  For spherical triangles
c  Output to a file every  triangle in the list to fort.2
      open (unit=11, file='fort.2')
      do 1500 it=1, ntria
        i1=iv(1,it)
        i2=iv(2,it)
        i3=iv(3,it)
        write(11,'(3p,9f8.0,a)') p(1,i1),p(2,i1),p(3,i1),
     $  p(1,i2),p(2,i2),p(3,i2), p(1,i3),p(2,i3),p(3,i3)
 1500 continue
c
      close (unit=11)
      close (unit=4)
      close (unit=5)
      write(iout,'(/a)') ' Spherical triangles written to file fort.2'
      write(iout,'(/a)') ' Triangle areas written to file fort.4'
      write(iout,'(/a)') ' Triangle angles written to file fort.5'
      return
      end
