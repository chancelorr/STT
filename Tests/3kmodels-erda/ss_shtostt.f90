program shtostt 
      implicit none 
!     Modified version of shtostt which produces only a SNAPSHOT 
!     Modified and translated into .f90 by Chance, Sept 2021
!                                                                       
!                                                                       
!     Modified version of fieldpred.f to read stt tesselation plus ipp z
! 	output is stt, br,ipp                                                
!     from models ARCH3k.1MAST, SED3k.1MAST or CALS3k.3MAST             
!     These three model files have to be in the same directory as the pr
!                                                                       
!	Uncertainty estimates are ignored                                     
!                                                                       
!     Cathy Constable, June 2009                                        
!                                                                       
!                                                                       
!     uses pieces of code by Monika Korte, Jeremy Bloxham, Andrew Jackso
!     Kathy Whaler, David Barraclough, Rick 0'Connel and Carl de Boor   
!---------------------------------------------------------------------  
!     for details of B-splines see:                                     
!            Carl de Boor "A Practical Guide to Splines"                
!            Springer-Verlag, 1978.                                     
!--------------------------------------------------------------------   
!---------------------------------------------------------------------  
!                                                                       
!     lmax   is maximum degree of spherical harmonics                   
!     nspl   is number of B-splines                                     
!                                                                       
!     gt     is the full array of coefficients in the expansion of      
!               the geomagnetic potential in spherical harmonics        
!               and cubic B-splines                                     
!                                                                       
!     g      is an array of geomagnetic main field coefficients at a    
!               particular point in time                                
!                                                                       
!     gd     is an array of geomagnetic secular variation coefficients  
!               at a particular point in time                           
!                                                                       
!     tknts  is an array of knot points                                 
!                                                                       
!     p      is an array of real associated Legendre polynomials        
!                                                                       
!     dp     is an array of derivatives of associated Legendre          
!               polynomials with respect to colatitude                  
!---------------------------------------------------------------------  
!                                                                       
!     CALLS:    interv   - calculates which knot lies immediately left  
!                          of the current time point                    
!                                                                       
!               bspline  - calculates B-splines at current time point   
!                                                                       
!               bspline1 - calculates first time derivative of          
!                          B-splines at current time point              
!                                                                       
!               b0, b1   - functions required by bspline1               
!                                                                       
!               plmbar   - calculates associated Legendre polynomials   
!                          and derivatives                              
!                                                                       
!               magfdz   - evaluates field model                        
!                                                                       
!               coords   - performs transformation between geodetic and 
!                          geocentric coords                            
!---------------------------------------------------------------------- 
!
      integer :: nspl,jord,npts,epoch
      real*8 :: fac
      real*8 :: alt,time,theta,phi,r,st,ct,sd,cd
      real*8 :: br,x,y,z,h,f,ainc,d,tStartEpoch,tEndEpoch
      character(len=30) :: outfile,modfile,tessel
      
      integer, parameter :: nxyz = 5000
      integer, parameter :: lmax=10
      integer, parameter :: nsplt=402
      
      integer, parameter :: n=lmax*(lmax+2)
      integer, parameter :: np=n*nsplt
      integer, parameter :: nl=(lmax+1)*(lmax+2)/2
      
      integer :: ipp(nxyz)
      
      real*8 :: gt(n,nsplt),dgt(n,nsplt)
      real*8 :: spl(nsplt),tknts(nsplt+4)
      real*8 :: g(n),dg(n)
      
      real*8 :: pt(4,nxyz),alat(nxyz), alon(nxyz)
      real*8 :: p(nl),dp(nl)
      real*8 :: dx(lmax*(lmax+2)),dy(lmax*(lmax+2)),dz(lmax*(lmax+2))
      
      integer :: lm,nm,k,j,i,nleft,flag,io
      
      data jord/4/
      data fac/1.74532925e-2/
      
      write(*,*) '            -Program shtostt- '
      write(*,*) 'converts continuous sh model to stt snapshots'
      write(*,*) 'reads tesselation x,y,z,ipp '
      write(*,*) 'with flag ipp for patch integration '
      write(*,*) 'outputs tesselation x,y,z,br, ipp '
      write(*,*) 'for processing by eflux'
      write(*,*) 'models CALS3k.3MAST, ARCH3k.1MAST or SED3k.1MAST'
      write(*,*) 'with our MAST estimation procedure for'
      write(*,*) 'the model coefficients.'
      write(*,*) 'Results are written to a plain text output file.'
      write(*,*)
      write(*,*) 'Choose model: 1 - CALS3k.3MAST'
      write(*,*) '              2 - ARCH3k.1MAST'
      write(*,*) '              3 - SED3k.1MAST'
      read(*,*) flag
      write(*,*) 'Give tesselation file name:'
      read(*,*) tessel
      write(*,*) 'Give output file name:'
      read(*,*) outfile
      write(*,*) 'Model Epoch:'
      read(*,*)  epoch
      if (flag.eq.1) then
            modfile='CALS3k.3MAST'
      else if (flag.eq.2) then
            modfile='ARCH3k.1MAST'
      else if (flag.eq.3) then
            modfile='SED3k.1MAST'
      else 
            write(*,*) 'ERROR: invalid model choice'
            stop        
      end if
      
      !********************************************************************
      !     read model, ignore block of uncertainties at the end
      open(7,file=modfile)
      
      read(7,*) tStartEpoch,tEndEpoch
      read(7,*) lm,nm,nspl,(tknts(i),i=1,nspl+4)
      read(7,*) gt
      !     read(7,*) dgt
      close(7)
      
      open(11,file=outfile)
      !
      !*********************************************************************
      !     Read teseslation and ipp patch configuration, convert x,y,z to lat,long
      !
      open(9, file=tessel)
      npts=0
      do j=1,nxyz
            read(9,*, IOSTAT=io) pt(1,j),pt(2,j),pt(3,j),ipp(j)
            if (io.ne.0) exit
            call xyzll(alat(j),alon(j),pt(1,j))
            npts=npts + 1
      enddo
      write (*,*) npts, ' read from tessel, ',tessel
      
      !********************************************************************
      
      i=epoch
      time = dfloat(i)
      ! not sure we need to print every single line
      write (*,*) 'Epoch ', time
      alt=0.0
      !-----
      !     calculate main field coefficients and uncertainties at time time
      call interv(tknts,time,nspl,nleft)
      call bspline(tknts,time,nspl,jord,nleft,spl(nleft-3))
      
      do  k=1,n
            g(k)=0.0
            dg(k)=0.0
            do j=1,4
                  g(k) = g(k) + spl(j+nleft-4)*gt(k,j+nleft-4)
                  dg(k)=dg(k) + spl(j+nleft-4)*dgt(k,j+nleft-4)
            enddo 
      enddo 
      
      ! Evaluate model at all tesselation points
      do j=1,npts
            theta = (90.0-alat(j))*fac
            phi = alon(j)*fac
            !     call coords(alt,theta,r,sd,cd)      
            !     Set radius to CMB
            cd=1.d0
            sd=0.d0
            r= 3485.d0
            
            st=sin(theta)
            ct=cos(theta)
            call plmbar(p,dp,ct,lmax)
            
            call magfdz(p,dp,theta,phi,r,lmax,g,dx,dy,dz,x,y,z,h,f, &
            & ainc,d,sd,cd)
                  br=-z
            ! call errfdz(p,dp,theta,phi,r,lmax,dg,dx,dy,dz,ex,ey,ez,sd,cd)
            
            ! eh=sqrt((1/h**2)*((x*ex)**2+(y*ey)**2))
            ! ef=sqrt((1/f**2)*((x*ex)**2+(y*ey)**2+(z*ez)**2))
            ! ei=sqrt((1/(1+(z/h)**2))**2*((ez/h)**2+((z*eh)/h**2)**2))
            ! ed=sqrt((1/(1+(y/x)**2))**2*((ey/x)**2+((y*ex)/x**2)**2))
            
            ! ainc=ainc/fac
            ! ei=ei/fac
            ! d=d/fac
            ! ed=ed/fac
            ! f=f/1000.
            ! ef=ef/1000.
            
            ! write(11,6200) i,d,ainc,f,ed,ei,ef
            write(11,*)(pt(k,j),k=1,3),br,ipp(j),time
      end do
      close(11)
end program shtostt
!
!
! _______________________________________________________________________________
!
!     
subroutine errfdz(p,dp,theta,phi,r,lmax,g,dx,dy,dz,x,y,z,sd,cd)
!     calculate uncertainties for field predictions from uncertainties
!     in coefficients g
!     
!     use error propagation rules
      implicit real*8 (a-h,o-z)
      dimension  :: g(lmax*(lmax+2))
      dimension  :: dx(lmax*(lmax+2)),dy(lmax*(lmax+2)),dz(lmax*(lmax+2))
      dimension  :: p((lmax+1)*(lmax+2)/2),dp((lmax+1)*(lmax+2)/2)
      
      
      b=6371.2/r
      x=0.
      y=0.
      z=0.
      st=sin(theta)
      if(abs(sinth).lt.1.e-10) sinth=1.e-10
      
      do l=1,lmax
            
            l1=l+1
            bb=b**(l+2)
            k=l*l
            k1=(l*l1)/2+1
            
            dx(k)=dp(k1)*bb
            dy(k)=0.
            dz(k)=-p(k1)*l1*bb
            x=x+(g(k)*dx(k))**2
            z=z+(g(k)*dz(k))**2
            
            do m=1,l
                  
                  t=float(m)*phi
                  k=l*l+2*m-1
                  k1=(l*l1)/2+m+1
                  sint=sin(t)
                  cost=cos(t)
                  
                  dxd = dp(k1)*bb
                  dx(k) = dxd*cost
                  dx(k+1) = dxd*sint
                  x = x + (g(k)*dx(k))**2 + (g(k+1)*dx(k+1))**2
                  
                  dxy = m*p(k1)*bb/sinth
                  dy(k) = dxy*sint
                  dy(k+1) = -dxy*cost
                  y = y + (g(k)*dy(k))**2 + (g(k+1)*dy(k+1))**2
                  
                  dzd = -l1*p(k1)*bb
                  dz(k) = dzd*cost
                  dz(k+1) = dzd*sint
                  z = z + (g(k)*dz(k))**2 + (g(k+1)*dz(k+1))**2
                  
            enddo
      enddo
                  
      xs = x
      x = x*(cd**2) + z*(sd**2)
      z = z*(cd**2) - xs*(sd**2)
                  
      x=sqrt(x)
      y=sqrt(y)
      z=sqrt(z)
                  
      do k=1,lmax*(lmax+2)
            dxk = dx(k)
            dzk = dz(k)
            dx(k) = dxk*cd + dzk*sd
            dz(k) = dzk*cd - dxk*sd
      enddo
                        
      return
end subroutine errfdz
                        
!--------------------------------------------------------------------------      
                              
subroutine interv(tknts,time,nspl,nleft)
      implicit real*8 (a-h,o-z)
      
      dimension tknts(nspl+4)
      
      if(time.lt.tknts(4).or.time.gt.tknts(nspl+1)) return
      
      do n=5,nspl+1
            if(time.le.tknts(n)) then
                  nleft=n-1
                  exit
            endif
      enddo 
            
      return
end subroutine interv
            
!-------------------------------------------------------------------

subroutine bspline(tknts,t,nspl,jorder,nleft,spl)
      
!     calculate splines of order jorder where 1 <= jorder <= 4
      implicit real*8 (a-h,o-z)
      dimension :: tknts(nspl+4)
      dimension :: spl(4)
      
      dimension :: deltal(4),deltar(4)
      
      spl(1)=1.0
      
      do j=1,jorder-1
            
            deltar(j) = tknts(nleft+j) - t
            deltal(j) = t - tknts(nleft+1-j)
            saved=0.0
            
            do i=1,j
                  term = spl(i)/(deltar(i)+deltal(j+1-i))
                  spl(i) = saved + deltar(i)*term
                  saved = deltal(j+1-i)*term
            enddo
                  spl(j+1) = saved
      enddo
      return
end subroutine bspline
                  
!-----------------------------------------------------------------

subroutine coords(h,theta,r,sd,cd)
      implicit real*8 (a-h,o-z)
      
      pi=3.14159265
      b1=40680925.
      b2=40408585.
      theta=pi/2-theta
      clat=cos(theta)
      slat=sin(theta)
      one=b1*clat*clat
      two=b2*slat*slat
      three=one+two
      four=sqrt(three)
      r=sqrt(h*(h+2.*four)+(b1*one+b2*two)/three)
      cd=(h+four)/r
      sd=(b1-b2)/four*slat*clat/r
      st=slat*cd-clat*sd
      ct=clat*cd+slat*sd
      theta=pi/2.-atan2(sinth,costh)
      return
end subroutine coords
      
!-----------------------------------------------------------------

subroutine magfdz(p,dp,theta,phi,r,lmax,g,dx,dy,dz,x,y,z,h,f,i,d,sd,cd)
! ***************************************************************
!      j bloxham  8 nov 1982 & 11 oct 1983
!      modified version of dg13.sv.fort:magfd & jb62.sv.progs:magfds
!      gives field components at radius r
! ***************************************************************
!      this version 16 jan 87
!      saves dx dy dz in computation
! ======================================================================
      
      implicit real*8 (a-h,o-z)
      dimension :: g(lmax*(lmax+2))
      dimension :: dx(lmax*(lmax+2)),dy(lmax*(lmax+2)),dz(lmax*(lmax+2))
      dimension :: p((lmax+1)*(lmax+2)/2),dp((lmax+1)*(lmax+2)/2)
      real*8 :: i
      
      
      b=6371.2/r
      x=0.
      y=0.
      z=0.
      sinth=sin(theta)
      if(abs(sinth).lt.1.e-10) sinth=1.e-10
      
      do l=1,lmax
            
            l1=l+1
            bb=b**(l+2)
            k=l*l
            k1=(l*l1)/2+1
            
            dx(k)=dp(k1)*bb
            dy(k)=0.
            dz(k)=-p(k1)*l1*bb
            x=x+g(k)*dx(k)
            z=z+g(k)*dz(k)
            
            do m=1,l
                  
                  t=float(m)*phi
                  k=l*l+2*m-1
                  k1=(l*l1)/2+m+1
                  sint=sin(t)
                  cost=cos(t)
                  
                  dxd = dp(k1)*bb
                  dx(k) = dxd*cost
                  dx(k+1) = dxd*sint
                  x = x + (g(k)*dx(k)) + (g(k+1)*dx(k+1))
                  
                  dxy = m*p(k1)*bb/sinth
                  dy(k) = dxy*sint
                  dy(k+1) = -dxy*cost
                  y = y + (g(k)*dy(k)) + (g(k+1)*dy(k+1))
                  
                  dzd = -l1*p(k1)*bb
                  dz(k) = dzd*cost
                  dz(k+1) = dzd*sint
                  z = z + (g(k)*dz(k)) + (g(k+1)*dz(k+1))
                  
            enddo
      enddo
                  
      xs = x
      x = x*cd + z*sd
      z = z*cd - xs*sd
                  
      do k=1,lmax*(lmax+2)
            dxk = dx(k)
            dzk = dz(k)
            dx(k) = dxk*cd + dzk*sd
            dz(k) = dzk*cd - dxk*sd
      enddo
                        
                        
      h=sqrt(x*x+y*y)
      f=sqrt(h*h+z*z)
      i=asin(z/f)
      d=atan2(y,x)
                        
      return
end subroutine magfdz
                        
!---------------------------------------------------------------------

subroutine plmbar(p,dp,z,lmax)
      
!    evaluates normalized associated legendre function p(l,m) as function of
!    z=cos(colatitude) using recurrence relation starting with p(l,l) 
!    and then increasing l keeping m fixed.  normalization is: 
!    integral(y(l,m)*y(l,m))=4.*pi, where y(l,m) = p(l,m)*exp(i*m*longitude),
!    which is incorporated into the recurrence relation. p(k) contains p(l,m)
!    with k=(l+1)*l/2+m+1; i.e. m increments through range 0 to l before 
!    incrementing l. routine is stable in single and double precision to
!    l,m = 511 at least; timing proportional to lmax**2
!    r.j.o'connell 7 sept. 1989

!    a.jackson 19 october 1989  code added at end:
!    (2) derivatives added and stored in dp(k)
!        using same arrangement as for p(k)
      
      implicit real*8(a-h,o-z)
      dimension :: p(*),dp(*)
!     dimension of p, dp must be (lmax+1)*(lmax+2)/2 in calling program
      
      if (lmax.lt.0.or.abs(z).gt.1.d0) then
            print *, 'bad arguments p, dp, z, lmax in plmbar' 
      endif
!           case for p(l,0)
      pm2=1.d0
      p(1)=1.d0
      dp(1)=0.d0
      if (lmax .eq. 0) return
      pm1=z
      p(2)=dsqrt(3.d0)*pm1
      k=2
      
      do l=2,lmax
            k=k+l
            plm=(dfloat(2*l-1)*z*pm1-dfloat(l-1)*pm2)/dfloat(l)
            p(k)=dsqrt(dfloat(2*l+1))*plm
            pm2=pm1
      enddo
            pm1=plm
!           case for m > 0
            pmm = 1.d0
            sintsq = (1.d0-z)*(1.d0+z)
            fnum = -1.0d0
            fden = 0.0d0
            kstart = 1
      do m =1 ,lmax
!           case for p(m,m) 
            kstart = kstart+m+1
            fnum = fnum+2.0d0
            fden = fden+2.0d0
            pmm = pmm*sintsq*fnum/fden
            pm2 = dsqrt(dfloat(4*m+2)*pmm)
            p(kstart) = pm2
            
            if (m .eq. lmax) exit

!           case for p(m+1,m)
            pm1=z*dsqrt(dfloat(2*m+3))*pm2
            k = kstart+m+1
            p(k) = pm1
!           case for p(l,m) with l > m+1
            if (m .lt. (lmax-1)) then
                  do l = m+2,lmax
                        k = k+l
                        f1=dsqrt(dfloat((2*l+1)*(2*l-1))/dfloat((l+m)*(l-m)))
                        f2=dsqrt(dfloat((2*l+1)*(l-m-1)*(l+m-1))/dfloat((2*l-3)*(l+m)*(l-m)))
                        plm=z*f1*pm1-f2*pm2
                        p(k) = plm
                        pm2 = pm1
                  enddo          
            pm1 = plm
            endif
      enddo
                        
!       Gauss-Schmidt normalisation:
      k=1
      do l=1,lmax
            fac=1.d0/dsqrt(dfloat(2*l+1))
            do m=0,l
                  k=k+1
                  p(k)=p(k)*fac
            enddo
      enddo
                  
!     now find derivatives of p(z) wrt theta, where z=cos(theta)
      dp(2)=-p(3)
      dp(3)=p(2)
      k=3
      do l=2,lmax
            
            k=k+1
!         treat m=0 and m=l separately
            dp(k)=-dsqrt(dfloat(l*(l+1))/2.d0)*p(k+1)
            dp(k+l)=dsqrt(dfloat(l)/2.d0)*p(k+l-1)
            do m=1,l-1
                  k=k+1
                  fac1=dsqrt( dfloat( (l-m)*(l+m+1) ) )
                  fac2=dsqrt( dfloat( (l+m)*(l-m+1) ) )
                  if(m.eq.1)fac2=fac2*dsqrt(2.d0)
                  dp(k)=0.5d0*( fac2*p(k-1) - fac1*p(k+1) )
            enddo
            k=k+1
                  
      enddo
      return
end subroutine plmbar
                              
!------------------------------------------------------------------

subroutine bspline1(tknts,t,nspl,nl,spl)
      
      implicit real*8(a-h,o-z)
      integer, parameter :: nsplt=402
      dimension :: tknts(nspl+4)
      dimension :: spl(nspl),spl1(nsplt)
      
      data jord/3/
      
      if(nspl.gt.nsplt)then
            write(6,*)' increase dimensions of spl1 in bspline1'
            stop
      endif
      
      do is=1,nspl
            spl(is)=0.0
      enddo
            
      call bspline(tknts,t,nspl,jord,nl,spl1(nl-2))
            
      spl(nl) = b0(tknts,nl)*spl1(nl)
      
      spl(nl-1) = b0(tknts,nl-1)*spl1(nl-1) + b1(tknts,nl-1)*spl1(nl)
      spl(nl-2) = b0(tknts,nl-2)*spl1(nl-2) + b1(tknts,nl-2)*spl1(nl-1)
      spl(nl-3) = b1(tknts,nl-3)*spl1(nl-2)
            
      return
end subroutine bspline1
            
function b0(tknts,i)
      implicit real*8(a-h,o-z)
      dimension :: tknts(*)
      
      b0 = 3.0/( tknts(i+3) - tknts(i) )
      
      return
end function b0      
      
function b1(tknts,i)
      implicit real*8(a-h,o-z)
      dimension :: tknts(*)
      
      b1 = -3.0/( tknts(i+4) - tknts(i+1) )
      
      return
end function b1
            
!____________________________________________________
subroutine xyzll(rlat,rlong,x)
!     calls no other routines
!     Maps unit vector x into lat, long (degrees)

      real*8  :: x(3),rlat,rlong,toDeg
      data toDeg/57.29577951/
      rlat=toDeg*atan2(x(3),sqrt(x(2)**2 +x(1)**2))
      rlong =toDeg*atan2(x(2),x(1))

      return
end subroutine xyzll