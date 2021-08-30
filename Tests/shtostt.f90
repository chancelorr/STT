!*==INTERV.spg  processed by SPAG 6.72Dc at 07:24 on 30 Aug 2021
!   by Chance, adapted from C. Constable
!  	reads snaphsot sh field in form of 1,m,glm,hlm and tesselation ifn form of x,y,z, ipp
!  	evlautes br on CMB and outputs to corefile for digestion by eflux
 
!       Modified version of fieldpred.f to read stt tesselation plus ipp zone
!   	output is stt, br,ipp
!       from models ARCH3k.1MAST, SED3k.1MAST or CALS3k.3MAST
!       These three model files have to be in the same directory as the program.
 
!  	Uncertainty estimates are ignored
 
!       Cathy Constable, June 2009
 
 
!       uses pieces of code by Monika Korte, Jeremy Bloxham, Andrew Jackson, David Gubbins
!       Kathy Whaler, David Barraclough, Rick 0'Connel and Carl de Boor
!  ---------------------------------------------------------------------
!       for details of B-splines see: 
!              Carl de Boor "A Practical Guide to Splines"
!              Springer-Verlag, 1978.
!  --------------------------------------------------------------------
!  ---------------------------------------------------------------------
 
!       lmax   is maximum degree of spherical harmonics
!       nspl   is number of B-splines
 
!       gt     is the full array of coefficients in the expansion of
!                 the geomagnetic potential in spherical harmonics
!                 and cubic B-splines
 
!       g      is an array of geomagnetic main field coefficients at a
!                 particular point in time
 
!       gd     is an array of geomagnetic secular variation coefficients
!                 at a particular point in time
 
!       tknts  is an array of knot points
 
!       p      is an array of real associated Legendre polynomials
 
!       dp     is an array of derivatives of associated Legendre
!                 polynomials with respect to colatitude
!  ---------------------------------------------------------------------
 
!       CALLS:    interv   - calculates which knot lies immediately left
!                            of the current time point
 
!                 bspline  - calculates B-splines at current time point
 
!                 bspline1 - calculates first time derivative of
!                            B-splines at current time point
 
!                 b0, b1   - functions required by bspline1
 
!                 plmbar   - calculates associated Legendre polynomials
!                            and derivatives
 
!                 magfdz   - evaluates field model
 
!                 coords   - performs transformation between geodetic and
!                            geocentric coords
 
!--------------------------------------------------------------------------
 
      SUBROUTINE INTERV(Tknts,Time,Nspl,Nleft)
        IMPLICIT NONE
  !*--INTERV10
  !*** Start of declarations inserted by SPAG
        INTEGER n , Nleft , Nspl
        REAL Time
  !*** End of declarations inserted by SPAG
  !      implicit real*8 (a-h,o-z)
   
        REAL*8 Tknts(Nspl+4)
   
        IF ( Time.LT.Tknts(4) .OR. Time.GT.Tknts(Nspl+1) ) RETURN
   
        DO n = 5 , Nspl + 1
           IF ( Time.LE.Tknts(n) ) THEN
              Nleft = n - 1
              GOTO 99999
           ENDIF
        ENDDO
   
   
  99999 END
  !*==BSPLINE.spg  processed by SPAG 6.72Dc at 07:24 on 30 Aug 2021
   
  !-------------------------------------------------------------------
   
        SUBROUTINE BSPLINE(Tknts,T,Nspl,Jorder,Nleft,Spl)
        IMPLICIT NONE
  !*--BSPLINE36
  !*** Start of declarations inserted by SPAG
        INTEGER i , j , Jorder , Nleft , Nspl
        REAL saved , term
        REAL*8 T
  !*** End of declarations inserted by SPAG
   
  ! calculate splines of order jorder where 1 <= jorder <= 4
  !       implicit real*8 (a-h,o-z)
        REAL*8 Tknts(Nspl+4)
        REAL*8 Spl(4)
   
        REAL*8 deltal(4) , deltar(4)
   
        Spl(1) = 1.0
   
        DO j = 1 , Jorder - 1
   
           deltar(j) = Tknts(Nleft+j) - T
           deltal(j) = T - Tknts(Nleft+1-j)
           saved = 0.0
   
           DO i = 1 , j
              term = Spl(i)/(deltar(i)+deltal(j+1-i))
              Spl(i) = saved + deltar(i)*term
              saved = deltal(j+1-i)*term
           ENDDO
   
           Spl(j+1) = saved
   
        ENDDO
   
   
   
        END
  !*==COORDS.spg  processed by SPAG 6.72Dc at 07:24 on 30 Aug 2021
   
  !-----------------------------------------------------------------
   
        SUBROUTINE COORDS(H,Theta,R,Sd,Cd)
        IMPLICIT NONE
  !*--COORDS76
  !*** Start of declarations inserted by SPAG
        REAL*8 b1 , b2 , Cd , clat , costh , four , H , one , pi , R ,    &
             & Sd , sinth , slat , Theta , three , two
  !*** End of declarations inserted by SPAG
   
        pi = 3.14159265
        b1 = 40680925.
        b2 = 40408585.
        Theta = pi/2 - Theta
        clat = COS(Theta)
        slat = SIN(Theta)
        one = b1*clat*clat
        two = b2*slat*slat
        three = one + two
        four = SQRT(three)
        R = SQRT(H*(H+2.*four)+(b1*one+b2*two)/three)
        Cd = (H+four)/R
        Sd = (b1-b2)/four*slat*clat/R
        sinth = slat*Cd - clat*Sd
        costh = clat*Cd + slat*Sd
        Theta = pi/2. - ATAN2(sinth,costh)
        END
  !*==MAGFDZ.spg  processed by SPAG 6.72Dc at 07:24 on 30 Aug 2021
   
  !-----------------------------------------------------------------
   
        SUBROUTINE MAGFDZ(P,Dp,Theta,Phi,R,Lmax,G,Dx,Dy,Dz,X,Y,Z,H,F,I,D, &
                        & Sd,Cd)
        IMPLICIT NONE
  !*--MAGFDZ106
  !*** Start of declarations inserted by SPAG
        REAL b , bb , Cd , cost , D , dxd , dxk , dxy , dzd , dzk , F ,   &
           & H , Phi , R , Sd , sint , sinth , t , Theta , X
        REAL xs , Y , Z
        INTEGER k , k1 , l , l1 , Lmax , m
  !*** End of declarations inserted by SPAG
   
  !
  !***************************************************************
  !
  !     j bloxham  8 nov 1982 & 11 oct 1983
  !
  !     modified version of dg13.sv.fort:magfd & jb62.sv.progs:magfds
  !
  !     gives field components at radius r
  !
  !***************************************************************
  !
  !
  !     this version 16 jan 87
  !
  !     saves dx dy dz in computation
  !
  !c======================================================================
   
  !      implicit real*8 (a-h,o-z)
        REAL*8 G(Lmax*(Lmax+2))
        REAL*8 Dx(Lmax*(Lmax+2)) , Dy(Lmax*(Lmax+2)) , Dz(Lmax*(Lmax+2))
        REAL*8 P((Lmax+1)*(Lmax+2)/2) , Dp((Lmax+1)*(Lmax+2)/2)
        REAL*8 I
   
   
        b = 6371.2/R
        X = 0.
        Y = 0.
        Z = 0.
        sinth = SIN(Theta)
        IF ( ABS(sinth).LT.1.E-10 ) sinth = 1.E-10
   
        DO l = 1 , Lmax
   
           l1 = l + 1
           bb = b**(l+2)
           k = l*l
           k1 = (l*l1)/2 + 1
   
           Dx(k) = Dp(k1)*bb
           Dy(k) = 0.
           Dz(k) = -P(k1)*l1*bb
           X = X + G(k)*Dx(k)
           Z = Z + G(k)*Dz(k)
   
           DO m = 1 , l
   
              t = FLOAT(m)*Phi
              k = l*l + 2*m - 1
              k1 = (l*l1)/2 + m + 1
              sint = SIN(t)
              cost = COS(t)
   
              dxd = Dp(k1)*bb
              Dx(k) = dxd*cost
              Dx(k+1) = dxd*sint
              X = X + (G(k)*Dx(k)) + (G(k+1)*Dx(k+1))
   
              dxy = m*P(k1)*bb/sinth
              Dy(k) = dxy*sint
              Dy(k+1) = -dxy*cost
              Y = Y + (G(k)*Dy(k)) + (G(k+1)*Dy(k+1))
   
              dzd = -l1*P(k1)*bb
              Dz(k) = dzd*cost
              Dz(k+1) = dzd*sint
              Z = Z + (G(k)*Dz(k)) + (G(k+1)*Dz(k+1))
   
           ENDDO
        ENDDO
   
        xs = X
        X = X*Cd + Z*Sd
        Z = Z*Cd - xs*Sd
   
        DO k = 1 , Lmax*(Lmax+2)
           dxk = Dx(k)
           dzk = Dz(k)
           Dx(k) = dxk*Cd + dzk*Sd
           Dz(k) = dzk*Cd - dxk*Sd
        ENDDO
   
   
        H = SQRT(X*X+Y*Y)
        F = SQRT(H*H+Z*Z)
        I = ASIN(Z/F)
        D = ATAN2(Y,X)
   
        END
  !*==PLMBAR.spg  processed by SPAG 6.72Dc at 07:24 on 30 Aug 2021
   
  !---------------------------------------------------------------------
   
        SUBROUTINE PLMBAR(P,Dp,Z,Lmax)
        IMPLICIT NONE
  !*--PLMBAR209
  !*** Start of declarations inserted by SPAG
        REAL f1 , f2 , fac , fac1 , fac2 , fden , fnum , plm , pm1 , pm2 ,&
           & pmm , sintsq , Z
        INTEGER k , kstart , l , Lmax , m
  !*** End of declarations inserted by SPAG
  !
  !  evaluates normalized associated legendre function p(l,m) as function of
  !   z=cos(colatitude) using recurrence relation starting with p(l,l)
  !   and then increasing l keeping m fixed.  normalization is:
  !   integral(y(l,m)*y(l,m))=4.*pi, where y(l,m) = p(l,m)*exp(i*m*longitude),
  !   which is incorporated into the recurrence relation. p(k) contains p(l,m)
  !   with k=(l+1)*l/2+m+1; i.e. m increments through range 0 to l before
  !   incrementing l. routine is stable in single and double precision to
  !   l,m = 511 at least; timing proportional to lmax**2
  !   r.j.o'connell 7 sept. 1989
   
  !   a.jackson 19 october 1989  code added at end:
  !   (2) derivatives added and stored in dp(k)
  !       using same arrangement as for p(k)
  !
  !      implicit real*8(a-h,o-z)
        REAL*8 P(*) , Dp(*)
  !     --dimension of p, dp must be (lmax+1)*(lmax+2)/2 in calling program
        IF ( Lmax.LT.0 .OR. ABS(Z).GT.1.D0 ) PAUSE 'bad arguments'
  !       --case for p(l,0)
        pm2 = 1.D0
        P(1) = 1.D0
        Dp(1) = 0.D0
        IF ( Lmax.EQ.0 ) RETURN
        pm1 = Z
        P(2) = DSQRT(3.D0)*pm1
        k = 2
        DO l = 2 , Lmax
           k = k + l
           plm = (DFLOAT(2*l-1)*Z*pm1-DFLOAT(l-1)*pm2)/DFLOAT(l)
           P(k) = DSQRT(DFLOAT(2*l+1))*plm
           pm2 = pm1
           pm1 = plm
        ENDDO
  !       --case for m > 0
        pmm = 1.D0
        sintsq = (1.D0-Z)*(1.D0+Z)
        fnum = -1.0D0
        fden = 0.0D0
        kstart = 1
        DO m = 1 , Lmax
  !         --case for p(m,m)
           kstart = kstart + m + 1
           fnum = fnum + 2.0D0
           fden = fden + 2.0D0
           pmm = pmm*sintsq*fnum/fden
           pm2 = DSQRT(DFLOAT(4*m+2)*pmm)
           P(kstart) = pm2
           IF ( m.EQ.Lmax ) GOTO 100
  !         --case for p(m+1,m)
           pm1 = Z*DSQRT(DFLOAT(2*m+3))*pm2
           k = kstart + m + 1
           P(k) = pm1
  !         --case for p(l,m) with l > m+1
           IF ( m.LT.(Lmax-1) ) THEN
              DO l = m + 2 , Lmax
                 k = k + l
                 f1 = DSQRT(DFLOAT((2*l+1)*(2*l-1))/DFLOAT((l+m)*(l-m)))
                 f2 = DSQRT(DFLOAT((2*l+1)*(l-m-1)*(l+m-1))               &
                    & /DFLOAT((2*l-3)*(l+m)*(l-m)))
                 plm = Z*f1*pm1 - f2*pm2
                 P(k) = plm
                 pm2 = pm1
                 pm1 = plm
              ENDDO
           ENDIF
        ENDDO
   
   
  !       Gauss-Schmidt normalisation:
   100  k = 1
        DO l = 1 , Lmax
           fac = 1.D0/DSQRT(DFLOAT(2*l+1))
           DO m = 0 , l
              k = k + 1
              P(k) = P(k)*fac
           ENDDO
        ENDDO
   
  !       now find derivatives of p(z) wrt theta, where z=cos(theta)
        Dp(2) = -P(3)
        Dp(3) = P(2)
        k = 3
        DO l = 2 , Lmax
   
           k = k + 1
  !         treat m=0 and m=l separately
           Dp(k) = -DSQRT(DFLOAT(l*(l+1))/2.D0)*P(k+1)
           Dp(k+l) = DSQRT(DFLOAT(l)/2.D0)*P(k+l-1)
           DO m = 1 , l - 1
              k = k + 1
              fac1 = DSQRT(DFLOAT((l-m)*(l+m+1)))
              fac2 = DSQRT(DFLOAT((l+m)*(l-m+1)))
              IF ( m.EQ.1 ) fac2 = fac2*DSQRT(2.D0)
              Dp(k) = 0.5D0*(fac2*P(k-1)-fac1*P(k+1))
           ENDDO
           k = k + 1
   
        ENDDO
        END
  !*==BSPLINE1.spg  processed by SPAG 6.72Dc at 07:24 on 30 Aug 2021
   
  !------------------------------------------------------------------
   
        SUBROUTINE BSPLINE1(Tknts,T,Nspl,Nl,Spl)
        IMPLICIT NONE
  !*--BSPLINE1321
  !*** Start of declarations inserted by SPAG
        REAL*8 B0 , B1 , T
        INTEGER is , jord , Nl , Nspl , NSPLT
  !*** End of declarations inserted by SPAG
   
  !      implicit real*8(a-h,o-z)
        PARAMETER (NSPLT=402)
        REAL*8 Tknts(Nspl+4)
        REAL*8 Spl(Nspl) , spl1(NSPLT)
   
        DATA jord/3/
   
        IF ( Nspl.GT.NSPLT ) THEN
           WRITE (6,*) ' increase dimensions of spl1 in bspline1'
           STOP
        ENDIF
   
        DO is = 1 , Nspl
           Spl(is) = 0.0
        ENDDO
   
        CALL BSPLINE(Tknts,T,Nspl,jord,Nl,spl1(Nl-2))
   
        Spl(Nl) = B0(Tknts,Nl)*spl1(Nl)
   
        Spl(Nl-1) = B0(Tknts,Nl-1)*spl1(Nl-1) + B1(Tknts,Nl-1)*spl1(Nl)
        Spl(Nl-2) = B0(Tknts,Nl-2)*spl1(Nl-2) + B1(Tknts,Nl-2)*spl1(Nl-1)
        Spl(Nl-3) = B1(Tknts,Nl-3)*spl1(Nl-2)
   
        END
  !*==B0.spg  processed by SPAG 6.72Dc at 07:24 on 30 Aug 2021
   
        FUNCTION B0(Tknts,I)
        IMPLICIT NONE
  !*--B0356
  !*** Start of declarations inserted by SPAG
        REAL*8 B0 , Tknts
        INTEGER I
  !*** End of declarations inserted by SPAG
        DIMENSION Tknts(*)
   
        B0 = 3.0/(Tknts(I+3)-Tknts(I))
   
        END
  !*==B1.spg  processed by SPAG 6.72Dc at 07:24 on 30 Aug 2021
   
        FUNCTION B1(Tknts,I)
        IMPLICIT NONE
  !*--B1370
  !*** Start of declarations inserted by SPAG
        REAL*8 B1 , Tknts
        INTEGER I
  !*** End of declarations inserted by SPAG
        DIMENSION Tknts(*)
   
        B1 = -3.0/(Tknts(I+4)-Tknts(I+1))
   
        END
  !*==XYZLL.spg  processed by SPAG 6.72Dc at 07:24 on 30 Aug 2021
   
  !____________________________________________________
        SUBROUTINE XYZLL(Rlat,Rlong,X)
        IMPLICIT NONE
  !*--XYZLL385
  !$$$$ calls no other routines
  !  Maps unit vector x into lat, long (degrees)
        REAL*8 X(3) , Rlat , Rlong , drad
        DATA drad/57.29577951/
        Rlat = drad*ATAN2(X(3),SQRT(X(2)**2+X(1)**2))
        Rlong = drad*ATAN2(X(2),X(1))
        END
  !_____________________________________________________________
  