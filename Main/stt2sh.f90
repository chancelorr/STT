PROGRAM STTTOSH
! Written by C. Constable
! Modified by Chancelor Roberts
! Updating July 2022 to produce a program for creating stt core field models then translating them to Schmidt normalized SH models at Earth Surface
! Uncompleted modification in June 2009 at GFZ- started form fflux7.f, stripping out forward and inverse parts
    IMPLICIT NONE
    
    DOUBLE PRECISION :: A , AMAx , AMIn , BC , BRO , BS , DLM , DMAx ,   &
                    & DMIn , fluxd0 , fluxt , FLUxv , GLM , GRAdm ,    &
                    & HLM , P , PLM , r , rc , REG
    DOUBLE PRECISION :: RO , SIGma , t, time , V , val , XYW
    double precision :: bColat1 , bColat2 , bLong1 , bLong2, r0, lat, lon
    INTEGER :: i , i1 , i2 , if , IFLag , INP , IOUt , ip , IPP ,        &
        & IPRint , IPTs , ISIgn , ISTack , ITRis , IV , j , k ,     &
        & kdim , kk , l
    INTEGER :: ldeg , LLPt , LMAX , LOBs , LPTr , LPTs , LTRis , m ,     &
        & MAXE , MAXOB , MAXP , MOBs , N , NEAr , NEPoch , nfound , &
        & NOBs , NP , NPAtch , NPTs
    INTEGER :: NQUad , NTRia , NTRis , NVERT , nwant1 , nwant2 , NXYZ ,  &
        & NXYZ2 , NXYZ6
! reads position cooordinates on the core, creates a tesselation as
! nearly equiangular as possible using dorito, and finds the nearest
! neighbours of every point
! Option exists to read predefined tesselation
!
!
! Next comments are dubious
! Evaluates flux integrals through patches on core surrounded by null flux curves
! Also will calculate flux (signed and unsigend?) through predefined regions on core (and perhaps Earth?)
! surface to look at hemispheric asymmetries, temporal evolution, etc.
!
! see readme for what is needed for conversion of sh models to input here
!
!
!******calls sbody, nneigh, triq, obs, bfield, gram, qr
!
!	nxyz:    maximum no. of points sampled on the core
!	nxyz2	 =2*nxyz
!	nvert	 maximum no. of nearest neighbours  allowed for any
!                point on the core
!	maxob    maximum no. of observations at surface
!
!	in /tessel/
!
!	n:       no. of points on the core
!	ntria:	no of triangles on core =2*n-4
!	p(i,j):  i=1,3 gives ith coordinate of jth point
!		 p(4,j) gives br at jth point
!	iv(i,k): point number of ith vertex of k-th triangle
!	near(i,j): i=1,nvert gives nearest neighbours of j-th point in
!		   p, ordered sequentially around p.
!
!	in /observ/
!	mobs:    total no. of observation points at Earth surface
!	ro(i,j,k): i=1,3 position of jth observation point, at kth
!		   epoch (r,theta,phi)
!	bro(J,k):  jth observation of a component of the magnetic field
!                  point, at kth epoch
!	iflag(j,k): flag identifying component type for jth observation
!	          at kth epoch
!	sigma(j,k):   rms misfit level for datum j at epoch k
!	nepoch:	  number of epochs
!	nobs(k):  number of observation sites at epoch k
!	lobs(k):  number of components at epoch k
!
!	in /calcul/
!
!	bc(j):  computed radial field at jth core point
!	bs(j,k):  jth computed field component at epoch k
!	a(i,j):  matrix for forward calculations, i=1,mobs,j=1,nn;
!		 also used to tore upper triangle in lowermost nxyz
!		 rows
!
!	in /rule/
!	integration rule for right triangle at origin
!	nquad:    no. of points to evaluate the function
!	xyw(i,j): i=1,2 gives x and y coordinate of point for function
!		  evaluation and i=3 gives weight of that point for
!		  integration rule.
!
!	in /patchy/
!	np(i): 	  total no. of patches of one sign for epoch i
!	npts(j,k):  #vertices (points) in patch j, for epoch k
!	isign(j,k)  sign of patch j, for epoch k
!	ipts(j):  starting index in lpts of the list of points in
!		  patch j
!  	lpts(ipts(j)+k-1):    kth vertex in patch j
!	ntris(j):             # triangles in patch j
!  	itris(j):             starting index in ltris of the list of
!                        triangles in patch j
!   	ltris(itris(j)+k-1):  kth triangle in patch j
!	lptr(j,k):	k=1:  positive patch (if any)associated with
!		  triangle j.  k=2:negative patch (if any)
!       ???????llpt(maxp,maxe),npatch,ipp(nxyz)??????
!
!	in /patchg/
!	fluxv(j,k):  flux through patch j, for epoch k
!	gradm(j,k,i): gradient of flux through patch j due to change in
!		field at vertex k, for epoch i
!
!
!	in /circle/
!
!	v(i,k):  i=1,4 contains coordinates of circumcenter
!	 	 and "radius" of circle circumscribing triangle k in
!		 list iv
!	istack:  is a work array of dimension at least nxyz2
!
!
!	in /shr/
!
!	glm(l+1,m+1), hlm(l+1,m+1):	partially normalized Schmidt
!	coeffs of degree l, order m
!	plm(l+1,m+1):  Legendre polynomial of degree l, order m
!
!
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ,NXYZ6=6*NXYZ,NVERT=10)
    PARAMETER (MAXOB=7000,MAXP=30,LMAX=101,MAXE=200)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    COMMON /OBSERV/ RO(3,MAXOB,MAXE) , BRO(MAXOB,MAXE) ,              &
                & SIGma(MAXOB,MAXE) , IFLag(MAXOB,MAXE) , MOBs ,    &
                & NOBs(MAXE) , NEPoch , LOBs(MAXE)
    COMMON /CALCUL/ BS(MAXOB,MAXE) , A(MAXOB,NXYZ) , BC(NXYZ2) ,      &
                & REG(NXYZ,NXYZ)
    COMMON /RULE  / XYW(3,100) , NQUad
    COMMON /PATCHY/ NP(MAXE) , NPTs(MAXP,MAXE) , ISIgn(MAXP,MAXE) ,   &
                & IPTs(MAXP,MAXE) , LPTs(2*NXYZ,MAXE) , NTRis(NXYZ) &
                & , ITRis(NXYZ+1) , LTRis(2*(2*NXYZ-4)) ,           &
                & LPTr(2*NXYZ-4,2) , LLPt(MAXP,MAXE) , NPAtch ,     &
                & IPP(NXYZ)
    COMMON /PATCHG/ FLUxv(MAXP,MAXE) , GRAdm(MAXP,NXYZ,MAXE)
    COMMON /SHR   / GLM(LMAX,LMAX) , HLM(LMAX,LMAX) , PLM(LMAX,LMAX) ,&
                & DLM(LMAX,LMAX)
    COMMON /CIRCLE/ V(4,NXYZ2) , ISTack(NXYZ2)
    COMMON /BOUND / DMIn , DMAx , AMIn , AMAx
    CHARACTER*80 :: filnam
    DIMENSION :: val(10)
!      data ((iv(i,j),i=1,3),j=1,8)/1,2,3,1,3,4,1,4,5,1,2,5,
!     $ 2,3,6,3,4,6,4,5,6,2,5,6/,
!     $ ((p(i,j),i=1,3),j=1,6)/0,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1/
    DATA INP , IOUt , IPRint/5 , 6 , 0/ , DMIn , DMAx/1.0 , -1.0/ ,   &
        & i1 , i2/ - 1 , 0/
    DATA rc/3485/
! number of parameters expected in calls to getval
    DATA nwant1 , nwant2/1 , 2/
!
!  Get parameters to run code
    CALL INCORE
!
!
!  Get required core points and create a tesselation using dorito
!
    CALL GETVAL('bodydim',val,nwant1,nfound)
    kdim = NINT(val(1))
!
    CALL SBODY(kdim,NXYZ,N,P,IV,NTRia)
    WRITE (IOUt,*) 'Returned from sbody'
    CALL FLUSH(IOUt)
!
!  Find the nearest neighbours of every point on the core and save in
!  near. This gives support region for basis functions.
!
    CALL NNEIGH
    WRITE (IOUt,*) 'Returned from nneigh'
    CALL FLUSH(IOUt)
!
!  Set up 7 point Stroud integration rule for forward calculations
!
    CALL TRIQ(i1,i2)
    WRITE (IOUt,*) 'Returned from triq'
    CALL FLUSH(IOUt)
    CALL GETVAL('problem',val,nwant1,nfound)
    if = NINT(val(1))
!	if (if.eq.1.or.if.eq.0)then
!
!   Read in core model with patch assignations
    CALL GETVAL('epochs',val,nwant1,nfound)
    NEPoch = NINT(val(1))
    CALL GETCHR('corefield',filnam,nfound)
    OPEN (UNIT=20,FILE=filnam)
    READ (20, *) bColat1 , bColat2 , bLong1 , bLong2
    READ (20, *) r0
    DO kk = 1 , NEPoch
        DO i = 1 , N
        READ (20,*) lat, lon , BC(i) , IPP(i), t
!,time
        P(4,i) = BC(i)
        ENDDO
!
!  Cycle core points into p(4,j)
!
        DO k = 1 , N
        P(4,k) = BC(k)
        ENDDO
!
!  Find SH representation of this core field
        CALL GETVAL('shrep',val,nwant2,nfound)
        ldeg = NINT(val(1))
        r = val(2)
        IF ( nfound.GT.0 .AND. ldeg.GT.0 ) THEN
        CALL SHEXP(ldeg,r)
        OPEN(UNIT=21, FILE='shOut/SHR'//filnam(8:))
        write (21, '(i2, i2, i3, a)', advance='no') ldeg, 0, NEpoch, NEW_LINE('')
        
        WRITE (IOUt,*)                                              &
            &' Writing Spherical Harmonic Representation to: '//'shOut/SHR'//filnam(8:)
        DO l = 0 , ldeg
            DO m = 0 , l
                write (21,'(i2, i10, S, f26.13, f26.13, a)', advance='no') l, m, GLM(l+1,m+1) , HLM(l+1,m+1), NEW_LINE('')
            ENDDO
        ENDDO
        ENDIF
!
!  For each epoch model, when ip=1, find patches on the core where Br has one sign
!
        CALL GETVAL('patches',val,nwant1,nfound)
        ip = NINT(val(1))
        IF ( ip.EQ.1 ) THEN
        WRITE (*,*) 'kk= ' , kk
        IF ( kk.EQ.1 ) THEN
            CALL PATCH(kk)
        ELSE
            NP(kk) = NP(kk-1)
        ENDIF
        WRITE (IOUt,*) 'no of patches = ' , NP(kk)
!  Compute flux through designated area, and axial dipole contribution
!
!  Cycle core points into p(4,j)
!
        DO k = 1 , N
            P(4,k) = BC(k)
        ENDDO
!
        CALL FLUX(kk)
        fluxd0 = 0.0
        fluxt = 0.0
        DO j = 1 , NP(kk)
            WRITE (IOUt,*) ' flux through patch ' , j , ' is ' ,     &
                        & FLUxv(j,kk)*rc**2
!	     fluxt=fluxt + abs(fluxv(j,kk))
            WRITE (IOUt,*) 'no of points' , NPTs(j,kk)
            WRITE (IOUt,*) 'points in this patch are ' ,             &
                        & (LPTs(IPTs(j,kk)+k-1,kk),k=1,NPTs(j,kk))
        ENDDO
!	   write(iout,*)'Total flux for this epoch is ',
!     $      fluxt*rc**2
        WRITE (9,*) time , FLUxv(1,kk)*rc**2 , FLUxv(2,kk)*rc**2
        CALL FLUSH(IOUt)
        ENDIF
    ENDDO
    CLOSE(unit=20)
    END
!_______________________________________________________________________
!_______________________________________________________________________
    FUNCTION P3AREA(A,B,C)
!  computes area of planar triangle specified by 3 vectors
    IMPLICIT NONE

    DOUBLE PRECISION A , ab , ac , ar , B , C , DOTTY , P3AREA
    INTEGER j
    
    DIMENSION A(3) , B(3) , C(3) , ab(3) , ac(3) , ar(3)
    DO j = 1 , 3
        ab(j) = B(j) - A(j)
        ac(j) = C(j) - A(j)
    ENDDO
    CALL AXB(ab,ac,ar)
    P3AREA = SQRT(DOTTY(ar,ar))/2.
    END

!_______________________________________________________________________
    FUNCTION TEGF(F,P1,P2,P3)
    IMPLICIT NONE

    DOUBLE PRECISION aa , ajac , b , bjac , c , F , gam1 , P1 , P2 ,  &
                    & P3 , PAREA , rg , rp , rpa , TEGF , ua , ub ,    &
                    & uc , vc , XYW
    INTEGER INP , IOUt , IPRint , k , NQUad

!$$$$$$$calls centr, gpole, gnomon, f, parea
!  integrates the function f over the spherical triangle defined by
!  vertices p1, p2, p3
!  core field br at vertices p, is in 4th component of arrays
!
!  The quadrature rule generated by 'triq' for integrals over
!  plane triangle with vertices  (0,0,0) (1,0,0) (0,1,0) is assumed
!  to be available in common /rule/.
!
!
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /RULE  / XYW(3,100) , NQUad
    DIMENSION rp(3) , vc(4) , aa(3) , b(3) , c(3) , P1(4) , P2(4) ,   &
            & P3(4) , ua(2) , ub(2) , uc(2) , rg(2) , rpa(3)
!
    TEGF = 0.
!
!  compute centroid vc of triangle
!
    CALL CENTR(P1,P2,P3,vc)
!
!  rotate so vertices are in coord system with vc as pole, p(1,j) at
!  zero longitude
    CALL GPOLE(1,vc,P1,P1,aa)
    CALL GPOLE(1,vc,P1,P2,b)
    CALL GPOLE(1,vc,P1,P3,c)
!
!  find gnomonic projection of vertices and jacobian of transformation
    CALL GNOMON(1,ua,aa)
    CALL GNOMON(1,ub,b)
    CALL GNOMON(1,uc,c)
    ajac = PAREA(ua,ub,uc)*2.
!
!  Run through the cubature rule for this triangle
!
    DO k = 1 , NQUad
!
!  find points on gnomonic surface corresponding to integration points
!
        rg(1) = ua(1) + (ub(1)-ua(1))*XYW(1,k) + (uc(1)-ua(1))*XYW(2,k)
        rg(2) = ua(2) + (ub(2)-ua(2))*XYW(1,k) + (uc(2)-ua(2))*XYW(2,k)
!
!  find points in physical space corres. to rg
!
        CALL GNOMON(-1,rg,rpa)
        CALL GPOLE(-1,vc,P1,rp,rpa)
!
!  evaluate function f at rp, rg, etc.
!
        gam1 = F(rg,rp,k,P1,P2,P3)
!
!  Weight integral according to xyw rule
!  Find Jacobian for gnomonic to spherical projection
!
        bjac = SQRT((1+rg(1)**2+rg(2)**2)**(-3))
        TEGF = TEGF + gam1*XYW(3,k)*ajac*bjac
    ENDDO
    END
!_______________________________________________________________________
    FUNCTION B2INT(Rg,Rp,K,P1,P2,P3)
    IMPLICIT NONE

    DOUBLE PRECISION B2INT , gamma , P1 , P2 , P3 , Rg , Rp , XYW
    INTEGER K , NQUad

!  interpolates b to a point rg on the gnomonic projection
    COMMON /RULE  / XYW(3,100) , NQUad
    DIMENSION Rp(3) , Rg(2) , gamma(3) , P1(4) , P2(4) , P3(4)
    gamma(1) = 1. - XYW(1,K) - XYW(2,K)
    gamma(2) = XYW(1,K)
    gamma(3) = XYW(2,K)
    B2INT = (gamma(1)*P1(4)+gamma(2)*P2(4)+gamma(3)*P3(4))**2
    END

!____________________________________________________________________________________________________
    SUBROUTINE FLUXPT(Llpt,Np,Fluxv)
    IMPLICIT NONE

    DOUBLE PRECISION Fluxv , rc , x
    INTEGER i , INP , IOUt , IPRint , k , Llpt , MAXE , MAXP , Np

!  get a ptr to order fluxv(np) in magnitude
!  ignore patches with no flux i.e. treat as null flux points
    COMMON /IO    / INP , IOUt , IPRint
    PARAMETER (MAXP=30,MAXE=200)
    DIMENSION x(MAXP,MAXE) , Fluxv(Np) , Llpt(Np)
    DATA rc/3485/
    k = 0
    DO i = 1 , Np
        IF ( ABS(Fluxv(i)).GT.1. ) THEN
        k = k + 1
        x(k,1) = (Fluxv(i))
        x(k,2) = FLOAT(i)
        ENDIF
    ENDDO
    Np = k
    CALL SORT22(MAXP,Np,x)
    DO i = 1 , Np
        Llpt(i) = NINT(x(i,2))
        WRITE (*,*) i , Fluxv(i)*rc**2 , Llpt(i) , Fluxv(Llpt(i))*rc**2
    ENDDO
    END

!____________________________________________________________________________________________________c_______________________________________________________________________
    SUBROUTINE SORT22(Maxp,No,X)
    IMPLICIT NONE

    INTEGER i , jo , ko , Maxp , mo , No
    DOUBLE PRECISION temp , temp2 , X

!$$$$$ calls no other routines
!  In-place rapid sorter.  Rearranges  x  into ascending order of x(i,1)
    DIMENSION X(Maxp,2)
!
    mo = No
100  IF ( mo.LE.15 ) THEN
        IF ( mo.LE.1 ) RETURN
        mo = 2*(mo/4) + 1
    ELSE
        mo = 2*(mo/8) + 1
    ENDIF
    ko = No - mo
    jo = 1
    i = jo
200  IF ( X(i,1).GT.X(i+mo,1) ) THEN
        temp = X(i,1)
        temp2 = X(i,2)
        X(i,1) = X(i+mo,1)
        X(i,2) = X(i+mo,2)
        X(i+mo,1) = temp
        X(i+mo,2) = temp2
        i = i - mo
        IF ( i.GE.1 ) GOTO 200
    ENDIF
    jo = 1 + jo
    IF ( jo.GT.ko ) GOTO 100
    i = jo
    GOTO 200
    END

!______________________________________________________________________
    SUBROUTINE DORITO(N,P,Ntria,Iv,V,Istack)
    IMPLICIT NONE

    DOUBLE PRECISION AMAx , AMIn , AREA , darea , DMAx , DMIn ,       &
                    & dotmax , DOTTY , fpi , P , p12 , p13 , p23 ,     &
                    & prism , rsq , sarea , V , x
    INTEGER i , i1 , i2 , i3 , id , INP , IOUt , IPRint , iprism ,    &
        & isp , Istack , itemp , Iv , j , jt , k , km , kmt , kt ,  &
        & kv
    INTEGER l1 , l2 , N , n1 , Ntria , nuc

!$$$$$ calls scirc, dotty, area
!  Routine to assign triangles to a set of  n  3-vectors  p,on the
!  surface of the unit sphere based upon
!  the method of Watson, Computers & Geosciences, v8,1, pp 97-101, 1982.
!  the final network is a set of Delaunay triangles with the property
!  that no vertex lies inside the circumcircle of any triangle.  This
!  system is as close to equiangular as possible.
!  n   is the number of input points.
!  p   is an array of data, dimensioned  p(4, n), where
!      the 4-th component in each vector is ignored.  It is
!      present to allow field data in the fourth element.
!  ntria on exit contains the total number of triangles generated.
!  iv  on exit contains the vertices of the tessellation -
!      iv(1,k),iv(2,k),iv(3,k)  defines the corners of the  k-th
!      triangle as the points  (p(1,iv(1,k), p(2,iv(1,k)),p(3,iv(1,k)),
!      (p(1,iv(2,k), p(2,iv(2,k)), p(3,iv(3,k)),
!      (p(1,iv(3,k)), p(2,iv(3,k),p(3,iv(3,k)).
!      iv  must be dimensioned at least as iv(3, 2*n).
!  v(i,k), i=1,,4 has coordinates of circumcenter & "radius" of circle
!      circumscribing triangle k in list iv.
!  istack  another work array of length 2*n
!
!
    DIMENSION P(4,*) , V(4,*) , Iv(3,*) , Istack(*) , prism(3,6) ,    &
            & iprism(6)
    DIMENSION itemp(3,2) , kv(2,500) , dotmax(6)
    COMMON /BOUND / DMIn , DMAx , AMIn , AMAx
    COMMON /IO    / INP , IOUt , IPRint
    DATA ((prism(i,j),i=1,3),j=1,6)/0 , 0 , 1 , 1 , 0 , 0 , 0 , 1 ,   &
        & 0 , -1 , 0 , 0 , 0 , -1 , 0 , 0 , 0 , -1/
!
!
    DATA itemp/1 , 1 , 2 , 2 , 3 , 3/ , fpi/12.56637062/
!
!  Find closest points in array p to those in prism and use these as a
!  starting configuration
!
!  Initialize dotmax
    DO j = 1 , 6
        dotmax(j) = 0.0
    ENDDO
!  Compare every point to prism points, keep the index of those closest
!
    DO j = 1 , N
        DO i = 1 , 6
        x = DOTTY(P(1,j),prism(1,i))
        IF ( x.GT.dotmax(i) ) THEN
            dotmax(i) = x
            iprism(i) = j
        ENDIF
        ENDDO
    ENDDO
!  Now make a tesselation from the points in iprism and save in iv(i,j)
!
    Iv(1,1) = iprism(1)
    Iv(2,1) = iprism(2)
    Iv(3,1) = iprism(3)
    Iv(1,2) = iprism(1)
    Iv(2,2) = iprism(3)
    Iv(3,2) = iprism(4)
    Iv(1,3) = iprism(1)
    Iv(2,3) = iprism(4)
    Iv(3,3) = iprism(5)
    Iv(1,4) = iprism(1)
    Iv(2,4) = iprism(2)
    Iv(3,4) = iprism(5)
    Iv(1,5) = iprism(2)
    Iv(2,5) = iprism(3)
    Iv(3,5) = iprism(6)
    Iv(1,6) = iprism(3)
    Iv(2,6) = iprism(4)
    Iv(3,6) = iprism(6)
    Iv(1,7) = iprism(4)
    Iv(2,7) = iprism(5)
    Iv(3,7) = iprism(6)
    Iv(1,8) = iprism(2)
    Iv(2,8) = iprism(5)
    Iv(3,8) = iprism(6)
!  put in starting pointers
!
    n1 = 2*N
    DO i = 1 , n1
        Istack(i) = i
    ENDDO
!
!  Starting configuration for the sphere is specified by
!  8 right triangles, in iv(i,j),j=1,...8
!  Compute the circumcentres and "radii" for these triangles.


    DO i = 1 , 8
        CALL SCIRC(P(1,Iv(1,i)),P(1,Iv(2,i)),P(1,Iv(3,i)),V(1,i))
    ENDDO
    isp = 8
    id = 9
    km = 1
!  Scan through the data forwards.
    DO nuc = 1 , N
!  Skip this point if it is part of the starting octahedron.
        DO i = 1 , 6
        IF ( nuc.EQ.iprism(i) ) GOTO 100
        ENDDO
        km = 0
!  Loop through the established 3-tuples.
!
        DO jt = 1 , isp
!  Test if the new data point is within the  jt  circumcircle.
        rsq = DOTTY(P(1,nuc),V(1,jt))
        IF ( rsq.GE.V(4,jt) ) THEN
!  The point is within.  Delete this 3-tuple but save its edges.
            id = id - 1
            Istack(id) = jt
!  Add edges to  kv  but delete if already present.
            DO i = 1 , 3
                l1 = itemp(i,1)
                l2 = itemp(i,2)
                IF ( km.GT.0 ) THEN
                    kmt = km
                    DO j = 1 , kmt
                    IF ( Iv(l1,jt).EQ.kv(1,j) ) THEN
                        IF ( Iv(l2,jt).EQ.kv(2,j) ) THEN
                            km = km - 1
                            IF ( j.LE.km ) THEN
                                DO k = j , km
                                kv(1,k) = kv(1,k+1)
                                kv(2,k) = kv(2,k+1)
                                ENDDO
                            ENDIF
                            GOTO 10
                        ENDIF
                    ENDIF
                    ENDDO
                ENDIF
                km = 1 + km
                kv(1,km) = Iv(l1,jt)
                kv(2,km) = Iv(l2,jt)
10            ENDDO
        ENDIF
        ENDDO
!
!  Form new 3-tuples.
        DO i = 1 , km
        kt = Istack(id)
        id = 1 + id
!  Calculate the circumcircle center and radius squared.
        i1 = kv(1,i)
        i2 = kv(2,i)
        CALL SCIRC(P(1,i1),P(1,i2),P(1,nuc),V(1,kt))
        Iv(1,kt) = kv(1,i)
        Iv(2,kt) = kv(2,i)
        Iv(3,kt) = nuc
        ENDDO
        isp = 2 + isp
100  ENDDO
!
!  check no. of triangles is right and surface area integrates to 4pi.
!
    Ntria = isp
    isp = 2*N - 4
    IF ( Ntria.NE.isp ) WRITE (IOUt,*) 'only ' , Ntria ,              &
                            &'triangles, should be ' , isp
    sarea = 0.0
    AMIn = fpi
    AMAx = 0.0
    DO i = 1 , Ntria
        i1 = Iv(1,i)
        i2 = Iv(2,i)
        i3 = Iv(3,i)
        x = AREA(P(1,i1),P(1,i2),P(1,i3))
        sarea = sarea + x
        IF ( x.LT.AMIn ) AMIn = x
        IF ( x.GT.AMAx ) AMAx = x
        p12 = DOTTY(P(1,i1),P(1,i2))
        DMAx = MAX(p12,DMAx)
        DMIn = MIN(p12,DMIn)
        p13 = DOTTY(P(1,i1),P(1,i3))
        DMAx = MAX(p13,DMAx)
        DMIn = MIN(DMIn,p13)
        p23 = DOTTY(P(1,i2),P(1,i3))
        DMAx = MAX(p23,DMAx)
        DMIn = MIN(DMIn,p23)
    ENDDO
    darea = 100.*sarea/fpi
    WRITE (IOUt,*) 'Surface area of tesselation is ' , sarea
    WRITE (IOUt,*) 'This is ' , darea , '% that of sphere'
    WRITE (IOUt,*) 'Minimum triangle area is ' , AMIn
    WRITE (IOUt,*) 'Maximum triangle area is ' , AMAx
    END

!
!
!_______________________________________________________________________
    SUBROUTINE SCIRC(X1,X2,X3,C)
    IMPLICIT NONE

    DOUBLE PRECISION a , b , C , cnorm , DOTTY , resq , X1 , X2 , X3
    INTEGER i , INP , IOUt , IPRint

!  for points x1,x2,x3,which are the vertices of a spherical triangle
!  calculates the circum circle origin c(1),c(2),c(3), and the dot
!  product of the centre with any vertex, c(4),
    DIMENSION X1(3) , X2(3) , X3(3) , C(4) , a(3,3) , b(3)
    COMMON /IO    / INP , IOUt , IPRint
    DO i = 1 , 3
        a(1,i) = X1(i)
        a(2,i) = X2(i)
        a(3,i) = X3(i)
        b(i) = 1.0
    ENDDO
!    write(6,*) resq,' resqin'
    CALL QR(3,3,3,a,b,C,resq)
!	write(6,*) resq,' resqout'
    IF ( resq.LT.0. ) THEN
        WRITE (IOUt,*) 'problem in qr'
        WRITE (IOUt,*) (X1(i),i=1,3)
        WRITE (IOUt,*) (X2(i),i=1,3)
        WRITE (IOUt,*) (X3(i),i=1,3)
        STOP
    ENDIF
!  normalize c
    cnorm = SQRT(DOTTY(C,C))
    DO i = 1 , 3
        C(i) = C(i)/cnorm
    ENDDO
!  compute "radius"
    C(4) = DOTTY(C,X1)
    END

!_______________________________________________________________________
    FUNCTION DOTTY(X,Y)
    IMPLICIT NONE

    DOUBLE PRECISION DOTTY , X , Y

!  computes the dot product of x and y
!
    DIMENSION X(3) , Y(3)
    DOTTY = X(1)*Y(1) + X(2)*Y(2) + X(3)*Y(3)
    END

!
!_______________________________________________________________________
    FUNCTION AREA(A,B,C)
    IMPLICIT NONE

    DOUBLE PRECISION A , aa , AREA , B , bb , C , cc , cosa , cosb ,  &
                    & cosc , sina , sinb , sinc
    INTEGER INP , IOUt , IPRint

!$$$$ calls no other routines
!  Given three unit vectors a, b, c  finds the area of the
!  SPHERICAL triangle formed at the tips.
    COMMON /IO    / INP , IOUt , IPRint
    DIMENSION A(3) , B(3) , C(3)
!
    cosa = B(1)*C(1) + B(2)*C(2) + B(3)*C(3)
    cosb = A(1)*C(1) + A(2)*C(2) + A(3)*C(3)
    cosc = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
    sina = SQRT(ABS(1.D0-cosa**2))
    sinb = SQRT(ABS(1.D0-cosb**2))
    sinc = SQRT(ABS(1.D0-cosc**2))
    aa = ACOS((cosa-cosb*cosc)/(sinb*sinc))
    bb = ACOS((cosb-cosa*cosc)/(sina*sinc))
    cc = ACOS((cosc-cosb*cosa)/(sinb*sina))
    AREA = aa + bb + cc - 3.14159265358979324
!      write(iout, '(9f8.4/a,f8.4)')a,b,c,' Area =',area
    END

    SUBROUTINE QR(Ndime,M,N,A,B,X,Resq)
    IMPLICIT NONE

    DOUBLE PRECISION A , B , const , qv1 , Resq , sq , u1 , X
    INTEGER i , i1 , ii , j , j1 , jj , M , N , Ndime

!$$$$  calls no other routines
!  solves over-determined least-squares problem  ax = b
!  where  a  is an  m by n  matrix,  b  is an m-vector .
!  resq  is the sum of squared residuals of optimal solution.  also used
!  to signal error conditions - if -2 , system is underdetermined,  if
!  -1,  system is singular.
!  method - successive householder rotations.  see lawson+hanson - solv
!  -ing least squares problems.
!  routine will also work when m=n.
!*****   caution -  a and b  are overwritten by this routine.
    DIMENSION A(Ndime,N+1) , B(Ndime) , X(N+1)
    DOUBLE PRECISION sum , dot
!
    Resq = -2.0D0
    IF ( M.LT.N ) RETURN
!   loop ending on 1800 rotates  a  into upper triangular form
    DO j = 1 , N
!  find constants for rotation and diagonal entry
        sq = 0.00
        DO i = j , M
        sq = A(i,j)**2 + sq
        ENDDO
        qv1 = -SIGN(SQRT(sq),A(j,j))
        u1 = A(j,j) - qv1
        A(j,j) = qv1
        j1 = j + 1
        IF ( j1.LE.N ) THEN
!  rotate remaining columns of sub-matrix
        DO jj = j1 , N
            dot = u1*A(j,jj)
            DO i = j1 , M
                dot = A(i,jj)*A(i,j) + dot
            ENDDO
            const = dot/ABS(qv1*u1)
            DO i = j1 , M
                A(i,jj) = A(i,jj) - const*A(i,j)
            ENDDO
            A(j,jj) = A(j,jj) - const*u1
        ENDDO
        ENDIF
!  rotate  b  vector
        dot = u1*B(j)
        IF ( j1.LE.M ) THEN
        DO i = j1 , M
            dot = B(i)*A(i,j) + dot
        ENDDO
        ENDIF
        const = dot/ABS(qv1*u1)
        B(j) = B(j) - const*u1
        IF ( j1.LE.M ) THEN
        DO i = j1 , M
            B(i) = B(i) - const*A(i,j)
        ENDDO
        ENDIF
    ENDDO
!  solve triangular system by back-substitution.
    Resq = -1.00
    DO ii = 1 , N
        i = N - ii + 1
        sum = B(i)
        IF ( ii.NE.1 ) THEN
        i1 = i + 1
        DO j = i1 , N
            sum = sum - A(i,j)*X(j)
        ENDDO
        ENDIF
        IF ( A(i,i).EQ.0.0D0 ) RETURN
        X(i) = sum/A(i,i)
    ENDDO
!  find residual in overdetermined case.
    Resq = 0.0D0
    IF ( M.EQ.N ) RETURN
    i1 = N + 1
    DO i = i1 , M
        Resq = B(i)**2 + Resq
    ENDDO
    END

!
!
    SUBROUTINE SSPQR(Mdim,M,N,A,B,X,Resq,Kk)
    IMPLICIT NONE

    DOUBLE PRECISION A , B , const , qv1 , Resq , sq , u1 , X
    INTEGER i , i1 , ii , j , j1 , jj , Kk , l , M , Mdim , N

!     calls no other routines
!     solves overdetermined least squares problem ax=b
!     for special matrices a, where a has dimension 2*n+kk by n, and
!     a(i,j)=0., if i.gt.2*j+kk.  this special form appears in the depleted
!     basis analysis of harmonic splines with kk=0, and elsewhere.
!     resq is the sum of squared residuals of optimal solution.  also
!     used to signal error conditions-  if -2 the system is underdetermined, if
!     -1, system is singular.
!     method - successive householder rotations. see lawson & hanson,
!     solving least squares problems.
!**********caution - a and b are overwritten by this routine********************

    DIMENSION A(Mdim,N) , B(Mdim) , X(N)
    DOUBLE PRECISION sum , dot
!
    Resq = -2.0
!      if(m.ne.2*n+kk)return
!    loop ending on 1800 rotates a into upper triangular form.
    DO j = 1 , N
!  find constants for rotation and diagonal entry
        sq = 0.0D0
!     dot products and other action computed only down to l=2*j+kk in each column
        l = 2*j + Kk
        DO i = j , l
        sq = A(i,j)**2 + sq
        ENDDO
        qv1 = -SIGN(SQRT(sq),A(j,j))
        u1 = A(j,j) - qv1
        A(j,j) = qv1
        j1 = j + 1
        IF ( j1.LE.N ) THEN
!  rotate remaining columns of sub-matrix
        DO jj = j1 , N
            dot = u1*A(j,jj)
            DO i = j1 , l
                dot = A(i,jj)*A(i,j) + dot
            ENDDO
            const = dot/ABS(qv1*u1)
            DO i = j1 , l
                A(i,jj) = A(i,jj) - const*A(i,j)
            ENDDO
            A(j,jj) = A(j,jj) - const*u1
        ENDDO
        ENDIF
!  rotate  b  vector
        dot = u1*B(j)
        IF ( j1.LE.l ) THEN
        DO i = j1 , l
            dot = B(i)*A(i,j) + dot
        ENDDO
        ENDIF
        const = dot/ABS(qv1*u1)
        B(j) = B(j) - const*u1
        IF ( j1.LE.l ) THEN
        DO i = j1 , l
            B(i) = B(i) - const*A(i,j)
        ENDDO
        ENDIF
    ENDDO
!  solve triangular system by back-substitution.
    Resq = -1.0D0
    DO ii = 1 , N
        i = N - ii + 1
        sum = B(i)
        IF ( ii.NE.1 ) THEN
        i1 = i + 1
        DO j = i1 , N
            sum = sum - A(i,j)*X(j)
        ENDDO
        ENDIF
        IF ( A(i,i).EQ.0.0D0 ) RETURN
        X(i) = sum/A(i,i)
    ENDDO
!  find residual in overdetermined case.
    Resq = 0.0D0
    IF ( M.EQ.N ) RETURN
    i1 = N + 1
    DO i = i1 , M
        Resq = B(i)**2 + Resq
    ENDDO
    END
!

!_____________________________________________________________________
    SUBROUTINE SBODY(Kdim,Maxmum,N,P,Iv,Ntria)
    IMPLICIT NONE

    DOUBLE PRECISION AMAx , AMIn , cnorm , DMAx , DMIn , DOTTY , P ,  &
                    & V , val
    double precision :: bColat1 , bColat2 , bLong1 , bLong2, r0, t
    INTEGER i , INP , IOUt , IPRint , ISTack , Iv , j , k , Kdim ,    &
        & Maxmum , maxv , N , new , nfound , niv , Ntria , nwant1 , &
        & NXYZ , NXYZ2

!$$$$ calls dorito,synops
!  If kdim=3, gets a series of x,y,z coords that define points on the
!  surface of a sphere; if kdim=4 x,y,z coordinates plus associated
!  field value are read in.
!  Note that the position vector is dimensioned p(4,*) IN EITHER CASE.
!  Computes the tessellation or reads the appropriate
!  pointers from a file.
!
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ)
    CHARACTER*80 name
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /BOUND / DMIn , DMAx , AMIn , AMAx
    COMMON /CIRCLE/ V(4,NXYZ2) , ISTack(NXYZ2)
    DIMENSION P(4,*), Iv(*)
    DIMENSION val(10)
!
!
    CALL GETVAL('tnew',val,nwant1,nfound)
    new = NINT(val(1))
!
!  Get x, y, z, triples from disk.
    IF ( Kdim.EQ.3 ) CALL GETCHR('corepoints',name,nfound)
    IF ( Kdim.EQ.4 ) CALL GETCHR('corefield',name,nfound)
    WRITE (IOUt,'(1x,a,a50)') 'Core points read from file: ' , name
    OPEN (UNIT=12,FILE=name)
    READ (12, *) bColat1 , bColat2 , bLong1 , bLong2
    READ (12, *) r0
!   Read lat and lon into the array sized for xyz coordinates
!   call lltoxyz to reassign the values accordingly
!   so p(3, i) is just a dummy variable when calling the function 
    DO i = 1 , NXYZ
        READ (12,*,ERR=100,END=200) (p(j, i), j=1, Kdim), t
        Call LLTOXYZ(p(1:3, i))
    ENDDO
100  WRITE (IOUt,'(1x,a,i6,a/a)') 'Core point file was truncated at ' ,&
                                & NXYZ , ' data' ,                     &
                            &' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
200  N = i - 1
!  check unit vectors
    DO i = 1 , N
        cnorm = SQRT(DOTTY(P(1,i),P(1,i)))
        P(1,i) = P(1,i)/cnorm
        P(2,i) = P(2,i)/cnorm
        P(3,i) = P(3,i)/cnorm
    ENDDO
    CLOSE (UNIT=12)
    WRITE (IOUt,'(1x,a,i6)') 'Number of coordinates read =' , N
!
!
!  Pointers on disk.  Get pointer vectors from disk.
    IF ( new.EQ.0 ) THEN
        CALL GETCHR('pointer',name,nfound)
        WRITE (*,'(a)') name
        OPEN (UNIT=13,FILE=name)
        maxv = 2*NXYZ - 4
        Ntria = 0
        niv = 1
        DO k = 1 , maxv
        READ (13,*,ERR=250,END=300) Iv(niv) , Iv(niv+1) , Iv(niv+2)
        niv = niv + 3
        Ntria = Ntria + 1
        ENDDO
250     WRITE (IOUt,'(//(1x,a/))')                                     &
                    &'Pointer file incompatible with core point data'&
                    & , '*******************************************'
        STOP
300     WRITE (IOUt,*) 'ntria = ' , Ntria
        IF ( Ntria.NE.2*N-4 ) GOTO 250
        WRITE (IOUt,'(1x,a,a50)') 'Pointers read from file: ' , name
        CLOSE (UNIT=13)
!
!  No pointers on disk.  Therefore generate pointer series with
!  Delaunay algorithm. First perform the calculations.
    ELSE
!
        CALL DORITO(N,P,Ntria,Iv,V,ISTack)
        WRITE (IOUt,'(1x,a)') 'Delaunay phase complete'
!
!  Save pointers on the named file.
        CALL GETCHR('pointer',name,nfound)
        OPEN (UNIT=14,FILE=name)
        WRITE (14,'(3I7)') (Iv(j),j=1,3*Ntria)
        CLOSE (UNIT=14)
        WRITE (IOUt,'(1x,a)') 'Pointers written to disk'
    ENDIF
    CALL SYNOPS(new,N,P,Ntria,Iv)
    END

!_______________________________________________________________________
    SUBROUTINE NNEIGH
    IMPLICIT NONE

    DOUBLE PRECISION AMAx , AMIn , cov , DMAx , DMIn , P , SAREA , sum
    INTEGER i , INP , IOUt , IPRint , IV , ivj , j , jv , k , N ,     &
        & nay , NEAr , neigh , nn , NTRia , NVERT , NXYZ , NXYZ2

!  From an arbitrarily ordered input list of unit vectors, the
!  program creates for each one a list of neighbors
!  such that they form a convex polygon around that vertex and the
!  polygon encloses no other vertex.  These polygons taken
!  together cover the sphere exactly three times.  The polygon
!  forms the support region for an interpolation function on the sphere
!  centered on the associated vertex.
!
!$$$$ Requires: isort, order
!
!  The parameter nxyz is the maximum number of unit vectors allowed.
!  The parameter nvert is the maximum number of faces at any vertex.
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ,NVERT=10)
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    DIMENSION neigh(100)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /BOUND / DMIn , DMAx , AMIn , AMAx
!
    DMIn = 1.0
    DMAx = -1.0
    WRITE (IOUt,*) ' Entering nneigh'
!
!  For all unit vectors discover the neighbor set.
!
    DO j = 1 , N
        NEAr(1,j) = 0
!  Scan through the triangle list, seeking the vertex  j.
        nay = 0
        DO jv = 1 , NTRia
        IF ( IV(1,jv).EQ.j .OR. IV(2,jv).EQ.j .OR. IV(3,jv).EQ.j )  &
            & THEN
            DO i = 1 , 3
                ivj = IV(i,jv)
                IF ( ivj.NE.j .AND. ivj.LE.N ) THEN
                    nay = 1 + nay
                    neigh(nay) = ivj
                ENDIF
            ENDDO
        ENDIF
        ENDDO
!  Sort the list, then discard duplicate vertices.
        CALL ISORT(nay,neigh)
        nn = 1
        DO i = 2 , nay
        IF ( neigh(i).NE.neigh(i-1) ) THEN
            nn = nn + 1
            IF ( nn.GT.NVERT ) WRITE (IOUt,'(a,i4,a,i4)')            &
                & ' WARNING: More than ' , NVERT ,                   &
                &' neighbors at vertex ' , j
            neigh(nn) = neigh(i)
        ENDIF
        ENDDO
!  Rank the list by longitude about the point p(.,j), so that vertices
!  progress sequentially about centre. neigh(1) is assumed zero
!  longitude
        CALL ORDER(j,nn,neigh)
!
        IF ( nn.LT.NVERT ) NEAr(nn+1,j) = 0
!        write(iout, '(i5,(4x,14i5/))')j,(near(i,j),i=1,nn)
!
!
    ENDDO
!  check completeness of covering polygons by computing total area.
!  should be 12pi
!
    sum = 0.0
    DO j = 1 , N
        DO i = 2 , NVERT
        IF ( NEAr(i,j).EQ.0 ) GOTO 50
        sum = sum + SAREA(P(1,j),P(1,NEAr(i-1,j)),P(1,NEAr(i,j)))
        k = i
        ENDDO
50      sum = sum + SAREA(P(1,j),P(1,NEAr(1,j)),P(1,NEAr(k,j)))
    ENDDO
    cov = sum/(12*3.1415926)
    WRITE (IOUt,'(/a,f7.4)') 'Covering area /12*pi =' , cov
    IF ( ABS(cov-1.0).GT.0.02 ) WRITE (IOUt,*)                        &
        & 'covering is inconsistent'
!
    END

!_______________________________________________________________________
    SUBROUTINE ISORT(No,X)
    IMPLICIT NONE

    INTEGER i , jo , ko , mo , No
    DOUBLE PRECISION temp

!$$$$$ calls no other routines
!  In-place rapid sorter.  Rearranges  x  into ascending order
    INTEGER X(No)
!
    mo = No
100  IF ( mo.LE.15 ) THEN
        IF ( mo.LE.1 ) RETURN
        mo = 2*(mo/4) + 1
    ELSE
        mo = 2*(mo/8) + 1
    ENDIF
    ko = No - mo
    jo = 1
    i = jo
200  IF ( X(i).GT.X(i+mo) ) THEN
        temp = X(i)
        X(i) = X(i+mo)
        X(i+mo) = temp
        i = i - mo
        IF ( i.GE.1 ) GOTO 200
    ENDIF
    jo = 1 + jo
    IF ( jo.GT.ko ) GOTO 100
    i = jo
    GOTO 200
    END

!
!_______________________________________________________________________
    SUBROUTINE ORDER(J,Nn,Neigh)
    IMPLICIT NONE

    INTEGER i , IV , J , N , NEAr , Neigh , Nn , NTRia , NVERT ,      &
        & NXYZ , NXYZ2
    DOUBLE PRECISION P , phi , tpi , x

!
!  orders the nn nearest neighbours of point j, ranking by longitude
!  from neigh(nn)
!
!*****requires sort2, gpole
!
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ,NVERT=10)
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    DIMENSION x(3) , phi(NVERT,2) , Neigh(Nn)
    DATA tpi/6.28318531/
!
!  first rotate so that p(.,j) is the pole, zero longitude is given by
!  neigh(1) and compute longitude phi for each neighbour.
    DO i = 1 , Nn
        CALL GPOLE(1,P(1,J),P(1,Neigh(1)),P(1,Neigh(i)),x)
        phi(i,2) = FLOAT(i)
        phi(i,1) = ATAN2(x(2),x(1))
        phi(i,1) = DMOD(phi(i,1)+tpi,tpi)
    ENDDO
!  rank neigh(i) by phi(i,1)
    CALL SORT2(NVERT,Nn,phi)
    DO i = 1 , Nn
        NEAr(i,J) = Neigh(NINT(phi(i,2)))
    ENDDO
    END

!
!_______________________________________________________________________
    SUBROUTINE SORT2(Nd,No,X)
    IMPLICIT NONE

    INTEGER i , jo , ko , mo , Nd , No
    DOUBLE PRECISION temp , temp2 , X

!$$$$$ calls no other routines
!  In-place rapid sorter.  Rearranges  x  into ascending order of first
!  column
    DIMENSION X(Nd,2)
!
    mo = No
100  IF ( mo.LE.15 ) THEN
        IF ( mo.LE.1 ) RETURN
        mo = 2*(mo/4) + 1
    ELSE
        mo = 2*(mo/8) + 1
    ENDIF
    ko = No - mo
    jo = 1
    i = jo
200  IF ( X(i,1).GT.X(i+mo,1) ) THEN
        temp = X(i,1)
        temp2 = X(i,2)
        X(i,1) = X(i+mo,1)
        X(i,2) = X(i+mo,2)
        X(i+mo,1) = temp
        X(i+mo,2) = temp2
        i = i - mo
        IF ( i.GE.1 ) GOTO 200
    ENDIF
    jo = 1 + jo
    IF ( jo.GT.ko ) GOTO 100
    i = jo
    GOTO 200
    END

!
!_______________________________________________________________________
    FUNCTION SAREA(A,B,C)
    IMPLICIT NONE

    DOUBLE PRECISION A , aa , B , bb , C , cc , cosa , cosb , cosc ,  &
                    & SAREA , SGN , sina , sinb , sinc
    INTEGER INP , IOUt , IPRint

!$$$$ calls sgn
!  Given three unit vectors a, b, c  finds the area of the
!  SPHERICAL triangle formed at the tips.
    DIMENSION A(3) , B(3) , C(3)
    COMMON /IO    / INP , IOUt , IPRint
!
!	write(iout,*)'Entering sarea'
    cosa = B(1)*C(1) + B(2)*C(2) + B(3)*C(3)
    cosb = A(1)*C(1) + A(2)*C(2) + A(3)*C(3)
    cosc = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
    sina = SQRT(ABS(1.0-cosa**2))
    sinb = SQRT(ABS(1.0-cosb**2))
    sinc = SQRT(ABS(1.0-cosc**2))
    aa = ((cosa-cosb*cosc)/(sinb*sinc))
    IF ( ABS(aa).GT.1.0 ) aa = SGN(aa)
    bb = ((cosb-cosa*cosc)/(sina*sinc))
    IF ( ABS(bb).GT.1.0 ) bb = SGN(bb)
    cc = ((cosc-cosb*cosa)/(sinb*sina))
    IF ( ABS(cc).GT.1.0 ) cc = SGN(cc)
    aa = ACOS(aa)
    bb = ACOS(bb)
    cc = ACOS(cc)
    SAREA = aa + bb + cc - 3.14159265358979324
!	write(iout, '(9f8.4/a,f8.4)')a,b,c,' Area =',sarea
    END

!
!_______________________________________________________________________
    FUNCTION SGN(X)
    IMPLICIT NONE

    DOUBLE PRECISION SGN , X

    IF ( X.LE.0. ) THEN
        SGN = -1.0
    ELSE
        SGN = 1.0
    ENDIF
    END

!_________________________________________________________________________
    SUBROUTINE TRIQ(Nx,Ny)
    IMPLICIT NONE

    INTEGER i , INP , IOUt , IPRint , ixp , iyp , j , m , mx , mx1 ,  &
        & mx2 , my , my1 , my2 , NQUad , ntwo , Nx , Ny
    DOUBLE PRECISION xw , XYW , yw , zw

!$$$$ calls no other routines
!  Constructs the conical product quadrature formula for integration
!  over the plane triangular region  0.le.x, 0.le.y, x+y.le.1. see
!  Stroud,1971, 'Approximate Calculation of Multiple Integrals',pp28-30.
!  nx  is the number of sample points desired along  x  lines
!  ny  the number along  y  lines
!  nquad   the number of sample points generated.  ideally this is
!      nx*ny, but may differ from this because of the limited
!      number of gauss formulas stored.
! xyw  triples of  x-y coordintes and weights generated for the
!      quadrature formula.  if n=min(nx,ny), the formula is
!      exact for polynomials in x**p*y**q  where p+q.le.2*n-1.
!  if  nx  is negative the routine supplies the special 7-point degree 5
!  formula of Stroud (*T2 5-1, page 315).
!
!  The results are loaded into the common /rule/
!
    COMMON /RULE  / XYW(3,100) , NQUad
    COMMON /IO    / INP , IOUt , IPRint
    DIMENSION ixp(15) , xw(316) , iyp(9) , yw(72) , zw(21)
!
!  Gauss formulas with w(x)=1, with orders 2,3,...10,12,16,20,24,32.
    DATA ixp/1 , 5 , 11 , 19 , 29 , 41 , 55 , 71 , 89 , 109 , 133 ,   &
        & 165 , 205 , 253 , 317/
    DATA (xw(j),j=1,54)/0.2113248654 , 0.5000000000 , 0.7886751346 ,  &
        & 0.5000000000 , 0.1127016654 , 0.2777777778 , 0.5000000000 , &
        & 0.4444444444 , 0.8872983346 , 0.2777777778 , 0.0694318442 , &
        & 0.1739274225 , 0.3300094782 , 0.3260725774 , 0.6699905218 , &
        & 0.3260725774 , 0.9305681558 , 0.1739274225 , 0.0469100771 , &
        & 0.1184634425 , 0.2307653450 , 0.2393143352 , 0.5000000000 , &
        & 0.2844444444 , 0.7692346550 , 0.2393143352 , 0.9530899229 , &
        & 0.1184634425 , 0.0337652429 , 0.0856622462 , 0.1693953068 , &
        & 0.1803807865 , 0.3806904070 , 0.2339569673 , 0.6193095930 , &
        & 0.2339569673 , 0.8306046932 , 0.1803807865 , 0.9662347571 , &
        & 0.0856622462 , 0.0254460439 , 0.0647424831 , 0.1292344072 , &
        & 0.1398526957 , 0.2970774243 , 0.1909150252 , 0.5000000000 , &
        & 0.2089795918 , 0.7029225757 , 0.1909150252 , 0.8707655928 , &
        & 0.1398526957 , 0.9745539561 , 0.0647424831/
    DATA (xw(j),j=55,108)/0.0198550718 , 0.0506142681 , 0.1016667613 ,&
        & 0.1111905172 , 0.2372337951 , 0.1568533229 , 0.4082826788 , &
        & 0.1813418917 , 0.5917173212 , 0.1813418917 , 0.7627662049 , &
        & 0.1568533229 , 0.8983332387 , 0.1111905172 , 0.9801449282 , &
        & 0.0506142681 , 0.0159198803 , 0.0406371942 , 0.0819844464 , &
        & 0.0903240803 , 0.1933142837 , 0.1303053482 , 0.3378732883 , &
        & 0.1561735385 , 0.5000000000 , 0.1651196775 , 0.6621267117 , &
        & 0.1561735385 , 0.8066857163 , 0.1303053482 , 0.9180155536 , &
        & 0.0903240803 , 0.9840801197 , 0.0406371942 , 0.0130467358 , &
        & 0.0333356722 , 0.0674683167 , 0.0747256746 , 0.1602952159 , &
        & 0.1095431812 , 0.2833023030 , 0.1346333596 , 0.4255628305 , &
        & 0.1477621123 , 0.5744371695 , 0.1477621123 , 0.7166976970 , &
        & 0.1346333596 , 0.8397047841 , 0.1095431812 , 0.9325316833 , &
        & 0.0747256746 , 0.9869532642 , 0.0333356722/
    DATA (xw(j),j=109,164)/0.0092196830 , 0.0235876680 ,              &
        & 0.0479413720 , 0.0534696630 , 0.1150486630 , 0.0800391645 , &
        & 0.2063410230 , 0.1015837135 , 0.3160842505 , 0.1167462685 , &
        & 0.4373832955 , 0.1245735230 , 0.5626167045 , 0.1245735230 , &
        & 0.6839157495 , 0.1167462685 , 0.7936589770 , 0.1015837135 , &
        & 0.8849513370 , 0.0800391645 , 0.9520586280 , 0.0534696630 , &
        & 0.9907803170 , 0.0235876680 , 0.0052995328 , 0.0135762297 , &
        & 0.0277124885 , 0.0311267620 , 0.0671843988 , 0.0475792558 , &
        & 0.1222977958 , 0.0623144856 , 0.1910618778 , 0.0747979944 , &
        & 0.2709916112 , 0.0845782597 , 0.3591982246 , 0.0913017075 , &
        & 0.4524937451 , 0.0947253052 , 0.5475062549 , 0.0947253052 , &
        & 0.6408017754 , 0.0913017075 , 0.7290083888 , 0.0845782597 , &
        & 0.8089381222 , 0.0747979944 , 0.8777022042 , 0.0623144856 , &
        & 0.9328156012 , 0.0475792558 , 0.9722875115 , 0.0311267620 , &
        & 0.9947004672 , 0.0135762297/
    DATA (xw(j),j=165,240)/0.0034357004 , 0.0088070035 ,              &
        & 0.0180140364 , 0.0203007149 , 0.0438827859 , 0.0313360241 , &
        & 0.0804415141 , 0.0416383708 , 0.1268340468 , 0.0509650599 , &
        & 0.1819731597 , 0.0590972660 , 0.2445664990 , 0.0658443192 , &
        & 0.3131469557 , 0.0710480546 , 0.3861070745 , 0.0745864933 , &
        & 0.4617367395 , 0.0763766935 , 0.5382632605 , 0.0763766935 , &
        & 0.6138929255 , 0.0745864933 , 0.6868530443 , 0.0710480546 , &
        & 0.7554335010 , 0.0658443192 , 0.8180268403 , 0.0590972660 , &
        & 0.8731659532 , 0.0509650599 , 0.9195584859 , 0.0416383708 , &
        & 0.9561172141 , 0.0313360241 , 0.9819859636 , 0.0203007149 , &
        & 0.9965642996 , 0.0088070035 , 0.0024063900 , 0.0061706150 , &
        & 0.0126357220 , 0.0142656945 , 0.0308627240 , 0.0221387195 , &
        & 0.0567922365 , 0.0296492925 , 0.0899990070 , 0.0366732405 , &
        & 0.1299379040 , 0.0430950807 , 0.1759531740 , 0.0488093260 , &
        & 0.2272892645 , 0.0537221350 , 0.2831032460 , 0.0577528340 , &
        & 0.3424786600 , 0.0608352365 , 0.4044405663 , 0.0629187280 , &
        & 0.4679715535 , 0.0639690975 , 0.5320284465 , 0.0639690975 , &
        & 0.5955594337 , 0.0629187280 , 0.6575213400 , 0.0608352365 , &
        & 0.7168967540 , 0.0577528340 , 0.7727107355 , 0.0537221350 , &
        & 0.8240468260 , 0.0488093260/
    DATA (xw(j),j=241,316)/0.8700620960 , 0.0430950807 ,              &
        & 0.9100009930 , 0.0366732405 , 0.9432077635 , 0.0296492925 , &
        & 0.9691372760 , 0.0221387195 , 0.9873642780 , 0.0142656945 , &
        & 0.9975936100 , 0.0061706150 , 0.0013680695 , 0.0035093050 , &
        & 0.0071942445 , 0.0081371970 , 0.0176188725 , 0.0126960325 , &
        & 0.0325469625 , 0.0171369310 , 0.0518394225 , 0.0214179490 , &
        & 0.0753161935 , 0.0254990295 , 0.1027581025 , 0.0293420465 , &
        & 0.1339088910 , 0.0329111110 , 0.1684778670 , 0.0361728970 , &
        & 0.2061421215 , 0.0390969475 , 0.2465500460 , 0.0416559620 , &
        & 0.2893243620 , 0.0438260465 , 0.3340656990 , 0.0455869390 , &
        & 0.3803563190 , 0.0469221995 , 0.4277640195 , 0.0478193600 , &
        & 0.4758461675 , 0.0482700440 , 0.5241538325 , 0.0482700440 , &
        & 0.5722359805 , 0.0478193600 , 0.6196436810 , 0.0469221995 , &
        & 0.6659343010 , 0.0455869390 , 0.7106756380 , 0.0438260465 , &
        & 0.7534499540 , 0.0416559620 , 0.7938578785 , 0.0390969475 , &
        & 0.8315221330 , 0.0361728970 , 0.8660911090 , 0.0329111110 , &
        & 0.8972418975 , 0.0293420465 , 0.9246838065 , 0.0254990295 , &
        & 0.9481605775 , 0.0214179490 , 0.9674530375 , 0.0171369310 , &
        & 0.9823811275 , 0.0126960325 , 0.9928057555 , 0.0081371970 , &
        & 0.9986319305 , 0.0035093050/
!
!  Gauss formulas w(x)=x, with orders 1,2 ... 8.
    DATA iyp/1 , 3 , 7 , 13 , 21 , 31 , 43 , 57 , 73/
    DATA yw/0.6666666667 , 0.5 , .3550510257 , .1819586183 ,          &
        & .8449489743 , .3180413817 , .2123405382 , .0698269799 ,      &
        & .5905331356 , .2292411064 , .9114120405 , .2009319137 ,      &
        & .1397598643 , .0311809710 , .4164095676 , .1298475476 ,      &
        & .7231569864 , .2034645680 , .9428958039 , .1355069134 ,      &
        & .0985350858 , .0157479145 , .3045357266 , .0739088701 ,      &
        & .5620251898 , .1463869871 , .8019865821 , .1671746381 ,      &
        & .9601901429 , .0967815902 , .0730543287 , .0087383018 ,      &
        & .2307661380 , .0439551656 , .4413284812 , .0986611509 ,      &
        & .6630153097 , .1407925538 , .8519214003 , .1355424972 ,      &
        & .9706835728 , .0723103307 , .0562625605 , .0052143622 ,      &
        & .1802406917 , .0274083567 , .3526247171 , .0663846965 ,      &
        & .5471536263 , .1071250657 , .7342101772 , .1273908973 ,      &
        & .8853209468 , .1105092582 , .9775206136 , .0559673634 ,      &
        & .0446339553 , .0032951914 , .1443662570 , .0178429027 ,      &
        & .2868247571 , .0454393195 , .4548133152 , .0791995995 ,      &
        & .6280678354 , .1060473594 , .7856915206 , .1125057995 ,      &
        & .9086763921 , .0911190236 , .9822200849 , .0445508044/
!  Special degree 5 formula.
    DATA zw/0.333333333 , 0.333333333 , 0.1125 , 0.101286507 ,        &
        & 0.101286507 , 0.062969590 , 0.797426985 , 0.101286507 ,      &
        & 0.062969590 , 0.101286507 , 0.797426985 , 0.062969590 ,      &
        & 0.470142064 , 0.470142064 , 0.066197076 , 0.059715871 ,      &
        & 0.470142064 , 0.066197076 , 0.470142064 , 0.059715871 ,      &
        & 0.066197076/
!
!
!  Find sample numbers closest to those input
    IF ( Nx.LT.0 ) THEN
!  Copy special formula into output array.
        m = 7
        DO j = 1 , 7
        DO i = 1 , 3
            XYW(i,j) = zw(3*(j-1)+i)
        ENDDO
        ENDDO
        NQUad = m
        IF ( IPRint.GE.3 ) WRITE (IOUt,99001)                          &
                                & ((XYW(i,j),i=1,3),j=1,NQUad)
        GOTO 99999
    ELSE
        ntwo = 2*Nx
        DO mx = 1 , 14
        IF ( ixp(mx+1)-ixp(mx).GE.ntwo ) GOTO 100
        ENDDO
        mx = 14
    ENDIF
100  my = MIN0(Ny,8)
!
    mx1 = ixp(mx)
    mx2 = ixp(mx+1) - 1
    my1 = iyp(my)
    my2 = iyp(my+1) - 1
!
!  Create the product formula.
    m = 0
    DO i = mx1 , mx2 , 2
        DO j = my1 , my2 , 2
        m = 1 + m
        XYW(1,m) = 1.0 - yw(j)
        XYW(2,m) = yw(j)*xw(i)
        XYW(3,m) = yw(j+1)*xw(i+1)
        ENDDO
    ENDDO
    NQUad = m
    IF ( IPRint.GE.3 ) WRITE (IOUt,99001) ((XYW(i,j),i=1,3),j=1,NQUad)
    RETURN
99001 FORMAT (/' Quadrature coordinates and weights'/(3F8.4,3x,3F8.4))
99999 END

!_______________________________________________________________________
!_______________________________________________________________________
    SUBROUTINE GNOMON(Iway,U,R)
    IMPLICIT NONE

    INTEGER INP , IOUt , IPRint , Iway
    DOUBLE PRECISION psi , R , rho , theta , U

!  if iway=1 performs the gnomonic projection of r(x,y,z) on the sphere
!  into u on the plane tangent at the north pole
!  if iway=-1 the inverse transformation is performed
!  All coordinates are cartesian
    COMMON /IO    / INP , IOUt , IPRint
    DIMENSION U(2) , R(3)
    IF ( Iway.EQ.1 ) THEN
        IF ( R(3).EQ.0 ) THEN
        WRITE (IOUt,*) 'gnomonic mapping undefined'
        STOP
        ENDIF
        rho = SQRT(R(1)**2+R(2)**2)/R(3)
        psi = ATAN2(R(2),R(1))
        U(1) = rho*COS(psi)
        U(2) = rho*SIN(psi)
    ELSEIF ( Iway.EQ.-1 ) THEN
        theta = ATAN(SQRT(U(1)**2+U(2)**2))
        psi = ATAN2(U(2),U(1))
        R(1) = SIN(theta)*COS(psi)
        R(2) = SIN(theta)*SIN(psi)
        R(3) = COS(theta)
    ENDIF
    END

!______________________________________________________________________
    SUBROUTINE GPOLE(Iway,V,A,X1,X2)
    IMPLICIT NONE

    DOUBLE PRECISION A , cth , DOTTY , V , vc , vxa , X1 , X2
    INTEGER i , Iway

!  if iway = 1
!  converts unit vector x1 from a geographic Cartesian coord system to
!  one whose N-pole is given by v, a is any arbitrary vector in the
!  geographic system and converts to a unit vector which is at zero
!  longitude and correct latitude in the new system
!
!  if iway=-1 the inverse transformation from x2 to x1 is performed
!  v,a are always expressed in the geographic system
!
!
    DIMENSION V(3) , A(3) , X1(3) , X2(3) , vxa(3) , vc(3)
!
    cth = DOTTY(A,V)
    IF ( Iway.EQ.1 ) THEN
!  new z
        X2(3) = DOTTY(V,X1)
!  compute new y direction
        CALL AXB(V,A,vxa)
        CALL UNIT(vxa,vxa)
        X2(2) = DOTTY(vxa,X1)
!  compute new x
        DO i = 1 , 3
        vc(i) = V(i)
        ENDDO
        CALL SCALEY(cth,vc)
        DO i = 1 , 3
        vc(i) = A(i) - vc(i)
        ENDDO
        CALL UNIT(vc,vc)
        X2(1) = DOTTY(vc,X1)
!
    ELSEIF ( Iway.EQ.-1 ) THEN
        CALL AXB(V,A,vxa)
        CALL UNIT(vxa,vxa)
        DO i = 1 , 3
        vc(i) = V(i)
        ENDDO
        CALL SCALEY(cth,vc)
        DO i = 1 , 3
        vc(i) = A(i) - vc(i)
        ENDDO
        CALL UNIT(vc,vc)
        X1(1) = vc(1)*X2(1) + vxa(1)*X2(2) + V(1)*X2(3)
        X1(2) = vc(2)*X2(1) + vxa(2)*X2(2) + V(2)*X2(3)
        X1(3) = vc(3)*X2(1) + vxa(3)*X2(2) + V(3)*X2(3)
    ENDIF
    END

!_______________________________________________________________________
    SUBROUTINE AXB(A,B,C)
    IMPLICIT NONE

    DOUBLE PRECISION A , B , C

!  returns in c the vector cross product of 3-vectors a and b
    DIMENSION A(3) , B(3) , C(3)
    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)
    END

!_______________________________________________________________________
    SUBROUTINE UNIT(X,Ux)
    IMPLICIT NONE

    INTEGER i
    DOUBLE PRECISION unorm , Ux , X

!  converts x into a unit vector ux
    DIMENSION X(3) , Ux(3)
    unorm = 0.
    DO i = 1 , 3
        unorm = unorm + X(i)**2
    ENDDO
    unorm = SQRT(unorm)
    DO i = 1 , 3
        Ux(i) = X(i)/unorm
    ENDDO
    END

!_______________________________________________________________________
    SUBROUTINE SCALEY(C,X)
    IMPLICIT NONE

    DOUBLE PRECISION C , X
    INTEGER i

!  scales a vector x by the constant c
    DIMENSION X(3)
    DO i = 1 , 3
        X(i) = X(i)*C
    ENDDO
    END

!_______________________________________________________________________
    SUBROUTINE CENTR(X1,X2,X3,C)
    IMPLICIT NONE

    DOUBLE PRECISION C , cnorm , DOTTY , X1 , X2 , X3
    INTEGER j

    DIMENSION X1(3) , X2(3) , X3(3) , C(3)
!  computes centroid of triangle vertices x1,x2,x3
    DO j = 1 , 3
        C(j) = X1(j) + X2(j) + X3(j)
    ENDDO
    cnorm = SQRT(DOTTY(C,C))
    DO j = 1 , 3
        C(j) = C(j)/cnorm
    ENDDO
    END

!______________________________________________________________________
    SUBROUTINE BFIELD
    IMPLICIT NONE

    DOUBLE PRECISION A , BC , BRO , BS , P , REG , RO , SIGma , V ,   &
                    & XYW
    INTEGER i , IFLag , INP , IOUt , IPRint , ISTack , IV , k , l ,   &
        & LOBs , m , MAXE , MAXOB , MOBs , N , NEAr , NEPoch ,      &
        & NOBs , NQUad , NTRia
    INTEGER NVERT , NXYZ , NXYZ2

!
!  A sphere is defined by spherical triangular facets with
!  corners at the points defined by the 3-vectors  p  and a pointer
!  array  iv  specifying the corners to be connected with  ntria
!  triangles. The 4th component of p specifies the radial field at the
!  core mantle boundary
!  This routine computes the magnetic field at any point
!  outside the core by using the gram matrix computed in gram
!  (this calculation is assumed performed already).
!
!
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ,NVERT=10)
    PARAMETER (MAXOB=7000,MAXE=200)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    COMMON /OBSERV/ RO(3,MAXOB,MAXE) , BRO(MAXOB,MAXE) ,              &
                & SIGma(MAXOB,MAXE) , IFLag(MAXOB,MAXE) , MOBs ,    &
                & NOBs(MAXE) , NEPoch , LOBs(MAXE)
    COMMON /CALCUL/ BS(MAXOB,MAXE) , A(MAXOB,NXYZ) , BC(NXYZ2) ,      &
                & REG(NXYZ,NXYZ)
    COMMON /RULE  / XYW(3,100) , NQUad
    COMMON /CIRCLE/ V(4,NXYZ2) , ISTack(NXYZ2)
!
!  Zero the field vector
    l = 0
    DO m = 1 , NEPoch
        DO i = 1 , NOBs(m)
        BS(i,m) = 0.0
        IF ( IFLag(i,m).GT.0 ) THEN
            l = l + 1
            DO k = 1 , N
                BS(i,m) = BS(i,m) + A(l,k+(m-1)*N)*BC((m-1)*N+k)      &
                        & *SIGma(i,m)
            ENDDO
        ELSE
            BS(i,m) = 99999.
        ENDIF
        ENDDO
        l = l + 1
    ENDDO
    END

!_______________________________________________________________________
    SUBROUTINE GRAM
    IMPLICIT NONE
    DOUBLE PRECISION A , aa , ajac , b , BC , bjac , BRO , BS , c ,   &
                    & ELT , g , gamma , P , PAREA , REG , rg , RO ,    &
                    & rp , rpa , SIGma
    DOUBLE PRECISION ua , ub , uc , V , vc , xx , xxx , XYW
    INTEGER i , IFLag , INP , iob , IOUt , IPRint , ISTack , IV ,     &
        & iv1 , iv2 , iv3 , j , k , ll , LOBs , MAXE , MAXOB ,      &
        & MOBs , N , NEAr
    INTEGER NEPoch , nn , NOBs , NQUad , NTRia , NVERT , NXYZ , NXYZ2
!$$$$$ calls bgreen and bint, to compute green functions and for
!  interpolation of b on core surface
!  For each triangle on the core,
!  performs the numerical cubature of three surface integrals for br,
!  btheta, bphi
!  This provides the matrix a for the computation of the observed
!  surface fields, given the field at source points p(i,j)
!  core field br at vertices p, is in 4th component of arrays
!
!  The quadrature rule generated by 'triq' for integrals over
!  plane triangle with vertices  (0,0,0) (1,0,0) (0,1,0) is assumed
!  to be available in common /rule/.
!
!  The list of observer coordinates is in /observ/.
!
!  The results are returned in common /calcul/
!
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ,NVERT=10)
    PARAMETER (MAXOB=7000,MAXE=200)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    COMMON /OBSERV/ RO(3,MAXOB,MAXE) , BRO(MAXOB,MAXE) ,              &
                & SIGma(MAXOB,MAXE) , IFLag(MAXOB,MAXE) , MOBs ,    &
                & NOBs(MAXE) , NEPoch , LOBs(MAXE)
    COMMON /CALCUL/ BS(MAXOB,MAXE) , A(MAXOB,NXYZ) , BC(NXYZ2) ,      &
                & REG(NXYZ,NXYZ)
    COMMON /RULE  / XYW(3,100) , NQUad
    COMMON /CIRCLE/ V(4,NXYZ2) , ISTack(NXYZ2)
    DIMENSION rp(3) , g(3) , vc(4) , aa(3) , b(3) , c(3) , ua(2) ,    &
            & ub(2) , uc(2) , rg(2) , rpa(3) , gamma(3)
!
    WRITE (IOUt,*) 'entering gram'
    xxx = 0.
    nn = N*NEPoch
!  initialize matrix to zero
    DO j = 1 , nn
        DO iob = 1 , MOBs + 3
        A(iob,j) = 0.
        ENDDO
    ENDDO
    DO i = 1 , NEPoch
        WRITE (IOUt,*) ' Epoch ' , i , ' has ' , NOBs(i) ,             &
                    &' observations'
    ENDDO
!
    DO j = 1 , NTRia
!
!  vertices of triangle are points iv1,iv2,iv3
        iv1 = IV(1,j)
        iv2 = IV(2,j)
        iv3 = IV(3,j)
!
!  compute centroid vc of current triangle
!
        CALL CENTR(P(1,iv1),P(1,iv2),P(1,iv3),vc)
!	xxx=xxx + area(p(1,iv1),p(1,iv2),p(1,iv3))
!
!  rotate so vertices are in coord system with vc as pole, p(1,iv1) at
!  zero longitude
        CALL GPOLE(1,vc,P(1,iv1),P(1,iv1),aa)
        CALL GPOLE(1,vc,P(1,iv1),P(1,iv2),b)
        CALL GPOLE(1,vc,P(1,iv1),P(1,iv3),c)
!
!  find gnomonic projection of vertices and jacobian of transformation
        CALL GNOMON(1,ua,aa)
        CALL GNOMON(1,ub,b)
        CALL GNOMON(1,uc,c)
        ajac = PAREA(ua,ub,uc)*2.
!
!  Perform facet integration for all observation points.
!  Run through the cubature rule for this triangle for
!  each observation station.
!
        DO k = 1 , NQUad
!
!  find points on gnomonic surface corresponding to integration points
!
        rg(1) = ua(1) + (ub(1)-ua(1))*XYW(1,k) + (uc(1)-ua(1))      &
                & *XYW(2,k)
        rg(2) = ua(2) + (ub(2)-ua(2))*XYW(1,k) + (uc(2)-ua(2))      &
                & *XYW(2,k)
!
!  Find interpolant
!
        gamma(1) = 1. - XYW(1,k) - XYW(2,k)
        gamma(2) = XYW(1,k)
        gamma(3) = XYW(2,k)
!
!  find points in physical space corres. to rg
!
        CALL GNOMON(-1,rg,rpa)
        CALL GPOLE(-1,vc,P(1,iv1),rp,rpa)
!
!  Find Jacobian for gnomonic to spherical projection
!
        bjac = SQRT((1+rg(1)**2+rg(2)**2)**(-3))
!
!  evaluate green functions at rp for all observer points
!
        ll = 0
        DO i = 1 , NEPoch
            LOBs(i) = 0
            DO iob = 1 , NOBs(i)
                CALL BGREEN(RO(1,iob,i),rp,g)
!
!  Weight integral according to xyw rule
!
                IF ( IFLag(iob,i).GT.0 ) THEN
                    LOBs(i) = LOBs(i) + 1
                    ll = ll + 1
                    xx = ELT(IFLag(iob,i),g,BRO(iob,i))
                    A(ll,iv1+(i-1)*N) = A(ll,iv1+(i-1)*N) + xx*gamma(1)&
                    & *XYW(3,k)*ajac*bjac/SIGma(iob,i)
                    A(ll,iv2+(i-1)*N) = A(ll,iv2+(i-1)*N) + xx*gamma(2)&
                    & *XYW(3,k)*ajac*bjac/SIGma(iob,i)
                    A(ll,iv3+(i-1)*N) = A(ll,iv3+(i-1)*N) + xx*gamma(3)&
                    & *XYW(3,k)*ajac*bjac/SIGma(iob,i)
                ENDIF
            ENDDO
            ll = ll + 1
!  No monopoles allowed
            A(ll,iv1+(i-1)*N) = A(ll,iv1+(i-1)*N) + gamma(1)*XYW(3,k)&
                                & *ajac*bjac
            A(ll,iv2+(i-1)*N) = A(ll,iv2+(i-1)*N) + gamma(2)*XYW(3,k)&
                                & *ajac*bjac
            A(ll,iv3+(i-1)*N) = A(ll,iv3+(i-1)*N) + gamma(3)*XYW(3,k)&
                                & *ajac*bjac
        ENDDO
        ENDDO
    ENDDO
!	write(iout,*)'sum of s. triangle areas',xxx
    END
!______________________________________________________________________
    SUBROUTINE BGREEN(Ro,So,G)
    IMPLICIT NONE
    DOUBLE PRECISION amu , c , fpi , G , p , r , rho , rhocu , rhosq ,&
                    & Ro , rrr , So , sounit , t , tt , uphi , uthet , &
                    & x
    INTEGER i
!  computes greens function for magnetic field components br,bthet, bphi
!  at observer point ro due to source point so.
!  ro specified as r, theta, phi coords
!  so specified as x, y, z , geocentric.
!
    DIMENSION Ro(3) , So(3) , G(3)
    DATA fpi/12.56637062/ , c/1.0/
    t = Ro(2)
    p = Ro(3)
    sounit = SQRT(So(1)**2+So(2)**2+So(3)**2)
    DO i = 1 , 3
        So(i) = So(i)/sounit
    ENDDO
    rho = c/Ro(1)
    rhosq = rho*rho
    rhocu = rhosq*rho
!  compute runit.runit'
    amu = SIN(t)*(COS(p)*So(1)+SIN(p)*So(2)) + COS(t)*So(3)
    r = SQRT(1.0-2.0*amu*rho+rhosq)
    rrr = r**3
    tt = 1 + r - amu*rho
!  compute runit'.thetaunit and runit'.phiunit
    uthet = COS(t)*(COS(p)*So(1)+SIN(p)*So(2)) - SIN(t)*So(3)
!
    uphi = -SIN(p)*So(1) + COS(p)*So(2)
!
!  compute g(1),g(2),g(3), green functions for r,theta,phi
    G(1) = rhosq*((1.0-rhosq)/rrr-1.)/fpi
    x = rhocu*(1.0+2.0*r-rhosq)/(rrr*tt*fpi)
    G(2) = -uthet*x
    G(3) = -uphi*x
    END
!_______________________________________________________________________
    FUNCTION PAREA(A,B,C)
    IMPLICIT NONE
    DOUBLE PRECISION A , ab , ac , ar , B , C , DOTTY , PAREA
    INTEGER j
!  computes area of planar triangle
    DIMENSION A(2) , B(2) , C(2) , ab(3) , ac(3) , ar(3)
    DO j = 1 , 2
        ab(j) = B(j) - A(j)
        ac(j) = C(j) - A(j)
    ENDDO
    ab(3) = 0.
    ac(3) = 0.
    CALL AXB(ab,ac,ar)
    PAREA = SQRT(DOTTY(ar,ar))/2.
    END
!______________________________________________________________________
    FUNCTION ELT(Iflag,G,Bro)
    IMPLICIT NONE
    DOUBLE PRECISION Bro , drad , ELT , G , h
    INTEGER Iflag , INP , IOUt , IPRint
!  converts green function for br, btheta, bphi given in g to other
!  elements fo the geomagnetic field
!	iflag=1  declination
!	      2  inclination
!	      3  horizontal intensity
!	      4  x or north component
!	      5  y or east component
!	      6  z or vertical down component
!	      7  f total intensity
!
    COMMON /IO    / INP , IOUt , IPRint
    DIMENSION G(3)
    DATA drad/57.2958/
!	write(iout,*)'iflag=',iflag
    IF ( Iflag.EQ.1 ) THEN
        ELT = G(2)*SIN(Bro/drad) + G(3)*COS(Bro/drad)
    ELSEIF ( Iflag.EQ.2 ) THEN
        h = SQRT(G(2)**2+G(3)**2)
        ELT = drad*ATAN2(-G(1),h)
    ELSEIF ( Iflag.EQ.3 ) THEN
        ELT = SQRT(G(2)**2+G(3)**2)
    ELSEIF ( Iflag.EQ.4 ) THEN
        ELT = -G(2)
    ELSEIF ( Iflag.EQ.5 ) THEN
        ELT = G(3)
    ELSEIF ( Iflag.EQ.6 ) THEN
        ELT = -G(1)
    ELSEIF ( Iflag.EQ.7 ) THEN
        ELT = SQRT(G(1)**2+G(2)**2+G(3)**2)
    ELSE
        WRITE (IOUt,*) 'Unknown measurement type'
    ENDIF
    END
!--------------------------------------------------------------
    SUBROUTINE FLUX(I)
    IMPLICIT NONE
    DOUBLE PRECISION aa , ajac , b , BINT , bjac , c , FLUxv , g ,    &
                    & gam1 , GRAdm , P , PAREA , rg , rp , rpa , SINT ,&
                    & ua , ub , uc , V
    DOUBLE PRECISION vc , XYW , yes
    INTEGER I , INP , IOUt , IPP , IPRint , IPTs , ISIgn , ISTack ,   &
        & ITRis , IV , j , k , kkk , KSIGNM , LLPt , LPTr , LPTs ,  &
        & LTRis , MAXE , MAXP
    INTEGER N , NEAr , NP , NPAtch , NPTs , NQUad , NTRia , NTRis ,   &
        & NVERT , NXYZ , NXYZ2
!
!  Computes the flux through the designated regions  on the core,
!  and the gradient matrix for the flux through the curves with
!  respect to changes in the values of the field at the vertices
!  on the core.
!
!  Assumes patch has already been called to define the variables
!  in /patchy/.
!
!
!  Outline:
!
!  Initialize the flux vector and flux gradient matrix to zero
!  Loop in triangles
!    for each triangle, compute rotation and gnomonic projection
!    loop in integration sample points
!       compute preimage of the sample point on the core
!       compute the jacobian for the gnomonic at that sample point
!       compute field value at that point
!       add the value, weighted by the int. wt. and the jacobians,
!           to the positive patch or negative patch, according
!           to its sign.
!       add interp. wt*int. wt.* jacobian to the appropriate elements
!           of the grad. matrix (according to signs of br and vertex)
!    next integration sample point
!  next triangle.
!
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ,MAXP=30,NVERT=10,MAXE=200)
    DIMENSION rp(3) , g(3) , vc(4) , aa(3) , b(3) , c(3) , ua(2) ,    &
            & ub(2) , uc(2) , rg(2) , rpa(3)
    COMMON /PATCHY/ NP(MAXE) , NPTs(MAXP,MAXE) , ISIgn(MAXP,MAXE) ,   &
                & IPTs(MAXP,MAXE) , LPTs(2*NXYZ,MAXE) , NTRis(NXYZ) &
                & , ITRis(NXYZ+1) , LTRis(2*(2*NXYZ-4)) ,           &
                & LPTr(2*NXYZ-4,2) , LLPt(MAXP,MAXE) , NPAtch ,     &
                & IPP(NXYZ)
    COMMON /PATCHG/ FLUxv(MAXP,MAXE) , GRAdm(MAXP,NXYZ,MAXE)
!
!  see subroutine patch for definitions of the variables in /patchy/
!
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /RULE  / XYW(3,100) , NQUad
    COMMON /CIRCLE/ V(4,NXYZ2) , ISTack(NXYZ2)
!
!  initialize flux vector and gradient matrix.
!******
    WRITE (IOUt,*) ' arrived in flux '
!      do 10 j=1,np(i)
!      fluxv(j,i)=0.
!      do 10 k=1,n
! 10   gradm(j,k,i)=0.
!
!  Loop in patches, then in triangles within each patch.
    DO kkk = 1 , NP(I)
! Flux in patch kkk at epoch i
        FLUxv(kkk,I) = 0.D0
        DO j = LTRis(ITRis(kkk)) , LTRis(ITRis(kkk)+NTRis(kkk)-1)
!
!
!  compute centroid vc of current triangle
        CALL CENTR(P(1,IV(1,j)),P(1,IV(2,j)),P(1,IV(3,j)),vc)
!
!  rotate so vertices are in coord system with vc as pole, p(1,iv1) at
!  zero longitude
        CALL GPOLE(1,vc,P(1,IV(1,j)),P(1,IV(1,j)),aa)
        CALL GPOLE(1,vc,P(1,IV(1,j)),P(1,IV(2,j)),b)
        CALL GPOLE(1,vc,P(1,IV(1,j)),P(1,IV(3,j)),c)
!
!  find gnomonic projection of vertices and jacobian of transformation
        CALL GNOMON(1,ua,aa)
        CALL GNOMON(1,ub,b)
        CALL GNOMON(1,uc,c)
        ajac = PAREA(ua,ub,uc)*2.
!*******
!      write(iout,*)' entering integration loop '
!
!  Loop in integration sample points.
        DO k = 1 , NQUad
!
!  find points on gnomonic surface corresponding to integration points
            rg(1) = ua(1) + (ub(1)-ua(1))*XYW(1,k) + (uc(1)-ua(1))   &
                    & *XYW(2,k)
            rg(2) = ua(2) + (ub(2)-ua(2))*XYW(1,k) + (uc(2)-ua(2))   &
                    & *XYW(2,k)
!
!  Find Jacobian for gnomonic to spherical projection
            bjac = SQRT((1+rg(1)**2+rg(2)**2)**(-3))
!
!  Find interpolating vector g(.)---linear interpolation in gnomonic
!  plane.
            g(1) = 1. - XYW(1,k) - XYW(2,k)
            g(2) = XYW(1,k)
            g(3) = XYW(2,k)
!
!  find the point rp on the core corresponding to rg.
            CALL GNOMON(-1,rg,rpa)
            CALL GPOLE(-1,vc,P(1,IV(1,j)),rp,rpa)
!
!  Find field value and sign of patch.
            gam1 = BINT(rp,rg,k,P(1,IV(1,j)),P(1,IV(2,j)),           &
                & P(1,IV(3,j)))
            yes = SINT(rp,rg,k,P(1,IV(1,j)),P(1,IV(2,j)),P(1,IV(3,j))&
                & ,IPP(IV(1,j)),IPP(IV(2,j)),IPP(IV(3,j)))
            IF ( yes.LE.0 .AND. KSIGNM(ISIgn(kkk,I)).LT.0 ) THEN
                yes = 1.D0
            ELSEIF ( yes.GT.0 .AND. KSIGNM(ISIgn(kkk,I)).EQ.1 ) THEN
                yes = 1.D0
            ELSE
                yes = 0.0
            ENDIF
!
!  Update the flux vector fluxv and the gradient matrix gradm.
!*******
!      write(iout,*)' tri ',j,' sample ', k,' br ',br, ' g ',
!     +  (g(ii),ii=1,3),' ndx ',ndx,' ajac ',ajac,' bjac ',bjac
!*******
            FLUxv(kkk,I) = FLUxv(kkk,I) + ajac*bjac*XYW(3,k)*gam1*yes
!	write(*,*) j, iv(1,j), ksignm(iv(1,j)),ndx, fluxv(ndx,i),i,br
!            do 40 l=1,3
!               if(signm(iv(l,j)).eq.sgnbr)
!              gradm(lptr(j,ndx),iv(l,j),i)=gradm(lptr(j,ndx),
!     $         iv(l,j),i)+g(l)*ajac*bjac*xyw(3,k)
! 40         continue
!
!  Next integration sample point.
        ENDDO
!
!  Next triangle.
        ENDDO
    ENDDO
    END
!-------------------------------------------------------------
    SUBROUTINE PATCH(I)
    IMPLICIT NONE
    DOUBLE PRECISION FLUxv , GRAdm , P
    INTEGER I , INP , IOUt , IPP , IPRint , IPTs , ISIgn , ITRis ,    &
        & IV , j , jj , k , kk , KSIGNM , l , LLPt , LPTr , LPTs ,  &
        & LTRis , m
    INTEGER MAXE , MAXP , mm , N , ndx , NEAr , next , now , NP ,     &
        & NPAtch , NPTs , NTRia , NTRis , NVERT , NXYZ , NXYZ2
!
!  constructs lists of triangles and vertices associated with
!  each patch on the core on which Br has one sign
!
!  constructs a list of the signs of the patches.
!  constructs an array indicating which points are in which
!  patches, and how many patches a point appears in.
!
!  assumes the  near  array has been previously constructed.
!
!* calls signm
!
    PARAMETER (NXYZ=20000,NVERT=10,NXYZ2=2*NXYZ)
    PARAMETER (MAXP=30,MAXE=200)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /PATCHY/ NP(MAXE) , NPTs(MAXP,MAXE) , ISIgn(MAXP,MAXE) ,   &
                & IPTs(MAXP,MAXE) , LPTs(2*NXYZ,MAXE) , NTRis(NXYZ) &
                & , ITRis(NXYZ+1) , LTRis(2*(2*NXYZ-4)) ,           &
                & LPTr(2*NXYZ-4,2) , LLPt(MAXP,MAXE) , NPAtch ,     &
                & IPP(NXYZ)
    COMMON /PATCHG/ FLUxv(MAXP,MAXE) , GRAdm(MAXP,NXYZ,MAXE)
!
!  np(i):                   total # patches of one sign for epoch i
!  npts(j,i):              # vertices (points) in patch j for epoch i
!  isign(j,i):             sign of patch j for epoch i
!  ipts(j,i):              starting index in lpts of the list of points
!                        in patch j, for epoch i
!  lpts(ipts(j)+k-1,i):    kth vertex in patch j for epoch i
!  ntris(j):             # triangles in patch j
!  itris(j):             starting index in ltris of the list of
!                        triangles in patch j
!  ltris(itris(j)+k-1):  kth triangle in patch j
!  lptr(j,k):            k=1:  positive patch (if any) associated with
!                        triangle j.  k=2: negative patch (if any)
!  fluxv(j,i):  		 flux through patch j for epoch i
!  gradm(j,k,i): 		 gradient of flux through patch j due to change
!			 in field at vertex k for epoch i
!
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
!
!  n:                    total # points on the core
!  ntria:                total # triangles on the core
!  p(k,j):               kth component of jth point, k=1,2,3;
!                        Br(jth pt), k=4
!  iv(k,j):              kth vertex in jth triangle
!  near(k,j):            kth nearest neighbor of jth point
!
!  initialize indices.
    ITRis(1) = 1
    NP(I) = 0
!
    WRITE (IOUt,*) ' Entering patch' , I
!	  do 1 k=1,maxp
!	  do 1 j=1,maxe
!1	  npts(k,i)=0
!
!  Are all vertices in a list?
    DO j = 1 , N
        DO k = 1 , NP(I)
        DO m = 1 , NPTs(k,I)
            IF ( LPTs(IPTs(k,I)+m-1,I).EQ.j ) GOTO 100
        ENDDO
        ENDDO
!  Vertex j is not on a list.  Make a new patch with j as its leading
!  point.
        NP(I) = NP(I) + 1
        IF ( NP(I).GT.30 ) THEN
        WRITE (IOUt,*) 'too many patches, cant go on'
        STOP
        ENDIF
        NPTs(NP(I),I) = 1
        ISIgn(NP(I),I) = KSIGNM(j)
        IF ( NP(I).EQ.1 ) THEN
        IPTs(NP(I),I) = 1
        ELSE
        IPTs(NP(I),I) = IPTs(NP(I)-1,I) + NPTs(NP(I)-1,I)
        ENDIF
        LPTs(IPTs(NP(I),I),I) = j
!
!  List all neighbors, neighbors of neighbors, etc. with same sign Br.
        DO jj = 1 , N
        IF ( jj.GT.NPTs(NP(I),I) ) GOTO 100
        now = LPTs(IPTs(NP(I),I)+jj-1,I)
!  now is the current point.  Examine  now 's neighbors.
        DO kk = 1 , NVERT
            next = NEAr(kk,now)
!  next is candidate neighbor.  If next=0,  now  has no more neighbors.
            IF ( next.EQ.0 ) GOTO 50
            IF ( KSIGNM(next).EQ.KSIGNM(now) ) THEN
!  next should be on the list---is it already?
                DO mm = 1 , NPTs(NP(I),I)
                    IF ( LPTs(IPTs(NP(I),I)+mm-1,I).EQ.next ) GOTO 20
                ENDDO
!  Add  next  to  lpts.
                NPTs(NP(I),I) = NPTs(NP(I),I) + 1
                LPTs(IPTs(NP(I),I)+NPTs(NP(I),I)-1,I) = next
            ENDIF
20         ENDDO
50      ENDDO
100  ENDDO
    DO j = 1 , NP(I)
        WRITE (IOUt,*) ' patch ' , j , ' number of points ' , NPTs(j,I)&
                    & , ' sign ' , ISIgn(j,I)
        WRITE (IOUt,*) ' points are ' ,                                &
                    & (LPTs(IPTs(j,I)+k-1,I),k=1,NPTs(j,I))
    ENDDO
!
!  Make a nonredundant list of all triangles including any point in
!  patch j.
    DO j = 1 , NP(I)
        NTRis(j) = 0
        DO k = 1 , NPTs(j,I)
        DO m = 1 , NTRia
            DO l = 1 , 3
                IF ( IV(l,m).EQ.LPTs(IPTs(j,I)+k-1,I) ) THEN
!  see if triangle m is already listed
                    DO jj = 1 , NTRis(j)
                    IF ( LTRis(ITRis(j)+jj-1).EQ.m ) GOTO 120
                    ENDDO
!  it isn't; add it
                    NTRis(j) = NTRis(j) + 1
                    LTRis(ITRis(j)+NTRis(j)-1) = m
                    GOTO 120
                ENDIF
            ENDDO
120        ENDDO
        ENDDO
!
!  find the starting index of the set of triangles for the next patch
        ITRis(j+1) = ITRis(j) + NTRis(j)
    ENDDO
!
    DO j = 1 , NP(I)
!	 write(iout,*)' patch ',j,' number of triangles ',ntris(j)
!	 write(iout,*)' triangles are ',(ltris(itris(j)+k-1),k=1,ntris(j))
    ENDDO
!  For each patch, make a nonredundant list of points in the
!  triangles that intersect it.
    DO j = 1 , NP(I)
        NPTs(j,I) = 0
        DO k = 1 , NTRis(j)
        DO m = 1 , 3
!  Check if the next point is already on the list
            DO jj = 1 , NPTs(j,I)
                IF ( LPTs(IPTs(j,I)+jj-1,I)                           &
                    & .EQ.IV(m,LTRis(ITRis(j)+k-1)) ) GOTO 140
            ENDDO
!  it isn't; add it
            NPTs(j,I) = NPTs(j,I) + 1
            LPTs(IPTs(j,I)+NPTs(j,I)-1,I) = IV(m,LTRis(ITRis(j)+k-1))
140        ENDDO
        ENDDO
!
!  Sort the list of points associated with each patch
        CALL ISORT(NPTs(j,I),LPTs(IPTs(j,I),I))
!
!  Update the point-patch index
        IPTs(j+1,I) = IPTs(j,I) + NPTs(j,I)
    ENDDO
!
!  Make lptr, the cross-reference array of which patches a triangle
!  appears in.
    DO j = 1 , NTRia
        LPTr(j,1) = 0
        LPTr(j,2) = 0
    ENDDO
    DO j = 1 , NP(I)
        IF ( ISIgn(j,I).GT.0 ) THEN
!  Patch j is positive and should be in the first column of lptr.
        ndx = 1
        ELSE
!  Patch j should be in the second column of lptr.
        ndx = 2
        ENDIF
        DO k = 1 , NTRis(j)
        l = LTRis(ITRis(j)+k-1)
        LPTr(l,ndx) = j
        ENDDO
    ENDDO
    END
!----------------------------------------------------------
    FUNCTION SIGNM(I)
    IMPLICIT NONE
    INTEGER I , IV , N , NEAr , NTRia , NVERT , NXYZ , NXYZ2
    DOUBLE PRECISION P , SIGNM
!
!  returns the sign of Br at the ith vertex.
!* Calls no other routines.
!
    PARAMETER (NXYZ=20000,NVERT=10,NXYZ2=2*NXYZ)
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    SIGNM = -1
    IF ( P(4,I).GE.0.0 ) SIGNM = 1
    END
!-----------------------------------------------------------
!----------------------------------------------------------
    INTEGER FUNCTION KSIGNM(I)
    IMPLICIT NONE
    INTEGER I , IPP , IPTs , ISIgn , ITRis , IV , LLPt , LPTr , LPTs ,&
        & LTRis , MAXE , MAXP , N , NEAr , NP , NPAtch , NPTs ,     &
        & NTRia , NTRis , NVERT
    INTEGER NXYZ , NXYZ2
    DOUBLE PRECISION P
!
!  returns the sign of ipp at the ith vertex.
!* Calls no other routines.
!
    PARAMETER (NXYZ=20000,NVERT=10,NXYZ2=2*NXYZ)
    PARAMETER (MAXP=30,MAXE=200)
    COMMON /PATCHY/ NP(MAXE) , NPTs(MAXP,MAXE) , ISIgn(MAXP,MAXE) ,   &
                & IPTs(MAXP,MAXE) , LPTs(2*NXYZ,MAXE) , NTRis(NXYZ) &
                & , ITRis(NXYZ+1) , LTRis(2*(2*NXYZ-4)) ,           &
                & LPTr(2*NXYZ-4,2) , LLPt(MAXP,MAXE) , NPAtch ,     &
                & IPP(NXYZ)
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    KSIGNM = -1
    IF ( IPP(I).GE.0 ) KSIGNM = 1
    END
!-----------------------------------------------------------
!_______________________________________________________________________
    SUBROUTINE SHEXP(Ldeg,Rr)
    IMPLICIT NONE
    DOUBLE PRECISION aa , ajac , b , BINT , bjac , c , const , ct ,   &
                    & DLM , ELF , gam1 , GLM , HLM , P , PAREA , part ,&
                    & phi , pi , PLM , r
    DOUBLE PRECISION rc , rg , rp , rpa , Rr , rrn , theta , ua , ub ,&
                    & uc , vc , xxx , XYW
    INTEGER i , il , im , INP , IOUt , IPRint , IV , iv1 , iv2 , iv3 ,&
        & j , k , Ldeg , LMAX , N , NEAr , NQUad , NTRia , NVERT ,  &
        & NXYZ
    INTEGER NXYZ2
!  For each triangle on the core,
!  performs the numerical cubature of br core model against partially
!  normalised spherical harmonics to get sh rep'n out to degree and
!  order ldeg at radius rr
!  input rr in km is normalized so that core radius =1
!  core field br at vertices p, is in 4th component of arrays
!
!  The quadrature rule generated by 'triq' for integrals over
!  plane triangle with vertices  (0,0,0) (1,0,0) (0,1,0) is assumed
!  to be available in common /rule/.
!
!  The list of observer coordinates is in /observ/.
!
!  The results are returned in common /calcul/
!
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ,NVERT=10)
    PARAMETER (LMAX=101)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    COMMON /RULE  / XYW(3,100) , NQUad
    COMMON /SHR   / GLM(LMAX,LMAX) , HLM(LMAX,LMAX) , PLM(LMAX,LMAX) ,&
                & DLM(LMAX,LMAX)
    DIMENSION rp(3) , vc(4) , aa(3) , b(3) , c(3) , ua(2) , ub(2) ,   &
            & uc(2) , rg(2) , rpa(3)
    DATA pi/3.14159265/ , rc/3486.D0/
    rrn = Rr/rc
!
    xxx = 0.
    IF ( Ldeg.GE.LMAX ) THEN
        WRITE (IOUt,*) 's.h.expansion truncated to ' , LMAX - 1
        Ldeg = LMAX - 1
    ENDIF
!  initialize she to zero
    DO j = 1 , Ldeg
        DO i = 1 , Ldeg
        GLM(i,j) = 0.
        HLM(i,j) = 0.
        ENDDO
    ENDDO
!
    DO j = 1 , NTRia
!
!  vertices of triangle are points iv1,iv2,iv3
        iv1 = IV(1,j)
        iv2 = IV(2,j)
        iv3 = IV(3,j)
!
!  compute centroid vc of current triangle
!
        CALL CENTR(P(1,iv1),P(1,iv2),P(1,iv3),vc)
!	xxx=xxx + area(p(1,iv1),p(1,iv2),p(1,iv3))
!
!  rotate so vertices are in coord system with vc as pole, p(1,iv1) at
!  zero longitude
        CALL GPOLE(1,vc,P(1,iv1),P(1,iv1),aa)
        CALL GPOLE(1,vc,P(1,iv1),P(1,iv2),b)
        CALL GPOLE(1,vc,P(1,iv1),P(1,iv3),c)
!
!  find gnomonic projection of vertices and jacobian of transformation
!  for plane to unit triangle
!
        CALL GNOMON(1,ua,aa)
        CALL GNOMON(1,ub,b)
        CALL GNOMON(1,uc,c)
        ajac = PAREA(ua,ub,uc)*2.
!
!  integrate over triangle
!
        DO k = 1 , NQUad
!
!  find points on gnomonic surface corresponding to integration points
!
        rg(1) = ua(1) + (ub(1)-ua(1))*XYW(1,k) + (uc(1)-ua(1))      &
                & *XYW(2,k)
        rg(2) = ua(2) + (ub(2)-ua(2))*XYW(1,k) + (uc(2)-ua(2))      &
                & *XYW(2,k)
!
!  find points in physical space corres. to rg
!
        CALL GNOMON(-1,rg,rpa)
        CALL GPOLE(-1,vc,P(1,iv1),rp,rpa)
!
!  convert to colat, longitude
!
        CALL RTHPHI(1,rp(1),rp(2),rp(3),r,theta,phi)
        ct = COS(theta)
!
!  Evaluate Legendre polynomials at this site, and field
!
        CALL PLMDLM(ct,Ldeg,LMAX,PLM,DLM)
        gam1 = BINT(rp,rg,k,P(1,iv1),P(1,iv2),P(1,iv3))
!
!  Find Jacobian for gnomonic to spherical projection
!
        bjac = SQRT((1+rg(1)**2+rg(2)**2)**(-3))
!
!  evaluate partially normalized spherical harmonic functions at rp
!
!  Weight integral according to xyw rule
!
        DO il = 0 , Ldeg
            im = 0
            const = FLOAT(2*il+1)/(rrn**(il+2)*FLOAT(il+1)*4.0*pi)
            GLM(il+1,im+1) = GLM(il+1,im+1)                          &
                            & + gam1*const*PLM(il+1,im+1)             &
                            & *COS(FLOAT(im)*phi)*ajac*bjac*XYW(3,k)
            HLM(il+1,im+1) = HLM(il+1,im+1)                          &
                            & + gam1*const*PLM(il+1,im+1)             &
                            & *SIN(FLOAT(im)*phi)*ajac*bjac*XYW(3,k)
            DO im = 1 , il
                part = SQRT(2.*ELF(il-im)/ELF(il+im))*(-1)**im
                GLM(il+1,im+1) = GLM(il+1,im+1)                       &
                                & + gam1*const*PLM(il+1,im+1)          &
                                & *COS(FLOAT(im)*phi)                  &
                                & *ajac*bjac*XYW(3,k)*part
                HLM(il+1,im+1) = HLM(il+1,im+1)                       &
                                & + gam1*const*PLM(il+1,im+1)          &
                                & *SIN(FLOAT(im)*phi)                  &
                                & *ajac*bjac*XYW(3,k)*part
            ENDDO
        ENDDO
        ENDDO
    ENDDO
!	print *,'sum of s. triangle areas',xxx
    END
!
!______________________________________________________________________
    FUNCTION ELF(L)
    IMPLICIT NONE
    DOUBLE PRECISION ELF
    INTEGER i , L
!  Finds l!
    ELF = 1.0
    DO i = 1 , L
        ELF = ELF*FLOAT(i)
    ENDDO
    END
!_______________________________________________________________________
    SUBROUTINE PLMDLM(C,Lmax,Ldim,Plm,Dlm)
    IMPLICIT NONE
    DOUBLE PRECISION C , Dlm , elms , Plm , s , twos
    INTEGER l , Ldim , lm , Lmax , m , mmax
!$$$$$ calls no other routines
!  finds associated legendre functions of all orders and degrees for
!  a fixed argument by recurrence.  the sign convention is that of
!  abramowitz & stegun, so that plm(l=1,m=1) = - sin theta.  precisely
!    plm = (-1)**m * (1-x**2)**(m/2) * (d/dx)**m pl(x),
!  where  pl(x)  is the legendre polynomial of degree  l.  special care
!  has been taken to avoid loss of precision or failure at or near  c=1.
!
!  c     is the (real) argument,  c .le. 1.0 in magnitude
!  lmax  is the maximum degree desired.
!  ldim  is the 1st dimension of both  plm and dlm  in the caller.  note
!        ldim .ge. lmax+1
!  plm()  array for results.  plm  of degree  l  and order  m  is held
!        in element  plm(l+1,m+1).  the elements  plm(l,l+1)  are used,
!        but other elements above the diagonal are not disturbed.
!  dlm() array for derivatives  (d/d theta) plm,  where  c=cos theta.
!
    DIMENSION Plm(Ldim,Lmax) , Dlm(Ldim,Lmax)
!
    Plm(1,1) = 1.0
    Dlm(1,1) = 0.0
    IF ( Lmax.EQ.0 ) RETURN
!  zero out super-diagonal.
    DO l = 1 , Lmax
        Plm(l,l+1) = 0.0
    ENDDO
    s = SQRT(DMAX1(0.D0,1.D0-C*C))
    twos = 2.0*s
!  recur diagonally.
    DO l = 1 , Lmax
        elms = l*s
        mmax = Lmax - l + 1
        DO m = 1 , mmax
        lm = l + m
        Plm(lm,m+1) = C*Plm(lm-1,m+1) - elms*Plm(lm-1,m)
        elms = twos + elms
        ENDDO
        Plm(l+1,1) = C*Plm(l,1) + s*Plm(l,2)/l
    ENDDO
!  find derivatives.
    DO l = 1 , Lmax
        Dlm(l+1,1) = Plm(l+1,2)
        DO m = 1 , l
        IF ( m.LT.l ) Dlm(l+1,m+1)                                  &
            & = 0.5*(Plm(l+1,m+2)-(l+m)*(l-m+1)*Plm(l+1,m))
        IF ( m.EQ.l ) Dlm(l+1,l+1) = -l*Plm(l+1,l)
        ENDDO
    ENDDO
    END
!
!_______________________________________________________________________
    SUBROUTINE RTHPHI(Iway,X,Y,Z,R,Theta,Phi)
    IMPLICIT NONE
    INTEGER iout , Iway
    DOUBLE PRECISION Phi , R , Theta , X , Y , Z
!
!  if iway=1 transforms x,y,z to r, theta, phi
!  if iway=-1 inverse transformation is performed
!
    IF ( Iway.EQ.1 ) THEN
        R = SQRT(X**2+Y**2+Z**2)
        Theta = ACOS(Z/R)
        Phi = ATAN2(Y,X)
    ELSEIF ( Iway.EQ.-1 ) THEN
        X = R*SIN(Theta)*COS(Phi)
        Y = R*SIN(Theta)*SIN(Phi)
        Z = R*COS(Theta)
    ELSE
        WRITE (iout,*) 'invalid subroutine call to rthphi'
    ENDIF
    END
!_______________________________________________________________________
    FUNCTION SINT(Rp,Rg,K,P1,P2,P3,Ipp1,Ipp2,Ipp3)
    IMPLICIT NONE
    DOUBLE PRECISION gamma , P1 , P2 , P3 , Rg , Rp , SINT , XYW
    INTEGER Ipp1 , Ipp2 , Ipp3 , K , NQUad
!  interpolates patch identity yes or no to a point rg on the gnomonic projection
    COMMON /RULE  / XYW(3,100) , NQUad
    DIMENSION Rp(3) , Rg(2) , gamma(3) , P1(4) , P2(4) , P3(4)
    gamma(1) = 1. - XYW(1,K) - XYW(2,K)
    gamma(2) = XYW(1,K)
    gamma(3) = XYW(2,K)
    SINT = gamma(1)*Ipp1 + gamma(2)*Ipp2 + gamma(3)*Ipp3
    END
!_______________________________________________________________________
    FUNCTION BINT(Rp,Rg,K,P1,P2,P3)
    IMPLICIT NONE
    DOUBLE PRECISION BINT , gamma , P1 , P2 , P3 , Rg , Rp , XYW
    INTEGER K , NQUad
!  interpolates b to a point rg on the gnomonic projection
    COMMON /RULE  / XYW(3,100) , NQUad
    DIMENSION Rp(3) , Rg(2) , gamma(3) , P1(4) , P2(4) , P3(4)
    gamma(1) = 1. - XYW(1,K) - XYW(2,K)
    gamma(2) = XYW(1,K)
    gamma(3) = XYW(2,K)
    BINT = gamma(1)*P1(4) + gamma(2)*P2(4) + gamma(3)*P3(4)
    END
!_______________________________________________________________________
    SUBROUTINE OBS(Kdim)
    IMPLICIT NONE
    DOUBLE PRECISION BRO , drad , h , RO , SIGma , val , x , y
    INTEGER i , IFLag , INP , IOUt , IPRint , j , k , Kdim , LOBs ,   &
        & m , MAXE , MAXOB , MOBs , NEPoch , nfound , NOBs , nwant1
!  if kdim=3 reads the observation points for the magnetic field into
!  array ro.
!  if kdim=6 reads magnetic field values at those points as well
!  ro(i,j,k) :  i=1,3 position of jth observation (r,theta,phi),kth
!               epoch
!  bro(i,j,k) :  i=1,3 magnetic field at jth observation point (br,
!  btheta,bphi in local coords), kth epoch
    PARAMETER (MAXOB=7000,MAXE=200)
    CHARACTER*80 filnam(2) , dform(2)
    COMMON /IO    / INP , IOUt , IPRint
!
    COMMON /OBSERV/ RO(3,MAXOB,MAXE) , BRO(MAXOB,MAXE) ,              &
                & SIGma(MAXOB,MAXE) , IFLag(MAXOB,MAXE) , MOBs ,    &
                & NOBs(MAXE) , NEPoch , LOBs(MAXE)
    DIMENSION x(3) , y(3)
    DIMENSION val(10)
    DATA drad/57.295780/
    MOBs = 0
    CALL GETVAL('epochs',val,nwant1,nfound)
    NEPoch = NINT(val(1))
    DO i = 1 , NEPoch
        NOBs(i) = 0
    ENDDO
    IF ( Kdim.EQ.3 ) THEN
        CALL GETCHR('e1file',filnam(1),nfound)
        OPEN (UNIT=16,FILE=filnam(1))
        CALL GETCHR('e1format',dform(1),nfound)
        CALL GETVAL('rmse',SIGma,NEPoch,nfound)
        DO j = 1 , MAXOB + 1
        READ (16,dform(1),ERR=50,END=50) (RO(i,j,1),i=1,3)
        MOBs = MOBs + 1
        ENDDO
50      WRITE (IOUt,*) MOBs , 'observation points read'
    ELSEIF ( Kdim.EQ.6 ) THEN
        CALL GETCHR('e1file',filnam(1),nfound)
        CALL GETCHR('e1format',dform(1),nfound)
        IF ( NEPoch.EQ.2 ) THEN
        CALL GETCHR('e2file',filnam(2),nfound)
        CALL GETCHR('e2format',dform(2),nfound)
        ENDIF
        DO k = 1 , NEPoch
        OPEN (UNIT=16,FILE=filnam(k))
        IF ( dform(k).EQ.'()' ) THEN
!  Assumed format is one component per line, observatory type
!  Ellipticity correction applied.
            DO j = 1 , MAXOB + 1
                READ (16,*,END=60,ERR=60) IFLag(j,k) , RO(2,j,k) ,    &
                    & RO(3,j,k) , RO(1,j,k) , BRO(j,k) , SIGma(j,k)
                IF ( IFLag(j,k).GE.0 ) THEN
!  manipulate to form required here
                    RO(3,j,k) = RO(3,j,k)/drad
                    h = RO(1,j,k)/1000.
                    RO(1,j,k) = 6378.139
                    RO(2,j,k) = RO(2,j,k)/drad
!  do ellipticity correction
                    CALL ELLIPT(RO(1,j,k),RO(2,j,k),RO(1,j,k),RO(2,j,k)&
                            & )
                    RO(1,j,k) = (RO(1,j,k)+h)/3485
                    NOBs(k) = NOBs(k) + 1
                ENDIF
            ENDDO
            WRITE (IOUt,'(a)') 'data file ' , filnam
        ELSE
!  Assumed format is satellite data, x,y,z, manipulated to above single
!  component form
!  Read measurement uncertainty from input file
            CALL GETVAL('rmse',SIGma(1,k),nwant1,nfound)
            DO j = 1 , MAXOB + 1
                READ (16,dform(k),END=60,ERR=60) y(2) , y(3) , y(1) , &
                    & x(1) , x(2) , x(3)
                DO m = 1 , 3
                    IF ( ABS(x(m)-99999.).GT.1.0E-05 ) THEN
                    NOBs(k) = NOBs(k) + 1
                    RO(3,NOBs(k),k) = y(3)/drad
                    RO(2,NOBs(k),k) = (90.-y(2))/drad
                    RO(1,NOBs(k),k) = (6371.2+y(1))/3485.7
                    BRO(NOBs(k),k) = x(m)
                    IFLag(NOBs(k),k) = m + 3
                    SIGma(NOBs(k),k) = SIGma(1,k)
                    ENDIF
                ENDDO
            ENDDO
        ENDIF
60         WRITE (IOUt,'(a)') 'format ' , dform(k)
        WRITE (IOUt,'(a)') 'data file ' , filnam(k)
        WRITE (IOUt,*) NOBs(k) ,                                    &
                    &' observation points and magnetic field data read'
        MOBs = MOBs + NOBs(k)
        ENDDO
    ENDIF
    WRITE (IOUt,*) 'Total observation points' , MOBs
    IF ( MOBs.GT.MAXOB ) THEN
        WRITE (IOUt,*) 'maxob = ' , MAXOB ,                            &
                    &' recompile with larger array    dimensions'
        STOP
    ENDIF
    END
!
!_______________________________________________________________________
    SUBROUTINE ELLIPT(R,Theta,Gr,Gtheta)
    IMPLICIT NONE
    DOUBLE PRECISION a , b , coesq , esq , glam , Gr , Gtheta , pi ,  &
                    & R , Theta
!
!******calls no other routines
!
!  Converts geodetic radius, gr, and colatitude, gtheta, into
!  geocentric, r, theta
!  All angles in radians.
!
    DATA pi/3.14159265/ , a/6378.139/ , b/6356.750/
!
    esq = 1.0 - (b/a)**2
    coesq = (b/a)**2
    glam = pi/2. - Gtheta
    Theta = ATAN(1.0/(coesq*TAN(glam)))
    IF ( Theta.LT.0. ) Theta = Theta + pi
    R = a*SQRT((1-esq*(2-esq)*(SIN(glam))**2)/(1-esq*(SIN(glam))**2))
    END
!
!______________________________________________________________________
!  UNIT 2:  CORE DATA ROUTINES
!_______________________________________________________________________
    SUBROUTINE INCORE
    IMPLICIT NONE
    INTEGER IECho , INMX , INP , IOUt , IPRint , ISTate , l , nflag , &
        & nfound , NIN , nwant1
    DOUBLE PRECISION val
!
!$$$$ calls no other routines
!  The input routine for the downward continuation of geomagnetic
!  field to the surface of the core
!  Reads from the standard input until eof or continue
!  Saves the lines in the character array input for later
!  parsing by getchr and getval.
!  Prints a list of commands for this group of routines
!  Option exists to read all inputs from a specified file
!
    PARAMETER (INMX=200)
    CHARACTER*80 INPut(INMX) , line , filnam
    COMMON /DICTA / INPut
    COMMON /NDICT / IECho , NIN , ISTate(INMX)
    COMMON /IO    / INP , IOUt , IPRint
    DIMENSION val(20)
!
    NIN = 0
    WRITE (IOUt,'(a)') ' ' , '                    =================' ,&
                    &' ' ,                                          &
                &'Enter commands for core field inversion (? for help)'
!
    DO l = NIN + 1 , INMX
        READ (*,'(80a)',END=200) line
        IF ( line(1:4).EQ.'cont' ) GOTO 200
        IF ( line(1:4).EQ.'cfil' ) GOTO 100
        IF ( line(1:1).EQ.'?' ) THEN
        WRITE (IOUt,'(a/a/(2x,a))')                                 &
                            &'Enter commands from the following list:'&
                            & ,                                        &
                    &'   (M) means mandatory information (O) optional'&
                    & , ' ' ,                                          &
                    &'?:        Remind me of the command list again' ,&
            &'continue: Quit reading commands and begin calculations'&
            & , 'cfile filnam: Read all commands from filnam' ,       &
        &'obsdim n: dimension of observation array, 3 for just points'&
        & , '       6 for field values too (M)' ,                      &
        &'e1file file: name of file with evaluation points for forward'&
    & , '       model or observations for inversion' ,              &
        &'e1format a: data format in e1file (M)' ,                     &
        &'e2file file: name of file with evaluation points for forward'&
    & ,                                                             &
    &'       model or observations for inversion at second epoch (O)'&
    & , 'e2format a: data format in e2file (O)' ,                     &
    &'design n:   Construct design matrix (1), read it (0) or read' ,&
    &'     upper triangular form after initial qr (-1) (M)' ,        &
    &'problem n:  Forward (1) or inverse (-1) or fflux (0)' ,        &
    &'     problem (M)' ,                                            &
    &'corefield file: file with core model for forward calculation (O)'
        WRITE (IOUt,'(/(2x,a))')                                    &
    &'shrep ldeg r: degree and radius at whish to evaluate spherical'&
    & , '       harmonic representation (O)' ,                        &
    &'patches 1: Find null flux curves and flux through patches (O)' &
    & , 'epochs n:  Number of time snapshots (M)' ,                   &
    &'rmse noise:  rms msifit level for obsfiles (M)' ,              &
    &'lambda rlam1, rlam2,...:  Lagrange multipliers (M)' ,          &
    &'tnew n: tesselation exists (0), needs creating (1,2), 2 saves a'&
    & , '       file for plotting' ,                                   &
    &'corepoints file: file with tesselation points (M)' ,                &
    &'pointer file: file with pointers for triangles (M)' ,           &
    &'bodydim n: dimension of tesselation points (M)' ,               &
    &'zero n1,n2,...: make the listed points in the tesselation have' &
    & , '      zero radial field' ,                                    &
    &'outfil file: file for output of core model(s)(O)' ,             &
    &'fiter n: no of fr flux iterations allowed' ,                    &
    &'invtype minbrsq (or mingrad or mindiff), type of' ,             &
    &'      regularization' , ' '
        ENDIF
        IF ( line(1:1).NE.' ' ) THEN
        NIN = NIN + 1
        INPut(NIN) = line
        ISTate(NIN) = 0
        ENDIF
    ENDDO
100  NIN = NIN + 1
    INPut(NIN) = line
    ISTate(NIN) = 0
    CALL GETCHR('cfil',filnam,nfound)
    WRITE (IOUt,'(a)') 'Commands read from file ' , filnam
    OPEN (UNIT=10,FILE=filnam)
    DO l = NIN + 1 , INMX
        READ (10,'(80a)',END=200) line
        IF ( line(1:4).EQ.'cont' ) GOTO 200
        IF ( line(1:1).EQ.'?' ) THEN
        WRITE (IOUt,'(a/a/(2x,a))')                                 &
                            &'Enter commands from the following list:'&
                            & ,                                        &
                    &'   (M) means mandatory information (O) optional'&
                    & , ' ' ,                                          &
                    &'?:        Remind me of the command list again' ,&
            &'continue: Quit reading commands and begin calculations'&
            & ,                                                       &
        &'cfil filnam: Read all subsequent commands from file filnam'&
        & ,                                                           &
        &'obsdim n: dimension of observation array, 3 for just points'&
        & , '6 for field values too (M)' ,                             &
        &'e1file file: name of file with evaluation points for forward'&
    & , 'model or observations for inversion' ,                     &
        &'e1format a: data format in e1file (M)' ,                     &
        &'e2file file: name of file with evaluation points for forward'&
    & , 'model or observations for inversion at second epoch (O)' , &
        &'e2format a: data format in e2file (O)' ,                     &
        &'design n:   Construct design matrix (1), read it (0) or read'&
    & , 'upper triangule form after initial qr (-1) (M)' ,          &
        &'problem n:  Forward (1) or inverse (-1) problem (M)' ,       &
    &'corefield file: file with core model for forward calculation (O)'&
    & ,                                                                &
    &'shrep ldeg r: degree and radius at whish to evaluate spherical' &
    & , 'harmonic representation (O)'
        WRITE (IOUt,'(/(2x,a))')                                    &
    &'patches 1: Find null flux curves and flux trhough patches (O)'&
    & , 'epochs n:  Number of time snapshots (M)' ,                  &
    &'rmse noise:  rms msifit level for obsfiles (M)' ,             &
    &'lambda rlam1, rlam2,...:  Lagrange multipliers (M)' ,         &
    &'tnew n: tesselation exists (0), needs creating (1,2), 2 saves a'&
    & , 'file for plotting' ,                                          &
    &'corepoints file: file with tesselation points (M)' ,                &
    &'pointer file: file with pointers for triangles (M)' ,           &
    &'bodydim n: dimension of tesselation points (M)' , ' '
        ENDIF
        IF ( line(1:1).NE.' ' ) THEN
        NIN = NIN + 1
        INPut(NIN) = line
        ISTate(NIN) = 0
        ENDIF
    ENDDO
!
!  Check that all mandatory parameters are available
!
200  nflag = 0
    CALL GETCHR('e1file',filnam,nfound)
!	if(nfound.lt.0)nflag=1
    CALL GETCHR('corepoints',filnam,nfound)
    IF ( nfound.LT.0 ) nflag = 1
    CALL GETCHR('pointer',filnam,nfound)
    IF ( nfound.LT.0 ) nflag = 1
!	call getval('obsdim',val,nwant1,nfound)
!	if(nfound.lt.0)nflag=1
    CALL GETVAL('bodydim',val,nwant1,nfound)
    IF ( nfound.LT.0 ) nflag = 1
!	call getval('design',val,nwant1,nfound)
!	if(nfound.lt.0)nflag=1
    CALL GETVAL('problem',val,nwant1,nfound)
    IF ( nfound.LT.0 ) nflag = 1
    CALL GETVAL('epochs',val,nwant1,nfound)
    IF ( nfound.LT.0 ) nflag = 1
!	call getval('rmse',val,nwant1,nfound)
!	if(nfound.eq.0)nflag=1
!	call getval('lambda',val,nwant1,nfound)
!	if(nfound.eq.0)nflag=1
!	call getchr('e1format',filnam,nfound)
!	if(nfound.lt.0)nflag=1
!
!  If not all there then stop now
    IF ( nflag.EQ.1 ) THEN
        WRITE (IOUt,*) 'Insufficient information supplied.'
        WRITE (IOUt,'(a)')                                             &
    &'Check list of mandatory parameters and your input for omissions o&
    &r errors'
        STOP
    ENDIF
    END
!
!
!=======================================================================
!==UNIT 6:  MISCELLANEOUS ROUTINES======================================
!=======================================================================
    SUBROUTINE GETVAL(Code,Values,Nwant,Nfound)
    IMPLICIT NONE
    
    INTEGER IECho , iex , INMX , INP , IOUt , IPRint , ISTate , k ,   &
        & k1 , k2 , kn , l , lin , m , nb , Nfound , NIN , nval ,   &
        & Nwant
    DOUBLE PRECISION Values
!$$$$ calls no other routines
!
!  Evaluates numbers in the string  char.  nwant  numbers are expected,
!  nfound  are actually found in the character string.  If an error is
!  discovered  nfound=-n,  where  n  is the number of numbers properly
!  decoded.  If there are no numbers after the codeword nfound=0.
    PARAMETER (INMX=200)
    CHARACTER*80 INPut(INMX) , line
    COMMON /DICTA / INPut
    COMMON /NDICT / IECho , NIN , ISTate(INMX)
    COMMON /IO    / INP , IOUt , IPRint
    DIMENSION Values(10)
    CHARACTER bore(2)*1 , form*8 , Code*4 , char*80
!
    DATA bore/'e' , ' '/
!
!e
    DO lin = NIN , 1 , -1
        line = INPut(lin)
        IF ( Code.EQ.line(1:4) ) THEN
        IF ( IECho.GE.1 ) WRITE (IOUt,'(2a)') '>>> ' , line
        ISTate(lin) = ISTate(lin) + 1
        nb = INDEX(line,' ') + 1
        char = line(nb:80)
        kn = 80 - nb + 1
        GOTO 100
        ENDIF
    ENDDO
!  code word not found
    Nfound = -99
    RETURN
!
100  k1 = 1
!  Up to  nwant  numbers are sought.
    DO nval = 1 , Nwant
        DO k = k1 , kn
        IF ( char(k:k).NE.' ' ) GOTO 150
        ENDDO
        Nfound = nval - 1
        RETURN
150     iex = 1
        DO l = k , kn
        IF ( char(l:l).EQ.',' ) GOTO 200
        IF ( char(l:l).EQ.' ' ) THEN
            IF ( char(l-1:l-1).NE.bore(iex) ) GOTO 200
            iex = 2
        ENDIF
        ENDDO
200     m = l - k
!  Augment unpredictable error trapping of some compilers e.g. Sun
        IF ( INDEX('Ee+-.0123456789',char(k:k)).EQ.0 ) GOTO 400
        k2 = l + 1
        WRITE (form,'(2h(f, i3, 3h.0) )') m
        READ (char(k:l-1),form,ERR=400) Values(nval)
        k1 = k2
    ENDDO
    nval = 1 - Nwant
300  Nfound = 1 - nval
    RETURN
!
400  WRITE (IOUt,*) ' ' ,                                              &
                &'***** Unreadable numbers in this input line:' ,   &
                & line
    GOTO 300
    END
!______________________________________________________________
    SUBROUTINE GETCHR(Code,Char,Nfound)
    IMPLICIT NONE
    INTEGER IECho , INMX , INP , IOUt , IPRint , ISTate , k , lin ,   &
        & nb , Nfound , NIN
!$$$$ calls no other routines
    PARAMETER (INMX=200)
    CHARACTER*80 INPut(INMX) , line
    CHARACTER Code*4
    CHARACTER*(*) Char
    COMMON /DICTA / INPut
    COMMON /NDICT / IECho , NIN , ISTate(INMX)
    COMMON /IO    / INP , IOUt , IPRint
!e
    DO lin = NIN , 1 , -1
        line = INPut(lin)
        IF ( Code.EQ.line(1:4) ) THEN
        IF ( IECho.GE.1 ) WRITE (IOUt,'(2a)') '>>> ' , line
        ISTate(lin) = ISTate(lin) + 1
        GOTO 100
        ENDIF
    ENDDO
!  code word not found
    Nfound = -99
    RETURN
!
100  nb = INDEX(line,' ') + 1
    DO k = nb , 80
        IF ( line(k:k).NE.' ' ) THEN
        Char = line(k:80)
        Nfound = 1
        RETURN
        ENDIF
    ENDDO
!
!  Blank field after code word
    Nfound = 0
    END
!______________________________________________________________
    SUBROUTINE QL(Nm,N,A,D,E,Ierr)
    IMPLICIT NONE
    DOUBLE PRECISION A , D , E
    INTEGER Ierr , N , Nm
!$$$$ calls tred2, tql2
!  using  eispack  routines tred2, tql2,  solves the symmetric
!  eigenvalue-eigenvector problem for a real matrix.
!  on input
!  nm   row dimension of the symmetric array  a  in the caller.
!  n    order of the array (<= nm)
!  a    the real symmetric array to be treated
!  e     a working array at least  n  long
!
!  on output
!  d    the array of eigenvalues ascending
!  a    the corresponding array of eigenvectors, with row
!       dimension  nm.  original  a  is overwritten.
!  ierr 0  if all's well.
!       j  if eigenvalue number j and above not found.
!
!
    DIMENSION A(Nm,*) , D(*) , E(*)
!
    CALL TRED2(Nm,N,A,D,E,A)
!
    CALL TQL2(Nm,N,D,E,A,Ierr)
    END
!______________________________________________________________________
    SUBROUTINE TQL2(Nm,N,D,E,Z,Ierr)
    IMPLICIT NONE
    DOUBLE PRECISION b , c , D , E , f , g , h , p , r , s , Z
!$$$$ calls no other routines
!
    INTEGER i , j , k , l , m , N , ii , l1 , Nm , mml , Ierr
    DIMENSION D(N) , E(N) , Z(Nm,N)
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input:
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement;
!
!        n is the order of the matrix;
!
!        d contains the diagonal elements of the input matrix;
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary;
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output:
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1;
!
!        e has been destroyed;
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues;
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     ------------------------------------------------------------------
!
!     :::::::::: machep is a machine dependent parameter specifying
!                the relative precision of floating point arithmetic.
!                machep = 16.0d0**(-13) for long form arithmetic
!                on s360 ::::::::::
    REAL*4 machep
    DATA machep/1.0E-06/
!
    Ierr = 0
    IF ( N.NE.1 ) THEN
!
        DO i = 2 , N
        E(i-1) = E(i)
        ENDDO
!
        f = 0.0D0
        b = 0.0D0
        E(N) = 0.0D0
!
        DO l = 1 , N
        j = 0
        h = machep*(ABS(D(l))+ABS(E(l)))
        IF ( b.LT.h ) b = h
!     :::::::::: look for small sub-diagonal element ::::::::::
        DO m = l , N
            IF ( ABS(E(m)).LE.b ) GOTO 20
!     :::::::::: e(n) is always zero, so there is no exit
!                through the bottom of the loop ::::::::::
        ENDDO
!
20         IF ( m.NE.l ) THEN
30            IF ( j.EQ.30 ) GOTO 100
            j = j + 1
!     :::::::::: form shift ::::::::::
            l1 = l + 1
            g = D(l)
            p = (D(l1)-g)/(2.0D0*E(l))
            r = SQRT(p*p+1.0D0)
            D(l) = E(l)/(p+SIGN(r,p))
            h = g - D(l)
!
            DO i = l1 , N
                D(i) = D(i) - h
            ENDDO
!
            f = f + h
!     :::::::::: ql transformation ::::::::::
            p = D(m)
            c = 1.0D0
            s = 0.0D0
            mml = m - l
!     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
            DO ii = 1 , mml
                i = m - ii
                g = c*E(i)
                h = c*p
                IF ( ABS(p).LT.ABS(E(i)) ) THEN
                    c = p/E(i)
                    r = SQRT(c*c+1.0D0)
                    E(i+1) = s*E(i)*r
                    s = 1.0D0/r
                    c = c*s
                ELSE
                    c = E(i)/p
                    r = SQRT(c*c+1.0D0)
                    E(i+1) = s*p*r
                    s = c/r
                    c = 1.0D0/r
                ENDIF
                p = c*D(i) - s*g
                D(i+1) = h + s*(c*g+s*D(i))
!     :::::::::: form vector ::::::::::
                DO k = 1 , N
                    h = Z(k,i+1)
                    Z(k,i+1) = s*Z(k,i) + c*h
                    Z(k,i) = c*Z(k,i) - s*h
                ENDDO
!
            ENDDO
!
            E(l) = s*p
            D(l) = c*p
            IF ( ABS(E(l)).GT.b ) GOTO 30
        ENDIF
        D(l) = D(l) + f
        ENDDO
!     :::::::::: order eigenvalues and eigenvectors ::::::::::
        DO ii = 2 , N
        i = ii - 1
        k = i
        p = D(i)
!
        DO j = ii , N
            IF ( D(j).LT.p ) THEN
                k = j
                p = D(j)
            ENDIF
        ENDDO
!
        IF ( k.NE.i ) THEN
            D(k) = D(i)
            D(i) = p
!
            DO j = 1 , N
                p = Z(j,i)
                Z(j,i) = Z(j,k)
                Z(j,k) = p
            ENDDO
        ENDIF
!
!
        ENDDO
    ENDIF
    GOTO 99999
!     :::::::::: set error -- no convergence to an
!                eigenvalue after 30 iterations ::::::::::
100  Ierr = l
99999 END
!______________________________________________________________________
    SUBROUTINE TRED2(Nm,N,A,D,E,Z)
    IMPLICIT NONE
    DOUBLE PRECISION A , D , E , f , g , h , hh , scale , Z
!$$$$ calls no other routines
!
    INTEGER i , j , k , l , N , ii , Nm , jp1
    DIMENSION A(Nm,N) , D(N) , E(N) , Z(Nm,N)
!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     on input:
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement;
!
!        n is the order of the matrix;
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output:
!
!        d contains the diagonal elements of the tridiagonal matrix;
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero;
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction;
!
!        a and z may coincide.  if distinct, a is unaltered.
!
!     questions and comments should be directed to b. s. garbow,
!     applied mathematics division, argonne national laboratory
!
!     ------------------------------------------------------------------
!
    DO i = 1 , N
!
        DO j = 1 , i
        Z(i,j) = A(i,j)
        ENDDO
    ENDDO
!
    IF ( N.NE.1 ) THEN
!     :::::::::: for i=n step -1 until 2 do -- ::::::::::
        DO ii = 2 , N
        i = N + 2 - ii
        l = i - 1
        h = 0.0D0
        scale = 0.0D0
        IF ( l.LT.2 ) THEN
            E(i) = Z(i,l)
        ELSE
!     :::::::::: scale row (algol tol then not needed) ::::::::::
            DO k = 1 , l
                scale = scale + ABS(Z(i,k))
            ENDDO
!
            IF ( scale.NE.0.0D0 ) THEN
!
                DO k = 1 , l
                    Z(i,k) = Z(i,k)/scale
                    h = h + Z(i,k)*Z(i,k)
                ENDDO
!
                f = Z(i,l)
                g = -SIGN(SQRT(h),f)
                E(i) = scale*g
                h = h - f*g
                Z(i,l) = f - g
                f = 0.0D0
!
                DO j = 1 , l
                    Z(j,i) = Z(i,j)/h
                    g = 0.0D0
!     :::::::::: form element of a*u ::::::::::
                    DO k = 1 , j
                    g = g + Z(j,k)*Z(i,k)
                    ENDDO
!
                    jp1 = j + 1
                    IF ( l.GE.jp1 ) THEN
!
                    DO k = jp1 , l
                        g = g + Z(k,j)*Z(i,k)
                    ENDDO
                    ENDIF
!     :::::::::: form element of p ::::::::::
                    E(j) = g/h
                    f = f + E(j)*Z(i,j)
                ENDDO
!
                hh = f/(h+h)
!     :::::::::: form reduced a ::::::::::
                DO j = 1 , l
                    f = Z(i,j)
                    g = E(j) - hh*f
                    E(j) = g
!
                    DO k = 1 , j
                    Z(j,k) = Z(j,k) - f*E(k) - g*Z(i,k)
                    ENDDO
                ENDDO
            ELSE
                E(i) = Z(i,l)
            ENDIF
        ENDIF
!
        D(i) = h
        ENDDO
    ENDIF
!
    D(1) = 0.0D0
    E(1) = 0.0D0
!     :::::::::: accumulation of transformation matrices ::::::::::
    DO i = 1 , N
        l = i - 1
        IF ( D(i).NE.0.0D0 ) THEN
!
        DO j = 1 , l
            g = 0.0D0
!
            DO k = 1 , l
                g = g + Z(i,k)*Z(k,j)
            ENDDO
!
            DO k = 1 , l
                Z(k,j) = Z(k,j) - g*Z(k,i)
            ENDDO
        ENDDO
        ENDIF
!
        D(i) = Z(i,i)
        Z(i,i) = 1.0D0
        IF ( l.GE.1 ) THEN
!
        DO j = 1 , l
            Z(i,j) = 0.0D0
            Z(j,i) = 0.0D0
        ENDDO
        ENDIF
!
    ENDDO
!
    END
!______________________________________________________________________
    SUBROUTINE EIGV(F,Ev,R,Tr)
    IMPLICIT NONE
    DOUBLE PRECISION a , a0 , d1 , d2 , d3 , da , eps , Ev , F , fm , &
                    & fm2 , g , g42 , g52 , g62 , R , ss , tol , Tr , x
    DOUBLE PRECISION x2 , y11 , y21 , y22 , y31 , y32 , y33 , yn ,    &
                    & ynor
    INTEGER i , is , k
!  Calculates the eigenvalues , ev(k), and eigenvectors, r(i,k), of a
!  real symmetric 3x3 matrix, with elements f1-f6.
!  tr=trace, f1=m11,f2=m22, f3=m33, f4=m12, f5=m13, f6=m23
    DIMENSION F(6) , Ev(3) , R(3,3) , g(6) , yn(3)
    DATA tol/2.D-14/
    Tr = (F(1)+F(2)+F(3))/3.D0
!***put traceless part into g
    DO i = 1 , 3
        g(i) = F(i) - Tr
        g(i+3) = F(i+3)
    ENDDO
!***cope with zero elements
    DO i = 1 , 6
        IF ( g(i).EQ.0.D0 ) g(i) = tol
    ENDDO
    g62 = g(6)*g(6)
    g52 = g(5)*g(5)
    g42 = g(4)*g(4)
    fm2 = 0.5D0*(g(1)*g(1)+g(2)*g(2)+g(3)*g(3)) + g42 + g52 + g62
    fm = SQRT(fm2)
    a0 = g(1)*g62 + g(2)*g52 + g(3)*g42 - g(1)*g(2)*g(3) - 2.D0*g(4)  &
        & *g(5)*g(6)
    a = a0/(fm*fm2)
    da = ABS(a)
!***iterate for largest root
    x = 1.D0
100  x2 = x*x
    eps = (da+x*(1.D0-x2))/(3.D0*x2-1.D0)
    x = x + eps
    IF ( ABS(eps).GT.tol ) GOTO 100
    IF ( a.GT.0.D0 ) x = -x
    ss = SQRT(DMAX1(1.0+4.0*a/(x**3),tol))
!***construct all eigenvals
    Ev(1) = x*fm
    Ev(2) = -0.5D0*Ev(1)*(1.D0+ss)
    Ev(3) = -0.5D0*Ev(1)*(1.D0-ss)
!***do eignevectors
    DO k = 1 , 3
        d1 = g(1) - Ev(k)
        d2 = g(2) - Ev(k)
        d3 = g(3) - Ev(k)
        y11 = d2*d3 - g62
        y21 = g(5)*g(6) - d3*g(4)
        y31 = g(4)*g(6) - d2*g(5)
        y22 = d1*d3 - g52
        y32 = g(4)*g(5) - d1*g(6)
        y33 = d1*d2 - g42
        yn(1) = y11*y11 + y21*y21 + y31*y31
        yn(2) = y21*y21 + y22*y22 + y32*y32
        yn(3) = y31*y31 + y32*y32 + y33*y33
        is = 1
        IF ( yn(2).GT.yn(1) ) is = 2
        IF ( yn(3).GT.yn(is) ) is = 3
        ynor = SQRT(yn(is))
        IF ( is.EQ.2 ) THEN
        R(1,k) = y21/ynor
        R(2,k) = y22/ynor
        R(3,k) = y32/ynor
        ELSEIF ( is.EQ.3 ) THEN
        R(1,k) = y31/ynor
        R(2,k) = y32/ynor
        R(3,k) = y33/ynor
        ELSE
        R(1,k) = y11/ynor
        R(2,k) = y21/ynor
        R(3,k) = y31/ynor
        ENDIF
    ENDDO
    Ev(1) = Ev(1) + Tr
    Ev(2) = Ev(2) + Tr
    Ev(3) = Ev(3) + Tr
    END
!______________________________________________________________________________________
    SUBROUTINE PATCEQ
    IMPLICIT NONE
    DOUBLE PRECISION EQUIV , FLUxv , GRAdm , P , perc , percmax , val
    INTEGER i , iff , ii , INP , IOUt , IPP , IPRint , IPTs , ISIgn , &
        & ITRis , IV , j , k , LLPt , LPTr , LPTs , LTRis , m1 ,    &
        & m2 , MAXE
    INTEGER MAXP , mnp , N , NEAr , nepoch , nfound , nllpt , nmr ,   &
        & nmrag , nnp , NP , NPAtch , npr , nprag , NPTs , NTRia ,  &
        & NTRis , NVERT , NXYZ , NXYZ2
!
!  Finds equivalent patches for different epochs
!  Assumes fluxpt has been called to give pointer to patches ranked
!  by magnitude of the flux through them
!  Starts with largest patches, checks sign and points in common,
!  works down in size.
!  More than one patch may be equivalent to another if less than 60%
!  of points are in common.  If no equivalent patch exists then add a
!  null patch to equivalence
!  In difficult cases option exists to specify patch equivalents in
!  run file and run one iteration at a time
!  llpt is reordered so that llpt(i,1) is a patch equivalent to llpt(i,2)
!
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ,NVERT=10)
    PARAMETER (MAXP=30,MAXE=200)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    COMMON /PATCHY/ NP(MAXE) , NPTs(MAXP,MAXE) , ISIgn(MAXP,MAXE) ,   &
                & IPTs(MAXP,MAXE) , LPTs(2*NXYZ,MAXE) , NTRis(NXYZ) &
                & , ITRis(NXYZ+1) , LTRis(2*(2*NXYZ-4)) ,           &
                & LPTr(2*NXYZ-4,2) , LLPt(MAXP,MAXE) , NPAtch ,     &
                & IPP(NXYZ)
    COMMON /PATCHG/ FLUxv(MAXP,MAXE) , GRAdm(MAXP,NXYZ,MAXE)
    DIMENSION val(30) , nllpt(MAXP,MAXE) , nprag(MAXP,MAXE) ,         &
            & nmrag(MAXP,MAXE) , npr(2) , nmr(2)
    DATA iff/0/
!
    nepoch = 2
    CALL GETVAL('equiv',val,30,nfound)
    IF ( nfound.GT.0 .AND. iff.EQ.0 ) THEN
        NPAtch = nfound/2
        WRITE (IOUt,*) NPAtch , ' equivalences read'
        DO i = 1 , nfound/2
        LLPt(i,1) = NINT(val(2*i-1))
        LLPt(i,2) = NINT(val(2*i))
        ENDDO
        iff = 1
        RETURN
    ELSE
        nnp = MIN(NP(1),NP(2))
        mnp = MAX(NP(1),NP(2))
        k = 0
        DO j = 1 , 2
        npr(j) = 0
        nmr(j) = 0
        ENDDO
!
!  Starting with biggest patches of positive and negative flux,
!  sort their arrays of points by index. If less
!  than 60% are in common save these patches for reordering.  Otherwise keep these
!  patches as equivalent. Goto next patch etc.
!
        DO i = 1 , nnp/2
!  Start with negative patch
        CALL ISORT(NPTs(LLPt(i,1),1),LPTs(IPTs(LLPt(i,1),1),1))
        CALL ISORT(NPTs(LLPt(i,2),2),LPTs(IPTs(LLPt(i,2),2),2))
        perc = EQUIV(LPTs(IPTs(LLPt(i,1),1),1),NPTs(LLPt(i,1),1),   &
                & LPTs(IPTs(LLPt(i,2),2),2),NPTs(LLPt(i,2),2))
        IF ( perc.LT.60.0 ) THEN
!  Add to record of unmatched patches
            DO j = 1 , 2
                IF ( FLUxv(LLPt(i,j),j).GE.0 ) THEN
                    npr(j) = npr(j) + 1
                    nprag(npr(j),j) = LLPt(i,j)
                ELSEIF ( FLUxv(LLPt(i,j),j).LE.0 ) THEN
                    nmr(j) = nmr(j) + 1
                    nmrag(nmr(j),j) = LLPt(i,j)
                ENDIF
            ENDDO
        ELSE
!  Otherwise keep the match gven by flux ordering
            k = k + 1
            nllpt(k,1) = LLPt(i,1)
            nllpt(k,2) = LLPt(i,2)
            WRITE (*,'(a,i2,a,i2,a,f6.2)') ' Epoch 1, patch ' ,      &
                & LLPt(i,1) , ' equivalenced to Epoch 2, patch ' ,  &
                & LLPt(i,2) ,                                       &
                & ' percentage of points in common is ' , perc
        ENDIF
!  Goto positive patch
        CALL ISORT(NPTs(LLPt(NP(1)-i+1,1),1),                       &
                    & LPTs(IPTs(LLPt(NP(1)-i+1,1),1),1))
        CALL ISORT(NPTs(LLPt(NP(2)-i+1,2),2),                       &
                    & LPTs(IPTs(LLPt(NP(2)-i+1,2),2),2))
        perc = EQUIV(LPTs(IPTs(LLPt(NP(1)-i+1,1),1),1),             &
                & NPTs(LLPt(NP(1)-i+1,1),1),                           &
                & LPTs(IPTs(LLPt(NP(2)-i+1,2),2),2),                   &
                & NPTs(LLPt(NP(2)-i+1,2),2))
        IF ( perc.LT.60.0 ) THEN
!  Add to record of unmatched patches
            DO j = 1 , 2
                IF ( FLUxv(LLPt(NP(j)-i+1,j),j).GE.0. ) THEN
                    npr(j) = npr(j) + 1
                    nprag(npr(j),j) = LLPt(NP(j)-i+1,j)
                ELSEIF ( FLUxv(LLPt(NP(j)-i+1,j),j).LE.0. ) THEN
                    nmr(j) = nmr(j) + 1
                    nmrag(nmr(j),j) = LLPt(NP(j)-i+1,j)
                ENDIF
            ENDDO
        ELSE
!  Otherwise keep the match gven by flux ordering
            k = k + 1
            nllpt(k,1) = LLPt(NP(1)-i+1,1)
            nllpt(k,2) = LLPt(NP(2)-i+1,2)
            WRITE (*,'(a,i2,a,i2,a,f6.2)') ' Epoch 1, patch ' ,      &
                & LLPt(NP(1)-i+1,1) ,                               &
                    &' equivalenced to Epoch 2, patch ' ,              &
                & LLPt(NP(2)-i+1,2) ,                               &
                    &' percentage of points in common is ' , perc
        ENDIF
        ENDDO
        IF ( MOD(nnp,2).EQ.1 ) THEN
        i = nnp/2 + 1
!  Do the last patch
        CALL ISORT(NPTs(LLPt(i,1),1),LPTs(IPTs(LLPt(i,1),1),1))
        CALL ISORT(NPTs(LLPt(i,2),2),LPTs(IPTs(LLPt(i,2),2),2))
        perc = EQUIV(LPTs(IPTs(LLPt(i,1),1),1),NPTs(LLPt(i,1),1),   &
                & LPTs(IPTs(LLPt(i,2),2),2),NPTs(LLPt(i,2),2))
        IF ( perc.LT.60.0 ) THEN
!  Add to record of unmatched patches
            DO j = 1 , 2
                IF ( FLUxv(LLPt(i,j),j).GE.0 ) THEN
                    npr(j) = npr(j) + 1
                    nprag(npr(j),j) = LLPt(i,j)
                ELSEIF ( FLUxv(LLPt(i,j),j).LE.0 ) THEN
                    nmr(j) = nmr(j) + 1
                    nmrag(nmr(j),j) = LLPt(i,j)
                ENDIF
            ENDDO
        ELSE
!  Otherwise keep the match gven by flux ordering
            k = k + 1
            nllpt(k,1) = LLPt(i,1)
            nllpt(k,2) = LLPt(i,2)
            WRITE (*,'(a,i2,a,i2,a,f6.2)') ' Epoch 1, patch ' ,      &
                & LLPt(i,1) , ' equivalenced to Epoch 2, patch ' ,  &
                & LLPt(i,2) ,                                       &
                & ' percentage of points in common is ' , perc
        ENDIF
        ENDIF
!  If more patches at 1 epoch than other add these to appropriate rag end list
        IF ( NP(1).GT.NP(2) ) THEN
        m1 = 1
        ELSE
        m1 = 2
        ENDIF
        DO i = nnp/2 + MOD(nnp,2) + 1 , mnp - nnp/2
        CALL ISORT(NPTs(LLPt(i,m1),m1),LPTs(IPTs(LLPt(i,m1),m1),m1))
!  Add to record of unmatched patches
        IF ( FLUxv(LLPt(i,m1),m1).GE.0 ) THEN
            npr(m1) = npr(m1) + 1
            nprag(npr(m1),m1) = LLPt(i,m1)
        ELSEIF ( FLUxv(LLPt(i,m1),m1).LE.0 ) THEN
            nmr(m1) = nmr(m1) + 1
            nmrag(nmr(m1),m1) = LLPt(i,m1)
        ENDIF
        ENDDO
        WRITE (*,*) npr(1) , npr(2) , nmr(1) , nmr(2)
        IF ( npr(1).GT.0 .OR. nmr(1).GT.0 .OR. npr(2).GT.0 .OR. nmr(2) &
        & .GT.0 ) THEN
!
!  Sort out rag ends by equivalencing all possible patches of same sign
!  to find those with most points in common
!  Find epoch with most patches and use this for reference
        WRITE (*,'(a)') ' Attempting to match +ve patches'
        IF ( npr(1).GT.npr(2) ) THEN
            WRITE (*,*) 'Epoch 1 has more patches'
            m1 = 1
            m2 = 2
        ELSEIF ( npr(2).GT.npr(1) ) THEN
            WRITE (*,*) 'Epoch 2 has more patches'
            m1 = 2
            m2 = 1
        ELSE
            WRITE (*,*) 'Epochs 1 and 2 have same no. of patches'
            m1 = 1
            m2 = 2
        ENDIF
        DO i = 1 , npr(m1)
            percmax = 0.
            k = k + 1
            nllpt(k,m1) = nprag(i,m1)
            DO j = 1 , npr(m2)
                perc = EQUIV(LPTs(IPTs(nprag(i,m1),m1),m1),           &
                    & NPTs(nprag(i,m1),m1),                          &
                    & LPTs(IPTs(nprag(j,m2),m2),m2),                 &
                    & NPTs(nprag(j,m2),m2))
                WRITE (*,'(a,i2,a,i2,a,f6.2,a)') ' Patches ' ,        &
                    & nprag(i,m1) , ' and ' , nprag(j,m2) ,          &
                    & ' have ' , perc ,                              &
                    &' percent of points in common'
!  Save match if better than before
                IF ( perc.GE.percmax ) THEN
                    percmax = perc
                    nllpt(k,m2) = nprag(j,m2)
                ENDIF
            ENDDO
!  If no match equivalence to dummy patch no. 99
            IF ( percmax.LT.1. ) nllpt(k,m2) = 99
            WRITE (*,'(a,i2,a,i2,a,i2,a,i2,a,f6.2)') ' Epoch ' , m1 ,&
                    &' patch ' , nllpt(k,m1) ,                         &
                    &' equivalenced to Epoch ' , m2 , ' patch ' ,      &
                & nllpt(k,m2) ,                                     &
                    &' percentage of points in common is ' , percmax
        ENDDO
!  Check all npr(m2) patches have been included
        DO j = 1 , npr(m2)
            ii = 0
            DO i = 1 , k
                IF ( nllpt(i,m2).EQ.nprag(j,m2) ) ii = 1
            ENDDO
!  Add another patch if not there
            IF ( ii.EQ.0 ) THEN
                k = k + 1
                nllpt(k,m2) = nprag(j,m2)
                nllpt(k,m1) = 99
            ENDIF
        ENDDO
        WRITE (*,'(a)') ' Attempting to match -ve patches'
        IF ( nmr(1).GT.nmr(2) ) THEN
            WRITE (*,*) 'Epoch 1 has more patches'
            m1 = 1
            m2 = 2
        ELSEIF ( nmr(2).GT.nmr(1) ) THEN
            WRITE (*,*) 'Epoch 2 has more patches'
            m1 = 2
            m2 = 1
        ELSE
            WRITE (*,*) 'Epochs 1 and 2 have same no. of patches'
            m1 = 1
            m2 = 2
        ENDIF
        DO i = 1 , nmr(m1)
            percmax = 0.
            k = k + 1
            nllpt(k,m1) = nmrag(i,m1)
            DO j = 1 , nmr(m2)
                perc = EQUIV(LPTs(IPTs(nmrag(i,m1),m1),m1),           &
                    & NPTs(nmrag(i,m1),m1),                          &
                    & LPTs(IPTs(nmrag(j,m2),m2),m2),                 &
                    & NPTs(nmrag(j,m2),m2))
                WRITE (*,'(a,i2,a,i2,a,f6.2,a)') ' Patches ' ,        &
                    & nmrag(i,m1) , ' and ' , nmrag(j,m2) ,          &
                    & ' have ' , perc ,                              &
                    &' percent of points in common'
!  Save match if better than before
                IF ( perc.GE.percmax ) THEN
                    percmax = perc
                    nllpt(k,m2) = nmrag(j,m2)
                ENDIF
            ENDDO
!  If no match equivalence to dummy patch no. 99
            IF ( percmax.LT.1.0 ) nllpt(k,m2) = 99
            WRITE (*,'(a,i2,a,i2,a,i2,a,i2,a,f6.2)') ' Epoch ' , m1 ,&
                    &' patch ' , nllpt(k,m1) ,                         &
                    &' equivalenced to Epoch ' , m2 , ' patch ' ,      &
                & nllpt(k,m2) ,                                     &
                    &' percentage of points in common is ' , percmax
        ENDDO
!  Check all nmr(m2) patches have been included
        DO j = 1 , nmr(m2)
            ii = 0
            DO i = 1 , k
                IF ( nllpt(i,m2).EQ.nmrag(j,m2) ) ii = 1
            ENDDO
!  Add another patch if not there
            IF ( ii.EQ.0 ) THEN
                k = k + 1
                nllpt(k,m2) = nmrag(j,m2)
                nllpt(k,m1) = 99
            ENDIF
        ENDDO
        ENDIF
        NPAtch = k
        WRITE (*,*) 'npatch = ' , NPAtch
!  copy nllpt into llpt
        DO i = 1 , NPAtch
        DO j = 1 , nepoch
            LLPt(i,j) = nllpt(i,j)
        ENDDO
        ENDDO
    ENDIF
    END
!______________________________________________________________________________________
    FUNCTION EQUIV(List1,N1,List2,N2)
    IMPLICIT NONE
    DOUBLE PRECISION count , EQUIV
    INTEGER i1 , i2 , List1 , List2 , N1 , N2
    DIMENSION List1(N1) , List2(N2)
!  Counts how many points list1 and list2 have in common and returns the
!  value 200*npts/(n1+n2).  Assumes list1 and list2 are sorted in ascending order
    count = 0.
    i1 = 1
    i2 = 1
100  IF ( List1(i1).EQ.List2(i2) ) THEN
        count = count + 1.
        i1 = i1 + 1
        i2 = i2 + 1
        IF ( i1.LE.N1 .AND. i2.LE.N2 ) GOTO 100
        EQUIV = 200.*count/FLOAT(N1+N2)
    ELSEIF ( List1(i1).GT.List2(i2) ) THEN
        i2 = i2 + 1
        IF ( i2.LE.N2 ) GOTO 100
        EQUIV = 200.*count/FLOAT(N1+N2)
    ELSE
        i1 = i1 + 1
        IF ( i1.LE.N1 ) GOTO 100
        EQUIV = 200.*count/FLOAT(N1+N2)
    ENDIF
    END
!______________________________________________________________________________________
!
!
!--------------------------------------------------------------
    SUBROUTINE PMAKER(Kk,I)
    IMPLICIT NONE
    DOUBLE PRECISION aa , ajac , b , bjac , br , c , FLUxv , g ,      &
                    & GRAdm , P , PAREA , rg , rp , rpa , ua , ub ,    &
                    & uc , V , vc , XYW
    INTEGER I , INP , IOUt , IPP , IPRint , IPTs , ISIgn , ISTack ,   &
        & ITRis , IV , j , k , Kk , l , LLPt , LPTr , LPTs , LTRis ,&
        & m , MAXE
    INTEGER MAXP , N , NEAr , NP , NPAtch , NPTs , NQUad , NTRia ,    &
        & NTRis , NVERT , NXYZ , NXYZ2
!
!  Makes a patch on the core out of points in patch kk at epoch i
!  Computes the gradient matrix for this patch
!  Note that this patch is not distinct in sign from
!  its surroundings, and may even be mixed in sign
!
!  Assumes patch  and flux have already been called to define the
!  variables in /patchy/ for epoch 2.
!
!
!  Outline:
!
!  Initialize the flux vector and flux gradient matrix to zero
!  Loop in triangles
!    for each triangle, compute rotation and gnomonic projection
!    loop in integration sample points
!       compute preimage of the sample point on the core
!       compute the jacobian for the gnomonic at that sample point
!       compute field value at that point
!       add the value, weighted by the int. wt. and the jacobians,
!           to the positive patch or negative patch, according
!           to its sign.
!       add interp. wt*int. wt.* jacobian to the appropriate elements
!           of the grad. matrix (according to signs of br and vertex)
!    next integration sample point
!  next triangle.
!
    PARAMETER (NXYZ=20000,NXYZ2=2*NXYZ,MAXP=30,NVERT=10,MAXE=200)
    DIMENSION rp(3) , g(3) , vc(4) , aa(3) , b(3) , c(3) , ua(2) ,    &
            & ub(2) , uc(2) , rg(2) , rpa(3)
    COMMON /PATCHY/ NP(MAXE) , NPTs(MAXP,MAXE) , ISIgn(MAXP,MAXE) ,   &
                & IPTs(MAXP,MAXE) , LPTs(2*NXYZ,MAXE) , NTRis(NXYZ) &
                & , ITRis(NXYZ+1) , LTRis(2*(2*NXYZ-4)) ,           &
                & LPTr(2*NXYZ-4,2) , LLPt(MAXP,MAXE) , NPAtch ,     &
                & IPP(NXYZ)
    COMMON /PATCHG/ FLUxv(MAXP,MAXE) , GRAdm(MAXP,NXYZ,MAXE)
!
!  see subroutine patch for definitions of the variables in /patchy/
!
    COMMON /TESSEL/ N , NTRia , P(4,NXYZ) , IV(3,NXYZ2) ,             &
                & NEAr(NVERT,NXYZ)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /RULE  / XYW(3,100) , NQUad
    COMMON /CIRCLE/ V(4,NXYZ2) , ISTack(NXYZ2)
!
!  initialize flux vector and gradient matrix.
!******
    WRITE (IOUt,*) ' arrived in pmaker '
    FLUxv(MAXP,1) = 0.
    DO j = 1 , N
        GRAdm(MAXP,j,I) = 0.
    ENDDO
!
!  Loop in triangles associated with the patch
    DO m = 1 , NTRis(Kk)
        j = LTRis(ITRis(Kk)+m-1)
!
!
!  compute centroid vc of current triangle
        CALL CENTR(P(1,IV(1,j)),P(1,IV(2,j)),P(1,IV(3,j)),vc)
!
!  rotate so vertices are in coord system with vc as pole, p(1,iv1) at
!  zero longitude
        CALL GPOLE(1,vc,P(1,IV(1,j)),P(1,IV(1,j)),aa)
        CALL GPOLE(1,vc,P(1,IV(1,j)),P(1,IV(2,j)),b)
        CALL GPOLE(1,vc,P(1,IV(1,j)),P(1,IV(3,j)),c)
!
!  find gnomonic projection of vertices and jacobian of transformation
        CALL GNOMON(1,ua,aa)
        CALL GNOMON(1,ub,b)
        CALL GNOMON(1,uc,c)
        ajac = PAREA(ua,ub,uc)*2.
!*******
!      write(iout,*)' entering integration loop '
!
!  Loop in integration sample points.
        DO k = 1 , NQUad
!
!  find points on gnomonic surface corresponding to integration points
        rg(1) = ua(1) + (ub(1)-ua(1))*XYW(1,k) + (uc(1)-ua(1))      &
                & *XYW(2,k)
        rg(2) = ua(2) + (ub(2)-ua(2))*XYW(1,k) + (uc(2)-ua(2))      &
                & *XYW(2,k)
!
!  Find Jacobian for gnomonic to spherical projection
        bjac = SQRT((1+rg(1)**2+rg(2)**2)**(-3))
!
!  Find interpolating vector g(.)---linear interpolation in gnomonic
!  plane.
        g(1) = 1. - XYW(1,k) - XYW(2,k)
        g(2) = XYW(1,k)
        g(3) = XYW(2,k)
!
!  find the point rp on the core corresponding to rg.
        CALL GNOMON(-1,rg,rpa)
        CALL GPOLE(-1,vc,P(1,IV(1,j)),rp,rpa)
!
!  Find field value and its sign.
        br = 0.
        DO l = 1 , 3
            br = br + g(l)*P(4,IV(l,j))
        ENDDO
!
!  Update the flux vector fluxv and the gradient matrix gradm.
        FLUxv(MAXP,I) = FLUxv(MAXP,I) + ajac*bjac*XYW(3,k)*br
        DO l = 1 , 3
            GRAdm(MAXP,IV(l,j),I) = GRAdm(MAXP,IV(l,j),I) + g(l)     &
                                    & *ajac*bjac*XYW(3,k)
        ENDDO
!
!  Next integration sample point.
        ENDDO
!
!  Next triangle.
    ENDDO
    END
!_______________________________________________________________________
    SUBROUTINE SYNOPS(New,N,P,Ntria,Iv)
    IMPLICIT NONE
    DOUBLE PRECISION AMAx , AMIn , AREA , DMAx , DMIn , DOTTY , P ,   &
                    & p12 , p13 , p23 , sarea , x
    INTEGER i , i1 , i2 , i3 , INP , IOUt , IPRint , it , Iv , N ,    &
        & New , Ntria
!$$$$ calls dotty
!  prints number of triangles in tesselation and outputs list of
!  triangles to fort.1 if new>=2
!
    DIMENSION P(4,*) , Iv(3,Ntria)
    COMMON /IO    / INP , IOUt , IPRint
    COMMON /BOUND / DMIn , DMAx , AMIn , AMAx
    sarea = 0.
    IF ( New.EQ.0 ) THEN
        AMIn = 12.
        AMAx = 0.
        DO i = 1 , Ntria
        i1 = Iv(1,i)
        i2 = Iv(2,i)
        i3 = Iv(3,i)
!	print *, i1,i2,i3
        x = AREA(P(1,i1),P(1,i2),P(1,i3))
        sarea = sarea + x
        IF ( x.LT.AMIn ) AMIn = x
        IF ( x.GT.AMAx ) AMAx = x
        p12 = DOTTY(P(1,i1),P(1,i2))
        DMAx = MAX(p12,DMAx)
        DMIn = MIN(p12,DMIn)
        p13 = DOTTY(P(1,i1),P(1,i3))
        DMAx = MAX(p13,DMAx)
        DMIn = MIN(DMIn,p13)
        p23 = DOTTY(P(1,i2),P(1,i3))
        DMAx = MAX(p23,DMAx)
        DMIn = MIN(DMIn,p23)
        ENDDO
    ELSEIF ( New.EQ.2 ) THEN
        OPEN (UNIT=4,FILE='fort.4')
        OPEN (UNIT=5,FILE='fort.5')
        DO i = 1 , Ntria
        i1 = Iv(1,i)
        i2 = Iv(2,i)
        i3 = Iv(3,i)
        x = AREA(P(1,i1),P(1,i2),P(1,i3))
        WRITE (4,*) x
        sarea = sarea + x
        IF ( x.LT.AMIn ) AMIn = x
        IF ( x.GT.AMAx ) AMAx = x
        p12 = DOTTY(P(1,i1),P(1,i2))
        DMAx = MAX(p12,DMAx)
        DMIn = MIN(p12,DMIn)
        p13 = DOTTY(P(1,i1),P(1,i3))
        DMAx = MAX(p13,DMAx)
        DMIn = MIN(DMIn,p13)
        p23 = DOTTY(P(1,i2),P(1,i3))
        DMAx = MAX(p23,DMAx)
        DMIn = MIN(DMIn,p23)
        WRITE (5,*) p12 , p13 , p23
        ENDDO
    ENDIF
!
!  Describe triangle model
    WRITE (IOUt,'(//a,i6)') ' Number of triangles in body ' , Ntria
!
    WRITE (IOUt,'(/a,2f8.2,a)') ' Triangle sides range between ' ,    &
                            & 57.296*ACOS(DMAx) , 57.296*ACOS(DMIn) &
                            & , ' degrees'
    WRITE (IOUt,'(/a,f8.4,a,f8.4)') ' Triangle areas range between ' ,&
                                & AMAx , ' and ' , AMIn
!
    IF ( New.LE.1 ) RETURN
!
!  For spherical triangles
!  Output to a file every  triangle in the list to fort.2
    OPEN (UNIT=11,FILE='fort.2')
    DO it = 1 , Ntria
        i1 = Iv(1,it)
        i2 = Iv(2,it)
        i3 = Iv(3,it)
        WRITE (11,'(3p,9f8.0,a)') P(1,i1) , P(2,i1) , P(3,i1) , P(1,i2)&
                                & , P(2,i2) , P(3,i2) , P(1,i3) ,      &
                                & P(2,i3) , P(3,i3)
    ENDDO
!
    CLOSE (UNIT=11)
    CLOSE (UNIT=4)
    CLOSE (UNIT=5)
    WRITE (IOUt,'(/a)') ' Spherical triangles written to file fort.2'
    WRITE (IOUt,'(/a)') ' Triangle areas written to file fort.4'
    WRITE (IOUt,'(/a)') ' Triangle angles written to file fort.5'
    END

!_________________________________________________________________________
!

    SUBROUTINE LLTOXYZ(x)
        IMPLICIT NONE

    ! Read in lat, lon, dummy var    

        double precision :: colat, lon
        double precision :: theta, toRad, ct, st, cp, sp, phi
        double precision :: x(3)

        data toRad/0.01745329/

    !   Store lat, lon

        colat=90 - x(1)
        lon=x(2)

        theta = toRad * colat
        ct = cos(theta)
        st = sin(theta)

        phi = toRad * lon
        cp = cos(phi)
        sp = sin(phi)
     
    !   Compute x, y, z 
        x(1) = (st*cp)
        x(2) = (st*sp)
        x(3) = ct

    END
!_________________________________________________________________


