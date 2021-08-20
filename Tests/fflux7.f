	program drive7
	implicit double precision (a-h,o-z)
c reads position cooordinates on the core, creates a tesselation as
c nearly equiangular as possible using dorito, and finds the nearest
c neighbours of every point
c Reads B field measurments on surface, inverts for field values on
c core surface with integrated br squared regularization.
c  Joint inversion for 2 epochs is possible, to minimize difference
c  between the two models
c
c******calls sbody, nneigh, triq, obs, bfield, gram, qr
c
c	nxyz:    maximum no. of points sampled on the core
c	nxyz2	 =2*nxyz
c	nvert	 maximum no. of nearest neighbours  allowed for any
c                point on the core
c	maxob    maximum no. of observations at surface
c
c	in /tessel/
c	
c	n:       no. of points on the core
c	ntria:	no of triangles on core =2*n-4
c	p(i,j):  i=1,3 gives ith coordinate of jth point
c		 p(4,j) gives br at jth point
c	iv(i,k): point number of ith vertex of k-th triangle
c	near(i,j): i=1,nvert gives nearest neighbours of j-th point in
c		   p, ordered sequentially around p.
c
c	in /observ/
c	mobs:    total no. of observation points at Earth surface
c	ro(i,j,k): i=1,3 position of jth observation point, at kth
c		   epoch (r,theta,phi)
c	bro(J,k):  jth observation of a component of the magnetic field 
c                  point, at kth epoch 
c	iflag(j,k): flag identifying component type for jth observation
c	          at kth epoch
c	sigma(j,k):   rms misfit level for datum j at epoch k
c	nepoch:	  number of epochs
c	nobs(k):  number of observation sites at epoch k
c	lobs(k):  number of components at epoch k
c
c	in /calcul/
c	
c	bc(j):  computed radial field at jth core point
c	bs(j,k):  jth computed field component at epoch k
c	a(i,j):  matrix for forward calculations, i=1,mobs,j=1,nn;
c		 also used to tore upper triangle in lowermost nxyz
c		 rows
c
c	in /rule/ 
c	integration rule for right triangle at origin
c	nquad:    no. of points to evaluate the function
c	xyw(i,j): i=1,2 gives x and y coordinate of point for function
c		  evaluation and i=3 gives weight of that point for
c		  integration rule.		
c
c	in /patchy/
c	np(i): 	  total no. of patches of one sign for epoch i
c	npts(j,k):  #vertices (points) in patch j, for epoch k
c	isign(j,k)  sign of patch j, for epoch k
c	ipts(j):  starting index in lpts of the list of points in
c		  patch j
c  	lpts(ipts(j)+k-1):    kth vertex in patch j
c	ntris(j):             # triangles in patch j
c  	itris(j):             starting index in ltris of the list of 
c                        triangles in patch j
c   	ltris(itris(j)+k-1):  kth triangle in patch j
c	lptr(j,k):	k=1:  positive patch (if any)associated with
c		  triangle j.  k=2:negative patch (if any)
c
c	in /patchg/
c	fluxv(j,k):  flux through patch j, for epoch k
c	gradm(j,k,i): gradient of flux through patch j due to change in
c		field at vertex k, for epoch i
c
c	
c	in /circle/
c
c	v(i,k):  i=1,4 contains coordinates of circumcenter
c	 	 and "radius" of circle circumscribing triangle k in
c		 list iv
c	istack:  is a work array of dimension at least nxyz2
c
c
c	in /shr/
c	
c	glm(l+1,m+1), hlm(l+1,m+1):	partially normalized Schmidt 
c	coeffs of degree l, order m
c	plm(l+1,m+1):  Legendre polynomial of degree l, order m
c
c
	parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,nvert=10)
	parameter (maxob=7000,maxp=30,lmax=101,maxe=2)
	common /io/inp,iout,iprint
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
	common /observ/ro(3,maxob,maxe),bro(maxob,maxe),sigma(maxob,
     $           maxe),iflag(maxob,maxe),mobs,nobs(maxe),nepoch,
     $           lobs(maxe)
	common /calcul/bs(maxob,maxe),a(maxob,nxyz),bc(nxyz2),
     $    reg(nxyz,nxyz)
	common /rule/xyw(3,100),nquad
        common/patchy/np(maxe),npts(maxp,maxe),isign(maxp,maxe),
     $   ipts(maxp,maxe),
     $   lpts(2*nxyz,maxe),ntris(nxyz),itris(nxyz+1),
     $   ltris(2*(2*nxyz-4)),lptr(2*nxyz-4,2),
     $   llpt(maxp,maxe),npatch
      common/patchg/  fluxv(maxp,maxe),gradm(maxp,nxyz,maxe)
	common /shr/glm(lmax,lmax),hlm(lmax,lmax),plm(lmax,lmax),
     $        dlm(lmax,lmax)
	common /circle/v(4,nxyz2),istack(nxyz2)
	common /bound/dmin,dmax,amin,amax
	character*80 filnam
	dimension val(10)
c      data ((iv(i,j),i=1,3),j=1,8)/1,2,3,1,3,4,1,4,5,1,2,5,
c     $ 2,3,6,3,4,6,4,5,6,2,5,6/,
c     $ ((p(i,j),i=1,3),j=1,6)/0,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1/
	data inp,iout,iprint/5,6,0/,dmin,dmax/1.0,-1.0/,i1,i2/-1,0/
	data rc/3485/
c
c  Get parameters to run code
	call incore
c
c
c  Get required core points and create a tesselation using dorito
c
	call getval('bodydim',val,1,nfound)
	kdim=nint(val(1))
c
	call sbody(kdim,nxyz,n,p,iv,ntria)
	write(iout,*)'Returned from sbody'
	call flush(iout)
c
c  Find the nearest neighbours of every point on the core and save in 
c  near. This gives support region for basis functions.
c
	call nneigh
		write(iout,*)'Returned from nneigh'
	call flush(iout)
c
c  Set up 7 point Stroud integration rule for forward calculations
c	
	call triq(i1,i2)
		write(iout,*)'Returned from triq'
	call flush(iout)
c
c  Read observation points 
c  kdim=3 for just points kdim=6 for associated field observations too.
c
	call getval('obsdim',val,1,nfound)
	kdim=nint(val(1))
	write(iout,*)'Observation dimension = ',kdim
	call obs(kdim)
c
c  Compute field values at observation points
c  First construct the design matrix
c
c
	call getval('design',val,1,nfound)
	ig=nint(val(1))
	if(ig.eq.1)then
	  call gram
c
c  write down design matrix for possible future use
c
	write(iout,*)'design matrix calculation complete'
	call flush(iout)
	if(nepoch.eq.1)lobs(2)=0 
	  open (unit=20,file='designmatrix',form='unformatted')
	  write(iout,*) 'Lobs(1,2)= ',lobs(1),lobs(2)
	  write(20)lobs(1),lobs(2)
	  do 1 j=1,n
1	  write(20)(a(i,j),i=1,lobs(1))
	if(nepoch.eq.2)then
	do 2 j=n+1,2*n
2        write(20)(a(i,j),i=lobs(1)+2,lobs(1)+lobs(2)+1)
	endif
	  close(unit=20)
	  write(iout,*)'design matrix written to disk'
	else if(ig.eq.0)then
	  open (unit=20,file='designmatrix',form='unformatted')
	  read(20)lobs(1),lobs(2)
	do 3 j=1,n
3	  read(20)(a(i,j),i=1,lobs(1))
	if(nepoch.eq.2)then
	do 4 j=n+1,2*n
4	  read(20)(a(i,j),i=lobs(1)+2,lobs(1)+lobs(2)+1)
	  do 5 j=1,n
	  do 6 i=1,lobs(1)+1
6	  a(i,j+n)=0.0
	  do 7 i=lobs(1)+1,lobs(2)+lobs(1)+2
7       a(i,j)=0.0
	  a(lobs(2)+lobs(1)+2,j)=0.
5	  continue
	endif
	  close(unit=20)
	  write(iout,*)'design matrix read in'
	  write(iout,*)lobs(1),lobs(2)
	else if (ig.ne.-1.and.ig.ne.-2)then
	  write(iout,*)'check input file, (1,0,-1) only valid options'
	endif
	call flush(iout)
	call getval('problem',val,1,nfound)
	if=nint(val(1))
	if (if.eq.1.or.if.eq.0)then
c
c  Forward problem
c   Read in core model
	call getchr('corefile',filnam,nfound)
	open (unit=20,file=filnam)
	do 11 i=1,n
	  read(20,*)x,y,z,bc(i)
c	  read(20,*)x,y,z,bc(i),bc(n+i)
11	continue
	close(unit=20)
	   if(if.eq.1) then
           call bfield
c
c  write out resulting field values at surface
	   do 10 k=1,nepoch
	   write(iout,*)'Epoch ',k
	   do 10 i=1,nobs(k)
	     write(iout,*)(ro(j,i,k),j=1,3),bs(i,k)
10	   continue
	   call flush(iout)
c
c Compute residual and rms misfit
	   do 19 k=1,nepoch
	   write(iout,*)'Epoch ',k
	   rms=0.
	   do 9 i=1,nobs(k)
	     r=-bs(i,k)+bro(i,k)
	     write(iout,*)iflag(i,k),(ro(j,i,k),j=1,3),r
8		rms=rms+ r**2
9	   continue
  	   write(iout,*)'Misfit for epoch ',k, ' is ',
     $  sqrt(rms/float(lobs(k)))
	   call flush(iout)
19	   continue
	   endif
	else if(if.eq.-1)then
c
c  Inverse problem
	call invert(ig)
	call getchr('outfile',filnam,nfound)
	if(nfound.eq.0)filnam='coremodel'
	write(iout,'(a)')' Core field coords and br written to ',filnam
	call flush(iout)
	open(unit=21,file=filnam)
	   do 20 i=1,n
	     p(4,i)=bc(i)
	     write(21,'(5(1x,g15.7))') (p(j,i),j=1,4),bc(n+i)
20	   continue
	endif
	call flush(21)
c
c  Find SH representation of this core field
	call getval('shrep',val,2,nfound)
	ldeg=nint(val(1))
	r=val(2)
	if (nfound.le.0.or.ldeg.le.0)goto 23
	do 21 i=1,nepoch
	if(if.eq.1)then
c
c  Cycle core points into p(4,j)
c
	   do 24 k=1,n
	     p(4,k)=bc(k+(i-1)*n)
24	   continue
	endif
	call shexp(ldeg,r)
	write(iout,*)' Spherical Harmonic Representation follows: '
	do  21 l=0,ldeg
	do 22 m=0,l
	write(iout,*)l,m,glm(l+1,m+1),hlm(l+1,m+1)
22	continue
21	continue
c
c  For each epoch model find patches on the core where Br has one sign
c
23	call getval('patches',val,1,nfound)
	ip=nint(val(1))
	if(ip.eq.1)then
	  do 101 i=1,nepoch
c
c  Cycle core points into p(4,j)
c
	   do 224 k=1,n
	     p(4,k)=bc(k+(i-1)*n)
224	   continue
	   call patch(i)
	   write(iout,*)'no of patches = ',np(i)
c  Compute flux through the null flux curves and gradient of flux
c  wrt changes in field values at the vertices on the core
	   call flux(i)
	   fluxd0=0.0
	   fluxt=0.0
	   do 100 j=1,np(i)
	     write(iout,*)' flux through patch ',j,' is ',
     $        fluxv(j,i)*rc**2
	     fluxt=fluxt + abs(fluxv(j,i))
	     write(iout,*)'no of points',npts(J,i)
	  write(iout,*)'points in this patch are ',
     $           (lpts(ipts(j,i)+k-1,i),k=1,npts(j,i))
100	   continue
	   write(iout,*)'Total flux for this epoch is ',
     $      fluxt*rc**2
	   call flush(iout)
c  Get a pointer llpt to allow reordering of patches and gradient matrix according to 
c  magnitude of fluxv in each patch (llpt(np(i),j) points to largest patch at epoch j).
	   call fluxpt(llpt(1,i),np(i),fluxv(1,i))
c	   	   do 103 j=1,np(i)
c	     fluxv(llpt(np(i)+1-j,i),i)=0.0
c		 do 105 k=1,n
c		   fluxv(llpt(np(i)+1-j,i),i)=fluxv(llpt(np(i)+1-j,i),i) + 
c    $   gradm(llpt(np(i)+1-j,i),k,i)*p(4,k)
c105		 continue
c	     write(iout,*)' flux through patch ',j,' is ',
c     $ fluxv(llpt(np(i)+1-j,i),i)
c	     write(iout,*)'no of points',npts(llpt(np(i)+1-J,i),i)
c	  write(iout,*)'points in this patch are ',
c     $   (lpts(ipts(llpt(np(i)+1-j,i),i)+k-1,i),
c     $       k=1,npts(llpt(np(i)+1-j,i),i))
c103	   continue
101	   continue
	endif
c  Try and equivalence patches from the two epochs, npatch is the number
c  of distinct patches, those with no equivalent at other epoch are flagged
c  by 99 in llpt(k,i)
	call getchr('outfile',filnam,nfound)
	if(nfound.eq.0)filnam='coremodel'
	write(iout,'(a)')' Core field coords and br written to ',filnam
	open(unit=21,file=filnam)
c read number of fr lfux iterations to be performed
	call getval('fiter',val,1,nfound)
	if (nfound.eq.0)then
	     iter=0
	else
	     iter=nint(val(1))
	endif
	do 2000 jj=1,iter
	write(*,'(a,i3)')'Iteration number ',jj
	write(21,'(a,i3)')'Iteration number ',jj
	call patceq
c Compute flux difference norm squared
	   fluxd0=fluxd
	   fluxd=0.0
	   do 1999 i=1,npatch
	     if(llpt(i,1).eq.99)fluxv(llpt(i,1),1)=0.
	     if(llpt(i,2).eq.99)fluxv(llpt(i,2),2)=0.
	     fluxd=fluxd + (fluxv(llpt(i,1),1) - fluxv(llpt(i,2),2))**2
1999	   continue
	   write(iout,*)'Norm squared of flux diff is ',fluxd
	   fratio=fluxd0/fluxd
	   if(jj.eq.1.or.fratio.lt.1.0)fratio=1
	call frflux(ig,fratio)
	   do 200 i=1,n
	     p(4,i)=bc(i)
	     write(21,'(5(1x,g15.7))') (p(j,i),j=1,4),bc(n+i)
200	   continue
c
c  For each epoch model find patches on the core where Br has one sign
c
230	call getval('patches',val,1,nfound)
	ip=nint(val(1))
	if(ip.eq.1)then
	  do 111 i=1,nepoch
c
c  Cycle core points into p(4,j)
c
	   do 240 k=1,n
	     p(4,k)=bc(k+(i-1)*n)
240	   continue
	   call patch(i)
	   write(iout,*)'no of patches = ',np(i)
c  Compute flux through the null flux curves and gradient of flux
c  wrt changes in field values at the vertices on the core
	   call flux(i)
	   fluxt=0.0
	   do 110 j=1,np(i)
	     write(iout,*)' flux through patch ',j,' is ',fluxv(j,i)*rc**2
	     fluxt=fluxt + abs(fluxv(j,i))
	     write(iout,*)'no of points',npts(J,i)
	  write(iout,*)'points in this patch are ',
     $           (lpts(ipts(j,i)+k-1,i),k=1,npts(j,i))
110	   continue
	   write(iout,*)'Total flux for this epoch is ',
     $          fluxt*rc**2
c  Get a pointer llpt to allow reordering of patches and gradient matrix according to 
c  magnitude of fluxv in each patch (llpt(np(i),j) points to largest patch at epoch j).
	   call fluxpt(llpt(1,i),np(i),fluxv(1,i))
111	   continue
	endif
2000	continue
	stop
	end
c_______________________________________________________________________
	subroutine invert(ig)
c
	implicit double precision (a-h,o-z)
	parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,nvert=10)
	parameter (maxob=7000,maxp=30,lmax=101,maxe=2)
	common /io/inp,iout,iprint
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
	common /observ/ro(3,maxob,maxe),bro(maxob,maxe),sigma(maxob,
     $                maxe), iflag(maxob,maxe),mobs,nobs(maxe),nepoch,
     $                lobs(maxe)
	common /calcul/bs(maxob,maxe),a(maxob,nxyz),bc(nxyz2),
     $         reg(nxyz,nxyz)
	common /rule/xyw(3,100),nquad
      common /patchy/np(maxe),npts(maxp,maxe),isign(maxp,maxe),
     $   ipts(maxp,maxe),
     $   lpts(2*nxyz,maxe),ntris(nxyz),itris(nxyz+1),
     $   ltris(2*(2*nxyz-4)),lptr(2*nxyz-4,2),
     $   llpt(maxp,maxe),npatch
      common/patchg/  fluxv(maxp,maxe),gradm(maxp,nxyz,maxe)
	common /shr/glm(lmax,lmax),hlm(lmax,lmax),plm(lmax,lmax),
     $        dlm(lmax,lmax)
	common /circle/v(4,nxyz2),istack(nxyz2)
	dimension ba(maxob),parm(maxe),rlam(10),at(nxyz,nxyz),bat(nxyz)
	character*20 problem
	external b2int
	write(*,*)'npatch = ',npatch
c	calculate model dimension
	nn=nepoch*n
	write(iout,*)'nepoch= ',nepoch, ' nn= ',nn
c  clear out lower part of a
	do 1 i=1,nxyz
	do 1 j=1,nxyz
1	at(i,j)=0.0
c
c  
c  for new problem read data and make gram matrix, otherwise read
c  grammatrix
	if(ig.eq.0.or.ig.eq.1)then
c  
c  make vector of surface field, absent data are flagged by 99999.
c
  	  l=0
	   do 222 i=1,nepoch
	   lobs(i)=0.
	   do 2 k=1,nobs(i)
	      if(iflag(k,i).gt.0)then
	          l=l+1
		      lobs(i)=lobs(i)+1
		      if(iflag(k,i).gt.1)then
   	             ba(l)=bro(k,i)/sigma(k,i)
			  else
			     ba(l)=0.0
			  endif
	      endif
3	     continue
2	   continue
	     l=l+1
222	   continue
	   write(iout,*)'A total of ',l-nepoch,' measurements were used'
c
c  is unregularized problem under or over-determined?
c
	resq1=0.0
	kk=1
	k2k=1
	do 4 k=1,nepoch
	if (lobs(k).lt.n)then
c
c  problem is underdetermined
c
	write(iout,*)'Epoch ',k,' LS problem is underdetermined'
	   resq=0.0
c
c  copy a, and ba into lower part of arrays for later use
c
	  do 5 i=1,lobs(k)
	  do 6 j=1,nn
6 	   at(i+k2k-1,j)=a(i+k2k-1,j)
 	  bat(i+k2k-1)=ba(i+k2k-1)
5	  continue
	  do 65 j=1,n
65	  at(lobs(k)+1,j)=1.0
	  bat(lobs(k)+1)=0.
	  kk=kk+lobs(k)+1
	  k2k=k2k+lobs(k)+1
	else
c
c  do least squares accummmulation phase for this overdetermined part
c
c  Move relevant part of a to upper corner
c
	write(iout,*)'Epoch ',k,' LS problem is overdetermined'
	write(iout,*)'lobs= ',lobs(k), 'k2k= ',k2k,'kk= ',kk
	  do 7 i=1,lobs(k)
	    do 8 j=1,n
8	    a(i,j) = a(kk+i-1,j+(k-1)*n)
7	  ba(i)=ba(kk+i-1)
c  Include no monopole constraint
c	   do 85 j=1,n
c85	   a(lobs(k)+1,j)=1.0
	   ba(lobs(k)+1)=0.
	   ll=lobs(k)+1
	   call qr(maxob,ll,n,a,ba,bc,resq)
	   if (resq.le.0)
     $   write(iout,*)resq,' singular matrix in initial qr'
	   write(iout,*)'Epoch ',k
	   write(iout,*)'minimum possible rms misfit is ',
     $       sqrt(resq/float(lobs(k)))
	   call flush(iout)
	   resq1=resq1+resq
c  copy upper triangle of a, and ba into at and bat for later
c  use
	write(iout,*)' Copying upper triangle to at,bt'
	call flush(iout)
	  do 9 i=1,n
	  do 10 j=i,n
10 	  at(i+k2k-1,j+(k-1)*n)=a(i,j)
	  bat(i+k2k-1)=ba(i)
9	  continue
	  kk=kk+lobs(k)+1
	  k2k=k2k+n
	endif
4	continue
c
c  Write triangular form of a  and ba to disk
c
	open(unit=20,file='lsp',form='unformatted')
	write(iout,*)' writing upper triangle of at,bt'
	write(iout,*)'lobs(1),lobs(2)',lobs(1),lobs(2)
	write(iout,*)'l, k2k, nn',l,k2k,nn
	write(iout,*)'resq1',resq1
	call flush(iout)
	write(20)lobs(1),lobs(2),l,k2k,resq1
	do 90 j=1,nn
	write(20)(at(i,j),i=1,k2k-1)
90	continue
	write(20)(bat(j),j=1,k2k-1)
	close(unit=20)
	write(iout,*)'Upper triangular form of design matrix and data
     $ written to disk'
	call flush(iout)
c
c  otherwise read upper triangle of a and ba from lsp
c  and save in at and bat
c	
	else if(ig.eq.-1.or.ig.eq.-2) then
  	  open(unit=20,file='lsp',form='unformatted')
	  read(20)lobs(1),lobs(2),l,k2k,resq1
	  do 91 j=1,nn
	  read(20)(at(i,j),i=1,k2k-1)
91	  continue
	  read(20)(bat(j),j=1,k2k-1)
	  close(unit=20)
	  write(*,*)l-nepoch,' measurements',k2k-1,' rows in at ',
     $      resq1,' = total misfit'
	endif
c
c  Decide what problem to solve
c 
	call getchr('invtype',problem,nfound)
c
c*********************************************************************
c  Start with single epoch case
c
	if(nepoch.eq.1)then
c  choose regularization type and solve
	if(problem(1:7).eq.'minbrsq')then
	   nf=2
	else if(problem(1:7).eq.'mingrad')then
	   nf=4
	else
	  write(iout,*)'Unknown regularization type, min brsq used'
	  nf=2
	endif
	call regul(nf)
	call getval('lambda',rlam,1,nfound)
	nlam=nfound
	write(iout,*)'rlam1 = ',rlam(1)
	rlam1=rlam(1)*float(l)
c  Construct new a matrix for solution by qr
	  if(nf.eq.4)then
             nnn=n-1
	  else
	     nnn=nn
	  endif
	  do 11 i=1,nnn
	    ba(i)=0.
	     do 12 j=1,n
12	     a(i,j)=reg(i,j)*sqrt(rlam1)
11	  continue
	  do 111 i=1,n
	    ba(nnn+i)=bat(i)
	  do 112 j=1,n
112	     a(nnn+i,j)=at(i,j)
	if(lobs(1).gt.n)then
	  do 113 j=1,i-1
113	   a(nnn+i,j)=0.
	endif
111	  continue
c  solve this regularized problem
c
	ll=nnn+k2k-1
	write(iout,*)'ll= ',ll,' nn= ',nn
	call qr(maxob,ll,nn,a,ba,bc,resqr)
	if (resqr.le.0)write(iout,*)resqr, ' singular matrix in regul qr'
	write(iout,*)'qr rms misfit =', sqrt(resqr/float(l))
	call flush(iout)
c
c  compute integrated regularization parameter 
c
	p1=0.0
	do 14 i=1,nnn
	 ba(i)=0.
	 do 15 j=1,n
	  ba(i)=ba(i) + reg(i,j)*bc(j)
15	 continue
14	continue
	parm(1)=0.
	do 16 j=1,nnn
	parm(1)=parm(1)+ ba(j)*ba(j)
16	continue
	p1=p1+parm(1)
	if(nf.eq.2)then
	write(iout,*)'integrated br**2 over core is ',parm(1)
	   write(iout,*)'regularized model rms misfit =',
     $      sqrt((resqr-rlam1*p1+resq1)/float(l))
	else if(nf.eq.4)then
	write(iout,*)'integrated (grad br)**2 over core is ',parm(1)
	   write(iout,*)'regularized model rms misfit =',
     $      sqrt((resqr-rlam1*p1+resq1)/float(l))
	endif
17	continue
	else if(nepoch.eq.2)then
c
c********************************************************************
c  Now have two epochs
c
c
	if(problem(1:7).eq.'minbrsq'.or.problem(1:7).eq.'mingrad')then
	if(problem(1:7).eq.'minbrsq')then
	   nf=2
	   nnn=n
	else if(problem(1:7).eq.'mingrad')then
	   nf=4
	   nnn=n-1
	endif
	call regul(nf)
	write(iout,*)'Return from regul'
	call flush(iout)
	call getval('lambda',rlam,2,nfound)
	nlam=nfound
	write(iout,*)'rlam1 = ',rlam(1)
	rlam1=rlam(1)*float(l)
	write(iout,*)'rlam2 = ',rlam(2)
	rlam2=rlam(2)*float(lobs(2))/float(lobs(1))
c  expand regularization matrix to handle two epochs
	do 20 i=1,nnn
	do 22 j=1,n
	  reg(i,j+n)=0.0
	  reg(i+nnn,j)=0.0
	  reg(i+nnn,j+n) = reg(i,j)
22	continue
20	continue
c  Construct new a matrix for solution by qr
	  do 23 i=1,nn
	    ba(2*nnn+i)=bat(i)
	     do 24 j=1,nn
24	     a(2*nnn+i,j)=at(i,j)
23	  continue
	  do 25 i=1,2*nnn
	   ba(i)=0.
	   do 25 j=1,nn
	   a(i,j)=reg(i,j)*sqrt(rlam1)
25	  continue
c
c
c Allow for individual epoch satisfaction of misfit criteria
	kk=k2k-n-1
	write(iout,*)'First ', kk,' data weighted by sqrt(rlam2)'
	do 26 i=1,kk
	  ba(2*nnn+i)=sqrt(rlam2)*ba(2*nnn+i)
	do 26 j=1,nn
	  a(2*nnn+i,j)=sqrt(rlam2)*a(2*nnn+i,j)
26	continue
	ll=2*nnn+k2k-1
c  
c  solve this regularized problem
c
	write(iout,*)'ll= ',ll
	call qr(maxob,ll,nn,a,ba,bc,resqr)
	if (resqr.le.0)write(iout,*) resqr, ' singular matrix in regul qr'
	write(iout,*)'qr rms misfit =', sqrt(resqr/float(l))
c
c  compute integrated regularization parameter for each epoch
c
	p1=0.0
	do 27 k=1,nepoch
	do 28 i=1,nnn
	 ba(i+(k-1)*n)=0.
	 do 29 j=1,n
	  ba(i+(k-1)*n)=ba(i+(k-1)*n) + reg(i,j)*bc(j+(k-1)*n)
29	 continue
28	continue
	parm(k)=0.
	do 30 j=1,nnn
	parm(k)=parm(k)+ ba(j+(k-1)*n)*ba(j+(k-1)*n)
30	continue
	p1=p1+parm(k)
	write(iout,*)'Epoch ',k
	if(nf.eq.2)then
	  write(iout,*)'integrated br**2 over core is ',parm(k)
	else if(nf.eq.4)then
	  write(iout,*)'integrated (grad br)**2 over core is ',parm(k)
	endif
27	continue
	   write(iout,*)'regularized model rms misfit =',
     $      sqrt((resqr-rlam1*p1+resq1)/float(l))
32	continue
	diff=0.
	do 325 j=1,nnn
325	diff =diff + (ba(j)-ba(j+n))**2
	write(iout,*)'Integrated squared difference in models is',diff	
c
	endif
	endif
	return 
	end
c
c_______________________________________________________________________
	subroutine regul(nf)
	implicit double precision (a-h,o-z)
c  finds the matrix reg(n,n) to do integration of (Br)**2 or 
c  (grad Br)**2 over the surface of the core
c  so that br.reg.br yields the integral
c  then cholesky factor of reg is computed for regularized solution by
c  QR in invert.
c  Polyhedral core is used to enable closed form of integral
c  nf is type of regularization
c  nf=2 for br, nf=4 for grad br.
	parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,nvert=10)
	parameter (maxob=7000,maxp=30,lmax=101,maxe=2)
	common /io/inp,iout,iprint
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
	common /calcul/bs(maxob,maxe),a(maxob,nxyz),bc(nxyz2),
     $         reg(nxyz,nxyz)
	common /rule/xyw(3,100),nquad
c
c  Initialize matrix to zero
c
	do 1 i=1,n
	do 1 j=1,n
1	reg(i,j)=0.
	area=0.
	areas=0.
c
c  Add contribution to matrix from each triangle
c
	do 100 j=1,ntria
c
c  Vertices of current triangle are points iv1,iv2,iv3
c
	iv1=iv(1,j)
	iv2=iv(2,j)
	iv3=iv(3,j)
c
c  find area of current triangle
c
	delta=p3area(p(1,iv1),p(1,iv2),p(1,iv3))
	sdelta=sarea(p(1,iv1),p(1,iv2),p(1,iv3))
	area=area +delta
	areas=areas+sdelta
c
c  for br regularization, construct a symmetric +ve definite matrix
c
	if (nf.eq.2)then
	  const=delta/6.
	  reg(iv1,iv1)=const + reg(iv1,iv1)
	  reg(iv2,iv2)=const + reg(iv2,iv2)
	  reg(iv3,iv3)=const + reg(iv3,iv3)
	  reg(iv1,iv2)=const/2. + reg(iv1,iv2)
	  reg(iv1,iv3)=const/2. + reg(iv1,iv3)
	  reg(iv2,iv3)=const/2. + reg(iv2,iv3)
	  reg(iv2,iv1)=const/2. + reg(iv2,iv1)
	  reg(iv3,iv1)=const/2. + reg(iv3,iv1)
	  reg(iv3,iv2)=const/2. + reg(iv3,iv2)
c
c
c  otherwise for grad br regularization, get a symmetric, semi+ve def
c
	else if (nf.eq.4)then
	  s12=2*(1-dotty(p(1,iv2),p(1,iv3)))
	  s22=2*(1-dotty(p(1,iv1),p(1,iv3)))
	  s32=2*(1-dotty(p(1,iv2),p(1,iv1)))
	  const=4.*sdelta
	  reg(iv1,iv1)=s12/const + reg(iv1,iv1)
	  reg(iv2,iv2)=s22/const + reg(iv2,iv2)
	  reg(iv3,iv3)=s32/const + reg(iv3,iv3)
	  reg(iv1,iv2)=(s32-s12-s22)/(2.*const) + reg(iv1,iv2)
	  reg(iv1,iv3)=(s22-s12-s32)/(2.*const) + reg(iv1,iv3)
	  reg(iv2,iv3)=(s12-s22-s32)/(2.*const) + reg(iv2,iv3)
	  reg(iv2,iv1)=(s32-s12-s22)/(2.*const) + reg(iv2,iv1)
	  reg(iv3,iv1)=(s22-s12-s32)/(2.*const) + reg(iv3,iv1)
	  reg(iv3,iv2)=(s12-s22-s32)/(2.*const) + reg(iv3,iv2)
c
c use eigenvalue decomposition to find square root of reg
c
	endif
100	continue
	if(nf.eq.2) then
c  Find Cholesky decomposition of reg so can solve by QR
c  Upper triangular form returned
	  call cholsk(n,reg,nxyz,ifbad)
	else if(nf.eq.4) then
	  call ql(nxyz,n,reg,bc,a,ierr)
	  if(ierr.ne.0)then
	   write(iout,*) 'Trouble in eigenvalue decomposition of reg'
	  endif
c  compute sqrt(lambda)Ts(l), discard two zero eigenvalues and vectors
	  do 10 i=1,n
	  do 11 j=1,n
11	  a(i,j)=reg(j,i)
10	  continue
	  do 12 i=1,n-1
	  do 13 j=1,n
13	  reg(i,j)=a(i+1,j)
	  reg(i,i)=sqrt(bc(i+1))*reg(i,i)
12	  continue
	endif
	write(iout,*)' area of regularization polyhedron = ',area
	write(iout,*)' area of regularization sphere = ',areas
	return
	end
C
C
C
      SUBROUTINE CHOLSK(N,A,LDA,IFBAD)
	implicit double precision (a-h,o-z)
C     COURTESY BOB PARKER
C     CALLS NO OTHER ROUTINES
C     FACTORIZATION OF POSITIVE DEFINITIVE SYMMETRIC MATRICES
C     EXPLANATION        
C     CHOLSK CALCULATES THE FACTOR L IN THE CHOLESKY FACTORIZATION
C     A=L L-TRANS. THE ORIGINAL ARRAY IS OVERWRITTEN.
C     ARGUMENTS
C     A AN ARRAY DIMENSION LDA,LDA
C     ON ENTRY THE INPUT MATRIX SHOULD OCCUPY THE 1ST N COLS, OR IF
C     DESIRED, JUST THE LOWER LEFT TRIANGLE.  ON EXIT THE ORIGINAL
C     ARRAY IS REPLACED WITH L THE LEFT CHOLESKY FACTOR, IN THE UPPER HALF.
C     N THE ORDER OF THE MATRIX
C     IFBAD ERROR FLAG, SET TO 1 IF MATRIX NOT POSITIVE DEFINITE. SET TO
C     0 IF OK.
      DIMENSION A(LDA,LDA)
C
      IFBAD=0
      DO 1500 I=1,N
      I1=I+1
      DO 1400 J=I,N
      J1=J+1
      X=A(J,I)
      IF(I.EQ.1)GOTO 1200
      DO 1100 KBACK=2,I
      K=I-KBACK+1
1100  X=X-A(K,J1)*A(K,I1)
1200  IF(I.NE.J)GOTO 1300
      IF(X.LE.0.0)GOTO 4000
      Y=1.0/SQRT(X)
      A(I,I1)=Y
      GOTO 1400
1300  A(I,J1)=X*Y
1400  CONTINUE
1500  CONTINUE
C
      DO 1600 I=1,N
1600  A(I,I+1)=1.0/A(I,I+1)
C
C     COPY L INTO LOWER HALF OF MATRIX AND ZERO UPPER HALF.
      DO 2100 I=1,N
      DO 2100 J=1,I
      A(I,J)=A(J,I+1)
2100  A(J,I+1)=0.0  
C     take transpose to yield upper triangular form   
      DO 220 I=1,N
      DO 220 J=1,I
      A(J,I)=A(I,J)
220   IF(I.NE.J)A(I,J)=0.D0
      RETURN
C
C     nonpositive matrix encountered
4000  ifbad=1
	write(iout,*)'non-positive matrix encountered in Cholsk'
      return
      end 
c_______________________________________________________________________
	function p3area(a,b,c)
c  computes area of planar triangle specified by 3 vectors
	implicit double precision (a-h,o-z)
	dimension a(3),b(3),c(3),ab(3),ac(3),ar(3)
	do 5 j=1,3
	  ab(j)=b(j)-a(j)
	  ac(j)=c(j)-a(j)
5	continue
	call axb(ab,ac,ar)
	p3area=sqrt(dotty(ar,ar))/2.
	return
	end
c_______________________________________________________________________
      function tegf(f,p1,p2,p3)
	implicit double precision (a-h,o-z)
c$$$$$$$calls centr, gpole, gnomon, f, parea
c  integrates the function f over the spherical triangle defined by
c  vertices p1, p2, p3
c  core field br at vertices p, is in 4th component of arrays
c
c  The quadrature rule generated by 'triq' for integrals over
c  plane triangle with vertices  (0,0,0) (1,0,0) (0,1,0) is assumed
c  to be available in common /rule/.
c
c
	common /io/inp,iout,iprint
	common /rule/xyw(3,100),nquad
      dimension rp(3),vc(4),aa(3),b(3),c(3),p1(4),p2(4),p3(4),
     $          ua(2),ub(2),uc(2),rg(2),rpa(3)
c
	tegf=0.
c
c  compute centroid vc of triangle
c
	call centr(p1,p2,p3,vc)
c
c  rotate so vertices are in coord system with vc as pole, p(1,j) at
c  zero longitude
	call gpole(1,vc,p1,p1,aa)
	call gpole(1,vc,p1,p2,b)
	call gpole(1,vc,p1,p3,c)
c
c  find gnomonic projection of vertices and jacobian of transformation
	call gnomon(1,ua,aa)
	call gnomon(1,ub,b)
	call gnomon(1,uc,c)
	ajac = parea (ua,ub,uc)*2.
c
c  Run through the cubature rule for this triangle 
c
        do 2400 k=1, nquad
c
c  find points on gnomonic surface corresponding to integration points
c
	rg(1)= ua(1) + (ub(1)-ua(1))*xyw(1,k)  +(uc(1)-ua(1))*xyw(2,k)
	rg(2)= ua(2) + (ub(2)-ua(2))*xyw(1,k)  +(uc(2)-ua(2))*xyw(2,k)
c
c  find points in physical space corres. to rg
c
	call gnomon(-1,rg,rpa)
	call gpole(-1,vc,p1,rp,rpa)
c
c  evaluate function f at rp, rg, etc.
c
	gam1=f(rg,rp,k,p1,p2,p3)
c
c  Weight integral according to xyw rule
c  Find Jacobian for gnomonic to spherical projection
c
	bjac= sqrt((1+rg(1)**2 + rg(2)**2)**-3)
	tegf=tegf + gam1*xyw(3,k)*ajac*bjac
 2400   continue
      return
      end
c_______________________________________________________________________
	function b2int(rg,rp,k,p1,p2,p3)
	implicit double precision (a-h,o-z)
c  interpolates b to a point rg on the gnomonic projection
	common /rule/xyw(3,100),nquad
	dimension rp(3),rg(2),gamma(3),p1(4),p2(4),p3(4)
	gamma(1)=1.-xyw(1,k)-xyw(2,k)
	gamma(2)=xyw(1,k)
	gamma(3)=xyw(2,k)
	b2int=(gamma(1)*p1(4) +gamma(2)*p2(4) +gamma(3)*p3(4))**2
	return
	end
c____________________________________________________________________________________________________
	subroutine fluxpt(llpt,np,fluxv)
	implicit double precision (a-h,o-z)
c  get a ptr to order fluxv(np) in magnitude
c  ignore patches with no flux i.e. treat as null flux points
	common /io/inp,iout,iprint
	parameter (maxob=7000,maxp=30,lmax=101,maxe=2)
	dimension x(maxp,maxe),fluxv(np),llpt(np)
	data rc/3485/
	k=0
	do 10 i=1,np
	if(abs(fluxv(i)).gt.1.) then
	  k=k+1
	  x(k,1)=(fluxv(i))
	  x(k,2)=float(i)
	endif
10	continue
        np=k
	call sort22(maxp,np,x)
	do 15 i=1,np
	  llpt(i)=nint(x(i,2))
	  write(*,*)i,fluxv(i)*rc**2,llpt(i),fluxv(llpt(i))*rc**2
15	continue
	return
	end
c____________________________________________________________________________________________________c_______________________________________________________________________
      subroutine sort22(maxp,no, x)
	implicit double precision (a-h,o-z)
c$$$$$ calls no other routines
c  In-place rapid sorter.  Rearranges  x  into ascending order of x(i,1)
      dimension x(maxp,2)
c
      mo=no
 2    if (mo-15) 21,21,23
 21   if (mo.le.1) return
      mo=2*(mo/4) + 1
      go to 24
 23   mo=2*(mo/8) + 1
 24   ko=no - mo
      jo=1
 25   i=jo
 26   if (x(i,1) - x(i+mo,1)) 28,28,27
 27   temp=x(i,1)
 	temp2=x(i,2)
      x(i,1)=x(i+mo,1)
	x(i,2)=x(i+mo,2)
      x(i+mo,1)=temp
	x(i+mo,2)=temp2
      i=i - mo
      if (i-1) 28,26,26
 28   jo=1 + jo
      if (jo-ko) 25,25,2
      end 
c______________________________________________________________________________________
 	subroutine frflux(ig,fratio)
c
	implicit double precision (a-h,o-z)
	parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,nvert=10)
	parameter (maxob=7000,maxp=30,lmax=101,maxe=2)
	common /io/inp,iout,iprint
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
	common /observ/ro(3,maxob,maxe),bro(maxob,maxe),sigma(maxob,
     $                maxe), iflag(maxob,maxe),mobs,nobs(maxe),nepoch,
     $                lobs(maxe)
	common /calcul/bs(maxob,maxe),a(maxob,nxyz),bc(nxyz2),
     $         reg(nxyz,nxyz)
	common /rule/xyw(3,100),nquad
      common /patchy/np(maxe),npts(maxp,maxe),isign(maxp,maxe),
     $   ipts(maxp,maxe),
     $   lpts(2*nxyz,maxe),ntris(nxyz),itris(nxyz+1),
     $   ltris(2*(2*nxyz-4)),lptr(2*nxyz-4,2),
     $   llpt(maxp,maxe),npatch
      common/patchg/  fluxv(maxp,maxe),gradm(maxp,nxyz,maxe)
	common /shr/glm(lmax,lmax),hlm(lmax,lmax),plm(lmax,lmax),
     $        dlm(lmax,lmax)
	common /circle/v(4,nxyz2),istack(nxyz2)
	dimension ba(maxob),parm(maxe),rlam(10),at(nxyz,nxyz),bat(nxyz)
	character*20 problem
	external b2int
	write(*,*)'npatch = ',npatch
c	calculate model dimension
	nn=nepoch*n
      write(*,*)'nn = ',nn
c  clear out lower part of a
	do 1 i=1,nn+3
	do 1 j=1,nn+3
1	at(i,j)=0.0
c
c  
  	  open(unit=20,file='lsp',form='unformatted')
	  read(20)lobs(1),lobs(2),l,k2k,resq1
	  do 2 j=1,nn
	  read(20)(at(i,j),i=1,k2k-1)
2	  continue
	  read(20)(bat(j),j=1,k2k-1)
	  close(unit=20)
	  write(*,*)'Epoch 1 ',lobs(1), ' measurements'
	  write(*,*)'Epoch 2 ',lobs(2), ' measurements'
	  write(*,*)l-2,' measurements',k2k-1,' rows in at ',
     $      resq1,' = total misfit'
c
c
c  Attempt to minimize difference in integrated flux over equivalent
c  patches.  patch, flux and patceq must have been called first
c
c  Put flux difference matrix in to regularize. It will have npatch rows.
c
	call getval('lambda',rlam,3,nfound)
	nlam=nfound
	write(iout,*)'rlam1 = ',rlam(1)
	rlam1=rlam(1)*float(lobs(1))
	write(iout,*)'rlam2 = ',rlam(2)
	rlam2=rlam(2)
	write(iout,*)'rlam3 = ',rlam(3)
	rlam3=rlam(3)*float(lobs(1))*fratio
c
c  Construct new a matrix for solution by qr. First part is flux difference mninimization
c
	write(*,*)npatch ,' patches acting as constraints'
	nnpatch=npatch
	lrow=0
	jflag1=0
	jflag2=0
	do 40 i=1,npatch
	   ba(i)=0.
	   do 40 j=1,nn
40	a(i,j)=0.
	do 36 i=1,npatch
c  For each row in flux constraint matrix check that each patch has
c  not appeared before and use this to decide which row  (nrow)
c  constraint gets added to
	    write(*,*)'LLPT(i,1),LLPT(i,2)',llpt(i,1),llpt(i,2)
	    do 1111 k=1,i-1
	       if(llpt(k,1).eq.llpt(i,1).and.llpt(k,1).ne.99) then
	         write(*,*)llpt(k,1),llpt(i,1)
	         nrow=k
	         nnpatch=nnpatch-1
	         jflag1=1
	         goto 3
	       else if(llpt(k,2).eq.llpt(i,2).and.llpt(k,2).ne.99) then
	         write(*,*)llpt(k,2),llpt(i,2)
	         nrow=k
	         nnpatch=nnpatch-1
	         jflag2=1
	         goto 3
	       endif
1111	    continue
	    lrow=lrow+1
	    nrow=lrow
3	    write(iout,*)nrow,llpt(i,1),llpt(i,2)
	    if(jflag1.eq.0.and.jflag2.eq.0)then
	     if(llpt(i,2).ne.99.and.llpt(i,1).ne.99) then
c  Try and equivalence flux at epoch 1 and 2
	       ba(nrow)=fluxv(llpt(i,2),2)*sqrt(rlam3)
	       do 37 j=1,n
37	       a(nrow,j)=a(nrow,j) + gradm(llpt(i,1),j,1)*sqrt(rlam3)
	     else if(llpt(i,2).eq.99) then
c  Try and reduce flux through epoch 1 patch to zero
	       ba(nrow)=0.
	       do 38 j=1,n
38	       a(nrow,j)=a(nrow,j) + gradm(llpt(i,1),j,1)*sqrt(rlam3)
	     else if(llpt(i,1).eq.99) then
c  Make a new patch at epoch 1 for compatibility with 2
	       ba(nrow)=fluxv(llpt(i,2),2)*sqrt(rlam3)
	       call pmaker(llpt(i,2),1)
	       do 34 j=1,n
34	       a(nrow,j)=a(nrow,j) + gradm(maxp,j,1)*sqrt(rlam3)
	       write(iout,*)'Flux through pmaker patch is ',
     $         fluxv(maxp,1)
	     endif
	   else if (jflag1.eq.1) then
c  Equate 1 patch at epoch 1 to 2 at epoch 2
	     ba(nrow)=ba(nrow) + fluxv(llpt(i,2),2)*sqrt(rlam3)
	   else if(jflag2.eq.1) then
c  Equate 2 patches at epoch 1 to 1 at epoch 2
	     do 39 j=1,n
39	     a(nrow,j)=a(nrow,j) + gradm(llpt(i,1),j,1) *sqrt(rlam3)
	   endif
	  jflag1=0
	  jflag2=0
36	continue
        write (*,*) 'patch part of matrix complete'
	write(*,*)nnpatch,' patches acting as constraints'
	npatch=nnpatch
	if(k2k-1.eq.nn)then
c  Original ls problem was overdetermined
c  Add in triangular part of data in alternating rows
c  Allow for individual epoch satisfaction of misfit critieria
	   kk=k2k-n-1
	   write(iout,*)' First  ',kk, ' data weighted by sqrt(rlam2)'
	   do 409 i=1,kk
	     ba(npatch+2*i-1)=sqrt(rlam2)*bat(i)
	     do 409 j=1,n
	     a(npatch+2*i-1,j)=sqrt(rlam2)*at(i,j)
409	   continue
c
c  Add regularization constraint
c
cc
	   if(problem(1:7).eq.'minbrsq')then
	     nf=2
	     nnn=n
	   else if(problem(1:7).eq.'mingrad')then
	     nf=4
	     nnn=n-1
	   else 
	     nf=2
	     nnn=n
	   endif
	   write(*,*)' Computing regularization matrix for epoch 1'
	   call regul(nf)
	   write(*,*)'returned from regul'
	  do 25 i=1, nnn
	   ba(npatch+2*i)=0.
	   do 25 j=1,n
	   a(npatch + 2*i,j)=reg(i,j)*sqrt(rlam1)
25	  continue
	  nkk=npatch
	else
c  Original ls problem was underdetermined
c  Put data in uppermost block with gradients matrix
c  Allow for individual epoch satisfaction of misfit critieria
	kk=k2k-n-1
	write(iout,*)' First  ',kk, ' data weighted by sqrt(rlam2)'
	do 419 i=1,kk
	  ba(npatch+i)=sqrt(rlam2)*bat(i)
	  do 419 j=1,n
	  a(npatch+i,j)=sqrt(rlam2)*at(i,j)
419	continue
c
c  Add regularization constraint on both epochs
c
cc
	if(problem(1:7).eq.'minbrsq')then
	   nf=2
	   nnn=n
	else if(problem(1:7).eq.'mingrad')then
	   nf=4
	   nnn=n-1
	else 
	   nf=2
	   nnn=n
	endif
	write(*,*)' Computing regularization matrix for epoch 1'
	call regul(nf)
	write(*,*)'returned from regul'
	  do 215 i=1, nnn
	   ba(npatch+k2k-1+i)=0.
	   do 215 j=1,n
	   a(npatch + k2k-1+i,j)=reg(i,j)*sqrt(rlam1)
215	  continue
	  nkk=npatch + kk
	endif
c  Save pathc part of matirx in at,bat
	do 35 i=1,npatch
	   bat(i)=ba(i)
	do 35 j=1,n
35	at(i,j)=a(i,j)
c
	ll=npatch + kk + nnn
c  
c  solve this regularized problem
c
	write(iout,*)'ll= ',ll, ' nkk= ',nkk
	call sspqr(maxob,ll,n,a,ba,bc,resqr,nkk)
	if (resqr.le.0)write(iout,*) resqr, ' singular matrix in regul qr'
	write(iout,*)'qr rms misfit =', sqrt(resqr/float(lobs(1)))
c
c  compute approx integrated difference parameter between epochs
c
	fdiff=0.0
	do 43 i=1,npatch
	 ba(i)=0.
	 do 44 j=1,n
	  ba(i)=ba(i) + at(i,j)*bc(j)
44	 continue
	 write(iout,*)'Nrow ',i, ' flux ', ba(i),
     $             ' target flux ',bat(i)
         fdiff=fdiff + (ba(i)-bat(i))**2
43	continue
	write(iout,*)'linearized approx to squared integrated diff. 
     $  in flux over core is ', fdiff
     c
c  compute integrated regularization parameter for each epoch
c
	p1=0.0
	do 27 k=1,nepoch
	do 28 i=1,nnn
	 ba(i+(k-1)*n)=0.
	 do 29 j=1,n
	  ba(i+(k-1)*n)=ba(i+(k-1)*n) + reg(i,j)*bc(j+(k-1)*n)
29	 continue
28	continue
	parm(k)=0.
	do 30 j=1,nnn
	parm(k)=parm(k)+ ba(j+(k-1)*n)*ba(j+(k-1)*n)
30	continue
	p1=p1+parm(k)
	write(iout,*)'Epoch ',k
	if(nf.eq.2)then
	write(iout,*)'integrated br**2 over core is ',parm(k)
	else if(nf.eq.4)then
	write(iout,*)'integrated (grad br)**2 over core is ',parm(k)
	endif
27	continue
	   write(iout,*)'regularized model rms misfit =',
     $      sqrt((resqr-rlam3*fdiff+resq1-rlam1*p1)/float(l))
      return
	end
c_______________________________________________________________________
c______________________________________________________________________
      subroutine dorito(n, p, ntria, iv, v, istack)
	implicit double precision (a-h,o-z)
c$$$$$ calls scirc, dotty, area
c  Routine to assign triangles to a set of  n  3-vectors  p,on the
c  surface of the unit sphere based upon
c  the method of Watson, Computers & Geosciences, v8,1, pp 97-101, 1982.
c  the final network is a set of Delaunay triangles with the property
c  that no vertex lies inside the circumcircle of any triangle.  This
c  system is as close to equiangular as possible.
c  n   is the number of input points.
c  p   is an array of data, dimensioned  p(4, n), where 
c      the 4-th component in each vector is ignored.  It is
c      present to allow field data in the fourth element.
c  ntria on exit contains the total number of triangles generated.
c  iv  on exit contains the vertices of the tessellation -
c      iv(1,k),iv(2,k),iv(3,k)  defines the corners of the  k-th
c      triangle as the points  (p(1,iv(1,k), p(2,iv(1,k)),p(3,iv(1,k)),
c      (p(1,iv(2,k), p(2,iv(2,k)), p(3,iv(3,k)),  
c      (p(1,iv(3,k)), p(2,iv(3,k),p(3,iv(3,k)).
c      iv  must be dimensioned at least as iv(3, 2*n).
c  v(i,k), i=1,,4 has coordinates of circumcenter & "radius" of circle
c      circumscribing triangle k in list iv.
c  istack  another work array of length 2*n
c
c
      dimension p(4,*), v(4,*), iv(3,*),istack(*),prism(3,6),iprism(6)
      dimension itemp(3,2),kv(2,500),dotmax(6)
      common /bound/dmin,dmax,amin,amax
	  common /io/ inp,iout,iprint
      data ((prism(i,j),i=1,3),j=1,6)/0,0,1,1,0,0,0,1,0,-1,0,0,0,
     $             -1,0,0,0,-1/
c
c
      data itemp/1,1,2,2,3,3/,fpi/12.56637062/
c
c  Find closest points in array p to those in prism and use these as a
c  starting configuration
c
c  Initialize dotmax
	do 1 j=1,6
1	dotmax(j)=0.0
c  Compare every point to prism points, keep the index of those closest
c
	do 10 j=1,n
	  do 11 i=1,6
	    x=dotty(p(1,j),prism(1,i))
	    if  (x.gt.dotmax(i)) then
	      dotmax(i)=x
	      iprism(i)=j
	    endif
11	  continue
10	continue
c  Now make a tesselation from the points in iprism and save in iv(i,j)
c
	iv(1,1)=iprism(1)
	iv(2,1)=iprism(2)
	iv(3,1)=iprism(3)
	iv(1,2)=iprism(1)
	iv(2,2)=iprism(3)
	iv(3,2)=iprism(4)
	iv(1,3)=iprism(1)
	iv(2,3)=iprism(4)
	iv(3,3)=iprism(5)
	iv(1,4)=iprism(1)
	iv(2,4)=iprism(2)
	iv(3,4)=iprism(5)
	iv(1,5)=iprism(2)
	iv(2,5)=iprism(3)
	iv(3,5)=iprism(6)
	iv(1,6)=iprism(3)
	iv(2,6)=iprism(4)
	iv(3,6)=iprism(6)
	iv(1,7)=iprism(4)
	iv(2,7)=iprism(5)
	iv(3,7)=iprism(6)
	iv(1,8)=iprism(2)
	iv(2,8)=iprism(5)
	iv(3,8)=iprism(6)
c  put in starting pointers
c
      n1=2*n 
      do 6 i=1, n1
        istack(i)=i
 6    continue
c
c  Starting configuration for the sphere is specified by
c  8 right triangles, in iv(i,j),j=1,...8
c  Compute the circumcentres and "radii" for these triangles.
	do 12 i=1,8
	 call scirc(p(1,iv(1,i)),p(1,iv(2,i)),p(1,iv(3,i)),v(1,i))
12	continue
      isp=8
      id=9
      km=1
c  Scan through the data forwards.
      do 50 nuc=1, n
c  Skip this point if it is part of the starting octahedron.
	do 29 i=1,6
	  if(nuc.eq.iprism(i))goto 50
29	continue
        km=0
c  Loop through the established 3-tuples.
c
        do 30 jt=1, isp
c  Test if the new data point is within the  jt  circumcircle.
          rsq=dotty(p(1,nuc),v(1,jt))
          if (rsq.lt.v(4,jt)) go to 30
c  The point is within.  Delete this 3-tuple but save its edges.
          id=id - 1
          istack(id)=jt
c  Add edges to  kv  but delete if already present.
          do 28 i=1, 3
            l1=itemp(i,1)
            l2=itemp(i,2)
            if (km.le.0) go to 26
            kmt=km
            do 24 j=1, kmt
              if (iv(l1,jt) .ne. kv(1,j)) goto 24
              if (iv(l2,jt) .ne. kv(2,j)) goto 24
              km=km - 1
              if (j .gt. km) goto 28
              do 20 k=j, km
                 kv(1,k)=kv(1,k+1)
                 kv(2,k)=kv(2,k+1)
 20           continue
              go to 28
 24         continue
 26         km=1 + km
            kv(1,km)=iv(l1,jt)
            kv(2,km)=iv(l2,jt)
 28       continue
 30     continue
c
c  Form new 3-tuples.
        do 48 i=1, km
          kt=istack(id)
          id=1 + id
c  Calculate the circumcircle center and radius squared.
          i1=kv(1,i)
          i2=kv(2,i)
	  call scirc(p(1,i1),p(1,i2),p(1,nuc),v(1,kt))
          iv(1,kt)=kv(1,i)
          iv(2,kt)=kv(2,i)
          iv(3,kt)=nuc
 48     continue
      isp=2 + isp
 50   continue
c
c  check no. of triangles is right and surface area integrates to 4pi.
c
      ntria=isp
	isp=2*n-4
	if(ntria.ne.isp) then
	  write(iout,*)'only ',ntria,'triangles, should be ',isp
	endif
      sarea=0.0
	amin=fpi
	amax=0.0
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
 60   continue
	darea=100.*sarea/fpi
        write(iout,*)'Surface area of tesselation is ',sarea
	write(iout,*)'This is ',darea,'% that of sphere'
	write(iout,*)'Minimum triangle area is ',amin
	write(iout,*)'Maximum triangle area is ',amax
      return
      end                                                               
c
c
c_______________________________________________________________________
	subroutine scirc(x1,x2,x3,c)
	implicit double precision (a-h,o-z)
c  for points x1,x2,x3,which are the vertices of a spherical triangle
c  calculates the circum circle origin c(1),c(2),c(3), and the dot
c  product of the centre with any vertex, c(4),
	dimension x1(3),x2(3),x3(3),c(4),a(3,3),b(3)
	      common /io/ inp,iout,iprint
	do 10 i=1,3
	  a(1,i)=x1(i)
	  a(2,i)=x2(i)
	  a(3,i)=x3(i)
	  b(i)=1.0
10	continue
c	write(6,*) resq,' resqin'
	call qr(3,3,3,a,b,c,resq)
c	write(6,*) resq,' resqout'
	if(resq.lt.0.)then
	  write(iout,*)'problem in qr'
	  write(iout,*) (x1(i),i=1,3)
	  write(iout,*) (x2(i),i=1,3)
	  write(iout,*) (x3(i),i=1,3)
	stop
	endif
c  normalize c
	cnorm=sqrt(dotty(c,c))
	do 11 i=1,3
	c(i)=c(i)/cnorm
11	continue
c  compute "radius"
	c(4)=dotty(c,x1)
	return
	end
c_______________________________________________________________________
	function dotty(x,y)
	implicit double precision (a-h,o-z)
c  computes the dot product of x and y
c
	dimension x(3),y(3)
	dotty=x(1)*y(1) + x(2)*y(2) +x(3)*y(3)
	return
	end
c
c_______________________________________________________________________
      function area(a, b, c)
	implicit double precision (a-h,o-z)
c$$$$ calls no other routines
c  Given three unit vectors a, b, c  finds the area of the
c  SPHERICAL triangle formed at the tips.
      common /io/ inp,iout,iprint
      dimension a(3),b(3),c(3)
c
      cosa=b(1)*c(1)+b(2)*c(2)+b(3)*c(3)
      cosb=a(1)*c(1)+a(2)*c(2)+a(3)*c(3)
      cosc=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      sina=sqrt(abs(1.D0 - cosa**2))
      sinb=sqrt(abs(1.D0 - cosb**2))
      sinc=sqrt(abs(1.D0 - cosc**2))
      AA=acos((cosa-cosb*cosc)/(sinb*sinc))
      BB=acos((cosb-cosa*cosc)/(sina*sinc))
      CC=acos((cosc-cosb*cosa)/(sinb*sina))
      area=AA+BB+CC - 3.14159265358979324
c      write(iout, '(9f8.4/a,f8.4)')a,b,c,' Area =',area
      return
      end
      subroutine qr(ndime,m, n, a, b, x, resq)
	implicit double precision (a-h,o-z)
c$$$$  calls no other routines
c  solves over-determined least-squares problem  ax = b
c  where  a  is an  m by n  matrix,  b  is an m-vector .
c  resq  is the sum of squared residuals of optimal solution.  also used
c  to signal error conditions - if -2 , system is underdetermined,  if
c  -1,  system is singular.
c  method - successive householder rotations.  see lawson+hanson - solv
c  -ing least squares problems.
c  routine will also work when m=n.
c*****   caution -  a and b  are overwritten by this routine.
      dimension a(ndime,n+1),b(ndime),x(n+1)
      double precision sum,dot
c
      resq=-2.0d0
      if (m.lt.n) return
c   loop ending on 1800 rotates  a  into upper triangular form
      do 1800 j=1,n
c  find constants for rotation and diagonal entry
      sq=0.00
      do 1100 i=j,m
 1100 sq=a(i,j)**2 + sq
      qv1=-sign(sqrt(sq),a(j,j))
      u1=a(j,j) - qv1
      a(j,j)=qv1
      j1=j + 1
      if (j1.gt.n) go to 1500
c  rotate remaining columns of sub-matrix
      do 1400 jj=j1,n
      dot=u1*a(j,jj)
      do 1200 i=j1,m
 1200 dot=a(i,jj)*a(i,j) + dot
      const=dot/abs(qv1*u1)
      do 1300 i=j1,m
 1300 a(i,jj)=a(i,jj) - const*a(i,j)
      a(j,jj)=a(j,jj) - const*u1
 1400 continue
c  rotate  b  vector
 1500 dot=u1*b(j)
      if (j1.gt.m) go to 1610
      do 1600 i=j1,m
 1600 dot=b(i)*a(i,j) + dot
 1610 const=dot/abs(qv1*u1)
      b(j)=b(j) - const*u1
      if (j1.gt.m) go to 1800
      do 1700 i=j1,m
 1700 b(i)=b(i) - const*a(i,j)
 1800 continue
c  solve triangular system by back-substitution.
      resq=-1.00
      do 2200 ii=1,n
      i=n-ii+1
      sum=b(i)
      if (ii.eq.1) go to 2110
      i1=i+1
      do 2100 j=i1,n
 2100 sum=sum - a(i,j)*x(j)
 2110 if (a(i,i).eq. 0.0d0) return
 2200 x(i)=sum/a(i,i)
c  find residual in overdetermined case.
      resq=0.0d0
      if (m.eq.n) return
      i1=n+1
      do 2300 i=i1,m
 2300 resq=b(i)**2 + resq
      return
      end
c
c
      subroutine sspqr(mdim,m,n,a,b,x,resq,kk)
	implicit double precision (a-h,o-z)
c     calls no other routines
c     solves overdetermined least squares problem ax=b
c     for special matrices a, where a has dimension 2*n+kk by n, and
c     a(i,j)=0., if i.gt.2*j+kk.  this special form appears in the depleted
c     basis analysis of harmonic splines with kk=0, and elsewhere.
c     resq is the sum of squared residuals of optimal solution.  also
c     used to signal error conditions-  if -2 the system is underdetermined, if
c     -1, system is singular.
c     method - successive householder rotations. see lawson & hanson,
c     solving least squares problems.
c**********caution - a and b are overwritten by this routine********************
      dimension a(mdim,n),b(mdim),x(n)
      double precision sum,dot
c
      resq=-2.0
c      if(m.ne.2*n+kk)return
c    loop ending on 1800 rotates a into upper triangular form.
      do 1800 j=1,n
c  find constants for rotation and diagonal entry
      sq=0.0d0                                                
c     dot products and other action computed only down to l=2*j+kk in each column
      l=2*j+kk
      do 1100 i=j,l
 1100 sq=a(i,j)**2 + sq
      qv1=-sign(sqrt(sq),a(j,j))
      u1=a(j,j) - qv1
      a(j,j)=qv1
      j1=j + 1
      if (j1.gt.n) go to 1500
c  rotate remaining columns of sub-matrix
      do 1400 jj=j1,n
      dot=u1*a(j,jj)
      do 1200 i=j1,l
 1200 dot=a(i,jj)*a(i,j) + dot
      const=dot/abs(qv1*u1)
      do 1300 i=j1,l
 1300 a(i,jj)=a(i,jj) - const*a(i,j)
      a(j,jj)=a(j,jj) - const*u1
 1400 continue
c  rotate  b  vector
 1500 dot=u1*b(j)
      if (j1.gt.l) go to 1610
      do 1600 i=j1,l
 1600 dot=b(i)*a(i,j) + dot
 1610 const=dot/abs(qv1*u1)
      b(j)=b(j) - const*u1
      if (j1.gt.l) go to 1800
      do 1700 i=j1,l
 1700 b(i)=b(i) - const*a(i,j)
 1800 continue
c  solve triangular system by back-substitution.
      resq=-1.0d0
      do 2200 ii=1,n
      i=n-ii+1
      sum=b(i)
      if (ii.eq.1) go to 2110
      i1=i+1
      do 2100 j=i1,n
 2100 sum=sum - a(i,j)*x(j)
 2110 if (a(i,i).eq. 0.0d0) return
 2200 x(i)=sum/a(i,i)
c  find residual in overdetermined case.
      resq=0.0d0
      if (m.eq.n) return
      i1=n+1
      do 2300 i=i1,m
 2300 resq=b(i)**2 + resq
      return
      end
      subroutine sbody(kdim, maxmum, n, p, iv, ntria)
	implicit double precision (a-h,o-z)
c$$$$ calls dorito,synops
c  If kdim=3, gets a series of x,y,z coords that define points on the 
c  surface of a sphere; if kdim=4 x,y,z coordinates plus associated
c  field value are read in.
c  Note that the position vector is dimensioned p(4,*) IN EITHER CASE.
c  Computes the tessellation or reads the appropriate
c  pointers from a file.
c
	parameter(nxyz=6500,nxyz2=2*nxyz)
      character*80 name
      common /io/ inp,iout,iprint
      common /bound/dmin,dmax,amin,amax
      common /circle/v(4,nxyz2),istack(nxyz2)
      dimension p(4,*), iv(*)
c
c
	call getval('tnew',val,1,nfound)
	new=nint(val)
c
c  Get x, y, z, triples from disk.
      if (kdim .eq. 3) call getchr('tessel',name,nfound)
      if (kdim .eq. 4) call getchr('corefile',name,nfound)
      write(iout,'(1x,a,a50)') 'Tesselation read from file: ',name
      open (unit=12, file=name)
      do 1100 i=1, nxyz
        read (12,*, err=1151,end=1150) (p(j,i), j=1, kdim)
 1100 continue
1151      write(iout,'(1x,a,i6,a/a)')
     $'Topography file was truncated at ',nxyz,' data',
     $' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
 1150 n=i-1
c  check unit vectors
	do 1200 i=1,n
	  cnorm=sqrt(dotty(p(1,i),p(1,i)))
	  p(1,i)=p(1,i)/cnorm
	  p(2,i)=p(2,i)/cnorm
	  p(3,i)=p(3,i)/cnorm
1200    continue
      close (unit=12)
      write(iout,'(1x,a,i6)') 'Number of coordinates read =',n
c
c
c  Pointers on disk.  Get pointer vectors from disk.
      if ( new .eq. 0 ) then
	call getchr('pointer',name,nfound)
	write(6,'(a)') name
        open (unit=13, file=name)
        maxv=2*nxyz-4
	ntria=0
	niv=1
	do 1999 k=1,maxv
         read (13,*,err=2001,end=2000)iv(niv),iv(niv+1),iv(niv+2)
	 niv= niv+3
	 ntria= ntria+1
1999	continue
2001        write (iout,'(//(1x,a/))')
     $  'Pointer file incompatible with topography data',
     $  '*******************************************'
        stop
 2000   write(iout,*) 'ntria = ', ntria
	if(ntria.ne.2*n-4)goto 2001
      write(iout,'(1x,a,a50)') 'Pointers read from file: ',name
        close (unit=13)
c
c  No pointers on disk.  Therefore generate pointer series with 
c  Delaunay algorithm. First perform the calculations.
      else
c 
        call dorito(n, p, ntria, iv, v, istack)
        write(iout,'(1x,a)') 'Delaunay phase complete'
c
c  Save pointers on the named file.
	  call getchr('pointer',name,nfound)
          open (unit=14, file=name)
          write(14,'(3I7)') (iv(j), j=1, 3*ntria)
          close (unit=14)
          write(iout,'(1x,a)') 'Pointers written to disk'
      endif
      call synops(new, n, p, ntria, iv)
      return
      end
c_______________________________________________________________________
	subroutine nneigh
	implicit double precision (a-h,o-z)
c  From an arbitrarily ordered input list of unit vectors, the
c  program creates for each one a list of neighbors
c  such that they form a convex polygon around that vertex and the
c  polygon encloses no other vertex.  These polygons taken 
c  together cover the sphere exactly three times.  The polygon
c  forms the support region for an interpolation function on the sphere
c  centered on the associated vertex.
c
c$$$$ Requires: isort, order
c
c  The parameter nxyz is the maximum number of unit vectors allowed.
c  The parameter nvert is the maximum number of faces at any vertex.
      parameter (nxyz=6500, nxyz2=2*nxyz, nxyz6=6*nxyz,nvert=10)
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
      dimension neigh(100)
      common /io/inp,iout,iprint
      common /bound/ dmin,dmax,amin,amax
c
      dmin=1.0
      dmax=-1.0
	  write(iout,*)' Entering nneigh'
c
c  For all unit vectors discover the neighbor set.
c
      do 1500 j=1,n
        near(1,j)=0
c  Scan through the triangle list, seeking the vertex  j.
        nay=0
        do 1400 jv=1, ntria
          if (iv(1,jv).eq.j .or. iv(2,jv).eq.j .or. iv(3,jv).eq.j) then
          do 1250 i=1, 3
            ivj=iv(i,jv)
            if (ivj.eq.j .or. ivj.gt.n) goto 1250
            nay=1+nay
            neigh(nay)=ivj
 1250     continue
          endif
 1400   continue
c  Sort the list, then discard duplicate vertices.
        call isort(nay, neigh)
        nn=1
        do 1450 i=2, nay
          if (neigh(i) .eq. neigh(i-1)) goto 1450
          nn=nn +1
          if (nn .gt. nvert) write(iout, '(a,i4,a,i4)')
     $    ' WARNING: More than ',nvert,' neighbors at vertex ',j
	  neigh(nn)=neigh(i)
 1450   continue
c  Rank the list by longitude about the point p(.,j), so that vertices
c  progress sequentially about centre. neigh(1) is assumed zero
c  longitude
	call order (j,nn,neigh)
c  
        if (nn .lt. nvert) near(nn+1,j)=0
c        write(iout, '(i5,(4x,14i5/))')j,(near(i,j),i=1,nn)
c
c
 1500 continue
c  check completeness of covering polygons by computing total area.
c  should be 12pi
c
	sum=0.0
	do 3500 j=1,n
	   do 3100 i=2,nvert
	     if(near(i,j).eq.0)goto 3200
	     sum=sum + sarea(p(1,j), p(1,near(i-1,j)), p(1,near(i,j)))
	     k=i
3100       continue
3200	   sum=sum + sarea(p(1,j), p(1,near(1,j)), p(1,near(k,j)))
3500    continue
	cov = sum/(12*3.1415926)
	write(iout, '(/a,f7.4)') 'Covering area /12*pi =',cov
	if  (abs(cov-1.0).gt.0.02) 
     $   write(iout,*)'covering is inconsistent'
c
      return
      end
c_______________________________________________________________________
      subroutine isort(no, x)
	implicit double precision (a-h,o-z)
c$$$$$ calls no other routines
c  In-place rapid sorter.  Rearranges  x  into ascending order
      integer x(no)
c
      mo=no
 2    if (mo-15) 21,21,23
 21   if (mo.le.1) return
      mo=2*(mo/4) + 1
      go to 24
 23   mo=2*(mo/8) + 1
 24   ko=no - mo
      jo=1
 25   i=jo
 26   if (x(i) - x(i+mo)) 28,28,27
 27   temp=x(i)
      x(i)=x(i+mo)
      x(i+mo)=temp
      i=i - mo
      if (i-1) 28,26,26
 28   jo=1 + jo
      if (jo-ko) 25,25,2
      end                                                               xsort
c
c_______________________________________________________________________
	subroutine order(j,nn,neigh)
	implicit double precision (a-h,o-z)
c
c  orders the nn nearest neighbours of point j, ranking by longitude
c  from neigh(nn)
c
c*****requires sort2, gpole
c
	parameter (nxyz=6500, nxyz2=2*nxyz, nxyz6=6*nxyz,nvert=10)
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
	dimension x(3),phi(nvert,2),neigh(nn)
	data tpi/6.28318531/
c
c  first rotate so that p(.,j) is the pole, zero longitude is given by
c  neigh(1) and compute longitude phi for each neighbour.
	do 50 i=1,nn
	   call gpole(1,p(1,j),p(1,neigh(1)),p(1,neigh(i)),x)
	   phi(i,2)=float(i)
	   phi(i,1)=atan2(x(2),x(1))
	   phi(i,1)=dmod(phi(i,1)+tpi,tpi)
50	continue
c  rank neigh(i) by phi(i,1)
	call sort2(nvert,nn,phi)
	do 60 i=1,nn
	  near(i,j)=neigh(nint(phi(i,2)))
60	continue
	return
	end
c
c_______________________________________________________________________
      subroutine sort2(nd,no, x)
	implicit double precision (a-h,o-z)
c$$$$$ calls no other routines
c  In-place rapid sorter.  Rearranges  x  into ascending order of first
c  column
      dimension x(nd,2)
c
      mo=no
 2    if (mo-15) 21,21,23
 21   if (mo.le.1) return
      mo=2*(mo/4) + 1
      go to 24
 23   mo=2*(mo/8) + 1
 24   ko=no - mo
      jo=1
 25   i=jo
 26   if (x(i,1) - x(i+mo,1)) 28,28,27
 27   temp=x(i,1)
	temp2=x(i,2)
      x(i,1)=x(i+mo,1)
      x(i,2)=x(i+mo,2)
      x(i+mo,1)=temp
      x(i+mo,2)=temp2
      i=i - mo
      if (i-1) 28,26,26
 28   jo=1 + jo
      if (jo-ko) 25,25,2
      end                                                               xsort
c
c_______________________________________________________________________
      function sarea(a, b, c)
	implicit double precision (a-h,o-z)
c$$$$ calls sgn
c  Given three unit vectors a, b, c  finds the area of the
c  SPHERICAL triangle formed at the tips.
      dimension a(3),b(3),c(3)
	 common /io/ inp,iout,iprint
c
c	write(iout,*)'Entering sarea'
      cosa=b(1)*c(1)+b(2)*c(2)+b(3)*c(3)
      cosb=a(1)*c(1)+a(2)*c(2)+a(3)*c(3)
      cosc=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      sina=sqrt(abs(1.0 - cosa**2))
      sinb=sqrt(abs(1.0 - cosb**2))
      sinc=sqrt(abs(1.0 - cosc**2))
      AA=((cosa-cosb*cosc)/(sinb*sinc))
	if(abs(aa).gt.1.0)aa=sgn(aa)
      BB=((cosb-cosa*cosc)/(sina*sinc))
	if(abs(bb).gt.1.0)bb=sgn(bb)
      CC=((cosc-cosb*cosa)/(sinb*sina))
	if(abs(cc).gt.1.0)cc=sgn(cc)
	aa=acos(aa)
	bb=acos(bb)
	cc=acos(cc)
      sarea=AA+BB+CC - 3.14159265358979324
c	write(iout, '(9f8.4/a,f8.4)')a,b,c,' Area =',sarea
      return
      end
c
c_______________________________________________________________________
	function sgn(x)
	implicit double precision (a-h,o-z)
	if(x.le.0.)then
	  sgn=-1.0
	else
	  sgn=1.0
	endif
	return
	end
      subroutine triq(nx, ny)
	implicit double precision (a-h,o-z)
c$$$$ calls no other routines
c  Constructs the conical product quadrature formula for integration
c  over the plane triangular region  0.le.x, 0.le.y, x+y.le.1. see
c  Stroud,1971, 'Approximate Calculation of Multiple Integrals',pp28-30.
c  nx  is the number of sample points desired along  x  lines
c  ny  the number along  y  lines
c  nquad   the number of sample points generated.  ideally this is
c      nx*ny, but may differ from this because of the limited
c      number of gauss formulas stored.
c xyw  triples of  x-y coordintes and weights generated for the
c      quadrature formula.  if n=min(nx,ny), the formula is
c      exact for polynomials in x**p*y**q  where p+q.le.2*n-1.
c  if  nx  is negative the routine supplies the special 7-point degree 5
c  formula of Stroud (*T2 5-1, page 315).
c
c  The results are loaded into the common /rule/
c
      common /rule/ xyw(3,100),nquad
      common /io/ inp,iout,iprint
      dimension ixp(15),xw(316),iyp(9),yw(72),zw(21)
c
c  Gauss formulas with w(x)=1, with orders 2,3,...10,12,16,20,24,32.
      data ixp/1,5,11,19,29,41,55,71,89,109,133,165,205,253,317/
      data (xw(j),j=1,54)/
     $ 0.2113248654, 0.5000000000, 0.7886751346, 0.5000000000,            ord 2
     $ 0.1127016654, 0.2777777778, 0.5000000000, 0.4444444444,            ord 3
     $ 0.8872983346, 0.2777777778,
     $ 0.0694318442, 0.1739274225, 0.3300094782, 0.3260725774,            ord 4
     $ 0.6699905218, 0.3260725774, 0.9305681558, 0.1739274225,
     $ 0.0469100771, 0.1184634425, 0.2307653450, 0.2393143352,            ord 5
     $ 0.5000000000, 0.2844444444, 0.7692346550, 0.2393143352,
     $ 0.9530899229, 0.1184634425,
     $ 0.0337652429, 0.0856622462, 0.1693953068, 0.1803807865,           ord 6
     $ 0.3806904070, 0.2339569673, 0.6193095930, 0.2339569673,
     $ 0.8306046932, 0.1803807865, 0.9662347571, 0.0856622462,
     $ 0.0254460439, 0.0647424831, 0.1292344072, 0.1398526957,           ord 7
     $ 0.2970774243, 0.1909150252, 0.5000000000, 0.2089795918,
     $ 0.7029225757, 0.1909150252, 0.8707655928, 0.1398526957,
     $ 0.9745539561, 0.0647424831/
      data (xw(j),j=55,108)/
     $ 0.0198550718, 0.0506142681, 0.1016667613, 0.1111905172,           ord 8
     $ 0.2372337951, 0.1568533229, 0.4082826788, 0.1813418917,
     $ 0.5917173212, 0.1813418917, 0.7627662049, 0.1568533229,
     $ 0.8983332387, 0.1111905172, 0.9801449282, 0.0506142681,
     $ 0.0159198803, 0.0406371942, 0.0819844464, 0.0903240803,           ord 9
     $ 0.1933142837, 0.1303053482, 0.3378732883, 0.1561735385,
     $ 0.5000000000, 0.1651196775, 0.6621267117, 0.1561735385,
     $ 0.8066857163, 0.1303053482, 0.9180155536, 0.0903240803,
     $ 0.9840801197, 0.0406371942,
     $ 0.0130467358, 0.0333356722, 0.0674683167, 0.0747256746,           ord 10
     $ 0.1602952159, 0.1095431812, 0.2833023030, 0.1346333596,
     $ 0.4255628305, 0.1477621123, 0.5744371695, 0.1477621123,
     $ 0.7166976970, 0.1346333596, 0.8397047841, 0.1095431812,
     $ 0.9325316833, 0.0747256746, 0.9869532642, 0.0333356722/
      data (xw(j),j=109,164)/
     $ 0.0092196830, 0.0235876680, 0.0479413720, 0.0534696630,           ord 12
     $ 0.1150486630, 0.0800391645, 0.2063410230, 0.1015837135,
     $ 0.3160842505, 0.1167462685, 0.4373832955, 0.1245735230,
     $ 0.5626167045, 0.1245735230, 0.6839157495, 0.1167462685,
     $ 0.7936589770, 0.1015837135, 0.8849513370, 0.0800391645,
     $ 0.9520586280, 0.0534696630, 0.9907803170, 0.0235876680,
     $ 0.0052995328, 0.0135762297, 0.0277124885, 0.0311267620,           ord 16
     $ 0.0671843988, 0.0475792558, 0.1222977958, 0.0623144856,
     $ 0.1910618778, 0.0747979944, 0.2709916112, 0.0845782597,
     $ 0.3591982246, 0.0913017075, 0.4524937451, 0.0947253052,
     $ 0.5475062549, 0.0947253052, 0.6408017754, 0.0913017075,
     $ 0.7290083888, 0.0845782597, 0.8089381222, 0.0747979944,
     $ 0.8777022042, 0.0623144856, 0.9328156012, 0.0475792558,
     $ 0.9722875115, 0.0311267620, 0.9947004672, 0.0135762297/
      data (xw(j),j=165,240)/
     $ 0.0034357004, 0.0088070035, 0.0180140364, 0.0203007149,           ord 20
     $ 0.0438827859, 0.0313360241, 0.0804415141, 0.0416383708,
     $ 0.1268340468, 0.0509650599, 0.1819731597, 0.0590972660,
     $ 0.2445664990, 0.0658443192, 0.3131469557, 0.0710480546,
     $ 0.3861070745, 0.0745864933, 0.4617367395, 0.0763766935,
     $ 0.5382632605, 0.0763766935, 0.6138929255, 0.0745864933,
     $ 0.6868530443, 0.0710480546, 0.7554335010, 0.0658443192,
     $ 0.8180268403, 0.0590972660, 0.8731659532, 0.0509650599,
     $ 0.9195584859, 0.0416383708, 0.9561172141, 0.0313360241,
     $ 0.9819859636, 0.0203007149, 0.9965642996, 0.0088070035,
     $ 0.0024063900, 0.0061706150, 0.0126357220, 0.0142656945,           ord 24
     $ 0.0308627240, 0.0221387195, 0.0567922365, 0.0296492925,
     $ 0.0899990070, 0.0366732405, 0.1299379040, 0.0430950807,
     $ 0.1759531740, 0.0488093260, 0.2272892645, 0.0537221350,
     $ 0.2831032460, 0.0577528340, 0.3424786600, 0.0608352365,
     $ 0.4044405663, 0.0629187280, 0.4679715535, 0.0639690975,
     $ 0.5320284465, 0.0639690975, 0.5955594337, 0.0629187280,
     $ 0.6575213400, 0.0608352365, 0.7168967540, 0.0577528340,
     $ 0.7727107355, 0.0537221350, 0.8240468260, 0.0488093260/
      data (xw(j),j=241,316)/
     $ 0.8700620960, 0.0430950807, 0.9100009930, 0.0366732405,
     $ 0.9432077635, 0.0296492925, 0.9691372760, 0.0221387195,
     $ 0.9873642780, 0.0142656945, 0.9975936100, 0.0061706150,
     $ 0.0013680695, 0.0035093050, 0.0071942445, 0.0081371970,           ord 32
     $ 0.0176188725, 0.0126960325, 0.0325469625, 0.0171369310,
     $ 0.0518394225, 0.0214179490, 0.0753161935, 0.0254990295,
     $ 0.1027581025, 0.0293420465, 0.1339088910, 0.0329111110,
     $ 0.1684778670, 0.0361728970, 0.2061421215, 0.0390969475,
     $ 0.2465500460, 0.0416559620, 0.2893243620, 0.0438260465,
     $ 0.3340656990, 0.0455869390, 0.3803563190, 0.0469221995,
     $ 0.4277640195, 0.0478193600, 0.4758461675, 0.0482700440,
     $ 0.5241538325, 0.0482700440, 0.5722359805, 0.0478193600,
     $ 0.6196436810, 0.0469221995, 0.6659343010, 0.0455869390,
     $ 0.7106756380, 0.0438260465, 0.7534499540, 0.0416559620,
     $ 0.7938578785, 0.0390969475, 0.8315221330, 0.0361728970,
     $ 0.8660911090, 0.0329111110, 0.8972418975, 0.0293420465,
     $ 0.9246838065, 0.0254990295, 0.9481605775, 0.0214179490,
     $ 0.9674530375, 0.0171369310, 0.9823811275, 0.0126960325,
     $ 0.9928057555, 0.0081371970, 0.9986319305, 0.0035093050/
c
c  Gauss formulas w(x)=x, with orders 1,2 ... 8.
      data iyp/1,3,7,13,21,31,43,57,73/
      data yw/
     $0.6666666667, 0.5,                                                 ord 1
     $.3550510257, .1819586183,.8449489743,.3180413817,                  ord 2
     $.2123405382,.0698269799,.5905331356,.2292411064,                   ord 3
     $.9114120405,.2009319137,                                           ord 3
     $.1397598643,.0311809710,.4164095676,.1298475476,                   ord 4
     $.7231569864,.2034645680,.9428958039,.1355069134,                   ord 4
     $.0985350858,.0157479145,.3045357266,.0739088701,.5620251898,       ord 5
     $.1463869871,.8019865821,.1671746381,.9601901429,.0967815902,       ord 5
     $.0730543287,.0087383018,.2307661380,.0439551656,.4413284812,       ord 6
     $.0986611509,.6630153097,.1407925538,.8519214003,.1355424972,       ord 6
     $.9706835728,.0723103307,                                           ord 6
     $.0562625605,.0052143622,.1802406917,.0274083567,.3526247171,       ord 7
     $.0663846965,.5471536263,.1071250657,.7342101772,.1273908973,       ord 7
     $.8853209468,.1105092582,.9775206136,.0559673634,                   ord 7
     $.0446339553,.0032951914,.1443662570,.0178429027,.2868247571,       ord 8
     $.0454393195,.4548133152,.0791995995,.6280678354,.1060473594,       ord 8
     $.7856915206,.1125057995,.9086763921,.0911190236,                   ord 8
     $.9822200849,.0445508044/                                           ord 8
c  Special degree 5 formula.
      data zw/0.333333333, 0.333333333, 0.1125,
     $0.101286507, 0.101286507, 0.062969590,
     $0.797426985, 0.101286507, 0.062969590,
     $0.101286507, 0.797426985, 0.062969590,
     $0.470142064, 0.470142064, 0.066197076,
     $0.059715871, 0.470142064, 0.066197076,
     $0.470142064, 0.059715871, 0.066197076/
c
c
c  Find sample numbers closest to those input
      if (nx.lt.0) goto 2000
      ntwo=2*nx
      do 1100 mx=1, 14
 1100 if (ixp(mx+1)-ixp(mx).ge.ntwo) goto 1110
      mx=14
 1110 my=min0(ny, 8)
c
      mx1=ixp(mx)
      mx2=ixp(mx+1) - 1
      my1=iyp(my)
      my2=iyp(my+1) - 1
c
c  Create the product formula.
      m=0
      do 1500 i=mx1, mx2, 2
      do 1500 j=my1, my2, 2
      m=1 + m
      xyw(1,m)=1.0 - yw(j)
      xyw(2,m)=yw(j)*xw(i)
      xyw(3,m)=yw(j+1)*xw(i+1)
 1500 continue
      nquad=m
      if (iprint .ge. 3) write(iout,150) ((xyw(i,j),i=1,3), j=1,nquad)
 150  format(/' Quadrature coordinates and weights'/(3f8.4,3x,3f8.4))
      return
c  Copy special formula into output array.
 2000 m=7
      do 2100 j=1, 7
      do 2100 i=1,3
 2100 xyw(i,j)=zw(3*(j-1)+i)
      nquad=m
      if (iprint .ge. 3) write(iout,150) ((xyw(i,j),i=1,3), j=1,nquad)
      return
      end
c_______________________________________________________________________
c_______________________________________________________________________
	subroutine gnomon(iway,u,r)
	implicit double precision (a-h,o-z)
c  if iway=1 performs the gnomonic projection of r(x,y,z) on the sphere
c  into u on the plane tangent at the north pole
c  if iway=-1 the inverse transformation is performed
c  All coordinates are cartesian
	 common /io/inp,iout,iprint
	dimension u(2),r(3)
	if (iway.eq.1) then
	   if (r(3).eq.0)then
	     write(iout,*)'gnomonic mapping undefined'
	     stop
	   endif
	   rho=sqrt(r(1)**2 + r(2)**2)/r(3)
	   psi=atan2(r(2),r(1))
	   u(1)=rho*cos(psi)
	   u(2)=rho*sin(psi)
	else if(iway.eq.-1) then
	   theta=atan(sqrt(u(1)**2 + u(2)**2))
	   psi=atan2(u(2),u(1))
	   r(1)=sin(theta)*cos(psi)
	   r(2)=sin(theta)*sin(psi)
	   r(3)=cos(theta)
	endif
	return
	end
c______________________________________________________________________
	subroutine gpole(iway,v,a,x1,x2)
	implicit double precision (a-h,o-z)
c  if iway = 1
c  converts unit vector x1 from a geographic Cartesian coord system to 
c  one whose N-pole is given by v, a is any arbitrary vector in the
c  geographic system and converts to a unit vector which is at zero 
c  longitude and correct latitude in the new system
c  
c  if iway=-1 the inverse transformation from x2 to x1 is performed
c  v,a are always expressed in the geographic system
c
c
	dimension v(3),a(3),x1(3),x2(3),vxa(3),vc(3)
c
	cth=dotty(a,v)
	if (iway.eq.1)then
c  new z
  	   x2(3)=dotty(v,x1)
c  compute new y direction
	   call axb(v,a,vxa)
	   call unit(vxa,vxa)
	   x2(2)=dotty(vxa,x1)
c  compute new x
	   do 1 i=1,3
1	   vc(i)=v(i)
	   call scaley(cth,vc)
	   do 2 i=1,3
2	   vc(i)=a(i)-vc(i)
	   call unit(vc,vc)
	   x2(1)=dotty(vc,x1)
c
	else if (iway.eq.-1) then
	   call axb(v,a,vxa)
	   call unit(vxa,vxa)
	   do 3 i=1,3
3	   vc(i)=v(i)
	   call scaley(cth,vc)
	   do 4 i=1,3
4	   vc(i)=a(i)-vc(i)
	   call unit(vc,vc)
	   x1(1)=vc(1)*x2(1) + vxa(1)*x2(2) + v(1)*x2(3)
	   x1(2)=vc(2)*x2(1) + vxa(2)*x2(2) + v(2)*x2(3)
	   x1(3)=vc(3)*x2(1) + vxa(3)*x2(2) + v(3)*x2(3)
	endif
	return
	end
c_______________________________________________________________________
	subroutine axb(a,b,c)
	implicit double precision (a-h,o-z)
c  returns in c the vector cross product of 3-vectors a and b
	dimension a(3),b(3),c(3)
	c(1)=a(2)*b(3) -a(3)*b(2)
	c(2)=a(3)*b(1)-a(1)*b(3)
	c(3)=a(1)*b(2)-a(2)*b(1)
	return
	end
c_______________________________________________________________________
	subroutine unit(x,ux)
	implicit double precision (a-h,o-z)
c  converts x into a unit vector ux
	dimension x(3),ux(3)
	unorm=0.
	do 1 i=1,3
1	unorm= unorm +x(i)**2
	unorm =sqrt(unorm)
	do 2 i=1,3
2	ux(i)= x(i)/unorm
	return
	end
c_______________________________________________________________________
	subroutine scaley(c,x)
	implicit double precision (a-h,o-z)
c  scales a vector x by the constant c
	dimension x(3)
	do 1 i=1,3
1	x(i)=x(i)*c
	return
	end
c_______________________________________________________________________
	subroutine centr(x1,x2,x3,c)
	implicit double precision (a-h,o-z)
	dimension x1(3),x2(3),x3(3),c(3)
c  computes centroid of triangle vertices x1,x2,x3
	do 11 j=1,3
	c(j)= x1(j) + x2(j) + x3(j)
11	continue
	cnorm =sqrt(dotty(c,c))
	do 12 j=1,3
	c(j)=c(j)/cnorm
12	continue
	return
	end
c______________________________________________________________________
      subroutine bfield
	implicit double precision (a-h,o-z)
c
c  A sphere is defined by spherical triangular facets with
c  corners at the points defined by the 3-vectors  p  and a pointer
c  array  iv  specifying the corners to be connected with  ntria
c  triangles. The 4th component of p specifies the radial field at the
c  core mantle boundary
c  This routine computes the magnetic field at any point
c  outside the core by using the gram matrix computed in gram
c  (this calculation is assumed performed already).
c
c
	parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,nvert=10)
	parameter (maxob=7000,maxp=30,maxe=2)
	common /io/inp,iout,iprint
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
	common /observ/ro(3,maxob,maxe),bro(maxob,maxe),sigma(maxob,
     $         maxe),iflag(maxob,maxe),mobs,nobs(maxe),nepoch,lobs(maxe)
	common /calcul/bs(maxob,maxe),a(maxob,nxyz),bc(nxyz2),
     $                 reg(nxyz,nxyz)
	common /rule/xyw(3,100),nquad
	common /circle/v(4,nxyz2),istack(nxyz2)
c
c  Zero the field vector
	l=0
	do 1150 m=1,nepoch
      do 1100 i=1, nobs(m)
 1200       bs(i,m)=0.0
	if (iflag(i,m).gt.0)then
	  l=l+1
	  do 1400 k=1,n
	  bs(i,m) =bs(i,m) + a(l,k+(m-1)*n)*bc((m-1)*n+k)*
     $                sigma(i,m)
1400	  continue
	else
	  bs(i,m)=99999.
	endif
1300	continue
 1100 continue
	 l=l+1
1150  continue
      return
      end
c_______________________________________________________________________
      subroutine gram
	implicit double precision (a-h,o-z)
c$$$$$ calls bgreen and bint, to compute green functions and for
c  interpolation of b on core surface
c  For each triangle on the core,
c  performs the numerical cubature of three surface integrals for br,
c  btheta, bphi
c  This provides the matrix a for the computation of the observed
c  surface fields, given the field at source points p(i,j)
c  core field br at vertices p, is in 4th component of arrays
c
c  The quadrature rule generated by 'triq' for integrals over
c  plane triangle with vertices  (0,0,0) (1,0,0) (0,1,0) is assumed
c  to be available in common /rule/.
c
c  The list of observer coordinates is in /observ/.
c
c  The results are returned in common /calcul/
c
	parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,nvert=10)
	parameter (maxob=7000,maxp=30,maxe=2)
	common /io/inp,iout,iprint
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
	common /observ/ro(3,maxob,maxe),bro(maxob,maxe),sigma(maxob,
     $	             maxe),iflag(maxob,maxe), mobs,nobs(maxe),nepoch,
     $               lobs(maxe)
	common /calcul/bs(maxob,maxe),a(maxob,nxyz),bc(nxyz2),
     $                 reg(nxyz,nxyz)
	common /rule/xyw(3,100),nquad
	common /circle/v(4,nxyz2),istack(nxyz2)
      dimension rp(3),g(3),vc(4),aa(3),b(3),c(3),
     $          ua(2),ub(2),uc(2),rg(2),rpa(3),gamma(3)
c
	write(iout,*)'entering gram'
	xxx=0.
	nn=n*nepoch
c  initialize matrix to zero
	do 1900 j=1,nn
  	  do 1700 iob=1,mobs+3
	  a(iob,j)=0.
1700	  continue
1900	continue
	do 1899 i=1,nepoch
	  write(iout,*)' Epoch ',i,' has ',nobs(i), ' observations'
1899	continue
c
	do 1500 j=1,ntria
c  
c  vertices of triangle are points iv1,iv2,iv3
	iv1=iv(1,j)
	iv2=iv(2,j)
	iv3=iv(3,j)
c
c  compute centroid vc of current triangle
c
	call centr(p(1,iv1),p(1,iv2),p(1,iv3),vc)
c	xxx=xxx + area(p(1,iv1),p(1,iv2),p(1,iv3))
c
c  rotate so vertices are in coord system with vc as pole, p(1,iv1) at
c  zero longitude
	call gpole(1,vc,p(1,iv1),p(1,iv1),aa)
	call gpole(1,vc,p(1,iv1),p(1,iv2),b)
	call gpole(1,vc,p(1,iv1),p(1,iv3),c)
c
c  find gnomonic projection of vertices and jacobian of transformation
	call gnomon(1,ua,aa)
	call gnomon(1,ub,b)
	call gnomon(1,uc,c)
	ajac = parea (ua,ub,uc)*2.
c
c  Perform facet integration for all observation points.
c  Run through the cubature rule for this triangle for
c  each observation station.
c
        do 2400 k=1, nquad
c
c  find points on gnomonic surface corresponding to integration points
c
	rg(1)= ua(1) + (ub(1)-ua(1))*xyw(1,k)  +(uc(1)-ua(1))*xyw(2,k)
	rg(2)= ua(2) + (ub(2)-ua(2))*xyw(1,k)  +(uc(2)-ua(2))*xyw(2,k)
c
c  Find interpolant 
c
	gamma(1)=1. -xyw(1,k)-xyw(2,k)
	gamma(2)=xyw(1,k)
	gamma(3)=xyw(2,k)
c
c  find points in physical space corres. to rg
c
	call gnomon(-1,rg,rpa)
	call gpole(-1,vc,p(1,iv1),rp,rpa)
c
c  Find Jacobian for gnomonic to spherical projection
c
	bjac= sqrt((1+rg(1)**2 + rg(2)**2)**-3)
c
c  evaluate green functions at rp for all observer points
c
      ll=0
      do 2550 i=1,nepoch
	lobs(i)=0
      do 2500 iob=1, nobs(i)
	call bgreen(ro(1,iob,i),rp,g)
c
c  Weight integral according to xyw rule
c
	if(iflag(iob,i).gt.0)then
	  lobs(i)=lobs(i)+1
	  ll=ll+1
	  xx=elt(iflag(iob,i),g,bro(iob,i))
  	  a(ll,iv1+(i-1)*n)=a(ll,iv1+(i-1)*n) + 
     $     xx*gamma(1)*xyw(3,k)*ajac*bjac/sigma(iob,i)
	  a(ll,iv2+(i-1)*n)=a(ll,iv2+(i-1)*n) + 
     $     xx*gamma(2)*xyw(3,k)*ajac*bjac/sigma(iob,i)
	  a(ll,iv3+(i-1)*n)=a(ll,iv3+(i-1)*n) + 
     $     xx*gamma(3)*xyw(3,k)*ajac*bjac/sigma(iob,i)
	endif
 2500 continue
	  ll=ll+1
c  No monopoles allowed
  	  a(ll,iv1+(i-1)*n)=a(ll,iv1+(i-1)*n) + 
     $     gamma(1)*xyw(3,k)*ajac*bjac
	  a(ll,iv2+(i-1)*n)=a(ll,iv2+(i-1)*n) + 
     $     gamma(2)*xyw(3,k)*ajac*bjac
	  a(ll,iv3+(i-1)*n)=a(ll,iv3+(i-1)*n) + 
     $     gamma(3)*xyw(3,k)*ajac*bjac
2550  continue
 2400   continue
1500	continue
c	write(iout,*)'sum of s. triangle areas',xxx
      return
      end
c______________________________________________________________________
	subroutine bgreen(ro,so,g)
	implicit double precision (a-h,o-z)
c  computes greens function for magnetic field components br,bthet, bphi
c  at observer point ro due to source point so.
c  ro specified as r, theta, phi coords
c  so specified as x, y, z , geocentric.
c
	dimension ro(3),so(3),g(3)
	data fpi/12.56637062/,c/1.0/
	t=ro(2)
	p=ro(3)
	sounit=sqrt(so(1)**2 + so(2)**2 + so(3)**2)
	do 1 i=1,3
1	so(i)=so(i)/sounit
	rho=c/ro(1)
	rhosq=rho*rho
	rhocu=rhosq*rho
c  compute runit.runit'
	amu=sin(t)*(cos(p)*so(1) + sin(p)*so(2)) + cos(t)*so(3)
	r=sqrt(1.0-2.0*amu*rho +rhosq)
	rrr=r**3
	tt=1+r-amu*rho
c  compute runit'.thetaunit and runit'.phiunit
	uthet=cos(t)*(cos(p)*so(1) + sin(p)*so(2)) - sin(t)*so(3)
c	
	uphi=-sin(p)*so(1) + cos(p)*so(2)
c
c  compute g(1),g(2),g(3), green functions for r,theta,phi
	g(1)=rhosq*((1.0-rhosq)/rrr -1.)/fpi
	x=rhocu*(1.0 + 2.0*r-rhosq)/(rrr*tt*fpi)
	g(2)=-uthet*x
	g(3)=-uphi*x
	return
	end
c_______________________________________________________________________
	function parea(a,b,c)
	implicit double precision (a-h,o-z)
c  computes area of planar triangle
	dimension a(2),b(2),c(2),ab(3),ac(3),ar(3)
	do 5 j=1,2
	  ab(j)=b(j)-a(j)
	  ac(j)=c(j)-a(j)
5	continue
	 ab(3)=0.
	ac(3)=0.
	call axb(ab,ac,ar)
	parea=sqrt(dotty(ar,ar))/2.
	return
	end
c______________________________________________________________________
	function elt(iflag,g,bro)
	implicit double precision (a-h,o-z)
c  converts green function for br, btheta, bphi given in g to other
c  elements fo the geomagnetic field
c	iflag=1  declination
c	      2  inclination
c	      3  horizontal intensity
c	      4  x or north component
c	      5  y or east component
c	      6  z or vertical down component
c	      7  f total intensity
c	
	common /io/inp,iout,iprint
	dimension g(3)
	data drad/57.2958/
c	write(iout,*)'iflag=',iflag
	if(iflag.eq.1)then
	  elt= g(2)*sin(bro/drad) + g(3)*cos(bro/drad)
	else if(iflag.eq.2)then
	  h=sqrt(g(2)**2 + g(3)**2)
	  elt= drad*atan2(-g(1),h)
	else if(iflag.eq.3)then
	  elt=sqrt(g(2)**2 + g(3)**2)
	else if(iflag.eq.4)then
	  elt=  -g(2)
	else if(iflag.eq.5)then
	  elt=  g(3)
	else if(iflag.eq.6)then
	  elt=  -g(1)
	else if(iflag.eq.7)then
	  elt=  sqrt(g(1)**2 + g(2)**2 + g(3)**2)
	else
	  write(iout,*)'Unknown measurement type'
	endif
	return
	end
c--------------------------------------------------------------
      subroutine flux(i)
	implicit double precision (a-h,o-z)
c
c  Computes the flux through the null-flux curves on the core,
c  and the gradient matrix for the flux through the curves with
c  respect to changes in the values of the field at the vertices
c  on the core.
c
c  Assumes patch has already been called to define the variables
c  in /patchy/.
c
c
c  Outline:
c
c  Initialize the flux vector and flux gradient matrix to zero
c  Loop in triangles
c    for each triangle, compute rotation and gnomonic projection
c    loop in integration sample points
c       compute preimage of the sample point on the core
c       compute the jacobian for the gnomonic at that sample point
c       compute field value at that point
c       add the value, weighted by the int. wt. and the jacobians,
c           to the positive patch or negative patch, according
c           to its sign.
c       add interp. wt*int. wt.* jacobian to the appropriate elements
c           of the grad. matrix (according to signs of br and vertex)
c    next integration sample point
c  next triangle.
c
      parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,maxp=30,
     $   nvert=10,maxe=2)
      dimension rp(3),g(3),vc(4),aa(3),b(3),c(3),
     +   ua(2),ub(2),uc(2),rg(2),rpa(3)
      common/patchy/np(maxe),npts(maxp,maxe),isign(maxp,maxe),
     $  ipts(maxp,maxe),lpts(2*nxyz,maxe),ntris(nxyz),
     $  itris(nxyz+1),
     $  ltris(2*(2*nxyz-4)),lptr(2*nxyz-4,2),
     $   llpt(maxp,maxe),npatch
      common/patchg/ fluxv(maxp,maxe),gradm(maxp,nxyz,maxe)
c
c  see subroutine patch for definitions of the variables in /patchy/
c
      common/tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
      common /io/inp,iout,iprint
      common /rule/xyw(3,100),nquad
      common /circle/v(4,nxyz2),istack(nxyz2)
c
c  initialize flux vector and gradient matrix.
c******
      write(iout,*)' arrived in flux '
      do 10 j=1,np(i)
      fluxv(j,i)=0.
      do 10 k=1,n
 10   gradm(j,k,i)=0.
c
c  Loop in triangles.
      do 100 j=1,ntria
c  
c
c  compute centroid vc of current triangle
         call centr(p(1,iv(1,j)),p(1,iv(2,j)),p(1,iv(3,j)),vc)
c
c  rotate so vertices are in coord system with vc as pole, p(1,iv1) at
c  zero longitude
         call gpole(1,vc,p(1,iv(1,j)),p(1,iv(1,j)),aa)
         call gpole(1,vc,p(1,iv(1,j)),p(1,iv(2,j)),b)
         call gpole(1,vc,p(1,iv(1,j)),p(1,iv(3,j)),c)
c
c  find gnomonic projection of vertices and jacobian of transformation
         call gnomon(1,ua,aa)
         call gnomon(1,ub,b)
         call gnomon(1,uc,c)
         ajac = parea (ua,ub,uc)*2.
c*******
c      write(iout,*)' entering integration loop '
c
c  Loop in integration sample points.
         do 50 k=1, nquad
c
c  find points on gnomonic surface corresponding to integration points
            rg(1)= ua(1) + (ub(1)-ua(1))*xyw(1,k)  +
     $        (uc(1)-ua(1))*xyw(2,k)
            rg(2)= ua(2) + (ub(2)-ua(2))*xyw(1,k)  +
     $        (uc(2)-ua(2))*xyw(2,k)
c
c  Find Jacobian for gnomonic to spherical projection
            bjac= sqrt((1+rg(1)**2 + rg(2)**2)**-3)
c
c  Find interpolating vector g(.)---linear interpolation in gnomonic
c  plane.
            g(1)=1.-xyw(1,k)-xyw(2,k)
            g(2)=xyw(1,k)
            g(3)=xyw(2,k)
c
c  find the point rp on the core corresponding to rg.
            call gnomon(-1,rg,rpa)
            call gpole(-1,vc,p(1,iv(1,j)),rp,rpa)
c
c  Find field value and its sign.
            br=0.
            do 30 l=1,3
 30         br=br+g(l)*p(4,iv(l,j))
            if (br.ge.0.) then
               ndx=1
               sgnbr=1.
            else
               ndx=2
               sgnbr=-1.
            endif
c
c  Update the flux vector fluxv and the gradient matrix gradm.
c*******
c      write(iout,*)' tri ',j,' sample ', k,' br ',br, ' g ',
c     +  (g(ii),ii=1,3),' ndx ',ndx,' ajac ',ajac,' bjac ',bjac
c*******
            fluxv(lptr(j,ndx),i)=fluxv(lptr(j,ndx),i)+
     $        ajac*bjac*xyw(3,k)*br
            do 40 l=1,3
c               if(signm(iv(l,j)).eq.sgnbr) 
              gradm(lptr(j,ndx),iv(l,j),i)=gradm(lptr(j,ndx),
     $         iv(l,j),i)+g(l)*ajac*bjac*xyw(3,k)
 40         continue
c
c  Next integration sample point.
 50      continue
c
c  Next triangle.
100   continue
      return
      end
c-------------------------------------------------------------
      subroutine patch(i)
	implicit double precision (a-h,o-z)
c
c  constructs lists of triangles and vertices associated with
c  each patch on the core on which Br has one sign
c  
c  constructs a list of the signs of the patches.
c  constructs an array indicating which points are in which
c  patches, and how many patches a point appears in.
c
c  assumes the  near  array has been previously constructed.
c
c* calls signm
c
      parameter (nxyz=6500, nvert=10, nxyz2=2*nxyz)
	parameter(maxp=30,maxe=2)
		common /io/inp,iout,iprint
      common/patchy/np(maxe),npts(maxp,maxe),isign(maxp,maxe),
     $   ipts(maxp,maxe),
     $   lpts(2*nxyz,maxe),ntris(nxyz),itris(nxyz+1),
     $   ltris(2*(2*nxyz-4)),lptr(2*nxyz-4,2),
     $   llpt(maxp,maxe),npatch
      common/patchg/  fluxv(maxp,maxe),gradm(maxp,nxyz,maxe)
c
c  np(i):                   total # patches of one sign for epoch i
c  npts(j,i):              # vertices (points) in patch j for epoch i
c  isign(j,i):             sign of patch j for epoch i
c  ipts(j,i):              starting index in lpts of the list of points 
c                        in patch j, for epoch i
c  lpts(ipts(j)+k-1,i):    kth vertex in patch j for epoch i
c  ntris(j):             # triangles in patch j
c  itris(j):             starting index in ltris of the list of 
c                        triangles in patch j
c  ltris(itris(j)+k-1):  kth triangle in patch j
c  lptr(j,k):            k=1:  positive patch (if any) associated with 
c                        triangle j.  k=2: negative patch (if any)
c  fluxv(j,i):  		 flux through patch j for epoch i
c  gradm(j,k,i): 		 gradient of flux through patch j due to change 
c			 in field at vertex k for epoch i
c
      common/tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
c
c  n:                    total # points on the core
c  ntria:                total # triangles on the core
c  p(k,j):               kth component of jth point, k=1,2,3; 
c                        Br(jth pt), k=4
c  iv(k,j):              kth vertex in jth triangle
c  near(k,j):            kth nearest neighbor of jth point
c
c  initialize indices.
      itris(1)=1
      np(i)=0
c 
	write(iout,*)' Entering patch',i
c	  do 1 k=1,maxp
c	  do 1 j=1,maxe
c1	  npts(k,i)=0
c
c  Are all vertices in a list?
      do 60 j=1, n
        do 10 k=1, np(i)
        do 10 m=1, npts(k,i)
          if ( lpts(ipts(k,i)+m-1,i) .eq. j ) go to 60
 10     continue
c  Vertex j is not on a list.  Make a new patch with j as its leading
c  point.
        np(i)=np(i)+1
	if(np(i).gt.30)then
	    write(iout,*)'too many patches, cant go on'
	    stop
	endif
        npts(np(i),i)=1
        isign(np(i),i)=signm(j)
	if(np(i).eq.1)then
	   ipts(np(i),i)=1
	else
           ipts(np(i),i)=ipts(np(i)-1,i)+npts(np(i)-1,i)
	endif
        lpts(ipts(np(i),i),i)=j
c
c  List all neighbors, neighbors of neighbors, etc. with same sign Br.
        do 50 jj=1,n
          if ( jj .gt. npts(np(i),i) ) go to 60
          now=lpts(ipts(np(i),i)+jj-1,i)
c  now is the current point.  Examine  now 's neighbors.
          do 40 kk=1, nvert
            next=near(kk, now)
c  next is candidate neighbor.  If next=0,  now  has no more neighbors.
            if ( next .eq. 0 ) go to 50
            if ( signm(next) .eq. signm(now) ) then
c  next should be on the list---is it already?
              do 30 mm=1, npts(np(i),i)
 30           if ( lpts(ipts(np(i),i)+mm-1,i) .eq. next ) go to 40
c  Add  next  to  lpts.
              npts(np(i),i)=npts(np(i),i)+1
              lpts(ipts(np(i),i)+npts(np(i),i)-1,i)=next
            endif
 40       continue
 50     continue
 60   continue
	do 61 j=1,np(i)
	 write(iout,*)' patch ',j,' number of points ',npts(j,i),
     $     ' sign ',isign(j,i)
c	 write(iout,*)' points are ',(lpts(ipts(j,i)+k-1,i),k=1,npts(j,i))
61	continue
c
c  Make a nonredundant list of all triangles including any point in 
c  patch j.
      do 90 j=1, np(i)
        ntris(j)=0
        do 80 k=1, npts(j,i)
        do 80 m=1, ntria
          do 70 l=1, 3
            if ( iv(l,m) .eq. lpts(ipts(j,i)+k-1,i) ) then
c  see if triangle m is already listed
              do 65 jj=1, ntris(j)
 65           if ( ltris(itris(j)+jj-1) .eq. m ) go to 80
c  it isn't; add it
              ntris(j)=ntris(j)+1
              ltris(itris(j)+ntris(j)-1)=m
              go to 80
            endif
 70       continue
 80     continue
c
c  find the starting index of the set of triangles for the next patch
        itris(j+1)=itris(j)+ntris(j)
 90   continue
c
	do 91 j=1,np(i)
c	 write(iout,*)' patch ',j,' number of triangles ',ntris(j)
c	 write(iout,*)' triangles are ',(ltris(itris(j)+k-1),k=1,ntris(j))
91	continue
c  For each patch, make a nonredundant list of points in the 
c  triangles that intersect it.
      do 110 j=1, np(i)
         npts(j,i)=0
         do 100 k=1, ntris(j)
         do 100 m=1, 3
c  Check if the next point is already on the list
            do 95 jj=1,npts(j,i)
 95         if ( lpts(ipts(j,i)+jj-1,i) .eq. 
     +        iv(m, ltris(itris(j)+k-1)) ) go to 100
c  it isn't; add it
            npts(j,i)=npts(j,i)+1
            lpts(ipts(j,i)+npts(j,i)-1,i) = iv(m,ltris(itris(j)+k-1))
 100     continue
c
c  Sort the list of points associated with each patch
	call isort(npts(j,i),lpts(ipts(j,i),i))
c
c  Update the point-patch index
         ipts(j+1,i)=ipts(j,i)+npts(j,i)
 110  continue
c
c  Make lptr, the cross-reference array of which patches a triangle
c  appears in.
      do 120 j=1,ntria
      lptr(j,1)=0
 120  lptr(j,2)=0
      do 130 j=1,np(i)
         if(isign(j,i).gt.0) then
c  Patch j is positive and should be in the first column of lptr.
            ndx=1
         else
c  Patch j should be in the second column of lptr.
            ndx=2
         endif 
         do 130 k=1,ntris(j)
            l=ltris(itris(j)+k-1)
            lptr(l,ndx)=j
 130  continue
      return
      end
c----------------------------------------------------------
      function signm(i)
	implicit double precision (a-h,o-z)
c
c  returns the sign of Br at the ith vertex.
c* Calls no other routines.
c
      parameter (nxyz=6500, nvert=10, nxyz6=6*nxyz, nxyz2=2*nxyz)
      common/tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
      signm=-1
      if ( p(4,i) .ge. 0.0 ) signm=1
      return
      end
c-----------------------------------------------------------
c_______________________________________________________________________
      subroutine shexp(ldeg,rr)
	implicit double precision (a-h,o-z)
c  For each triangle on the core,
c  performs the numerical cubature of br core model against partially
c  normalised spherical harmonics to get sh rep'n out to degree and
c  order ldeg at radius rr
c  input rr in km is normalized so that core radius =1
c  core field br at vertices p, is in 4th component of arrays
c
c  The quadrature rule generated by 'triq' for integrals over
c  plane triangle with vertices  (0,0,0) (1,0,0) (0,1,0) is assumed
c  to be available in common /rule/.
c
c  The list of observer coordinates is in /observ/.
c
c  The results are returned in common /calcul/
c
	parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,nvert=10)
	parameter (lmax=101)
	common /io/inp,iout,iprint
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
	common /rule/xyw(3,100),nquad
	common /shr/glm(lmax,lmax),hlm(lmax,lmax),plm(lmax,lmax),
     $        dlm(lmax,lmax)
      dimension rp(3),vc(4),aa(3),b(3),c(3),
     $          ua(2),ub(2),uc(2),rg(2),rpa(3)
	data pi/3.14159265/,rc/3486.d0/
	rrn=rr/rc
c
	xxx=0.
	if (ldeg.ge.lmax)then
	write(iout,*)'s.h.expansion truncated to ',lmax-1
	ldeg=lmax-1
	endif
c  initialize she to zero
	do 1900 j=1,ldeg
  	  do 1700 i=1,ldeg
	  glm(i,j)=0.
	  hlm(i,j)=0.
1700	  continue
1900	continue
c
	do 1500 j=1,ntria
c  
c  vertices of triangle are points iv1,iv2,iv3
	iv1=iv(1,j)
	iv2=iv(2,j)
	iv3=iv(3,j)
c
c  compute centroid vc of current triangle
c
	call centr(p(1,iv1),p(1,iv2),p(1,iv3),vc)
c	xxx=xxx + area(p(1,iv1),p(1,iv2),p(1,iv3))
c
c  rotate so vertices are in coord system with vc as pole, p(1,iv1) at
c  zero longitude
	call gpole(1,vc,p(1,iv1),p(1,iv1),aa)
	call gpole(1,vc,p(1,iv1),p(1,iv2),b)
	call gpole(1,vc,p(1,iv1),p(1,iv3),c)
c
c  find gnomonic projection of vertices and jacobian of transformation
c  for plane to unit triangle
c
	call gnomon(1,ua,aa)
	call gnomon(1,ub,b)
	call gnomon(1,uc,c)
	ajac = parea (ua,ub,uc)*2.
c
c  integrate over triangle
c
        do 2400 k=1, nquad
c
c  find points on gnomonic surface corresponding to integration points
c
	rg(1)= ua(1) + (ub(1)-ua(1))*xyw(1,k)  +(uc(1)-ua(1))*xyw(2,k)
	rg(2)= ua(2) + (ub(2)-ua(2))*xyw(1,k)  +(uc(2)-ua(2))*xyw(2,k)
c
c  find points in physical space corres. to rg
c
	call gnomon(-1,rg,rpa)
	call gpole(-1,vc,p(1,iv1),rp,rpa)
c
c  convert to colat, longitude
c
	call rthphi(1,rp(1),rp(2),rp(3),r, theta,phi)
	ct=cos(theta)
c
c  Evaluate Legendre polynomials at this site, and field
c
	call plmdlm(ct,ldeg,lmax,plm,dlm)
	gam1=bint(rp,rg,k,p(1,iv1),p(1,iv2),p(1,iv3))
c
c  Find Jacobian for gnomonic to spherical projection
c
	bjac= sqrt((1+rg(1)**2 + rg(2)**2)**-3)
c
c  evaluate partially normalized spherical harmonic functions at rp 
c
c  Weight integral according to xyw rule
c
      do 2500 il=0,ldeg
	  im=0
	  const=float(2*il+1)/(rrn**(il+2)*float(il+1)*4.0*pi)
	  glm(il+1,im+1)=glm(il+1,im+1) + gam1*const*
     $       plm(il+1,im+1)*cos(float(im)*phi)*ajac*bjac*xyw(3,k)
	  hlm(il+1,im+1)=hlm(il+1,im+1) + gam1*const*
     $       plm(il+1,im+1)*sin(float(im)*phi)*ajac*bjac*xyw(3,k)
	 do 2600 im=1,il
	  part=sqrt(2.*elf(il-im)/elf(il+im))*(-1)**im
	  glm(il+1,im+1)=glm(il+1,im+1) + gam1*const*
     $       plm(il+1,im+1)*cos(float(im)*phi)*ajac*bjac*xyw(3,k)*part
	  hlm(il+1,im+1)=hlm(il+1,im+1) + gam1*const*
     $       plm(il+1,im+1)*sin(float(im)*phi)*ajac*bjac*xyw(3,k)*part
2600	 continue	
2500 	continue
2400   continue
1500	continue
c	print *,'sum of s. triangle areas',xxx
      return
      end
c
c______________________________________________________________________
	function elf(l)
	implicit double precision (a-h,o-z)
c  Finds l!
	elf=1.0
	do 1 i=1,l
1	elf=elf*float(i)
	return
	end
c_______________________________________________________________________
      subroutine plmdlm(c, lmax, ldim, plm, dlm)
	implicit double precision (a-h,o-z)
c$$$$$ calls no other routines
c  finds associated legendre functions of all orders and degrees for
c  a fixed argument by recurrence.  the sign convention is that of
c  abramowitz & stegun, so that plm(l=1,m=1) = - sin theta.  precisely
c    plm = (-1)**m * (1-x**2)**(m/2) * (d/dx)**m pl(x),
c  where  pl(x)  is the legendre polynomial of degree  l.  special care
c  has been taken to avoid loss of precision or failure at or near  c=1.
c
c  c     is the (real) argument,  c .le. 1.0 in magnitude
c  lmax  is the maximum degree desired.
c  ldim  is the 1st dimension of both  plm and dlm  in the caller.  note
c        ldim .ge. lmax+1
c  plm()  array for results.  plm  of degree  l  and order  m  is held
c        in element  plm(l+1,m+1).  the elements  plm(l,l+1)  are used,
c        but other elements above the diagonal are not disturbed.
c  dlm() array for derivatives  (d/d theta) plm,  where  c=cos theta.
c
      dimension plm(ldim,lmax), dlm(ldim,lmax)
c
      plm(1,1)=1.0
      dlm(1,1)=0.0
      if (lmax.eq.0) return
c  zero out super-diagonal.
      do 1100 l=1, lmax
 1100 plm(l,l+1)=0.0
      s=sqrt(dmax1(0.D0, 1.D0 - c*c))
      twos=2.0*s
c  recur diagonally.
      do 1300 l=1, lmax
      elms=l*s
      mmax=lmax - l + 1
      do 1200 m=1, mmax
      lm=l + m
      plm(lm,m+1)=c*plm(lm-1,m+1) - elms*plm(lm-1,m)
 1200 elms=twos + elms
 1300 plm(l+1,1)=c*plm(l,1) + s*plm(l,2)/l
c  find derivatives.
      do 1500 l=1, lmax
      dlm(l+1,1)=plm(l+1,2)
      do 1500 m=1, l
      if (m .lt. l)
     $ dlm(l+1,m+1)=0.5*(plm(l+1,m+2) - (l+m)*(l-m+1)*plm(l+1,m))
 1500 if (m .eq. l) dlm(l+1,l+1)=-l*plm(l+1,l)
      return
      end                                                               plmdlm
c
c_______________________________________________________________________
	subroutine rthphi(iway,x,y,z,r,theta,phi)
	implicit double precision (a-h,o-z)
c
c  if iway=1 transforms x,y,z to r, theta, phi
c  if iway=-1 inverse transformation is performed
c  
	if (iway.eq.1)then
	  r=sqrt(x**2 + y**2 +z**2)
	  theta=acos(z/r)
	  phi=atan2(y,x)
	else if (iway.eq.-1)then
	  x=r*sin(theta)*cos(phi)
	  y=r*sin(theta)*sin(phi)
	  z=r*cos(theta)
	else
	  write(iout,*)'invalid subroutine call to rthphi'
	endif
	return
	end
c_______________________________________________________________________
	function bint(rp,rg,k,p1,p2,p3)
	implicit double precision (a-h,o-z)
c  interpolates b to a point rg on the gnomonic projection
	common /rule/xyw(3,100),nquad
	dimension rp(3),rg(2),gamma(3),p1(4),p2(4),p3(4)
	gamma(1)=1.-xyw(1,k)-xyw(2,k)
	gamma(2)=xyw(1,k)
	gamma(3)=xyw(2,k)
	bint=gamma(1)*p1(4) +gamma(2)*p2(4) +gamma(3)*p3(4)
	return
	end
c_______________________________________________________________________
	subroutine obs(kdim)
	implicit double precision (a-h,o-z)
c  if kdim=3 reads the observation points for the magnetic field into
c  array ro.
c  if kdim=6 reads magnetic field values at those points as well
c  ro(i,j,k) :  i=1,3 position of jth observation (r,theta,phi),kth
c               epoch
c  bro(i,j,k) :  i=1,3 magnetic field at jth observation point (br,
c  btheta,bphi in local coords), kth epoch
	parameter (maxob=7000,maxe=2)
	character *80 filnam(2),dform(2)
	common /io/inp,iout,iprint
	
	common /observ/ro(3,maxob,maxe),bro(maxob,maxe),sigma(maxob,
     $        maxe),iflag(maxob,maxe),mobs,nobs(maxe),nepoch,
     $	     lobs(maxe)
	dimension x(3),y(3)
	data drad/57.295780/
	mobs=0
	call getval('epochs',val,1,nfound)
	nepoch=nint(val)
	do 1 i=1,nepoch
1	nobs(i)=0
	if (kdim.eq.3) then
	call getchr('e1file',filnam(1),nfound)
	    open (unit=16, file =filnam(1))
	call getchr('e1format',dform(1),nfound)
	call getval('rmse',sigma,nepoch,nfound)
	    do 20 j=1,maxob+1
	      read(16,dform(1),err=101,end=101)(ro(i,j,1),i=1,3) 
	      mobs=mobs+1
20	    continue
101	write(iout,*) mobs, 'observation points read'
	else if (kdim.eq.6) then
	    call getchr('e1file',filnam(1),nfound)
	    call getchr('e1format',dform(1),nfound)
	if (nepoch.eq.2)then
	    call getchr('e2file',filnam(2),nfound)
	    call getchr('e2format',dform(2),nfound)
 	endif
	    do 102 k=1,nepoch
	    open (unit=16, file =filnam(k))
	if (dform(k).eq.'()')then
c  Assumed format is one component per line, observatory type
c  Ellipticity correction applied.
	    do 30 j=1,maxob+1
	         read(16,*,end=100,err=100)iflag(j,k),ro(2,j,k),
     $                 ro(3,j,k),ro(1,j,k),bro(j,k),sigma(j,k)
	         if(iflag(j,k).lt.0)goto 30
c  manipulate to form required here
                ro(3,j,k)=ro(3,j,k)/drad
		  h=ro(1,j,k)/1000.
	         ro(1,j,k)=6378.139
                ro(2,j,k)=ro(2,j,k)/drad
c  do ellipticity correction
	        call ellipt(ro(1,j,k),ro(2,j,k),ro(1,j,k),ro(2,j,k))
		 ro(1,j,k)=(ro(1,j,k) + h)/3485
	        nobs(k)=nobs(k)+1
30	    continue
	    write(iout, '(a)')'data file ', filnam
	else
c  Assumed format is satellite data, x,y,z, manipulated to above single 
c  component form
c  Read measurement uncertainty from input file
	    call getval('rmse',sigma(1,k),1,nfound)
	    do 40 j=1,maxob+1
	      read(16,dform(k),end=100,err=100)y(2),y(3),y(1),x(1),
     $                                     x(2),x(3)
	      do 39 m=1,3
	      if(abs(x(m)-99999.).gt.1.0e-05)then
	        nobs(k)=nobs(k)+1
	        ro(3,nobs(k),k)=y(3)/drad
	        ro(2,nobs(k),k)=(90.-y(2))/drad
	        ro(1,nobs(k),k)=(6371.2 + y(1))/3485.7
		 bro(nobs(k),k)=x(m)
		 iflag(nobs(k),k)=m+3
		 sigma(nobs(k),k)=sigma(1,k)
	      endif
39	      continue
40	    continue
	endif
100	write(iout, '(a)') 'format ',dform(k)
	write(iout,'(a)')'data file ', filnam(k)
	write(iout, *)
     $   nobs(k), ' observation points and magnetic field data read'
	mobs=mobs +nobs(k)
102	continue
	endif
	write(iout, *) 'Total observation points',mobs
	if (mobs.gt.maxob)then
	  write(iout, *)'maxob = ',maxob, ' recompile with larger array 
     $ dimensions'
	stop
	endif
	return 
	end
c
c_______________________________________________________________________
	subroutine ellipt(r,theta,gr,gtheta)
	implicit double precision (a-h,o-z)
c
c******calls no other routines
c  
c  Converts geodetic radius, gr, and colatitude, gtheta, into
c  geocentric, r, theta
c  All angles in radians.
c
	data pi/3.14159265/,a/6378.139/,b/6356.750/
c
	esq=1.0 - (b/a)**2
	coesq=(b/a)**2
	glam=pi/2. - gtheta
	theta=atan(1.0/(coesq*tan(glam)))
	if (theta.lt.0.)theta=theta + pi
	r=a*sqrt((1-esq*(2-esq)*(sin(glam))**2)/(1-esq*(sin(glam))**2))
	return
	end
c
c______________________________________________________________________
c  UNIT 2:  CORE DATA ROUTINES
c_______________________________________________________________________
      subroutine incore
	implicit double precision (a-h,o-z)
c
c$$$$ calls no other routines
c  The input routine for the downward continuation of geomagnetic
c  field to the surface of the core
c  Reads from the standard input until eof or continue
c  Saves the lines in the character array input for later
c  parsing by getchr and getval.
c  Prints a list of commands for this group of routines
c  Option exists to read all inputs from a specified file
c
      parameter (inmx=200)
      character *80 input(inmx),line,filnam
      common /dicta/ input
      common /ndict/ iecho,nin,istate(inmx)
      common /io/ inp,iout,iprint
	dimension val(20)
c
      nin=0
      write(iout,'(a)') ' ',
     $'                    =================',' ',
     $'Enter commands for core field inversion (? for help)'
c
      do 1500 l=nin+1,inmx
        read(*,'(80a)', end=2000) line
        if (line(1:4) .eq. 'cont') goto 2000
		if (line(1:4).eq. 'cfil') goto 1900
        if (line(1:1) .eq. '?') then
       write(iout, '(a/a/(2x,a))')
     $'Enter commands from the following list:',
     $'   (M) means mandatory information (O) optional',' ',
     $'?:        Remind me of the command list again',
     $'continue: Quit reading commands and begin calculations',
     $'cfile filnam: Read all commands from filnam',
     $'obsdim n: dimension of observation array, 3 for just points',
     $'       6 for field values too (M)',
     $'e1file file: name of file with evaluation points for forward',
     $ '       model or observations for inversion',
     $'e1format a: data format in e1file (M)',
     $'e2file file: name of file with evaluation points for forward',
     $ '       model or observations for inversion at second epoch (O)',
     $'e2format a: data format in e2file (O)',
     $'design n:   Construct design matrix (1), read it (0) or read',
     $'     upper triangular form after initial qr (-1) (M)',
     $'problem n:  Forward (1) or inverse (-1) or fflux (0)',
     $'     problem (M)',
     $'corefile file: file with core model for forward calculation (O)'
       write(iout, '(/(2x,a))')
     $'shrep ldeg r: degree and radius at whish to evaluate spherical',
     $'       harmonic representation (O)',
     $'patches 1: Find null flux curves and flux through patches (O)',
     $'epochs n:  Number of time snapshots (M)',
     $'rmse noise:  rms msifit level for obsfiles (M)',
     $'lambda rlam1, rlam2,...:  Lagrange multipliers (M)',
     $'tnew n: tesselation exists (0), needs creating (1,2), 2 saves a',
     $'       file for plotting',
     $'tessel file: file with tesselation points (M)',
     $'pointer file: file with pointers for triangles (M)',
     $'bodydim n: dimension of tesselation points (M)',
     $'zero n1,n2,...: make the listed points in the tesselation have',
     $'      zero radial field',
     $'outfil file: file for output of core model(s)(O)',
     $'fiter n: no of fr flux iterations allowed',
     $'invtype minbrsq (or mingrad or mindiff), type of',
     $'      regularization',
     $' '
	endif
        if (line(1:1) .ne. ' ') then
          nin=nin + 1
          input(nin)=line
          istate(nin)=0
        endif
 1500 continue
 1900  nin=nin+1
 	   input(nin)=line
	   istate(nin)=0
       call getchr('cfil',filnam,nfound)
	   write(iout,'(a)') 'Commands read from file ',filnam
       open(unit=10,file=filnam)
      do 1600 l=nin+1,inmx
        read(10,'(80a)', end=2000) line
        if (line(1:4) .eq. 'cont') goto 2000
        if (line(1:1) .eq. '?') then
       write(iout, '(a/a/(2x,a))')
     $'Enter commands from the following list:',
     $'   (M) means mandatory information (O) optional',' ',
     $'?:        Remind me of the command list again',
     $'continue: Quit reading commands and begin calculations',
     $'cfil filnam: Read all subsequent commands from file filnam',
     $'obsdim n: dimension of observation array, 3 for just points',
     $'6 for field values too (M)',
     $'e1file file: name of file with evaluation points for forward',
     $ 'model or observations for inversion',
     $'e1format a: data format in e1file (M)',
     $'e2file file: name of file with evaluation points for forward',
     $ 'model or observations for inversion at second epoch (O)',
     $'e2format a: data format in e2file (O)',
     $'design n:   Construct design matrix (1), read it (0) or read',
     $'upper triangule form after initial qr (-1) (M)',
     $'problem n:  Forward (1) or inverse (-1) problem (M)',
     $'corefile file: file with core model for forward calculation (O)',
     $'shrep ldeg r: degree and radius at whish to evaluate spherical',
     $'harmonic representation (O)'
       write(iout, '(/(2x,a))')
     $'patches 1: Find null flux curves and flux trhough patches (O)',
     $'epochs n:  Number of time snapshots (M)',
     $'rmse noise:  rms msifit level for obsfiles (M)',
     $'lambda rlam1, rlam2,...:  Lagrange multipliers (M)',
     $'tnew n: tesselation exists (0), needs creating (1,2), 2 saves a',
     $'file for plotting',
     $'tessel file: file with tesselation points (M)',
     $'pointer file: file with pointers for triangles (M)',
     $'bodydim n: dimension of tesselation points (M)',
     $' '
	endif
        if (line(1:1) .ne. ' ') then
          nin=nin + 1
          input(nin)=line
          istate(nin)=0
        endif
 1600 continue
c
c  Check that all mandatory parameters are available
c
2000	nflag=0
	call getchr('e1file',filnam,nfound)
	if(nfound.lt.0)nflag=1
	call getchr('tessel',filnam,nfound)
	if(nfound.lt.0)nflag=1
	call getchr('pointer',filnam,nfound)
	if(nfound.lt.0)nflag=1
	call getval('obsdim',val,1,nfound)
	if(nfound.lt.0)nflag=1
	call getval('bodydim',val,1,nfound)
	if(nfound.lt.0)nflag=1
	call getval('design',val,1,nfound)
	if(nfound.lt.0)nflag=1
	call getval('problem',val,1,nfound)
	if(nfound.lt.0)nflag=1
	call getval('epochs',val,1,nfound)
	if(nfound.lt.0)nflag=1
	call getval('rmse',val,1,nfound)
	if(nfound.eq.0)nflag=1
	call getval('lambda',val,1,nfound)
	if(nfound.eq.0)nflag=1
	call getchr('e1format',filnam,nfound)
	if(nfound.lt.0)nflag=1
c
c  If not all there then stop now
	if (nflag.eq.1) then
	  write(iout, *)'Insufficient information supplied.'
	  write(iout, '(a)')'Check list of mandatory parameters and your input 
     $for omissions or errors'
	  stop
	endif
 	return
      end
c
c
c=======================================================================
c==UNIT 6:  MISCELLANEOUS ROUTINES======================================
c=======================================================================
      subroutine getval(code, values, nwant, nfound)
	implicit double precision (a-h,o-z)
c$$$$ calls no other routines
c
c  Evaluates numbers in the string  char.  nwant  numbers are expected,
c  nfound  are actually found in the character string.  If an error is
c  discovered  nfound=-n,  where  n  is the number of numbers properly
c  decoded.  If there are no numbers after the codeword nfound=0.
      parameter (inmx=200)
      character *80 input(inmx),line
      common /dicta/ input
      common /ndict/ iecho,nin,istate(inmx)
      common /io/ inp,iout,iprint
      dimension values(*)
      character bore(2)*1, form*8, code*4, char*80
c
      data bore/'e',' '/
c
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
        if (code .eq. line(1:4) ) then
          if (iecho .ge. 1) write(iout,'(2a)')'>>> ',line
          istate(lin)=istate(lin) + 1
          nb=index(line, ' ')+1
          char=line(nb:80)
          kn=80 - nb + 1
          goto 1020
        endif
 1010 continue
c  code word not found
      nfound=-99
      return
c
 1020 continue
      k1=1
c  Up to  nwant  numbers are sought.
      do 1800 nval=1, nwant
        do 1100 k=k1, kn
          if (char(k:k) .ne. ' ') goto 1200
 1100   continue
        nfound=nval-1
        return
 1200   iex=1
        do 1300 l=k, kn
          if (char(l:l).eq. ',') goto 1500
          if (char(l:l).ne. ' ') goto 1300
          if (char(l-1:l-1).ne.bore(iex)) goto 1500
          iex=2
 1300   continue
 1500   m=l - k
c  Augment unpredictable error trapping of some compilers e.g. Sun
        if (index('Ee+-.0123456789', char(k:k)) .eq. 0) goto 2000
        k2=l+1
        write(form, '(2h(f, i3, 3h.0) )') m
        read (char(k:l-1), form, err=2000) values(nval)
 1800 k1=k2
      nval=1 - nwant
 1900 nfound=1 - nval
      return
c
 2000 write(iout, *)' ',
     $'***** Unreadable numbers in this input line:',line
      goto 1900
      end
c______________________________________________________________
      subroutine getchr(code, char, nfound)
	implicit double precision (a-h,o-z)
c$$$$ calls no other routines
      parameter (inmx=200)
      character *80 input(inmx),  line,  code*4
      character *(*)char
      common /dicta/ input
      common /ndict/ iecho,nin,istate(inmx)
      common /io/ inp,iout,iprint
ce
      do 1010 lin=nin, 1, -1
        line=input(lin)
        if (code .eq. line(1:4) ) then
          if (iecho .ge. 1) write(iout,'(2a)')'>>> ',line
          istate(lin)=istate(lin) + 1
          goto 1020
        endif
 1010 continue
c  code word not found
      nfound=-99
      return
c
 1020 continue
      nb=index(line, ' ')+1
      do 1200 k=nb, 80
        if (line(k:k) .ne. ' ') then
          char=line(k:80)
          nfound=1
          return
        endif
 1200 continue
c
c  Blank field after code word
      nfound=0
      return
      end
c______________________________________________________________
      subroutine ql(nm, n, a, d, e, ierr)
	implicit double precision (a-h,o-z)
c$$$$ calls tred2, tql2
c  using  eispack  routines tred2, tql2,  solves the symmetric
c  eigenvalue-eigenvector problem for a real matrix.
c  on input
c  nm   row dimension of the symmetric array  a  in the caller.
c  n    order of the array (<= nm)
c  a    the real symmetric array to be treated
c  e     a working array at least  n  long
c
c  on output
c  d    the array of eigenvalues ascending
c  a    the corresponding array of eigenvectors, with row
c       dimension  nm.  original  a  is overwritten.
c  ierr 0  if all's well.
c       j  if eigenvalue number j and above not found.
c
c
      dimension a(nm,*), d(*), e(*)
c
      call tred2(nm, n, a, d, e, a)
c
      call tql2(nm, n, d, e, a, ierr)
      return
      end
c______________________________________________________________________
      subroutine tql2(nm, n, d, e, z, ierr)
	implicit double precision (a-h,o-z)
c$$$$ calls no other routines
c
      integer i,j,k,l,m,n,ii,l1,nm,mml,ierr
      dimension d(n),e(n),z(nm,n)                                
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input:
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement;
c
c        n is the order of the matrix;
c
c        d contains the diagonal elements of the input matrix;
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary;
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output:
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1;
c
c        e has been destroyed;
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues;
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     ------------------------------------------------------------------
c
c     :::::::::: machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c                machep = 16.0d0**(-13) for long form arithmetic
c                on s360 ::::::::::
	real*4 machep
      data machep/1.0e-06/                                              *doub
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      b = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = machep * (abs(d(l)) + abs(e(l)))
         if (b .lt. h) b = h
c     :::::::::: look for small sub-diagonal element ::::::::::
         do 110 m = l, n
            if (abs(e(m)) .le. b) go to 120
c     :::::::::: e(n) is always zero, so there is no exit
c                through the bottom of the loop ::::::::::
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     :::::::::: form shift ::::::::::
         l1 = l + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = sqrt(p*p+1.0d0)
         d(l) = e(l) / (p + sign(r,p))
         h = g - d(l)
c
         do 140 i = l1, n
  140    d(i) = d(i) - h
c
         f = f + h
c     :::::::::: ql transformation ::::::::::
         p = d(m)
         c = 1.0d0
         s = 0.0d0
         mml = m - l
c     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
         do 200 ii = 1, mml
            i = m - ii
            g = c * e(i)
            h = c * p
            if (abs(p) .lt. abs(e(i))) go to 150
            c = e(i) / p
            r = sqrt(c*c+1.0d0)
            e(i+1) = s * p * r
            s = c / r
            c = 1.0d0 / r
            go to 160
  150       c = p / e(i)
            r = sqrt(c*c+1.0d0)
            e(i+1) = s * e(i) * r
            s = 1.0d0 / r
            c = c * s
  160       p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     :::::::::: form vector ::::::::::
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         e(l) = s * p
         d(l) = c * p
         if (abs(e(l)) .gt. b) go to 130
  220    d(l) = d(l) + f
  240 continue
c     :::::::::: order eigenvalues and eigenvectors ::::::::::
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     :::::::::: set error -- no convergence to an
c                eigenvalue after 30 iterations ::::::::::
 1000 ierr = l
 1001 return
      end
c______________________________________________________________________
      subroutine tred2(nm, n, a, d, e, z)
	implicit double precision (a-h,o-z)
c$$$$ calls no other routines
c
      integer i,j,k,l,n,ii,nm,jp1
      dimension a(nm,n),d(n),e(n),z(nm,n)                        
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input:
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement;
c
c        n is the order of the matrix;
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output:
c
c        d contains the diagonal elements of the tridiagonal matrix;
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero;
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction;
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 100 j = 1, i
            z(i,j) = a(i,j)
  100 continue
c
      if (n .eq. 1) go to 320
c     :::::::::: for i=n step -1 until 2 do -- ::::::::::
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     :::::::::: scale row (algol tol then not needed) ::::::::::
         do 120 k = 1, l
  120    scale = scale + abs(z(i,k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = z(i,l)
         go to 290
c
  140    do 150 k = 1, l
            z(i,k) = z(i,k) / scale
            h = h + z(i,k) * z(i,k)
  150    continue
c
         f = z(i,l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         z(i,l) = f - g
         f = 0.0d0
c
         do 240 j = 1, l
            z(j,i) = z(i,j) / h
            g = 0.0d0
c     :::::::::: form element of a*u ::::::::::
            do 180 k = 1, j
  180       g = g + z(j,k) * z(i,k)
c
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
  200       g = g + z(k,j) * z(i,k)
c     :::::::::: form element of p ::::::::::
  220       e(j) = g / h
            f = f + e(j) * z(i,j)
  240    continue
c
         hh = f / (h + h)
c     :::::::::: form reduced a ::::::::::
         do 260 j = 1, l
            f = z(i,j)
            g = e(j) - hh * f
            e(j) = g
c
            do 260 k = 1, j
               z(j,k) = z(j,k) - f * e(k) - g * z(i,k)
  260    continue
c
  290    d(i) = h
  300 continue
c
  320 d(1) = 0.0d0
      e(1) = 0.0d0
c     :::::::::: accumulation of transformation matrices ::::::::::
      do 500 i = 1, n
         l = i - 1
         if (d(i) .eq. 0.0d0) go to 380
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(i,k) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * z(k,i)
  360    continue
c
  380    d(i) = z(i,i)
         z(i,i) = 1.0d0
         if (l .lt. 1) go to 500
c
         do 400 j = 1, l
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  400    continue
c
  500 continue
c
      return
      end
c______________________________________________________________________
	subroutine eigv(f,ev,r,tr)
	implicit double precision (a-h,o-z)
c  Calculates the eigenvalues , ev(k), and eigenvectors, r(i,k), of a
c  real symmetric 3x3 matrix, with elements f1-f6.
c  tr=trace, f1=m11,f2=m22, f3=m33, f4=m12, f5=m13, f6=m23
	dimension f(6),ev(3),r(3,3),g(6),yn(3)
	data tol/2.d-14/
	tr=(f(1) + f(2) + f(3))/3.d0
c***put traceless part into g
	 do 1 i=1,3
	   g(i)=f(i)-tr
1	g(i+3)=f(i+3)
c***cope with zero elements
	do 3 i=1,6
 	  if(g(i).eq.0.d0)g(i)=tol
3	continue
	g62=g(6)*g(6)
	g52=g(5)*g(5)
	g42=g(4)*g(4)
	fm2=0.5d0*(g(1)*g(1) + g(2)*g(2) + g(3)*g(3)) + g42 +g52 + g62
	fm=sqrt(fm2)
	a0=g(1)*g62 + g(2)*g52 + g(3)*g42
     $                  - g(1)*g(2)*g(3) - 2.d0*g(4)*g(5)*g(6)
	a=a0/(fm*fm2)
	da=abs(a)
c***iterate for largest root
	x=1.d0
5	x2=x*x
	  eps=(da +x*(1.d0-x2))/(3.d0*x2-1.d0)
	  x=x+eps
	if(abs(eps).gt.tol)goto5
	if(a.gt.0.d0)x=-x
	ss=sqrt(dmax1(1.0 + 4.0*a/(x**3),tol))
c***construct all eigenvals
	ev(1)=x*fm 
	ev(2)=-0.5d0*ev(1)*(1.d0 +ss)
	ev(3)=-0.5d0*ev(1)*(1.d0 -ss)
c***do eignevectors
	do 50 k=1,3
	  d1=g(1)- ev(k)
	  d2=g(2)- ev(k)
	  d3=g(3)- ev(k)
	  y11=d2*d3-g62
	  y21=g(5)*g(6) -d3*g(4)
	  y31=g(4)*g(6) -d2*g(5)
	  y22=d1*d3-g52
	  y32=g(4)*g(5) -d1*g(6)
	  y33=d1*d2-g42
	  yn(1)=y11*y11 +y21*y21 + y31*y31
	  yn(2)=y21*y21 +y22*y22 + y32*y32
	  yn(3)=y31*y31 +y32*y32 + y33*y33
	  is=1
	  if(yn(2).gt.yn(1)) is=2
	  if(yn(3).gt.yn(is)) is=3
	  ynor=sqrt(yn(is))
	  goto(31,32,33),is
31	    r(1,k)=y11/ynor
	    r(2,k)=y21/ynor
	    r(3,k)=y31/ynor
	goto 50
32	    r(1,k)=y21/ynor
	    r(2,k)=y22/ynor
	    r(3,k)=y32/ynor
	goto 50
33	    r(1,k)=y31/ynor
	    r(2,k)=y32/ynor
	    r(3,k)=y33/ynor
50	continue
	ev(1)=ev(1) +tr
	ev(2)=ev(2) + tr
	ev(3) = ev(3) + tr
	return
	end
c______________________________________________________________________________________
	subroutine patceq
	implicit double precision (a-h,o-z)
c
c  Finds equivalent patches for different epochs
c  Assumes fluxpt has been called to give pointer to patches ranked
c  by magnitude of the flux through them
c  Starts with largest patches, checks sign and points in common, 
c  works down in size.
c  More than one patch may be equivalent to another if less than 60%
c  of points are in common.  If no equivalent patch exists then add a
c  null patch to equivalence
c  In difficult cases option exists to specify patch equivalents in
c  run file and run one iteration at a time
c  llpt is reordered so that llpt(i,1) is a patch equivalent to llpt(i,2)
c  
	parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,nvert=10)
	parameter (maxp=30,lmax=101,maxe=2)
	common /io/inp,iout,iprint
	common /tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
        common/patchy/np(maxe),npts(maxp,maxe),isign(maxp,maxe),
     $   ipts(maxp,maxe),
     $   lpts(2*nxyz,maxe),ntris(nxyz),itris(nxyz+1),
     $   ltris(2*(2*nxyz-4)),lptr(2*nxyz-4,2),
     $   llpt(maxp,maxe),npatch
        common/patchg/  fluxv(maxp,maxe),gradm(maxp,nxyz,maxe)
	dimension val(30),nllpt(maxp,maxe),
     $     nprag(maxp,maxe),nmrag(maxp,maxe),npr(2),nmr(2)
	data iff/0/
c
	nepoch=2
	call getval('equiv',val,30,nfound)
	if (nfound.gt.0.and.iff.eq.0)then
	     npatch =nfound/2
	      write(iout,*)npatch,' equivalences read'
	  do 1 i=1,nfound/2
	   llpt(i,1)=nint(val(2*i-1))
	   llpt(i,2)=nint(val(2*i))
1	   continue
	iff=1
	return
	else
	nnp=min(np(1),np(2))
	mnp=max(np(1),np(2))
	k=0
	do 3 j=1,2
	npr(j)=0
3	nmr(j)=0
c
c  Starting with biggest patches of positive and negative flux, 
c  sort their arrays of points by index. If less
c  than 60% are in common save these patches for reordering.  Otherwise keep these 
c  patches as equivalent. Goto next patch etc. 
c
	do 100 i=1,nnp/2
c  Start with negative patch
	  call isort(npts(llpt(i,1),1),lpts(ipts(llpt(i,1),1),1))
	  call isort(npts(llpt(i,2),2),lpts(ipts(llpt(i,2),2),2))
	  perc=equiv(lpts(ipts(llpt(i,1),1),1),npts(llpt(i,1),1),
     $      lpts(ipts(llpt(i,2),2),2),npts(llpt(i,2),2))
         if(perc.lt.60.0)then
c  Add to record of unmatched patches
	    do 997 j=1,2
	    if(fluxv(llpt(i,j),j).ge.0)then
	     npr(j)=npr(j)+1
	     nprag(npr(j),j)=llpt(i,j)
	    else if(fluxv(llpt(i,j),j).le.0)then
	     nmr(j)=nmr(j)+1
	     nmrag(nmr(j),j)=llpt(i,j)
	    endif
997	    continue
	  else
c  Otherwise keep the match gven by flux ordering
	     k=k+1
	     nllpt(k,1)=llpt(i,1)
	     nllpt(k,2)=llpt(i,2)
	     write(*,'(a,i2,a,i2,a,f6.2)')' Epoch 1, patch ',llpt(i,1),
     $           ' equivalenced to Epoch 2, patch ',llpt(i,2),
     $            ' percentage of points in common is ',perc
         endif
c  Goto positive patch
     	  call isort(npts(llpt(np(1)-i+1,1),1),
     $              lpts(ipts(llpt(np(1)-i+1,1),1),1))
	  call isort(npts(llpt(np(2)-i+1,2),2),
     $              lpts(ipts(llpt(np(2)-i+1,2),2),2))
	  perc=equiv(lpts(ipts(llpt(np(1)-i+1,1),1),1),
     $      npts(llpt(np(1)-i+1,1),1),lpts(ipts(llpt(np(2)-i+1,2),2),2),
     $      npts(llpt(np(2)-i+1,2),2))
         if(perc.lt.60.0)then
c  Add to record of unmatched patches
	    do 998 j=1,2
           if(fluxv(llpt(np(j)-i+1,j),j).ge.0.)then
	     npr(j)=npr(j)+1
	     nprag(npr(j),j)=llpt(np(j)-i+1,j)
	    else if(fluxv(llpt(np(j)-i+1,j),j).le.0.)then
	     nmr(j)=nmr(j)+1
	     nmrag(nmr(j),j)=llpt(np(j)-i+1,j)
	    endif
998	    continue
	  else
c  Otherwise keep the match gven by flux ordering
	    k=k+1
	    nllpt(k,1)=llpt(np(1)-i+1,1)
	    nllpt(k,2)=llpt(np(2)-i+1,2)
	    write(*,'(a,i2,a,i2,a,f6.2)')' Epoch 1, patch ',
     $           llpt(np(1)-i+1,1),
     $           ' equivalenced to Epoch 2, patch ',llpt(np(2)-i+1,2),
     $            ' percentage of points in common is ',perc
         endif
100	continue
	if(mod(nnp,2).eq.1)then
	   i=nnp/2+1
c  Do the last patch
	  call isort(npts(llpt(i,1),1),lpts(ipts(llpt(i,1),1),1))
	  call isort(npts(llpt(i,2),2),lpts(ipts(llpt(i,2),2),2))
	  perc=equiv(lpts(ipts(llpt(i,1),1),1),npts(llpt(i,1),1),
     $      lpts(ipts(llpt(i,2),2),2),npts(llpt(i,2),2))
         if(perc.lt.60.0)then
c  Add to record of unmatched patches
	    do 999 j=1,2
	    if(fluxv(llpt(i,j),j).ge.0)then
	     npr(j)=npr(j)+1
	     nprag(npr(j),j)=llpt(i,j)
	    else if(fluxv(llpt(i,j),j).le.0)then
	     nmr(j)=nmr(j)+1
	     nmrag(nmr(j),j)=llpt(i,j)
	    endif
999	    continue
	  else
c  Otherwise keep the match gven by flux ordering
	     k=k+1
	     nllpt(k,1)=llpt(i,1)
	     nllpt(k,2)=llpt(i,2)
	     write(*,'(a,i2,a,i2,a,f6.2)')' Epoch 1, patch ',llpt(i,1),
     $           ' equivalenced to Epoch 2, patch ',llpt(i,2),
     $            ' percentage of points in common is ',perc
         endif
	endif
c  If more patches at 1 epoch than other add these to appropriate rag end list
	if(np(1).gt.np(2))then
	  m1=1
	else
	  m1=2
	endif
	do 105 i=nnp/2+mod(nnp,2)+1,mnp-nnp/2
	   call isort(npts(llpt(i,m1),m1),lpts(ipts(llpt(i,m1),m1),m1))
c  Add to record of unmatched patches
	    if(fluxv(llpt(i,m1),m1).ge.0)then
	     npr(m1)=npr(m1)+1
	     nprag(npr(m1),m1)=llpt(i,m1)
	    else if(fluxv(llpt(i,m1),m1).le.0)then
	     nmr(m1)=nmr(m1)+1
	     nmrag(nmr(m1),m1)=llpt(i,m1)
	    endif
105    continue
	write(*,*)npr(1),npr(2),nmr(1),nmr(2)
	if (npr(1).gt.0.or.nmr(1).gt.0.or.npr(2).gt.0
     $  .or.nmr(2).gt.0) then
c  
c  Sort out rag ends by equivalencing all possible patches of same sign
c  to find those with most points in common
c  Find epoch with most patches and use this for reference
	write(*,'(a)')' Attempting to match +ve patches'
	if(npr(1).gt.npr(2))then
	   write(*,*)'Epoch 1 has more patches'
	   m1=1
	   m2=2
	else if(npr(2).gt.npr(1)) then
	   write(*,*)'Epoch 2 has more patches'
	   m1=2
	   m2=1
	else
	   write(*,*)'Epochs 1 and 2 have same no. of patches'
	   m1=1
	   m2=2
	endif
	do 200 i=1,npr(m1)
	  percmax=0.
	  k=k+1
	  nllpt(k,m1)=nprag(i,m1)
	  do 201 j=1,npr(m2)
	    perc=equiv(lpts(ipts(nprag(i,m1),m1),m1),npts(nprag(i,m1),m1),
     $            lpts(ipts(nprag(j,m2),m2),m2),npts(nprag(j,m2),m2))
           write(*,'(a,i2,a,i2,a,f6.2,a)')' Patches ',nprag(i,m1),
     $     ' and ', nprag(j,m2), 
     $       ' have ', perc,' percent of points in common'
c  Save match if better than before
	    if(perc.ge.percmax)then
	      percmax=perc
	      nllpt(k,m2)=nprag(j,m2)
	    endif
201	  continue
c  If no match equivalence to dummy patch no. 99
	     if(percmax.lt.1.) nllpt(k,m2)=99
	     write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.2)')
     $           ' Epoch ', m1,' patch ',nllpt(k,m1),
     $           ' equivalenced to Epoch ', m2,' patch ',nllpt(k,m2),
     $            ' percentage of points in common is ',percmax
200	continue
c  Check all npr(m2) patches have been included
	   do 199 j=1,npr(m2)
	   ii=0
	    do 198 i=1,k
	     if(nllpt(i,m2).eq.nprag(j,m2))ii=1
198	    continue
c  Add another patch if not there
	   if(ii.eq.0)then
	    k=k+1
	    nllpt(k,m2)=nprag(j,m2)
	     nllpt(k,m1)=99
	   endif
199	   continue
	write(*,'(a)')' Attempting to match -ve patches'
		if(nmr(1).gt.nmr(2))then
	   write(*,*)'Epoch 1 has more patches'
	   m1=1
	   m2=2
	else if(nmr(2).gt.nmr(1)) then
	   write(*,*)'Epoch 2 has more patches'
	   m1=2
	   m2=1
	else
	   write(*,*)'Epochs 1 and 2 have same no. of patches'
	   m1=1
	   m2=2
	endif
	do 203 i=1,nmr(m1)
	  percmax=0.
	  k=k+1
	  nllpt(k,m1)=nmrag(i,m1)
	  do 202 j=1,nmr(m2)
	    perc=equiv(lpts(ipts(nmrag(i,m1),m1),m1),npts(nmrag(i,m1),m1),
     $            lpts(ipts(nmrag(j,m2),m2),m2),npts(nmrag(j,m2),m2))
           write(*,'(a,i2,a,i2,a,f6.2,a)')' Patches ',nmrag(i,m1),
     $     ' and ', nmrag(j,m2), 
     $       ' have ', perc,' percent of points in common'
c  Save match if better than before
	    if(perc.ge.percmax)then
	      percmax=perc
	      nllpt(k,m2)=nmrag(j,m2)
	    endif
202	  continue
c  If no match equivalence to dummy patch no. 99
	     if(percmax.lt.1.0) nllpt(k,m2)=99
	     write(*,'(a,i2,a,i2,a,i2,a,i2,a,f6.2)')
     $           ' Epoch ', m1,' patch ',nllpt(k,m1),
     $           ' equivalenced to Epoch ', m2,' patch ',nllpt(k,m2),
     $            ' percentage of points in common is ',percmax
203	continue
c  Check all nmr(m2) patches have been included
	   do 196 j=1,nmr(m2)
	   ii=0
	    do 197 i=1,k
	     if(nllpt(i,m2).eq.nmrag(j,m2))ii=1
197	    continue
c  Add another patch if not there
	   if(ii.eq.0)then
	    k=k+1
	    nllpt(k,m2)=nmrag(j,m2)
	     nllpt(k,m1)=99
	   endif
196	   continue
	endif
	npatch=k
	write(*,*)'npatch = ',npatch
c  copy nllpt into llpt
	do 101 i=1,npatch
	do 101 j=1,nepoch
101	llpt(i,j)=nllpt(i,j)
	endif
	return
	end
c______________________________________________________________________________________
	function equiv(list1,n1,list2,n2)
	implicit double precision (a-h,o-z)
	dimension list1(n1),list2(n2)
c  Counts how many points list1 and list2 have in common and returns the
c  value 200*npts/(n1+n2).  Assumes list1 and list2 are sorted in ascending order
	count=0.
	i1=1
	i2=1
1	if(list1(i1).eq.list2(i2))then
	    count=count+1.
	    i1=i1+1
	    i2=i2+1
	    if(i1.gt.n1.or.i2.gt.n2)goto10
	else if(list1(i1).gt.list2(i2))then
	    i2=i2+1
	    if(i2.gt.n2)goto 10
	else
	    i1=i1+1
	    if(i1.gt.n1)goto 10
	endif
	goto 1
10	equiv=200.*count/float(n1+n2)
	return
	end
c______________________________________________________________________________________
	
	
c--------------------------------------------------------------
      subroutine pmaker(kk,i)
	implicit double precision (a-h,o-z)
c
c  Makes a patch on the core out of points in patch kk at epoch i
c  Computes the gradient matrix for this patch
c  Note that this patch is not distinct in sign from
c  its surroundings, and may even be mixed in sign
c
c  Assumes patch  and flux have already been called to define the 
c  variables in /patchy/ for epoch 2.
c
c
c  Outline:
c
c  Initialize the flux vector and flux gradient matrix to zero
c  Loop in triangles
c    for each triangle, compute rotation and gnomonic projection
c    loop in integration sample points
c       compute preimage of the sample point on the core
c       compute the jacobian for the gnomonic at that sample point
c       compute field value at that point
c       add the value, weighted by the int. wt. and the jacobians,
c           to the positive patch or negative patch, according
c           to its sign.
c       add interp. wt*int. wt.* jacobian to the appropriate elements
c           of the grad. matrix (according to signs of br and vertex)
c    next integration sample point
c  next triangle.
c
      parameter(nxyz=6500,nxyz2=2*nxyz,nxyz6=6*nxyz,maxp=30,
     $   nvert=10,maxe=2)
      dimension rp(3),g(3),vc(4),aa(3),b(3),c(3),
     +   ua(2),ub(2),uc(2),rg(2),rpa(3)
      common/patchy/np(maxe),npts(maxp,maxe),isign(maxp,maxe),
     $  ipts(maxp,maxe),lpts(2*nxyz,maxe),ntris(nxyz),
     $  itris(nxyz+1),
     $  ltris(2*(2*nxyz-4)),lptr(2*nxyz-4,2),
     $   llpt(maxp,maxe),npatch
      common/patchg/ fluxv(maxp,maxe),gradm(maxp,nxyz,maxe)
c
c  see subroutine patch for definitions of the variables in /patchy/
c
      common/tessel/n,ntria,p(4,nxyz),iv(3,nxyz2),near(nvert,nxyz)
      common /io/inp,iout,iprint
      common /rule/xyw(3,100),nquad
      common /circle/v(4,nxyz2),istack(nxyz2)
c
c  initialize flux vector and gradient matrix.
c******
      write(iout,*)' arrived in pmaker '
	fluxv(maxp,1)=0.
      do 10 j=1,n
 10   gradm(maxp,j,i)=0.
c
c  Loop in triangles associated with the patch
      do 100 m=1,ntris(kk)
	j=ltris(itris(kk) +m-1)
c  
c
c  compute centroid vc of current triangle
         call centr(p(1,iv(1,j)),p(1,iv(2,j)),p(1,iv(3,j)),vc)
c
c  rotate so vertices are in coord system with vc as pole, p(1,iv1) at
c  zero longitude
         call gpole(1,vc,p(1,iv(1,j)),p(1,iv(1,j)),aa)
         call gpole(1,vc,p(1,iv(1,j)),p(1,iv(2,j)),b)
         call gpole(1,vc,p(1,iv(1,j)),p(1,iv(3,j)),c)
c
c  find gnomonic projection of vertices and jacobian of transformation
         call gnomon(1,ua,aa)
         call gnomon(1,ub,b)
         call gnomon(1,uc,c)
         ajac = parea (ua,ub,uc)*2.
c*******
c      write(iout,*)' entering integration loop '
c
c  Loop in integration sample points.
         do 50 k=1, nquad
c
c  find points on gnomonic surface corresponding to integration points
            rg(1)= ua(1) + (ub(1)-ua(1))*xyw(1,k)  +
     $        (uc(1)-ua(1))*xyw(2,k)
            rg(2)= ua(2) + (ub(2)-ua(2))*xyw(1,k)  +
     $        (uc(2)-ua(2))*xyw(2,k)
c
c  Find Jacobian for gnomonic to spherical projection
            bjac= sqrt((1+rg(1)**2 + rg(2)**2)**-3)
c
c  Find interpolating vector g(.)---linear interpolation in gnomonic
c  plane.
            g(1)=1.-xyw(1,k)-xyw(2,k)
            g(2)=xyw(1,k)
            g(3)=xyw(2,k)
c
c  find the point rp on the core corresponding to rg.
            call gnomon(-1,rg,rpa)
            call gpole(-1,vc,p(1,iv(1,j)),rp,rpa)
c
c  Find field value and its sign.
            br=0.
            do 30 l=1,3
 30         br=br+g(l)*p(4,iv(l,j))
c
c  Update the flux vector fluxv and the gradient matrix gradm.
            fluxv(maxp,i)=fluxv(maxp,i)+
     $        ajac*bjac*xyw(3,k)*br
            do 40 l=1,3
              gradm(maxp,iv(l,j),i)=gradm(maxp,
     $         iv(l,j),i)+g(l)*ajac*bjac*xyw(3,k)
 40         continue
c
c  Next integration sample point.
 50      continue
c
c  Next triangle.
100   continue
      return
      end
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
c	print *, i1,i2,i3
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
