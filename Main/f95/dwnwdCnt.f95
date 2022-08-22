program dwnwdCont
    implicit none

    !   Reads a snapshot spherical harmonic model (at earth's surface)
    !   Downward continues the gauss coefficients to a new specified radius
    !   Scaling is performed as coefficients are read in
    !
    !   by Chancelor Roberts, 08/22
    !
    !   input: infile, outfile, r, flag
    !   flag: 0 for potential, 1 for Br
    !
    !   Please insert full filepaths from calling directory
    !
    !   infile: sh model to be read
    !   outfile: sh model to be output (for mapping in gmt)
    !   r: new radius (as ratio to surface radius)
    !   a: surface radius = 6371e3 m
    !   g, h: gauss coefficients
    !   pref: prefactor to coefficients determined by flag
    !       0: pref = a (Psi)
    !       1: pref = a*(l+1)/r (Br)

    integer :: ll, mm, l, m, lmax
    integer :: nm, nspl, flag
    real*8 :: g, h, a, r, t, pref
    character :: infile*30, outfile*30

    data a/6371000/

    read (*, '(a)') infile
    write (*, *) infile
    read (*, '(a)') outfile
    write (*, *) outfile
    read (*, '(f5.3)') r
    read (*, *) flag
    if (flag.eq.0) print *, 'building Psi'
    if (flag.eq.1) print *, 'building Br'


    print *, 'opening ', infile
    open (unit=10, file=infile)
    open (unit=20, file=outfile)

    read (10, *) lmax, nm, nspl
    if (nspl.gt.1) print *, 'Multiple epochs provided, only reading first available model'
    write (20, '(i2, i2, i6, a)', advance='no') lmax , nm, 1, NEW_LINE('')
    read (10, *) t
    write (20, '(S, f19.13)') t

    print *, 'Downward continuing ', infile, ' to r=', r

    pref=a
!   read monopole
    do ll=0, lmax
        if (flag.eq.1) pref=(ll+1)/r
        do mm=0, ll
            read (10, *) l, m, g, h
!            print *, 'l, ', l, 'll, ', ll
            write (20,'(i2, i10, S, f30.13, f30.13, a)', advance='no') &
                l, m, pref*((1/r)**(l+1))*g, pref*((1/r)**(l+1))*h, NEW_LINE('')
        ENDDO
    enddo

    print *, 'New model written to ', outfile

    close (10)
    close (20)

end program dwnwdCont
