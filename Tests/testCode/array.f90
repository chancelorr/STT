! test arrays
program arrays

integer :: p(3, 6), q(3, 6)

DATA ((p(i,j),i=1,3),j=1,6)/0 , 0 , 1 , 0 , 0 , -1 , -1 , -1 , 0 ,&
			& 1 , -1 , 0 , 1 , 1 , 0 , -1 , 1 , 0/

open (unit=10, file='arr')
do j=1, 6
    write (10, *) (p(i, j), i=1, 3)
end do

print *, 'Values written to "arr"'

! read (10, *) q

! print *, q

end program arrays
