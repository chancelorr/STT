! test program for testing what fortran spits out

program tester
implicit none

integer :: numBig
real :: num, n

print *, 'enter num, n'
read *, num, n
print *, 'num: ', num

! numBig is an integer, this is the truncation step
numBig = num * (10.0 ** n)
num = numBig/(10.0 ** n)

print *, 'numBig: ', numBig
print *, 'num: ', num

end
