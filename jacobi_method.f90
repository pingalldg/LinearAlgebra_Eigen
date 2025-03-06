program main
!this program evaluates eigen values using jacobi techniques.
real::x,beta
integer::i,j,flag,n1,n2
real,allocatable,dimension(:,:)::a,a1,q
write(*,*)"Square matrix dimension; >= 2"
read(*,*)n
allocate(a(n,n))
allocate(a1(n,n))
allocate(q(n,n))
beta=0.
write(*,*)"Now define your aray elements,give always symmetric matrix!"
do i=1,n
do j=1,n
read(*,*)x
a(i,j)=x
a1(i,j)=0.
q(i,j)=0.
end do
end do
!!!!printing the matrix!!!!
write(*,*)"you have given this array"
do i=1,n
write(*,*)(a(i,j),j=1,n)
end do
!!!!!!!!main body!!!!!!!!!!!!
flag=0
do while(flag.eq.0)
 call findgreat(a,n,n1,n2)

if(abs(a(n1,n2)).lt. 0.0001)then
exit
end if
beta=(a(n1,n1)-a(n2,n2))/(2*a(n1,n2))
!write(*,*)beta
!stop
call createmat(q,n,n1,n2,beta)
call matmult(a,q,a1,n)
call matmult(transpose(q),a1,a,n)
flag2=1
do i=1,n-1
do j=i+1,n
if(abs(a(i,j)).lt.0.0001)then
flag2=flag2+1
end if
end do
end do
iter=int((n*n-n)/2.)
if(flag2.eq.iter)then
 exit
end if
end do
!!!!printing final matrix!!!!
write(*,*)"**************"
do i=1,n
write(*,*)(a(i,j),j=1,n)
end do
!write(*,*)"**************"
!do i=1,n
!write(*,*)(a1(i,j),j=1,n)
!end do
!write(*,*)"**************"
!do i=1,n
!write(*,*)(q(i,j),j=1,n)
!end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program main
!!!!!!!!!!!!!!!!!!
subroutine matmult(b,c,d,n)
real,dimension(n,n)::b,c,d
integer::i,j,k
do i=1,n
do j=1,n
d(i,j)=0.
do k=1,n
d(i,j)=d(i,j)+b(i,k)*c(k,j)
end do

if(abs(d(i,j)).lt.0.0001)then
d(i,j)=0.
end if
end do
end do

end subroutine matmult

subroutine findgreat(a,n,n1,n2)
integer::n,i,j,n1,n2
real,dimension(n,n)::a
real::x
n1=1
n2=2
x=abs(a(n1,n2))
do i=1,n-1
do j=i+1,n
if(abs(a(i,j)).gt.x)then
x=abs(a(i,j))
n1=i
n2=j
end if
end do
end do

end subroutine findgreat


subroutine createmat(q,n,n1,n2,beta)
real::x,beta,c,s
integer::i,j,n1,n2,n
real,dimension(n,n)::q

 s=sqrt(1./2.-beta/(2.*sqrt(1+beta*beta)))
 c=sqrt(1./2.+beta/(2.*sqrt(1+beta*beta)))
! write(*,*)"******",c,s
do i=1,n
do j=1,n
q(i,j)=0.
if(i.eq.j)then
q(i,j)=1.
else
q(i,j)=0.
end if
end do
end do
q(n1,n1)=c
q(n2,n2)=c
q(n1,n2)=-s
q(n2,n1)=s
end subroutine createmat
