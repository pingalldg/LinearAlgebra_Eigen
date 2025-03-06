program main
!this program evaluates eigen values using qr decomposition techniques.
integer::i,j,flag,flag2,iter
real,allocatable,dimension(:,:)::a,a1,q,r,test
write(*,*)"Square matrix dimension"
read(*,*)n
allocate(a(n,n))
allocate(a1(n,n))
allocate(q(n,n))
allocate(r(n,n))
allocate(test(n,n))

write(*,*)"Now define your aray elements,give always symmetric matrix!"
do i=1,n
do j=1,n
read(*,*)x
a(i,j)=x
end do
end do
!!!!printing the matrix!!!!
write(*,*)"you have given this array"
do i=1,n
write(*,*)(a(i,j),j=1,n)
end do
!!!!!!!!!main body!!!!!!!!!!!!!!!
flag=0
do while (flag.eq.0)
!do iter=1,100
!!!!initialiation!!!
do i=1,n
do j=1,n
a1(i,j)=0.
q(i,j)=0.
r(i,j)=0.
test(i,j)=0.0
end do 
end do
!!!!!!!!!qr decomposition!!!!!!!!!!!!!!!

do j=1,n
call gramschmitt(a,a1,j,n)
end do

do j=1,n
!pg=0.
do i=1,n
q(i,j)=a1(i,j)
!pg=pg+q(i,j)*q(i,j)
end do 
!write(*,*)"write norm",pg
end do
call matmult(transpose(q),a,r,n)
!!!!!!!!!!!! new rq form!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call matmult(r,q,test,n)
!call matmult(r,q,a,n)
flag2=0
do i=1,n
if(abs(a(i,i)-test(i,i)).lt.0.001)then
flag2=1
else
flag2=0
exit
end if
end do

if(flag2.eq.0)then
do l=1,n
do m=1,n
a(l,m)=test(l,m)
end do 
end do
else
write(*,*)(test(i,i),i=1,n)
flag=1
end if
end do
write(*,*)"***********"
!do j=1,n
!write(*,*)(a(i,j),i=1,n)
!end do
!!!!!!!!!!!!!!!!!!!!debug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*)"*******************"
!do i=1,n
!write(*,*)(a(i,j),j=1,n)
!end do
!write(*,*)"*******************"
!do i=1,n
!write(*,*)(q(i,j),j=1,n)
!end do
!write(*,*)"*******************"
!do i=1,n
!write(*,*)(r(i,j),j=1,n)
!end do
!write(*,*)"*******************"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate(a)
deallocate(a1)
deallocate(q)
deallocate(r)
end

subroutine gramschmitt(ar,ar1,j,n)
real,dimension(n,n)::ar,ar1
real,dimension(n)::br
integer::n,j,k,m
real::s,coeff
!expression
do i=1,n
br(i)=ar(i,j)
k=1
do while((j-k).ge.1)
coeff=0.
do m=1,n
coeff=coeff+ar(m,j)*ar1(m,j-k)
end do
br(i)=br(i)-coeff*ar1(i,j-k)
k=k+1
end do
end do
s=0.
!norm
do i=1,n
s=s+br(i)*br(i)
end do
!new vector a1
do i=1,n
br(i)=br(i)/sqrt(s)
ar1(i,j)=br(i)
if(abs(ar1(i,j)).lt.0.0001)then
ar1(i,j)=0.
end if
end do 
end subroutine gramschmitt


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


