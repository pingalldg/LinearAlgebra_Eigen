program main
!this program evaluates largest eigen value of a matrix using power iteration method and its corresponding eigen value.....Can you find any similarity  with this and the variational princile in quantum mechanics where you determine the ground state energy.
integer::i,j,flag,iter
real::x,s,s1,mu,mu1
real,allocatable,dimension(:,:)::a
real,allocatable,dimension(:)::b,b1
write(*,*)"Square matrix dimension"
read(*,*)n
allocate(a(n,n))
allocate(b(n))
allocate(b1(n))
write(*,*)"Now define your aray elements"
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
!!initialization!!
do i=1,n
b1(i)=0.
b(i)=0.
end do
!!!start with  n dimensional test column vector-----b!!!
s=0.0
do i=1,n
! call change_seed()
 call init_random_seed()
 call random_number(x)
b(i)=x
s=s+x*x
end do
do i=1,n
b(i)=b(i)/sqrt(s)
end do
write(*,*)"you have made your starting test vector"
!!!Now we will start our iteration process!!!
flag=0
iter=0
mu1=200.
do while(flag.eq.0)
iter=iter+1
s=0.
s1=0.
mu=0.
!write(*,*)mu1,iter
do i=1,n
do j=1,n
b1(i)=b1(i)+a(i,j)*b(j)
end do
s1=s1+b1(i)*b1(i)
s=s+b1(i)*b(i)
mu=mu+b(i)*b(i)
end do
mu=s/mu
if(abs(mu-mu1).gt.0.0001)then
do i=1,n
b(i)=b1(i)/sqrt(s1)
b1(i)=0.
end do
mu1=mu
else
write(*,*)"your largest eigen value",mu1,mu
flag=1
end if
end do
end program main

 subroutine change_seed()
 integer :: i, n, clock
 integer, dimension(:), allocatable :: seed
 call random_seed(size = n)
 allocate(seed(n))
 call system_clock(count=clock)
 seed = clock + 37 * (/ (i - 1, i = 1, n) /)
 call random_seed(put = seed)
 deallocate(seed)
 end subroutine change_seed

  subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t
  call random_seed(size = n)
  allocate(seed(n))
  open(newunit=un, file='/dev/urandom', access='stream', status='old', action='read', form='unformatted', iostat=istat)
  if (istat == 0) then
    read(un) seed
    close(un)
  else
    call system_clock(t)
    if (t == 0) then
      call date_and_time(values = dt)
      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 + dt(3) * 24_int64 * 60 * 60   &
      * 1000 + dt(5) * 60 * 60 * 1000 + dt(6) * 60 * 1000 + dt(7) * 1000 + dt(8)
    end if
    pid = getpid()
    t = ieor(t, int(pid, kind(t)))
    do i = 1, n
      seed(i) = lcg(t)
    end do
  end if
  call random_seed(put = seed)
  deallocate(seed)
contains
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
      s = 104729
    else
      s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed

