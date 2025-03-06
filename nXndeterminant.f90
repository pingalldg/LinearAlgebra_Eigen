!###### defining a matrix and  finding its determinant ############
	  program main
	  implicit none
	  real,allocatable,dimension(:,:)::arr
	  integer:: m,n,j,i,k,l,ind
          real::coeff,det

	  write(*,*)"write down the dimension of the matrix (should be an integer)"
	  read(*,*) m
	  allocate (arr(m,m))
	  call arrdef(arr,m)
          j=1
          ind=0
          do while(j.lt.m)
	  call findpivot(arr,m,j,i)
          call findpivot(arr,m,j+1,k)
!           write(*,*)i,k
           if (i.gt.k) then
           call interchange(arr,m,j,j+1)
           j=0
           ind=ind+1
           end if

           if(i.eq.k)then
           if(i.lt.m)then
           coeff=-(arr(j,i+1)/arr(j+1,i+1))
           call replacement(arr,m,j,j+1,coeff)
           write(*,*)"j,j+1 replacement",j,j+1,coeff,i
           j=0
           end if
          end if
           j=j+1
           end do
           call diagmult(arr,m,det)
           do n=1,m
           write(*,*)(arr(n,l),l=1,m)
           end do
           write(*,*)"determinant =",((-1)**ind)*det
          deallocate(arr)
	  end program main  

!#####################################################  
	!#######this subroutine will take the 

	subroutine arrdef(arr1,m)
	real::x
	integer::m
	real,dimension(m,m)::arr1
	write(*,*)"array declaration!!give a integer matrix"
	do j=1,m
	do k=1,m
	read(*,*)x
	arr1(j,k)=x
	end do 
	end do
        write(*,*)"array is defined"
	end subroutine arrdef 

!#####################################################         
 
       subroutine findpivot(arr1,m,j,i)
       integer::m,j,i,flag
       real,dimension(m,m)::arr1
       real::x
       flag=0
       i=0
       do while(flag.eq.0)
       i=i+1
       x=arr1(j,i)
       if((abs(x).gt.0.000000).or.(i.gt.m))then
       flag=1
       i=i-1
       else
       arr1(j,i)=0
       end if
       end do
!       write(*,*)"pivot position",i,"of ",j," th row"
       end subroutine findpivot
!#####################################################  

       subroutine interchange(arr1,m,l,k)
       integer::l,k,m,j
       real,dimension(m,m)::arr1
       real::x,y
       do j=1,m
       x=arr1(l,j)
       y=arr1(k,j)
       arr1(k,j)=x
       arr1(l,j)=y
       end do
!       write(*,*)"interchange between",l,"th and  ",k,"th row"
       end subroutine interchange
!#####################################################  

       subroutine multiply(arr1,m,l,coeff)
       integer::l,m,j
       real,dimension(m,m)::arr1
       real::coeff
       do j=1,m
       arr1(l,j)=coeff*arr1(l,j)
       end do
!      write(*,*)"multiply constant to the ",l,"th row"
       end subroutine multiply
!#####################################################  

       subroutine replacement(arr1,m,l,k,coeff)
       integer::l,k,m,j
       real,dimension(m,m)::arr1
       real::coeff
       do j=1,m
       arr1(l,j)=arr1(l,j)+coeff*arr1(k,j)
       end do
!       write(*,*)l,"th row +   ",coeff,"x  ",k,"th row"
       end subroutine replacement

!#####################################################  
 
       subroutine diagmult(arr1,m,d)
       integer::m
       real,dimension(m,m)::arr1
       real::d
       d=1.
       do i=1,m
       d=d*arr1(i,i)
       end do
!      write(*,*)d
       end subroutine diagmult
     
       
