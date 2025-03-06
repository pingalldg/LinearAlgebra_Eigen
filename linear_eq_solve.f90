!###### solving linear equation ############
	  program main
	  implicit none
	  real,allocatable,dimension(:,:)::arr
	  integer:: m,n,j,i,k,l,o
          real::coeff,det

	  write(*,*)"write down the dimension of the matrix (should be an integer)"
	  read(*,*) m,o
	  allocate (arr(m,o))
	  call arrdef(arr,m,o)
          j=1

           do while(j.lt.m)
	   call findpivot(arr,m,o,j,i)
           call findpivot(arr,m,o,j+1,k)
           coeff=1./arr(j,i+1)
           call multiply(arr,m,o,j,coeff)
           coeff=1./arr(j+1,k+1)
           call multiply(arr,m,o,j+1,coeff)
!           write(*,*)i,k
           if (i.gt.k) then
            call interchange(arr,m,o,j,j+1)
            coeff=-(arr(j,i+1)/arr(j+1,i+1))
            call replacement(arr,m,o,j,j+1,coeff)
            j=0
           end if
           if(i.eq.k)then
           if(i.lt.o)then
           coeff=-1
           call replacement(arr,m,o,j,j+1,coeff)
!           write(*,*)"j,j+1 replacement",j,j+1,coeff,i
           j=0
           end if
           end if
           j=j+1
           end do
           
           do n=1,m
           write(*,*)"the reduced row echolon form",(arr(n,l),l=1,o)
           end do
           deallocate(arr)
	   end program main  

!#####################################################!

	subroutine arrdef(arr1,m,o)
	real::x
	integer::m,o
	real,dimension(m,o)::arr1
	write(*,*)"array declaration!!give a integer matrix"
	do j=1,m
	do k=1,o
	read(*,*)x
	arr1(j,k)=x
	end do 
	end do
        write(*,*)"array is defined"
	end subroutine arrdef 

!#####################################################         
 
       subroutine findpivot(arr1,m,o,j,i)
       integer::m,j,i,flag,o
       real,dimension(m,o)::arr1
       real::x
       flag=0
       i=0
       do while(flag.eq.0)
       i=i+1
       x=arr1(j,i)
       if((abs(x).gt.0.000000).or.(i.gt.o))then
       flag=1
       i=i-1
       else
       arr1(j,i)=0
       end if
       end do
!       write(*,*)"pivot position",i,"of ",j," th row"
       end subroutine findpivot

!#####################################################  

       subroutine interchange(arr1,m,o,l,k)
       integer::l,k,m,j,o
       real,dimension(m,o)::arr1
       real::x,y
       do j=1,o
       x=arr1(l,j)
       y=arr1(k,j)
       arr1(k,j)=x
       arr1(l,j)=y
       end do
!       write(*,*)"interchange between",l,"th and  ",k,"th row"
       end subroutine interchange

!#####################################################  

       subroutine multiply(arr1,m,o,l,coeff)
       integer::l,m,j,o
       real,dimension(m,o)::arr1
       real::coeff
       do j=1,o
       arr1(l,j)=coeff*arr1(l,j)
       end do
!      write(*,*)"multiply constant to the ",l,"th row"
       end subroutine multiply

!#####################################################  

       subroutine replacement(arr1,m,o,l,k,coeff)
       integer::l,k,m,j,o
       real,dimension(m,o)::arr1
       real::coeff
       do j=1,o
       arr1(l,j)=arr1(l,j)+coeff*arr1(k,j)
       end do
!       write(*,*)l,"th row +   ",coeff,"x  ",k,"th row"
       end subroutine replacement

!#####################################################  
     
       
