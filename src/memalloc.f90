! --------------------------------------------------------------------------
! -                                                                        -
! -  MEMALLOCATE : A storage manage package for finite element code        -
! -                                                                        -
! -     Xiong Zhang, (2013)                                                -
! -     Computational Dynamics Group, School of Aerospace                  -
! -     Tsinghua Univerity                                                 -
! -                                                                        -
! -  List of subroutine                                                    -
! -                                                                        -
! -     memalloca - allocate an array in the shared storage                -
! -     memfree   - deallocate the specified array                         -
! -     memfreefrom   - deallocate all arrays from the specified array     -
! -     memfreefromto - deallocate all arrays between the specified arrays -
! -     memprint      - print the contents of the specified array          -
! -     memprintptr   - print a subset of the storage in given format      -
! -     meminfo   - list all allocated arrays                              -
! -                                                                        -
! --------------------------------------------------------------------------


module memAllocate

   integer, parameter :: MTOT = 10000  ! Speed storage available for execution
   integer, parameter :: ITWO = 2      ! Double precision indicator
                                        !   1 - Single precision arithmetic
                                        !   2 - Double precision arithmetic
   real(4) :: A(MTOT)
   real(8) :: DA(MTOT/ITWO)
   integer :: IA(MTOT)

   equivalence (A,IA), (A,DA)  ! A, DA, and IA share the same storage units

   integer, parameter :: amax = 200    ! Maximum number of arrays allowed

   integer :: np(amax) = 0    ! Pointer to each array
   integer :: alen(amax) = 0  ! Length of each array
   integer :: aprec(amax) = 0 ! Precision of each array
   character*8 :: aname(amax) = ""

   integer :: nplast = 0      ! Pointer to the last allocated element in A
                               ! nplast is in the unit of single precision

contains

   subroutine memalloc(num, name, len, prec)
! -----------------------------------------------------------------------------
! -  Purpose                                                                  -
! -     Allocate an array in the storage of A                                 -
! -                                                                           -
! -  Input                                                                    -
! -     num  - Number of the array allocated                                  -
! -     name - Name of the array                                              -
! -     len  - Length of the array (total number of elements of the array)    -
! -     prec - Precision of the array                                         -
! -            1: Single precision                                            -
! -            2 : Double precesion                                           -
! -                                                                           -
! -----------------------------------------------------------------------------
      implicit none
      integer :: num, len, prec
      character*5 name
      
      integer :: i, npfirst

      if (num < 1 .or. num > amax) then
         write(*,'("*** Error *** Invalid array number:  ",I3)') num
         stop
      end if

      if (prec < 1 .or. prec > 2) then
         write(*,'("*** Error *** Invalid array type:  ",I3)') prec
         stop
      end if

      if (np(num) > 0) call memfree(num)  ! array num exists

      if (nplast+len*prec > MTOT) then
         write(*,'("*** Error *** No adequate storage available in A",/, &
                   "    Required  :", I10, /,  &
                   "    Available :", I10)') len*prec, MTOT - nplast
         stop
      end if

      npfirst = nplast + 1
      np(num) = nplast/prec + 1   ! In the unit of allocated array
      aname(num) = name
      alen(num)  = len
      aprec(num) = prec

      nplast = nplast + len*prec
      if (mod(nplast,2) == 1) nplast = nplast+1  ! Make nplast an even number
      
      do i = npfirst, nplast
         A(i) = 0
      end do

   end subroutine memalloc


   subroutine memfree(num)
! -----------------------------------------------------------------------------
! -  Purpose                                                                  -
! -     Free the array num and compact the storage if necessary               -
! -                                                                           -
! -  Input                                                                    -
! -     num  - Number of the array to be deallocated                          -
! -                                                                           -
! -----------------------------------------------------------------------------
      implicit none
      integer :: i, num, npbase, nplen

      if (np(num) <= 0) return  ! The array has not been allocated

!     Base address of the array num in the single precision unit
      npbase = (np(num)-1)*aprec(num)

!     Length of the array num in the single precision unit
      nplen  = ceiling(alen(num)*aprec(num)/2.0)*2  ! Make nplen an even number

!     Compact the storage if neccessary
      if (npbase+nplen < nplast) then
!        Move arrays behind the array num forward to reuse its storage
         do i = npbase+nplen+1, nplast
            A(i-nplen) = A(i)
         end do

!        Update the pointer of arrays behind the array num
         do i = 1, amax
            if ((np(i)-1)*aprec(i) > npbase) np(i) = np(i) - nplen/aprec(i)
         end do
      end if

      np(num)    = 0
      aname(num) = ""
      alen(num)  = 0
      aprec(num) = 0

      nplast = nplast - nplen
   end subroutine memfree


   subroutine memfreefrom(num)
! -----------------------------------------------------------------------------
! -  Purpose                                                                  -
! -     Free all arrays from num to the end                                   -
! -                                                                           -
! -  Input                                                                    -
! -     num  - Number of the array to be deallocated from                     -
! -                                                                           -
! -----------------------------------------------------------------------------
      implicit none
      integer :: i, num
      
      do i=amax,num,-1
         call memfree(i)
      end do
      
   end subroutine memfreefrom


   subroutine memfreefromto(n1,n2)
! -----------------------------------------------------------------------------
! -  Purpose                                                                  -
! -     Free all arrays from n1 to n2                                         -
! -                                                                           -
! -  Input                                                                    -
! -     n1  - Number of the array to be deallocated from                      -
! -     n2  - Number of the array to be deallocated to                        -
! -                                                                           -
! -----------------------------------------------------------------------------
      implicit none
      integer :: i, n1, n2
      
      do i=n2,n1,-1
         call memfree(i)
      end do
      
   end subroutine memfreefromto


   subroutine memprint(num)
! -----------------------------------------------------------------------------
! -  Purpose                                                                  -
! -     Print the contents of the array num                                   -
! -                                                                           -
! -  Input                                                                    -
! -     num  - Number of the array to be printed                              -
! -                                                                           -
! -----------------------------------------------------------------------------
      implicit none
      integer :: num,i

      if (np(num) <= 0) then
         write(*,'("*** Error *** Array ", I3, " has not been allocated.")') num
         return
      end if

      write(*,'("Contents of Array ", A5, ":")') aname(num)
      if (aprec(num) == 1) then
         write(*,'(8I10)') (IA(i), i=np(num),np(num)+alen(num)-1)
      else
         write(*,'(8E10.2)') (DA(i), i=np(num),np(num)+alen(num)-1)
      end if

   end subroutine memprint


   subroutine memprintptr(ptr, len, atype)
! -----------------------------------------------------------------------------
! -  Purpose                                                                  -
! -     Print the contents of the stroage starting from ptr                   -
! -                                                                           -
! -  Input                                                                    -
! -     ptr   - Pointer to the first entry (in single precision unit)         -
! -     len   - Total number of entries to be printed                         -
! -     atype - Type of the entries (0 - integer; 1 - float;  2 - double)     -
! -                                                                           -
! -----------------------------------------------------------------------------
      implicit none
      integer :: i, ptr, len, atype
      character*8 dtype(3)
      data dtype/"integer","real","double"/

      write(*,'("Contents of storage starting from ", I5, " in ", A8, ":")') ptr, dtype(atype+1)
      if (atype == 0) then
         write(*,'(8I10)') (IA(i), i=ptr,ptr+len-1)
      else if (atype == 1) then
         write(*,'(8E10.2)') (A(i), i=ptr,ptr+len-1)
      else if (atype == 2) then
         write(*,'(8E10.2)') (DA(i), i=(ptr-1)/ITWO+1, (ptr-1)/ITWO+len)
      end if

   end subroutine memprintptr


   subroutine meminfo
! -----------------------------------------------------------------------------
! -  Purpose                                                                  -
! -     Print the information of the storage                                  -
! -                                                                           -
! -----------------------------------------------------------------------------
      implicit none
      integer :: i

      write(*,'("List of all arrays:")')
      write(*,'("   Number   Name   Length   Pointer   Precision")')
      do i=1,amax
         if (np(i) == 0) cycle
         write(*,'(I7, 4X, A5, I9, I10, I12)') i, aname(i), alen(i), np(i), aprec(i)
      end do
   end subroutine meminfo

end module memAllocate
