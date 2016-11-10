program mpmTest
   use memAllocate

   integer i, n1, n2, n3, n4, n5

   write(*,'(/, "Create array ", I3, ":", A5)') 5, "X"
   call memalloc(5, "X    ", 20, 2)
   do i=np(5), np(5)+alen(5)-1
      DA(i) = i
   end do
   call memprint(5)

   write(*,'(/, "Create array ", I3, ":", A5)') 1, "Y"
   call memalloc(1, "Y    ", 10, 2)
   do i=np(1), np(1)+alen(1)-1
      DA(i) = 20+i
   end do
   call memprint(1)

   write(*,'(/, "Create array ", I3, ":", A5)') 11, "IX"
   call memalloc(11, "IX   ", 30, 1)
   do i=np(11), np(11)+alen(11)-1
      IA(i) = 30+i
   end do
   call memprint(11)

   write(*,'(/, "Create array ", I3, ":", A5)') 15, "IY"
   call memalloc(15, "IY   ", 5, 1)
   do i=np(15), np(15)+alen(15)-1
      IA(i) = 15+i
   end do
   call memprint(15)

   call meminfo

   write(*,'(/, "Delete array ", I3, ":", A5)') 5, "X"
   call memfree(5)
   call meminfo

   call memprint(1)
   call memprint(11)
   call memprint(15)

   write(*,'(/, "Delete array ", I3, ":", A5)') 11, "IX"
   call memfree(11)
   call meminfo

   call memprint(1)
   call memprint(15)

   write(*,'(/, "Create array ", I3, ":", A5)') 20, "D"
   call memalloc(20, "D    ", 7, 2)
   call meminfo
   do i=np(20), np(20)+alen(20)-1
      DA(i) = 15+i
   end do
   call memprint(20)

   write(*,'(/, "Delete array ", I3, ":", A5)') 15, "IY"
   call memfree(15)
   call meminfo
   call memprint(1)
   call memprint(20)

   write(*,'(/, "Delete array ", I3, ":", A5)') 1, "Y"
   call memfree(1)
   call meminfo
   call memprint(20)

   write(*,'(/, "Create array ", I3, ":", A5)') 5, "X"
   call memalloc(5, "X    ", 20, 2)
   do i=np(5), np(5)+alen(5)-1
      DA(i) = i
   end do
   call memprint(5)

   write(*,'(/, "Create array ", I3, ":", A5)') 1, "Y"
   call memalloc(1, "Y    ", 10, 2)
   do i=np(1), np(1)+alen(1)-1
      DA(i) = 20+i
   end do
   call memprint(1)

   write(*,'(/, "Create array ", I3, ":", A5)') 11, "IX"
   call memalloc(11, "IX   ", 30, 1)
   do i=np(11), np(11)+alen(11)-1
      IA(i) = 30+i
   end do
   call memprint(11)

   write(*,'(/, "Create array ", I3, ":", A5)') 15, "IY"
   call memalloc(15, "IY   ", 5, 1)
   do i=np(15), np(15)+alen(15)-1
      IA(i) = 15+i
   end do
   call memprint(15)

   call meminfo

   write(*,'(/, "Delete array ", I3, ":", A5)') 11, "IX"
   call memfree(11)
   call meminfo
   call memprint(1)
   call memprint(5)
   call memprint(15)
   call memprint(20)

   write(*,'(/, "Delete array ", I3, ":", A5)') 20, "D"
   call memfree(20)

   write(*,'(/, "Delete array ", I3, ":", A5)') 5, "X"
   call memfree(5)

   write(*,'(/, "Delete array ", I3, ":", A5)') 1, "Y"
   call memfree(1)

   call meminfo
   call memprint(15)

end