subroutine bdopt(ID)
    USE LIST_CLASS
    USE NODE_CLASS
    USE GLOBALS, ONLY : IIN, IOUT, NEQ, NUMNP, NWK
    
    implicit none
    TYPE(LIST) :: LISTS(NEQ)
    TYPE(NODE), POINTER :: P_NODE
    INTEGER :: ID(6, NUMNP), ID_NEW
    INTEGER :: NUME ! 全部总共的节点个数
    INTEGER :: N    ! 节点中的函数
    INTEGER :: TEMP_NODE(8) ! 用来暂存
    INTEGER :: TEMP_ID(48) ! 用来暂存节点信息
    INTEGER :: I, J, K ! 循环变量
    INTEGER :: temp_length, TEMP_INDEX
    INTEGER :: FREE_DOF ! 单元中自由的自由度
    character (len = 10) :: cFmt
    DO I = 1,NUMNP
        DO J = 1,6
            IF(ID(J,I) .NE. 0) THEN
                !CALL SET(lists(ID(J,I)), I, J)
            END IF
        END DO
    END DO
    ! input phase

     READ (IIN,"(I10)") NUME ! 总element数
     DO I = 1,NUME
         READ (IIN,"(I10)") N ! 这个element对应的节点数
         READ (IIN, "(8I10)") (TEMP_NODE(J), J = 1,N)
         FREE_DOF = 0
         DO J = 1,N
             DO K = 1,6
                 IF (ID(K, TEMP_NODE(J)) .NE. 0) THEN
                     FREE_DOF = FREE_DOF + 1
                     TEMP_ID(FREE_DOF) = ID(K, TEMP_NODE(J))
                 END IF
             END DO
         END DO
         DO J = 1, FREE_DOF
             DO K = 1, FREE_DOF
                 CALL NEW(P_NODE, TEMP_ID(J))
                 call ADD_WITH_SEARCH(lists(TEMP_ID(K)), P_NODE)
             END DO
         END DO
     END DO
        ! CM algorithm
        ! 先找到最小的节点
        ID_NEW = 1
        temp_length = LISTS(1)%length_
        TEMP_INDEX = 1
        DO I = 1,NEQ
             IF(temp_length > lists(I)%length_) THEN
                temp_length = lists(I)%length_
                TEMP_INDEX = I
            END IF
        END DO
        !ID(lists(TEMP_INDEX)%column_, lists(TEMP_INDEX)%row_) = ID_NEW
        ID_NEW = ID_NEW + 1
        ! 后续的搜索循环
        p_node => lists(temp_index)%head_
        do while(associated(p_node))
            lists(p_node%index_)%sign_ = 1
            p_node => p_node%next_
        end do
        lists(temp_index)%sign_ = -1
        do while( id_new .le. neq)
            temp_index = 0
            temp_length = 100000 ! 应该要放最大整数，但是忘了是多少了
            do i = 1, neq
              if(lists(i)%sign_ .eq. 1) then
                  if(temp_length > lists(i)%length_) then
                    temp_length = lists(i)%length_
                    temp_index  = i
                  end if
             end if
           end do
           !id(lists(temp_index)%column_, lists(temp_index)%row_) = id_new
           id_new = id_new+1
           p_node => lists(temp_index)%head_
           do while(associated(p_node))
              if(lists(p_node%index_)%sign_ .ne. -1) then
                  lists(p_node%index_)%sign_ = 1
              end if
              p_node => p_node%next_
           end do
           lists(temp_index)%sign_ = -1
        end do
        ! Write equation numbers
        WRITE (IOUT,"(//,' EQUATION NUMBERS',//,'   NODE',9X,  &
                     'DEGREES OF FREEDOM',/,'  NUMBER',/,  &
                     '     N',13X,'X    Y    Z   RX   RY   RZ',/,(1X,I10,9X,6I10))") (N,(ID(I,N),I=1,6),N=1,NUMNP)
     !  清空list，释放内存
     do i = 1, neq
         call delete_all(lists(i))
     end do
end subroutine bdopt
    
subroutine pardiso_input(ID)
    USE vector_CLASS
    USE GLOBALS, ONLY : IIN, IOUT, NEQ, NUMNP, NWK, pardisodoor, TIM, huge
    USE MEMALLOCATE
    
    implicit none
    TYPE(vector) :: vectors(NEQ)
    INTEGER :: ID(6, NUMNP), ID_NEW
    INTEGER :: NUME ! 全部总共的节点个数
    INTEGER :: N    ! 节点中的函数
    INTEGER :: TEMP_NODE(8) ! 用来暂存
    INTEGER :: TEMP_ID(48) ! 用来暂存节点信息
    INTEGER :: I, J, K ! 循环变量
    INTEGER :: temp_length, TEMP_INDEX
    INTEGER :: FREE_DOF ! 单元中自由的自由度
    ! input phase
     READ (IIN,"(I10)") NUME ! 总element数
     ! 为了节省空间，这里用sign_表示序号好了
     WRITE(*,'("Begin create vectors ")')
     DO I = 1,NUME
         READ (IIN,"(I10)") N ! 这个element对应的节点数
         READ (IIN, "(<N>I10)") (TEMP_NODE(J), J = 1,N)
         FREE_DOF = 0
         ! 首先确认vector的长度
         DO J = 1,N
             DO K = 1,6
                 if(ID(K, TEMP_NODE(J)) .ne. 0) then
                     select case(N)
                     case(2)
                        vectors(ID(K, TEMP_NODE(J)))%length_ = vectors(ID(K, TEMP_NODE(J)))%length_+11
                     case(4)
                        vectors(ID(K, TEMP_NODE(J)))%length_ = vectors(ID(K, TEMP_NODE(J)))%length_+19
                     case(8)
                        vectors(ID(K, TEMP_NODE(J)))%length_ = vectors(ID(K, TEMP_NODE(J)))%length_+23
                    end select
                 end if
             END DO
        END DO
     END DO
    WRITE(*,'("Begin allocate vectors ")')
    do i = 1,neq
        call set(vectors(i), i)
    end do
    WRITE(*,'("End allocate vectors ")')
    ! 向vector中加入节点
     DO I = 1,NUME
         READ (IIN,"(I10)") N ! 这个element对应的节点数
         READ (IIN, "(<N>I10)") (TEMP_NODE(J), J = 1,N)
         FREE_DOF = 0
         DO J = 1,N
             DO K = 1,6
                     FREE_DOF = FREE_DOF + 1
                     TEMP_ID(FREE_DOF) = ID(K, TEMP_NODE(J))
             END DO
         END DO
         DO K = 1, FREE_DOF
             DO J = 1, FREE_DOF
                 if((temp_id(K) .ne. 0).and.(temp_id(J) .ne. 0)) then
                    call add(vectors(TEMP_ID(K)), temp_id(j))
                end if
             END DO
         END DO
     END DO    
     WRITE(*,'("End create vectors ")')
     CALL SECOND (TIM(2))
     ! 注意这里分配了rowIndex
     WRITE(*,'("Begin calculating rowIndex ")')
     CALL MEMALLOC(2,"rowIndex",NEQ+1,1)
     CALL assign_rowIndex(vectors, IA(NP(2)))
     WRITE(*,'("End calculating rowIndex ")')
     if( MTOT < nplast + nwk*2 + neq) then ! Job-4
         huge = .true.
     end if
     WRITE(*,'("Begin calculating columns ")')
     if(huge) then
        CALL MEMALLOC(3,"STFF ",1,ITWO)
        CALL MEMALLOC(4,"R    ",NEQ,ITWO)
        CALL MEMALLOC(5,"columns",1,1)
        allocate(stff(nwk))
        allocate(columns(nwk))
        CALL assign_columns(vectors, columns)
     else
         !huge = .false.
        CALL MEMALLOC(3,"STFF ",NWK,ITWO)
        CALL MEMALLOC(4,"R    ",NEQ,ITWO)
        CALL MEMALLOC(5,"columns",NWK,1)
        CALL assign_columns(vectors, IA(NP(5)))
     end if
     WRITE(*,'("End calculating columns ")')
     WRITE(*,'("Begin delete vectors ")')
     do i = 1, neq
         call delete(vectors(i))
     end do
     WRITE(*,'("End delete vectors ")')
    end subroutine pardiso_input

subroutine assign_rowIndex(vectors, rowIndex)
    use GLOBALS, only : neq, nwk
    use vector_class
    
    implicit none
    type(vector) :: vectors(neq)
    integer :: rowIndex(neq+1)
    integer :: i, j
    
    nwk = 0
    do i=1, neq
        rowIndex(i) = 1 + nwk
        nwk = nwk + vectors(i)%length_real_
    end do
    rowIndex(neq+1) = nwk+1
    end subroutine assign_rowIndex
    
subroutine assign_columns(vectors, columns)
    use GLOBALS, only : neq, nwk
    use vector_class
    
    implicit none
    type(vector) :: vectors(neq)
    integer :: columns(nwk)
    integer :: i, j, k
    
    k = 1
    do i = 1, neq
        do j = 1, vectors(i)%length_real_
            columns(k) = vectors(i)%nodes_(j)
            k = k+1
        end do
    end do
end subroutine assign_columns