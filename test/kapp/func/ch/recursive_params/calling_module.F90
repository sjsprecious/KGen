module calling_module

        USE kernel, only : add

        public calling_subroutine

        contains

        subroutine calling_subroutine()

                integer(kind=4), parameter :: real_kind = 4

                real(kind=real_kind), allocatable :: array_1(:,:,:,:), &
                                                        array_2(:,:,:,:), &
                                                        array_3(:,:,:,:)

                integer(kind=4), parameter :: new_size=9
                integer(kind=4), parameter :: temp_size=new_size+1, array_size = temp_size


                allocate(array_1(array_size,array_size,array_size,array_size))
                allocate(array_2(array_size,array_size,array_size,array_size))
                allocate(array_3(array_size,array_size,array_size,array_size))

                array_1(:,:,:,:) = 1.0
                array_2(:,:,:,:) = 1.0

                call add(array_1, array_2, array_3)

                !print *, "a3 = ", array_3

                if(allocated(array_1)) deallocate(array_1)
                if(allocated(array_2)) deallocate(array_2)
                if(allocated(array_3)) deallocate(array_3)

        end subroutine

end module

                
