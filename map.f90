program map
use random_gen
implicit none

    integer, parameter :: LSIZE=1000
    real, dimension(LSIZE,LSIZE) :: random_terrain, elevation
    real :: corr_length, variance, distance
    integer :: i,j,l,m

    call init_random_seed

    corr_length=100
    variance = 10


    do i = 1,LSIZE
        do j = 1,LSIZE

            random_terrain(i,j) = rnd_gauss(0.0, variance)
            elevation(i,j) = 0

        enddo
    enddo


    do i = 1,LSIZE
        do j=1,LSIZE

            !write(6,*) i,j

             do l = i - 2*corr_length, i + 2*corr_length
                do m = j - 2*corr_length, j + 2*corr_length
                    
                    distance = sqrt( real((l-i)*(l-i) + (j-m)*(j-m)) )

                    if(distance.lt.2*corr_length) then 
                        elevation(i,j) = elevation(i,j) + random_terrain(l,m)*exp(-distance/corr_length)
                    endif

                enddo
            enddo
        enddo
    enddo


    do i = 1,LSIZE
        do j = 1,LSIZE

            write(6,*) i,j,random_terrain(i,j),elevation(i,j)            

        enddo
    
        write(6,*) ''
    
    enddo


end program

    
