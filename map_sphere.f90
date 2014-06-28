program map
use random_gen
use constants
implicit none

    integer, parameter :: LSIZE=100
    real, dimension(LSIZE,LSIZE) :: random_terrain, elevation, normalization
    real :: corr_length, variance, distance
    integer :: i,j,l,m
    integer :: rad, rad_m
    real :: phi,x,y,z !real-space positions of lattice points of point (i,j)
    real :: phi_m,x_m,y_m,z_m !same for points (l,m)

    call init_random_seed

    corr_length=10
    variance = 10


!initialize arrays
    random_terrain = 0
    elevation = 0
    normalization = 0


!fill the necessary part of the lattice with random numbers
    do j = 1,LSIZE
    
        rad = sqrt(real(LSIZE*LSIZE/4 - (j-LSIZE/2)*(j-LSIZE/2)))

        do i = int(0.5*LSIZE-rad),int(0.5*LSIZE+rad)

            random_terrain(i,j) = rnd_gauss(0.0, variance)          

        enddo
    enddo



!integrate over a kernel to introduce correlations
    do j = 1,LSIZE

        rad = int(sqrt(real(LSIZE*LSIZE/4 - (j-LSIZE/2)*(j-LSIZE/2))))
   
        do i = int(0.5*LSIZE-rad),int(0.5*LSIZE+rad)

            if (rad.eq.0) rad=1
  
            phi = 2*pi*real(i-0.5*LSIZE+rad)/real(2*rad)
            z = (j-LSIZE/2)  
            x = rad*cos(phi)  
            y = rad*sin(phi)

            do m = 1,LSIZE
 
                rad_m = sqrt(real(LSIZE*LSIZE/4 - (m-LSIZE/2)*(m-LSIZE/2))) 

                if (rad.eq.0) rad=1

                do l = int(0.5*LSIZE-rad),int(0.5*LSIZE+rad)
                
                    phi_m = 2*pi*real(l-0.5*LSIZE+rad_m)/real(2*rad_m)
                    z_m = (m-LSIZE/2)
                    x_m = rad_m*cos(phi_m)
                    y_m = rad_m*sin(phi_m) 
                    distance = sqrt( real((x-x_m)*(x-x_m) + (y-y_m)*(y-y_m) + (z-z_m)*(z-z_m)) )

 
                    if(distance<2*corr_length) then
                        elevation(i,j) = elevation(i,j) + random_terrain(l,m)*exp(-distance/corr_length)  
                        normalization(i,j) = normalization(i,j) + exp(-distance/corr_length)      
                    endif   

                 enddo

            enddo
           

        enddo
    enddo


!print out the grid
    do j = 1,LSIZE
        rad = sqrt(real(LSIZE*LSIZE/4 - (j-LSIZE/2)*(j-LSIZE/2)))

        if (rad.eq.0) rad=1

        do i = int(0.5*LSIZE-rad),int(0.5*LSIZE+rad)

            phi = 2*pi*real(i-0.5*LSIZE+rad)/real(2*rad)
            z = (j-LSIZE/2)
            x = rad*cos(phi)
            y = rad*sin(phi)

            write(6,*) i,j,random_terrain(i,j),elevation(i,j),normalization(i,j),x,y,z,phi        

        enddo
    
        write(6,*) ''
    
    enddo


end program

    
