        program transfer

        integer :: i,j
        real :: h(401,401),h1(401,401),h2(401,401)

        open(11,file='height.txt',form='formatted')

        read(11,*) h

        close(11)

        open(12,file='h.dat',form='unformatted',access='direct',recl=4*401*401)

        write(12,rec=1) h

        close(12) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        open(11,file='h1.txt',form='formatted')

        read(11,*) h1

        close(11)

        open(12,file='h1.dat',form='unformatted',access='direct',recl=4*401*401)

        write(12,rec=1) h1

        close(12) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(11,file='h2.txt',form='formatted')

        read(11,*) h2

        close(11)

        open(12,file='h2.dat',form='unformatted',access='direct',recl=4*401*401)

        write(12,rec=1) h2

        close(12) 

        stop

        end
