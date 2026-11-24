module donnees

    implicit none

    integer, parameter                      ::  PR = 8

    integer                                 ::  nx, schx, choix_f, nr, nx0

    real(PR)                                ::  Lx, dx, dt, Tg, Td, D, cfl

    real(PR), allocatable, DIMENSION(:)     ::  T, x, b
    real(PR), allocatable, DIMENSION(:,:)   ::  A

    contains

        subroutine Initialisation()
            integer :: i = 0

            open(unit=10, file="parametres.dat", status="old", action="read")
            read(10,*) Lx, nx0, schx, choix_f, nr, Tg, Td, D, cfl 
            close(10)

            nx = nx0
            dx = Lx / real(nx+1, PR)
            
            ! Calcul correct du dt avec le facteur 2.0
            dt = cfl * (dx**2 / (2.0_PR * D))

            allocate(T(0:nx+1), x(0:nx+1), b(nx), A(nx,nx))

            T(0) = Tg
            T(nx+1) = Td   

            x = [(i*dx, i=0, nx+1)]

        end subroutine Initialisation
        
        subroutine Fin()
            if (allocated(T)) deallocate(T, x, b, A)
        end subroutine Fin

        subroutine raffinement(j)
            integer, intent(in)     :: j
            integer                 :: i
            
            if(allocated(T) .or. allocated(x)) then
                deallocate(T)
                deallocate(x)
                deallocate(A)
                deallocate(b)
            end if

            nx = (2**j) * nx0
            
            dx = Lx / real(nx+1, PR)
            
            ! dt dépend de dx
            dt = cfl * (dx**2 / (2.0_PR * D))

            allocate(T(0:nx+1), x(0:nx+1), b(1:nx) ,A(1:nx, 1:nx))

            T(0) = Tg
            T(nx+1) = Td 

            x = [(i*dx, i=0,nx+1)]
            
        end subroutine raffinement

end module donnees