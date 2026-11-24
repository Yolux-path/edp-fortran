module donnees

    implicit none

    integer, parameter                  ::  PR = 8

    integer                             ::  nx, schx, choix_f, nr, nx0

    real(PR)                            ::  Lx, dx

    real(PR), allocatable, DIMENSION(:) ::  u, x

    contains

        subroutine Initialisation() !Validé

            integer :: i = 0

            open(unit=10, file="parametres.dat", status="old", action="read")

            read(10,*) Lx, nx0, schx, choix_f, nr !Longueur du domaine, nombre de points, schéma de discrétisation, fonction choisie, nombre de raffinements
    
            nx = nx0

            dx = Lx / real(nx, PR)

            allocate(u(0:nx-1), x(0:nx-1))

            x = [(i*dx, i=0,nx-1)]
    
            close(10)


        end subroutine Initialisation
        
        subroutine Fin()

            deallocate(u, x)

        end subroutine Fin

        subroutine raffinement(j)
            
            integer, intent(in)     :: j
            integer                 :: i
            
            if(allocated(u) .or. allocated(x)) then
                deallocate(u)
                deallocate(x)
            end if

            nx = (2**j) * nx0
            dx = Lx / real(nx, PR)

            allocate(u(0:nx-1), x(0:nx-1))

            do i = 0, nx-1
                x(i) = i*dx
            end do
            
        end subroutine raffinement

end module donnees
