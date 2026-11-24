module donnees
    implicit none
    integer, parameter                      :: PR = 8
    integer                                 :: nx, schx, choix_f, nr, nx0
    real(PR)                                :: Lx, dx, dt, a, cfl, c
    real(PR), allocatable, DIMENSION(:)     :: Ue, x, b
    real(PR), allocatable, DIMENSION(:,:)   :: M
    real(PR), parameter :: D_diff = 1.0e-3_PR

contains

    subroutine Initialisation()
        integer :: i
        open(unit=10, file="parametres.dat", status="old", action="read")
        read(10,*) Lx, nx0, schx, choix_f, nr, a, cfl
        close(10)

        nx = nx0
        dx = Lx / real(nx+1, PR)
        dt = cfl * (dx**2) / (abs(a)*dx + 2.0_PR*D_diff)
        c = (a * dt) / dx

        allocate(Ue(0:nx+1), x(0:nx+1), b(nx), M(nx,nx))
        x = [(real(i, PR)*dx, i=0, nx+1)]
    end subroutine Initialisation

    subroutine raffinement(j)
        integer, intent(in) :: j
        integer             :: i
        if(allocated(x)) deallocate(Ue, x, b, M)
        
        nx = nx0 * (2**j) 
        dx = Lx / real(nx+1, PR)
        dt = cfl * (dx**2) / (abs(a)*dx + 2.0_PR*D_diff)
        c = (a * dt) / dx
    
        allocate(Ue(0:nx+1), x(0:nx+1), b(nx), M(nx,nx))
        do i = 0, nx+1
            x(i) = real(i, PR) * dx
        end do
    end subroutine raffinement

    subroutine Fin()
        if(allocated(Ue)) deallocate(Ue, x, b, M)
    end subroutine Fin
end module donnees