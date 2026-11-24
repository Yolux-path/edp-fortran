module post_traitement

    use fonctions

    implicit none
    

    contains

        subroutine Calcul_erreurs(e)

            real(PR), intent(out) :: e(3)
            real(PR)              :: e_inf, e_1, e_2
            integer               :: i

            e_inf = 0.0_PR
            do i = 0, nx-1
                e_inf = max(e_inf, abs(u(i) - fp(x(i))))
            end do

            e_1 = 0.0_PR
            do i = 0, nx-1
                e_1 = e_1 + (dx * abs(u(i) - fp(x(i))) ) 
            end do
            e_1 = e_1 / Lx

            e_2 = 0.0_PR
            do i = 0, nx-1
                e_2 = e_2 + (dx * (u(i) - fp(x(i)))**2)
            end do
            e_2 = sqrt(e_2)
            e_2 = e_2 / Lx

            e(1) = e_inf
            e(2) = e_1
            e(3) = e_2

        end subroutine Calcul_erreurs
        
        subroutine ecrit_erreurs(tab_dx,err)
            real(PR), allocatable, dimension(:), intent(in)     :: tab_dx
            real(PR), allocatable, dimension(:,:), intent(in)   :: err
            integer :: j

            open(unit=20, file="err.dat", action="write")

            do j = 1, size(tab_dx)
                write(20,*) tab_dx(j), err(:,j)
            end do

            close(20)
        
        end subroutine ecrit_erreurs

        subroutine affiche_pentes(tab_dx, err)
            real(PR), dimension(:), intent(in)     :: tab_dx
            real(PR), dimension(:,:), intent(in)   :: err
            real(PR) :: pente
            integer  :: j
        
            print*, "--- Analyse des pentes locales (Convergence) ---"
            print*, " j  |  Delta x  |  Erreur (L-inf) |  Pente locale"
            print*, "------------------------------------------------"
            
            do j = 1, size(tab_dx) - 1
                ! Calcul de la pente locale pour l'erreur L-inf (err(1,j))
                pente = (log(err(1, j+1)) - log(err(1, j))) / &
                        (log(tab_dx(j+1)) - log(tab_dx(j)))
                
                write(*, '(I3, " | ", E10.3, " | ", E12.5, " | ", F6.2)') &
                      j, tab_dx(j), err(1,j), pente
            end do
        end subroutine affiche_pentes
    
end module post_traitement