module post_traitement
    use fonctions
    implicit none
    
contains

    subroutine Calcul_erreurs(e, approx)
        real(PR), dimension(:), allocatable, intent(in) :: approx
        real(PR), intent(inout) :: e(3)
        real(PR) :: e_inf, e_1, e_2
        integer  :: i

        e_inf = 0.0_PR
        do i = 0, nx+1
            e_inf = max(e_inf, abs(Ue(i) - approx(i)))
        end do
   
        e_1 = 0.0_PR
        do i = 0, nx+1
            e_1 = e_1 + (dx * abs(Ue(i) - approx(i))) 
        end do
        
        e_2 = 0.0_PR
        do i = 0, nx+1
            e_2 = e_2 + (dx * (Ue(i) - approx(i))**2)
        end do
        e_2 = sqrt(e_2) / Lx

        e(1) = e_inf
        e(2) = e_1
        e(3) = e_2
    end subroutine Calcul_erreurs
    
    subroutine ecrit_erreurs(tab_dx, err)
        real(PR), allocatable, dimension(:), intent(in)     :: tab_dx
        real(PR), allocatable, dimension(:,:), intent(in)   :: err
        integer  :: j, k
        real(PR) :: p
    
        open(unit=30, file="erreurs.dat", status="replace")
        write(30,*) "# dx | E_inf | Pente | E_L1 | Pente | E_L2 | Pente"
    
        do j = 1, size(tab_dx)
            write(30, '(E15.8)', advance='no') tab_dx(j)
            do k = 1, 3
                if (j > 1) then
                    p = log(err(k, j-1) / err(k, j)) / log(tab_dx(j-1) / tab_dx(j))
                    write(30, '(2X, E15.8, 2X, F8.3)', advance='no') err(k, j), p
                else
                    write(30, '(2X, E15.8, 2X, A8)', advance='no') err(k, j), "----"
                end if
            end do
            write(30, *) ""
        end do
        close(30)
    end subroutine ecrit_erreurs

end module post_traitement