module schemas_DF

    use fonctions
    use mod_algebre

    implicit none

    contains

        ! subroutine test_f_fp()
        !     real(PR) :: z
        !     write(*,*) "z", "f(z)", "fp(z)"
        !     do z = 0.0_PR, 1.0_PR, 0.1_PR
        !         write(*,*) z, f(z), fp(z)
        !     end do
        ! end subroutine test_f_fp

    subroutine Calc_DF()
        real(PR), allocatable, dimension(:,:) :: L, M_LU
        logical :: test
        integer :: nxa
    
        call fill_mat(A)
        call fill_vect(b)

        nxa = size(A,1)
        test = .true.
    
        if (schx == 1) then
            ! Le schéma à 3 points est symétrique : Cholesky est ok
            allocate(L(nxa,nxa))
            call chol_decomposition(A, L, test)
            if (test) then
                call chol_res(L, b, u(1:nx), test)
            end if
            deallocate(L)
        else if (schx == 2) then
            ! Le schéma à 5 points est asymétrique : Utiliser LU
            allocate(M_LU(nxa,nxa))
            call lu_decomposition(A, M_LU, test)
            if (test) then
                call lu_res(M_LU, b, u(1:nx), test)
            end if
            deallocate(M_LU)
        end if
    
        if (.not. test) print*, "Erreur lors de la résolution du système !"
    end subroutine Calc_DF
    
end module schemas_DF
