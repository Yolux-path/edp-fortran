module schemas_DF


    use fonctions
    
    implicit none

    contains

        ! subroutine test_f_fp()
        !     real(PR) :: z
        !     write(*,*) "z", "f(z)", "fp(z)"
        !     do z = 0.0_PR, 1.0_PR, 0.1_PR
        !         write(*,*) z, f(z), fp(z)
        !     end do
        ! end subroutine test_f_fp

        subroutine Calc_DF() ! Foncion de calcul du schéma de diff finie (Validé)
            
            integer                 :: i
            
            if(schx == 1) then
                u(0) = fp(x(0))
                u(nx-1) = fp(x(nx-1))
                do i = 1, nx-2
                    u(i) = ( f(x(i+1)) - f(x(i-1)) ) / (2.0_PR*dx) ! schéma centré
                end do
            else if (schx == 2) then
                u(0) = fp(x(0))
                u(nx-1) = fp(x(nx-1))
                do i = 1, nx-2
                    u(i) = ( f(x(i+1)) - f(x(i)) ) / dx 
                end do
            else
                u(0) = fp(x(0))
                u(1) = fp(x(1))
                u(nx-1) = fp(x(nx-1))
                do i = 2, nx-2
                    u(i) = ( 2*f(x(i+1)) + 3*f(x(i)) - 6*f(x(i-1)) + f(x(i-2)) ) / (6.0_PR*dx)
                end do
            end if

            ! print *, "u", "fp(x)", "dfférence entre u et fp(x)"
            ! do i = 0, nx-1
            !     print *, u(i), fp(x(i)), abs(u(i) - fp(x(i)))
            ! end do
        end subroutine Calc_DF
    
end module schemas_DF
