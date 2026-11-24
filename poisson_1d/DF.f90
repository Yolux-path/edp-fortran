program main
    
    use donnees
    use mod_algebre
    use schemas_DF
    use fonctions
    use post_traitement

    implicit none

    integer                                 :: j
    real(PR), allocatable, dimension(:)     :: tab_dx
    real(PR), allocatable, dimension(:,:)   :: err
    
    call Initialisation()
    
    allocate(tab_dx(nr), err(3, nr))

    do j = 1, nr
        !print*, "Raffinement n°", j
        call Calc_DF()
        call Calcul_erreurs( err(:,j) )
        tab_dx(j) = dx
        if (j < nr) then
            call raffinement(j)
        end if
    end do

    open(unit=10, file="solution_DF.dat", status="replace")

    do j= 0, nx+1
        write(10,*) x(j), u(j), sol_exacte(x(j))
    end do

    close(10)

    call ecrit_erreurs(tab_dx, err)

    call Fin()

    deallocate(tab_dx, err)

    !TEST DU MODULE ALGEBRE

    ! real(PR), dimension(:,:), allocatable   :: A,L ,Lt
    ! real(PR), dimension(:), allocatable     :: p, b
    ! integer                                 :: n
    ! logical                                 :: test
    ! real(PR)                                :: h

    ! n = 3
    ! h = 1/real(n+1)
    ! test = .true.

    ! allocate(A(n,n), L(n,n), Lt(n,n), p(n), b(n))

    ! call fill_mat(A,h)

    ! call print_mat(A,"A")

    ! call chol_decomposition(A,L,test)

    ! Lt = transpose(L)

    ! call print_mat(L,"L")
    ! call print_mat(Lt,"Lt")

    ! ! différence relative entre A et LLT (max :  10-15) OK
    ! call print_mat((ABS(A - MATMUL(L,Lt))),"A-LLt")

    ! call random_number(b)

    ! call chol_res(L,b,p,test)

    ! ! différence relative entre Ap et b (max :  10-16) OK
    ! print*, MAXVAL(abs(MATMUL(A,p) - b))

    ! deallocate(A, L, p, b)

    

end program main