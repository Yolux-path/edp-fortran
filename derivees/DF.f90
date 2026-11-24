program main
    
    use donnees
    use schemas_DF
    use post_traitement

    implicit none

    integer                                 :: j
    real(PR), allocatable, dimension(:)     :: tab_dx
    real(PR), allocatable, dimension(:,:)   :: err
    
    call Initialisation()

    allocate(tab_dx(nr), err(3, nr))
    
    ! call test_f_fp()

    ! print*, "schema", schx
    ! print*, "choix_f", choix_f

    do j = 1, nr
        call Calc_DF()
        call Calcul_erreurs( err(:,j) )
        tab_dx(j) = dx
        call raffinement(j)
    end do

    call ecrit_erreurs(tab_dx, err)
    
    call Fin
    
    deallocate(tab_dx, err)

end program main