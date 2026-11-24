program advection
    use donnees
    use fonctions
    use post_traitement
    use mod_schemas
    implicit none 

    integer                                 :: imax, i, i_ref, n_it
    real(PR)                                :: t_final, tnn, dt_init
    real(PR), dimension(:), allocatable     :: Un, Unp1, tab_dx
    real(PR), allocatable, dimension(:,:)   :: err

    call Initialisation()
    allocate(err(3, nr), tab_dx(nr))

    do i_ref = 1, nr
        if (i_ref > 1) call raffinement(i_ref-1)
        
        imax = nx
        tab_dx(i_ref) = dx
        dt_init = dt
        
        if (allocated(Un)) deallocate(Un, Unp1)
        allocate(Un(0:imax+1), Unp1(0:imax+1))

        t_final = 25.0_PR 
        tnn = 0.0_PR 
        n_it = 0
        
        do i = 0, imax+1
            Un(i) = Uinit(0.0_PR, x(i))
        end do
   
        ! Ouverture unique du fichier pour toute la simulation
        open(unit=40, file="solution.dat", status="replace")
        
        do while (tnn < t_final)
            if (tnn + dt_init > t_final) dt_init = t_final - tnn
            
            n_it = n_it + 1
            c = (a * dt_init) / dx

            Un(0) = Ug(tnn)

            do i = 1, imax
                select case (schx)
                case (1)
                    Unp1(i) = amont(Un(i), Un(i-1))
                case (2)
                    Unp1(i) = centre(Un(i), Un(i-1), Un(i+1))
                case (3)
                    Unp1(i) = lax_wendroff(Un(i), Un(i-1), Un(i+1))
                case default
                    stop "schx inconnu (1=amont, 2=centre, 3=Lax-Wendroff)"
                end select
                Unp1(i) = Unp1(i) + (dt_init * D_diff / dx**2) * (Un(i+1) - 2.0_PR*Un(i) + Un(i-1))
            end do
            
            Unp1(0) = Ug(tnn + dt_init)
            Unp1(imax+1) = Unp1(imax)
            Un = Unp1
            tnn = tnn + dt_init

            ! Sauvegarde toutes les 100 itérations
            if (mod(n_it, 100) == 0) then 
                write(40,*) "" ! Ligne vide 1
                write(40,*) "" ! Ligne vide 2 (Standard Gnuplot pour 'index')
                do i = 0, imax+1
                    write(40, '(3ES17.8E3)') x(i), Un(i), tnn
                end do
            end if
        end do
        close(40)

        do i = 0, imax+1
            Ue(i) = Uexacte(tnn, x(i))
        end do
        call Calcul_erreurs(err(:, i_ref), Un)
    end do 

    call ecrit_erreurs(tab_dx, err)
    call Fin()
end program advection