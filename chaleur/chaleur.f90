program chaleur
    use donnees
    use mod_algebre
    use schemas_chaleur
    use fonctions
    use post_traitement

    implicit none 

    integer                               :: imax, i, j, l
    real(PR)                              :: t1, tnn = 0
    real(PR), dimension(:), allocatable   :: tab_t, Tn, Tnp1
    
    ! Variables pour le schéma implicite
    real(PR), dimension(:,:), allocatable :: I_mat, M, L_chol
    logical                               :: check_chol

    ! Variables pour l'étude de convergence
    real(PR), dimension(:), allocatable   :: tab_dx
    real(PR), dimension(:,:), allocatable :: err_tab
    real(PR), dimension(3)                :: e

    character(len=40)                     :: ct
    character(len=2)                      :: prefix_sch

    call Initialisation()

    ! Choix du préfixe pour les fichiers de sortie
    if (schx == 1) then
        prefix_sch = "EE"
        print*, "=== Lancement avec le schéma EXPLICITE ==="
    else if (schx == 2) then
        prefix_sch = "EI"
        print*, "=== Lancement avec le schéma IMPLICITE ==="
    else
        print*, "ERREUR : schx inconnu (1=Explicite, 2=Implicite)"
        stop
    end if

    ! =========================================================================
    ! PARTIE 1 : Génération des profils temporels sur le maillage de base (nx0)
    ! =========================================================================
    print*, "-> Partie 1 : Génération des profils temporels..."
    imax = nx
    allocate(tab_t(0:5))
    allocate(Tn(0:imax+1), Tnp1(1:imax))
    tab_t = [0._PR, 0.02_PR, 0.04_PR, 0.06_PR, 0.08_PR, 0.10_PR]
    
    ! Préparation Implicite (si sélectionné)
    if (schx == 2) then
        allocate(I_mat(imax, imax), M(imax, imax), L_chol(imax, imax))
        I_mat = 0.0_PR
        do i = 1, imax
            I_mat(i,i) = 1.0_PR
        end do
        call fill_mat(A)
        M = I_mat + dt * A
        call chol_decomposition(M, L_chol, check_chol)
    end if

    ! Boucle sur les temps d'observation
    do j=0, 5
        t1 = tab_t(j)
        tnn = 0.0_PR
       
        ! Condition initiale
        do i=0,imax+1
            Tn(i) = Tinit(x(i))
        end do
        
        ! Boucle en temps
        do while (tnn < t1)
            if (schx == 1) then
                do i=1,imax 
                    Tnp1(i) = Tn(i) + D*dt*(Tn(i+1) - 2._PR*Tn(i) + Tn(i-1))/dx**2
                end do
            else if (schx == 2) then
                call chol_res(L_chol, Tn(1:imax), Tnp1(1:imax), check_chol)
            end if
            tnn = tnn + dt
            Tn(1:imax) = Tnp1(1:imax)
        end do
        
        write(ct,*) t1
       
        ! Écriture numérique
        open(unit=11, file="sol_"//prefix_sch//"_1_"//trim(adjustl(ct))//".dat")
        do i=0,imax+1
            write(11,*) x(i), Tn(i)
        end do
        close(11)

        ! Écriture exacte au temps réel d'arrêt (tnn au lieu de t1)
        open(unit=10, file="sol_exacte_1_"//trim(adjustl(ct))//".dat")
        do i=0,imax+1
            write(10,*) x(i), Texacte(tnn, x(i))
        end do
        close(10)
    end do

    ! Nettoyage de la partie 1
    deallocate(tab_t, Tn, Tnp1)
    if (schx == 2) deallocate(I_mat, M, L_chol)


    ! =========================================================================
    ! PARTIE 2 : Étude de convergence (Boucle sur les raffinements)
    ! =========================================================================
    print*, "-> Partie 2 : Calcul de la convergence (err.dat)..."
    
    allocate(tab_dx(nr), err_tab(3, nr))
    t1 = 0.1_PR ! Temps final cible
    
    do l = 1, nr
        
        call raffinement(l) ! Met à jour nx, dx et réalloue le domaine
        imax = nx
        allocate(Tn(0:imax+1), Tnp1(1:imax))
        
        ! Condition initiale
        do i=0,imax+1
            Tn(i) = Tinit(x(i))
        end do

        ! Préparation Implicite pour le nouveau maillage
        if (schx == 2) then
            allocate(I_mat(imax, imax), M(imax, imax), L_chol(imax, imax))
            I_mat = 0.0_PR
            do i = 1, imax
                I_mat(i,i) = 1.0_PR
            end do
            call fill_mat(A)
            M = I_mat + dt * A
            call chol_decomposition(M, L_chol, check_chol)
        end if

        tnn = 0.0_PR
        
        ! Boucle en temps jusqu'à t1 = 0.1
        do while (tnn < t1)
            if (schx == 1) then
                do i=1,imax 
                    Tnp1(i) = Tn(i) + D*dt*(Tn(i+1) - 2._PR*Tn(i) + Tn(i-1))/dx**2
                end do
            else if (schx == 2) then
                call chol_res(L_chol, Tn(1:imax), Tnp1(1:imax), check_chol)
            end if
            tnn = tnn + dt
            Tn(1:imax) = Tnp1(1:imax)
        end do

        ! Calcul et stockage de l'erreur avec tnn (temps réel d'arrêt)
        
        call Calcul_erreurs(Tn, tnn, e)
        tab_dx(l) = dx
        err_tab(:, l) = e

        print*, "   Raffinement", l, "| nx =", nx, "| Erreur L_inf =", e(1)

        ! Nettoyage itération
        deallocate(Tn, Tnp1)
        if (schx == 2) deallocate(I_mat, M, L_chol)

    end do

    ! Écriture du fichier final des erreurs
    call ecrit_erreurs(tab_dx, err_tab)
    print*, "=== Terminé ! Fichier err.dat généré. ==="

    deallocate(tab_dx, err_tab)
    call Fin()

end program chaleur