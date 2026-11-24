module mod_algebre

    use donnees
    use fonctions

    implicit none
    
contains

    subroutine print_mat(tab, c)
        real(PR), dimension(:,:), intent(in)    :: tab
        character(len=*), intent(in)            ::  c
        
        integer :: i, j

        ! print*, size(tab(1,:))
        ! print*, size(tab(:,1))

        j = size(tab(:,1))

        print*, c
        do i = 1, j
                print*, tab(i,:)
        end do
    
        
    end subroutine print_mat
    
    subroutine m_to_lu(M, L, U)
        real(PR), dimension(:,:), intent(in)    :: M
        real(PR), dimension(:,:), intent(out)   :: L
        real(PR), dimension(:,:), intent(out)   :: U

        integer :: n, i, j

        n = size(M(1,:))

        ! print*, n

        do i = 1, n
            do j = 1 , i
                L(i,j) = M(i,j)
            end do
            L(i,i) = 1._PR
        end do

        ! call print_mat(L, "L=")

        do i = 1, n
            do j = i , n
                U(i,j) = M(i,j)
            end do
        end do

        ! call print_mat(U, "U=")

    end subroutine m_to_lu

    subroutine lu_check(A, M, check)
        real(PR), dimension(:,:), intent(in)    :: A
        real(PR), dimension(:,:), intent(in)    :: M
        logical, intent(out)                    :: check


        real(PR) :: tolerance
        integer :: n, i
        real(PR), dimension(:,:), allocatable :: L
        real(PR), dimension(:,:), allocatable :: U


        tolerance = 1e-3
        check = .true.
        n = size(M(1,:))

        allocate(L(n,n),U(n,n))

        call m_to_lu(M,L,U)

        ! call print_mat(L, "L")
        ! call print_mat( MATMUL(L,U),"LU =")
        ! print*, A(1,:) - L(1,:)
        ! print*,    shape(A(1,:) - L(1,:))


        L = MATMUL(L,U)

        do i = 1, n
            if ( MAXVAL( A(i,:) - L(i,:))  > tolerance ) then
                check = .false.
            end if
        end do

        deallocate(L,U)
        
    end subroutine lu_check

    subroutine lu_decomposition(A, M, check)
        real(PR), dimension(:,:), intent(in)    :: A
        real(PR), dimension(:,:), intent(out)   :: M
        logical, intent(out) :: check

        integer     :: n, i, j
        real(PR)    :: tolerance
        
        n = size(A(1,:))
        tolerance = 1e-3
        check = .true.

        if ( size(M(1,:)) /= n .or. size(M(:,1)) /= n) then 
            print*, "ERROR : Taille incompatible ! dans lu_decomposition "
            check = .false.
            return
        end if

        if( size(A(1,:)) /= size(A(:,1))) then
            print*, "ERROR : La matrice n'est pas carré dans lu_decomposition"
            check = .false.
            return
        end if

        ! print*, check

        M = A
        do i = 1, n-1
            do j = i+1, n  
                !check de la divison par 0 !
                if ( abs( M(i,i) ) < tolerance  ) then
                    check = .false.
                    print*, "ERROR : Division par 0 ! dans lu_decomposition"
                    return
                end if
                M (j,i) = M(j,i) / M(i, i)
                M (j,i+1:n)= M(j,i+1:n) - M(j,i) * M(i,i+1:n)
            end do 
        end do 

        ! call print_mat(M, "M=")
        
    end subroutine lu_decomposition

    subroutine lu_res(M, b, x, check)

        real(PR), dimension(:,:), intent(in)    :: M
        real(PR), dimension(:), intent(in)      :: b
        logical ,               intent(out)     :: check 
        real(PR), dimension(:), intent(out)     :: x

        real(PR), dimension(:,:), allocatable   :: L, U
        real(PR), dimension(:), allocatable     :: y
        real(PR)                                :: sum, tolerance
        integer                                 :: n, i, j
        

        n = size(b)
        tolerance = 1e-3
        check = .true.
        allocate(L(n,n), U(n,n), y(n))

        if(size(M(1,:)) /= n .or. size(M(1,:)) /= size(x)  ) then
            print*, "ERROR1 : Taille incompatible dans lu_res!! "
            check = .false.
            return
        end if

        y = 0.

        call m_to_lu(M, L, U)

        do i = 1, n
            sum = 0.
            do j = 1, i-1
                sum = sum + L(i,j) * y(j)
            end do
            ! if ( abs(L(i,j)) < tolerance ) then
            !     print*, "ERROR : Division par 0  dans lu_res!"
            !     check = .false.
            !     return
            ! end if
            y(i) = (b(i) - sum) / L(i,j)
        end do

        ! print*, "y =", y

        do i = n, 1, -1
            sum = 0.
            do j = i+1, n
                sum = sum + U(i,j) * x(j)
            end do
            if ( abs(U(i,i)) < tolerance ) then
                print*, "ERROR2 : Division par 0  dans lu_res!"
                check = .false.
                return
            end if
            x(i) = (y(i) - sum) / U(i,i)
        end do
        
        ! print*, "x = ", x

        deallocate(L, U, y)
    end subroutine lu_res

    subroutine lu_test(A, b, x, check)
        real(PR), dimension(:,:), intent(in)    :: A
        real(PR), dimension(:), intent(in)      :: b
        real(PR), dimension(:), intent(out)     :: x
        logical ,               intent(out)     :: check 

        real(PR)                                :: tolerance
        integer                                 :: n

        n = size(b)
        tolerance = 1e-3
        check = .true.

        !check de la taille
        if(size(A(1,:)) /= n .or. size(A(1,:)) /= size(x)  ) then
            print*, "ERROR : Taille incompatible dans lu_test!! "
            check = .false.
            return
        end if

        if( size(A(1,:)) /= size(A(:,1))) then
            print*, "ERROR : La matrice n'est pas carré dans lu_test"
            check = .false.
            return
        end if

        !Ax = b ?

        if ( MAXVAL( MATMUL(A,x) - b)  > tolerance ) then
            check = .false.
            print*, "ax /= b "
        end if


    end subroutine lu_test

    subroutine chol_check(A, L, check)
        real(PR), dimension(:,:), intent(in)    :: A
        real(PR), dimension(:,:), intent(in)    :: L
        logical, intent(out)                    :: check

        integer                                 :: n, i
        real(PR), dimension(:,:), allocatable   :: Lt, LLt
        real(PR)                                :: tolerance

        n = size(A(1,:))
        tolerance = 1e-3
        check = .true.

        allocate(Lt(n,n), LLt(n,n))

        if ( size(L(1,:)) /= n .or. size(L(:,1)) /= n) then 
            print*, "ERROR : Taille incompatible ! dans chol "
            check = .false.
            return
        end if

        if( size(A(1,:)) /= size(A(:,1))) then
            print*, "ERROR : La matrice n'est pas carré dans chol"
            check = .false.
            return
        end if

        
        Lt = TRANSPOSE(L)
        
        LLt = MATMUL(L, Lt)

        ! call print_mat(LLt, "L*Lt = ")

        do i = 1, n
            if ( MAXVAL( A(i,:) - LLt(i,:))  > tolerance ) then
                check = .false.
            end if
        end do

        deallocate(Lt,LLt)

    end subroutine chol_check

    subroutine chol_decomposition(A, L, check)
        real(PR), dimension(:,:), intent(in)    :: A
        real(PR), dimension(:,:), intent(out)   :: L
        logical, intent(out)                    :: check

        integer                                 :: n, i, j, k
        real(PR)                                :: sum, tolerance
        
        n = size(A(1,:))
        tolerance = 1e-3
        check = .true.

        !check de la taille
        if ( size(L(1,:)) /= n .or. size(L(:,1)) /= n) then 
            print*, "ERROR : Taille incompatible ! dans chol "
            check = .false.
            return
        end if

        !check carré
        if( size(A(1,:)) /= size(A(:,1))) then
            print*, "ERROR : La matrice n'est pas carré dans chol"
            check = .false.
            return
        end if

        if ( any( A /= TRANSPOSE(A) ) ) then
            print*, "ERROR : La matrice n'est pas symétrique dans chol_decomposition"
            check = .false.
            return
        end if

        !Cholesky decomposition
        L = 0.
        do i = 1, n
            sum = 0.
            do k = 1, i-1
                sum = sum + L(i,k) ** 2
            end do

            if( A(i,i) - sum < 0 ) then
                print*, "ERROR : La matrice n'est pas définie positive dans chol_decomposition"
                check = .false.
                return
            end if

            L(i,i) = sqrt( A(i,i) - sum )

            do j = i+1, n 
                sum = 0.
                do k = 1, i-1
                    sum = sum + L(i,k) * L(j,k)
                end do
                if ( abs(L(i,i)) < tolerance ) then
                    print*, "ERROR : Division par 0  dans chol_decomposition!"
                    check = .false.
                    return
                end if
                L(j,i) = (A(j,i) - sum ) / L(i,i)
            end do
        end do

        
    end subroutine chol_decomposition

    subroutine chol_res(L, b, x, check)
        real(PR), dimension(:,:), intent(in)        :: L
        real(PR), dimension(:), intent(in)          :: b
        real(PR), dimension(:), intent(inout)       :: x
        logical ,               intent(inout)       :: check

        real(PR), dimension(:), allocatable     :: y
        real(PR)                                :: sum, tolerance
        integer                                 :: n, i, j

        n = size(b)
        tolerance = 1e-3

        allocate(y(n))

        if(size(L(1,:)) /= n .or. size(L(1,:)) /= size(x)  ) then
            print*, "ERROR : Taille incompatible dans chol_res!! "
            check = .false.
            return
        end if

        y = 0.

        do i = 1, n
            sum = 0.
            do j = 1, i-1
                sum = sum + L(i,j) * y(j)
            end do
            if ( abs(L(i,j)) < tolerance ) then
                print*, "ERROR : Division par 0  dans chol_res!"
                check = .false.
                return
            end if
            y(i) = (b(i) - sum) / L(i,j)
        end do 

        do i = n, 1, -1
            sum = 0.
            do j = i+1, n
                sum = sum + L(j,i) * x(j)
            end do
            if ( abs(L(i,i)) < tolerance ) then
                print*, "ERROR : Division par 0  dans chol_res!"
                check = .false.
                return
            end if
            x(i) = (y(i) - sum) / L(i,i)
        end do


        deallocate(y)
    end subroutine chol_res

    subroutine chol_test(A, b, x, check)
        real(PR), dimension(:,:), intent(in)    :: A
        real(PR), dimension(:), intent(in)      :: b
        real(PR), dimension(:), intent(out)     :: x
        logical ,               intent(out)     :: check 

        real(PR)                                :: tolerance
        integer                                 :: n

        n = size(b)
        tolerance = 1e-3
        check = .true.

        !check de la taille
        if(size(A(1,:)) /= n .or. size(A(1,:)) /= size(x)  ) then
            print*, "ERROR : Taille incompatible dans lu_test!! "
            check = .false.
            return
        end if

        if( size(A(1,:)) /= size(A(:,1))) then
            print*, "ERROR : La matrice n'est pas carré dans lu_test"
            check = .false.
            return
        end if

        !Ax = b ?

        if ( MAXVAL(abs( MATMUL(A,x) - b))  > tolerance ) then
            check = .false.
            print*, "ax /= b "
        end if


    end subroutine chol_test

    subroutine fill_mat(A)
        real(PR), dimension(:,:), intent(inout)   :: A
        integer                                   :: i, j 

        if (size(A(1,:)) /= size(A(:,1))) then
            print*, "ERROR : La matrice n'est pas carré dans fill_mat"
            return
        else
            print*, "Remplissage de la matrice A pour le schéma centré 3pts"
            do i = 1, size(A(1,:))
                do j = 1, size(A(:,1))
                    if (i == j) then
                        A(i,j) = -2._PR
                    else if ( abs(i-j) == 1 ) then
                        A(i,j) = 1._PR
                    else 
                        A(i,j) = 0._PR
                    end if
                end do
            end do
            A = A*(-D/dx**2)
        end if

    end subroutine fill_mat


end module mod_algebre