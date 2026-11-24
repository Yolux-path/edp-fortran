module fonctions
    use donnees 


    real(PR), parameter                 :: Pi = 4.0_PR*atan(1.0_PR)
    integer                             :: k  = 1
    

    contains 

    
    function Tinit(x) result(y)
        real(PR), intent(in)    :: x 
        real(PR)                :: y
        real(PR)                :: sum
        integer                 :: n

        if (choix_f == 1) then
            y = sin(2*k*Pi*x)
        else if (choix_f == 2) then
            sum = 0.25
            n = 1
            do while (n <=10)
                sum = sum + (((-1)**n - 1)/(Pi*n)**2)*cos(2*n*Pi*x)
                n = n + 1
            end do
            y = sum
        end if
    end function

    function Texacte(t,x) result(temp) 
        real(PR), intent(in)    :: t, x
        real(PR)                :: temp, sum
        integer                 :: n = 1
        
        if (choix_f == 1) then
            temp = exp(-D*(2*k*Pi)**2 * t) * sin(2*k*Pi*x)
        else if (choix_f == 2) then
            sum = 0.25
            n = 1
            do while (n <=100)
                sum = sum + exp(-D*(2*n*Pi)**2 * t) *(((-1)**n - 1)/(Pi*n)**2)*cos(2*n*Pi*x)
                n = n + 1
            end do
            temp = sum
        end if
    end function

end module fonctions