module fonctions

    use donnees
    
    implicit none
    
    real(PR), parameter         :: Pi = 4.0_PR*atan(1.0_PR)
    integer                     :: k  = 2
    
    contains

            function f(z) result(y) !Verifié avec les 3 choix
        
                real(PR), intent(in)    :: z
                real(PR)                :: y
                
                if(choix_f == 1) then
                    y = cos(30*Pi*z)

                else if(choix_f == 2) then
                    y = z**k

                else
                    y = 1
                end if
                
        
            end function f

            function fp(z) result(y)
        
                real(PR), intent(in)    :: z
                real(PR)                :: y
                
                if(choix_f == 1) then
                    y = -30*Pi*sin(30*Pi*z)

                else if(choix_f == 2) then
                    y = k*z**(k-1)

                else
                    y = 0
                end if
        
            end function fp
    
end module fonctions
