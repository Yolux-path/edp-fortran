module fonctions

    use donnees
    implicit none
    
    real(PR), parameter :: Pi = 4.0_PR*atan(1.0_PR), epsilon = 5*1e-3_PR
    integer             :: k = 1
    real(PR), parameter :: G_grav = 6.67430E-11_PR
    real(PR), parameter :: h_grav = 1000.0_PR
    
contains

    function f(z) result(y)
        real(PR), intent(in) :: z
        real(PR)             :: y
        
        if(choix_f == 1) then
            ! Cas test régulier
            y = 64*(Pi**2)*sin(8*Pi*z) + 90*(Pi**2)*cos(30*Pi*z)
        else if(choix_f == 2) then
            ! Cas f(x) = x^k
            y = z**k
        else if(choix_f == 3) then
            ! Cas f(x) = 1
            y = 1.0_PR
        else if(choix_f == 4) then
            y = 64*(Pi**2)*sin(8*Pi*z) + 90*(Pi**2)*cos(30*Pi*z)
        else if(choix_f == 5) then
            !Cas raide
            y= (1.0_PR/10) * (12*Pi)**2 * sin(12*Pi*z) - ((exp(z/epsilon))  / ((exp(1/epsilon)-1)*(epsilon**2)))
        end if
    end function f

    function sol_exacte(z) result(y)
        real(PR), intent(in) :: z
        real(PR)             :: y
        real(PR)             :: c1, c2
        
        if(choix_f == 1) then
            ! Solution exacte du cas régulier
            y = sin(8*Pi*z) + (1.0_PR/10.0_PR)*cos(30*Pi*z)

        else if(choix_f == 2) then
            ! -u'' = z^k => u'' = -z^k => u' = -z^(k+1)/(k+1) + C1
            ! u(z) = -z^(k+2)/((k+1)*(k+2)) + C1*z + C2
            c2 = ug
            c1 = (ud - ug + (Lx**(k+2))/((k+1.0_PR)*(k+2.0_PR))) / Lx
            y = -(z**(k+2))/((k+1.0_PR)*(k+2.0_PR)) + c1*z + c2
            
        else if (choix_f == 3) then
            ! -u'' = 1 => u'' = -1 => u' = -z + C1
            ! u(z) = -z^2/2 + C1*z + C2
            c2 = ug
            c1 = (ud - ug + (Lx**2)/2.0_PR) / Lx
            y = -(z**2)/2.0_PR + c1*z + c2

        else if (choix_f == 4) then
            ! Solution exacte du cas régulier
            y = sin(8*Pi*z) +(0.1_PR)*cos(30*Pi*z)
        else if (choix_f == 5) then
            ! Solution exacte du cas raide
            y = (1.0_PR/10) * sin(12*Pi*z) +( (exp(z/epsilon)-1) / (exp(1/epsilon)-1) ) 
        end if


    end function sol_exacte

    function f_gravite(z) result(y)
        real(PR), intent(in) :: z
        real(PR)             :: y, x_h
        
        x_h = abs(z) / h_grav
        
        if (x_h <= 1.0_PR) then
            y = (2.0_PR/3.0_PR) - x_h**2 + 0.5_PR * x_h**3
        else if (x_h > 1.0_PR .and. x_h <= 2.0_PR) then
            y = (1.0_PR/6.0_PR) * (2.0_PR - x_h)**3
        else
            y = 0.0_PR
        end if
        
        ! La source f(x) pour Poisson est 4 * Pi * G * rho
        ! On multiplie par le facteur donne dans l'enonce (2.4e4 / h)
        y = y * (2.4E4_PR / h_grav)
        y = 4.0_PR * Pi * G_grav * y
    end function f_gravite

end module fonctions