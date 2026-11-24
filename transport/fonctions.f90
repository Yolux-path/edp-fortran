module fonctions

    use donnees
    
    implicit none
    
    real(PR), parameter         :: Pi = 4.0_PR*atan(1.0_PR)
    integer                     :: k  = 3
    
    contains

        function Uinit(t,x) result(u_i)
            
            real(PR), intent(in) :: x,t
            real(PR)            :: u_i, m

            m = x-(a*t)

            if(choix_f == 0) then
                u_i = cos(2.0_PR*Pi*m)
            else if(choix_f == 1) then
                ! On centre à 0.3 et on divise par 0.01 pour l'affiner
                u_i = exp(-(m-0.7_PR)**2 / 0.01_PR)
            else if(choix_f == 2) then
                u_i = exp(-(m-5)**2) + exp(-20*(m-3)**2)
            else if(choix_f == 3) then
                if(m >= 0.4 .and. m <= 0.6) then
                    u_i = 1
                else
                    u_i = 0
                end if
            else if(choix_f == 4) then
                u_i = m**k
            end if

        end function Uinit

        function Ug(t) result(u_g)
            real(PR), intent(in) :: t
            real(PR) :: u_g
            if (t >= 2.0_PR .and. t <= 3.0_PR) then
                u_g = 30.0_PR
            else
                u_g = 0.0_PR
            end if
        end function Ug

        function Uexacte(t,x) result(u_e)
            real(PR), intent(in) :: x,t
            real(PR)             :: u_e, t_retard
            
            t_retard = t - (x / a)
            
            if (t_retard >= 2.0_PR .and. t_retard <= 3.0_PR) then
                u_e = 30.0_PR
            else
                u_e = 0.0_PR
            end if
        end function Uexacte

end module fonctions
