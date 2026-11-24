module mod_schemas
    use donnees
    implicit none
contains

    function amont(Ui, Uim1) result(Unp1)
        real(PR), intent(in)    :: Ui, Uim1
        real(PR)                :: Unp1
        Unp1 = Ui - c * (Ui - Uim1)
    end function amont

    function centre(Ui, Uim1, Uip1) result(Unp1)
        real(PR), intent(in)      :: Ui, Uim1, Uip1
        real(PR)                  :: Unp1
        Unp1 = Ui - (c / 2.0_PR) * (Uip1 - Uim1)
    end function centre

    function lax_wendroff(Ui, Uim1, Uip1) result(Unp1)
        real(PR), intent(in)    :: Ui, Uim1, Uip1
        real(PR)                :: Unp1
        Unp1 = Ui - (c / 2.0_PR) * (Uip1 - Uim1) + ((c**2) / 2.0_PR) * (Uip1 - 2.0_PR * Ui + Uim1)
    end function lax_wendroff

end module mod_schemas