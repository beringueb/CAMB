module interfaces
    implicit none
    ! Declare interfaces for functions
    interface
        function Recombination_rayleigh_eff(a)
        real(dl), intent(in) :: a
        real(dl) zst,z,az,bz,Recombination_rayleigh_eff
        real(dl), parameter ::  zinitial = 1e4_dl !highest redshift
        real(dl), parameter ::  zfinal=0._dl
        integer,  parameter :: Nz=10000
        real(dl), parameter :: delta_z = (zinitial-zfinal)/Nz
        integer ilo,ihi
        z=1/a-1
        if (z.ge.zrec(1)) then
            Recombination_rayleigh_eff=x_rayleigh_eff(1)
        else
            if (z.le.zrec(nz)) then
                Recombination_rayleigh_eff=x_rayleigh_eff(nz)
            else
                zst=(zinitial-z)/delta_z
                ihi= int(zst)
                ilo = ihi+1
                az=zst - int(zst)
                bz=1-az     
                Recombination_rayleigh_eff=az*x_rayleigh_eff(ilo)+bz*x_rayleigh_eff(ihi)+ &
                ((az**3-az)*dx_rayleigh_eff(ilo)+(bz**3-bz)*dx_rayleigh_eff(ihi))/6._dl
            endif
        endif
        end function Recombination_rayleigh_eff    
    end interface
end module interfaces
