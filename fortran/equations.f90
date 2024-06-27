    ! Equations module for background and ! To avoid circular module issues, some things are not part of module

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    ! Background evolution, return d tau/ d a, where tau is the conformal time
    function dtauda(this,a)
    !use results
    !and
    use Recombination
    !and
    use DarkEnergyInterface
    implicit none
    class(CAMBdata) :: this
    real(dl), intent(in) :: a
    real(dl) :: dtauda, grhoa2, grhov_t
    integer :: f_i = 1
    !integer, parameter :: num_cmb_freq = 6 !!!
    !integer, parameter :: nscatter = num_cmb_freq+1

    call this%CP%DarkEnergy%BackgroundDensityAndPressure(this%grhov, a, grhov_t)

    !  8*pi*G*rho*a**4.
    grhoa2 = this%grho_no_de(a) +  grhov_t * a**2
    if (grhoa2 <= 0) then
        call GlobalError('Universe stops expanding before today (recollapse not supported)', error_unsupported_params)
        dtauda = 0
    else
        dtauda = sqrt(3 / grhoa2)
    end if

    end function dtauda

    !Gauge-dependent perturbation equations

    module GaugeInterface
    use precision
    !use results
    !and
    use Recombination
    !and
    use MassiveNu
    use DarkEnergyInterface
    use Transfer
    implicit none
    public

    ! Equations and relation to synchronous gauge variables documented in the notes:
    ! https://cosmologist.info/notes/CAMB.pdf

    !Description of this file. Change if you make modifications.
    ! andrea
    character(LEN=*), parameter :: Eqns_name = 'cdm_rayleigh_backeffect'
    logical, parameter :: rayleigh_back = .false.
    ! andrea
    logical, parameter :: plot_evolve = .false. !for outputing time evolution

    integer, parameter :: basic_num_eqns = 4
    integer, parameter :: ix_etak=1, ix_clxc=2, ix_clxb=3, ix_vb=4 !Scalar array indices for each quantity
    integer, parameter :: ixt_H = 1, ixt_shear = 2 !tensor indices

    logical :: DoTensorNeutrinos = .true.

    logical, parameter :: second_order_tightcoupling = .true.

    real(dl) :: Magnetic = 0._dl
    !Vector mode anisotropic stress in units of rho_gamma
    real(dl) :: vec_sig0 = 1._dl
    !Vector mode shear
    integer, parameter :: max_l_evolve = 256 !Maximum l we are ever likely to propagate
    !Note higher values increase size of Evolution vars, hence memory

    !Supported scalar initial condition flags
    integer, parameter :: initial_adiabatic=1, initial_iso_CDM=2, &
        initial_iso_baryon=3,  initial_iso_neutrino=4, initial_iso_neutrino_vel=5, initial_vector = 0
    integer, parameter :: initial_nummodes =  initial_iso_neutrino_vel

    type EvolutionVars

        ! andrea
        real(dl) RayleighSwitchOnTime
        logical Rayleigh
        integer g_ix_freq, polind_freq, freq_neq
        ! andrea
        real(dl) q, q2
        real(dl) k_buf,k2_buf ! set in initial
        logical :: is_cosmological_constant

        integer w_ix !Index of two quintessence equations
        integer Tg_ix !index of matter temerature perturbation
        integer reion_line_ix !index of matter temerature perturbation

        integer xe_ix !index of x_e perturbation
        integer Ts_ix !index of Delta_{T_s}

        integer r_ix !Index of the massless neutrino hierarchy
        integer g_ix !Index of the photon neutrino hierarchy

        integer q_ix !index into q_evolve array that gives the value q
        logical TransferOnly

        !       nvar  - number of scalar (tensor) equations for this k
        integer nvar,nvart, nvarv

        !Max_l for the various hierarchies
        integer lmaxg,lmaxnr,lmaxnu,lmaxgpol,MaxlNeeded
        integer lmaxnrt, lmaxnut, lmaxt, lmaxpolt, MaxlNeededt
        logical EvolveTensorMassiveNu(max_nu)
        integer lmaxnrv, lmaxv, lmaxpolv
        integer lmaxline !21cm multipoles for getting reionization effect

        integer polind  !index into scalar array of polarization hierarchy

        !array indices for massive neutrino equations
        integer nu_ix(max_nu), nu_pert_ix
        integer nq(max_nu), lmaxnu_pert
        logical has_nu_relativistic

        !Initial values for massive neutrino v*3 variables calculated when switching
        !to non-relativistic approx
        real(dl) G11(max_nu),G30(max_nu)
        !True when using non-relativistic approximation
        logical MassiveNuApprox(max_nu)
        real(dl) MassiveNuApproxTime(max_nu)

        !True when truncating at l=2,3 when k*tau>>1 (see arXiv:1104.2933)
        logical high_ktau_neutrino_approx

        !Massive neutrino scheme being used at the moment
        integer NuMethod

        !True when using tight-coupling approximation (required for stability at early times)
        logical TightCoupling, TensTightCoupling
        real(dl) TightSwitchoffTime

        !Numer of scalar equations we are propagating
        integer ScalEqsToPropagate
        integer TensEqsToPropagate
        !beta > l for closed models
        integer FirstZerolForBeta
        !Tensor vars
        real(dl) aux_buf

        real(dl) pig, pigdot
        real(dl) poltruncfac

        logical no_nu_multpoles, no_phot_multpoles
        integer lmaxnu_tau(max_nu)  !lmax for massive neutinos at time being integrated
        logical nu_nonrelativistic(max_nu)

        real(dl) denlk(max_l_evolve),denlk2(max_l_evolve), polfack(max_l_evolve)
        real(dl) Kf(max_l_evolve)

        integer E_ix, B_ix !tensor polarizatisdon indices
        real(dl) denlkt(4,max_l_evolve),Kft(max_l_evolve)

        logical :: saha !still high x_e
        logical :: evolve_TM !\delta T_g evolved separately
        logical :: evolve_baryon_cs !evolving sound speed rather than using background approx

        !Workaround for ifort, gives class pointer to avoid creating temps and huge slow down
        class(TThermoData), pointer :: ThermoData => null()

        real, pointer :: OutputTransfer(:) => null()
        real(dl), pointer :: OutputSources(:) => null()
        real(dl), pointer :: CustomSources(:) => null()
        integer :: OutputStep = 0

    end type EvolutionVars

    ABSTRACT INTERFACE
    SUBROUTINE TSource_func(sources, tau, a, adotoa, grho, gpres,w_lam, cs2_lam,  &
        grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
        k,etak, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc,clxr, clxnu, clxde, delta_p_b, &
        dgq, qg, qr, qde, vb, qgdot, qrdot, vbdot, &
        dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E, Edot, &
        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
        tau0, tau_maxvis, Kf, f_K)
    use precision
    integer, parameter :: num_cmb_freq = 6 !!!
    integer, parameter :: nscatter = num_cmb_freq+1
    real(dl), intent(out) :: sources(:)
    real(dl), intent(in) :: tau, a, adotoa, grho, gpres,w_lam, cs2_lam,  &
                            grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
                            k,etak, etakdot, phi, phidot, sigma, sigmadot, &
                            dgrho, clxg,clxb,clxc, clxr, clxnu, clxde, delta_p_b, &
                            dgq, qg, qr, qde, vb, qgdot, qrdot, vbdot, &
                            dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
                            polter, polterdot, polterddot, octg, octgdot, E(2:3), Edot(2:3), &
                            tau0, tau_maxvis
    real(dl), intent(in) :: opacity(nscatter), dopacity(nscatter), ddopacity(nscatter), exptau(nscatter), &
                            visibility(nscatter), dvisibility(nscatter), ddvisibility(nscatter)
    REAL(dl), intent(in) :: Kf(*)
    real(dl), external :: f_K
    END SUBROUTINE TSource_func
    END INTERFACE

    class(CAMBdata), pointer :: State !Current state.
    !(Due to ifort bug State needs to be class pointer to avoid making temp class pointers when calling functions)
    type(CAMBParams), pointer :: CP   !Current parameters (part of state)

    !precalculated arrays
    real(dl) polfac(max_l_evolve),denl(max_l_evolve),vecfac(max_l_evolve),vecfacpol(max_l_evolve)

    real(dl), parameter :: ep0=1.0d-2
    integer, parameter :: lmaxnu_high_ktau=4 !Jan2015, increased from 3 to fix mpk for mnu~6eV

    real(dl) epsw
    real(dl), allocatable :: nu_tau_notmassless(:,:)
    real(dl) nu_tau_nonrelativistic(max_nu), nu_tau_massive(max_nu)

    procedure(state_function), private :: dtauda
    contains

    ! andrea
    subroutine Phot_Integrate_L012(EV,y,drhonu,fnu,pinu,opacity, opac_qgphot, opac_phot)
        type(EvolutionVars) EV
!  Compute the perturbations of density and energy flux
!  of photons in units of the mean density for Blackbody
        real(dl), intent(in) :: y(EV%nvar)
        real(dl), intent(OUT) ::  drhonu,fnu, pinu
        real(dl), optional :: opacity(nscatter), opac_qgphot, opac_phot
        integer iq, ind

!  q is the comoving momentum in units of k_B*T_gamma0/c.

        drhonu=0
        fnu=0
        pinu=0
        if (present(opacity)) then
        opac_qgphot=0
        opac_phot=0
        end if
        do iq=1,num_cmb_freq
            ind = EV%g_ix_freq + (iq-1)*EV%freq_neq
            drhonu=drhonu+ phot_int_kernel(iq)* y(ind)
            fnu=fnu+phot_int_kernel(iq)* y(ind+1)
            pinu=pinu+ phot_int_kernel(iq)*y(ind+2)
            if (present(opacity)) then
            opac_qgphot = opac_qgphot + phot_int_kernel(iq)* y(ind+1) * opacity(iq+1)
            opac_phot = opac_phot + phot_int_kernel(iq)*opacity(iq+1)
            end if
        end do

    end subroutine  Phot_Integrate_L012
    ! andrea

    subroutine SetActiveState(P)
    class(CAMBdata), target :: P

    State => P
    CP => P%CP

    end subroutine SetActiveState


    subroutine GaugeInterface_ScalEv(EV, this, etat, y, tau, tauend, tol1, ind, c1,w)
    type(EvolutionVars) EV
    class(TThermoData) :: this
    class(TRecfast) :: etat
    real(dl) c1(24),w(EV%nvar,9), y(EV%nvar), tol1, tau, tauend
    integer ind
    call dverk_derivs(EV, this, etat, EV%ScalEqsToPropagate, tau,y,tauend,tol1,ind,c1,EV%nvar,w)
    if (ind==-3) then
        call GlobalError('Dverk error -3: the subroutine was unable  to  satisfy  the  error ' &
            //'requirement  with a particular step-size that is less than or * ' &
            //'equal to hmin, which may mean that tol is too small' &
            //'--- but most likely you''ve messed up the y array indexing; ' &
            //'compiling with bounds checking may (or may not) help find the problem.',error_evolution)
    end if
    end subroutine GaugeInterface_ScalEv

    function f_K(x)
    real(dl) :: f_K
    real(dl), intent(in) :: x

    f_K = State%curvature_radius*State%rofChi(x/State%curvature_radius)

    end function f_K

    function next_nu_nq(nq) result (next_nq)
    integer, intent(in) :: nq
    integer q, next_nq

    if (nq==0) then
        next_nq=1
    else
        q = int(State%NuPerturbations%nu_q(nq))
        if (q>=10) then
            next_nq = State%NuPerturbations%nqmax
        else
            next_nq = nq+1
        end if
    end if

    end function next_nu_nq

    recursive subroutine GaugeInterface_EvolveScal(EV, this, etat, tau,y,tauend,tol1,ind,c1,w)
    use Recombination, only : CB1
    type(EvolutionVars) EV, EVout
    class(TThermoData) :: this
    class(TRecfast) :: etat
    real(dl) c1(24),w(EV%nvar,9), y(EV%nvar), yout(EV%nvar), tol1, tau, tauend
    integer ind, nu_i
    real(dl) cs2, opacity(nscatter), dopacity(nscatter)
    real(dl) tau_switch_ktau, tau_switch_nu_massless, tau_switch_nu_massive, next_switch
    real(dl) tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles,tau_switch_nu_nonrel
    real(dl) noSwitch, smallTime
    !Sources
    real(dl) tau_switch_saha, Delta_TM, xe,a,tau_switch_evolve_TM

    noSwitch= State%tau0+1
    smallTime =  min(tau, 1/EV%k_buf)/100

    tau_switch_ktau = noSwitch
    tau_switch_no_nu_multpoles= noSwitch
    tau_switch_no_phot_multpoles= noSwitch

    !Massive neutrino switches
    tau_switch_nu_massless = noSwitch
    tau_switch_nu_nonrel = noSwitch
    tau_switch_nu_massive= noSwitch

    !Sources
    tau_switch_saha=noSwitch
    ! and
    if (CP%Evolve_delta_xe .and. EV%saha)  tau_switch_saha = this%recombination_saha_tau
    tau_switch_evolve_TM=noSwitch
    if (EV%Evolve_baryon_cs .and. .not. EV%Evolve_tm) tau_switch_evolve_TM = this%recombination_Tgas_tau
    ! and

    !Evolve equations from tau to tauend, performing switches in equations if necessary.

    if (.not. EV%high_ktau_neutrino_approx .and. .not. EV%no_nu_multpoles ) then
        tau_switch_ktau=  max(20, EV%lmaxnr-4)/EV%k_buf
    end if

    if (CP%Num_Nu_massive /= 0) then
        do nu_i = 1, CP%Nu_mass_eigenstates
            if (EV%nq(nu_i) /= State%NuPerturbations%nqmax) then
                tau_switch_nu_massless = min(tau_switch_nu_massless,nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i))
            else if (.not. EV%nu_nonrelativistic(nu_i)) then
                tau_switch_nu_nonrel = min(nu_tau_nonrelativistic(nu_i),tau_switch_nu_nonrel)
            else if (EV%NuMethod==Nu_trunc .and..not. EV%MassiveNuApprox(nu_i)) then
                tau_switch_nu_massive = min(tau_switch_nu_massive,EV%MassiveNuApproxTime(nu_i))
            end if
        end do
    end if
    if (CP%DoLateRadTruncation) then
        if (.not. EV%no_nu_multpoles) & !!.and. .not. EV%has_nu_relativistic .and. tau_switch_nu_massless ==noSwitch)  &
            tau_switch_no_nu_multpoles= &
            max(15/EV%k_buf*CP%Accuracy%AccuracyBoost,min(State%taurend,this%matter_verydom_tau))

        if (.not. EV%no_phot_multpoles .and. (.not.CP%WantCls .or. EV%k_buf>0.03*CP%Accuracy%AccuracyBoost)) &
            tau_switch_no_phot_multpoles =max(15/EV%k_buf,State%taurend)*CP%Accuracy%AccuracyBoost
    end if
    ! andrea
    next_switch = min(tau_switch_ktau, tau_switch_nu_massless, EV%TightSwitchoffTime, tau_switch_nu_massive, &
                  tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles, tau_switch_nu_nonrel, &
                  EV%RayleighSwitchOnTime, noSwitch, tau_switch_saha, tau_switch_evolve_TM)
    !andrea
    if (next_switch < tauend) then
        if (next_switch > tau+smallTime) then
            call GaugeInterface_ScalEv(EV, this, etat, y, tau,next_switch,tol1,ind,c1,w)
            if (global_error_flag/=0) return
        end if
        EVout=EV
        if (next_switch == EV%TightSwitchoffTime) then
            !TightCoupling
            EVout%TightCoupling=.false.
            EVout%TightSwitchoffTime = noSwitch
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            EV=EVout
            y=yout
            ind=1
            !Set up variables with their tight coupling values
            y(EV%g_ix+2) = EV%pig
            call this%values(tau,a, cs2,opacity,dopacity)
            if (second_order_tightcoupling) then
                ! Francis-Yan Cyr-Racine November 2010

                y(EV%g_ix+3) = (3._dl/7._dl)*y(EV%g_ix+2)*(EV%k_buf/opacity(1))*(1._dl+dopacity(1)/opacity(1)**2) + &
                    (3._dl/7._dl)*EV%pigdot*(EV%k_buf/opacity(1)**2)*(-1._dl)

                y(EV%polind+2) = EV%pig/4 + EV%pigdot*(1._dl/opacity(1))*(-5._dl/8._dl- &
                    (25._dl/16._dl)*dopacity(1)/opacity(1)**2) + &
                    EV%pig*(EV%k_buf/opacity(1))**2*(-5._dl/56._dl)
                y(EV%polind+3) = (3._dl/7._dl)*(EV%k_buf/opacity(1))*y(EV%polind+2)*(1._dl + &
                    dopacity(1)/opacity(1)**2) + (3._dl/7._dl)*(EV%k_buf/opacity(1)**2)*((EV%pigdot/4._dl)* &
                    (1._dl+(5._dl/2._dl)*dopacity(1)/opacity(1)**2))*(-1._dl)
            else
                y(EV%g_ix+3) = 3./7*y(EV%g_ix+2)*EV%k_buf/opacity(1)
                y(EV%polind+2) = EV%pig/4
                y(EV%polind+3) =y(EV%g_ix+3)/4
            end if
        ! andrea
        else if (next_switch == EV%RayleighSwitchOnTime) then
            EVout%RayleighSwitchOnTime = noSwitch
            EVout%Rayleigh = .false.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            EV=EVout
            y=yout
            ind=1
            ! andrea
            !k tau >> 1, evolve massless neutrino effective fluid up to l=2
        else if (next_switch==tau_switch_ktau) then
            EVout%high_ktau_neutrino_approx=.true.
            EVout%nq(1:CP%Nu_mass_eigenstates) = State%NuPerturbations%nqmax
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch == tau_switch_nu_massless) then
            !Mass starts to become important, start evolving next momentum mode
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (EV%nq(nu_i) /= State%NuPerturbations%nqmax) then
                    if (next_switch == nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i)) then
                        EVOut%nq(nu_i) = next_nu_nq(EV%nq(nu_i))
                        call SetupScalarArrayIndices(EVout)
                        call CopyScalarVariableArray(y,yout, EV, EVout)
                        EV=EVout
                        y=yout
                        exit
                    endif
                end if
            end do
        else if (next_switch == tau_switch_nu_nonrel) then
            !Neutrino becomes non-relativistic, don't need high L
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (.not. EV%nu_nonrelativistic(nu_i) .and.  next_switch==nu_tau_nonrelativistic(nu_i) ) then
                    EVout%nu_nonrelativistic(nu_i)=.true.
                    call SetupScalarArrayIndices(EVout)
                    call CopyScalarVariableArray(y,yout, EV, EVout)
                    EV=EVout
                    y=yout
                    exit
                end if
            end do
        else if (next_switch == tau_switch_nu_massive) then
            !Very non-relativistic neutrinos, switch to truncated velocity-weight hierarchy
            ! and modification DU CAMB INITIAL HORS RAYLEIGH
            call this%values(tau,a,cs2,opacity,dopacity)
            ! and
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (.not. EV%MassiveNuApprox(nu_i) .and.  next_switch== EV%MassiveNuApproxTime(nu_i) ) then
                    call SwitchToMassiveNuApprox(EV, a, y, nu_i)
                    exit
                end if
            end do
        else if (next_switch==tau_switch_no_nu_multpoles) then
            !Turn off neutrino hierarchies at late time where slow and not needed.
            ind=1
            EVout%no_nu_multpoles=.true.
            EVOut%nq(1:CP%Nu_mass_eigenstates ) = State%NuPerturbations%nqmax
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch==tau_switch_no_phot_multpoles) then
            !Turn off photon hierarchies at late time where slow and not needed.
            ind=1
            EVout%no_phot_multpoles=.true.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch==tau_switch_saha) then
            !Sources
            ind=1
            EVout%saha = .false.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            !and modification DU CAMB INITIAL HORS RAYLEIGH
            call this%Values(tau,a, cs2,opacity, dopacity)
            ! and
            y=yout
            EV=EVout
            Delta_Tm = y(EV%g_ix)/4 ! assume delta_TM = delta_T_gamma
            xe= CP%Recomb%x_e(a)
            y(EV%xe_ix) = (1-xe)/(2-xe)*(-y(ix_clxb) + (3./2+  CB1/(CP%TCMB/a))*Delta_TM)
        else if (next_switch==tau_switch_evolve_TM) then
            !Sources
            ind=1
            EVout%evolve_TM = .true.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
            y(EV%Tg_ix) =y(EV%g_ix)/4 ! assume delta_TM = delta_T_gamma
        end if
        call GaugeInterface_EvolveScal(EV,this, etat, tau,y,tauend,tol1,ind,c1,w)
        return
    end if

    call GaugeInterface_ScalEv(EV,this,etat,y,tau,tauend,tol1,ind,c1,w)

    end subroutine GaugeInterface_EvolveScal

    subroutine GaugeInterface_EvolveTens(EV, this, etat, tau,y,tauend,tol1,ind,c1,w)
    use recombination
    type(EvolutionVars) EV, EVOut
    class(TThermoData) :: this
    class(TRecfast) :: etat
    real(dl) c1(24),w(EV%nvart,9), y(EV%nvart),yout(EV%nvart), tol1, tau, tauend
    integer ind
    real(dl) opacity(nscatter), cs2, a

    if (EV%TensTightCoupling .and. tauend > EV%TightSwitchoffTime) then
        if (EV%TightSwitchoffTime > tau) then
            call dverk_derivst(EV, this, etat, EV%TensEqsToPropagate,tau,y,EV%TightSwitchoffTime,tol1,ind,c1,EV%nvart,w)
        end if
        EVOut=EV
        EVOut%TensTightCoupling = .false.
        call SetupTensorArrayIndices(EVout)
        call CopyTensorVariableArray(y,yout,Ev, Evout)
        Ev = EvOut
        y=yout
        ! and modification DU CAMB INITIAL HORS RAYLEIGH
        call this%values(tau,a,cs2,opacity)
        ! and
        y(EV%g_ix+2)= 32._dl/45._dl*EV%k_buf/opacity(1)*y(ixt_shear)
        y(EV%E_ix+2) = y(EV%g_ix+2)/4
        ! andrea
        EVout%Rayleigh = .false.
        call SetupTensorArrayIndices(EVout)
        call CopyTensorVariableArray(y,yout, EV, EVout)
        EV=EVout
        y=yout
        ind=1
        ! andrea
    end if

    call dverk_derivst(EV, this, etat, EV%TensEqsToPropagate,tau,y,tauend,tol1,ind,c1,EV%nvart,w)

    end subroutine GaugeInterface_EvolveTens

    function DeltaTimeMaxed(a1,a2, tol) result(t)
    real(dl) a1,a2,t
    real(dl), optional :: tol
    if (a1>1._dl) then
        t=0
    elseif (a2 > 1._dl) then
        t = State%DeltaTime(a1,1.01_dl, tol)
    else
        t = State%DeltaTime(a1,a2, tol)
    end if
    end function DeltaTimeMaxed

    subroutine GaugeInterface_Init
    !Precompute various arrays and other things independent of wavenumber
    integer j, nu_i
    real(dl) a_nonrel, a_mass,a_massive, time, nu_mass

    epsw = 100/State%tau0

    if (CP%WantScalars) then
        do j=2,max_l_evolve
            polfac(j)=real((j+3)*(j-1),dl)/(j+1)
        end do
    end if

    if (CP%WantVectors) then
        do j=2,max_l_evolve
            vecfac(j)=real((j+2),dl)/(j+1)
            vecfacpol(j)=real((j+3)*j,dl)*(j-1)*vecfac(j)/(j+1)**2
        end do
    end if

    do j=1,max_l_evolve
        denl(j)=1._dl/(2*j+1)
    end do

    if (CP%Nu_Mass_eigenstates>0) then
        associate(nqmax => State%NuPerturbations%nqmax, nu_q => State%NuPerturbations%nu_q)
            if (allocated(nu_tau_notmassless)) deallocate(nu_tau_notmassless)
            allocate(nu_tau_notmassless(nqmax,max_nu))
            do nu_i=1, CP%Nu_Mass_eigenstates
                nu_mass = max(0.1_dl,State%nu_masses(nu_i))
                a_mass =  1.e-1_dl/nu_mass/CP%Accuracy%lAccuracyBoost
                time=State%DeltaTime(0._dl,State%NuPerturbations%nu_q(1)*a_mass)
                nu_tau_notmassless(1, nu_i) = time
                do j=2,nqmax
                    !times when each momentum mode becomes signficantly nonrelativistic
                    time= time + DeltaTimeMaxed(nu_q(j-1)*a_mass,nu_q(j)*a_mass, 0.01_dl)
                    nu_tau_notmassless(j, nu_i) = time
                end do

                a_nonrel =  2.5d0/nu_mass*CP%Accuracy%AccuracyBoost
                nu_tau_nonrelativistic(nu_i) =DeltaTimeMaxed(0._dl,a_nonrel)
                a_massive =  17.d0/nu_mass*CP%Accuracy%AccuracyBoost
                nu_tau_massive(nu_i) =nu_tau_nonrelativistic(nu_i) + DeltaTimeMaxed(a_nonrel,a_massive)
            end do
        end associate
    end if

    end subroutine GaugeInterface_Init


    subroutine SetupScalarArrayIndices(EV, max_num_eqns)
    !Set up array indices after the lmax have been decided
    use MassiveNu
    !andrea
    use Recombination
    ! andrea
    !Set the numer of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    integer, intent(out), optional :: max_num_eqns
    integer neq, maxeq, nu_i

    neq=basic_num_eqns
    maxeq=neq
    if (.not. EV%no_phot_multpoles) then
        !Photon multipoles
        EV%g_ix=basic_num_eqns+1
        if (EV%TightCoupling) then
            neq=neq+2
        else
            neq = neq+ (EV%lmaxg+1)
            !Polarization multipoles
            EV%polind = neq -1 !polind+2 is L=2, for polarizationthe first calculated
            neq=neq + EV%lmaxgpol-1
        end if
        ! andrea
        if (EV%Rayleigh) then
            EV%g_ix_freq=neq+1
            EV%polind_freq = neq+ (EV%lmaxg+1) -1
            EV%freq_neq=(EV%lmaxg+1+EV%lmaxgpol-1)
            neq=neq+EV%freq_neq*num_cmb_freq
        end if
        ! andrea
    end if
    ! andrea
    maxeq=maxeq+(EV%lmaxg+1+EV%lmaxgpol-1)*num_cmb_freq
    ! andrea
    if (.not. EV%no_nu_multpoles) then
        !Massless neutrino multipoles
        EV%r_ix= neq+1
        if (EV%high_ktau_neutrino_approx) then
            neq=neq + 3
        else
            neq=neq + (EV%lmaxnr+1)
        end if
    end if
    maxeq = maxeq +  (EV%lmaxg+1)+(EV%lmaxnr+1)+EV%lmaxgpol-1

    !Dark energy
    if (.not. CP%DarkEnergy%is_cosmological_constant) then
        EV%w_ix = neq + 1
        neq = neq + CP%DarkEnergy%num_perturb_equations
        maxeq = maxeq + CP%DarkEnergy%num_perturb_equations
    else
        EV%w_ix = 0
    end if

    !Sources
    if (CP%Evolve_delta_xe) then
        if (.not. EV%saha) then
            EV%xe_ix = neq+1
            neq=neq+1
        end if
        maxeq=maxeq+1
    end if

    if (EV%Evolve_baryon_cs) then
        if (EV%Evolve_TM) then
            EV%Tg_ix = neq+1
            neq=neq+1
        end if
        maxeq=maxeq+1
        if (CP%Do21cm .and. CP%SourceTerms%line_reionization) then
            EV%reion_line_ix = neq+1
            neq=neq+ EV%lmaxline+1 +  EV%lmaxline-1
            maxeq=maxeq+EV%lmaxline+1 +  EV%lmaxline-1
        end if
    end if

    if (CP%Evolve_delta_Ts) then
        EV%Ts_ix = neq+1
        neq=neq+1
        maxeq=maxeq+1
    end if

    !Massive neutrinos
    if (CP%Num_Nu_massive /= 0) then
        EV%has_nu_relativistic = any(EV%nq(1:CP%Nu_Mass_eigenstates)/=State%NuPerturbations%nqmax)
        if (EV%has_nu_relativistic) then
            EV%lmaxnu_pert=EV%lmaxnu
            EV%nu_pert_ix=neq+1
            neq = neq+ EV%lmaxnu_pert+1
            maxeq=maxeq+ EV%lmaxnu_pert+1
        else
            EV%lmaxnu_pert=0
        end if

        do nu_i=1, CP%Nu_Mass_eigenstates
            if (EV%high_ktau_neutrino_approx) then
                EV%lmaxnu_tau(nu_i) = int(lmaxnu_high_ktau *CP%Accuracy%lAccuracyBoost)
                if (CP%Transfer%accurate_massive_neutrinos) EV%lmaxnu_tau(nu_i) = EV%lmaxnu_tau(nu_i) *3
            else
                EV%lmaxnu_tau(nu_i) =max(min(nint(0.8_dl*EV%q*nu_tau_nonrelativistic(nu_i) &
                    *CP%Accuracy%lAccuracyBoost),EV%lmaxnu),3)
                !!!Feb13tweak
                if (EV%nu_nonrelativistic(nu_i)) EV%lmaxnu_tau(nu_i)= &
                    min(EV%lmaxnu_tau(nu_i),nint(4*CP%Accuracy%lAccuracyBoost))
            end if
            if (State%nu_masses(nu_i) > 5000 .and. CP%Transfer%high_precision) &
                EV%lmaxnu_tau(nu_i) = EV%lmaxnu_tau(nu_i)*2 !megadamping
            EV%lmaxnu_tau(nu_i)=min(EV%lmaxnu,EV%lmaxnu_tau(nu_i))

            EV%nu_ix(nu_i)=neq+1
            if (EV%MassiveNuApprox(nu_i)) then
                neq = neq+4
            else
                neq = neq+ EV%nq(nu_i)*(EV%lmaxnu_tau(nu_i)+1)
            endif
            maxeq = maxeq + State%NuPerturbations%nqmax*(EV%lmaxnu+1)
        end do
    else
        EV%has_nu_relativistic = .false.
    end if

    EV%ScalEqsToPropagate = neq
    if (present(max_num_eqns)) then
        max_num_eqns=maxeq
    end if

    end subroutine SetupScalarArrayIndices

    subroutine CopyScalarVariableArray(y,yout, EV, EVout)
    ! andrea
    use Recombination
    ! andrea
    type(EvolutionVars) EV, EVOut
    real(dl), intent(in) :: y(EV%nvar)
    real(dl), intent(out) :: yout(EVout%nvar)
    integer lmax,i, nq
    integer nnueq,nu_i, ix_off, ix_off2, ind, ind2
    real(dl) q, pert_scale

    yout=0
    yout(1:basic_num_eqns) = y(1:basic_num_eqns)

    ! DarkEnergy
    if (CP%DarkEnergy%num_perturb_equations > 0) &
        yout(EVOut%w_ix:EVOut%w_ix + CP%DarkEnergy%num_perturb_equations - 1) = &
        y(EV%w_ix:EV%w_ix + CP%DarkEnergy%num_perturb_equations - 1)

    if (.not. EV%no_phot_multpoles .and. .not. EVout%no_phot_multpoles) then
        if (EV%TightCoupling .or. EVOut%TightCoupling) then
            lmax=1
        else
            lmax = min(EV%lmaxg,EVout%lmaxg)
        end if
        yout(EVout%g_ix:EVout%g_ix+lmax)=y(EV%g_ix:EV%g_ix+lmax)
        if (.not. EV%TightCoupling .and. .not. EVOut%TightCoupling) then
            lmax = min(EV%lmaxgpol,EVout%lmaxgpol)
            yout(EVout%polind+2:EVout%polind+lmax)=y(EV%polind+2:EV%polind+lmax)
        end if
        ! andrea
        if (EV%Rayleigh .and. EVout%Rayleigh) then
             do i=1,num_cmb_freq
                 !assume number of frequencies is fixed
                 ix_off2 = EVOut%g_ix_freq + (i-1)*EVOut%freq_neq
                 ix_off = EV%g_ix_freq + (i-1)*EV%freq_neq
                 lmax = min(EV%lmaxg,EVout%lmaxg)
                 yout(ix_off2:ix_off2+lmax)=y(ix_off:ix_off+lmax)
                 lmax = min(EV%lmaxgpol,EVout%lmaxgpol)
                 ix_off2 = EVOut%polind_freq + (i-1)*EVOut%freq_neq
                 ix_off = EV%polind_freq + (i-1)*EV%freq_neq
                 yout(ix_off2+2:ix_off2+lmax)=y(ix_off+2:ix_off+lmax)
             end do
        end if
        ! andrea
    end if

    if (.not. EV%no_nu_multpoles .and. .not. EVout%no_nu_multpoles) then
        if (EV%high_ktau_neutrino_approx .or. EVout%high_ktau_neutrino_approx) then
            lmax=2
        else
            lmax = min(EV%lmaxnr,EVout%lmaxnr)
        end if
        yout(EVout%r_ix:EVout%r_ix+lmax)=y(EV%r_ix:EV%r_ix+lmax)
    end if

    if (CP%Num_Nu_massive /= 0) then
        do nu_i=1,CP%Nu_mass_eigenstates
            ix_off=EV%nu_ix(nu_i)
            ix_off2=EVOut%nu_ix(nu_i)
            if (EV%MassiveNuApprox(nu_i) .and. EVout%MassiveNuApprox(nu_i)) then
                nnueq=4
                yout(ix_off2:ix_off2+nnueq-1)=y(ix_off:ix_off+nnueq-1)
            else if (.not. EV%MassiveNuApprox(nu_i) .and. .not. EVout%MassiveNuApprox(nu_i)) then
                lmax=min(EV%lmaxnu_tau(nu_i),EVOut%lmaxnu_tau(nu_i))
                nq = min(EV%nq(nu_i), EVOut%nq(nu_i))
                do i=1,nq
                    ind= ix_off + (i-1)*(EV%lmaxnu_tau(nu_i)+1)
                    ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                    yout(ind2:ind2+lmax) = y(ind:ind+lmax)
                end do
                do i=nq+1, EVOut%nq(nu_i)
                    lmax = min(EVOut%lmaxnu_tau(nu_i), EV%lmaxnr)
                    ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                    yout(ind2:ind2+lmax) = y(EV%r_ix:EV%r_ix+lmax)

                    !Add leading correction for the mass
                    q=State%NuPerturbations%nu_q(i)
                    pert_scale=(State%nu_masses(nu_i)/q)**2/2
                    lmax = min(lmax,EV%lmaxnu_pert)
                    yout(ind2:ind2+lmax) = yout(ind2:ind2+lmax) &
                        + y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)*pert_scale
                end do
            end if
        end do

        if (EVOut%has_nu_relativistic .and. EV%has_nu_relativistic) then
            lmax = min(EVOut%lmaxnu_pert, EV%lmaxnu_pert)
            yout(EVout%nu_pert_ix:EVout%nu_pert_ix+lmax)=  y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)
        end if
    end if
    !Sources
    if (.not. EV%saha .and. .not. EVOut%saha) then
        yout(EVOut%xe_ix) =y(EV%xe_ix)
    end if
    if (EV%Evolve_baryon_cs) then
        if (EV%Evolve_TM .and. EVout%Evolve_TM) yout(EVOut%Tg_ix) = y(EV%Tg_ix)
        if (CP%Do21cm .and. CP%SourceTerms%line_reionization) then
            yout(EVOut%reion_line_ix:EVOut%reion_line_ix+EVout%lmaxline +  EVout%lmaxline-1) = &
                y(EV%reion_line_ix:EV%reion_line_ix+EV%lmaxline +  EV%lmaxline-1)
        end if
    end if
    if (CP%Evolve_delta_Ts) then
        yout(EVOut%Ts_ix) = y(EV%Ts_ix)
    end if

    end subroutine CopyScalarVariableArray


    subroutine SetupTensorArrayIndices(EV, maxeq)
    type(EvolutionVars) EV
    integer nu_i, neq
    integer, optional, intent(out) :: maxeq
    neq=2
    EV%g_ix = neq-1 !EV%g_ix+2 is quadrupole
    if (.not. EV%TensTightCoupling) then
        EV%E_ix = EV%g_ix + (EV%lmaxt-1)
        EV%B_ix = EV%E_ix + (EV%lmaxpolt-1)
        neq = neq+ (EV%lmaxt-1)+(EV%lmaxpolt-1)*2
        ! andrea
        if (EV%Rayleigh) then
               EV%g_ix_freq=neq-1
               EV%polind_freq = neq+ (EV%lmaxt-1) -1
               EV%freq_neq=EV%lmaxt-1+(EV%lmaxpolt-1)*2
               neq=neq+EV%freq_neq*num_cmb_freq
        end if
    end if
    if (present(maxeq)) then
        maxeq = 2 + ((EV%lmaxt-1)+(EV%lmaxpolt-1)*2) * (num_cmb_freq+1)
        ! andrea
    end if
    EV%r_ix = neq -1
    if (DoTensorNeutrinos) then
        neq = neq + EV%lmaxnrt-1
        if (present(maxeq)) maxeq = maxeq+EV%lmaxnrt-1
        if (CP%Num_Nu_massive /= 0 ) then
            do nu_i=1, CP%nu_mass_eigenstates
                EV%EvolveTensorMassiveNu(nu_i) = &
                    nu_tau_nonrelativistic(nu_i) < 0.8*State%tau_maxvis*CP%Accuracy%AccuracyBoost
                if (EV%EvolveTensorMassiveNu(nu_i)) then
                    EV%nu_ix(nu_i)=neq-1
                    neq = neq+ State%NuPerturbations%nqmax*(EV%lmaxnut-1)
                    if (present(maxeq)) maxeq = maxeq + State%NuPerturbations%nqmax*(EV%lmaxnut-1)
                end if
            end do
        end if
    end if

    EV%TensEqsToPropagate = neq

    end  subroutine SetupTensorArrayIndices

    subroutine CopyTensorVariableArray(y,yout, EV, EVout)
    type(EvolutionVars) EV, EVOut
    real(dl), intent(in) :: y(EV%nvart)
    real(dl), intent(out) :: yout(EVout%nvart)
    integer lmaxpolt, lmaxt, nu_i, ind, ind2, i
    ! andrea
    integer ix_off, ix_off2, lmax
    ! andrea
    yout=0
    yout(1:2) = y(1:2)
    if (.not. EVOut%TensTightCoupling .and. .not.EV%TensTightCoupling) then
        lmaxt = min(EVOut%lmaxt,EV%lmaxt)
        yout(EVout%g_ix+2:EVout%g_ix+lmaxt)=y(EV%g_ix+2:EV%g_ix+lmaxt)
        lmaxpolt = min(EV%lmaxpolt, EVOut%lmaxpolt)
        yout(EVout%E_ix+2:EVout%E_ix+lmaxpolt)=y(EV%E_ix+2:EV%E_ix+lmaxpolt)
        yout(EVout%B_ix+2:EVout%B_ix+lmaxpolt)=y(EV%B_ix+2:EV%B_ix+lmaxpolt)

        ! andrea
        if (EV%Rayleigh .and. EVout%Rayleigh) then
             do i=1,num_cmb_freq
             ix_off2 = EVOut%g_ix_freq + (i-1)*EVOut%freq_neq
             ix_off = EV%g_ix_freq + (i-1)*EV%freq_neq
             lmax = min(EV%lmaxt,EVout%lmaxt)
             yout(ix_off2+2:ix_off2+lmax)=y(ix_off+2:ix_off+lmax)
             lmax = min(EV%lmaxgpol,EVout%lmaxgpol)
             ix_off2 = EVOut%polind_freq + (i-1)*EVOut%freq_neq
             ix_off = EV%polind_freq + (i-1)*EV%freq_neq
             yout(ix_off2+2:ix_off2+lmax)=y(ix_off+2:ix_off+lmax)
             yout(ix_off2+2+(lmax-1):ix_off2+(lmax-1)+lmax)=y(ix_off+2+(lmax-1):ix_off+(lmax-1)+lmax)
             end do
        end if
        ! andrea

    end if
    if (DoTensorNeutrinos) then
        lmaxt=min(EV%lmaxnrt,EVOut%lmaxnrt)
        yout(EVout%r_ix+2:EVout%r_ix+lmaxt)=y(EV%r_ix+2:EV%r_ix+lmaxt)
        do nu_i =1, CP%nu_mass_eigenstates
            if (EV%EvolveTensorMassiveNu(nu_i)) then
                lmaxt=min(EV%lmaxnut,EVOut%lmaxnut)
                do i=1,State%NuPerturbations%nqmax
                    ind= EV%nu_ix(nu_i) + (i-1)*(EV%lmaxnut-1)
                    ind2=EVOut%nu_ix(nu_i)+ (i-1)*(EVOut%lmaxnut-1)
                    yout(ind2+2:ind2+lmaxt) = y(ind+2:ind+lmaxt)
                end do
            end if
        end do
    end if

    end subroutine CopyTensorVariableArray

    subroutine GetNumEqns(EV)
    use MassiveNu
    !Set the numer of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    real(dl) scal, max_nu_mass
    integer nu_i,q_rel,j

    if (CP%Num_Nu_massive == 0) then
        EV%lmaxnu=0
        max_nu_mass=0
    else
        max_nu_mass = maxval(State%nu_masses(1:CP%Nu_mass_eigenstates))
        do nu_i = 1, CP%Nu_mass_eigenstates
            !Start with momentum modes for which t_k ~ time at which mode becomes non-relativistic
            q_rel=0
            do j=1, State%NuPerturbations%nqmax
                !two different q's here EV%q ~k
                if (State%NuPerturbations%nu_q(j) > State%nu_masses(nu_i)*State%adotrad/EV%q) exit
                q_rel = q_rel + 1
            end do

            if (q_rel>= State%NuPerturbations%nqmax-2 .or. CP%WantTensors) then
                EV%nq(nu_i)=State%NuPerturbations%nqmax
            else
                EV%nq(nu_i)=q_rel
            end if
            !q_rel = nint(nu_masses(nu_i)*adotrad/EV%q) !two dffierent q's here EV%q ~k
            !EV%nq(nu_i)=max(0,min(nqmax0,q_rel)) !number of momentum modes to evolve intitially
            EV%nu_nonrelativistic(nu_i) = .false.
        end do

        EV%NuMethod = CP%MassiveNuMethod
        if (EV%NuMethod == Nu_Best) EV%NuMethod = Nu_Trunc
        !l_max for massive neutrinos
        if (CP%Transfer%high_precision) then
            EV%lmaxnu=nint(25*CP%Accuracy%lAccuracyBoost)
        else
            EV%lmaxnu=max(3,nint(10*CP%Accuracy%lAccuracyBoost))
            if (max_nu_mass>700) EV%lmaxnu=max(3,nint(15*CP%Accuracy%lAccuracyBoost)) !Feb13 tweak
        endif
    end if

    if (State%closed) then
        EV%FirstZerolForBeta = nint(EV%q*State%curvature_radius)
    else
        EV%FirstZerolForBeta= 100000 !a large number
    end if

    EV%high_ktau_neutrino_approx = .false.
    if (CP%WantScalars) then
        EV%TightCoupling=.true.
        EV%no_phot_multpoles =.false.
        EV%no_nu_multpoles =.false.
        EV%MassiveNuApprox=.false.
        !Sources
        EV%saha = .true.
        EV%Evolve_TM = .false.

        if (CP%Accuracy%AccuratePolarization) then
            EV%lmaxg  = max(nint(11*CP%Accuracy%lAccuracyBoost),3)
        else
            EV%lmaxg  = max(nint(8*CP%Accuracy%lAccuracyBoost),3)
        end if
        EV%lmaxnr = max(nint(14*CP%Accuracy%lAccuracyBoost),3)
        if (max_nu_mass>700) EV%lmaxnr = max(nint(32*CP%Accuracy%lAccuracyBoost),3) !Feb13 tweak

        EV%lmaxgpol = EV%lmaxg
        if (.not.CP%Accuracy%AccuratePolarization) EV%lmaxgpol=max(nint(4*CP%Accuracy%lAccuracyBoost),3)

        if (EV%q < 0.05) then
            !Large scales need fewer equations
            scal  = 1
            if (CP%Accuracy%AccuratePolarization) scal = 4  !But need more to get polarization right
            EV%lmaxgpol=max(3,nint(min(8,nint(scal* 150* EV%q))*CP%Accuracy%lAccuracyBoost))
            EV%lmaxnr=max(3,nint(min(7,nint(sqrt(scal)* 150 * EV%q))*CP%Accuracy%lAccuracyBoost))
            if (EV%lmaxnr < EV%lmaxnu) then
                ! Nov 2020 change following Pavel Motloch report
                EV%lmaxnr = EV%lmaxnu
                !EV%lmaxnu = min(EV%lmaxnu, EV%lmaxnr) ! may be better but have not tested and makes small result changes
            endif
            EV%lmaxg=max(3,nint(min(8,nint(sqrt(scal) *300 * EV%q))*CP%Accuracy%lAccuracyBoost))
            !Sources
            if (CP%SourceTerms%line_phot_quadrupole) then
                EV%lmaxg=EV%lmaxg*8
                EV%lmaxgpol=EV%lmaxgpol*4
            elseif (CP%Accuracy%AccurateReionization) then
                EV%lmaxg=EV%lmaxg*4
                EV%lmaxgpol=EV%lmaxgpol*2
            end if
        end if

        if (EV%TransferOnly) then
            EV%lmaxgpol = min(EV%lmaxgpol,nint(5*CP%Accuracy%lAccuracyBoost))
            EV%lmaxg = min(EV%lmaxg,nint(6*CP%Accuracy%lAccuracyBoost))
        end if
        if (CP%Transfer%high_precision .or. CP%Do21cm) then
            EV%lmaxnr=max(nint(45*CP%Accuracy%lAccuracyBoost),3)
            if (EV%q > 0.04 .and. EV%q < 0.5) then !baryon oscillation scales
                EV%lmaxg=max(EV%lmaxg,10)
            end if
        end if

        if (CP%Do21cm .and. CP%SourceTerms%line_reionization) then
            EV%lmaxg =  EV%lmaxg*8
            EV%lmaxgpol = EV%lmaxgpol*3
        end if

        EV%Evolve_baryon_cs = CP%Do21cm .or.CP%Evolve_delta_xe .or. CP%Evolve_delta_Ts

        if (CP%Do21cm .and. CP%SourceTerms%line_reionization) then
            EV%lmaxline  = EV%lmaxg
        end if

        if (State%closed) then
            EV%lmaxnu=min(EV%lmaxnu, EV%FirstZerolForBeta-1)
            EV%lmaxnr=min(EV%lmaxnr, EV%FirstZerolForBeta-1)
            EV%lmaxg=min(EV%lmaxg, EV%FirstZerolForBeta-1)
            EV%lmaxgpol=min(EV%lmaxgpol, EV%FirstZerolForBeta-1)
        end if

        EV%poltruncfac=real(EV%lmaxgpol,dl)/max(1,(EV%lmaxgpol-2))
        EV%MaxlNeeded=max(EV%lmaxg,EV%lmaxnr,EV%lmaxgpol,EV%lmaxnu)
        if (EV%MaxlNeeded > max_l_evolve) call MpiStop('Need to increase max_l_evolve')
        call SetupScalarArrayIndices(EV,EV%nvar)
        if (State%closed) EV%nvar=EV%nvar+1 !so can reference lmax+1 with zero coefficient
        EV%lmaxt=0
    else
        EV%nvar=0
    end if

    if (CP%WantTensors) then
        EV%TensTightCoupling = .true.
        EV%lmaxt=max(3,nint(8*CP%Accuracy%lAccuracyBoost))
        EV%lmaxpolt = max(3,nint(4*CP%Accuracy%lAccuracyBoost))
        ! if (EV%q < 1e-3) EV%lmaxpolt=EV%lmaxpolt+1
        if (DoTensorNeutrinos) then
            EV%lmaxnrt=nint(6*CP%Accuracy%lAccuracyBoost)
            EV%lmaxnut=EV%lmaxnrt
        else
            EV%lmaxnut=0
            EV%lmaxnrt=0
        end if
        if (State%closed) then
            EV%lmaxt=min(EV%FirstZerolForBeta-1,EV%lmaxt)
            EV%lmaxpolt=min(EV%FirstZerolForBeta-1,EV%lmaxpolt)
            EV%lmaxnrt=min(EV%FirstZerolForBeta-1,EV%lmaxnrt)
            EV%lmaxnut=min(EV%FirstZerolForBeta-1,EV%lmaxnut)
        end if
        EV%MaxlNeededt=max(EV%lmaxpolt,EV%lmaxt, EV%lmaxnrt, EV%lmaxnut)
        if (EV%MaxlNeededt > max_l_evolve) call MpiStop('Need to increase max_l_evolve')
        call SetupTensorArrayIndices(EV, EV%nvart)
    else
        EV%nvart=0
    end if


    if (CP%WantVectors) then
        EV%lmaxv=max(10,nint(8*CP%Accuracy%lAccuracyBoost))
        EV%lmaxpolv = max(5,nint(5*CP%Accuracy%lAccuracyBoost))

        EV%nvarv=(EV%lmaxv)+(EV%lmaxpolv-1)*2+3

        EV%lmaxnrv=nint(30*CP%Accuracy%lAccuracyBoost)

        EV%nvarv=EV%nvarv+EV%lmaxnrv
        if (CP%Num_Nu_massive /= 0 ) then
            call MpiStop('massive neutrinos not supported for vector modes')
        end if
    else
        EV%nvarv=0
    end if

    end subroutine GetNumEqns

    subroutine SwitchToMassiveNuApprox(EV,a, y, nu_i)
    !When the neutrinos are no longer highly relativistic we use a truncated
    !energy-integrated hierarchy going up to third order in velocity dispersion
    type(EvolutionVars) EV, EVout
    integer, intent(in) :: nu_i
    real(dl) a,a2,pnu,clxnu,dpnu,pinu,rhonu
    real(dl) qnu
    real(dl) y(EV%nvar), yout(EV%nvar)

    a2=a*a
    EVout=EV
    EVout%MassiveNuApprox(nu_i)=.true.
    call SetupScalarArrayIndices(EVout)
    call CopyScalarVariableArray(y,yout, EV, EVout)

    !Get density and pressure as ratio to massles by interpolation from table
    call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),rhonu,pnu)

    !Integrate over q
    call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
    !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
    dpnu=dpnu/rhonu
    qnu=qnu/rhonu
    clxnu = clxnu/rhonu
    pinu=pinu/rhonu

    yout(EVout%nu_ix(nu_i))=clxnu
    yout(EVout%nu_ix(nu_i)+1)=dpnu
    yout(EVout%nu_ix(nu_i)+2)=qnu
    yout(EVout%nu_ix(nu_i)+3)=pinu

    call Nu_Intvsq(EV,y, a, nu_i, EVout%G11(nu_i),EVout%G30(nu_i))
    !Analytic solution for higher moments, proportional to a^{-3}
    EVout%G11(nu_i)=EVout%G11(nu_i)*a2*a/rhonu
    EVout%G30(nu_i)=EVout%G30(nu_i)*a2*a/rhonu

    EV=EVout
    y=yout

    end subroutine SwitchToMassiveNuApprox

    subroutine MassiveNuVarsOut(EV,y,yprime,a,adotoa,grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum,clxnu_all)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), yprime(EV%nvar),a, adotoa
    real(dl), optional :: grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum,clxnu_all
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    !dgpi = a^2 kappa pi (anisotropic stress)
    !dgpi_diff = a^2 kappa (3*p -rho)*pi

    integer nu_i
    real(dl) pinudot,grhormass_t, rhonu, pnu,  rhonudot
    real(dl) grhonu_t,gpnu_t
    real(dl) clxnu, qnu, pinu, dpnu, grhonu, dgrhonu

    grhonu=0
    dgrhonu=0
    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass_t=State%grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),rhonu,pnu)

        if (EV%MassiveNuApprox(nu_i)) then
            clxnu=y(EV%nu_ix(nu_i))
            !dpnu = y(EV%iq0+1+off_ix)
            qnu=y(EV%nu_ix(nu_i)+2)
            pinu=y(EV%nu_ix(nu_i)+3)
            pinudot=yprime(EV%nu_ix(nu_i)+3)
        else
            !Integrate over q
            call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            !dpnu=dpnu/rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
            pinu=pinu/rhonu
            rhonudot = ThermalNuBack%drho(a*State%nu_masses(nu_i),adotoa)

            call Nu_pinudot(EV,y, yprime, a,adotoa, nu_i,pinudot)
            pinudot=pinudot/rhonu - rhonudot/rhonu*pinu
        endif

        grhonu_t=grhormass_t*rhonu
        gpnu_t=grhormass_t*pnu

        grhonu = grhonu  + grhonu_t
        if (present(gpres)) gpres= gpres + gpnu_t

        dgrhonu= dgrhonu + grhonu_t*clxnu
        if (present(dgq)) dgq  = dgq   + grhonu_t*qnu
        if (present(dgpi)) dgpi = dgpi  + grhonu_t*pinu
        if (present(dgpi_diff)) dgpi_diff = dgpi_diff + pinu*(3*gpnu_t-grhonu_t)
        if (present(pidot_sum)) pidot_sum = pidot_sum + grhonu_t*pinudot
    end do
    if (present(grho)) grho = grho  + grhonu
    if (present(dgrho)) dgrho= dgrho + dgrhonu
    if (present(clxnu_all)) clxnu_all = dgrhonu/grhonu

    end subroutine MassiveNuVarsOut

    subroutine Nu_Integrate_L012(EV,y,a,nu_i,drhonu,fnu,dpnu,pinu)
    type(EvolutionVars) EV
    !  Compute the perturbations of density and energy flux
    !  of one eigenstate of massive neutrinos, in units of the mean
    !  density of one eigenstate of massless neutrinos, by integrating over
    !  momentum.
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  drhonu,fnu
    real(dl), optional, intent(OUT) :: dpnu,pinu
    real(dl) tmp, am, aq,v, pert_scale
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.

    drhonu=0
    fnu=0
    if (present(dpnu)) then
        dpnu=0
        pinu=0
    end if
    am=a*State%nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    associate(nu_q=>State%NuPerturbations%nu_q, nu_int_kernel=>State%NuPerturbations%nu_int_kernel)
        do iq=1,EV%nq(nu_i)
            aq=am/nu_q(iq)
            v=1._dl/sqrt(1._dl+aq*aq)
            drhonu=drhonu+ nu_int_kernel(iq)* y(ind)/v
            fnu=fnu+nu_int_kernel(iq)* y(ind+1)
            if (present(dpnu)) then
                dpnu=dpnu+  nu_int_kernel(iq)* y(ind)*v
                pinu=pinu+ nu_int_kernel(iq)*y(ind+2)*v
            end if
            ind=ind+EV%lmaxnu_tau(nu_i)+1
        end do
        ind = EV%nu_pert_ix
        do iq=EV%nq(nu_i)+1,State%NuPerturbations%nqmax
            !Get the rest from perturbatively relativistic expansion
            aq=am/nu_q(iq)
            v=1._dl/sqrt(1._dl+aq*aq)
            pert_scale=(State%nu_masses(nu_i)/nu_q(iq))**2/2
            tmp = nu_int_kernel(iq)*(y(EV%r_ix)  + pert_scale*y(ind))
            drhonu=drhonu+ tmp/v
            fnu=fnu+nu_int_kernel(iq)*(y(EV%r_ix+1)+ pert_scale*y(ind+1))
            if (present(dpnu)) then
                dpnu=dpnu+ tmp*v
                pinu = pinu+ nu_int_kernel(iq)*(y(EV%r_ix+2)+ pert_scale*y(ind+2))*v
            end if
        end do
    end associate
    if (present(dpnu)) then
        dpnu = dpnu/3
    end if

    end subroutine Nu_Integrate_L012

    subroutine Nu_pinudot(EV,y, ydot, a,adotoa, nu_i,pinudot)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a,adotoa, y(EV%nvar), ydot(EV%nvar)

    !  Compute the time derivative of the mean density in massive neutrinos
    !  and the shear perturbation.
    real(dl) pinudot
    real(dl) aq,q,v,aqdot,vdot
    real(dl) psi2,psi2dot
    real(dl) am, pert_scale
    integer iq,ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    pinudot=0._dl
    ind=EV%nu_ix(nu_i)+2
    am=a*State%nu_masses(nu_i)
    do iq=1,EV%nq(nu_i)
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        pinudot=pinudot+State%NuPerturbations%nu_int_kernel(iq)*(ydot(ind)*v+y(ind)*vdot)
        ind=ind+EV%lmaxnu_tau(nu_i)+1
    end do
    ind = EV%nu_pert_ix+2
    do iq=EV%nq(nu_i)+1,State%NuPerturbations%nqmax
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        pert_scale=(State%nu_masses(nu_i)/q)**2/2
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        psi2dot=ydot(EV%r_ix+2)  + pert_scale*ydot(ind)
        psi2=y(EV%r_ix+2)  + pert_scale*y(ind)
        pinudot=pinudot+State%NuPerturbations%nu_int_kernel(iq)*(psi2dot*v+psi2*vdot)
    end do

    end subroutine Nu_pinudot

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function Nu_pi(EV, y, a, nu_i) result(pinu)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvart)
    real(dl) :: am
    real(dl) pinu,q,aq,v
    integer iq, ind

    if (EV%nq(nu_i)/=State%NuPerturbations%nqmax) call MpiStop('Nu_pi: nq/=nqmax')
    pinu=0
    ind=EV%nu_ix(nu_i)+2
    am=a*State%nu_masses(nu_i)
    do iq=1, EV%nq(nu_i)
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        pinu=pinu+State%NuPerturbations%nu_int_kernel(iq)*y(ind)*v
        ind =ind+EV%lmaxnut+1
    end do

    end function Nu_pi

    !cccccccccccccccccccccccccccccccccccccccccccccc
    subroutine Nu_Intvsq(EV,y, a, nu_i, G11,G30)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  G11,G30

    !  Compute the third order variables (in velocity dispersion)
    !by integrating over momentum.
    real(dl) aq,q,v, am
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    am=a*State%nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    G11=0._dl
    G30=0._dl
    if (EV%nq(nu_i)/=State%NuPerturbations%nqmax) call MpiStop('Nu_Intvsq nq/=nqmax0')
    do iq=1, EV%nq(nu_i)
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        G11=G11+State%NuPerturbations%nu_int_kernel(iq)*y(ind+1)*v**2
        if (EV%lmaxnu_tau(nu_i)>2) then
            G30=G30+State%NuPerturbations%nu_int_kernel(iq)*y(ind+3)*v**2
        end if
        ind = ind+EV%lmaxnu_tau(nu_i)+1
    end do

    end subroutine Nu_Intvsq


    subroutine MassiveNuVars(EV,y,a,grho,gpres,dgrho,dgq, wnu_arr)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), a, grho,gpres,dgrho,dgq
    real(dl), intent(out), optional :: wnu_arr(max_nu)
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    integer nu_i
    real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu

    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass_t=State%grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),rhonu,pnu)

        if (EV%MassiveNuApprox(nu_i)) then
            clxnu=y(EV%nu_ix(nu_i))
            qnu=y(EV%nu_ix(nu_i)+2)
        else
            !Integrate over q
            call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu)
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
        endif

        grhonu_t=grhormass_t*rhonu
        gpnu_t=grhormass_t*pnu

        grho = grho  + grhonu_t
        gpres= gpres + gpnu_t
        dgrho= dgrho + grhonu_t*clxnu
        dgq  = dgq   + grhonu_t*qnu

        if (present(wnu_arr)) then
            wnu_arr(nu_i) =pnu/rhonu
        end if
    end do

    end subroutine MassiveNuVars


    function Get21cm_source2(a,Delta_source,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe, vterm )
    !Delta_Tspin - Delta_TCMB
    use constants
    use DarkAge21cm
    !vterm = hdot/clh + k*n/3/clh
    real(dl), intent(in) :: a,Delta_source,Delta_TCMB,Delta_Tm,Tmat,Trad,xe, Delta_xe, vterm
    real(dl) :: Get21cm_source2
    real(dl) Rgamma,Rm
    real(dl) dC10, n_H,C10, C10_HH, C10_eH
    real(dl) kappa_HH,kappa_eH
    real(dl) tau_eps
    real(dl) H
    real(dl) TSpin
    n_H = State%NNow/a**3
    kappa_HH = kappa_HH_21cm(Tmat, .false.)
    kappa_eH = kappa_eH_21cm(Tmat, .false.)
    C10_HH = n_H*kappa_HH* (1- xe)
    C10_eH = n_H*kappa_eH*xe
    C10 = C10_HH + C10_eH    !only relevant when He ionization is negligible
    Rgamma = 1._dl/(C10+A10*Trad/T_21cm)
    Rm = 1._dl/(C10+A10*Tmat/T_21cm)

    !          TSpin=TsRecfast(a)
    !          write(*,'(9e15.5)') 1/a-1,Tmat,Tspin, Trad,C10_HH,C10_eH,A10*Trad/T_21cm,xe,&
    !                  n_H*kappa_pH_21cm(Tmat, .false.)*xe
    !          if (a>0.5) stop

    dC10 = (C10*Delta_source + (C10_HH*kappa_HH_21cm(Tmat, .true.) &
        + C10_eH*kappa_eH_21cm(Tmat, .true.)) * &
        Delta_Tm + (kappa_eH-kappa_HH)*xe*n_H*Delta_xe)

    Get21cm_source2 =  dC10*(Rgamma-Rm) +  C10*(Rm*Delta_tm - Delta_TCMB*Rgamma)

    TSpin=CP%Recomb%T_s(a)
    H = (1/(a*dtauda(State,a)))
    tau_eps = a*line21_const*State%NNow/a**3/H/Tspin/1000

    Get21cm_source2 = Get21cm_source2 + &
        tau_eps/2*A10*( 1/(C10*T_21cm/Tmat+A10) -  1/(C10*T_21cm/Trad+A10) ) * &
        (Delta_source -vterm + dC10/C10 + 2*( - Rgamma*dC10 + Delta_TCMB*(C10*Rgamma-1)) &
        + Trad/(Tmat-Trad)*(Delta_tm-Delta_TCMB)   )

    end function Get21cm_source2


    function Get21cm_dTs(a,Delta_n,Delta_Ts,Delta_TCMB,Delta_Tm,Tmat,Trad,xe )
    !d Delta T_s / d eta dropping small \Delta_xe terms
    use constants
    use DarkAge21cm
    real(dl), intent(in) :: a,Delta_n,Delta_Ts,Delta_TCMB,Delta_Tm,Tmat,Trad,xe
    real(dl) :: Get21cm_dTs
    real(dl) n_H,C10, C10_HH, C10_eH, delta_C10
    real(dl) kappa_HH,kappa_eH, TSpin

    n_H = State%NNow/a**3
    kappa_HH = kappa_HH_21cm(Tmat, .false.)
    kappa_eH = kappa_eH_21cm(Tmat, .false.)
    C10_HH = n_H*kappa_HH* (1- xe)
    C10_eH = n_H*kappa_eH*xe
    C10 = C10_HH + C10_eH    !only relevant when He ionization is negligible
    TSpin=CP%Recomb%T_s(a)
    delta_C10 = C10*Delta_n + (C10_HH*kappa_HH_21cm(Tmat, .true.) &
        +C10_eH*kappa_eH_21cm(Tmat, .true.))*Delta_Tm

    !          write(*,'(9e15.5)') 1/a-1,Tmat,Tspin, Trad,C10_HH,C10_eH,A10*Trad/T_21cm,xe,&
    !                  n_H*kappa_pH_21cm(Tmat, .false.)*xe

    Get21cm_dTs =  4*a*( TSpin/TMat*(Delta_Tm-Delta_ts)*C10 + (1-TSpin/TMat)*delta_C10 + &
        (Trad*Delta_TCMB - Tspin*Delta_Ts)*A10/T_21cm ) * MPC_in_sec

    end function Get21cm_dTs


    subroutine output_window_sources(EV, this, sources, y, yprime, &
        tau, a, adotoa, grho, gpres, &
        k, etak, z, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc,clxnu, Delta_TM, Delta_xe,  &
        dgq, qg,  vb, qgdot, vbdot, &
        dgpi, pig, pigdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E, Edot, &
        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau)
    !Line of sight sources for number counts, lensing and 21cm redshift windows
    type(EvolutionVars) EV
    class(TThermoData) :: this
    real(dl) y(EV%nvar), yprime(EV%nvar)
    real(dL), intent(out) :: sources(:)
    real(dL), intent(in) :: tau, a, adotoa, grho, gpres, &
        k,etak, z, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc,clxnu,  &
        dgq, qg, vb, qgdot, vbdot, &
        dgpi, pig, pigdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E(2:3), Edot(2:3)
    ! and
    real(dl), intent(in) :: opacity(nscatter), dopacity(nscatter), ddopacity(nscatter), &
        visibility(nscatter), dvisibility(nscatter), ddvisibility(nscatter), exptau(nscatter)
    ! and
    real(dl), intent(in) :: Delta_TM, Delta_xe
    real(dl) s(0:10), t(0:10)
    real(dl) counts_radial_source, counts_velocity_source, counts_density_source, counts_ISW_source, &
        counts_redshift_source, counts_timedelay_source, counts_potential_source
    integer w_ix, lineoff,lineoffpol
    real(dl) Delta_TCMB
    integer j
    real(dl) Tmat,Trad, Delta_source, Delta_source2
    real(dl) xe, chi, polter_line

    j = EV%OutputStep
    if (CP%SourceTerms%line_reionization) sources(2)=0

    ! and
    if (tau <= this%tau_start_redshiftwindows) return
    ! and

    !There are line of sight contributions...
    if (CP%Do21cm) then
        Delta_TCMB = clxg/4
        Delta_source = clxb
        Trad = CP%TCMB/a

        xe = CP%Recomb%x_e(a)
        Tmat = CP%Recomb%T_m(a)

        Delta_source2 = Get21cm_source2(a,Delta_source,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe, &
            k*(z+vb)/adotoa/3)
    end if

    do w_ix = 1, State%num_redshiftwindows
        associate (W => State%Redshift_W(w_ix))

            if (W%kind == window_lensing) then
                sources(3+w_ix) =-2*phi*W%win_lens(j)
            elseif (W%kind == window_counts) then
                !assume zero velocity bias and relevant tracer is CDM perturbation
                !neglect anisotropic stress in some places

                !Main density source
                if (CP%SourceTerms%counts_density) then
                    counts_density_source= W%wing(j)*(clxc*W%Window%GetBias(k,a) + (W%comoving_density_ev(j) - 3*adotoa)*sigma/k)
                    !Newtonian gauge count density; bias assumed to be on synchronous gauge CDM density
                else
                    counts_density_source= 0
                endif

                if (CP%SourceTerms%counts_redshift) then
                    !Main redshift distortion from kV_N/H j'' integrated by parts twice (V_N = sigma in synch gauge)
                    counts_redshift_source = ((4.D0*adotoa**2+gpres+grho/3.D0)/k*W%wing2(j) + &
                        (-4.D0*W%dwing2(j)*adotoa+W%ddwing2(j))/k)*sigma+(-etak/adotoa*k/3.D0-dgrho/ &
                        adotoa/6.D0+(etak/adotoa*k/3.D0+dgrho/adotoa/6.D0+(dgq/2.D0-2.D0*etak*adotoa)/k) &
                        /EV%Kf(1))*W%wing2(j)+2.D0*W%dwing2(j)*etak/k/EV%Kf(1)
                else
                    counts_redshift_source= 0
                end if

                ! 2v j'/(H\chi) geometric term
                if (State%tau0-tau > 0.1_dl .and. CP%SourceTerms%counts_radial) then
                    chi =State%tau0-tau
                    counts_radial_source= (1-2.5*W%Window%dlog10Ndm)*((-4.D0*W%wing2(j)/chi*adotoa &
                        -2.D0*(-W%dwing2(j)*chi-W%wing2(j))/chi**2)/ &
                        k*sigma+2.D0*W%wing2(j)*etak/chi/k/EV%Kf(1))
                else
                    counts_radial_source = 0
                end if

                if (CP%SourceTerms%counts_timedelay) then
                    !time delay; WinV is int g/chi
                    counts_timedelay_source= 2*(1-2.5*W%Window%dlog10Ndm)*W%WinV(j)*2*phi
                else
                    counts_timedelay_source = 0
                end if

                if (CP%SourceTerms%counts_ISW) then
                    !WinF is int wingtau
                    counts_ISW_source = W%WinF(j)*2*phidot
                else
                    counts_ISW_source = 0
                end if

                if (CP%SourceTerms%counts_potential) then
                    !approx phi = psi
                    counts_potential_source = ( phidot/adotoa + phi +(5*W%Window%dlog10Ndm-2)*phi ) * W%wing(j) &
                        + phi * W%wingtau(j)
                else
                    counts_potential_source = 0
                end if

                if (CP%SourceTerms%counts_velocity) then
                    counts_velocity_source =  (-2.D0*W%wingtau(j)*adotoa+W%dwingtau(j))/k*sigma &
                        +W%wingtau(j)*etak/k/EV%Kf(1) &
                        - counts_radial_source  !don't double count terms; counts_radial is part of counts_velocity with 1/H/chi
                else
                    counts_velocity_source = 0
                end if

                sources(3+w_ix)=  counts_radial_source +  counts_density_source + counts_redshift_source &
                    + counts_timedelay_source + counts_potential_source &
                    + counts_ISW_source + counts_velocity_source

                sources(3+w_ix)=sources(3+w_ix)/W%Fq

                if (CP%SourceTerms%counts_lensing) &
                    sources(3+W%mag_index+State%num_redshiftwindows) = phi*W%win_lens(j)*(2-5*W%Window%dlog10Ndm)
            elseif (W%kind == window_21cm) then
                if (CP%SourceTerms%line_basic) then
                    sources(3+w_ix)= exptau(1)*(W%wing(j)*Delta_source + W%wing2(j)*Delta_source2 &
                        - W%Wingtau(j)*(clxb - (Delta_source2+clxg/4)))
                    !!    sources(3+w_ix)= exptau*W%wing(j)*phi
                else
                    sources(3+w_ix)= 0
                end if

                if (CP%SourceTerms%line_distortions ) then
                    !With baryon velocity, dropping small terms
                    s(1) =  (sigma/adotoa/3.D0-etak/adotoa**2/3.D0)*W%wing(j)*exptau(1)*k
                    s(2) =  -1.D0/adotoa**2*exptau(1)*W%wing(j)*dgrho/6.D0+((((4.D0*sigma+ &
                        vb)*adotoa+(-grho*sigma/2.D0-vb*grho/3.D0)/adotoa+(sigma*grho**2/18.D0 + &
                        vb*grho**2/18.D0)/adotoa**3)*W%wing(j)-4.D0*W%dwing(j)*sigma+(W%ddwing(j)*sigma + &
                        W%ddwing(j)*vb)/adotoa+(W%dwing(j)*sigma*grho/3.D0+W%dwing(j)*vb*grho/3.D0)/ &
                        adotoa**2-2.D0*W%dwing(j)*vb+((-2.D0*etak+etak*grho/adotoa**2/3.D0)*W%wing(j) &
                        + 2.D0*W%dwing(j)*etak/adotoa)/EV%Kf(1))*exptau(1) + &
                        (-4.D0*visibility(1)*sigma- 2.D0*visibility(1)*vb + &
                        (dvisibility(1)*sigma+dvisibility(1)*vb)/adotoa+(visibility(1)*grho*sigma/3.D0 + &
                        visibility(1)*vb*grho/3.D0)/adotoa**2)*W%wing(j)+2.D0*visibility(1)*etak/adotoa*W%wing(j)/ &
                        EV%Kf(1)+(2.D0*visibility(1)*W%dwing(j)*sigma+2.D0*visibility(1)*W%dwing(j)*vb)/adotoa)/k
                    t(0) =  s(1)+s(2)

                    sources(3+w_ix)= sources(3+w_ix) + t(0)
                end if


                if (CP%SourceTerms%line_extra) then
                    !All sources except below
                    if (CP%SourceTerms%line_basic .and. CP%SourceTerms%line_distortions) then
                        sources(3+w_ix) =  (-2.D0/3.D0*sigma+2.D0/3.D0*etak/adotoa)*W%winV(j)*exptau(1)*k+ &
                            (W%wing2(j)*Delta_source2+W%wing(j)*Delta_source+1.D0/adotoa*W%winV(j)*dgrho/3.D0)* &
                            exptau(1)+((-W%dwing(j)*vb+(-(3.D0*gpres+grho)*sigma/3.D0 &
                            - 4.D0*adotoa**2*sigma)*W%winV(j)+4.D0*adotoa*W%dwinV(j)*sigma+(-sigma- &
                            vb)*W%ddWinV(j)-vbdot*W%wing(j)-W%dwinV(j)*vbdot+(-2.D0*W%dwinV(j)*etak &
                            + 2.D0*etak*adotoa*W%winV(j))/EV%Kf(1))*exptau(1)-2.D0*visibility(1)*sigma*W%dwinV(j)+ &
                            (4.D0*visibility(1)*sigma*adotoa-dvisibility(1)*sigma)*W%winV(j)-2.D0*visibility(1)*W%winV(j)*etak/ &
                            EV%Kf(1)-visibility(1)*W%dwinV(j)*vb-visibility(1)*W%wing(j)*vb)/k+((2.D0*W%dwinV(j)*dgpi+ &
                            diff_rhopi*W%winV(j))*exptau(1)+2.D0*visibility(1)*W%winV(j)*dgpi)/k**2
                    else
                        s(1) =  ((-2.D0/3.D0*sigma+2.D0/3.D0*etak/adotoa)*W%winV(j)+(-sigma/adotoa/3.D0 + &
                            etak/adotoa**2/3.D0)*W%wing(j))*exptau(1)*k+(1.D0/adotoa*W%winV(j)*dgrho/3.D0 &
                            + 1.D0/adotoa**2*W%wing(j)*dgrho/6.D0)*exptau(1)
                        s(2) =  s(1)
                        s(6) =  ((-vb-sigma)*W%ddWinV(j)+(-4.D0*adotoa**2*sigma-&
                            (18.D0*gpres+ 6.D0*grho)*sigma/18.D0)*W%winV(j)+((-4.D0*sigma-vb)*adotoa-vbdot + &
                            (grho*sigma/ 2.D0+vb*grho/3.D0)/adotoa+(-grho**2*sigma/18.D0-vb*grho**2/18.D0)/ &
                            adotoa**3)*W%wing(j)+W%dwing(j)*vb+(-W%ddwing(j)*sigma-W%ddwing(j)*vb)/adotoa &
                            + 4.D0*W%dwinV(j)*sigma*adotoa+4.D0*W%dwing(j)*sigma+(-W%dwing(j)*grho*sigma/3.D0 - &
                            W%dwing(j)*vb*grho/3.D0)/adotoa**2-W%dwinV(j)*vbdot+((2.D0*etak-etak*grho/ &
                            adotoa**2/3.D0)*W%wing(j)-2.D0*W%dwing(j)*etak/adotoa-2.D0*W%dwinV(j)*etak &
                            + 2.D0*etak*adotoa*W%winV(j))/EV%Kf(1))*exptau(1)-visibility(1)*W%dwinV(j)*vb + &
                            (4.D0*visibility(1)*sigma*adotoa-dvisibility(1)*sigma)*W%winV(j)
                        s(5) =  s(6)+(-2.D0*visibility(1)*etak/adotoa*W%wing(j)-2.D0*visibility(1)*W%winV(j)*etak)/ &
                            EV%Kf(1)+(4.D0*visibility(1)*sigma+(-visibility(1)*grho*sigma/3.D0-visibility(1)*vb*grho/3.D0)/ &
                            adotoa**2+visibility(1)*vb+(-dvisibility(1)*sigma-dvisibility(1)*vb)/adotoa)*W%wing(j)+ &
                            (-2.D0*visibility(1)*W%dwing(j)*sigma-2.D0*visibility(1)*W%dwing(j)*vb)/adotoa &
                            - 2.D0*visibility(1)*W%dwinV(j)*sigma
                        s(6) =  1.D0/k
                        s(4) =  s(5)*s(6)
                        s(5) =  ((diff_rhopi*W%winV(j)+2.D0*W%dwinV(j)*dgpi)*exptau(1) &
                            + 2.D0*visibility(1)*dgpi*W%winV(j))/k**2
                        s(3) =  s(4)+s(5)
                        t(0) =  s(2)+s(3)

                        sources(3+w_ix) =   sources(3+w_ix) + t(0)
                    end if
                end if


                if (CP%SourceTerms%line_reionization) then
                    if (State%num_redshiftwindows>1) stop 'reionization only for one window at the mo'
                    lineoff=EV%reion_line_ix
                    lineoffpol = lineoff+EV%lmaxline-1

                    if (tau  < State%tau0) then
                        polter_line = 0.1_dl*y(lineoff+2)+9._dl/15._dl*y(lineoffpol+2)
                        sources(2)=visibility(1)*polter_line*(15._dl/2._dl)/(f_K(State%tau0-tau)*k)**2
                    else
                        sources(2)=0
                    end if

                    if (.not. CP%SourceTerms%use_21cm_mK) sources(2)= sources(2) /W%Fq

                    s(1) =  visibility(1)*y(lineoff+2)/4.D0+visibility(1)*y(lineoff)
                    s(2) =  s(1)
                    s(4) =  (-1.D0/EV%Kf(1)*visibility(1)*W%winV(j)*etak/10.D0-visibility(1)*sigma*W%dwinV(j)/10.D0 &
                        - 9.D0/20.D0*visibility(1)*yprime(lineoff+2)-27.D0/100.D0*visibility(1)*opacity(1)*y(lineoff+1) &
                        - 9.D0/10.D0*dvisibility(1)*y(lineoff+3)-3.D0/20.D0*visibility(1)*opacity(1)*EV%Kf(2)*y(lineoffpol+3)+ &
                        visibility(1)*W%dwinV(j)*vb+81.D0/200.D0*visibility(1)*opacity(1)*y(lineoff+3) &
                        +3.D0/5.D0*dvisibility(1)*y(lineoff+1)+3.D0/10.D0*visibility(1)*yprime(lineoff+1)+ &
                        (visibility(1)*adotoa*sigma/5.D0+(36.D0*visibility(1)*opacity(1)-80.D0*dvisibility(1))*sigma/400.D0+ &
                        dvisibility(1)*vb+visibility(1)*vbdot)*W%winV(j))/k
                    s(5) =  (visibility(1)*W%winV(j)*dgpi/10.D0+9.D0/20.D0*visibility(1)*dopacity(1)*y(lineoffpol+2) &
                        + 261.D0/400.D0*visibility(1)*opacity(1)**2.D0*y(lineoff+2)&
                        -117.D0/200.D0*visibility(1)*opacity(1)**2.D0*y(lineoffpol+2)+3.D0/4.D0*ddvisibility(1)*y(lineoff+2) &
                        - 27.D0/20.D0*dvisibility(1)*opacity(1)*y(lineoff+2)+9.D0/10.D0*dvisibility(1)*opacity(1)*y(lineoffpol+2)&
                        -27.D0/40.D0*visibility(1)*dopacity(1)*y(lineoff+2))/k**2
                    s(3) =  s(4)+s(5)
                    t(0) =  s(2)+s(3)

                    sources(3+w_ix)= sources(3+w_ix) + t(0)
                end if

                if (CP%SourceTerms%line_phot_quadrupole) then
                    s(1) =  (EV%kf(1)*W%wing2(j)*pig/2.D0+(-clxg/4.D0-5.D0/8.D0*pig)*W%wing2(j))*exptau(1)
                    s(3) =  ((-1.D0/EV%kf(1)*W%wing2(j)*etak+(-sigma+9.D0/8.D0*EV%kf(2)*y(9) &
                        -3.D0/4.D0*qg)*W%dwing2(j)+(-opacity(1)*vb+2.D0*adotoa*sigma+9.D0/8.D0*EV%kf(2)*yprime(9) &
                        + 3.D0/8.D0*opacity(1)*EV%kf(2)*E(3)+3.D0/4.D0*opacity(1)*qg)*W%wing2(j))*exptau(1)+ &
                        (-3.D0/4.D0*visibility(1)*qg-visibility(1)*sigma+9.D0/8.D0*visibility(1)*EV%kf(2)*y(9))*W%wing2(j))/k
                    s(4) =  (((27.D0/16.D0*opacity(1)*pig-15.D0/8.D0*pigdot &
                        -9.D0/8.D0*opacity(1)*E(2))*W%dwing2(j)+(27.D0/16.D0*dopacity(1)*pig &
                        +9.D0/8.D0*opacity(1)**2.D0*E(2)-9.D0/8.D0*opacity(1)**2.D0*polter &
                        +27.D0/16.D0*opacity(1)*pigdot+dgpi-9.D0/8.D0*dopacity(1)*E(2))*W%wing2(j)&
                        -15.D0/8.D0*W%ddwing2(j)*pig)*exptau(1)-15.D0/4.D0*visibility(1)*W%dwing2(j)*pig+(- &
                        (-27.D0*visibility(1)*opacity(1)+30.D0*dvisibility(1))*pig/16.D0-9.D0/8.D0*visibility(1)*opacity(1)*E(2) &
                        - 15.D0/8.D0*visibility(1)*pigdot)*W%wing2(j))/k**2
                    s(2) =  s(3)+s(4)
                    t(0) =  s(1)+s(2)

                    sources(3+w_ix)= sources(3+w_ix)+ t(0)
                end if


                if (CP%SourceTerms%line_phot_dipole) then
                    sources(3+w_ix)=sources(3+w_ix) + (EV%kf(1)*W%wing2(j)*pig/2.D0-W%wing2(j)*clxg/4.D0)*exptau(1) &
                        +(((vbdot- opacity(1)*vb+3.D0/4.D0*opacity(1)*qg)*&
                        W%wing2(j)+(vb-3.D0/4.D0*qg)*W%dwing2(j))*exptau(1)+&
                        (visibility(1)*vb-3.D0/4.D0*visibility(1)*qg)*W%wing2(j))/k
                end if

                if (.not. CP%SourceTerms%use_21cm_mK) sources(3+w_ix)= sources(3+w_ix) /W%Fq
            end if
        end associate
    end do
    end subroutine output_window_sources
    ! andrea
    subroutine output(EV, this, etat, CLdata, yin, j, tau, sources, num_custom_sources)
    ! andrea
    Type(EvolutionVars) EV
    class(TThermoData) :: this
    class(TRecfast) :: etat
    class(TCLData) :: CLData
    ! andrea
    real(dl) yin(EV%nvar),yprimein(EV%nvar)
    ! andrea
    real(dl) y(EV%nvar), yprime(EV%nvar)
    integer, intent(in) :: j
    real(dl) tau
    real(dl), target :: sources(:)
    integer, intent(in) :: num_custom_sources
    
    ! andrea
    y(1:EV%nvar)=yin(1:EV%nvar)
    yprimein = 0
    EV%OutputSources => Sources
    EV%OutputStep = j
    if (num_custom_sources>0) &
        EV%CustomSources => sources(CLdata%CTransScal%NumSources - num_custom_sources+1:)
    call derivs(EV, this, etat, EV%ScalEqsToPropagate,tau,y,yprimein)
    yprime(1:EV%nvar) = yprimein(1:EV%nvar)
    ! andrea
    !write(*,*) 'Equations Sources', Sources
    nullify(EV%OutputSources, EV%CustomSources)

    end subroutine output

    ! andrea
    subroutine outputt(EV, this, etat ,ytin,n,tau,dt,dte,dtb)
    ! andrea
    !calculate the tensor sources for open and closed case
    implicit none
    integer n
    type(EvolutionVars) :: EV
    class(TThermoData) :: this
    class(TRecfast) :: etat
    real(dl), target :: yt(n), ytprime(n)
    ! andrea
    real(dl)  :: ytin(n), ytprimein(n)
    real(dl) tau,dt(nscatter),dte(nscatter),dtb(nscatter),x,polterdot,polterddot,prefac
    ! andrea
    real(dl) pig, pigdot, octg, aux, polter, shear, adotoa,a
    real(dl) sinhxr,cothxor
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Bprime,Eprime
    real(dl), target :: pol(3),polEprime(3), polBprime(3)
    ! andrea
    real(dl) :: lenswindow
    real(dl) :: opacity(nscatter), dopacity(nscatter), ddopacity(nscatter), &
        visibility(nscatter), dvisibility(nscatter), ddvisibility(nscatter), exptau(nscatter)
    ! and
    integer :: f_i, ix_off
    yt(1:EV%nvart)=ytin(1:EV%nvart)
    ytprimein = 0
    call derivst(EV, this, etat,EV%nvart,tau,ytin,ytprimein)
    ytprime(1:EV%nvart)=ytprimein(1:EV%nvart)
    ! andrea
    k2=EV%k2_buf
    k=EV%k_buf
    aux=EV%aux_buf
    shear = yt(ixt_shear)

    x=(State%tau0-tau)/State%curvature_radius
    call this%IonizationFunctionsAtTime(tau, a, opacity, dopacity, ddopacity, &
        visibility, dvisibility, ddvisibility, exptau, lenswindow)

    !  And the electric part of the Weyl.
    if (.not. EV%TensTightCoupling) then
        !  Use the full expression for pigdt
        pig=yt(EV%g_ix+2)
        pigdot=ytprime(EV%g_ix+2)
        E => yt(EV%E_ix+1:)
        Eprime=> ytprime(EV%E_ix+1:)
        Bprime => ytprime(EV%B_ix+1:)
        octg=ytprime(EV%g_ix+3)
    else
        !  Use the tight-coupling approximation
        adotoa = 1/(a*dtauda(State,a))
        ! and
        pigdot=32._dl/45._dl*k/opacity(1)*(2._dl*adotoa*shear+ytprime(ixt_shear))
        pig = 32._dl/45._dl*k/opacity(1)*shear
        ! and
        pol=0
        polEprime=0
        polBprime=0
        E=>pol
        EPrime=>polEPrime
        BPrime=>polBPrime
        E(2)=pig/4._dl
        EPrime(2)=pigdot/4
        octg=0
    endif

    sinhxr=State%rofChi(x)*State%curvature_radius

    if (EV%q*sinhxr > 1.e-8_dl) then
        ! andrea
        do f_i = 1, nscatter
            if (EV%Rayleigh .and. f_i>1) then
                ix_off=EV%g_ix_freq+2+(f_i-2)*EV%freq_neq
                yt(EV%g_ix+2:EV%g_ix+2+EV%freq_neq-1) = ytin(EV%g_ix+2:EV%g_ix+2+EV%freq_neq-1)+ &
                ytin(ix_off:ix_off+EV%freq_neq-1)
                ytprime(EV%g_ix+2:EV%g_ix+2+EV%freq_neq-1) = ytprimein(EV%g_ix+2:EV%g_ix+2+EV%freq_neq-1)+ &
                ytprimein(ix_off:ix_off+EV%freq_neq-1)
                pig =yt(EV%g_ix+2)
                pigdot=ytprime(EV%g_ix+2)
                octg=ytprime(EV%g_ix+3)
            end if
            ! andrea
            prefac=sqrt(EV%q2*State%curvature_radius*State%curvature_radius-State%Ksign)
            cothxor=State%cosfunc(x)/sinhxr
            ! and freestyle
            polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
            polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*pigdot
            polterddot = 9._dl/15._dl*(-dopacity(f_i)*(E(2)-polter)-opacity(f_i)*(  &
            Eprime(2)-polterdot) + k*(2._dl/3._dl*Bprime(2)*aux - 5._dl/27._dl*Eprime(3)*EV%Kft(2))) &
            +0.1_dl*(k*(-octg*EV%Kft(2)/3._dl + 8._dl/15._dl*ytprime(ixt_shear)) - &
            dopacity(f_i)*(pig - polter) - opacity(f_i)*(pigdot-polterdot))
            ! and freestyle

            ! and intuition ?
            dt(f_i)=(shear*exptau(f_i) + (15._dl/8._dl)*polter*visibility(f_i)/k)*State%curvature_radius/sinhxr**2/prefac
            dte(f_i)=State%curvature_radius*15._dl/8._dl/k/prefac* &
            ((ddvisibility(f_i)*polter + 2._dl*dvisibility(f_i)*polterdot + visibility(f_i)*polterddot)  &
            + 4._dl*cothxor*(dvisibility(f_i)*polter + visibility(f_i)*polterdot) - &
            visibility(f_i)*polter*(k2 -6*cothxor**2))
            dtb(f_i)=15._dl/4._dl*EV%q*State%curvature_radius/k/prefac*(visibility(f_i)*(2._dl*cothxor*polter + polterdot) + &
 dvisibility(f_i)*polter)

            if (rayleigh_diff .and. f_i>2) then
                dt(f_i)= dt(f_i) - dt(1)
                dte(f_i)= dte(f_i) - dte(1)
                dtb(f_i)= dtb(f_i) - dtb(1)
            end if
        end do
        ! andrea
    else
        ! and fix ?
        dt(1)=0._dl
        dte(1)=0._dl
        dtb(1)=0._dl
        ! and
    end if

    end subroutine outputt


    subroutine outputv(EV,this,etat,yv,n,tau,dt,dte,dtb)
    !calculate the vector sources
    implicit none
    integer n
    type(EvolutionVars) :: EV
    class(TThermoData) :: this
    class(TRecfast) :: etat
    real(dl), target :: yv(n), yvprime(n)
    ! and
    real(dl) tau,x,polterdot
    real(dl) dt, dte, dtb
    ! and
    real(dl) vb,qg, pig, polter, sigma
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Eprime
    real(dl) a, lenswindow
    integer, parameter :: num_cmb_freq = 6 !!!
    integer, parameter :: nscatter = num_cmb_freq+1
    real(dl) :: opacity(nscatter), dopacity(nscatter), ddopacity(nscatter), exptau(nscatter), &
                        visibility(nscatter), dvisibility(nscatter), ddvisibility(nscatter)


    call derivsv(EV,this,etat,EV%nvarv,tau,yv,yvprime)

    k2=EV%k2_buf
    k=EV%k_buf
    sigma = yv(2)
    vb  = yv(3)
    qg  = yv(4)
    pig = yv(5)


    x=(State%tau0-tau)*k

    if (x > 1.e-8_dl) then
        E => yv(EV%lmaxv+3:)
        Eprime=> yvprime(EV%lmaxv+3:)

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
        polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*yvprime(5)

        call this%IonizationFunctionsAtTime(tau, a, opacity, dopacity, ddopacity, &
            visibility, dvisibility, ddvisibility, exptau, lenswindow)


        if (yv(1) < 1e-3) then
            dt = 1
        else
            dt =0
        end if
        dt= (4*(vb+sigma)*visibility(1) + 15._dl/2/k*( visibility(1)*polterdot + dvisibility(1)*polter) &
            + 4*(exptau(1)*yvprime(2)) )/x

        dte= 15._dl/2*2*polter/x**2*visibility(1) + 15._dl/2/k*(dvisibility(1)*polter + visibility(1)*polterdot)/x

        dtb= -15._dl/2*polter/x*visibility(1)
    else
        dt=0
        dte=0
        dtb=0
    end if

    end subroutine outputv

    subroutine initial(EV, this, y, tau)
    !  Scalar initial conditions.
    implicit none

    type(EvolutionVars) EV
    class(TThermoData) :: this
    real(dl) y(EV%nvar)
    real(dl) Rp15,tau,x,x2,x3,om,omtau, &
        Rc,Rb,Rv,Rg,grhonu,chi
    real(dl) k,k2
    real(dl) a,a2, iqg, rhomass,a_massive, ep
    integer l,i, nu_i, j, ind
    integer, parameter :: i_clxg=1,i_clxr=2,i_clxc=3, i_clxb=4, &
        i_qg=5,i_qr=6,i_vb=7,i_pir=8, i_eta=9, i_aj3r=10,i_clxde=11,i_vde=12
    integer, parameter :: i_max = i_vde
    real(dl) initv(6,1:i_max), initvec(1:i_max)
    integer f_i

    nullify(EV%OutputTransfer) !Should not be needed, but avoids issues in ifort 14
    nullify(EV%OutputSources)
    nullify(EV%CustomSources)

    EV%is_cosmological_constant = State%CP%DarkEnergy%is_cosmological_constant

    if (State%flat) then
        EV%k_buf=EV%q
        EV%k2_buf=EV%q2
        EV%Kf(1:EV%MaxlNeeded)=1._dl
    else
        EV%k2_buf=EV%q2-State%curv
        EV%k_buf=sqrt(EV%k2_buf)

        do l=1,EV%MaxlNeeded
            EV%Kf(l)=1._dl-State%curv*(l*(l+2))/EV%k2_buf
        end do
    end if

    k=EV%k_buf
    k2=EV%k2_buf

    do j=1,EV%MaxlNeeded
        EV%denlk(j)=denl(j)*k*j
        EV%denlk2(j)=denl(j)*k*EV%Kf(j)*(j+1)
        EV%polfack(j)=polfac(j)*k*EV%Kf(j)*denl(j)
    end do

    !Get time to switch off tight coupling
    !The numbers here are a bit of guesswork
    !The high k increase saves time for very small loss of accuracy
    !The lower k ones are more delicate. Nead to avoid instabilities at same time
    !as ensuring tight coupling is accurate enough
    if (EV%k_buf > epsw) then
        if (EV%k_buf > epsw*5) then
            ep=ep0*5/CP%Accuracy%AccuracyBoost*0.65
        else
            ep=ep0
        end if
    else
        ep=ep0
    end if
    if (second_order_tightcoupling) ep=ep*2
    EV%TightSwitchoffTime = min(this%tight_tau, this%OpacityToTime(EV%k_buf/ep))
    ! andrea
    if (num_cmb_freq>0) then
        EV%RayleighSwitchOnTime = EV%TightSwitchoffTime
    else
        EV%RayleighSwitchOnTime=State%tau0+1
    end if
    EV%Rayleigh = .false.
    ! andrea

    y=0

    !  k*tau, (k*tau)**2, (k*tau)**3
    x=k*tau
    x2=x*x
    x3=x2*x
    rhomass =  sum(State%grhormass(1:CP%Nu_mass_eigenstates))
    grhonu=rhomass+State%grhornomass

    om = (State%grhob+State%grhoc)/sqrt(3*(State%grhog+grhonu))
    omtau=om*tau
    Rv=grhonu/(grhonu+State%grhog)

    Rg = 1-Rv
    Rc=CP%omch2/(CP%omch2+CP%ombh2)
    Rb=1-Rc
    Rp15=4*Rv+15

    if (CP%Scalar_initial_condition > initial_nummodes) &
        call MpiStop('Invalid initial condition for scalar modes')

    a=tau*State%adotrad*(1+omtau/4)
    a2=a*a

    initv=0

    !  Set adiabatic initial conditions

    chi=1  !Get transfer function for chi
    initv(1,i_clxg)=-chi*EV%Kf(1)/3*x2*(1-omtau/5)
    initv(1,i_clxr)= initv(1,i_clxg)
    initv(1,i_clxb)=0.75_dl*initv(1,i_clxg)
    initv(1,i_clxc)=initv(1,i_clxb)
    initv(1,i_qg)=initv(1,i_clxg)*x/9._dl
    initv(1,i_qr)=-chi*EV%Kf(1)*(4*Rv+23)/Rp15*x3/27
    initv(1,i_vb)=0.75_dl*initv(1,i_qg)
    initv(1,i_pir)=chi*4._dl/3*x2/Rp15*(1+omtau/4*(4*Rv-5)/(2*Rv+15))
    initv(1,i_aj3r)=chi*4/21._dl/Rp15*x3
    initv(1,i_eta)=-chi*2*EV%Kf(1)*(1 - x2/12*(-10._dl/Rp15 + EV%Kf(1)))

    if (CP%Scalar_initial_condition/= initial_adiabatic) then
        !CDM isocurvature

        initv(2,i_clxg)= Rc*omtau*(-2._dl/3 + omtau/4)
        initv(2,i_clxr)=initv(2,i_clxg)
        initv(2,i_clxb)=initv(2,i_clxg)*0.75_dl
        initv(2,i_clxc)=1+initv(2,i_clxb)
        initv(2,i_qg)=-Rc/9*omtau*x
        initv(2,i_qr)=initv(2,i_qg)
        initv(2,i_vb)=0.75_dl*initv(2,i_qg)
        initv(2,i_pir)=-Rc*omtau*x2/3/(2*Rv+15._dl)
        initv(2,i_eta)= Rc*omtau*(1._dl/3 - omtau/8)*EV%Kf(1)
        initv(2,i_aj3r)=0
        !Baryon isocurvature
        if (Rc==0) call MpiStop('Isocurvature initial conditions assume non-zero dark matter')

        initv(3,:) = initv(2,:)*(Rb/Rc)
        initv(3,i_clxc) = initv(3,i_clxb)
        initv(3,i_clxb)= initv(3,i_clxb)+1

        !neutrino isocurvature density mode

        initv(4,i_clxg)=Rv/Rg*(-1 + x2/6)
        initv(4,i_clxr)=1-x2/6
        initv(4,i_clxc)=-omtau*x2/80*Rv*Rb/Rg
        initv(4,i_clxb)= Rv/Rg/8*x2
        iqg = - Rv/Rg*(x/3 - Rb/4/Rg*omtau*x)
        initv(4,i_qg) =iqg
        initv(4,i_qr) = x/3
        initv(4,i_vb)=0.75_dl*iqg
        initv(4,i_pir)=x2/Rp15
        initv(4,i_eta)=EV%Kf(1)*Rv/Rp15/3*x2

        !neutrino isocurvature velocity mode

        initv(5,i_clxg)=Rv/Rg*x - 2*x*omtau/16*Rb*(2+Rg)/Rg**2
        initv(5,i_clxr)=-x -3*x*omtau*Rb/16/Rg
        initv(5,i_clxc)=-9*omtau*x/64*Rv*Rb/Rg
        initv(5,i_clxb)= 3*Rv/4/Rg*x - 9*omtau*x/64*Rb*(2+Rg)/Rg**2
        iqg = Rv/Rg*(-1 + 3*Rb/4/Rg*omtau+x2/6 +3*omtau**2/16*Rb/Rg**2*(Rg-3*Rb))
        initv(5,i_qg) =iqg
        initv(5,i_qr) = 1 - x2/6*(1+4*EV%Kf(1)/(4*Rv+5))
        initv(5,i_vb)=0.75_dl*iqg
        initv(5,i_pir)=2*x/(4*Rv+5)+omtau*x*6/Rp15/(4*Rv+5)
        initv(5,i_eta)=2*EV%Kf(1)*x*Rv/(4*Rv+5) + omtau*x*3*EV%Kf(1)*Rv/32*(Rb/Rg - 80/Rp15/(4*Rv+5))
        initv(5,i_aj3r) = 3._dl/7*x2/(4*Rv+5)

        !quintessence isocurvature mode
    end if

    if (CP%Scalar_initial_condition==initial_vector) then
        InitVec = 0
        do i=1,initial_nummodes
            InitVec = InitVec+ initv(i,:)*CP%InitialConditionVector(i)
        end do
    else
        InitVec = initv(CP%Scalar_initial_condition,:)
        if (CP%Scalar_initial_condition==initial_adiabatic) InitVec = -InitVec
        !So we start with chi=-1 as before
    end if

    y(ix_etak)= -InitVec(i_eta)*k/2
    !get eta_s*k, where eta_s is synchronous gauge variable

    !  CDM
    y(ix_clxc)=InitVec(i_clxc)

    !  Baryons
    y(ix_clxb)=InitVec(i_clxb)
    y(ix_vb)=InitVec(i_vb)

    !  Photons
    y(EV%g_ix)=InitVec(i_clxg)
    y(EV%g_ix+1)=InitVec(i_qg)

    ! DarkEnergy: This initializes also i_vq, when num_perturb_equations is set
    !             to 2.
    if (CP%DarkEnergy%num_perturb_equations > 0) then
        call CP%DarkEnergy%PerturbationInitial(InitVec(i_clxde:i_clxde + CP%DarkEnergy%num_perturb_equations - 1), &
            a, tau,  k)
        y(EV%w_ix:EV%w_ix + CP%DarkEnergy%num_perturb_equations - 1) = &
            InitVec(i_clxde:i_clxde + CP%DarkEnergy%num_perturb_equations - 1)
    end if

    if (CP%Evolve_delta_Ts) then
        y(EV%Ts_ix) = y(EV%g_ix)/4
    end if

    !  Neutrinos
    y(EV%r_ix)=InitVec(i_clxr)
    y(EV%r_ix+1)=InitVec(i_qr)
    y(EV%r_ix+2)=InitVec(i_pir)

    if (EV%lmaxnr>2) then
        y(EV%r_ix+3)=InitVec(i_aj3r)
    endif

    if (CP%Num_Nu_massive == 0) return

    do nu_i = 1, CP%Nu_mass_eigenstates
        EV%MassiveNuApproxTime(nu_i) = Nu_tau_massive(nu_i)
        a_massive =  20000*k/State%nu_masses(nu_i)*CP%Accuracy%AccuracyBoost*CP%Accuracy%lAccuracyBoost
        if (a_massive >=0.99) then
            EV%MassiveNuApproxTime(nu_i)=State%tau0+1
        else if (a_massive > 17.d0/State%nu_masses(nu_i)*CP%Accuracy%AccuracyBoost) then
            EV%MassiveNuApproxTime(nu_i)=max(EV%MassiveNuApproxTime(nu_i),State%DeltaTime(0._dl,a_massive, 0.01_dl))
        end if
        ind = EV%nu_ix(nu_i)
        do  i=1,EV%nq(nu_i)
            y(ind:ind+2)=y(EV%r_ix:EV%r_ix+2)
            if (EV%lmaxnu_tau(nu_i)>2) y(ind+3)=InitVec(i_aj3r)
            ind = ind + EV%lmaxnu_tau(nu_i)+1
        end do
    end do

    end subroutine initial


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initialt(EV, this, yt, tau)
    !  Initial conditions for tensors
    implicit none
    type(EvolutionVars) EV
    ! and
    class(TThermoData)  :: this
    ! and
    real(dl) bigR,tau,x,aj3r,elec, pir, rhomass
    integer l
    real(dl) k,k2 ,a, omtau
    real(dl) yt(EV%nvart)
    real(dl) tens0, ep, tensfac

    if (State%flat) then
        EV%aux_buf=1._dl
        EV%k2_buf=EV%q2
        EV%k_buf=EV%q
        EV%Kft(1:EV%MaxlNeededt)=1._dl !initialize for flat case
    else
        EV%k2_buf=EV%q2-3*State%curv
        EV%k_buf=sqrt(EV%k2_buf)
        EV%aux_buf=sqrt(1._dl+3*State%curv/EV%k2_buf)
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    do l=1,EV%MaxlNeededt
        if (.not. State%flat) EV%Kft(l)=1._dl-State%curv*((l+1)**2-3)/k2
        EV%denlkt(1,l)=k*denl(l)*l !term for L-1
        tensfac=real((l+3)*(l-1),dl)/(l+1)
        EV%denlkt(2,l)=k*denl(l)*tensfac*EV%Kft(l) !term for L+1
        EV%denlkt(3,l)=k*denl(l)*tensfac**2/(l+1)*EV%Kft(l) !term for polarization
        EV%denlkt(4,l)=k*4._dl/(l*(l+1))*EV%aux_buf !other for polarization
    end do

    if (k > 0.06_dl*epsw) then
        ep=ep0
    else
        ep=0.2_dl*ep0
    end if

    !    finished_tightcoupling = ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep))
    ! and
    EV%TightSwitchoffTime = min(this%tight_tau,this%OpacityToTime(EV%k_buf/ep))
    ! and

    rhomass =  sum(State%grhormass(1:CP%Nu_mass_eigenstates))
    omtau = tau*(State%grhob+State%grhoc)/sqrt(3*(State%grhog+rhomass+State%grhornomass))
    a=tau*State%adotrad*(1+omtau/4)

    if (DoTensorNeutrinos) then
        bigR = (rhomass+State%grhornomass)/(rhomass+State%grhornomass+State%grhog)
    else
        bigR = 0._dl
    end if

    x=k*tau

    tens0 = 1

    yt(ixt_H)= tens0
    !commented things are for the compensated mode with magnetic fields; can be neglected
    !-15/28._dl*x**2*(bigR-1)/(15+4*bigR)*Magnetic*(1-5./2*omtau/(2*bigR+15))

    elec=-tens0*(1+2*State%curv/k2)*(2*bigR+10)/(4*bigR+15) !elec, with H=1

    !shear
    yt(ixt_shear)=-5._dl/2/(bigR+5)*x*elec
    !          + 15._dl/14*x*(bigR-1)/(4*bigR+15)*Magnetic*(1 - 15./2*omtau/(2*bigR+15))

    yt(ixt_shear+1:EV%nvart)=0._dl

    !  Neutrinos
    if (DoTensorNeutrinos) then
        pir=-2._dl/3._dl/(bigR+5)*x**2*elec
        !           + (bigR-1)/bigR*Magnetic*(1-15./14*x**2/(15+4*bigR))
        aj3r=  -2._dl/21._dl/(bigR+5)*x**3*elec !&
        !           + 3._dl/7*x*(bigR-1)/bigR*Magnetic
        yt(EV%r_ix+2)=pir
        yt(EV%r_ix+3)=aj3r
        !Should set up massive too, but small anyway..
    end if

    end subroutine initialt

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initialv(EV,yv,tau)
    !  Initial conditions for vectors

    implicit none
    real(dl) bigR,Rc,tau,x,pir
    type(EvolutionVars) EV
    real(dl) k,k2 ,a, omtau
    real(dl) yv(EV%nvarv)

    if (State%flat) then
        EV%k2_buf=EV%q2
        EV%k_buf=EV%q
    else
        call MpiStop('Vectors not supported in non-flat models')
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    omtau = tau*(State%grhob+State%grhoc)/sqrt(3*(State%grhog+State%grhornomass))

    a=tau*State%adotrad*(1+omtau/4)

    x=k*tau

    bigR = (State%grhornomass)/(State%grhornomass+State%grhog)
    Rc=CP%omch2/(CP%omch2+CP%ombh2)

    yv(1)=a !Could eliminate this, but rarely used anyway

    yv(2)= vec_sig0*(1- 15._dl/2*omtau/(4*bigR+15)) + 45._dl/14*x*Magnetic*(BigR-1)/(4*BigR+15)
    !qg
    yv(4)= vec_sig0/3* (4*bigR + 5)/(1-BigR)*(1  -0.75_dl*omtau*(Rc-1)/(bigR-1)* &
        (1 - 0.25_dl*omtau*(3*Rc-2-bigR)/(BigR-1))) &
        -x/2*Magnetic
    yv(3)= 3._dl/4*yv(4)

    yv(5:EV%nvarv) = 0

    !        if (.false.) then
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = vec_sig0/6/bigR*x**2*(1+2*bigR*omtau/(4*bigR+15))
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+2) = -2/3._dl*vec_sig0/bigR*x*(1 +3*omtau*bigR/(4*bigR+15))
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+3) = 1/4._dl*vec_sig0/bigR*(5+4*BigR)
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+4) =1/9.*x*vec_sig0*(5+4*bigR)/bigR
    !         yv(4) = 0
    !         yv(3)= 3._dl/4*yv(4)
    !          return
    !        end if

    !  Neutrinos
    !q_r
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = -1._dl/3*vec_sig0*(4*BigR+5)/bigR &
        + x**2*vec_sig0/6/BigR +0.5_dl*x*(1/bigR-1)*Magnetic
    !pi_r
    pir=-2._dl/3._dl*x*vec_sig0/BigR - (1/bigR-1)*Magnetic
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +1)=pir
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +2)=3._dl/7*x*Magnetic*(1-1/BigR)

    end subroutine initialv


    subroutine outtransf(EV, this, etat, y, tau, Arr)
    !write out clxc, clxb, clxg, clxn
    use Transfer
    implicit none
    type(EvolutionVars) EV
    ! and au pif
    class(TThermoData) :: this
    class(TRecfast) :: etat
    ! and
    real(dl), intent(in) :: tau
    real, target :: Arr(:)
    real(dl) y(EV%nvar),yprime(EV%nvar)

    yprime = 0
    EV%OutputTransfer =>  Arr
    call derivs(EV, this, etat, EV%ScalEqsToPropagate,tau,y,yprime)
    nullify(EV%OutputTransfer)
    Arr(Transfer_kh+1:Transfer_max) = Arr(Transfer_kh+1:Transfer_max)/EV%k2_buf

    end subroutine outtransf

    subroutine derivs(EV, this, etat, n,tau,ay,ayprime)
    !  Evaluate the time derivatives of the scalar perturbations
    use constants, only : barssc0, Compton_CT, line21_const
    use MassiveNu
    use Recombination
    implicit none
    type(EvolutionVars) EV
    ! and au pif
    type(TThermoData) :: this
    type(TRecfast) :: etat
    ! and 
    integer n,nu_i
    real(dl) ay(n),ayprime(n)
    real(dl) w
    real(dl), intent(in) :: tau
    real(dl) k,k2
    ! andrea
    real(dl) y(EV%nvar), yin(EV%nvar), yprime(EV%nvar), yprimein(EV%nvar)
    ! andrea
    real(dl) photbar,cs2,pb43,grho,slip,clxgdot, &
        clxcdot,clxbdot,adotdota,gpres,clxrdot,etak
    real(dl) q,aq,v
    real(dl) G11_t,G30_t, wnu_arr(max_nu)

    real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t,sigma,polter
    real(dl) w_dark_energy_t !equation of state of dark energy
    real(dl) gpres_noDE !Pressure with matter and radiation, no dark energy
    real(dl) qgdot,qrdot,pigdot,pirdot,vbdot,dgrho,adotoa
    real(dl) a,a2,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
    ! andrea delete dopacity
    real(dl) E2
    ! andrea
    integer l,i,ind, ind2, off_ix, ix
    real(dl) dgs,sigmadot,dz
    real(dl) dgpi,dgrho_matter,grho_matter, clxnu, gpres_nu
    !non-flat vars
    real(dl) cothxor !1/tau in flat case
    real(dl) xe,Trad, Delta_TM, Tmat, Delta_TCMB
    real(dl) delta_p_b, wing_t, wing2_t,winv_t
    real(dl) Delta_source2, polter_line
    real(dl) Delta_xe, Tspin, tau_eps, tau_fac, Tb
    integer lineoff,lineoffpol
    !Variables for source calculation
    real(dl) diff_rhopi, pidot_sum, dgpi_diff, phi
    real(dl) E(2:3), Edot(2:3)
    real(dl) phidot, polterdot, polterddot, octg, octgdot
    real(dl) lenswindow
    !integer, parameter :: num_cmb_freq = 6 !!!
    !integer, parameter :: nscatter = num_cmb_freq+1
    real(dl) :: opacity(nscatter), dopacity(nscatter), &
                            ddopacity(nscatter), exptau(nscatter), &
                            visibility(nscatter), dvisibility(nscatter), &
                            ddvisibility(nscatter)
    real(dl) ISW, quadrupole_source, doppler, monopole_source, tau0, ang_dist
    real(dl) dgrho_de, dgq_de, cs2_de
    ! andrea
    integer f_i, ix_off, s_ix, kays
    real(dl) drhophot,qgphot,piphot,opac_qgphot, opac_phot, octgprime
    real(dl) polter_freq, opac_rayleigh, opac_tot
    ! andrea

    !write(*,*) 'is rayleigh :', EV%Rayleigh
    k=EV%k_buf
    k2=EV%k2_buf
    
    !  Get background scale factor, sound speed and ionisation fraction.
    if (EV%TightCoupling) then
    ! and
        call this%Values_array(tau,a,cs2,opacity,dopacity)
    else
        call this%Values_array(tau,a,cs2,opacity)
    ! and
    end if
    a2=a*a

    etak=ay(ix_etak)

    !  CDM variables
    clxc=ay(ix_clxc)

    !  Baryon variables
    clxb=ay(ix_clxb)
    vb=ay(ix_vb)
    !  Compute expansion rate from: grho 8*pi*rho*a**2

    grhob_t=State%grhob/a
    grhoc_t=State%grhoc/a
    grhor_t=State%grhornomass/a2
    grhog_t=State%grhog/a2

    if (EV%is_cosmological_constant) then
        grhov_t = State%grhov * a2
        w_dark_energy_t = -1_dl
    else
        call State%CP%DarkEnergy%BackgroundDensityAndPressure(State%grhov, a, grhov_t, w_dark_energy_t)
    end if

    !total perturbations: matter terms first, then add massive nu, de and radiation
    !  8*pi*a*a*SUM[rho_i*clx_i]
    dgrho_matter=grhob_t*clxb+grhoc_t*clxc
    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=grhob_t*vb

    gpres_nu=0
    grhonu_t=0

    if (State%CP%Num_Nu_Massive > 0) then
        call MassiveNuVars(EV,ay,a,grhonu_t,gpres_nu,dgrho_matter,dgq, wnu_arr)
    end if

    grho_matter=grhonu_t+grhob_t+grhoc_t
    grho = grho_matter+grhor_t+grhog_t+grhov_t
    gpres_noDE = gpres_nu + (grhor_t + grhog_t)/3

    if (State%flat) then
        adotoa=sqrt(grho/3)
        cothxor=1._dl/tau
    else
        adotoa=sqrt((grho+State%grhok)/3._dl)
        cothxor=1._dl/State%tanfunc(tau/State%curvature_radius)/State%curvature_radius
    end if

    dgrho = dgrho_matter

    if (EV%no_nu_multpoles) then
        !RSA approximation of arXiv:1104.2933, dropping opactity terms in the velocity
        !Approximate total density variables with just matter terms
        z=(0.5_dl*dgrho/k + etak)/adotoa
        dz= -adotoa*z - 0.5_dl*dgrho/k
        clxr=-4*dz/k
        qr=-4._dl/3*z
        pir=0
    else
        !  Massless neutrinos
        clxr=ay(EV%r_ix)
        qr  =ay(EV%r_ix+1)
        pir =ay(EV%r_ix+2)
    endif

    pig=0
    if (EV%no_phot_multpoles) then
        if (.not. EV%no_nu_multpoles) then
            z=(0.5_dl*dgrho/k + etak)/adotoa
            dz= -adotoa*z - 0.5_dl*dgrho/k
            clxg=-4*dz/k-4/k*opacity(1)*(vb+z)
            qg=-4._dl/3*z
        else
            clxg=clxr-4/k*opacity(1)*(vb+z)
            qg=qr
        end if
    else
        !  Photons
        clxg=ay(EV%g_ix)
        qg=ay(EV%g_ix+1)
        if (.not. EV%TightCoupling) pig=ay(EV%g_ix+2)
    end if

    !  8*pi*a*a*SUM[rho_i*clx_i] - radiation terms
    dgrho=dgrho + grhog_t*clxg+grhor_t*clxr

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=dgq + grhog_t*qg+grhor_t*qr

    if (.not. EV%TightCoupling .and. .not.EV%no_phot_multpoles )  then
            call Phot_Integrate_L012(EV,ay,drhophot,qgphot,piphot, opacity, opac_qgphot,opac_phot)
            if (rayleigh_back) then
            dgrho=dgrho + grhog_t*drhophot
            dgq=dgq + grhog_t*qgphot
            else
            opac_qgphot = 0
            opac_phot = opacity(1)
            end if
        end if

    !  Photon mass density over baryon mass density
    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    if (.not. EV%is_cosmological_constant) then
        call State%CP%DarkEnergy%PerturbedStressEnergy(dgrho_de, dgq_de, &
            a, dgq, dgrho, grho, grhov_t, w_dark_energy_t, gpres_noDE, etak, &
            adotoa, k, EV%Kf(1), ay, ayprime, EV%w_ix)
        dgrho = dgrho + dgrho_de
        dgq = dgq + dgq_de
    end if

    !  Get sigma (shear) and z from the constraints
    ! have to get z from eta for numerical stability
    z=(0.5_dl*dgrho/k + etak)/adotoa
    if (State%flat) then
        !eta*k equation
        sigma=(z+1.5_dl*dgq/k2)
        ayprime(ix_etak)=0.5_dl*dgq
    else
        sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)
        ayprime(ix_etak)=0.5_dl*dgq + State%curv*z
    end if

    if (.not. EV%is_cosmological_constant) &
        call State%CP%DarkEnergy%PerturbationEvolve(ayprime, w_dark_energy_t, &
        EV%w_ix, a, adotoa, k, z, ay)

    !  CDM equation of motion
    clxcdot=-k*z
    ayprime(ix_clxc)=clxcdot

    !  Baryon equation of motion.
    clxbdot=-k*(z+vb)
    ayprime(ix_clxb)=clxbdot
    !  Photon equation of motion
    clxgdot=-k*(4._dl/3._dl*z+qg)

    !Sources
    if (EV%Evolve_baryon_cs) then
        if (a > State%CP%Recomb%min_a_evolve_Tm) then
            Tmat = State%CP%Recomb%T_m(a)
        else
            Tmat = State%CP%TCMB/a
        end if
        if (EV%Evolve_TM) then
            Delta_TM = ay(EV%Tg_ix)
        else
            Delta_TM = clxg/4
        end if
        delta_p_b = barssc0*(1._dl-0.75d0*State%CP%yhe+(1._dl-State%CP%yhe)*opacity(1)*a2/State%akthom)*Tmat*(clxb + delta_tm)
    else
        Delta_TM = clxg/4
        delta_p_b = cs2*clxb
    end if

    if (State%CP%Evolve_delta_xe) then
        if (EV%saha) then
            xe=State%CP%Recomb%x_e(a)
            Delta_xe = (1-xe)/(2-xe)*(-clxb + (3._dl/2+  CB1/Tmat)*Delta_TM)
        else
            Delta_xe = ay(EV%xe_ix)
        end if
    else
        Delta_xe = 0
    end if

    ! Easy to see instability in k \sim 1e-3 by tracking evolution of vb

    !  Use explicit equation for vb if appropriate

    if (EV%TightCoupling) then
        !  ddota/a
        gpres = gpres_noDE + w_dark_energy_t*grhov_t
        adotdota=(adotoa*adotoa-gpres)/2

        pig = 32._dl/45/opacity(1)*k*(sigma+vb)

        !  First-order approximation to baryon-photon splip
        slip = - (2*adotoa/(1+pb43) + dopacity(1)/opacity(1))* (vb-3._dl/4*qg) &
            +(-adotdota*vb-k/2*adotoa*clxg +k*(cs2*clxbdot-clxgdot/4))/(opacity(1)*(1+pb43))

        if (second_order_tightcoupling) then
            ! by Francis-Yan Cyr-Racine simplified (inconsistently) by AL assuming flat
            !AL: First order slip seems to be fine here to 2e-4

            !  8*pi*G*a*a*SUM[rho_i*sigma_i]
            dgs = grhog_t*pig+grhor_t*pir

            ! Define shear derivative to first order
            sigmadot = -2*adotoa*sigma-dgs/k+etak

            !Once know slip, recompute qgdot, pig, pigdot
            qgdot = k*(clxg/4._dl-pig/2._dl) +opacity(1)*slip

            pig = 32._dl/45/opacity(1)*k*(sigma+3._dl*qg/4._dl)*(1+(dopacity(1)*11._dl/6._dl/opacity(1)**2)) &
                + (32._dl/45._dl/opacity(1)**2)*k*(sigmadot+3._dl*qgdot/4._dl)*(-11._dl/6._dl)

            pigdot = -(32._dl/45._dl)*(dopacity(1)/opacity(1)**2)*k*(sigma+3._dl*qg/4._dl)*(1 + &
                dopacity(1)*11._dl/6._dl/opacity(1)**2 ) &
                + (32._dl/45._dl/opacity(1))*k*(sigmadot+3._dl*qgdot/4._dl)*(1+(11._dl/6._dl) &
                *(dopacity(1)/opacity(1)**2))

            EV%pigdot = pigdot

        end if

        !  Use tight-coupling approximation for vb
        !  zeroth order approximation to vbdot + the pig term
        vbdot=(-adotoa*vb+cs2*k*clxb + k/4*pb43*(clxg-2*EV%Kf(1)*pig))/(1+pb43)

        vbdot=vbdot+pb43/(1+pb43)*slip
        EV%pig = pig
    else
        if (.not. EV%TightCoupling .and. .not.EV%no_phot_multpoles)  then 
            vbdot = -adotoa*vb+cs2*k*clxb -photbar*(opac_phot*(4._dl/3*vb-qg) -opac_qgphot)
        else
            vbdot=-adotoa*vb+cs2*k*clxb-photbar*opacity(1)*(4._dl/3*vb-qg)
        end if
    end if

    ayprime(ix_vb)=vbdot

    if (.not. EV%no_phot_multpoles) then
        !  Photon equations of motion
        ayprime(EV%g_ix)=clxgdot
        !qgdot=4._dl/3*(-vbdot-adotoa*vb+delta_p_b*k)/pb43 &
        if (.not. EV%TightCoupling .and. .not.EV%no_phot_multpoles)  then 
        qgdot=4._dl/3*(photbar*opacity(1)*(4._dl/3*vb-qg))/pb43 &
            +EV%denlk(1)*clxg-EV%denlk2(1)*pig
        else
        qgdot=4._dl/3*(-vbdot-adotoa*vb+cs2*k*clxb)/pb43 &
            +EV%denlk(1)*clxg-EV%denlk2(1)*pig  
        end if
        ayprime(EV%g_ix+1)=qgdot

        !  Use explicit equations for photon moments if appropriate
        if (.not. EV%tightcoupling) then
            E2=ay(EV%polind+2)
            polter = pig/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)
            ix= EV%g_ix+2
            if (EV%lmaxg>2) then
                pigdot=EV%denlk(2)*qg-EV%denlk2(2)*ay(ix+1)-opacity(1)*(pig - polter) &
                    +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
                do  l=3,EV%lmaxg-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1)-EV%denlk2(l)*ay(ix+1))-opacity(1)*ay(ix)
                end do
                ix=ix+1
                !  Truncate the photon moment expansion
                ayprime(ix)=k*ay(ix-1)-(EV%lmaxg+1)*cothxor*ay(ix) -opacity(1)*ay(ix)
            else !closed case
                pigdot=EV%denlk(2)*qg-opacity(1)*(pig - polter) +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
            endif
            !  Polarization
            !l=2
            ix=EV%polind+2
            if (EV%lmaxgpol>2) then
                ayprime(ix) = -opacity(1)*(ay(ix) - polter) - k/3._dl*ay(ix+1)
                do l=3,EV%lmaxgpol-1
                    ix=ix+1
                    ayprime(ix)=-opacity(1)*ay(ix) + (EV%denlk(l)*ay(ix-1)-EV%polfack(l)*ay(ix+1))
                end do
                ix=ix+1
                !truncate
                ayprime(ix)=-opacity(1)*ay(ix) + &
                    k*EV%poltruncfac*ay(ix-1)-(EV%lmaxgpol+3)*cothxor*ay(ix)
            else !closed case
                ayprime(ix) = -opacity(1)*(ay(ix) - polter)
            endif
        end if
    end if

    if (EV%Rayleigh) then
         !assume after tight coupling ended
         do f_i = 1, num_cmb_freq
                !write(*,*) 'freq_factors = ', freq_factors
                opac_rayleigh = etat%Recombination_rayleigh_eff(a)*State%akthom/a2*(min(1._dl,&
                  freq_factors(f_i,1)/a2**2 + freq_factors(f_i,2)/a2**3  + freq_factors(f_i,3)/a2**4))
                opac_tot=opac_rayleigh+opacity(1)
                ind = EV%g_ix_freq + (f_i-1)*EV%freq_neq
                ayprime(ind)=-k*ay(ind+1)
                ayprime(ind+1)= EV%denlk(1)*ay(ind)-EV%denlk2(1)*ay(ind+2) + &
                    4._dl/3*photbar*(-opacity(1)*ay(ind+1) + opac_rayleigh*(4._dl/3*vb-qg-ay(ind+1)))/pb43
                E2=ay(EV%polind_freq + (f_i-1)*EV%freq_neq +2)
                polter_freq = ay(ind+2)/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)
                if (EV%lmaxg>2) then
                     ayprime(ind+2)=EV%denlk(2)*ay(ind+1)-EV%denlk2(2)*ay(ind+3)-opac_tot*(ay(ind+2) - polter_freq) &
                       -opac_rayleigh*(ay(EV%g_ix+2) - polter) 
                     ix= ind+2
                     do  l=3,EV%lmaxg-1
                        ix=ix+1
                        ayprime(ix)=(EV%denlk(l)*ay(ix-1)-EV%denlk2(l)*ay(ix+1))-opac_tot*ay(ix)-(opac_rayleigh)*ay(EV%g_ix+l)
                     end do
                     ix=ix+1
                  !  Truncate the photon moment expansion
                     ayprime(ix)=k*ay(ix-1)-(EV%lmaxg+1)*cothxor*ay(ix) -opac_tot*ay(ix)-opac_rayleigh*ay(EV%g_ix+EV%lmaxg)
                else 
                     ayprime(ind+2)=EV%denlk(2)*ay(ind+1)-opac_tot*(ay(ind+2) - polter_freq)-opac_rayleigh*(ay(EV%g_ix+2) - polter) 
                endif
        !  Polarization
                    !l=2
                    ix=EV%polind_freq+(f_i-1)*EV%freq_neq + 2
                    if (EV%lmaxgpol>2) then
                      ayprime(ix) = -opac_tot*(ay(ix) - polter_freq) - k/3._dl*ay(ix+1)- (opac_rayleigh)*(ay(EV%polind+2) - polter)
                      do l=3,EV%lmaxgpol-1
                       ix=ix+1
                       ayprime(ix)=-opac_tot*ay(ix)-opac_rayleigh*ay(EV%polind+l) + (EV%denlk(l)*ay(ix-1)-EV%polfack(l)*ay(ix+1))
                      end do
                      ix=ix+1
                      !truncate
                      ayprime(ix)=-opac_tot*ay(ix)-opac_rayleigh*ay(EV%polind+EV%lmaxgpol) + &
                        k*EV%poltruncfac*ay(ix-1)-(EV%lmaxgpol+3)*cothxor*ay(ix)
                   else !closed case
                      ayprime(ix) = -opac_tot*(ay(ix) - polter_freq)-opac_rayleigh*ay(EV%polind+2)
                   endif
         end do
     end if

    if (.not. EV%no_nu_multpoles) then
        !  Massless neutrino equations of motion.
        clxrdot=-k*(4._dl/3._dl*z+qr)
        ayprime(EV%r_ix)=clxrdot
        qrdot=EV%denlk(1)*clxr-EV%denlk2(1)*pir
        ayprime(EV%r_ix+1)=qrdot
        if (EV%high_ktau_neutrino_approx) then
            !ufa approximation for k*tau>>1, more accurate when there are reflections from lmax
            !Method from arXiv:1104.2933
            !                if (.not. EV%TightCoupling) then
            !                 gpres=gpres+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
            !                 adotdota=(adotoa*adotoa-gpres)/2
            !                end if
            !                ddz=(2*adotoa**2 - adotdota)*z  &
            !                  + adotoa/(2*k)*( 6*(grhog_t*clxg+grhor_t*clxr) + 2*(grhoc_t*clxc+grhob_t*clxb) ) &
            !                   - 1._dl/(2*k)*( 2*(grhog_t*clxgdot+grhor_t*clxrdot) + grhoc_t*clxcdot + grhob_t*clxbdot )
            !                dz= -adotoa*z - 0.5_dl*dgrho/k
            !                pirdot= -3*pir*cothxor + k*(qr+4._dl/3*z)
            pirdot= -3*pir*cothxor - clxrdot
            ayprime(EV%r_ix+2)=pirdot

            !                pirdot=k*(0.4_dl*qr-0.6_dl*ay(EV%lmaxg+10)+8._dl/15._dl*sigma)
            !                ayprime(EV%lmaxg+9)=pirdot
            !                ayprime(3+EV%lmaxg+7)=k*ay(3+EV%lmaxg+6)- &
            !                                      (3+1)*cothxor*ay(3+EV%lmaxg+7)
            !               ayprime(3+EV%lmaxg+7+1:EV%lmaxnr+EV%lmaxg+7)=0
        else
            ix=EV%r_ix+2
            if (EV%lmaxnr>2) then
                pirdot=EV%denlk(2)*qr- EV%denlk2(2)*ay(ix+1)+8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
                do l=3,EV%lmaxnr-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1) - EV%denlk2(l)*ay(ix+1))
                end do
                !  Truncate the neutrino expansion
                ix=ix+1
                ayprime(ix)=k*ay(ix-1)- (EV%lmaxnr+1)*cothxor*ay(ix)
            else
                pirdot=EV%denlk(2)*qr +8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
            end if
        end if
    end if ! no_nu_multpoles

    if (EV%Evolve_baryon_cs) then
        if (EV%Evolve_TM) then
            Delta_TCMB = clxg/4
            xe = State%CP%Recomb%x_e(a)
            Trad = State%CP%TCMB/a

            !Matter temperature
            !Recfast_CT = (8./3.)*(sigma_T/(m_e*C))*a_R in Mpc [a_R = radiation constant]
            ayprime(EV%Tg_ix) = -2*k*(z+vb)/3 - a*  Compton_CT * (Trad**4) * xe / (1._dl+xe+State%fHe) * &
                ((1- Trad/Tmat)*(Delta_TCMB*4 + Delta_xe/(1+xe/(1+State%fHe))) &
                + Trad/Tmat*(Delta_Tm - Delta_TCMB)  )

            if (State%CP%Evolve_delta_Ts) then
                ayprime(EV%Ts_ix) =  Get21cm_dTs(a,clxb,ay(EV%Ts_ix),Delta_TCMB,Delta_Tm,Tmat,Trad,xe )
            end if
        else
            if (State%CP%Evolve_delta_Ts) then
                ayprime(EV%Ts_ix) = -k*(4._dl/3._dl*z+qg)/4  !Assume follows Delta_TM which follows clxg
            end if
        end if
    end if

    if (State%CP%Evolve_delta_xe .and. .not. EV%saha) then
        ayprime(EV%xe_ix) = &
            State%CP%Recomb%dDeltaxe_dtau(a, Delta_xe,clxb, Delta_Tm, k*z/3,k*vb, adotoa)
    end if

    if (State%CP%Do21cm) then
        if (a > State%CP%Recomb%min_a_evolve_Tm) then
            if (State%CP%SourceTerms%line_reionization) then
                lineoff = EV%reion_line_ix+1
                lineoffpol = lineoff+EV%lmaxline-1

                if (tau> this%tau_start_redshiftwindows) then
                    !Multipoles of 21cm

                    polter_line = ay(lineoff+2)/10+9._dl/15*ay(lineoffpol+2)

                    call interp_window(State%TimeSteps,State%Redshift_W(1),tau,wing_t,wing2_t,winv_t)

                    delta_source2 = Get21cm_source2(a,clxb,Delta_TCMB,Delta_Tm,Delta_xe,Tmat,Trad,xe,k*(z+vb)/adotoa/3)


                    !Drop some small terms since mulipoles only enter into reionzation anyway
                    !monopole
                    ayprime(lineoff) = -k*ay(lineoff+1) +  wing_t * clxb + wing2_t*delta_source2 + k*z/3*winV_t


                    !dipole
                    ayprime(lineoff+1)= EV%denlk(1)*ay(lineoff)-EV%denlk2(1)*ay(lineoff+2) - opacity(1)*ay(lineoff+1) &
                        -wing2_t * ( qg/4 - vb/3)   ! vb/3*WinV_t)

                    !quadrupole
                    ayprime(lineoff+2)= EV%denlk(2)*ay(lineoff+1)-EV%denlk2(2)*ay(lineoff+3) &
                        +opacity(1)*(polter_line -ay(lineoff+2) ) -   2._dl/15*k*sigma*winV_t &
                        - wing2_t * ay(EV%g_ix+2)/4

                    do  l=3,EV%lmaxline-1
                        ayprime(lineoff+l)=EV%denlk(l)*ay(lineoff+l-1)-EV%denlk2(l)*ay(lineoff+l+1)-opacity(1)*ay(lineoff+l) &
                            - wing2_t * ay(EV%g_ix+l)/4
                    end do
                    !truncate
                    ayprime(lineoff+EV%lmaxline)=k*ay(lineoff+EV%lmaxline-1)-(EV%lmaxline+1)*cothxor*ay(lineoff+EV%lmaxline)  &
                        -opacity(1)*ay(lineoff+EV%lmaxline) - wing2_t * ay(EV%g_ix+EV%lmaxline)/4

                    !  21cm Polarization
                    !l=2
                    ayprime(lineoffpol+2) = -opacity(1)*(ay(lineoffpol+2) - polter_line) - k/3._dl*ay(lineoffpol+3)
                    !and the rest
                    do l=3,EV%lmaxline-1
                        ayprime(lineoffpol+l)=-opacity(1)*ay(lineoffpol+l) + EV%denlk(l)*ay(lineoffpol+l-1) -&
                            EV%polfack(l)*ay(lineoffpol+l+1)
                    end do

                    !truncate
                    ayprime(lineoffpol+EV%lmaxline)=-opacity(1)*ay(lineoffpol+EV%lmaxline) + &
                        k*EV%poltruncfac*ay(lineoffpol+EV%lmaxline-1)-(EV%lmaxline+3)*cothxor*ay(lineoffpol+EV%lmaxline)
                else
                    ayprime(lineoff:lineoffpol+EV%lmaxline)=0
                end if
            end if
        end if
    end if



    !  Massive neutrino equations of motion.
    if (State%CP%Num_Nu_massive >0) then
        !DIR$ LOOP COUNT MIN(1), AVG(1)
        do nu_i = 1, State%CP%Nu_mass_eigenstates
            if (EV%MassiveNuApprox(nu_i)) then
                !Now EV%iq0 = clx, EV%iq0+1 = clxp, EV%iq0+2 = G_1, EV%iq0+3=G_2=pinu
                !see astro-ph/0203507
                G11_t=EV%G11(nu_i)/a/a2
                G30_t=EV%G30(nu_i)/a/a2
                off_ix = EV%nu_ix(nu_i)
                w=wnu_arr(nu_i)
                ayprime(off_ix)=-k*z*(w+1) + 3*adotoa*(w*ay(off_ix) - ay(off_ix+1))-k*ay(off_ix+2)
                ayprime(off_ix+1)=(3*w-2)*adotoa*ay(off_ix+1) - 5._dl/3*k*z*w - k/3*G11_t
                ayprime(off_ix+2)=(3*w-1)*adotoa*ay(off_ix+2) - k*(2._dl/3*EV%Kf(1)*ay(off_ix+3)-ay(off_ix+1))
                ayprime(off_ix+3)=(3*w-2)*adotoa*ay(off_ix+3) + 2*w*k*sigma - k/5*(3*EV%Kf(2)*G30_t-2*G11_t)
            else
                ind=EV%nu_ix(nu_i)
                !DIR$ LOOP COUNT MIN(3), AVG(3)
                do i=1,EV%nq(nu_i)
                    q=State%NuPerturbations%nu_q(i)
                    aq=a*State%nu_masses(nu_i)/q
                    v=1._dl/sqrt(1._dl+aq*aq)

                    ayprime(ind)=-k*(4._dl/3._dl*z + v*ay(ind+1))
                    ind=ind+1
                    ayprime(ind)=v*(EV%denlk(1)*ay(ind-1)-EV%denlk2(1)*ay(ind+1))
                    ind=ind+1
                    if (EV%lmaxnu_tau(nu_i)==2) then
                        ayprime(ind)=-ayprime(ind-2) -3*cothxor*ay(ind)
                    else
                        ayprime(ind)=v*(EV%denlk(2)*ay(ind-1)-EV%denlk2(2)*ay(ind+1)) &
                            +k*8._dl/15._dl*sigma
                        do l=3,EV%lmaxnu_tau(nu_i)-1
                            ind=ind+1
                            ayprime(ind)=v*(EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
                        end do
                        !  Truncate moment expansion.
                        ind = ind+1
                        ayprime(ind)=k*v*ay(ind-1)-(EV%lmaxnu_tau(nu_i)+1)*cothxor*ay(ind)
                    end if
                    ind = ind+1
                end do
            end if
        end do

        if (EV%has_nu_relativistic) then
            ind=EV%nu_pert_ix
            ayprime(ind)=+k*a2*qr -k*ay(ind+1)
            ind2= EV%r_ix
            do l=1,EV%lmaxnu_pert-1
                ind=ind+1
                ind2=ind2+1
                ayprime(ind)= -a2*(EV%denlk(l)*ay(ind2-1)-EV%denlk2(l)*ay(ind2+1)) &
                    +   (EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
            end do
            ind=ind+1
            ind2=ind2+1
            ayprime(ind)= k*(ay(ind-1) -a2*ay(ind2-1)) -(EV%lmaxnu_pert+1)*cothxor*ay(ind)
        end if
    end if

    if (associated(EV%OutputTransfer) .or. associated(EV%OutputSources)) then
        if (EV%TightCoupling .or. EV%no_phot_multpoles) then
            E=0
            Edot=0
        else
            E = ay(EV%polind+2:EV%polind+3)
            Edot = ayprime(EV%polind+2:EV%polind+3)
        end if
        if (EV%no_nu_multpoles) then
            pirdot=0
            qrdot = -4*dz/3
        end if
        if (EV%no_phot_multpoles) then
            pigdot=0
            octg=0
            octgdot=0
            qgdot = -4*dz/3
        else
            if (EV%TightCoupling) then
                if (second_order_tightcoupling) then
                    octg = (3._dl/7._dl)*pig*(EV%k_buf/opacity(1))
                    E(2) = pig/4 + pigdot*(1._dl/opacity(1))*(-5._dl/8._dl)
                    E(3) = (3._dl/7._dl)*(EV%k_buf/opacity(1))*E(2)
                    Edot(2)= (pigdot/4._dl)*(1+(5._dl/2._dl)*(dopacity(1)/opacity(1)**2))
                else
                    pigdot = -dopacity(1)/opacity(1)*pig + 32._dl/45*k/opacity(1)*(-2*adotoa*sigma  &
                        +etak/EV%Kf(1)-  dgpi/k +vbdot )
                    Edot(2) = pigdot/4
                    E(2) = pig/4
                    octg=0
                end if
                octgdot=0
            else
                octg=ay(EV%g_ix+3)
                octgdot=ayprime(EV%g_ix+3)
            end if
        end if
        if (EV%is_cosmological_constant) then
            dgrho_de=0
            dgq_de=0
        end if

        dgpi  = grhor_t*pir + grhog_t*pig
        dgpi_diff = 0  !sum (3*p_nu -rho_nu)*pi_nu
        pidot_sum = grhog_t*pigdot + grhor_t*pirdot
        clxnu =0
        if (State%CP%Num_Nu_Massive /= 0) then
            call MassiveNuVarsOut(EV,ay,ayprime,a, adotoa, dgpi=dgpi, clxnu_all=clxnu, &
                dgpi_diff=dgpi_diff, pidot_sum=pidot_sum)
        end if
        gpres = gpres_noDE + w_dark_energy_t*grhov_t
        diff_rhopi = pidot_sum - (4*dgpi+ dgpi_diff)*adotoa + &
            State%CP%DarkEnergy%diff_rhopi_Add_Term(dgrho_de, dgq_de, grho, &
            gpres, w_dark_energy_t, State%grhok, adotoa, &
            EV%kf(1), k, grhov_t, z, k2, ayprime, ay, EV%w_ix)
        phi = -((dgrho +3*dgq*adotoa/k)/EV%Kf(1) + dgpi)/(2*k2)

        if (associated(EV%OutputTransfer)) then
            EV%OutputTransfer(Transfer_kh) = k/(State%CP%h0/100._dl)
            EV%OutputTransfer(Transfer_cdm) = clxc
            EV%OutputTransfer(Transfer_b) = clxb
            EV%OutputTransfer(Transfer_g) = clxg
            EV%OutputTransfer(Transfer_r) = clxr
            EV%OutputTransfer(Transfer_nu) = clxnu
            EV%OutputTransfer(Transfer_tot) =  dgrho_matter/grho_matter !includes neutrinos
            EV%OutputTransfer(Transfer_nonu) = (grhob_t*clxb+grhoc_t*clxc)/(grhob_t + grhoc_t)
            EV%OutputTransfer(Transfer_tot_de) =  dgrho/grho_matter
            !Transfer_Weyl is k^2Phi, where Phi is the Weyl potential
            EV%OutputTransfer(Transfer_Weyl) = k2*phi
            EV%OutputTransfer(Transfer_Newt_vel_cdm)=  -k*sigma/adotoa
            EV%OutputTransfer(Transfer_Newt_vel_baryon) = -k*(vb + sigma)/adotoa
            EV%OutputTransfer(Transfer_vel_baryon_cdm) = vb
            if (State%CP%do21cm) then
                Tspin = State%CP%Recomb%T_s(a)
                xe = State%CP%Recomb%x_e(a)

                tau_eps = a*line21_const*State%NNow/a**3/adotoa/Tspin/1000
                delta_source2 = Get21cm_source2(a,clxb,clxg/4,Delta_Tm,Delta_xe,Tmat,&
                    State%CP%TCMB/a,xe,k*(z+vb)/adotoa/3)
                tau_fac = tau_eps/(exp(tau_eps)-1)
                EV%OutputTransfer(Transfer_monopole) = ( clxb + Trad/(Tspin-Trad)*delta_source2 )  &
                    + (tau_fac-1)*(clxb - (delta_source2 + clxg/4)  )

                EV%OutputTransfer(Transfer_vnewt) = tau_fac*k*(vb+sigma)/adotoa
                EV%OutputTransfer(Transfer_Tmat) =  delta_TM
                if (State%CP%SourceTerms%use_21cm_mK) then
                    Tb = (1-exp(-tau_eps))*a*(Tspin-Trad)*1000

                    EV%OutputTransfer(Transfer_monopole) = EV%OutputTransfer(Transfer_monopole)*Tb
                    EV%OutputTransfer(Transfer_vnewt) = EV%OutputTransfer(Transfer_vnewt)*Tb
                    EV%OutputTransfer(Transfer_Tmat) = EV%OutputTransfer(Transfer_Tmat)*Tb
                end if
            end if
        end if

        if (associated(EV%OutputSources)) then
          
            if (State%num_redshiftwindows > 0) then
                call output_window_sources(EV, this, EV%OutputSources, ay, ayprime, &
                tau, a, adotoa, grho, gpres, &
                k, etak, z, ayprime(ix_etak), phi, phidot, sigma, sigmadot, &
                dgrho, clxg,clxb,clxc,clxnu, Delta_TM, Delta_xe, &
                dgq, qg, vb, qgdot, vbdot, &
                dgpi, pig, pigdot, diff_rhopi, &
                polter, polterdot, polterddot, octg, octgdot, E, Edot, &
                opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau)
            end if

            if (associated(EV%CustomSources)) then
                select type(DE=>State%CP%DarkEnergy)
                class is (TDarkEnergyEqnOfState)
                    cs2_de = DE%cs2_lam
                class default
                    cs2_de=1
                end select
                block
                    procedure(TSource_func), pointer :: custom_sources_func

                    call c_f_procpointer(CP%CustomSources%c_source_func,custom_sources_func)

                    call custom_sources_func(EV%CustomSources, tau, a, adotoa, grho, gpres,w_dark_energy_t, cs2_de, &
                        grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
                        k, etak, ayprime(ix_etak), phi, phidot, sigma, sigmadot, &
                        dgrho, clxg,clxb,clxc,clxr,clxnu, dgrho_de/grhov_t, delta_p_b, &
                        dgq, qg, qr, dgq_de/grhov_t, vb, qgdot, qrdot, vbdot, &
                        dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
                        polter, polterdot, polterddot, octg, octgdot, E, Edot, &
                        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
                        tau0, State%tau_maxvis, EV%Kf,f_K)
                end block
            end if

            EV%OutputSources = 0

            call IonizationFunctionsAtTime(this, tau, a, opacity, dopacity, ddopacity, &
                 visibility, dvisibility, ddvisibility, exptau, lenswindow)

            !do kays = 1, nscatter
            !    if (visibility(kays)  > 0.2) then
                    !write(*,*) 'visibility : ', vis(scat), scat
            !        write(*,*) 'visibility(kays) : ', visibility(kays)
            !    end if
            !end do

            tau0 = State%tau0
            phidot = (1.0d0/2.0d0)*(adotoa*(-dgpi - 2*k2*phi) + dgq*k - &
                     diff_rhopi+ k*sigma*(gpres + grho))/k2
            !time derivative of shear
            sigmadot = -adotoa*sigma - 1.0d0/2.0d0*dgpi/k + k*phi
            !quadrupole source derivatives; polter = pi_g/10 + 3/5 E_2
            polter = pig/10+9._dl/15*E(2)
            polterdot = (1.0d0/10.0d0)*pigdot + (3.0d0/5.0d0)*Edot(2)
            polterddot = -2.0d0/25.0d0*adotoa*dgq/(k*EV%Kf(1)) - 4.0d0/75.0d0*adotoa* &
                         k*sigma - 4.0d0/75.0d0*dgpi - 2.0d0/75.0d0*dgrho/EV%Kf(1) - 3.0d0/ &
                         50.0d0*k*octgdot*EV%Kf(2) + (1.0d0/25.0d0)*k*qgdot - 1.0d0/5.0d0 &
                         *k*EV%Kf(2)*Edot(3) + (-1.0d0/10.0d0*pig + (7.0d0/10.0d0)* &
                         polter - 3.0d0/5.0d0*E(2))*dopacity(1) + (-1.0d0/10.0d0*pigdot &
                         + (7.0d0/10.0d0)*polterdot - 3.0d0/5.0d0*Edot(2))*opacity(1)
            !Temperature source terms, after integrating by parts in conformal time

            ! and fix
            s_ix = 0
            !write(*,*) 'nscatter', nscatter
            do f_i = 1, nscatter
                if (.not. EV%no_phot_multpoles .and. EV%Rayleigh .and. f_i>1) then
                    ix_off=EV%g_ix_freq+(f_i-2)*EV%freq_neq 
                    y(EV%g_ix:EV%g_ix+EV%freq_neq-1) = yin(EV%g_ix:EV%g_ix+EV%freq_neq-1)+ &
                                                       yin(ix_off:ix_off+EV%freq_neq-1)
                    yprime(EV%g_ix:EV%g_ix+EV%freq_neq-1) = yprimein(EV%g_ix:EV%g_ix+EV%freq_neq-1)+ &
                                                            yprimein(ix_off:ix_off+EV%freq_neq-1)

                    pig = y(EV%g_ix+2)
                    polter = 0.1_dl*pig+9._dl/15._dl*E(2)
                    pigdot = yprime(EV%g_ix+2)
                    octg = y(EV%g_ix+3)
                    octgprime = yprime(EV%g_ix+3)
                    clxg = y(EV%g_ix)
                    qg  = y(EV%g_ix+1)
                    qgdot = yprime(EV%g_ix+1)
                    !exptau(f_i) = exptau(1)
                end if
                !2phi' term (\phi' + \psi' in Newtonian gauge), phi is the Weyl potential
                ISW = 2*phidot*exptau(f_i)
                monopole_source =  (-etak/(k*EV%Kf(1)) + 2*phi + clxg/4)*visibility(f_i)
                doppler = ((sigma + vb)*dvisibility(f_i) + (sigmadot + vbdot)*visibility(f_i))/k
                quadrupole_source = (5.0d0/8.0d0)*(3*polter*ddvisibility(f_i) + 6*polterdot*dvisibility(f_i) &
                                    + (k**2*polter + 3*polterddot)*visibility(f_i))/k**2

                EV%OutputSources(s_ix+1) = ISW + doppler + monopole_source + quadrupole_source
                ang_dist = f_K(tau0-tau)
                !if (EV%OutputSources(4) < -0.5) then
                    !write(*,*) 'ISW : ', ISW + doppler + monopole_source + quadrupole_source
                    !write(*,*) 'ISW', monopole_source, doppler, quadrupole_source
                    !write(*,*) 'visibility', visibility(f_i)
                !end if
                if (tau < tau0) then
                    !E polarization source
                    EV%OutputSources(s_ix+2)=visibility(1)*polter*(15._dl/8._dl)/(ang_dist**2*k2)
                    !factor of four because no 1/16 later
                end if

                !if (rayleigh_diff .and. f_i>2) then
                !     EV%OutputSources(s_ix+1:s_ix+2)= EV%OutputSources(s_ix+1:s_ix+2) - EV%OutputSources(1:2)
                !end if

                s_ix = s_ix+2
            
                if (size(EV%OutputSources) > 2 .and. f_i==1) then
                    ! and
                    !Get lensing sources
                    if (tau>State%tau_maxvis .and. tau0-tau > 0.1_dl) then
                        EV%OutputSources(3) = -2*phi*f_K(tau-State%tau_maxvis)/(f_K(tau0-State%tau_maxvis)*ang_dist)
                        !We include the lensing factor of two here
                    end if
                    s_ix=s_ix+1
                end if
            end do !f_1
        !    and
        !if (EV%OutputSources(1) < -0.5) then
        !    EV%OutputSources(1) = 0.0005
        !end if
        !if (EV%OutputSources(2) < -0.5) then
        !    EV%OutputSources(2) = 0.0005
        !end if
        !if (EV%OutputSources(3) < -0.5) then
        !    EV%OutputSources(3) = 0.0005
        !end if
        !if (EV%OutputSources(4) < -0.5) then
        !    EV%OutputSources(4) = 0.0005
        !end if
        !if (EV%OutputSources(5) < -0.5) then
        !    EV%OutputSources(5) = 0.0005
        !end if
        !if (EV%OutputSources(6) < -0.5) then
        !    EV%OutputSources(6) = 0.0005
        !end if
        !if (EV%OutputSources(7) < -0.5) then
        !    EV%OutputSources(7) = 0.0005
        !end if
        !if (EV%OutputSources(8) < -0.5) then
        !    EV%OutputSources(8) = 0.0005
        !end if
        !if (EV%OutputSources(9) < -0.5) then
        !    EV%OutputSources(9) = 0.0005
        !end if
        !if (EV%OutputSources(10) < -0.5) then
        !    EV%OutputSources(10) = 0.0005
        !end if
        !if (EV%OutputSources(11) < -0.5) then
        !    EV%OutputSources(11) = 0.0005
        !end if
        !if (EV%OutputSources(12) < -0.5) then
        !    EV%OutputSources(12) = 0.0005
        !end if
        !if (EV%OutputSources(13) < -0.5) then
        !    EV%OutputSources(13) = 0.0005
        !end if
        !if (EV%OutputSources(14) < -0.5) then
        !    EV%OutputSources(14) = 0.0005
        !end if
        !if (EV%OutputSources(15) < -0.5) then
        !    EV%OutputSources(15) = 0.0005
        !end if
        
            !write(*,*) 'EV%OutputSources(4)', EV%OutputSources(4)
            !write(*,*) 'test', (tau>State%tau_maxvis .and. tau0-tau > 0.1_dl)
            !write(*,*) 'EV%OutputSources(1)', EV%OutputSources(1)
            !write(*,*) 'EV%OutputSources(2)', EV%OutputSources(2)
            !write(*,*) 'EV%OutputSources(3)', EV%OutputSources(3)
        !end if
        !if (tau < 210) then
        !write(*,*) EV%OutputSources(12)
        !write(*,*) '1', ISW + doppler + monopole_source + quadrupole_source
        !write(*,*) '2', visibility(1)*polter*(15._dl/8._dl)/(ang_dist**2*k2)
        !write(*,*) '3', 2*phi*f_K(tau-State%tau_maxvis)/(f_K(tau0-State%tau_maxvis)*ang_dist)
        !write(*,*) '1', visibility(:)
        !end if
        end if
    end if

    !if (EV%OutputSources(4) < -1.0) 
    !write(*,*) 'EV%OutputSources', EV%OutputSources
    !write(*,*) 'EV%OutputSources(4)', EV%OutputSources(4)
    
    end subroutine derivs

    subroutine derivsv(EV,this,etat,n,tau,yv,yvprime)
    !  Evaluate the time derivatives of the vector perturbations, flat case
    use MassiveNu
    implicit none
    type(EvolutionVars) EV
    class(TThermoData) :: this
    class(TRecfast) :: etat
    integer n,l
    real(dl), target ::  yv(n),yvprime(n)
    real(dl) ep,tau,grho,rhopi,cs2,opacity(nscatter),gpres
    logical finished_tightcoupling
    real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
    real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
    real(dl) sigma, qg,pig, qr, vb, rhoq, vbdot, photbar, pb43
    real(dl) k,k2,a,a2, adotdota
    real(dl) pir,adotoa
    real(dl) w_dark_energy_t

    k2=EV%k2_buf
    k=EV%k_buf

    !E and B start at l=2. Set up pointers accordingly to fill in y arrays
    E => yv(EV%lmaxv+3:)
    Eprime=> yvprime(EV%lmaxv+3:)
    B => E(EV%lmaxpolv:)
    Bprime => Eprime(EV%lmaxpolv:)
    neutprime => Bprime(EV%lmaxpolv+1:)
    neut => B(EV%lmaxpolv+1:)

    a=yv(1)

    sigma=yv(2)


    !  Get sound speed and opacity, and see if should use tight-coupling
    
    ! and modification DU CAMB INITIAL HORS RAYLEIGH
    call this%values(tau,a,cs2,opacity)
    ! and
    ! andrea not sure
    !call this%values_array(tau,a,cs2,opacity)
    ! andrea
    if (k > 0.06_dl*epsw) then
        ep=ep0
    else
        ep=0.2_dl*ep0
    end if
    a2=a*a

    finished_tightcoupling = &
        ((k/opacity(1) > ep).or.(1._dl/(opacity(1)*tau) > ep .and. k/opacity(1) > 1d-4))


    ! Compute expansion rate from: grho=8*pi*rho*a**2
    ! Also calculate gpres: 8*pi*p*a**2
    grhob_t=State%grhob/a
    grhoc_t=State%grhoc/a
    grhor_t=State%grhornomass/a2
    grhog_t=State%grhog/a2
    call CP%DarkEnergy%BackgroundDensityAndPressure(State%grhov, a, grhov_t, w_dark_energy_t)

    grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
    gpres=(grhog_t+grhor_t)/3._dl+grhov_t*w_dark_energy_t

    adotoa=sqrt(grho/3._dl)
    adotdota=(adotoa*adotoa-gpres)/2

    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    yvprime(1)=adotoa*a

    vb = yv(3)
    qg = yv(4)
    qr = neut(1)

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    rhoq=grhob_t*vb+grhog_t*qg+grhor_t*qr
    !  sigma = 2*rhoq/k**2
    !for non-large k this expression for sigma is unstable at early times
    !so propagate sigma equation separately (near total cancellation in rhoq)

    if (finished_tightcoupling) then
        !  Use explicit equations:

        pig = yv(5)

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

        vbdot = -adotoa*vb-photbar*opacity(1)*(4._dl/3*vb-qg) - 0.5_dl*k*photbar*Magnetic

        !  Equation for the photon heat flux stress

        yvprime(4)=-0.5_dl*k*pig + opacity(1)*(4._dl/3*vb-qg)

        !  Equation for the photon anisotropic stress
        yvprime(5)=k*(2._dl/5*qg -8/15._dl*yv(6))+8._dl/15._dl*k*sigma  &
            -opacity(1)*(pig - polter)
        ! And for the moments
        do  l=3,EV%lmaxv-1
            yvprime(l+3)=k*denl(l)*l*(yv(l+2)-   &
                vecfac(l)*yv(l+4))-opacity(1)*yv(l+3)
        end do
        !  Truncate the hierarchy
        yvprime(EV%lmaxv+3)=k*EV%lmaxv/(EV%lmaxv-1._dl)*yv(EV%lmaxv+2)- &
            (EV%lmaxv+2._dl)*yv(EV%lmaxv+3)/tau-opacity(1)*yv(EV%lmaxv+3)

        !E equations

        Eprime(2) = - opacity(1)*(E(2) - polter) + k*(1/3._dl*B(2) - 8._dl/27._dl*E(3))
        do l=3,EV%lmaxpolv-1
            Eprime(l) =-opacity(1)*E(l) + k*(denl(l)*(l*E(l-1) - &
                vecfacpol(l)*E(l+1)) + 2._dl/(l*(l+1))*B(l))
        end do
        !truncate
        Eprime(EV%lmaxpolv)=0._dl

        !B-bar equations

        do l=2,EV%lmaxpolv-1
            Bprime(l) =-opacity(1)*B(l) + k*(denl(l)*(l*B(l-1) - &
                vecfacpol(l)*B(l+1)) - 2._dl/(l*(l+1))*E(l))
        end do
        !truncate
        Bprime(EV%lmaxpolv)=0._dl
    else
        !Tight coupling expansion results

        pig = 32._dl/45._dl*k/opacity(1)*(vb + sigma)

        EV%pig = pig

        vbdot=(-adotoa*vb  -3._dl/8*pb43*k*Magnetic  -3._dl/8*k*pb43*pig &
            - pb43/(1+pb43)/opacity(1)*(0.75_dl*k*adotoa*pb43**2/(pb43+1)*Magnetic + vb*&
            ( 2*pb43*adotoa**2/(1+pb43) + adotdota)) &
            )/(1+pb43)

        !  Equation for the photon heat flux
        ! Get drag from vbdot expression
        yvprime(4)=-0.5_dl*k*pig - &
            (vbdot+adotoa*vb)/photbar - 0.5_dl*k*Magnetic

        !  Set the derivatives to zero
        yvprime(5:n)=0._dl
        yv(5)=pig
        E(2)=  pig/4
    endif

    yvprime(3) = vbdot

    !  Neutrino equations:

    !  Massless neutrino anisotropic stress
    pir=neut(2)
    neutprime(1)= -0.5_dl*k*pir
    neutprime(2)=2._dl/5*k*qr -8._dl/15._dl*k*neut(3)+ 8._dl/15._dl*k*sigma
    !  And for the moments
    do  l=3,EV%lmaxnrv-1
        neutprime(l)=k*denl(l)*l*(neut(l-1)- vecfac(l)*neut(l+1))
    end do

    !  Truncate the hierarchy
    neutprime(EV%lmaxnrv)=k*EV%lmaxnrv/(EV%lmaxnrv-1._dl)*neut(EV%lmaxnrv-1)-  &
        (EV%lmaxnrv+2._dl)*neut(EV%lmaxnrv)/tau


    !  Get the propagation equation for the shear

    rhopi=grhog_t*pig+grhor_t*pir+ grhog_t*Magnetic

    yvprime(2)=-2*adotoa*sigma -rhopi/k

    end subroutine derivsv


    subroutine derivst(EV, this, etat, n,tau,ayt,aytprime)
    !  Evaluate the time derivatives of the tensor perturbations.
    use MassiveNu
    implicit none
    type(EvolutionVars) EV
    ! and au pif
    class(TThermoData) this
    class(TRecfast) etat
    ! and
    integer n,l,i,ind, nu_i
    real(dl), target ::  ayt(n),aytprime(n)
    ! andrea
    real(dl) opac_rayleigh,opac_tot,polter_freq
    integer f_i
    integer, parameter :: num_cmb_freq = 6 !!!
    integer, parameter :: nscatter = num_cmb_freq+1
    real(dl) tau, rhopi, opacity(nscatter), pirdt
    real(dl), dimension(:),pointer :: Ef,Bf,Efprime,Bfprime
    ! andrea
    real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
    real(dl) q,aq,v
    real(dl) Hchi,pinu, pig, polter
    real(dl) k,k2,a,a2,grhog_t, grhor_t
    real(dl) pir, adot, adotoa, rhonu, shear
    real(dl) cothxor
    k2=EV%k2_buf
    k= EV%k_buf

    call this%Expansion_Values(tau,a, adot,opacity)

    Hchi=ayt(ixt_H)
    shear=ayt(ixt_shear)
    a2=a*a
    adotoa = adot/a

    if (State%flat) then
        cothxor=1._dl/tau
    else
        cothxor=1._dl/State%tanfunc(tau/State%curvature_radius)/State%curvature_radius
    end if

    if (.not. EV%TensTightCoupling) then
        !  Don't use tight coupling approx - use explicit equations:
        !  Equation for the photon anisotropic stress


        !E and B start at l=2. Set up pointers accordingly to fill in ayt arrays
        E => ayt(EV%E_ix+1:)
        B => ayt(EV%B_ix+1:)
        Eprime=> aytprime(EV%E_ix+1:)
        Bprime => aytprime(EV%B_ix+1:)

        ind = EV%g_ix+2

        !  Photon anisotropic stress
        pig=ayt(ind)
        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

        if (EV%lmaxt > 2) then
            ! andrea
            aytprime(ind)=-EV%denlkt(2,2)*ayt(ind+1)+k*8._dl/15._dl*shear  &
                -opacity(1)*(pig - polter)
            ! andrea

            do l=3, EV%lmaxt -1
                ind = ind+1
                ! andrea
                aytprime(ind)=EV%denlkt(1,L)*ayt(ind-1)-EV%denlkt(2,L)*ayt(ind+1)-opacity(1)*ayt(ind)
                ! andrea
            end do

            !Truncate the hierarchy
            ind=ind+1
            ! andrea
            aytprime(ind)=k*EV%lmaxt/(EV%lmaxt-2._dl)*ayt(ind-1)- &
                (EV%lmaxt+3._dl)*cothxor*ayt(ind)-opacity(1)*ayt(ind)

            !E and B-bar equations

            Eprime(2) = - opacity(1)*(E(2) - polter) + EV%denlkt(4,2)*B(2) - &
                EV%denlkt(3,2)*E(3)
            do l=3, EV%lmaxpolt-1
                Eprime(l) =(EV%denlkt(1,L)*E(l-1)-EV%denlkt(3,L)*E(l+1) + EV%denlkt(4,L)*B(l)) &
                    -opacity(1)*E(l)
            end do
            l= EV%lmaxpolt
            !truncate: difficult, but setting l+1 to zero seems to work OK
            Eprime(l) = (EV%denlkt(1,L)*E(l-1) + EV%denlkt(4,L)*B(l)) -opacity(1)*E(l)

            Bprime(2) =-EV%denlkt(3,2)*B(3) - EV%denlkt(4,2)*E(2)  -opacity(1)*B(2)
            do l=3, EV%lmaxpolt-1
                Bprime(l) =(EV%denlkt(1,L)*B(l-1) -EV%denlkt(3,L)*B(l+1) - EV%denlkt(4,L)*E(l)) &
                    -opacity(1)*B(l)
            end do
            l=EV%lmaxpolt
            !truncate
            Bprime(l) =(EV%denlkt(1,L)*B(l-1) - EV%denlkt(4,L)*E(l))  -opacity(1)*B(l)
            if (EV%Rayleigh) then
            !assume after tight coupling ended
                do f_i = 1, num_cmb_freq
                    opac_rayleigh = etat%Recombination_rayleigh_eff(a)*State%akthom/a2*(min(1._dl,&
                    freq_factors(f_i,1)/a2**2 + freq_factors(f_i,2)/a2**3  + freq_factors(f_i,3)/a2**4)) 
                    opac_tot=opac_rayleigh+opacity(1)
                    ind = EV%g_ix_freq + (f_i-1)*EV%freq_neq
                    Ef => ayt(ind + (EV%lmaxt-1)+1:)
                    Bf => ayt(ind + (EV%lmaxt-1)+(EV%lmaxpolt-1)+1:)
                    Efprime=> aytprime(ind + (EV%lmaxt-1)+1:)
                    Bfprime => aytprime(ind + (EV%lmaxt-1)+(EV%lmaxpolt-1)+1:)
                    ind =ind + 2
                    !  Photon anisotropic stres
                    polter_freq = 0.1_dl*ayt(ind) + 9._dl/15._dl*Ef(2)
                    aytprime(ind)=-EV%denlkt(2,2)*ayt(ind+1) -opac_rayleigh*(pig - polter) &
                    - opac_tot*(ayt(ind)-polter_freq)
                    do l=3, EV%lmaxt -1
                        ind = ind+1
                        aytprime(ind)=EV%denlkt(1,L)*ayt(ind-1)-EV%denlkt(2,L)*ayt(ind+1)-opac_rayleigh*ayt(EV%g_ix+l) &
                        - opac_tot*ayt(ind)
                    end do
                    !Truncate the hierarchy
                    ind=ind+1
                    aytprime(ind)=k*EV%lmaxt/(EV%lmaxt-2._dl)*ayt(ind-1)- &
                    (EV%lmaxt+3._dl)*cothxor*ayt(ind)-opac_tot*ayt(ind)-opac_rayleigh*ayt(EV%g_ix+EV%lmaxt)
                    !E and B-bar equations
                    Efprime(2) = - opac_rayleigh*(E(2) - polter) - opac_tot*(Ef(2)-polter_freq) + EV%denlkt(4,2)*Bf(2) - &
                    EV%denlkt(3,2)*Ef(3)
                    do l=3, EV%lmaxpolt-1
                        Efprime(l) =(EV%denlkt(1,L)*Ef(l-1)-EV%denlkt(3,L)*Ef(l+1) + EV%denlkt(4,L)*Bf(l)) &
                        -opac_rayleigh*E(l) - opac_tot*Ef(l)
                    end do
                    l= EV%lmaxpolt
                    !truncate: difficult, but setting l+1 to zero seems to work OK
                    Efprime(l) = (EV%denlkt(1,L)*Ef(l-1) + EV%denlkt(4,L)*Bf(l)) -opac_rayleigh*E(l) - opac_tot*Ef(l)
                    Bfprime(2) =-EV%denlkt(3,2)*Bf(3) - EV%denlkt(4,2)*Ef(2) -opac_rayleigh*B(2) -opac_tot*Bf(2)
                    do l=3, EV%lmaxpolt-1
                        Bfprime(l) =(EV%denlkt(1,L)*Bf(l-1) -EV%denlkt(3,L)*Bf(l+1) - EV%denlkt(4,L)*Ef(l)) &
                        -opac_rayleigh*B(l) - opac_tot*Bf(l)
                    end do
                    l=EV%lmaxpolt
                    !truncate
                    Bfprime(l) =(EV%denlkt(1,L)*Bf(l-1) - EV%denlkt(4,L)*Ef(l))  -opac_rayleigh*B(l) - opac_tot*Bf(l)
                end do
            end if
        else !lmax=2
            ! andrea pas sûr de l'indentation ...)'
            aytprime(ind)=k*8._dl/15._dl*shear-opacity(1)*(pig - polter) 
            Eprime(2) = - opacity(1)*(E(2) - polter) + EV%denlkt(4,2)*B(2)
            Bprime(2) = - EV%denlkt(4,2)*E(2)  -opacity(1)*B(2)
            ! andrea
        end if
    else  !Tight coupling
        pig = 32._dl/45._dl*k/opacity(1)*shear
    endif
    grhor_t=State%grhornomass/a2
    grhog_t=State%grhog/a2
    rhopi=grhog_t*pig
    !  Neutrino equations:
    !  Anisotropic stress
    if (DoTensorNeutrinos) then
        neutprime => aytprime(EV%r_ix+1:)
        neut => ayt(EV%r_ix+1:)
        !  Massless neutrino anisotropic stress
        pir=neut(2)
        rhopi=rhopi+grhor_t*pir
        if (EV%lmaxnrt>2) then
            pirdt=-EV%denlkt(2,2)*neut(3) + 8._dl/15._dl*k*shear
            neutprime(2)=pirdt
            !  And for the moments
            do  l=3, EV%lmaxnrt-1
                neutprime(l)= EV%denlkt(1,L)*neut(l-1) -EV%denlkt(2,L)*neut(l+1)
            end do
            !  Truncate the hierarchy
            neutprime(EV%lmaxnrt)=k*EV%lmaxnrt/(EV%lmaxnrt-2._dl)*neut(EV%lmaxnrt-1)-  &
                (EV%lmaxnrt+3._dl)*cothxor*neut(EV%lmaxnrt)
        else
            pirdt= 8._dl/15._dl*k*shear
            neutprime(2)=pirdt
        end if

        !  Massive neutrino equations of motion and contributions to anisotropic stress.
        if (State%CP%Num_Nu_massive > 0) then
            do nu_i=1,State%CP%Nu_mass_eigenstates
                if (.not. EV%EvolveTensorMassiveNu(nu_i)) then
                    rhopi=rhopi+ State%grhormass(nu_i)/a2*pir !- good approx, note no rhonu weighting
                else
                    ind=EV%nu_ix(nu_i)+2

                    pinu= Nu_pi(EV, ayt, a, nu_i)
                    rhopi=rhopi+ State%grhormass(nu_i)/a2*pinu

                    do i=1,State%NuPerturbations%nqmax
                        q=State%NuPerturbations%nu_q(i)
                        aq=a*State%nu_masses(nu_i)/q
                        v=1._dl/sqrt(1._dl+aq*aq)
                        if (EV%lmaxnut>2) then
                            aytprime(ind)=-v*EV%denlkt(2,2)*ayt(ind+1)+8._dl/15._dl*k*shear
                            do l=3,EV%lmaxnut-1
                                ind=ind+1
                                aytprime(ind)=v*(EV%denlkt(1,L)*ayt(ind-1)-EV%denlkt(2,L)*ayt(ind+1))
                            end do
                            ind = ind+1
                            !Truncate moment expansion.
                            aytprime(ind)=k*v*EV%lmaxnut/(EV%lmaxnut-2._dl)*ayt(ind-1)-(EV%lmaxnut+3)*cothxor*ayt(ind)
                        else
                            aytprime(ind)=8._dl/15._dl*k*shear
                        end if
                        ind=ind+1
                    end do
                end if
            end do
        end if
    end if

    !  Get the propagation equation for the shear

    if (State%flat) then
        aytprime(ixt_shear)=-2*adotoa*shear+k*Hchi-rhopi/k
    else
        aytprime(ixt_shear)=-2*adotoa*shear+k*Hchi*(1+2*State%curv/k2)-rhopi/k
    endif

    aytprime(ixt_H)=-k*shear

    end subroutine derivst

    subroutine dverk_derivs(EV, this, etat, n, x, y, xend, tol, ind, c, nw, w)
    use Precision
    use MpiUtils
    use Config, only : GlobalError, error_evolution
    integer n, ind, nw, k
    real(dl) x, y(n), xend, tol, c(*), w(nw,9), temp
    type(EvolutionVars) EV
    class(TThermoData) :: this
    class(TRecfast) :: etat
     !it isn't, but as long as it maintains it as a pointer we are OK
    !
    !***********************************************************************
    !                                                                      *
    ! note added 11/14/85.                                                 *
    !                                                                      *
    ! if you discover any errors in this subroutine, please contact        *
    !                                                                      *
    !        kenneth r. jackson                                            *
    !        department of computer science                                *
    !        university of toronto                                         *
    !        toronto, ontario,                                             *
    !        canada   m5s 1a4                                              *
    !                                                                      *
    !        phone: 416-978-7075                                           *
    !                                                                      *
    !        electronic mail:                                              *
    !        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
    !        csnet:  krj@toronto                                           *
    !        arpa:   krj.toronto@csnet-relay                               *
    !        bitnet: krj%toronto@csnet-relay.arpa                          *
    !                                                                      *
    ! dverk is written in fortran 66.                                      *
    !                                                                      *
    ! the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
    ! set for a  vax  in  double  precision.  they  should  be  reset,  as *
    ! described below, if this program is run on another machine.          *
    !                                                                      *
    ! the c array is declared in this subroutine to have one element only, *
    ! although  more  elements  are  referenced  in this subroutine.  this *
    ! causes some compilers to issue warning messages.  there is,  though, *
    ! no  error  provided  c is declared sufficiently large in the calling *
    ! program, as described below.                                         *
    !                                                                      *
    ! the following external statement  for  fcn  was  added  to  avoid  a *
    ! warning  message  from  the  unix  f77 compiler.  the original dverk *
    ! comments and code follow it.                                         *
    !                                                                      *
    !***********************************************************************
    !
    !external fcn
    !
    !***********************************************************************
    !                                                                      *
    !     purpose - this is a runge-kutta  subroutine  based  on  verner's *
    ! fifth and sixth order pair of formulas for finding approximations to *
    ! the solution of  a  system  of  first  order  ordinary  differential *
    ! equations  with  initial  conditions. it attempts to keep the global *
    ! error proportional to  a  tolerance  specified  by  the  user.  (the *
    ! proportionality  depends  on the kind of error control that is used, *
    ! as well as the differential equation and the range of integration.)  *
    !                                                                      *
    !     various options are available to the user,  including  different *
    ! kinds  of  error control, restrictions on step sizes, and interrupts *
    ! which permit the user to examine the state of the  calculation  (and *
    ! perhaps make modifications) during intermediate stages.              *
    !                                                                      *
    !     the program is efficient for non-stiff systems.  however, a good *
    ! variable-order-adams  method  will probably be more efficient if the *
    ! function evaluations are very costly.  such a method would  also  be *
    ! more suitable if one wanted to obtain a large number of intermediate *
    ! solution values by interpolation, as might be the case  for  example *
    ! with graphical output.                                               *
    !                                                                      *
    !                                    hull-enright-jackson   1/10/76    *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !     use - the user must specify each of the following                *
    !                                                                      *
    !     n  number of equations                                           *
    !                                                                      *
    !   fcn  name of subroutine for evaluating functions - the  subroutine *
    !           itself must also be provided by the user - it should be of *
    !           the following form                                         *
    !              subroutine fcn(n, x, y, yprime)                         *
    !              integer n                                               *
    !              real(dl) x, y(n), yprime(n)                     *
    !                      *** etc ***                                     *
    !           and it should evaluate yprime, given n, x and y            *
    !                                                                      *
    !     x  independent variable - initial value supplied by user         *
    !                                                                      *
    !     y  dependent variable - initial values of components y(1), y(2), *
    !           ..., y(n) supplied by user                                 *
    !                                                                      *
    !  xend  value of x to which integration is to be carried out - it may *
    !           be less than the initial value of x                        *
    !                                                                      *
    !   tol  tolerance - the subroutine attempts to control a norm of  the *
    !           local  error  in  such  a  way  that  the  global error is *
    !           proportional to tol. in some problems there will be enough *
    !           damping  of  errors, as well as some cancellation, so that *
    !           the global error will be less than tol. alternatively, the *
    !           control   can   be  viewed  as  attempting  to  provide  a *
    !           calculated value of y at xend which is the exact  solution *
    !           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
    !           is proportional to tol.  (the norm  is  a  max  norm  with *
    !           weights  that  depend on the error control strategy chosen *
    !           by the user.  the default weight for the k-th component is *
    !           1/max(1,abs(y(k))),  which therefore provides a mixture of *
    !           absolute and relative error control.)                      *
    !                                                                      *
    !   ind  indicator - on initial entry ind must be set equal to  either *
    !           1  or  2. if the user does not wish to use any options, he *
    !           should set ind to 1 - all that remains for the user to  do *
    !           then  is  to  declare c and w, and to specify nw. the user *
    !           may also  select  various  options  on  initial  entry  by *
    !           setting ind = 2 and initializing the first 9 components of *
    !           c as described in the next section.  he may also  re-enter *
    !           the  subroutine  with ind = 3 as mentioned again below. in *
    !           any event, the subroutine returns with ind equal to        *
    !              3 after a normal return                                 *
    !              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
    !              -1, -2, or -3 after an error condition (see below)      *
    !                                                                      *
    !     c  communications vector - the dimension must be greater than or *
    !           equal to 24, unless option c(1) = 4 or 5 is used, in which *
    !           case the dimension must be greater than or equal to n+30   *
    !                                                                      *
    !    nw  first dimension of workspace w -  must  be  greater  than  or *
    !           equal to n                                                 *
    !                                                                      *
    !     w  workspace matrix - first dimension must be nw and second must *
    !           be greater than or equal to 9                              *
    !                                                                      *
    !     the subroutine  will  normally  return  with  ind  =  3,  having *
    ! replaced the initial values of x and y with, respectively, the value *
    ! of xend and an approximation to y at xend.  the  subroutine  can  be *
    ! called  repeatedly  with new values of xend without having to change *
    ! any other argument.  however, changes in tol, or any of the  options *
    ! described below, may also be made on such a re-entry if desired.     *
    !                                                                      *
    !     three error returns are also possible, in which  case  x  and  y *
    ! will be the most recently accepted values -                          *
    !     with ind = -3 the subroutine was unable  to  satisfy  the  error *
    !        requirement  with a particular step-size that is less than or *
    !        equal to hmin, which may mean that tol is too small           *
    !     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
    !        probably  means  that the requested tol (which is used in the *
    !        calculation of hmin) is too small                             *
    !     with ind = -1 the allowed maximum number of fcn evaluations  has *
    !        been  exceeded,  but  this  can only occur if option c(7), as *
    !        described in the next section, has been used                  *
    !                                                                      *
    !     there are several circumstances that will cause the calculations *
    ! to  be  terminated,  along with output of information that will help *
    ! the user determine the cause of  the  trouble.  these  circumstances *
    ! involve  entry with illegal or inconsistent values of the arguments, *
    ! such as attempting a normal  re-entry  without  first  changing  the *
    ! value of xend, or attempting to re-enter with ind less than zero.    *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !     options - if the subroutine is entered with ind = 1, the first 9 *
    ! components of the communications vector are initialized to zero, and *
    ! the subroutine uses only default values  for  each  option.  if  the *
    ! subroutine  is  entered  with ind = 2, the user must specify each of *
    ! these 9 components - normally he would first set them all  to  zero, *
    ! and  then  make  non-zero  those  that  correspond to the particular *
    ! options he wishes to select. in any event, options may be changed on *
    ! re-entry  to  the  subroutine  -  but if the user changes any of the *
    ! options, or tol, in the course of a calculation he should be careful *
    ! about  how  such changes affect the subroutine - it may be better to *
    ! restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
    ! program  -  the information is available to the user, but should not *
    ! normally be changed by him.)                                         *
    !                                                                      *
    !  c(1)  error control indicator - the norm of the local error is  the *
    !           max  norm  of  the  weighted  error  estimate  vector, the *
    !           weights being determined according to the value of c(1) -  *
    !              if c(1)=1 the weights are 1 (absolute error control)    *
    !              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
    !                 control)                                             *
    !              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
    !                 (relative  error  control,  unless abs(y(k)) is less *
    !                 than the floor value, abs(c(2)) )                    *
    !              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
    !                 (here individual floor values are used)              *
    !              if c(1)=5 the weights are 1/abs(c(k+30))                *
    !              for all other values of c(1), including  c(1) = 0,  the *
    !                 default  values  of  the  weights  are  taken  to be *
    !                 1/max(1,abs(y(k))), as mentioned earlier             *
    !           (in the two cases c(1) = 4 or 5 the user must declare  the *
    !           dimension of c to be at least n+30 and must initialize the *
    !           components c(31), c(32), ..., c(n+30).)                    *
    !                                                                      *
    !  c(2)  floor value - used when the indicator c(1) has the value 3    *
    !                                                                      *
    !  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
    !           to be abs(c(3)) - otherwise it uses the default value      *
    !              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
    !           where dwarf is a very small positive  machine  number  and *
    !           rreb is the relative roundoff error bound                  *
    !                                                                      *
    !  c(4)  hstart specification - if not zero, the subroutine  will  use *
    !           an  initial  hmag equal to abs(c(4)), except of course for *
    !           the restrictions imposed by hmin and hmax  -  otherwise it *
    !           uses the default value of hmax*(tol)**(1/6)                *
    !                                                                      *
    !  c(5)  scale specification - this is intended to be a measure of the *
    !           scale of the problem - larger values of scale tend to make *
    !           the method more reliable, first  by  possibly  restricting *
    !           hmax  (as  described  below) and second, by tightening the *
    !           acceptance requirement - if c(5) is zero, a default  value *
    !           of  1  is  used.  for  linear  homogeneous  problems  with *
    !           constant coefficients, an appropriate value for scale is a *
    !           norm  of  the  associated  matrix.  for other problems, an *
    !           approximation to  an  average  value  of  a  norm  of  the *
    !           jacobian along the trajectory may be appropriate           *
    !                                                                      *
    !  c(6)  hmax specification - four cases are possible                  *
    !           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
    !              min(abs(c(6)),2/abs(c(5)))                              *
    !           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
    !           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
    !              2/abs(c(5))                                             *
    !           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
    !              of 2                                                    *
    !                                                                      *
    !  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
    !           error  return with ind = -1 will be caused when the number *
    !           of function evaluations exceeds abs(c(7))                  *
    !                                                                      *
    !  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
    !           interrupt   the  calculations  after  it  has  chosen  its *
    !           preliminary value of hmag, and just before choosing htrial *
    !           and  xtrial  in  preparation for taking a step (htrial may *
    !           differ from hmag in sign, and may  require  adjustment  if *
    !           xend  is  near) - the subroutine returns with ind = 4, and *
    !           will resume calculation at the point  of  interruption  if *
    !           re-entered with ind = 4                                    *
    !                                                                      *
    !  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
    !           interrupt   the  calculations  immediately  after  it  has *
    !           decided whether or not to accept the result  of  the  most *
    !           recent  trial step, with ind = 5 if it plans to accept, or *
    !           ind = 6 if it plans to reject -  y(*)  is  the  previously *
    !           accepted  result, while w(*,9) is the newly computed trial *
    !           value, and w(*,2) is the unweighted error estimate vector. *
    !           the  subroutine  will  resume calculations at the point of *
    !           interruption on re-entry with ind = 5 or 6. (the user  may *
    !           change ind in this case if he wishes, for example to force *
    !           acceptance of a step that would otherwise be rejected,  or *
    !           vice versa. he can also restart with ind = 1 or 2.)        *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !  summary of the components of the communications vector              *
    !                                                                      *
    !     prescribed at the option       determined by the program         *
    !           of the user                                                *
    !                                                                      *
    !                                    c(10) rreb(rel roundoff err bnd)  *
    !     c(1) error control indicator   c(11) dwarf (very small mach no)  *
    !     c(2) floor value               c(12) weighted norm y             *
    !     c(3) hmin specification        c(13) hmin                        *
    !     c(4) hstart specification      c(14) hmag                        *
    !     c(5) scale specification       c(15) scale                       *
    !     c(6) hmax specification        c(16) hmax                        *
    !     c(7) max no of fcn evals       c(17) xtrial                      *
    !     c(8) interrupt no 1            c(18) htrial                      *
    !     c(9) interrupt no 2            c(19) est                         *
    !                                    c(20) previous xend               *
    !                                    c(21) flag for xend               *
    !                                    c(22) no of successful steps      *
    !                                    c(23) no of successive failures   *
    !                                    c(24) no of fcn evals             *
    !                                                                      *
    !  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !  an overview of the program                                          *
    !                                                                      *
    !     begin initialization, parameter checking, interrupt re-entries   *
    !  ......abort if ind out of range 1 to 6                              *
    !  .     cases - initial entry, normal re-entry, interrupt re-entries  *
    !  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
    !  v........abort if n.gt.nw or tol.le.0                               *
    !  .        if initial entry without options (ind .eq. 1)              *
    !  .           set c(1) to c(9) equal to zero                          *
    !  .        else initial entry with options (ind .eq. 2)               *
    !  .           make c(1) to c(9) non-negative                          *
    !  .           make floor values non-negative if they are to be used   *
    !  .        end if                                                     *
    !  .        initialize rreb, dwarf, prev xend, flag, counts            *
    !  .     case 2 - normal re-entry (ind .eq. 3)                         *
    !  .........abort if xend reached, and either x changed or xend not    *
    !  .        re-initialize flag                                         *
    !  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
    !  v        transfer control to the appropriate re-entry point.......  *
    !  .     end cases                                                  .  *
    !  .  end initialization, etc.                                      .  *
    !  .                                                                v  *
    !  .  loop through the following 4 stages, once for each trial step .  *
    !  .     stage 1 - prepare                                          .  *
    !***********error return (with ind=-1) if no of fcn evals too great .  *
    !  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
    !  .        calc hmin, scale, hmax                                  .  *
    !***********error return (with ind=-2) if hmin .gt. hmax            .  *
    !  .        calc preliminary hmag                                   .  *
    !***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
    !  .        calc hmag, xtrial and htrial                            .  *
    !  .     end stage 1                                                .  *
    !  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
    !  .     stage 3 - calc the error estimate                          .  *
    !  .     stage 4 - make decisions                                   .  *
    !  .        set ind=5 if step acceptable, else set ind=6            .  *
    !***********interrupt no 2 if requested....................re-entry.v  *
    !  .        if step accepted (ind .eq. 5)                              *
    !  .           update x, y from xtrial, ytrial                         *
    !  .           add 1 to no of successful steps                         *
    !  .           set no of successive failures to zero                   *
    !**************return(with ind=3, xend saved, flag set) if x .eq. xend *
    !  .        else step not accepted (ind .eq. 6)                        *
    !  .           add 1 to no of successive failures                      *
    !**************error return (with ind=-3) if hmag .le. hmin            *
    !  .        end if                                                     *
    !  .     end stage 4                                                   *
    !  .  end loop                                                         *
    !  .                                                                   *
    !  begin abort action                                                  *
    !     output appropriate  message  about  stopping  the  calculations, *
    !        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
    !        previous xend,  no of  successful  steps,  no  of  successive *
    !        failures, no of fcn evals, and the components of y            *
    !     stop                                                             *
    !  end abort action                                                    *
    !                                                                      *
    !***********************************************************************
    !
    !     ******************************************************************
    !     * begin initialization, parameter checking, interrupt re-entries *
    !     ******************************************************************
    !
    !  ......abort if ind out of range 1 to 6
    if (ind.lt.1 .or. ind.gt.6) go to 500
    !
    !        cases - initial entry, normal re-entry, interrupt re-entries
    !         go to (5, 5, 45, 1111, 2222, 2222), ind
    if (ind==3) goto 45
    if (ind==4) goto 1111
    if (ind==5 .or. ind==6) goto 2222

    !        case 1 - initial entry (ind .eq. 1 or 2)
    !  .........abort if n.gt.nw or tol.le.0
    if (n.gt.nw .or. tol.le.0._dl) go to 500
    if (ind.eq. 2) go to 15
    !              initial entry without options (ind .eq. 1)
    !              set c(1) to c(9) equal to 0
    do k = 1, 9
        c(k) = 0._dl
    end do
    go to 35
15  continue
    !              initial entry with options (ind .eq. 2)
    !              make c(1) to c(9) non-negative
    do k = 1, 9
        c(k) = dabs(c(k))
    end do
    !              make floor values non-negative if they are to be used
    if (c(1).ne.4._dl .and. c(1).ne.5._dl) go to 30
    do k = 1, n
        c(k+30) = dabs(c(k+30))
    end do
30  continue
35  continue
    !           initialize rreb, dwarf, prev xend, flag, counts
    c(10) = 2._dl**(-56)
    c(11) = 1.d-35
    !           set previous xend initially to initial value of x
    c(20) = x
    do k = 21, 24
        c(k) = 0._dl
    end do
    go to 50
    !        case 2 - normal re-entry (ind .eq. 3)
    !  .........abort if xend reached, and either x changed or xend not
45  if (c(21).ne.0._dl .and. &
        (x.ne.c(20) .or. xend.eq.c(20))) go to 500
    !           re-initialize flag
    c(21) = 0._dl
    go to 50
    !        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
    !           transfer control to the appropriate re-entry point..........
    !           this has already been handled by the computed go to        .
    !        end cases                                                     v
50  continue
    !
    !     end initialization, etc.
    !
    !     ******************************************************************
    !     * loop through the following 4 stages, once for each trial  step *
    !     * until the occurrence of one of the following                   *
    !     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
    !     *        stage 4                                                 *
    !     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
    !     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
    !     *        requested, in stage 1 or stage 4                        *
    !     ******************************************************************
    !
99999 continue
    !
    !        ***************************************************************
    !        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
    !        * and some parameter  checking,  and  end  up  with  suitable *
    !        * values of hmag, xtrial and htrial in preparation for taking *
    !        * an integration step.                                        *
    !        ***************************************************************
    !
    !***********error return (with ind=-1) if no of fcn evals too great
    if (c(7).eq.0._dl .or. c(24).lt.c(7)) go to 100
    ind = -1
    return
100 continue
    !
    !           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    if (ind .eq. 6) go to 105
    call derivs(EV, this, etat, n, x, y, w(1,1))
    c(24) = c(24) + 1._dl
105 continue
    !
    !           calculate hmin - use default unless value prescribed
    c(13) = c(3)
    if (c(3) .ne. 0._dl) go to 165
    !              calculate default value of hmin
    !              first calculate weighted norm y - c(12) - as specified
    !              by the error control indicator c(1)
    temp = 0._dl
    if (c(1) .ne. 1._dl) go to 115
    !                 absolute error control - weights are 1
    do 110 k = 1, n
        temp = dmax1(temp, dabs(y(k)))
110 continue
    c(12) = temp
    go to 160
115 if (c(1) .ne. 2._dl) go to 120
    !                 relative error control - weights are 1/dabs(y(k)) so
    !                 weighted norm y is 1
    c(12) = 1._dl
    go to 160
120 if (c(1) .ne. 3._dl) go to 130
    !                 weights are 1/max(c(2),abs(y(k)))
    do 125 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(2))
125 continue
    c(12) = dmin1(temp, 1._dl)
    go to 160
130 if (c(1) .ne. 4._dl) go to 140
    !                 weights are 1/max(c(k+30),abs(y(k)))
    do 135 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(k+30))
135 continue
    c(12) = dmin1(temp, 1._dl)
    go to 160
140 if (c(1) .ne. 5._dl) go to 150
    !                 weights are 1/c(k+30)
    do 145 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(k+30))
145 continue
    c(12) = temp
    go to 160
150 continue
    !                 default case - weights are 1/max(1,abs(y(k)))
    do 155 k = 1, n
        temp = dmax1(temp, dabs(y(k)))
155 continue
    c(12) = dmin1(temp, 1._dl)
160 continue
    c(13) = 10._dl*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
    !
    !           calculate scale - use default unless value prescribed
    c(15) = c(5)
    if (c(5) .eq. 0._dl) c(15) = 1._dl
    !
    !           calculate hmax - consider 4 cases
    !           case 1 both hmax and scale prescribed
    if (c(6).ne.0._dl .and. c(5).ne.0._dl) &
        c(16) = dmin1(c(6), 2._dl/c(5))
    !           case 2 - hmax prescribed, but scale not
    if (c(6).ne.0._dl .and. c(5).eq.0._dl) c(16) = c(6)
    !           case 3 - hmax not prescribed, but scale is
    if (c(6).eq.0._dl .and. c(5).ne.0._dl) c(16) = 2._dl/c(5)
    !           case 4 - neither hmax nor scale is provided
    if (c(6).eq.0._dl .and. c(5).eq.0._dl) c(16) = 2._dl
    !
    !***********error return (with ind=-2) if hmin .gt. hmax
    if (c(13) .le. c(16)) go to 170
    ind = -2
    return
170 continue
    !
    !           calculate preliminary hmag - consider 3 cases
    if (ind .gt. 2) go to 175
    !           case 1 - initial entry - use prescribed value of hstart, if
    !              any, else default
    c(14) = c(4)
    if (c(4) .eq. 0._dl) c(14) = c(16)*tol**(1._dl/6._dl)
    go to 185
175 if (c(23) .gt. 1._dl) go to 180
    !           case 2 - after a successful step, or at most  one  failure,
    !              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
    !              overflow. then avoid reduction by more than half.
    temp = 2._dl*c(14)
    if (tol .lt. (2._dl/.9d0)**6*c(19)) &
        temp = .9d0*(tol/c(19))**(1._dl/6._dl)*c(14)
    c(14) = dmax1(temp, .5d0*c(14))
    go to 185
180 continue
    !           case 3 - after two or more successive failures
    c(14) = .5d0*c(14)
185 continue
    !
    !           check against hmax
    c(14) = dmin1(c(14), c(16))
    !
    !           check against hmin
    c(14) = dmax1(c(14), c(13))
    !
    !***********interrupt no 1 (with ind=4) if requested
    if (c(8) .eq. 0._dl) go to 1111
    ind = 4
    return
    !           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
    !
    !           calculate hmag, xtrial - depending on preliminary hmag, xend
    if (c(14) .ge. dabs(xend - x)) go to 190
    !              do not step more than half way to xend
    c(14) = dmin1(c(14), .5d0*dabs(xend - x))
    c(17) = x + dsign(c(14), xend - x)
    go to 195
190 continue
    !              hit xend exactly
    c(14) = dabs(xend - x)
    c(17) = xend
195 continue
    !
    !           calculate htrial
    c(18) = c(17) - x
    !
    !        end stage 1
    !
    !        ***************************************************************
    !        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
    !        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
    !        * stage 3. w(*,9) is temporary storage until finally it holds *
    !        * ytrial.                                                     *
    !        ***************************************************************
    !
    temp = c(18)/1398169080000._dl
    !
    do 200 k = 1, n
        w(k,9) = y(k) + temp*w(k,1)*233028180000._dl
200 continue
    call derivs(EV, this, etat, n, x + c(18)/6._dl, w(1,9), w(1,2))
    !
    do 205 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*74569017600._dl &
            + w(k,2)*298276070400._dl  )
205 continue
    call derivs(EV, this, etat, n, x + c(18)*(4._dl/15._dl), w(1,9), w(1,3))
    !
    do 210 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*1165140900000._dl &
            - w(k,2)*3728450880000._dl &
            + w(k,3)*3495422700000._dl )
210 continue
    call derivs(EV, this, etat, n, x + c(18)*(2._dl/3._dl), w(1,9), w(1,4))
    !
    do 215 k = 1, n
        w(k,9) = y(k) + temp*( - w(k,1)*3604654659375._dl &
            + w(k,2)*12816549900000._dl &
            - w(k,3)*9284716546875._dl &
            + w(k,4)*1237962206250._dl )
215 continue
    call derivs(EV, this, etat, n, x + c(18)*(5._dl/6._dl), w(1,9), w(1,5))
    !
    do 220 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*3355605792000._dl &
            - w(k,2)*11185352640000._dl &
            + w(k,3)*9172628850000._dl &
            - w(k,4)*427218330000._dl &
            + w(k,5)*482505408000._dl  )
220 continue
    call derivs(EV, this, etat, n, x + c(18), w(1,9), w(1,6))
    !
    do 225 k = 1, n
        w(k,9) = y(k) + temp*( - w(k,1)*770204740536._dl &
            + w(k,2)*2311639545600._dl &
            - w(k,3)*1322092233000._dl &
            - w(k,4)*453006781920._dl &
            + w(k,5)*326875481856._dl  )
225 continue
    call derivs(EV, this, etat, n, x + c(18)/15._dl, w(1,9), w(1,7))
    !
    do 230 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*2845924389000._dl &
            - w(k,2)*9754668000000._dl &
            + w(k,3)*7897110375000._dl &
            - w(k,4)*192082660000._dl &
            + w(k,5)*400298976000._dl &
            + w(k,7)*201586000000._dl  )
230 continue
    call derivs(EV, this, etat, n, x + c(18), w(1,9), w(1,8))
    !
    !           calculate ytrial, the extrapolated approximation and store
    !              in w(*,9)
    do 235 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*104862681000._dl &
            + w(k,3)*545186250000._dl &
            + w(k,4)*446637345000._dl &
            + w(k,5)*188806464000._dl &
            + w(k,7)*15076875000._dl &
            + w(k,8)*97599465000._dl   )
235 continue
    !
    !           add 7 to the no of fcn evals
    c(24) = c(24) + 7._dl
    !
    !        end stage 2
    !
    !        ***************************************************************
    !        * stage 3 - calculate the error estimate est. first calculate *
    !        * the  unweighted  absolute  error  estimate vector (per unit *
    !        * step) for the unextrapolated approximation and store it  in *
    !        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
    !        * specified by the error  control  indicator  c(1).  finally, *
    !        * modify  this result to produce est, the error estimate (per *
    !        * unit step) for the extrapolated approximation ytrial.       *
    !        ***************************************************************
    !
    !           calculate the unweighted absolute error estimate vector
    do 300 k = 1, n
        w(k,2) = (   w(k,1)*8738556750._dl &
            + w(k,3)*9735468750._dl &
            - w(k,4)*9709507500._dl &
            + w(k,5)*8582112000._dl &
            + w(k,6)*95329710000._dl &
            - w(k,7)*15076875000._dl &
            - w(k,8)*97599465000._dl)/1398169080000._dl
300 continue
    !
    !           calculate the weighted max norm of w(*,2) as specified by
    !           the error control indicator c(1)
    temp = 0._dl
    if (c(1) .ne. 1._dl) go to 310
    !              absolute error control
    do 305 k = 1, n
        temp = dmax1(temp,dabs(w(k,2)))
305 continue
    go to 360
310 if (c(1) .ne. 2._dl) go to 320
    !              relative error control
    do 315 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)/y(k)))
315 continue
    go to 360
320 if (c(1) .ne. 3._dl) go to 330
    !              weights are 1/max(c(2),abs(y(k)))
    do 325 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(c(2), dabs(y(k))) )
325 continue
    go to 360
330 if (c(1) .ne. 4._dl) go to 340
    !              weights are 1/max(c(k+30),abs(y(k)))
    do 335 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(c(k+30), dabs(y(k))) )
335 continue
    go to 360
340 if (c(1) .ne. 5._dl) go to 350
    !              weights are 1/c(k+30)
    do 345 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
345 continue
    go to 360
350 continue
    !              default case - weights are 1/max(1,abs(y(k)))
    do 355 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(1._dl, dabs(y(k))) )
355 continue
360 continue
    !
    !           calculate est - (the weighted max norm of w(*,2))*hmag*scale
    !              - est is intended to be a measure of the error  per  unit
    !              step in ytrial
    c(19) = temp*c(14)*c(15)
    !
    !        end stage 3
    !
    !        ***************************************************************
    !        * stage 4 - make decisions.                                   *
    !        ***************************************************************
    !
    !           set ind=5 if step acceptable, else set ind=6
    ind = 5
    if (c(19) .gt. tol) ind = 6
    !
    !***********interrupt no 2 if requested
    if (c(9) .eq. 0._dl) go to 2222
    return
    !           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
2222 continue
    !
    if (ind .eq. 6) go to 410
    !              step accepted (ind .eq. 5), so update x, y from xtrial,
    !                 ytrial, add 1 to the no of successful steps, and set
    !                 the no of successive failures to zero
    x = c(17)
    do 400 k = 1, n
        y(k) = w(k,9)
400 continue
    c(22) = c(22) + 1._dl
    c(23) = 0._dl
    !**************return(with ind=3, xend saved, flag set) if x .eq. xend
    if (x .ne. xend) go to 405
    ind = 3
    c(20) = xend
    c(21) = 1._dl
    return
405 continue
    go to 420
410 continue
    !              step not accepted (ind .eq. 6), so add 1 to the no of
    !                 successive failures
    c(23) = c(23) + 1._dl
    !**************error return (with ind=-3) if hmag .le. hmin
    if (c(14) .gt. c(13)) go to 415
    ind = -3
    return
415 continue
420 continue
    !
    !        end stage 4
    !
    go to 99999
    !     end loop
    !
    !  begin abort action
500 continue
    !

    write (*,*) 'Error in dverk_derivs, x =',x, 'xend=', xend
    call GlobalError('DVERK_DERIVS error', error_evolution)
    end subroutine dverk_derivs

    subroutine dverk_derivst(EV, this, etat, n, x, y, xend, tol, ind, c, nw, w)
    use Precision
    use MpiUtils
    use Config, only : GlobalError, error_evolution
    integer n, ind, nw, k
    real(dl) x, y(n), xend, tol, c(*), w(nw,9), temp
    type(EvolutionVars) EV
    class(TThermoData) :: this
    class(TRecfast) :: etat
     !it isn't, but as long as it maintains it as a pointer we are OK
    !
    !***********************************************************************
    !                                                                      *
    ! note added 11/14/85.                                                 *
    !                                                                      *
    ! if you discover any errors in this subroutine, please contact        *
    !                                                                      *
    !        kenneth r. jackson                                            *
    !        department of computer science                                *
    !        university of toronto                                         *
    !        toronto, ontario,                                             *
    !        canada   m5s 1a4                                              *
    !                                                                      *
    !        phone: 416-978-7075                                           *
    !                                                                      *
    !        electronic mail:                                              *
    !        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
    !        csnet:  krj@toronto                                           *
    !        arpa:   krj.toronto@csnet-relay                               *
    !        bitnet: krj%toronto@csnet-relay.arpa                          *
    !                                                                      *
    ! dverk is written in fortran 66.                                      *
    !                                                                      *
    ! the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
    ! set for a  vax  in  double  precision.  they  should  be  reset,  as *
    ! described below, if this program is run on another machine.          *
    !                                                                      *
    ! the c array is declared in this subroutine to have one element only, *
    ! although  more  elements  are  referenced  in this subroutine.  this *
    ! causes some compilers to issue warning messages.  there is,  though, *
    ! no  error  provided  c is declared sufficiently large in the calling *
    ! program, as described below.                                         *
    !                                                                      *
    ! the following external statement  for  fcn  was  added  to  avoid  a *
    ! warning  message  from  the  unix  f77 compiler.  the original dverk *
    ! comments and code follow it.                                         *
    !                                                                      *
    !***********************************************************************
    !
    !external fcn
    !
    !***********************************************************************
    !                                                                      *
    !     purpose - this is a runge-kutta  subroutine  based  on  verner's *
    ! fifth and sixth order pair of formulas for finding approximations to *
    ! the solution of  a  system  of  first  order  ordinary  differential *
    ! equations  with  initial  conditions. it attempts to keep the global *
    ! error proportional to  a  tolerance  specified  by  the  user.  (the *
    ! proportionality  depends  on the kind of error control that is used, *
    ! as well as the differential equation and the range of integration.)  *
    !                                                                      *
    !     various options are available to the user,  including  different *
    ! kinds  of  error control, restrictions on step sizes, and interrupts *
    ! which permit the user to examine the state of the  calculation  (and *
    ! perhaps make modifications) during intermediate stages.              *
    !                                                                      *
    !     the program is efficient for non-stiff systems.  however, a good *
    ! variable-order-adams  method  will probably be more efficient if the *
    ! function evaluations are very costly.  such a method would  also  be *
    ! more suitable if one wanted to obtain a large number of intermediate *
    ! solution values by interpolation, as might be the case  for  example *
    ! with graphical output.                                               *
    !                                                                      *
    !                                    hull-enright-jackson   1/10/76    *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !     use - the user must specify each of the following                *
    !                                                                      *
    !     n  number of equations                                           *
    !                                                                      *
    !   fcn  name of subroutine for evaluating functions - the  subroutine *
    !           itself must also be provided by the user - it should be of *
    !           the following form                                         *
    !              subroutine fcn(n, x, y, yprime)                         *
    !              integer n                                               *
    !              real(dl) x, y(n), yprime(n)                     *
    !                      *** etc ***                                     *
    !           and it should evaluate yprime, given n, x and y            *
    !                                                                      *
    !     x  independent variable - initial value supplied by user         *
    !                                                                      *
    !     y  dependent variable - initial values of components y(1), y(2), *
    !           ..., y(n) supplied by user                                 *
    !                                                                      *
    !  xend  value of x to which integration is to be carried out - it may *
    !           be less than the initial value of x                        *
    !                                                                      *
    !   tol  tolerance - the subroutine attempts to control a norm of  the *
    !           local  error  in  such  a  way  that  the  global error is *
    !           proportional to tol. in some problems there will be enough *
    !           damping  of  errors, as well as some cancellation, so that *
    !           the global error will be less than tol. alternatively, the *
    !           control   can   be  viewed  as  attempting  to  provide  a *
    !           calculated value of y at xend which is the exact  solution *
    !           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
    !           is proportional to tol.  (the norm  is  a  max  norm  with *
    !           weights  that  depend on the error control strategy chosen *
    !           by the user.  the default weight for the k-th component is *
    !           1/max(1,abs(y(k))),  which therefore provides a mixture of *
    !           absolute and relative error control.)                      *
    !                                                                      *
    !   ind  indicator - on initial entry ind must be set equal to  either *
    !           1  or  2. if the user does not wish to use any options, he *
    !           should set ind to 1 - all that remains for the user to  do *
    !           then  is  to  declare c and w, and to specify nw. the user *
    !           may also  select  various  options  on  initial  entry  by *
    !           setting ind = 2 and initializing the first 9 components of *
    !           c as described in the next section.  he may also  re-enter *
    !           the  subroutine  with ind = 3 as mentioned again below. in *
    !           any event, the subroutine returns with ind equal to        *
    !              3 after a normal return                                 *
    !              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
    !              -1, -2, or -3 after an error condition (see below)      *
    !                                                                      *
    !     c  communications vector - the dimension must be greater than or *
    !           equal to 24, unless option c(1) = 4 or 5 is used, in which *
    !           case the dimension must be greater than or equal to n+30   *
    !                                                                      *
    !    nw  first dimension of workspace w -  must  be  greater  than  or *
    !           equal to n                                                 *
    !                                                                      *
    !     w  workspace matrix - first dimension must be nw and second must *
    !           be greater than or equal to 9                              *
    !                                                                      *
    !     the subroutine  will  normally  return  with  ind  =  3,  having *
    ! replaced the initial values of x and y with, respectively, the value *
    ! of xend and an approximation to y at xend.  the  subroutine  can  be *
    ! called  repeatedly  with new values of xend without having to change *
    ! any other argument.  however, changes in tol, or any of the  options *
    ! described below, may also be made on such a re-entry if desired.     *
    !                                                                      *
    !     three error returns are also possible, in which  case  x  and  y *
    ! will be the most recently accepted values -                          *
    !     with ind = -3 the subroutine was unable  to  satisfy  the  error *
    !        requirement  with a particular step-size that is less than or *
    !        equal to hmin, which may mean that tol is too small           *
    !     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
    !        probably  means  that the requested tol (which is used in the *
    !        calculation of hmin) is too small                             *
    !     with ind = -1 the allowed maximum number of fcn evaluations  has *
    !        been  exceeded,  but  this  can only occur if option c(7), as *
    !        described in the next section, has been used                  *
    !                                                                      *
    !     there are several circumstances that will cause the calculations *
    ! to  be  terminated,  along with output of information that will help *
    ! the user determine the cause of  the  trouble.  these  circumstances *
    ! involve  entry with illegal or inconsistent values of the arguments, *
    ! such as attempting a normal  re-entry  without  first  changing  the *
    ! value of xend, or attempting to re-enter with ind less than zero.    *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !     options - if the subroutine is entered with ind = 1, the first 9 *
    ! components of the communications vector are initialized to zero, and *
    ! the subroutine uses only default values  for  each  option.  if  the *
    ! subroutine  is  entered  with ind = 2, the user must specify each of *
    ! these 9 components - normally he would first set them all  to  zero, *
    ! and  then  make  non-zero  those  that  correspond to the particular *
    ! options he wishes to select. in any event, options may be changed on *
    ! re-entry  to  the  subroutine  -  but if the user changes any of the *
    ! options, or tol, in the course of a calculation he should be careful *
    ! about  how  such changes affect the subroutine - it may be better to *
    ! restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
    ! program  -  the information is available to the user, but should not *
    ! normally be changed by him.)                                         *
    !                                                                      *
    !  c(1)  error control indicator - the norm of the local error is  the *
    !           max  norm  of  the  weighted  error  estimate  vector, the *
    !           weights being determined according to the value of c(1) -  *
    !              if c(1)=1 the weights are 1 (absolute error control)    *
    !              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
    !                 control)                                             *
    !              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
    !                 (relative  error  control,  unless abs(y(k)) is less *
    !                 than the floor value, abs(c(2)) )                    *
    !              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
    !                 (here individual floor values are used)              *
    !              if c(1)=5 the weights are 1/abs(c(k+30))                *
    !              for all other values of c(1), including  c(1) = 0,  the *
    !                 default  values  of  the  weights  are  taken  to be *
    !                 1/max(1,abs(y(k))), as mentioned earlier             *
    !           (in the two cases c(1) = 4 or 5 the user must declare  the *
    !           dimension of c to be at least n+30 and must initialize the *
    !           components c(31), c(32), ..., c(n+30).)                    *
    !                                                                      *
    !  c(2)  floor value - used when the indicator c(1) has the value 3    *
    !                                                                      *
    !  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
    !           to be abs(c(3)) - otherwise it uses the default value      *
    !              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
    !           where dwarf is a very small positive  machine  number  and *
    !           rreb is the relative roundoff error bound                  *
    !                                                                      *
    !  c(4)  hstart specification - if not zero, the subroutine  will  use *
    !           an  initial  hmag equal to abs(c(4)), except of course for *
    !           the restrictions imposed by hmin and hmax  -  otherwise it *
    !           uses the default value of hmax*(tol)**(1/6)                *
    !                                                                      *
    !  c(5)  scale specification - this is intended to be a measure of the *
    !           scale of the problem - larger values of scale tend to make *
    !           the method more reliable, first  by  possibly  restricting *
    !           hmax  (as  described  below) and second, by tightening the *
    !           acceptance requirement - if c(5) is zero, a default  value *
    !           of  1  is  used.  for  linear  homogeneous  problems  with *
    !           constant coefficients, an appropriate value for scale is a *
    !           norm  of  the  associated  matrix.  for other problems, an *
    !           approximation to  an  average  value  of  a  norm  of  the *
    !           jacobian along the trajectory may be appropriate           *
    !                                                                      *
    !  c(6)  hmax specification - four cases are possible                  *
    !           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
    !              min(abs(c(6)),2/abs(c(5)))                              *
    !           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
    !           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
    !              2/abs(c(5))                                             *
    !           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
    !              of 2                                                    *
    !                                                                      *
    !  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
    !           error  return with ind = -1 will be caused when the number *
    !           of function evaluations exceeds abs(c(7))                  *
    !                                                                      *
    !  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
    !           interrupt   the  calculations  after  it  has  chosen  its *
    !           preliminary value of hmag, and just before choosing htrial *
    !           and  xtrial  in  preparation for taking a step (htrial may *
    !           differ from hmag in sign, and may  require  adjustment  if *
    !           xend  is  near) - the subroutine returns with ind = 4, and *
    !           will resume calculation at the point  of  interruption  if *
    !           re-entered with ind = 4                                    *
    !                                                                      *
    !  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
    !           interrupt   the  calculations  immediately  after  it  has *
    !           decided whether or not to accept the result  of  the  most *
    !           recent  trial step, with ind = 5 if it plans to accept, or *
    !           ind = 6 if it plans to reject -  y(*)  is  the  previously *
    !           accepted  result, while w(*,9) is the newly computed trial *
    !           value, and w(*,2) is the unweighted error estimate vector. *
    !           the  subroutine  will  resume calculations at the point of *
    !           interruption on re-entry with ind = 5 or 6. (the user  may *
    !           change ind in this case if he wishes, for example to force *
    !           acceptance of a step that would otherwise be rejected,  or *
    !           vice versa. he can also restart with ind = 1 or 2.)        *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !  summary of the components of the communications vector              *
    !                                                                      *
    !     prescribed at the option       determined by the program         *
    !           of the user                                                *
    !                                                                      *
    !                                    c(10) rreb(rel roundoff err bnd)  *
    !     c(1) error control indicator   c(11) dwarf (very small mach no)  *
    !     c(2) floor value               c(12) weighted norm y             *
    !     c(3) hmin specification        c(13) hmin                        *
    !     c(4) hstart specification      c(14) hmag                        *
    !     c(5) scale specification       c(15) scale                       *
    !     c(6) hmax specification        c(16) hmax                        *
    !     c(7) max no of fcn evals       c(17) xtrial                      *
    !     c(8) interrupt no 1            c(18) htrial                      *
    !     c(9) interrupt no 2            c(19) est                         *
    !                                    c(20) previous xend               *
    !                                    c(21) flag for xend               *
    !                                    c(22) no of successful steps      *
    !                                    c(23) no of successive failures   *
    !                                    c(24) no of fcn evals             *
    !                                                                      *
    !  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !  an overview of the program                                          *
    !                                                                      *
    !     begin initialization, parameter checking, interrupt re-entries   *
    !  ......abort if ind out of range 1 to 6                              *
    !  .     cases - initial entry, normal re-entry, interrupt re-entries  *
    !  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
    !  v........abort if n.gt.nw or tol.le.0                               *
    !  .        if initial entry without options (ind .eq. 1)              *
    !  .           set c(1) to c(9) equal to zero                          *
    !  .        else initial entry with options (ind .eq. 2)               *
    !  .           make c(1) to c(9) non-negative                          *
    !  .           make floor values non-negative if they are to be used   *
    !  .        end if                                                     *
    !  .        initialize rreb, dwarf, prev xend, flag, counts            *
    !  .     case 2 - normal re-entry (ind .eq. 3)                         *
    !  .........abort if xend reached, and either x changed or xend not    *
    !  .        re-initialize flag                                         *
    !  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
    !  v        transfer control to the appropriate re-entry point.......  *
    !  .     end cases                                                  .  *
    !  .  end initialization, etc.                                      .  *
    !  .                                                                v  *
    !  .  loop through the following 4 stages, once for each trial step .  *
    !  .     stage 1 - prepare                                          .  *
    !***********error return (with ind=-1) if no of fcn evals too great .  *
    !  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
    !  .        calc hmin, scale, hmax                                  .  *
    !***********error return (with ind=-2) if hmin .gt. hmax            .  *
    !  .        calc preliminary hmag                                   .  *
    !***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
    !  .        calc hmag, xtrial and htrial                            .  *
    !  .     end stage 1                                                .  *
    !  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
    !  .     stage 3 - calc the error estimate                          .  *
    !  .     stage 4 - make decisions                                   .  *
    !  .        set ind=5 if step acceptable, else set ind=6            .  *
    !***********interrupt no 2 if requested....................re-entry.v  *
    !  .        if step accepted (ind .eq. 5)                              *
    !  .           update x, y from xtrial, ytrial                         *
    !  .           add 1 to no of successful steps                         *
    !  .           set no of successive failures to zero                   *
    !**************return(with ind=3, xend saved, flag set) if x .eq. xend *
    !  .        else step not accepted (ind .eq. 6)                        *
    !  .           add 1 to no of successive failures                      *
    !**************error return (with ind=-3) if hmag .le. hmin            *
    !  .        end if                                                     *
    !  .     end stage 4                                                   *
    !  .  end loop                                                         *
    !  .                                                                   *
    !  begin abort action                                                  *
    !     output appropriate  message  about  stopping  the  calculations, *
    !        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
    !        previous xend,  no of  successful  steps,  no  of  successive *
    !        failures, no of fcn evals, and the components of y            *
    !     stop                                                             *
    !  end abort action                                                    *
    !                                                                      *
    !***********************************************************************
    !
    !     ******************************************************************
    !     * begin initialization, parameter checking, interrupt re-entries *
    !     ******************************************************************
    !
    !  ......abort if ind out of range 1 to 6
    if (ind.lt.1 .or. ind.gt.6) go to 500
    !
    !        cases - initial entry, normal re-entry, interrupt re-entries
    !         go to (5, 5, 45, 1111, 2222, 2222), ind
    if (ind==3) goto 45
    if (ind==4) goto 1111
    if (ind==5 .or. ind==6) goto 2222

    !        case 1 - initial entry (ind .eq. 1 or 2)
    !  .........abort if n.gt.nw or tol.le.0
    if (n.gt.nw .or. tol.le.0._dl) go to 500
    if (ind.eq. 2) go to 15
    !              initial entry without options (ind .eq. 1)
    !              set c(1) to c(9) equal to 0
    do k = 1, 9
        c(k) = 0._dl
    end do
    go to 35
15  continue
    !              initial entry with options (ind .eq. 2)
    !              make c(1) to c(9) non-negative
    do k = 1, 9
        c(k) = dabs(c(k))
    end do
    !              make floor values non-negative if they are to be used
    if (c(1).ne.4._dl .and. c(1).ne.5._dl) go to 30
    do k = 1, n
        c(k+30) = dabs(c(k+30))
    end do
30  continue
35  continue
    !           initialize rreb, dwarf, prev xend, flag, counts
    c(10) = 2._dl**(-56)
    c(11) = 1.d-35
    !           set previous xend initially to initial value of x
    c(20) = x
    do k = 21, 24
        c(k) = 0._dl
    end do
    go to 50
    !        case 2 - normal re-entry (ind .eq. 3)
    !  .........abort if xend reached, and either x changed or xend not
45  if (c(21).ne.0._dl .and. &
        (x.ne.c(20) .or. xend.eq.c(20))) go to 500
    !           re-initialize flag
    c(21) = 0._dl
    go to 50
    !        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
    !           transfer control to the appropriate re-entry point..........
    !           this has already been handled by the computed go to        .
    !        end cases                                                     v
50  continue
    !
    !     end initialization, etc.
    !
    !     ******************************************************************
    !     * loop through the following 4 stages, once for each trial  step *
    !     * until the occurrence of one of the following                   *
    !     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
    !     *        stage 4                                                 *
    !     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
    !     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
    !     *        requested, in stage 1 or stage 4                        *
    !     ******************************************************************
    !
99999 continue
    !
    !        ***************************************************************
    !        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
    !        * and some parameter  checking,  and  end  up  with  suitable *
    !        * values of hmag, xtrial and htrial in preparation for taking *
    !        * an integration step.                                        *
    !        ***************************************************************
    !
    !***********error return (with ind=-1) if no of fcn evals too great
    if (c(7).eq.0._dl .or. c(24).lt.c(7)) go to 100
    ind = -1
    return
100 continue
    !
    !           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    if (ind .eq. 6) go to 105
    call derivst(EV, this, etat, n, x, y, w(1,1))
    c(24) = c(24) + 1._dl
105 continue
    !
    !           calculate hmin - use default unless value prescribed
    c(13) = c(3)
    if (c(3) .ne. 0._dl) go to 165
    !              calculate default value of hmin
    !              first calculate weighted norm y - c(12) - as specified
    !              by the error control indicator c(1)
    temp = 0._dl
    if (c(1) .ne. 1._dl) go to 115
    !                 absolute error control - weights are 1
    do 110 k = 1, n
        temp = dmax1(temp, dabs(y(k)))
110 continue
    c(12) = temp
    go to 160
115 if (c(1) .ne. 2._dl) go to 120
    !                 relative error control - weights are 1/dabs(y(k)) so
    !                 weighted norm y is 1
    c(12) = 1._dl
    go to 160
120 if (c(1) .ne. 3._dl) go to 130
    !                 weights are 1/max(c(2),abs(y(k)))
    do 125 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(2))
125 continue
    c(12) = dmin1(temp, 1._dl)
    go to 160
130 if (c(1) .ne. 4._dl) go to 140
    !                 weights are 1/max(c(k+30),abs(y(k)))
    do 135 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(k+30))
135 continue
    c(12) = dmin1(temp, 1._dl)
    go to 160
140 if (c(1) .ne. 5._dl) go to 150
    !                 weights are 1/c(k+30)
    do 145 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(k+30))
145 continue
    c(12) = temp
    go to 160
150 continue
    !                 default case - weights are 1/max(1,abs(y(k)))
    do 155 k = 1, n
        temp = dmax1(temp, dabs(y(k)))
155 continue
    c(12) = dmin1(temp, 1._dl)
160 continue
    c(13) = 10._dl*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
    !
    !           calculate scale - use default unless value prescribed
    c(15) = c(5)
    if (c(5) .eq. 0._dl) c(15) = 1._dl
    !
    !           calculate hmax - consider 4 cases
    !           case 1 both hmax and scale prescribed
    if (c(6).ne.0._dl .and. c(5).ne.0._dl) &
        c(16) = dmin1(c(6), 2._dl/c(5))
    !           case 2 - hmax prescribed, but scale not
    if (c(6).ne.0._dl .and. c(5).eq.0._dl) c(16) = c(6)
    !           case 3 - hmax not prescribed, but scale is
    if (c(6).eq.0._dl .and. c(5).ne.0._dl) c(16) = 2._dl/c(5)
    !           case 4 - neither hmax nor scale is provided
    if (c(6).eq.0._dl .and. c(5).eq.0._dl) c(16) = 2._dl
    !
    !***********error return (with ind=-2) if hmin .gt. hmax
    if (c(13) .le. c(16)) go to 170
    ind = -2
    return
170 continue
    !
    !           calculate preliminary hmag - consider 3 cases
    if (ind .gt. 2) go to 175
    !           case 1 - initial entry - use prescribed value of hstart, if
    !              any, else default
    c(14) = c(4)
    if (c(4) .eq. 0._dl) c(14) = c(16)*tol**(1._dl/6._dl)
    go to 185
175 if (c(23) .gt. 1._dl) go to 180
    !           case 2 - after a successful step, or at most  one  failure,
    !              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
    !              overflow. then avoid reduction by more than half.
    temp = 2._dl*c(14)
    if (tol .lt. (2._dl/.9d0)**6*c(19)) &
        temp = .9d0*(tol/c(19))**(1._dl/6._dl)*c(14)
    c(14) = dmax1(temp, .5d0*c(14))
    go to 185
180 continue
    !           case 3 - after two or more successive failures
    c(14) = .5d0*c(14)
185 continue
    !
    !           check against hmax
    c(14) = dmin1(c(14), c(16))
    !
    !           check against hmin
    c(14) = dmax1(c(14), c(13))
    !
    !***********interrupt no 1 (with ind=4) if requested
    if (c(8) .eq. 0._dl) go to 1111
    ind = 4
    return
    !           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
    !
    !           calculate hmag, xtrial - depending on preliminary hmag, xend
    if (c(14) .ge. dabs(xend - x)) go to 190
    !              do not step more than half way to xend
    c(14) = dmin1(c(14), .5d0*dabs(xend - x))
    c(17) = x + dsign(c(14), xend - x)
    go to 195
190 continue
    !              hit xend exactly
    c(14) = dabs(xend - x)
    c(17) = xend
195 continue
    !
    !           calculate htrial
    c(18) = c(17) - x
    !
    !        end stage 1
    !
    !        ***************************************************************
    !        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
    !        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
    !        * stage 3. w(*,9) is temporary storage until finally it holds *
    !        * ytrial.                                                     *
    !        ***************************************************************
    !
    temp = c(18)/1398169080000._dl
    !
    do 200 k = 1, n
        w(k,9) = y(k) + temp*w(k,1)*233028180000._dl
200 continue
    call derivst(EV, this, etat, n, x + c(18)/6._dl, w(1,9), w(1,2))
    !
    do 205 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*74569017600._dl &
            + w(k,2)*298276070400._dl  )
205 continue
    call derivst(EV, this, etat, n, x + c(18)*(4._dl/15._dl), w(1,9), w(1,3))
    !
    do 210 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*1165140900000._dl &
            - w(k,2)*3728450880000._dl &
            + w(k,3)*3495422700000._dl )
210 continue
    call derivst(EV, this, etat, n, x + c(18)*(2._dl/3._dl), w(1,9), w(1,4))
    !
    do 215 k = 1, n
        w(k,9) = y(k) + temp*( - w(k,1)*3604654659375._dl &
            + w(k,2)*12816549900000._dl &
            - w(k,3)*9284716546875._dl &
            + w(k,4)*1237962206250._dl )
215 continue
    call derivst(EV, this, etat, n, x + c(18)*(5._dl/6._dl), w(1,9), w(1,5))
    !
    do 220 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*3355605792000._dl &
            - w(k,2)*11185352640000._dl &
            + w(k,3)*9172628850000._dl &
            - w(k,4)*427218330000._dl &
            + w(k,5)*482505408000._dl  )
220 continue
    call derivst(EV, this, etat, n, x + c(18), w(1,9), w(1,6))
    !
    do 225 k = 1, n
        w(k,9) = y(k) + temp*( - w(k,1)*770204740536._dl &
            + w(k,2)*2311639545600._dl &
            - w(k,3)*1322092233000._dl &
            - w(k,4)*453006781920._dl &
            + w(k,5)*326875481856._dl  )
225 continue
    call derivst(EV, this, etat, n, x + c(18)/15._dl, w(1,9), w(1,7))
    !
    do 230 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*2845924389000._dl &
            - w(k,2)*9754668000000._dl &
            + w(k,3)*7897110375000._dl &
            - w(k,4)*192082660000._dl &
            + w(k,5)*400298976000._dl &
            + w(k,7)*201586000000._dl  )
230 continue
    call derivst(EV, this, etat, n, x + c(18), w(1,9), w(1,8))
    !
    !           calculate ytrial, the extrapolated approximation and store
    !              in w(*,9)
    do 235 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*104862681000._dl &
            + w(k,3)*545186250000._dl &
            + w(k,4)*446637345000._dl &
            + w(k,5)*188806464000._dl &
            + w(k,7)*15076875000._dl &
            + w(k,8)*97599465000._dl   )
235 continue
    !
    !           add 7 to the no of fcn evals
    c(24) = c(24) + 7._dl
    !
    !        end stage 2
    !
    !        ***************************************************************
    !        * stage 3 - calculate the error estimate est. first calculate *
    !        * the  unweighted  absolute  error  estimate vector (per unit *
    !        * step) for the unextrapolated approximation and store it  in *
    !        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
    !        * specified by the error  control  indicator  c(1).  finally, *
    !        * modify  this result to produce est, the error estimate (per *
    !        * unit step) for the extrapolated approximation ytrial.       *
    !        ***************************************************************
    !
    !           calculate the unweighted absolute error estimate vector
    do 300 k = 1, n
        w(k,2) = (   w(k,1)*8738556750._dl &
            + w(k,3)*9735468750._dl &
            - w(k,4)*9709507500._dl &
            + w(k,5)*8582112000._dl &
            + w(k,6)*95329710000._dl &
            - w(k,7)*15076875000._dl &
            - w(k,8)*97599465000._dl)/1398169080000._dl
300 continue
    !
    !           calculate the weighted max norm of w(*,2) as specified by
    !           the error control indicator c(1)
    temp = 0._dl
    if (c(1) .ne. 1._dl) go to 310
    !              absolute error control
    do 305 k = 1, n
        temp = dmax1(temp,dabs(w(k,2)))
305 continue
    go to 360
310 if (c(1) .ne. 2._dl) go to 320
    !              relative error control
    do 315 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)/y(k)))
315 continue
    go to 360
320 if (c(1) .ne. 3._dl) go to 330
    !              weights are 1/max(c(2),abs(y(k)))
    do 325 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(c(2), dabs(y(k))) )
325 continue
    go to 360
330 if (c(1) .ne. 4._dl) go to 340
    !              weights are 1/max(c(k+30),abs(y(k)))
    do 335 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(c(k+30), dabs(y(k))) )
335 continue
    go to 360
340 if (c(1) .ne. 5._dl) go to 350
    !              weights are 1/c(k+30)
    do 345 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
345 continue
    go to 360
350 continue
    !              default case - weights are 1/max(1,abs(y(k)))
    do 355 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(1._dl, dabs(y(k))) )
355 continue
360 continue
    !
    !           calculate est - (the weighted max norm of w(*,2))*hmag*scale
    !              - est is intended to be a measure of the error  per  unit
    !              step in ytrial
    c(19) = temp*c(14)*c(15)
    !
    !        end stage 3
    !
    !        ***************************************************************
    !        * stage 4 - make decisions.                                   *
    !        ***************************************************************
    !
    !           set ind=5 if step acceptable, else set ind=6
    ind = 5
    if (c(19) .gt. tol) ind = 6
    !
    !***********interrupt no 2 if requested
    if (c(9) .eq. 0._dl) go to 2222
    return
    !           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
2222 continue
    !
    if (ind .eq. 6) go to 410
    !              step accepted (ind .eq. 5), so update x, y from xtrial,
    !                 ytrial, add 1 to the no of successful steps, and set
    !                 the no of successive failures to zero
    x = c(17)
    do 400 k = 1, n
        y(k) = w(k,9)
400 continue
    c(22) = c(22) + 1._dl
    c(23) = 0._dl
    !**************return(with ind=3, xend saved, flag set) if x .eq. xend
    if (x .ne. xend) go to 405
    ind = 3
    c(20) = xend
    c(21) = 1._dl
    return
405 continue
    go to 420
410 continue
    !              step not accepted (ind .eq. 6), so add 1 to the no of
    !                 successive failures
    c(23) = c(23) + 1._dl
    !**************error return (with ind=-3) if hmag .le. hmin
    if (c(14) .gt. c(13)) go to 415
    ind = -3
    return
415 continue
420 continue
    !
    !        end stage 4
    !
    go to 99999
    !     end loop
    !
    !  begin abort action
500 continue
    !

    write (*,*) 'Error in dverk_derivst, x =',x, 'xend=', xend
    call GlobalError('DVERK_DERIVST error', error_evolution)
    end subroutine dverk_derivst

    subroutine dverk_derivsv(EV, this, etat, n, x, y, xend, tol, ind, c, nw, w)
    use Precision
    use MpiUtils
    use Config, only : GlobalError, error_evolution
    integer n, ind, nw, k
    real(dl) x, y(n), xend, tol, c(*), w(nw,9), temp
    type(EvolutionVars) EV
    class(TThermoData) :: this
    class(TRecfast) :: etat
     !it isn't, but as long as it maintains it as a pointer we are OK
    !
    !***********************************************************************
    !                                                                      *
    ! note added 11/14/85.                                                 *
    !                                                                      *
    ! if you discover any errors in this subroutine, please contact        *
    !                                                                      *
    !        kenneth r. jackson                                            *
    !        department of computer science                                *
    !        university of toronto                                         *
    !        toronto, ontario,                                             *
    !        canada   m5s 1a4                                              *
    !                                                                      *
    !        phone: 416-978-7075                                           *
    !                                                                      *
    !        electronic mail:                                              *
    !        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
    !        csnet:  krj@toronto                                           *
    !        arpa:   krj.toronto@csnet-relay                               *
    !        bitnet: krj%toronto@csnet-relay.arpa                          *
    !                                                                      *
    ! dverk is written in fortran 66.                                      *
    !                                                                      *
    ! the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
    ! set for a  vax  in  double  precision.  they  should  be  reset,  as *
    ! described below, if this program is run on another machine.          *
    !                                                                      *
    ! the c array is declared in this subroutine to have one element only, *
    ! although  more  elements  are  referenced  in this subroutine.  this *
    ! causes some compilers to issue warning messages.  there is,  though, *
    ! no  error  provided  c is declared sufficiently large in the calling *
    ! program, as described below.                                         *
    !                                                                      *
    ! the following external statement  for  fcn  was  added  to  avoid  a *
    ! warning  message  from  the  unix  f77 compiler.  the original dverk *
    ! comments and code follow it.                                         *
    !                                                                      *
    !***********************************************************************
    !
    !external fcn
    !
    !***********************************************************************
    !                                                                      *
    !     purpose - this is a runge-kutta  subroutine  based  on  verner's *
    ! fifth and sixth order pair of formulas for finding approximations to *
    ! the solution of  a  system  of  first  order  ordinary  differential *
    ! equations  with  initial  conditions. it attempts to keep the global *
    ! error proportional to  a  tolerance  specified  by  the  user.  (the *
    ! proportionality  depends  on the kind of error control that is used, *
    ! as well as the differential equation and the range of integration.)  *
    !                                                                      *
    !     various options are available to the user,  including  different *
    ! kinds  of  error control, restrictions on step sizes, and interrupts *
    ! which permit the user to examine the state of the  calculation  (and *
    ! perhaps make modifications) during intermediate stages.              *
    !                                                                      *
    !     the program is efficient for non-stiff systems.  however, a good *
    ! variable-order-adams  method  will probably be more efficient if the *
    ! function evaluations are very costly.  such a method would  also  be *
    ! more suitable if one wanted to obtain a large number of intermediate *
    ! solution values by interpolation, as might be the case  for  example *
    ! with graphical output.                                               *
    !                                                                      *
    !                                    hull-enright-jackson   1/10/76    *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !     use - the user must specify each of the following                *
    !                                                                      *
    !     n  number of equations                                           *
    !                                                                      *
    !   fcn  name of subroutine for evaluating functions - the  subroutine *
    !           itself must also be provided by the user - it should be of *
    !           the following form                                         *
    !              subroutine fcn(n, x, y, yprime)                         *
    !              integer n                                               *
    !              real(dl) x, y(n), yprime(n)                     *
    !                      *** etc ***                                     *
    !           and it should evaluate yprime, given n, x and y            *
    !                                                                      *
    !     x  independent variable - initial value supplied by user         *
    !                                                                      *
    !     y  dependent variable - initial values of components y(1), y(2), *
    !           ..., y(n) supplied by user                                 *
    !                                                                      *
    !  xend  value of x to which integration is to be carried out - it may *
    !           be less than the initial value of x                        *
    !                                                                      *
    !   tol  tolerance - the subroutine attempts to control a norm of  the *
    !           local  error  in  such  a  way  that  the  global error is *
    !           proportional to tol. in some problems there will be enough *
    !           damping  of  errors, as well as some cancellation, so that *
    !           the global error will be less than tol. alternatively, the *
    !           control   can   be  viewed  as  attempting  to  provide  a *
    !           calculated value of y at xend which is the exact  solution *
    !           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
    !           is proportional to tol.  (the norm  is  a  max  norm  with *
    !           weights  that  depend on the error control strategy chosen *
    !           by the user.  the default weight for the k-th component is *
    !           1/max(1,abs(y(k))),  which therefore provides a mixture of *
    !           absolute and relative error control.)                      *
    !                                                                      *
    !   ind  indicator - on initial entry ind must be set equal to  either *
    !           1  or  2. if the user does not wish to use any options, he *
    !           should set ind to 1 - all that remains for the user to  do *
    !           then  is  to  declare c and w, and to specify nw. the user *
    !           may also  select  various  options  on  initial  entry  by *
    !           setting ind = 2 and initializing the first 9 components of *
    !           c as described in the next section.  he may also  re-enter *
    !           the  subroutine  with ind = 3 as mentioned again below. in *
    !           any event, the subroutine returns with ind equal to        *
    !              3 after a normal return                                 *
    !              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
    !              -1, -2, or -3 after an error condition (see below)      *
    !                                                                      *
    !     c  communications vector - the dimension must be greater than or *
    !           equal to 24, unless option c(1) = 4 or 5 is used, in which *
    !           case the dimension must be greater than or equal to n+30   *
    !                                                                      *
    !    nw  first dimension of workspace w -  must  be  greater  than  or *
    !           equal to n                                                 *
    !                                                                      *
    !     w  workspace matrix - first dimension must be nw and second must *
    !           be greater than or equal to 9                              *
    !                                                                      *
    !     the subroutine  will  normally  return  with  ind  =  3,  having *
    ! replaced the initial values of x and y with, respectively, the value *
    ! of xend and an approximation to y at xend.  the  subroutine  can  be *
    ! called  repeatedly  with new values of xend without having to change *
    ! any other argument.  however, changes in tol, or any of the  options *
    ! described below, may also be made on such a re-entry if desired.     *
    !                                                                      *
    !     three error returns are also possible, in which  case  x  and  y *
    ! will be the most recently accepted values -                          *
    !     with ind = -3 the subroutine was unable  to  satisfy  the  error *
    !        requirement  with a particular step-size that is less than or *
    !        equal to hmin, which may mean that tol is too small           *
    !     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
    !        probably  means  that the requested tol (which is used in the *
    !        calculation of hmin) is too small                             *
    !     with ind = -1 the allowed maximum number of fcn evaluations  has *
    !        been  exceeded,  but  this  can only occur if option c(7), as *
    !        described in the next section, has been used                  *
    !                                                                      *
    !     there are several circumstances that will cause the calculations *
    ! to  be  terminated,  along with output of information that will help *
    ! the user determine the cause of  the  trouble.  these  circumstances *
    ! involve  entry with illegal or inconsistent values of the arguments, *
    ! such as attempting a normal  re-entry  without  first  changing  the *
    ! value of xend, or attempting to re-enter with ind less than zero.    *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !     options - if the subroutine is entered with ind = 1, the first 9 *
    ! components of the communications vector are initialized to zero, and *
    ! the subroutine uses only default values  for  each  option.  if  the *
    ! subroutine  is  entered  with ind = 2, the user must specify each of *
    ! these 9 components - normally he would first set them all  to  zero, *
    ! and  then  make  non-zero  those  that  correspond to the particular *
    ! options he wishes to select. in any event, options may be changed on *
    ! re-entry  to  the  subroutine  -  but if the user changes any of the *
    ! options, or tol, in the course of a calculation he should be careful *
    ! about  how  such changes affect the subroutine - it may be better to *
    ! restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
    ! program  -  the information is available to the user, but should not *
    ! normally be changed by him.)                                         *
    !                                                                      *
    !  c(1)  error control indicator - the norm of the local error is  the *
    !           max  norm  of  the  weighted  error  estimate  vector, the *
    !           weights being determined according to the value of c(1) -  *
    !              if c(1)=1 the weights are 1 (absolute error control)    *
    !              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
    !                 control)                                             *
    !              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
    !                 (relative  error  control,  unless abs(y(k)) is less *
    !                 than the floor value, abs(c(2)) )                    *
    !              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
    !                 (here individual floor values are used)              *
    !              if c(1)=5 the weights are 1/abs(c(k+30))                *
    !              for all other values of c(1), including  c(1) = 0,  the *
    !                 default  values  of  the  weights  are  taken  to be *
    !                 1/max(1,abs(y(k))), as mentioned earlier             *
    !           (in the two cases c(1) = 4 or 5 the user must declare  the *
    !           dimension of c to be at least n+30 and must initialize the *
    !           components c(31), c(32), ..., c(n+30).)                    *
    !                                                                      *
    !  c(2)  floor value - used when the indicator c(1) has the value 3    *
    !                                                                      *
    !  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
    !           to be abs(c(3)) - otherwise it uses the default value      *
    !              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
    !           where dwarf is a very small positive  machine  number  and *
    !           rreb is the relative roundoff error bound                  *
    !                                                                      *
    !  c(4)  hstart specification - if not zero, the subroutine  will  use *
    !           an  initial  hmag equal to abs(c(4)), except of course for *
    !           the restrictions imposed by hmin and hmax  -  otherwise it *
    !           uses the default value of hmax*(tol)**(1/6)                *
    !                                                                      *
    !  c(5)  scale specification - this is intended to be a measure of the *
    !           scale of the problem - larger values of scale tend to make *
    !           the method more reliable, first  by  possibly  restricting *
    !           hmax  (as  described  below) and second, by tightening the *
    !           acceptance requirement - if c(5) is zero, a default  value *
    !           of  1  is  used.  for  linear  homogeneous  problems  with *
    !           constant coefficients, an appropriate value for scale is a *
    !           norm  of  the  associated  matrix.  for other problems, an *
    !           approximation to  an  average  value  of  a  norm  of  the *
    !           jacobian along the trajectory may be appropriate           *
    !                                                                      *
    !  c(6)  hmax specification - four cases are possible                  *
    !           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
    !              min(abs(c(6)),2/abs(c(5)))                              *
    !           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
    !           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
    !              2/abs(c(5))                                             *
    !           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
    !              of 2                                                    *
    !                                                                      *
    !  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
    !           error  return with ind = -1 will be caused when the number *
    !           of function evaluations exceeds abs(c(7))                  *
    !                                                                      *
    !  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
    !           interrupt   the  calculations  after  it  has  chosen  its *
    !           preliminary value of hmag, and just before choosing htrial *
    !           and  xtrial  in  preparation for taking a step (htrial may *
    !           differ from hmag in sign, and may  require  adjustment  if *
    !           xend  is  near) - the subroutine returns with ind = 4, and *
    !           will resume calculation at the point  of  interruption  if *
    !           re-entered with ind = 4                                    *
    !                                                                      *
    !  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
    !           interrupt   the  calculations  immediately  after  it  has *
    !           decided whether or not to accept the result  of  the  most *
    !           recent  trial step, with ind = 5 if it plans to accept, or *
    !           ind = 6 if it plans to reject -  y(*)  is  the  previously *
    !           accepted  result, while w(*,9) is the newly computed trial *
    !           value, and w(*,2) is the unweighted error estimate vector. *
    !           the  subroutine  will  resume calculations at the point of *
    !           interruption on re-entry with ind = 5 or 6. (the user  may *
    !           change ind in this case if he wishes, for example to force *
    !           acceptance of a step that would otherwise be rejected,  or *
    !           vice versa. he can also restart with ind = 1 or 2.)        *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !  summary of the components of the communications vector              *
    !                                                                      *
    !     prescribed at the option       determined by the program         *
    !           of the user                                                *
    !                                                                      *
    !                                    c(10) rreb(rel roundoff err bnd)  *
    !     c(1) error control indicator   c(11) dwarf (very small mach no)  *
    !     c(2) floor value               c(12) weighted norm y             *
    !     c(3) hmin specification        c(13) hmin                        *
    !     c(4) hstart specification      c(14) hmag                        *
    !     c(5) scale specification       c(15) scale                       *
    !     c(6) hmax specification        c(16) hmax                        *
    !     c(7) max no of fcn evals       c(17) xtrial                      *
    !     c(8) interrupt no 1            c(18) htrial                      *
    !     c(9) interrupt no 2            c(19) est                         *
    !                                    c(20) previous xend               *
    !                                    c(21) flag for xend               *
    !                                    c(22) no of successful steps      *
    !                                    c(23) no of successive failures   *
    !                                    c(24) no of fcn evals             *
    !                                                                      *
    !  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
    !                                                                      *
    !***********************************************************************
    !                                                                      *
    !  an overview of the program                                          *
    !                                                                      *
    !     begin initialization, parameter checking, interrupt re-entries   *
    !  ......abort if ind out of range 1 to 6                              *
    !  .     cases - initial entry, normal re-entry, interrupt re-entries  *
    !  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
    !  v........abort if n.gt.nw or tol.le.0                               *
    !  .        if initial entry without options (ind .eq. 1)              *
    !  .           set c(1) to c(9) equal to zero                          *
    !  .        else initial entry with options (ind .eq. 2)               *
    !  .           make c(1) to c(9) non-negative                          *
    !  .           make floor values non-negative if they are to be used   *
    !  .        end if                                                     *
    !  .        initialize rreb, dwarf, prev xend, flag, counts            *
    !  .     case 2 - normal re-entry (ind .eq. 3)                         *
    !  .........abort if xend reached, and either x changed or xend not    *
    !  .        re-initialize flag                                         *
    !  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
    !  v        transfer control to the appropriate re-entry point.......  *
    !  .     end cases                                                  .  *
    !  .  end initialization, etc.                                      .  *
    !  .                                                                v  *
    !  .  loop through the following 4 stages, once for each trial step .  *
    !  .     stage 1 - prepare                                          .  *
    !***********error return (with ind=-1) if no of fcn evals too great .  *
    !  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
    !  .        calc hmin, scale, hmax                                  .  *
    !***********error return (with ind=-2) if hmin .gt. hmax            .  *
    !  .        calc preliminary hmag                                   .  *
    !***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
    !  .        calc hmag, xtrial and htrial                            .  *
    !  .     end stage 1                                                .  *
    !  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
    !  .     stage 3 - calc the error estimate                          .  *
    !  .     stage 4 - make decisions                                   .  *
    !  .        set ind=5 if step acceptable, else set ind=6            .  *
    !***********interrupt no 2 if requested....................re-entry.v  *
    !  .        if step accepted (ind .eq. 5)                              *
    !  .           update x, y from xtrial, ytrial                         *
    !  .           add 1 to no of successful steps                         *
    !  .           set no of successive failures to zero                   *
    !**************return(with ind=3, xend saved, flag set) if x .eq. xend *
    !  .        else step not accepted (ind .eq. 6)                        *
    !  .           add 1 to no of successive failures                      *
    !**************error return (with ind=-3) if hmag .le. hmin            *
    !  .        end if                                                     *
    !  .     end stage 4                                                   *
    !  .  end loop                                                         *
    !  .                                                                   *
    !  begin abort action                                                  *
    !     output appropriate  message  about  stopping  the  calculations, *
    !        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
    !        previous xend,  no of  successful  steps,  no  of  successive *
    !        failures, no of fcn evals, and the components of y            *
    !     stop                                                             *
    !  end abort action                                                    *
    !                                                                      *
    !***********************************************************************
    !
    !     ******************************************************************
    !     * begin initialization, parameter checking, interrupt re-entries *
    !     ******************************************************************
    !
    !  ......abort if ind out of range 1 to 6
    if (ind.lt.1 .or. ind.gt.6) go to 500
    !
    !        cases - initial entry, normal re-entry, interrupt re-entries
    !         go to (5, 5, 45, 1111, 2222, 2222), ind
    if (ind==3) goto 45
    if (ind==4) goto 1111
    if (ind==5 .or. ind==6) goto 2222

    !        case 1 - initial entry (ind .eq. 1 or 2)
    !  .........abort if n.gt.nw or tol.le.0
    if (n.gt.nw .or. tol.le.0._dl) go to 500
    if (ind.eq. 2) go to 15
    !              initial entry without options (ind .eq. 1)
    !              set c(1) to c(9) equal to 0
    do k = 1, 9
        c(k) = 0._dl
    end do
    go to 35
15  continue
    !              initial entry with options (ind .eq. 2)
    !              make c(1) to c(9) non-negative
    do k = 1, 9
        c(k) = dabs(c(k))
    end do
    !              make floor values non-negative if they are to be used
    if (c(1).ne.4._dl .and. c(1).ne.5._dl) go to 30
    do k = 1, n
        c(k+30) = dabs(c(k+30))
    end do
30  continue
35  continue
    !           initialize rreb, dwarf, prev xend, flag, counts
    c(10) = 2._dl**(-56)
    c(11) = 1.d-35
    !           set previous xend initially to initial value of x
    c(20) = x
    do k = 21, 24
        c(k) = 0._dl
    end do
    go to 50
    !        case 2 - normal re-entry (ind .eq. 3)
    !  .........abort if xend reached, and either x changed or xend not
45  if (c(21).ne.0._dl .and. &
        (x.ne.c(20) .or. xend.eq.c(20))) go to 500
    !           re-initialize flag
    c(21) = 0._dl
    go to 50
    !        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
    !           transfer control to the appropriate re-entry point..........
    !           this has already been handled by the computed go to        .
    !        end cases                                                     v
50  continue
    !
    !     end initialization, etc.
    !
    !     ******************************************************************
    !     * loop through the following 4 stages, once for each trial  step *
    !     * until the occurrence of one of the following                   *
    !     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
    !     *        stage 4                                                 *
    !     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
    !     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
    !     *        requested, in stage 1 or stage 4                        *
    !     ******************************************************************
    !
99999 continue
    !
    !        ***************************************************************
    !        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
    !        * and some parameter  checking,  and  end  up  with  suitable *
    !        * values of hmag, xtrial and htrial in preparation for taking *
    !        * an integration step.                                        *
    !        ***************************************************************
    !
    !***********error return (with ind=-1) if no of fcn evals too great
    if (c(7).eq.0._dl .or. c(24).lt.c(7)) go to 100
    ind = -1
    return
100 continue
    !
    !           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    if (ind .eq. 6) go to 105
    call derivsv(EV, this, etat, n, x, y, w(1,1))
    c(24) = c(24) + 1._dl
105 continue
    !
    !           calculate hmin - use default unless value prescribed
    c(13) = c(3)
    if (c(3) .ne. 0._dl) go to 165
    !              calculate default value of hmin
    !              first calculate weighted norm y - c(12) - as specified
    !              by the error control indicator c(1)
    temp = 0._dl
    if (c(1) .ne. 1._dl) go to 115
    !                 absolute error control - weights are 1
    do 110 k = 1, n
        temp = dmax1(temp, dabs(y(k)))
110 continue
    c(12) = temp
    go to 160
115 if (c(1) .ne. 2._dl) go to 120
    !                 relative error control - weights are 1/dabs(y(k)) so
    !                 weighted norm y is 1
    c(12) = 1._dl
    go to 160
120 if (c(1) .ne. 3._dl) go to 130
    !                 weights are 1/max(c(2),abs(y(k)))
    do 125 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(2))
125 continue
    c(12) = dmin1(temp, 1._dl)
    go to 160
130 if (c(1) .ne. 4._dl) go to 140
    !                 weights are 1/max(c(k+30),abs(y(k)))
    do 135 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(k+30))
135 continue
    c(12) = dmin1(temp, 1._dl)
    go to 160
140 if (c(1) .ne. 5._dl) go to 150
    !                 weights are 1/c(k+30)
    do 145 k = 1, n
        temp = dmax1(temp, dabs(y(k))/c(k+30))
145 continue
    c(12) = temp
    go to 160
150 continue
    !                 default case - weights are 1/max(1,abs(y(k)))
    do 155 k = 1, n
        temp = dmax1(temp, dabs(y(k)))
155 continue
    c(12) = dmin1(temp, 1._dl)
160 continue
    c(13) = 10._dl*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
    !
    !           calculate scale - use default unless value prescribed
    c(15) = c(5)
    if (c(5) .eq. 0._dl) c(15) = 1._dl
    !
    !           calculate hmax - consider 4 cases
    !           case 1 both hmax and scale prescribed
    if (c(6).ne.0._dl .and. c(5).ne.0._dl) &
        c(16) = dmin1(c(6), 2._dl/c(5))
    !           case 2 - hmax prescribed, but scale not
    if (c(6).ne.0._dl .and. c(5).eq.0._dl) c(16) = c(6)
    !           case 3 - hmax not prescribed, but scale is
    if (c(6).eq.0._dl .and. c(5).ne.0._dl) c(16) = 2._dl/c(5)
    !           case 4 - neither hmax nor scale is provided
    if (c(6).eq.0._dl .and. c(5).eq.0._dl) c(16) = 2._dl
    !
    !***********error return (with ind=-2) if hmin .gt. hmax
    if (c(13) .le. c(16)) go to 170
    ind = -2
    return
170 continue
    !
    !           calculate preliminary hmag - consider 3 cases
    if (ind .gt. 2) go to 175
    !           case 1 - initial entry - use prescribed value of hstart, if
    !              any, else default
    c(14) = c(4)
    if (c(4) .eq. 0._dl) c(14) = c(16)*tol**(1._dl/6._dl)
    go to 185
175 if (c(23) .gt. 1._dl) go to 180
    !           case 2 - after a successful step, or at most  one  failure,
    !              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
    !              overflow. then avoid reduction by more than half.
    temp = 2._dl*c(14)
    if (tol .lt. (2._dl/.9d0)**6*c(19)) &
        temp = .9d0*(tol/c(19))**(1._dl/6._dl)*c(14)
    c(14) = dmax1(temp, .5d0*c(14))
    go to 185
180 continue
    !           case 3 - after two or more successive failures
    c(14) = .5d0*c(14)
185 continue
    !
    !           check against hmax
    c(14) = dmin1(c(14), c(16))
    !
    !           check against hmin
    c(14) = dmax1(c(14), c(13))
    !
    !***********interrupt no 1 (with ind=4) if requested
    if (c(8) .eq. 0._dl) go to 1111
    ind = 4
    return
    !           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
    !
    !           calculate hmag, xtrial - depending on preliminary hmag, xend
    if (c(14) .ge. dabs(xend - x)) go to 190
    !              do not step more than half way to xend
    c(14) = dmin1(c(14), .5d0*dabs(xend - x))
    c(17) = x + dsign(c(14), xend - x)
    go to 195
190 continue
    !              hit xend exactly
    c(14) = dabs(xend - x)
    c(17) = xend
195 continue
    !
    !           calculate htrial
    c(18) = c(17) - x
    !
    !        end stage 1
    !
    !        ***************************************************************
    !        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
    !        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
    !        * stage 3. w(*,9) is temporary storage until finally it holds *
    !        * ytrial.                                                     *
    !        ***************************************************************
    !
    temp = c(18)/1398169080000._dl
    !
    do 200 k = 1, n
        w(k,9) = y(k) + temp*w(k,1)*233028180000._dl
200 continue
    call derivsv(EV, this, etat, n, x + c(18)/6._dl, w(1,9), w(1,2))
    !
    do 205 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*74569017600._dl &
            + w(k,2)*298276070400._dl  )
205 continue
    call derivsv(EV, this, etat, n, x + c(18)*(4._dl/15._dl), w(1,9), w(1,3))
    !
    do 210 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*1165140900000._dl &
            - w(k,2)*3728450880000._dl &
            + w(k,3)*3495422700000._dl )
210 continue
    call derivsv(EV, this, etat, n, x + c(18)*(2._dl/3._dl), w(1,9), w(1,4))
    !
    do 215 k = 1, n
        w(k,9) = y(k) + temp*( - w(k,1)*3604654659375._dl &
            + w(k,2)*12816549900000._dl &
            - w(k,3)*9284716546875._dl &
            + w(k,4)*1237962206250._dl )
215 continue
    call derivsv(EV, this, etat, n, x + c(18)*(5._dl/6._dl), w(1,9), w(1,5))
    !
    do 220 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*3355605792000._dl &
            - w(k,2)*11185352640000._dl &
            + w(k,3)*9172628850000._dl &
            - w(k,4)*427218330000._dl &
            + w(k,5)*482505408000._dl  )
220 continue
    call derivsv(EV, this, etat, n, x + c(18), w(1,9), w(1,6))
    !
    do 225 k = 1, n
        w(k,9) = y(k) + temp*( - w(k,1)*770204740536._dl &
            + w(k,2)*2311639545600._dl &
            - w(k,3)*1322092233000._dl &
            - w(k,4)*453006781920._dl &
            + w(k,5)*326875481856._dl  )
225 continue
    call derivsv(EV, this, etat, n, x + c(18)/15._dl, w(1,9), w(1,7))
    !
    do 230 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*2845924389000._dl &
            - w(k,2)*9754668000000._dl &
            + w(k,3)*7897110375000._dl &
            - w(k,4)*192082660000._dl &
            + w(k,5)*400298976000._dl &
            + w(k,7)*201586000000._dl  )
230 continue
    call derivsv(EV, this, etat, n, x + c(18), w(1,9), w(1,8))
    !
    !           calculate ytrial, the extrapolated approximation and store
    !              in w(*,9)
    do 235 k = 1, n
        w(k,9) = y(k) + temp*(   w(k,1)*104862681000._dl &
            + w(k,3)*545186250000._dl &
            + w(k,4)*446637345000._dl &
            + w(k,5)*188806464000._dl &
            + w(k,7)*15076875000._dl &
            + w(k,8)*97599465000._dl   )
235 continue
    !
    !           add 7 to the no of fcn evals
    c(24) = c(24) + 7._dl
    !
    !        end stage 2
    !
    !        ***************************************************************
    !        * stage 3 - calculate the error estimate est. first calculate *
    !        * the  unweighted  absolute  error  estimate vector (per unit *
    !        * step) for the unextrapolated approximation and store it  in *
    !        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
    !        * specified by the error  control  indicator  c(1).  finally, *
    !        * modify  this result to produce est, the error estimate (per *
    !        * unit step) for the extrapolated approximation ytrial.       *
    !        ***************************************************************
    !
    !           calculate the unweighted absolute error estimate vector
    do 300 k = 1, n
        w(k,2) = (   w(k,1)*8738556750._dl &
            + w(k,3)*9735468750._dl &
            - w(k,4)*9709507500._dl &
            + w(k,5)*8582112000._dl &
            + w(k,6)*95329710000._dl &
            - w(k,7)*15076875000._dl &
            - w(k,8)*97599465000._dl)/1398169080000._dl
300 continue
    !
    !           calculate the weighted max norm of w(*,2) as specified by
    !           the error control indicator c(1)
    temp = 0._dl
    if (c(1) .ne. 1._dl) go to 310
    !              absolute error control
    do 305 k = 1, n
        temp = dmax1(temp,dabs(w(k,2)))
305 continue
    go to 360
310 if (c(1) .ne. 2._dl) go to 320
    !              relative error control
    do 315 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)/y(k)))
315 continue
    go to 360
320 if (c(1) .ne. 3._dl) go to 330
    !              weights are 1/max(c(2),abs(y(k)))
    do 325 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(c(2), dabs(y(k))) )
325 continue
    go to 360
330 if (c(1) .ne. 4._dl) go to 340
    !              weights are 1/max(c(k+30),abs(y(k)))
    do 335 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(c(k+30), dabs(y(k))) )
335 continue
    go to 360
340 if (c(1) .ne. 5._dl) go to 350
    !              weights are 1/c(k+30)
    do 345 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
345 continue
    go to 360
350 continue
    !              default case - weights are 1/max(1,abs(y(k)))
    do 355 k = 1, n
        temp = dmax1(temp, dabs(w(k,2)) &
            / dmax1(1._dl, dabs(y(k))) )
355 continue
360 continue
    !
    !           calculate est - (the weighted max norm of w(*,2))*hmag*scale
    !              - est is intended to be a measure of the error  per  unit
    !              step in ytrial
    c(19) = temp*c(14)*c(15)
    !
    !        end stage 3
    !
    !        ***************************************************************
    !        * stage 4 - make decisions.                                   *
    !        ***************************************************************
    !
    !           set ind=5 if step acceptable, else set ind=6
    ind = 5
    if (c(19) .gt. tol) ind = 6
    !
    !***********interrupt no 2 if requested
    if (c(9) .eq. 0._dl) go to 2222
    return
    !           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
2222 continue
    !
    if (ind .eq. 6) go to 410
    !              step accepted (ind .eq. 5), so update x, y from xtrial,
    !                 ytrial, add 1 to the no of successful steps, and set
    !                 the no of successive failures to zero
    x = c(17)
    do 400 k = 1, n
        y(k) = w(k,9)
400 continue
    c(22) = c(22) + 1._dl
    c(23) = 0._dl
    !**************return(with ind=3, xend saved, flag set) if x .eq. xend
    if (x .ne. xend) go to 405
    ind = 3
    c(20) = xend
    c(21) = 1._dl
    return
405 continue
    go to 420
410 continue
    !              step not accepted (ind .eq. 6), so add 1 to the no of
    !                 successive failures
    c(23) = c(23) + 1._dl
    !**************error return (with ind=-3) if hmag .le. hmin
    if (c(14) .gt. c(13)) go to 415
    ind = -3
    return
415 continue
420 continue
    !
    !        end stage 4
    !
    go to 99999
    !     end loop
    !
    !  begin abort action
500 continue
    !

    write (*,*) 'Error in dverk_derivsv, x =',x, 'xend=', xend
    call GlobalError('DVERK_DERIVSV error', error_evolution)
    end subroutine dverk_derivsv

    end module GaugeInterface
