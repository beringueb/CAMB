    ! Modules used by cmbmain and other routines.

    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info) and Anthony Challinor
    !     See readme.html for documentation.
    !
    !     Based on CMBFAST  by  Uros Seljak and Matias Zaldarriaga, itself based
    !     on Boltzmann code written by Edmund Bertschinger, Chung-Pei Ma and Paul Bode.
    !     Original CMBFAST copyright and disclaimer:
    !
    !     Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
    !     the Massachusetts Institute of Technology.  All rights reserved.
    !
    !     THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. OR C.f.A. MAKE NO
    !     REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
    !     By way of example, but not limitation,
    !     M.I.T. AND C.f.A MAKE NO REPRESENTATIONS OR WARRANTIES OF
    !     MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
    !     THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
    !     ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
    !
    !     portions of this software are based on the COSMICS package of
    !     E. Bertschinger.  See the LICENSE file of the COSMICS distribution
    !     for restrictions on the modification and distribution of this software.

    module Transfer
    !Module for matter transfer function/matter power spectrum
    use Recombination
    use constants, only : const_pi, const_twopi
    use MiscUtils
    use RangeUtils
    use StringUtils
    use MathUtils
    use config
    use model
    use splines
    implicit none
    public

    integer, parameter :: Transfer_kh =1, Transfer_cdm=2,Transfer_b=3,Transfer_g=4, &
        Transfer_r=5, Transfer_nu = 6,  & !massless and massive neutrino
        Transfer_tot=7, Transfer_nonu=8, Transfer_tot_de=9,  &
        ! total perturbations with and without neutrinos, with neutrinos+dark energy in the numerator
        Transfer_Weyl = 10, & ! the Weyl potential, for lensing and ISW
        Transfer_Newt_vel_cdm=11, Transfer_Newt_vel_baryon=12,   & ! -k v_Newtonian/H
        Transfer_vel_baryon_cdm = 13 !relative velocity of baryons and CDM
    !Sources
    !Alternatively for 21cm
    integer, parameter :: Transfer_monopole=4, Transfer_vnewt=5, Transfer_Tmat = 6

    integer, parameter :: Transfer_max = Transfer_vel_baryon_cdm
    character(LEN=name_tag_len) :: Transfer_name_tags(Transfer_max-1) = &
        ['CDM     ', 'baryon  ', 'photon  ', 'nu      ', 'mass_nu ', 'total   ', &
        'no_nu   ', 'total_de', 'Weyl    ', 'v_CDM   ', 'v_b     ', 'v_b-v_c ']
    character(LEN=name_tag_len) :: Transfer21cm_name_tags(Transfer_max-1) = &
        ['CDM      ', 'baryon   ', 'photon   ', 'monopole ', 'v_newt   ', 'delta_T_g', &
        'no_nu    ', 'total_de ', 'Weyl     ', 'v_CDM    ', 'v_b      ', 'v_b-v_c  ']

    logical :: transfer_interp_matterpower  = .true. !output regular grid in log k
    !set to false to output calculated values for later interpolation

    integer :: transfer_power_var = Transfer_tot
    !What to use to calulcate the output matter power spectrum and sigma_8
    !Transfer_tot uses total matter perturbation

    !Sources
    Type Cl21cmVars
        Type(MatterPowerData), pointer :: PK
        integer l, itf
        logical logs
        real(dl) chi
    end Type Cl21cmVars

    interface Transfer_GetMatterPower
    module procedure Transfer_GetMatterPowerD,Transfer_GetMatterPowerS
    end interface

    contains

    subroutine Transfer_GetUnsplinedPower(State, M, PK,var1,var2, hubble_units)
    !Get 2pi^2/k^3 T_1 T_2 P_R(k)
    Type(MatterTransferData) :: M
    Type(CAMBdata) :: State
    real(dl), intent(inout):: PK(:,:)
    integer, optional, intent(in) :: var1
    integer, optional, intent(in) :: var2
    logical, optional, intent(in) :: hubble_units
    real(dl) :: h, k
    integer :: nz, nk, zix, ik, s1, s2
    logical :: hnorm

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)
    hnorm = DefaultTrue (hubble_units)

    nk=M%num_q_trans
    nz=State%CP%Transfer%PK_num_redshifts
    if (nk/= size(PK,1) .or. nz/=size(PK,2)) call MpiStop('Trasfer_GetUnsplinedPower wrong size')
    h = State%CP%H0/100

    if (s1==transfer_power_var .and. s2==transfer_power_var .and. allocated(State%CAMB_Pk)) then
        !Already computed
        do zix=1,nz
            PK(:,zix) = exp(State%CAMB_pk%matpower(:, State%PK_redshifts_index(nz-zix+1)))
        end do
        if (.not. hnorm) PK=  PK / h**3
    else

        do ik=1,nk
            k = M%TransferData(Transfer_kh,ik,1)*h
            do zix=1,nz
                PK(ik,zix) = M%TransferData(s1,ik,State%PK_redshifts_index(nz-zix+1))*&
                    M%TransferData(s2,ik,State%PK_redshifts_index(nz-zix+1))*k*&
                    const_pi*const_twopi*State%CP%InitPower%ScalarPower(k)
            end do
        end do

        if (hnorm) PK=  PK * h**3
    end if

    end subroutine Transfer_GetUnsplinedPower

    subroutine Transfer_GetNonLinRatio_index(State,M, ratio, itf)
    Type(MatterTransferData), intent(in) :: M
    Type(CAMBdata) :: State
    real(dl), allocatable, intent(out) :: ratio(:)
    integer, intent(in) :: itf
    Type(MatterPowerData) :: PKdata

    if (allocated(State%CAMB_PK)) then
        allocate(ratio, source = State%CAMB_PK%nonlin_ratio(:,itf))
    else
        call Transfer_GetMatterPowerData(State, M, PKdata,itf)
        call State%CP%NonLinearModel%GetNonLinRatios(State,PKdata)
        allocate(ratio, source = PKdata%nonlin_ratio(:,1))
    end if
    end subroutine Transfer_GetNonLinRatio_index


    subroutine Transfer_GetUnsplinedNonlinearPower(State,M, PK,var1,var2, hubble_units)
    !Get 2pi^2/k^3 T_1 T_2 P_R(k) after re-scaling for non-linear evolution (if turned on)
    Type(MatterTransferData), intent(in) :: M
    Type(CAMBdata) :: State
    real(dl), intent(inout):: PK(:,:)
    integer, optional, intent(in) :: var1
    integer, optional, intent(in) :: var2
    logical, optional, intent(in) :: hubble_units
    integer zix
    real(dl), allocatable :: ratio(:)

    if (.not. allocated(State%CAMB_Pk) .and. State%CP%Transfer%PK_num_redshifts == State%num_transfer_redshifts &
        .and. .not. State%OnlyTransfer) then
        allocate(State%CAMB_PK)
        call Transfer_GetMatterPowerData(State, State%MT, State%CAMB_PK)
        call State%CP%NonLinearModel%GetNonLinRatios(State, State%CAMB_PK)
    end if

    call Transfer_GetUnsplinedPower(State,M,PK,var1,var2, hubble_units)
    do zix=1, State%CP%Transfer%PK_num_redshifts
        call Transfer_GetNonLinRatio_index(State, M, ratio,State%PK_redshifts_index(State%CP%Transfer%PK_num_redshifts-zix+1))
        PK(:,zix) =  PK(:,zix) *ratio**2
    end do

    end subroutine Transfer_GetUnsplinedNonlinearPower

    subroutine Transfer_GetMatterPowerData(State, MTrans, PK_data, itf_only, var1, var2)
    !Does *NOT* include non-linear corrections
    !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
    !Here the definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
    !spectrum is generated to beyond the CMB k_max
    class(CAMBdata) :: State
    Type(MatterTransferData), intent(in) :: MTrans
    Type(MatterPowerData) :: PK_data
    integer, intent(in), optional :: itf_only
    integer, intent(in), optional :: var1, var2
    real(dl) :: h, kh, k, power
    integer :: ik, nz, itf, itf_start, itf_end, s1, s2

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)

    if (present(itf_only)) then
        itf_start=itf_only
        itf_end = itf_only
        nz = 1
    else
        itf_start=1
        nz= size(MTrans%TransferData,3)
        itf_end = nz
    end if
    PK_data%num_k = MTrans%num_q_trans
    PK_Data%num_z = nz

    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
    allocate(PK_data%log_kh(PK_data%num_k))
    allocate(PK_data%redshifts(nz))
    PK_data%redshifts = State%Transfer_Redshifts(itf_start:itf_end)

    h = State%CP%H0/100

    do ik=1,MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,1)
        k = kh*h
        PK_data%log_kh(ik) = log(kh)
        power = State%CP%InitPower%ScalarPower(k)
        if (global_error_flag/=0) then
            call MatterPowerdata_Free(PK_data)
            return
        end if
        do itf = 1, nz
            PK_data%matpower(ik,itf) = &
                log(MTrans%TransferData(s1,ik,itf_start+itf-1)*&
                MTrans%TransferData(s2,ik,itf_start+itf-1)*k &
                *const_pi*const_twopi*h**3*power)
        end do
    end do

    call MatterPowerdata_getsplines(PK_data)

    end subroutine Transfer_GetMatterPowerData

    subroutine MatterPowerData_Load(PK_data,fname)
    !Loads in kh, P_k from file for one redshiftr and one initial power spectrum
    !Not redshift is not stored in file, so not set correctly
    !Also note that output _matterpower file is already interpolated, so re-interpolating is probs not a good idea

    !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    use FileUtils
    character(LEN=*) :: fname
    type(MatterPowerData) :: PK_data
    type(TTextFile) :: F
    real(dl)kh, Pk
    integer ik
    integer nz

    nz = 1
    call F%Open(fname)

    PK_data%num_k = F%Lines()
    PK_Data%num_z = 1

    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
    allocate(PK_data%log_kh(PK_data%num_k))

    allocate(PK_data%redshifts(nz))
    PK_data%redshifts = 0

    do ik=1,PK_data%num_k
        read (F%unit, *) kh, Pk
        PK_data%matpower(ik,1) = log(Pk)
        PK_data%log_kh(ik) = log(kh)
    end do

    call MatterPowerdata_getsplines(PK_data)
    call F%Close()

    end subroutine MatterPowerData_Load


    subroutine MatterPowerdata_getsplines(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i

    do i = 1,PK_Data%num_z
        call spline_def(PK_data%log_kh,PK_data%matpower(:,i),PK_data%num_k,&
            PK_data%ddmat(:,i))
    end do

    end subroutine MatterPowerdata_getsplines

    !Sources
    subroutine MatterPowerdata_getsplines21cm(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i

    do i = 1,PK_Data%num_z
        call spline_def(PK_data%log_k,PK_data%matpower(:,i),PK_data%num_k,&
            PK_data%ddmat(:,i))
        call spline_def(PK_data%log_k,PK_data%vvpower(:,i),PK_data%num_k,&
            PK_data%ddvvpower(:,i))
        call spline_def(PK_data%log_k,PK_data%vdpower(:,i),PK_data%num_k,&
            PK_data%ddvdpower(:,i))
    end do

    end subroutine MatterPowerdata_getsplines21cm

    subroutine MatterPowerdata_Free(PK_data)
    Type(MatterPowerData) :: PK_data
    integer i
    !this shouldn't be needed when releasing the object.
    deallocate(PK_data%log_kh,stat=i)
    deallocate(PK_data%matpower,stat=i)
    deallocate(PK_data%ddmat,stat=i)
    deallocate(PK_data%nonlin_ratio,stat=i)
    deallocate(PK_data%redshifts,stat=i)
    !Sources
    deallocate(PK_data%log_k,stat=i)
    deallocate(PK_data%nonlin_ratio_vv,stat=i)
    deallocate(PK_data%nonlin_ratio_vd,stat=i)
    deallocate(PK_data%vvpower,stat=i)
    deallocate(PK_data%ddvvpower,stat=i)
    deallocate(PK_data%vdpower,stat=i)
    deallocate(PK_data%ddvdpower,stat=i)

    end subroutine MatterPowerdata_Free

    function MatterPowerData_k(PK, kh, itf, index_cache) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    Type(MatterPowerData) :: PK
    integer, intent(in) :: itf
    real (dl), intent(in) :: kh
    real(dl) :: logk
    integer llo,lhi
    real(dl) outpower, dp
    real(dl) ho,a0,b0
    integer, optional :: index_cache
    integer, save :: i_last = 1

    logk = log(kh)
    if (logk < PK%log_kh(1)) then
        dp = (PK%matpower(2,itf) -  PK%matpower(1,itf)) / &
            ( PK%log_kh(2)-PK%log_kh(1) )
        outpower = PK%matpower(1,itf) + dp*(logk - PK%log_kh(1))
    else if (logk > PK%log_kh(PK%num_k)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter

        dp = (PK%matpower(PK%num_k,itf) -  PK%matpower(PK%num_k-1,itf)) / &
            ( PK%log_kh(PK%num_k)-PK%log_kh(PK%num_k-1) )
        outpower = PK%matpower(PK%num_k,itf) + dp*(logk - PK%log_kh(PK%num_k))
    else
        llo=min(PresentDefault(i_last, index_cache),PK%num_k)
        do while (PK%log_kh(llo) > logk)
            llo=llo-1
        end do
        do while (PK%log_kh(llo+1) < logk)
            llo=llo+1
        end do
        if (present(index_cache)) then
            index_cache = llo
        else
            i_last =llo
        end if
        lhi=llo+1
        ho=PK%log_kh(lhi)-PK%log_kh(llo)
        a0=(PK%log_kh(lhi)-logk)/ho
        b0=1-a0

        outpower = a0*PK%matpower(llo,itf)+ b0*PK%matpower(lhi,itf)+&
            ((a0**3-a0)* PK%ddmat(llo,itf) &
            +(b0**3-b0)*PK%ddmat(lhi,itf))*ho**2/6
    end if

    outpower = exp(outpower)

    end function MatterPowerData_k

    !Sources
    subroutine MatterPower21cm_k(PK, k, itf, monopole, vv, vd)
    !Get monopole and velocity power at particular k by interpolation
    Type(MatterPowerData) :: PK
    integer, intent(in) :: itf
    real (dl), intent(in) :: k
    real(dl), intent(out) :: monopole, vv, vd
    real(dl) :: logk
    integer llo,lhi
    real(dl) ho,a0,b0,f1,f2,f3
    integer, save :: i_last = 1

    logk = log(k)
    if (logk < PK%log_k(1)) then
        monopole = 0
        vv=0
        return
    end if
    if (logk > PK%log_k(PK%num_k)) then
        monopole=0
        vv=0
        return
        !stop 'MatterPower21cm_k: out of bounds'
    else
        llo=min(i_last,PK%num_k)
        do while (PK%log_k(llo) > logk)
            llo=llo-1
        end do
        do while (PK%log_k(llo+1)< logk)
            llo=llo+1
        end do
        i_last =llo
        lhi=llo+1
        ho=PK%log_k(lhi)-PK%log_k(llo)
        a0=(PK%log_k(lhi)-logk)/ho
        b0=1-a0
        f1= (a0**3-a0)
        f2= (b0**3-b0)
        f3= ho**2/6


        monopole = a0*PK%matpower(llo,itf)+ b0*PK%matpower(lhi,itf)+&
            (f1* PK%ddmat(llo,itf) &
            +f2*PK%ddmat(lhi,itf))*f3
        vv = a0*PK%vvpower(llo,itf)+ b0*PK%vvpower(lhi,itf)+&
            (f1* PK%ddvvpower(llo,itf) &
            +f2*PK%ddvvpower(lhi,itf))*f3

        vd = a0*PK%vdpower(llo,itf)+ b0*PK%vdpower(lhi,itf)+&
            (f1* PK%ddvdpower(llo,itf) &
            +f2*PK%ddvdpower(lhi,itf))*f3
    end if

    monopole = exp(max(-30._dl,monopole))
    vv = exp(max(-30._dl,vv))
    vd = exp(max(-30._dl,vd))

    end subroutine MatterPower21cm_k


    subroutine Transfer_GetMatterPowerS(State, MTrans, outpower, itf, minkh, dlnkh, npoints, var1, var2)
    class(CAMBdata) :: state
    Type(MatterTransferData), intent(in) :: MTrans
    integer, intent(in) :: itf, npoints
    integer, intent(in), optional :: var1, var2
    real, intent(out) :: outpower(*)
    real, intent(in) :: minkh, dlnkh
    real(dl) :: outpowerd(npoints)
    real(dl):: minkhd, dlnkhd

    minkhd = minkh; dlnkhd = dlnkh
    call Transfer_GetMatterPowerD(State, MTrans, outpowerd, itf, minkhd, dlnkhd, npoints,var1, var2)
    outpower(1:npoints) = outpowerd(1:npoints)

    end subroutine Transfer_GetMatterPowerS

    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
    !Changed input variable from itf to itf_PK because we are looking for the itf_PK'th
    !redshift in the PK_redshifts array.  The position of this redshift in the master redshift
    !array, itf, is given by itf = CP%Transfer%Pk_redshifts_index(itf_PK)
    !Also changed (CP%NonLinear/=NonLinear_None) to
    !CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
    subroutine Transfer_GetMatterPowerD(State, MTrans, outpower, itf_PK, minkh, dlnkh, npoints, var1, var2)
    !Allows for non-smooth priordial spectra
    !if CP%Nonlinear/ = NonLinear_none includes non-linear evolution
    !Get total matter power spectrum at logarithmically equal intervals dlnkh of k/h starting at minkh
    !in units of (h Mpc^{-1})^3.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
    !sepctrum is generated to beyond the CMB k_max
    class(CAMBdata) :: state
    Type(MatterTransferData), intent(in) :: MTrans

    integer, intent(in) :: itf_PK, npoints
    real(dl), intent(out) :: outpower(npoints)
    real(dl), intent(in) :: minkh, dlnkh
    integer, intent(in), optional :: var1, var2

    integer ik, llo,il,lhi,lastix
    real(dl) matpower(MTrans%num_q_trans), kh, kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
    real(dl) atransfer,xi, a0, b0, ho, logmink,k, h
    integer itf
    integer :: s1,s2, sign
    logical log_interp
    real(dl), allocatable :: ratio(:)

    s1 = PresentDefault (transfer_power_var, var1)
    s2 = PresentDefault (transfer_power_var, var2)

    itf = State%PK_redshifts_index(itf_PK)

    if (npoints < 2) call MpiStop('Need at least 2 points in Transfer_GetMatterPower')

    if (minkh*exp((npoints-1)*dlnkh) > MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) &
        .and. FeedbackLevel > 0 ) &
        write(*,*) 'Warning: extrapolating matter power in Transfer_GetMatterPower'


    if (State%CP%NonLinear/=NonLinear_none .and. State%CP%NonLinear/=NonLinear_Lens) then
        call Transfer_GetNonLinRatio_index(State, Mtrans, ratio, itf)
    end if

    h = State%CP%H0/100
    logmink = log(minkh)
    do ik=1,MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,itf)
        k = kh*h
        kvals(ik) = log(kh)
        atransfer=MTrans%TransferData(s1,ik,itf)*MTrans%TransferData(s2,ik,itf)
        if (State%CP%NonLinear/=NonLinear_none .and. State%CP%NonLinear/=NonLinear_Lens) &
            atransfer = atransfer* ratio(ik)**2 !only one element, this itf
        matpower(ik) = atransfer*k*const_pi*const_twopi*h**3
        !Put in power spectrum later: transfer functions should be smooth, initial power may not be
    end do
    sign = 1
    log_interp = .true.
    if (any(matpower <= 0)) then
        if (all(matpower < 0)) then
            sign = -1
        else
            log_interp = .false.
        endif
    endif
    if (log_interp) matpower = log(sign*matpower)
    call spline_def(kvals,matpower,MTrans%num_q_trans, ddmat)

    llo=1
    lastix = npoints + 1
    do il=1, npoints
        xi=logmink + dlnkh*(il-1)
        if (xi < kvals(1)) then
            outpower(il)=-30.
            cycle
        end if
        do while ((xi > kvals(llo+1)).and.(llo < MTrans%num_q_trans))
            llo=llo+1
            if (llo >= MTrans%num_q_trans) exit
        end do
        if (llo == MTrans%num_q_trans) then
            lastix = il
            exit
        end if
        lhi=llo+1
        ho=kvals(lhi)-kvals(llo)
        a0=(kvals(lhi)-xi)/ho
        b0=(xi-kvals(llo))/ho

        outpower(il) = a0*matpower(llo)+ b0*matpower(lhi)+((a0**3-a0)* ddmat(llo) &
            +(b0**3-b0)*ddmat(lhi))*ho**2/6
    end do

    do while (lastix <= npoints)
        !Do linear extrapolation in the log
        !Obviouly inaccurate, non-linear etc, but OK if only using in tails of window functions
        outpower(lastix) = 2*outpower(lastix-1) - outpower(lastix-2)
        lastix = lastix+1
    end do

    if (log_interp) then
        outpower = sign*exp(max(-30.d0,outpower))
    end if
    associate(InitPower => State%CP%InitPower)
        do il = 1, npoints
            k = exp(logmink + dlnkh*(il-1))*h
            outpower(il) = outpower(il) * InitPower%ScalarPower(k)
            if (global_error_flag /= 0) exit
        end do
    end associate

    end subroutine Transfer_GetMatterPowerD

    subroutine Transfer_Get_SigmaR(State, MTrans, R, outvals, var1, var2, root)
    !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
    !of radius R h^{-1} Mpc, for all requested redshifts
    !set va1, var2 e.g. to get the value from some combination of transfer functions rather than total
    class(CAMBdata) :: State
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in) :: R
    integer, intent(in), optional :: var1, var2
    logical, intent(in), optional :: root !if true, give sigma8, otherwise sigma8^2
    real(dl), intent(out) :: outvals(:)
    real(dl) :: kh, k, h, x, win, lnk, dlnk, lnko, powers
    real(dl), dimension(State%CP%Transfer%PK_num_redshifts) :: dsig8, dsig8o, sig8, sig8o
    integer :: s1, s2, ik

    s1 = PresentDefault(transfer_power_var, var1)
    s2 = PresentDefault(transfer_power_var, var2)
    H=State%CP%h0/100._dl
    lnko=0
    dsig8o=0
    sig8=0
    sig8o=0
    do ik=1, MTrans%num_q_trans
        kh = MTrans%TransferData(Transfer_kh,ik,1)
        if (kh==0) cycle
        k = kh*H

        dsig8 = MTrans%TransferData(s1,ik, State%PK_redshifts_index(1:State%CP%Transfer%PK_num_redshifts))
        if (s1==s2) then
            dsig8 = dsig8**2
        else
            dsig8 = dsig8*MTrans%TransferData(s2,ik, State%PK_redshifts_index(1:State%CP%Transfer%PK_num_redshifts))
        end if
        x= kh *R
        if (x < 1e-2_dl) then
            win = 1._dl - x**2/10
        else
            win = 3*(sin(x)-x*cos(x))/x**3
        end if
        lnk=log(k)
        if (ik==1) then
            dlnk=0.5_dl
            !Approx for 2._dl/(Params%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
            !Contribution should be very small in any case
        else
            dlnk=lnk-lnko
        end if
        powers = State%CP%InitPower%ScalarPower(k)
        dsig8=(win*k**2)**2*powers*dsig8
        sig8=sig8+(dsig8+dsig8o)*dlnk/2
        dsig8o=dsig8
        lnko=lnk
    end do

    if (present(root)) then
        if (root) sig8 =sqrt(sig8)
    else
        sig8 =sqrt(sig8)
    end if
    outvals(1:State%CP%Transfer%PK_num_redshifts) = sig8

    end subroutine Transfer_Get_SigmaR

    subroutine Transfer_GetSigmaRArray(State, MTrans, R, sigmaR, redshift_ix, var1, var2)
    !Get array of SigmaR at (by default) redshift zero, for all values of R (in h^{-1}Mpc units)
    class(CAMBdata), target :: State
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in) :: R(:)
    real(dl), intent(out) :: SigmaR(:)
    integer, intent(in), optional :: redshift_ix, var1, var2
    integer red_ix, ik, subk
    real(dl) kh, k, h, dkh, k_step
    real(dl) lnk, dlnk, lnko, minR, maxR
    real(dl), dimension(size(R)) ::  x, win, dsig8, dsig8o, sig8, sig8o
    type(MatterPowerData), target:: PKspline
    type(MatterPowerData), pointer :: PK
    integer PK_ix
    integer :: nsub
    integer index_cache, nextra

    index_cache = 1
    nsub = 5 !interpolation steps
    h=State%CP%h0/100._dl
    minR = minval(R)/ h
    maxR = maxval(R)
    red_ix = PresentDefault(State%PK_redshifts_index(State%CP%Transfer%PK_num_redshifts), redshift_ix)

    if (allocated(State%CAMB_PK) .and. var1==transfer_power_var .and. var2==transfer_power_var) then
        PK => State%CAMB_PK
        PK_ix = red_ix
    else
        call Transfer_GetMatterPowerData(State, MTrans, PKspline, red_ix, var1, var2)
        PK => PKspline
        PK_ix = 1
    end if
    dkh = 0._dl
    lnko=0
    dsig8o=0
    sig8=0
    sig8o=0
    if (MTrans%TransferData(Transfer_kh,1,1)==0) call MpiStop('Transfer_GetSigmaRArray kh zero')

    !Steps to extrapolate beyond kmax for tail [could do analytically]
    nextra = (4 /minR - State%CP%Transfer%kmax)/h/ &
        (MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,1)- MTrans%TransferData(Transfer_kh,MTrans%num_q_trans-1,1))
    do ik=1, MTrans%num_q_trans + max(2, nextra)
        if (ik < MTrans%num_q_trans) then
            k_step= MTrans%TransferData(Transfer_kh,ik+1,1)-MTrans%TransferData(Transfer_kh,ik,1)
            kh = MTrans%TransferData(Transfer_kh,ik,1)
            nsub = max(1, nint(k_step*min(maxR,4/(kh*h))/0.5))
        else
            !after last step just extrapolate with previous size
            nsub = max(1, nint(k_step*min(maxR,4/(kh*h))))
        end if
        nsub = nint(nsub*State%CP%Accuracy%AccuracyBoost)
        dkh = k_step/nsub
        do subk = 1, nsub
            k = kh*h
            lnk=log(k)

            x= kh *R
            where (x < 1e-2_dl)
                win = 1._dl - x**2/5
            elsewhere
                win =(3*(sin(x)-x*cos(x))/x**3)**2
            end where
            if (ik==1 .and. subk==1) then
                dlnk=0.5_dl
                !Approx for 2._dl/(Params%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
                !Contribution should be very small in any case
            else
                dlnk=lnk-lnko
            end if
            dsig8=win*(MatterPowerData_k(PK, kh,PK_ix,index_cache)*k**3)
            sig8=sig8+(dsig8+dsig8o)*dlnk/2
            dsig8o=dsig8
            lnko=lnk
            kh = kh + dkh
        end do
    end do

    SigmaR=sqrt(sig8/(const_pi*const_twopi*h**3))

    end subroutine Transfer_GetSigmaRArray

    subroutine Transfer_Get_sigma8(State, MTrans, R, var1, var2)
    !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
    !of radius R h^{-1} Mpc
    ! set va1, var2 e.g. to get the value from some combination of transfer functions rather than total
    class(CAMBdata) :: State
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in), optional :: R
    integer, intent(in), optional :: var1, var2
    real(dl) :: radius

    if (global_error_flag /= 0) return

    radius = PresentDefault (8._dl, R)

    call Transfer_Get_SigmaR(State, MTrans, radius, MTrans%sigma_8, var1,var2)

    end subroutine Transfer_Get_sigma8

    subroutine Transfer_Get_sigmas(State, MTrans, R, var_delta, var_v)
    !Get sigma8 and sigma_{delta v} (for growth, like f sigma8 in LCDM)
    class(CAMBdata) :: State
    Type(MatterTransferData) :: MTrans
    real(dl), intent(in), optional :: R
    integer, intent(in), optional :: var_delta, var_v
    real(dl) :: radius
    integer :: s1, s2

    if (global_error_flag /= 0) return

    radius = PresentDefault (8._dl, R)
    s1 = PresentDefault (transfer_power_var, var_delta)
    s2 = PresentDefault (Transfer_Newt_vel_cdm, var_v)

    call Transfer_Get_SigmaR(State, MTrans, radius, MTrans%sigma_8, s1,s1)
    if (State%get_growth_sigma8) call Transfer_Get_SigmaR(State, MTrans, radius, &
        MTrans%sigma2_vdelta_8(:), s1, s2, root=.false.)

    end subroutine Transfer_Get_sigmas

    subroutine Transfer_output_Sig8(MTrans, State)
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CAMBdata), intent(in) :: State
    integer j_PK

    do j_PK=1, State%CP%Transfer%PK_num_redshifts
        write(*,'("at z =",f7.3," sigma8 (all matter) = ",f7.4)') &
            State%CP%Transfer%PK_redshifts(j_PK), MTrans%sigma_8(j_PK)
    end do
    if (State%get_growth_sigma8) then
        do j_PK=1, State%CP%Transfer%PK_num_redshifts
            write(*,'("at z =",f7.3," sigma8^2_vd/sigma8  = ",f7.4)') &
                State%CP%Transfer%PK_redshifts(j_PK), MTrans%sigma2_vdelta_8(j_PK)/MTrans%sigma_8(j_PK)
        end do
    end if

    end subroutine Transfer_output_Sig8


    subroutine Transfer_Allocate(MTrans, State)
    Type(MatterTransferData) :: MTrans
    class(CAMBdata) :: State

    call MTrans%Free()
    allocate(MTrans%q_trans(MTrans%num_q_trans))
    allocate(MTrans%TransferData(Transfer_max,MTrans%num_q_trans,State%num_transfer_redshifts))
    allocate(MTrans%sigma_8(State%CP%Transfer%PK_num_redshifts))
    if (State%get_growth_sigma8) allocate(MTrans%sigma2_vdelta_8(State%CP%Transfer%PK_num_redshifts))

    end subroutine Transfer_Allocate

    subroutine Transfer_SaveToFiles(MTrans,State,FileNames)
    use constants
    Type(MatterTransferData), intent(in) :: MTrans
    class(CAMBdata) :: State
    integer i,ik
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    integer i_PK
    integer unit

    do i_PK=1, State%CP%Transfer%PK_num_redshifts
        if (FileNames(i_PK) /= '') then
            i = State%PK_redshifts_index(i_PK)
            if (State%CP%do21cm) then
                unit = open_file_header(FileNames(i_PK), 'k/h', Transfer21cm_name_tags, 14)
            else
                unit = open_file_header(FileNames(i_PK), 'k/h', transfer_name_tags, 14)
            end if
            do ik=1,MTrans%num_q_trans
                if (MTrans%TransferData(Transfer_kh,ik,i)/=0) then
                    write(unit,'(*(E15.6))') MTrans%TransferData(Transfer_kh:Transfer_max,ik,i)
                end if
            end do
            close(unit)
        end if
    end do

    end subroutine Transfer_SaveToFiles

    subroutine Transfer_SaveMatterPower(MTrans, State,FileNames, all21cm)
    use constants
    !Export files of total  matter power spectra in h^{-1} Mpc units, against k/h.
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CAMBdata) :: State
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    character(LEN=name_tag_len) :: columns(3)
    integer itf, i, unit
    integer points
    real, dimension(:,:), allocatable :: outpower
    real minkh,dlnkh
    Type(MatterPowerData) :: PK_data
    integer ncol
    logical, intent(in), optional :: all21cm
    logical all21
    !JD 08/13 Changes in here to PK arrays and variables
    integer itf_PK

    all21 = DefaultFalse(all21cm)
    if (all21) then
        ncol = 3
    else
        ncol = 1
    end if

    do itf=1, State%CP%Transfer%PK_num_redshifts
        if (FileNames(itf) /= '') then
            if (.not. transfer_interp_matterpower ) then
                itf_PK = State%PK_redshifts_index(itf)

                points = MTrans%num_q_trans
                allocate(outpower(points,ncol))

                !Sources
                if (all21) then
                    call Transfer_Get21cmPowerData(MTrans, State, PK_data, itf_PK)
                else
                    call Transfer_GetMatterPowerData(State, MTrans, PK_data, itf_PK)
                    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
                    !Changed (CP%NonLinear/=NonLinear_None) to CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
                    if(State%CP%NonLinear/=NonLinear_none .and. State%CP%NonLinear/=NonLinear_Lens) then
                        call State%CP%NonLinearModel%GetNonLinRatios(State, PK_data)
                        PK_data%matpower = PK_data%matpower +  2*log(PK_data%nonlin_ratio)
                        call MatterPowerdata_getsplines(PK_data)
                    end if
                end if

                outpower(:,1) = exp(PK_data%matpower(:,1))
                !Sources
                if (all21) then
                    outpower(:,3) = exp(PK_data%vvpower(:,1))
                    outpower(:,2) = exp(PK_data%vdpower(:,1))

                    outpower(:,1) = outpower(:,1)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,2) = outpower(:,2)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,3) = outpower(:,3)/1d10*const_pi*const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                end if

                call MatterPowerdata_Free(PK_Data)
                columns = ['P   ', 'P_vd','P_vv']
                unit = open_file_header(FileNames(itf), 'k/h', columns(:ncol), 15)
                do i=1,points
                    write (unit, '(*(E15.6))') MTrans%TransferData(Transfer_kh,i,1),outpower(i,:)
                end do
                close(unit)
            else
                if (all21) stop 'Transfer_SaveMatterPower: if output all assume not interpolated'
                minkh = 1e-4
                dlnkh = 0.02
                points = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/dlnkh+1
                !             dlnkh = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/(points-0.999)
                allocate(outpower(points,1))
                call Transfer_GetMatterPowerS(State, MTrans, outpower(1,1), itf, minkh,dlnkh, points)

                columns(1) = 'P'
                unit = open_file_header(FileNames(itf), 'k/h', columns(:1), 15)

                do i=1,points
                    write (unit, '(*(E15.6))') minkh*exp((i-1)*dlnkh),outpower(i,1)
                end do
                close(unit)
            end if

            deallocate(outpower)
        end if
    end do

    end subroutine Transfer_SaveMatterPower


    subroutine Transfer_Get21cmPowerData(MTrans, State, PK_data, z_ix)
    !In terms of k, not k/h, and k^3 P_k /2pi rather than P_k
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CAMBdata) :: State
    Type(MatterPowerData) :: PK_data, PK_cdm
    real(dl) h, k, pow
    integer ik
    integer z_ix,nz

    nz = 1
    PK_data%num_k = MTrans%num_q_trans
    PK_Data%num_z = nz

    allocate(PK_data%matpower(PK_data%num_k,nz))
    allocate(PK_data%ddmat(PK_data%num_k,nz))
    allocate(PK_data%vvpower(PK_data%num_k,nz))
    allocate(PK_data%ddvvpower(PK_data%num_k,nz))
    allocate(PK_data%vdpower(PK_data%num_k,nz))
    allocate(PK_data%ddvdpower(PK_data%num_k,nz))
    allocate(PK_data%log_k(PK_data%num_k))
    allocate(PK_data%redshifts(nz))

    PK_data%redshifts = State%Transfer_Redshifts(z_ix)

    h = State%CP%H0/100

    if (State%CP%NonLinear/=NonLinear_None .and. State%CP%NonLinear/=NonLinear_Lens) then
        if (z_ix>1) stop 'not tested more than one redshift with Nonlinear 21cm'
        call Transfer_GetMatterPowerData(State, MTrans, PK_cdm, z_ix)
        call State%CP%NonLinearModel%GetNonLinRatios_All(State,PK_cdm)
    end if

    do ik=1,MTrans%num_q_trans
        k = MTrans%TransferData(Transfer_kh,ik,z_ix)*h
        PK_data%log_k(ik) = log(k)
        pow = State%CP%InitPower%ScalarPower(k)*1d10

        PK_data%matpower(ik,1) = &
            log( (MTrans%TransferData(Transfer_monopole,ik,z_ix)*k**2)**2 * pow)
        PK_data%vvpower(ik,1) = &
            log( (MTrans%TransferData(Transfer_vnewt ,ik,z_ix)*k**2)**2 * pow)
        PK_data%vdpower(ik,1) = &
            log( abs((MTrans%TransferData(Transfer_vnewt ,ik,z_ix)*k**2)*&
            (MTrans%TransferData(Transfer_monopole,ik,z_ix)*k**2))* pow)

        if (State%CP%NonLinear/=NonLinear_None) then
            PK_data%matpower(ik,1) = PK_data%matpower(ik,1) + 2*log(PK_cdm%nonlin_ratio(ik,z_ix))
            PK_data%vvpower(ik,1) = PK_data%vvpower(ik,1) + 2*log(PK_cdm%nonlin_ratio_vv(ik,z_ix))
            PK_data%vdpower(ik,1) = PK_data%vdpower(ik,1) + 2*log(PK_cdm%nonlin_ratio_vd(ik,z_ix))
        end if

    end do

    if (State%CP%NonLinear/=NonLinear_None)  call MatterPowerdata_Free(PK_cdm)


    call MatterPowerdata_getsplines21cm(PK_data)

    end subroutine Transfer_Get21cmPowerData

    function Get21cmCl_l(Vars,kin)
    !Direct integration with j^2/r, etc.
    class(Cl21cmVars) Vars
    real(dl) kin, x, jl,ddJl,k, jlm1
    real(dl) Get21cmCl_l
    real(dl) monopole, vv , vd
    external BJL_EXTERNAL
    if (Vars%logs) then
        k = exp(kin)
    else
        k = kin
    end if
    x= Vars%chi*k

    call MatterPower21cm_k(Vars%PK, k, Vars%itf, monopole, vv, vd)
    call bjl_external(Vars%l, x, jl)
    call bjl_external(Vars%l-1, x, jlm1)
    ddjl = -( 2/x*jlm1-(Vars%l+2)*real(Vars%l+1,dl)/x**2*jl + jl)

    Get21cmCl_l = jl**2*monopole + ddjl**2*vv - 2._dl *ddjl*jl*vd
    if (.not. Vars%logs)  Get21cmCl_l =  Get21cmCl_l / k

    end function Get21cmCl_l


    function Get21cmCl_l_avg(Vars,kin)
    !Asymptotic results where we take <cos^2>=1/2 assuming smooth power spectrum
    class(Cl21cmVars) Vars
    real(dl) kin, x, jl,ddJl,cross,k
    real(dl) Get21cmCl_l_avg
    real(dl) monopole, vv , vd,lphalf
    external BJL_EXTERNAL

    if (Vars%logs) then
        k = exp(kin)
    else
        k = kin
    end if
    x= Vars%chi*k

    call MatterPower21cm_k(Vars%PK, k, Vars%itf, monopole, vv, vd)
    lphalf=Vars%l+0.5_dl

    jl = 1/(2*x**2) /sqrt(1-(lphalf/x)**2)

    !  ddjl = (4/x**4+1)/(2*x**2)
    !
    ddjl = (x**4-2*x**2*lphalf**2+lphalf**4)/(x**4*sqrt(x**2-lphalf**2)*x)/2

    !    cross = (2-x**2)/(2*x**4)

    cross = (-x**2+lphalf**2)/(x**2*sqrt(x**2-lphalf**2)*x)/2

    Get21cmCl_l_avg = jl*monopole + ddjl*vv - 2._dl *cross*vd
    if (.not. Vars%logs)  Get21cmCl_l_avg =  Get21cmCl_l_avg / k

    !       Get21cmCl_l_avg=Get21cmCl_l_avg
    end function Get21cmCl_l_avg


    subroutine Transfer_Get21cmCls(MTrans, State,FileNames)
    use constants
    !Get 21cm C_l from sharp shell, using only monopole source and redshift distortions
    Type(MatterTransferData), intent(in) :: MTrans
    Type(CAMBdata), target :: State
    character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
    integer itf,ik, itf_PK
    integer points
    character(LEN=name_tag_len), dimension(3), parameter :: Transfer_21cm_name_tags = &
        ['CL  ','P   ','P_vv']
    Type(MatterPowerData), target ::PK_data
    real(dl)  tol,atol, chi, Cl
    integer l, lastl, unit
    real(dl) k_min, k_max,k, avg_fac
    Type(Cl21cmVars) vars
    Type(CAMBParams), pointer :: CP

    CP=>State%CP

    tol = 1e-5/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)
    do itf_PK=1, CP%Transfer%PK_num_redshifts
        itf = State%PK_redshifts_index(itf_PK)
        if (FileNames(itf) /= '') then
            chi = State%tau0-State%TimeOfz(CP%Transfer%PK_redshifts(itf_PK))

            points = MTrans%num_q_trans

            lastl=0

            call Transfer_Get21cmPowerData(MTrans, State, PK_data, itf)

            unit = open_file_header(FileNames(itf_PK), 'L', Transfer_21cm_name_tags, 8)

            do ik =1, points-1
                k =exp(PK_data%log_k(ik))
                l=nint(k*chi)
                !This is not an approximation, we are just chosing to sample l at values around (k samples)*chi

                if (l>1 .and. l/= lastl) then
                    lastl=l
                    Vars%l=l
                    Vars%chi = chi
                    Vars%PK => PK_data
                    Vars%itf = 1
                    Cl=0
                    atol = tol
                    avg_fac = 200
                    k_min = max(exp(PK_data%log_k(1)), k*(1-20*CP%Accuracy%AccuracyBoost/chi))
                    k_max = CP%Accuracy%AccuracyBoost*max(k*(1+avg_fac/chi), k*(1._dl+real(l,dl)**(-0.666_dl)))

                    if (k_max*chi < l+10) k_max = (l+10)/chi

                    Vars%logs = .false.
                    if (k_max < exp(PK_data%log_k(points))) then
                        !Integrate bessels properly
                        Cl = Integrate_Romberg(Vars,Get21cmCl_l,k_min,k_max, atol, 25)

                        Vars%logs = .true.
                        k_min = log(k_max)


                        if (l>2e6) then
                            !In baryon damping
                            !                     Vars%logs = .false.
                            !                     atol = tol/10
                            !                     k_min = max(exp(PK_data%log_k(1)), k*(1-10/chi) )
                            !                     k_max = k*(1+100/chi)

                            k_max = min(log(5*k), PK_data%log_k(points))

                        elseif (l>1e4) then
                            Vars%logs = .false.
                            k_min = k_max

                            k_max = min(k*35*CP%Accuracy%AccuracyBoost, exp(PK_data%log_k(points)))
                        else
                            !In white noise regime
                            k_max = min(log(max(0.3_dl,k)*18*CP%Accuracy%AccuracyBoost), PK_data%log_k(points))
                        end if

                        Cl = Cl+Integrate_Romberg(Vars,Get21cmCl_l_avg,k_min,k_max, atol, 25)
                    else
                        k_max = exp(PK_data%log_k(points))
                        Cl = Integrate_Romberg(Vars,Get21cmCl_l,k_min,k_max, atol, 25)
                    end if


                    Cl=exp(-2*State%optical_depths_for21cm(itf_PK))*const_fourpi*Cl* &
                        real(l,dl)*(l+1)/const_twopi/1d10

                    write (unit, '(1I8,3E15.5)') l, Cl, exp(PK_data%matpower(ik,1)/1d10), exp(PK_data%vvpower(ik,1)/1d10)
                end if
            end do

            close(unit)

            call MatterPowerdata_Free(PK_Data)
        end if
    end do

    end subroutine Transfer_Get21cmCls


    end module Transfer
