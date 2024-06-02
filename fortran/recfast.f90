    !Recombination module for CAMB, using RECFAST

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !C Integrator for Cosmic Recombination of Hydrogen and Helium,
    !C developed by Douglas Scott (dscott@astro.ubc.ca)
    !C based on calculations in the paper Seager, Sasselov & Scott
    !C (ApJ, 523, L1, 1999).
    !and "fudge" updates in Wong, Moss & Scott (2008).
    !C
    !C Permission to use, copy, modify and distribute without fee or royalty at
    !C any tier, this software and its documentation, for any purpose and without
    !C fee or royalty is hereby granted, provided that you agree to comply with
    !C the following copyright notice and statements, including the disclaimer,
    !C and that the same appear on ALL copies of the software and documentation,
    !C including modifications that you make for internal use or for distribution:
    !C
    !C Copyright 1999-2010 by University of British Columbia.  All rights reserved.
    !C
    !C THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO
    !C REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
    !C BY WAY OF EXAMPLE, BUT NOT LIMITATION,
    !c U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF
    !C MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
    !C THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
    !C ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
    !C
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    !CN     Name:        RECFAST
    !CV     Version: 1.5.2
    !C
    !CP     Purpose:  Calculate ionised fraction as a function of redshift.
    !CP            Solves for H and He simultaneously, and includes
    !CP           H "fudge factor" for low z effect, as well as
    !CP           HeI fudge factor.
    !C
    !CD     Description: Solves for ionisation history since recombination
    !CD     using the equations in Seager, Sasselov & Scott (ApJ, 1999).
    !CD     The Cosmological model can be flat or open.
    !CD  The matter temperature is also followed, with an update from
    !CD  Scott & Scott (2009).
    !CD  The values for \alpha_B for H are from Hummer (1994).
    !CD  The singlet HeI coefficient is a fit from the full code.
    !CD  Additional He "fudge factors" are as described in Wong, Moss
    !CD  and Scott (2008).
    !CD  Extra fitting function included (in optical depth) to account
    !CD  for extra H physics described in Rubino-Martin et al. (2010).
    !CD  Care is taken to use the most accurate constants.
    !C
    !CA     Arguments:
    !CA     Name, Description
    !CA     real(dl) throughout
    !CA
    !CA     z is redshift - W is sqrt(1+z), like conformal time
    !CA     x is total ionised fraction, relative to H
    !CA     x_H is ionized fraction of H - y(1) in R-K routine
    !CA     x_He is ionized fraction of He - y(2) in R-K routine
    !CA       (note that x_He=n_He+/n_He here and not n_He+/n_H)
    !CA     Tmat is matter temperature - y(3) in R-K routine
    !CA     f's are the derivatives of the Y's
    !CA     alphaB is case B recombination rate
    !CA     alpHe is the singlet only HeII recombination rate
    !CA     a_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
    !CA     b_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
    !CA     c_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
    !CA     d_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
    !CA     a_VF is Verner and Ferland type fitting parameter for Helium
    !CA     b_VF is Verner and Ferland type fitting parameter for Helium
    !CA     T_0 is Verner and Ferland type fitting parameter for Helium
    !CA     T_1 is Verner and Ferland type fitting parameter for Helium
    !CA     Tnow is the observed CMB temperature today
    !CA     Yp is the primordial helium abundace
    !CA     fHe is He/H number ratio = Yp/4(1-Yp)
    !CA     Trad and Tmat are radiation and matter temperatures
    !CA     epsilon is the approximate difference (=Trad-Tmat) at high z
    !CA     OmegaB is Omega in baryons today
    !CA     H is Hubble constant in units of 100 km/s/Mpc
    !CA     HO is Hubble constant in SI units
    !CA     bigH is 100 km/s/Mpc in SI units
    !CA     Hz is the value of H at the specific z (in ION)
    !CA     G is grvitational constant
    !CA     n is number density of hydrogen
    !CA     Nnow is number density today
    !CA     x0 is initial ionized fraction
    !CA     x_H0 is initial ionized fraction of Hydrogen
    !CA     x_He0 is initial ionized fraction of Helium
    !CA     rhs is dummy for calculating x0
    !CA     zinitial and zfinal are starting and ending redshifts
    !CA     zeq is the redshift of matter-radiation equality
    !CA     zstart and zend are for each pass to the integrator
    !CA     C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
    !CA     m_e,m_H: electron mass and mass of H atom in SI
    !CA     not4: ratio of 4He atomic mass to 1H atomic mass
    !CA     sigma: Thomson cross-section
    !CA     a_rad: radiation constant for u=aT^4
    !CA     Lambda: 2s-1s two photon rate for Hydrogen
    !CA     Lambda_He: 2s-1s two photon rate for Helium
    !CA     DeltaB: energy of first excited state from continuum = 3.4eV
    !CA     DeltaB_He: energy of first excited state from cont. for He = 3.4eV
    !CA     L_H_ion: level for H ionization in m^-1
    !CA     L_H_alpha: level for H Ly alpha in m^-1
    !CA     L_He1_ion: level for HeI ionization
    !CA     L_He2_ion: level for HeII ionization
    !CA     L_He_2s: level for HeI 2s
    !CA     L_He_2p: level for HeI 2p (21P1-11S0) in m^-1
    !CA     Lalpha: Ly alpha wavelength in SI
    !CA     Lalpha_He: Helium I 2p-1s wavelength in SI
    !CA     mu_H,mu_T: mass per H atom and mass per particle
    !CA     H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
    !CA     CDB=DeltaB/k_B                     Constants derived from B1,B2,R
    !CA     CDB_He=DeltaB_He/k_B  n=2-infinity for He in Kelvin
    !CA     CB1=CDB*4.         Lalpha and sigma_Th, calculated
    !CA     CB1_He1: CB1 for HeI ionization potential
    !CA     CB1_He2: CB1 for HeII ionization potential
    !CA     CR=2*Pi*(m_e/h_P)*(k_B/h_P)  once and passed in a common block
    !CA     CK=Lalpha**3/(8.*Pi)
    !CA     CK_He=Lalpha_He**3/(8.*Pi)
    !CA     CL=C*h_P/(k_B*Lalpha)
    !CA     CL_He=C*h_P/(k_B*Lalpha_He)
    !CA     CT=(8./3.)*(sigma/(m_e*C))*a
    !CA     Bfact=exp((E_2p-E_2s)/kT)    Extra Boltzmann factor
    !CA b_He= "fudge factor" for HeI, to approximate higher z behaviour
    !CA Heswitch=integer for modifying HeI recombination
    !CA Parameters and quantities to describe the extra triplet states
    !CA  and also the continuum opacity of H, with a fitting function
    !CA  suggested by KIV, astro-ph/0703438
    !CA a_trip: used to fit HeI triplet recombination rate
    !CA b_trip: used to fit HeI triplet recombination rate
    !CA L_He_2Pt: level for 23P012-11S0 in m^-1
    !CA L_He_2St: level for 23S1-11S0 in m^-1
    !CA L_He2St_ion: level for 23S1-continuum in m^-1
    !CA A2P_s: Einstein A coefficient for He 21P1-11S0
    !CA A2P_t: Einstein A coefficient for He 23P1-11S0
    !CA sigma_He_2Ps: H ionization x-section at HeI 21P1-11S0 freq. in m^2
    !CA sigma_He_2Pt: H ionization x-section at HeI 23P1-11S0 freq. in m^2
    !CA CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
    !CA CfHe_t: triplet statistical correction
    !CA Hswitch is an boolean for modifying the H recombination
    !CA AGauss1 is the amplitude of the 1st Gaussian for the H fudging
    !CA AGauss2 is the amplitude of the 2nd Gaussian for the H fudging
    !CA zGauss1 is the ln(1+z) central value of the 1st Gaussian
    !CA zGauss2 is the ln(1+z) central value of the 2nd Gaussian
    !CA wGauss1 is the width of the 1st Gaussian
    !CA wGauss2 is the width of the 2nd Gaussian


    !CA     tol: tolerance for the integrator
    !CA     cw(24),w(3,9): work space for DVERK
    !CA     Ndim: number of d.e.'s to solve (integer)
    !CA     Nz: number of output redshitf (integer)
    !CA     I: loop index (integer)
    !CA     ind,nw: work-space for DVERK (integer)
    !C
    !CF     File & device access:
    !CF     Unit /I,IO,O  /Name (if known)
    !C
    !CM     Modules called:
    !CM     DVERK (numerical integrator)
    !CM     GET_INIT (initial values for ionization fractions)
    !CM     ION (ionization and Temp derivatives)
    !C
    !CC     Comments:
    !CC     none
    !C
    !CH     History:
    !CH     CREATED            (simplest version) 19th March 1989
    !CH     RECREATED    11th January 1995
    !CH               includes variable Cosmology
    !CH               uses DVERK integrator
    !CH               initial conditions are Saha
    !CH     TESTED              a bunch, well, OK, not really
    !CH     MODIFIED     January 1995 (include Hummer's 1994 alpha table)
    !CH               January 1995 (include new value for 2s-1s rate)
    !CH               January 1995 (expand comments)
    !CH               March 1995 (add Saha for Helium)
    !CH               August 1997 (add HeII alpha table)
    !CH               July 1998 (include OmegaT correction and H fudge factor)
    !CH               Nov 1998 (change Trad to Tmat in Rup)
    !CH               Jan 1999 (tidied up for public consumption)
    !CH               Sept 1999 (switch to formula for alpha's, fix glitch)
    !CH                  Sept 1999 modified to CMBFAST by US & MZ
    !CH                     Nov 1999 modified for F90 and CAMB (AML)
    !CH                     Aug 2000 modified to prevent overflow erorr in He_Boltz (AML)
    !CH                     Feb 2001 corrected fix of Aug 2000 (AML)
    !CH                     Oct 2001 fixed error in hubble parameter, now uses global function (AML)
    !                       March 2003 fixed bugs reported by savita gahlaut
    !                       March 2005 added option for corrections from astro-ph/0501672.
    !                                  thanks to V.K.Dubrovich, S.I.Grachev
    !                       June 2006 defined RECFAST_fudge as free parameter (AML)
    !                       October 2006 (included new value for G)
    !                       October 2006 (improved m_He/m_H to be "not4")
    !                       October 2006 (fixed error, x for x_H in part of f(1))
    !CH              January 2008 (improved HeI recombination effects,
    !CH                       including HeI rec. fudge factor)
    !                Feb 2008   Recfast 1.4 changes above added (AML)
    !                           removed Dubrovich option (wrong anyway)
    !CH              Sept 2008 (added extra term to make transition, smoother for Tmat evolution)
    !                Sept 2008 Recfast 1.4.2 changes above added (AML)
    !                          General recombination module structure, fix to make He x_e smooth also in recfast (AML)
    !CH      Jan 2010 (added fitting function to modify K
    !CH              to match x_e(z) for new H physics)
    !AL             June 2012 updated fudge parameters to match HyRec and CosmoRec (AML)
    !AL             Sept 2012 changes now in public recfast, version number changed to match Recfast 1.5.2.



    module Recombination
    use constants
    use classes
    use DarkAge21cm
    use MathUtils
    use MpiUtils, only : MpiStop
    use MiscUtils
    use RangeUtils
    use StringUtils
    use config
    use model
    use splines
    !use iso_c_binding
    implicit none
    public

    character(LEN=*), parameter :: Recfast_Version = 'Recfast_1.5.2'

    real(dl), parameter ::  zinitial = 1e4_dl !highest redshift
    real(dl), parameter ::  zfinal=0._dl
    integer,  parameter :: Nz=10000
    real(dl), parameter :: delta_z = (zinitial-zfinal)/Nz
    !real(dl), private :: zrec(Nz),xrec(Nz),dxrec(Nz), Tsrec(Nz) ,dTsrec(Nz), tmrec(Nz),dtmrec(Nz)
    integer, parameter :: RECFAST_Heswitch_default = 6
    real(dl), parameter :: RECFAST_fudge_He_default = 0.86_dl ! Helium fudge parameter
    logical, parameter  :: RECFAST_Hswitch_default = .true.    ! Include H corrections (v1.5, 2010)
    real(dl), parameter :: RECFAST_fudge_default = 1.14_dl     ! 1.14_dl
    real(dl), parameter :: RECFAST_fudge_default2 = 1.105d0 + 0.02d0

    logical, parameter :: evolve_Ts = .false. !local equilibrium is very accurate
    real(dl), parameter :: Do21cm_minev = 1/(1+400.) !at which to evolve T_s

    real(dl), parameter :: bigH=100.0D3/Mpc !Ho in s-1
    ! and
    real(dl), parameter :: sigma1 = sigma_thomson
    ! and
    real(dl), parameter :: not4  = mass_ratio_He_H    !mass He/H atom

    real(dl), parameter :: B01 = 3*B10
    !Fundamental constants in SI units
    !      ("not4" pointed out by Gary Steigman)

    real(dl), parameter ::  Lambda = 8.2245809d0
    real(dl), parameter :: Lambda_He = 51.3d0    !new value from Dalgarno
    real(dl), parameter :: L_H_ion   = 1.096787737D7 !level for H ion. (in m^-1)
    real(dl), parameter :: L_H_alpha = 8.225916453D6 !averaged over 2 levels
    real(dl), parameter :: L_He1_ion = 1.98310772D7  !from Drake (1993)
    real(dl), parameter :: L_He2_ion = 4.389088863D7 !from JPhysChemRefData (1987)
    real(dl), parameter :: L_He_2s   = 1.66277434D7  !from Drake (1993)
    real(dl), parameter :: L_He_2p   = 1.71134891D7  !from Drake (1993)
    !   2 photon rates and atomic levels in SI units

    real(dl), parameter :: A2P_s     = 1.798287D9    !Morton, Wu & Drake (2006)
    real(dl), parameter :: A2P_t     = 177.58D0      !Lach & Pachuski (2001)
    real(dl), parameter :: L_He_2Pt  = 1.690871466D7 !Drake & Morton (2007)
    real(dl), parameter :: L_He_2St  = 1.5985597526D7 !Drake & Morton (2007)
    real(dl), parameter :: L_He2St_ion  =3.8454693845D6 !Drake & Morton (2007)
    real(dl), parameter :: sigma_He_2Ps  = 1.436289D-22  !Hummer & Storey (1998)
    real(dl), parameter :: sigma_He_2Pt  = 1.484872D-22  !Hummer & Storey (1998)
    !    Atomic data for HeI

    ! andrea
    real(dl), parameter :: HeRayleighFac = 0.1_dl !Rayleigh neutral He scattering cross section as ratio to H 
    ! andrea

    !       Set up some constants so they don't have to be calculated later
    real(dl), parameter :: Lalpha = 1.d0/L_H_alpha
    real(dl), parameter :: Lalpha_He = 1.d0/L_He_2p
    real(dl), parameter :: DeltaB = h_P*C*(L_H_ion-L_H_alpha)
    real(dl), parameter :: CDB = DeltaB/k_B
    real(dl), parameter :: DeltaB_He = h_P*C*(L_He1_ion-L_He_2s)   !2s, not 2p
    real(dl), parameter :: CDB_He = DeltaB_He/k_B
    real(dl), parameter :: CB1 = h_P*C*L_H_ion/k_B
    real(dl), parameter :: CB1_He1 = h_P*C*L_He1_ion/k_B   !ionization for HeI
    real(dl), parameter :: CB1_He2 = h_P*C*L_He2_ion/k_B   !ionization for HeII
    real(dl), parameter :: CR = const_twopi*(m_e/h_P)*(k_B/h_P)
    real(dl), parameter :: CK = Lalpha**3/(const_eightpi)
    real(dl), parameter :: CK_He = Lalpha_He**3/(const_eightpi)
    real(dl), parameter :: CL = C*h_P/(k_B*Lalpha)
    real(dl), parameter :: CL_He = C*h_P/(k_B/L_He_2s) !comes from det.bal. of 2s-1s
    real(dl), parameter :: CT = Compton_CT / MPC_in_sec

    real(dl), parameter :: Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B

    !       Matter departs from radiation when t(Th) > H_frac * t(H)
    !       choose some safely small number
    real(dl), parameter :: H_frac = 1D-3

    procedure(obj_function), private :: dtauda
    
    ! andrea
    ! public TRecfast, CB1
    ! andrea 
    
    Type RecombinationData
        real(dl) :: Recombination_saha_z !Redshift at which saha OK
        real(dl), private :: NNow, fHe
        real(dl), private :: zrec(Nz)
        real(dl), private :: xrec(Nz),dxrec(Nz), Tsrec(Nz) ,dTsrec(Nz), tmrec(Nz),dtmrec(Nz)
        ! andrea
        real(dl) x_rayleigh_eff(Nz),dx_rayleigh_eff(Nz)
        ! andrea
        real(dl), private :: DeltaB,DeltaB_He,Lalpha,mu_H,mu_T
        real(dl), private :: HO, Tnow, fu
        integer, private :: n_eq = 3
        logical :: doTspin = .false.

        !The following only used for approximations where small effect
        real(dl) :: OmegaK, OmegaT, z_eq
        class(CAMBdata), pointer :: State
        contains
            !procedure :: Recombination_rayleigh_eff
            !procedure :: total_scattering_eff
    end Type RecombinationData

    type, extends(TRecombinationModel) :: TRecfast
        real(dl) :: RECFAST_fudge  = RECFAST_fudge_default2
        real(dl) :: RECFAST_fudge_He = RECFAST_fudge_He_default
        integer  :: RECFAST_Heswitch = RECFAST_Heswitch_default
        logical  :: RECFAST_Hswitch = RECFAST_Hswitch_default
        !0) no change from old Recfast'
        !1) full expression for escape probability for singlet'
        !'   1P-1S transition'
        !2) also including effect of contiuum opacity of H on HeI'
        !'   singlet (based in fitting formula suggested by'
        !'   Kholupenko, Ivanchik & Varshalovich, 2007)'
        !3) only including recombination through the triplets'
        !4) including 3 and the effect of the contiuum '
        !'   (although this is probably negligible)'
        !5) including only 1, 2 and 3'
        !6) including all of 1 to 4'

        !fudge parameter if RECFAST_Hswitch
        !Gaussian fits for extra H physics (fit by Adam Moss , modified by Antony Lewis)
        real(dl) :: AGauss1 =      -0.14D0  !Amplitude of 1st Gaussian
        real(dl) :: AGauss2 =       0.079D0 ! 0.05D0  !Amplitude of 2nd Gaussian
        real(dl) :: zGauss1 =       7.28D0  !ln(1+z) of 1st Gaussian
        real(dl) :: zGauss2=        6.73D0  !ln(1+z) of 2nd Gaussian
        real(dl) :: wGauss1=        0.18D0  !Width of 1st Gaussian
        real(dl) :: wGauss2=        0.33D0  !Width of 2nd Gaussian
        Type(RecombinationData), allocatable :: Calc
    contains
    procedure :: ReadParams => TRecfast_ReadParams
    procedure :: Validate => TRecfast_Validate
    procedure :: Init => TRecfast_init
    procedure :: x_e => TRecfast_xe
    procedure :: xe_Tm => TRecfast_xe_Tm !ionization fraction and baryon temperature
    procedure :: T_m => TRecfast_tm !baryon temperature
    procedure :: T_s => TRecfast_ts !Spin temperature
    procedure :: Version => TRecfast_version
    procedure :: dDeltaxe_dtau => TRecfast_dDeltaxe_dtau
    procedure :: get_Saha_z => TRecfast_Get_Saha_z
    procedure, nopass :: SelfPointer => TRecfast_SelfPointer
    procedure :: Recombination_rayleigh_eff
    procedure :: total_scattering_eff

    end type TRecfast

    type :: TBackgroundOutputs
        real(dl), allocatable :: H(:)
        real(dl), allocatable :: DA(:)
        real(dl), allocatable :: rs_by_D_v(:)
    end type TBackgroundOutputs

    integer, parameter :: derived_age=1, derived_zstar=2, derived_rstar=3, derived_thetastar=4, derived_DAstar = 5, &
        derived_zdrag=6, derived_rdrag=7,derived_kD=8,derived_thetaD=9, derived_zEQ =10, derived_keq =11, &
        derived_thetaEQ=12, derived_theta_rs_EQ = 13
    integer, parameter :: nthermo_derived = 13
    ! andrea
    integer, parameter :: num_cmb_freq = 6 !!!
    logical :: rayleigh_diff = .true.
    logical :: rayleigh_pows(3) = [.true.,.true.,.true.]
    logical :: rayleigh_back_approx = .false.
    integer, parameter :: nscatter = num_cmb_freq+1
    real(dl) :: phot_freqs(num_cmb_freq)  !set in equations _Init
    real(dl) :: phot_int_kernel(num_cmb_freq)
    real(dl) :: freq_factors(num_cmb_freq,3)
    real(dl) :: av_freq_factors(3) 
    real(dl) :: ALens = 1._dl
    ! andrea

    Type lSamples
        integer :: nl = 0
        integer :: lmin = 2
        integer, allocatable :: l(:)
        logical :: use_spline_template = .true.
    contains
    procedure :: Init => lSamples_init
    procedure :: IndexOf => lSamples_indexOf
    procedure :: InterpolateClArr
    procedure :: InterpolateClArrTemplated
    end Type lSamples

    logical, parameter :: dowinlens = .false. !not used, test getting CMB lensing using visibility
    ! andrea
    integer, parameter :: thermal_history_def_timesteps = 40000
    ! andrea

    Type TThermoData
        logical :: HasThermoData = .false. !Has it been computed yet for current parameters?
        !Background thermal history, interpolated from precomputed tables
        integer :: nthermo !Number of table steps
        !baryon temperature, sound speed, ionization fractions, and opacity
        ! andrea
        real(dl), dimension(:), allocatable :: tb, cs2, xe
        ! e^(-tau) and derivatives
        real(dl), dimension(:), allocatable :: dcs2
        real(dl), dimension(:,:), allocatable :: dotmu, ddotmu, sdotmu, emmu, demmu, dddotmu, ddddotmu ! BB24
        ! real(dl) dotmu(nthermo,nscatter), ddotmu(nthermo,nscatter)
        ! real(dl) sdotmu(nthermo,nscatter),emmu(nthermo,nscatter)
        ! real(dl) demmu(nthermo,nscatter)
        ! real(dl) dddotmu(nthermo,nscatter),ddddotmu(nthermo,nscatter)
        real(dl), dimension(:), allocatable :: ScaleFactor, dScaleFactor, adot, dadot
        real(dl), dimension(:), allocatable :: winlens, dwinlens
        real(dl) tauminn,dlntau
        real(dl) :: tight_tau, actual_opt_depth
        !Times when 1/(opacity*tau) = 0.01, for use switching tight coupling approximation
        real(dl) :: matter_verydom_tau
        real(dl) :: recombination_saha_tau
        !sound horizon and recombination redshifts
        real(dl) :: r_drag0, z_star, z_drag  !!JH for updated BAO likelihood.
        real(dl), dimension(:), allocatable :: step_redshift, rhos_fac, drhos_fac
        real(dl) :: tau_start_redshiftwindows,tau_end_redshiftwindows
        logical :: has_lensing_windows = .false.
        real(dl) recombination_Tgas_tau
        Type(TCubicSpline) :: ScaleFactorAtTime
        !Mapping between redshift and time
        real(dl), private,dimension(:), allocatable :: redshift_time, dredshift_time
        real(dl), private, dimension(:), allocatable :: arhos_fac, darhos_fac, ddarhos_fac
    contains
    procedure :: Init => Thermo_Init
    procedure :: OpacityToTime => Thermo_OpacityToTime
    !procedure :: total_scattering_eff => Thermo_Total_Scattering_Eff
    !procedure :: values => Thermo_values
    procedure :: values => Thermo_values 
    procedure :: values_array => Thermo_values_array ! BB24
    procedure :: expansion_values => Thermo_expansion_values
    procedure :: expansion_values_array => Thermo_expansion_values_array !BB24
    procedure :: IonizationFunctionsAtTime
    procedure, private :: DoWindowSpline
    procedure, private :: SetTimeSteps
    procedure, private :: SetTimeStepWindows
    end type TThermoData

    !Sources
    Type CalWins
        real(dl), allocatable :: awin_lens(:), dawin_lens(:)
    end Type CalWins

    Type LimberRec
        integer n1,n2 !corresponding time step array indices
        real(dl), dimension(:), allocatable :: k
        real(dl), dimension(:), allocatable :: Source
    end Type LimberRec

    Type ClTransferData
        !Cl transfer function variables
        !values of q for integration over q to get C_ls
        Type (lSamples) :: ls ! l that are computed
        integer :: NumSources
        !Changes -scalars:  2 for just CMB, 3 for lensing
        !- tensors: T and E and phi (for lensing), and T, E, B respectively

        type (TRanges) :: q
        real(dl), dimension(:,:,:), allocatable :: Delta_p_l_k

        !The L index of the lowest L to use for Limber
        integer, dimension(:), allocatable :: Limber_l_min
        !For each l, the set of k in each limber window
        !indices LimberWindow(SourceNum,l)
        Type(LimberRec), dimension(:,:), allocatable :: Limber_windows

        !The maximum L needed for non-Limber
        integer max_index_nonlimber

    end Type ClTransferData

    type TCLdata
        Type(ClTransferData) :: CTransScal, CTransTens, CTransVec

        real(dl), dimension (:,:), allocatable :: Cl_scalar, Cl_tensor, Cl_vector
        !Indices are Cl_xxx( l , Cl_type)
        !where Cl_type is one of the above constants

        real(dl), dimension (:,:,:), allocatable :: Cl_Scalar_Array
        !Indices are Cl_xxx( l , field1,field2)
        !where ordering of fields is T, E, \psi (CMB lensing potential),T_freq_1,E_freq_1,.... window_1, window_2...

        !The following are set only if doing lensing
        integer lmax_lensed !Only accurate to rather less than this
        real(dl) , dimension (:,:), allocatable :: Cl_lensed
        !Cl_lensed(l, Cl_type) are the interpolated Cls
        ! andrea
        real(dl) , dimension (:,:,:,:,:), allocatable :: Cl_lensed_freqs, Cl_tensor_freqs
        ! andrea
    contains
    procedure :: InitCls => TCLdata_InitCls
    procedure :: output_cl_files => TCLdata_output_cl_files
    procedure :: output_lens_pot_files => TCLdata_output_lens_pot_files
    procedure :: NormalizeClsAtL => TCLdata_NormalizeClsAtL
    procedure :: output_veccl_files => TCLdata_output_veccl_files
    end type TCLdata

    Type TTimeSources
        ! values of q to evolve the propagation equations to compute the sources
        type(TRanges) :: Evolve_q
        real(dl), dimension(:,:,:), allocatable :: LinearSrc !Sources and second derivs
        !LinearSrc indices  ( k_index, source_index, time_step_index )
        integer SourceNum, NonCustomSourceNum
        !SourceNum is total number sources (2 or 3 for scalars, 3 for tensors).
    end type TTimeSources

    type, extends(TCAMBdata) :: CAMBdata

        type(CAMBparams) :: CP

        real(dl) ThermoDerivedParams(nthermo_derived)

        logical flat,closed

        !     grhocrit =kappa*a^2*rho_crit(0)
        !     grhornomass=grhor*number of massless neutrino species
        !     taurst,taurend - time at start/end of recombination
        !     dtaurec - dtau during recombination
        !     adotrad - a(tau) in radiation era
        real(dl) grhocrit,grhog,grhor,grhob,grhoc,grhov,grhornomass,grhok
        real(dl) taurst,dtaurec,taurend,tau_maxvis,adotrad

        real(dl) Omega_de
        real(dl) curv, curvature_radius, Ksign !curvature_radius = 1/sqrt(|curv|), Ksign = 1,0 or -1
        real(dl) tau0,chi0 !time today and rofChi(tau0/curvature_radius)
        real(dl) scale !relative to flat. e.g. for scaling lSamp%l sampling.

        real(dl) akthom !sigma_T * (number density of protons now)
        real(dl) fHe !n_He_tot / n_H_tot
        real(dl) Nnow
        real(dl) z_eq !defined assuming all neutrinos massless
        !Neutrinos
        real(dl) grhormass(max_nu)
        !     nu_masses=m_nu*c**2/(k_B*T_nu0)
        real(dl) nu_masses(max_nu)
        integer ::  num_transfer_redshifts = 1
        real(dl), allocatable  ::  transfer_redshifts(:)
        integer  ::  PK_redshifts_index(max_transfer_redshifts)

        logical :: OnlyTransfer = .false. !C_L/PK not computed; initial power spectrum data, instead get Delta_q_l array
        !If true, sigma_8 is not calculated either]]
        logical :: HasScalarTimeSources = .false. !No power spectra, only time transfer functions

        logical :: get_growth_sigma8 = .true.
        !gets sigma_vdelta, like sigma8 but using velocity-density cross power,
        !in late LCDM f*sigma8 = sigma_vdelta^2/sigma8

        logical :: needs_good_pk_sampling = .false.

        logical ::call_again = .false.
        !if being called again with same parameters to get different thing

        real(dl) :: reion_tau_start, reion_tau_complete
        integer :: reion_n_steps

        Type(TNuPerturbations) :: NuPerturbations

        Type(TBackgroundOutputs) :: BackgroundOutputs

        !Time steps for sampling sources
        Type(TRanges) :: TimeSteps
        !Background interpolation tables for thermal history etc.
        Type(TThermoData) :: ThermoData

        real(dl), allocatable :: transfer_times(:)

        !Matter transfer data
        Type (MatterTransferData):: MT

        !Matter power spectrum for default variable (used for non-linear corrections)
        Type(MatterPowerData), allocatable :: CAMB_PK

        Type(TClData) :: CLdata

        integer :: num_redshiftwindows = 0
        integer :: num_extra_redshiftwindows = 0
        Type(TRedWin), allocatable :: Redshift_W(:)
        real(dl), dimension(:), allocatable :: optical_depths_for21cm

        Type(TTimeSources), allocatable :: ScalarTimeSources
        integer :: Scalar_C_last = C_PhiE
    contains
    procedure :: DeltaTime => CAMBdata_DeltaTime
    procedure :: DeltaTimeArr => CAMBdata_DeltaTimeArr
    procedure :: TimeOfz => CAMBdata_TimeOfz
    procedure :: TimeOfzArr => CAMBdata_TimeOfzArr
    procedure :: DeltaPhysicalTimeGyr => CAMBdata_DeltaPhysicalTimeGyr
    procedure :: DeltaPhysicalTimeGyrArr => CAMBdata_DeltaPhysicalTimeGyrArr
    procedure :: AngularDiameterDistance => CAMBdata_AngularDiameterDistance
    procedure :: AngularDiameterDistanceArr => CAMBdata_AngularDiameterDistanceArr
    procedure :: AngularDiameterDistance2 => CAMBdata_AngularDiameterDistance2
    procedure :: AngularDiameterDistance2Arr => CAMBdata_AngularDiameterDistance2Arr
    procedure :: LuminosityDistance => CAMBdata_LuminosityDistance
    procedure :: ComovingRadialDistance => CAMBdata_ComovingRadialDistance
    procedure :: ComovingRadialDistanceArr => CAMBdata_ComovingRadialDistanceArr
    procedure :: GetBackgroundDensities => CAMBdata_GetBackgroundDensities
    procedure :: Hofz => CAMBdata_Hofz
    procedure :: HofzArr => CAMBdata_HofzArr
    procedure :: sound_horizon => CAMBdata_sound_horizon
    procedure :: sound_horizon_zArr => CAMBdata_sound_horizon_zArr
    procedure :: RedshiftAtTimeArr => CAMBdata_RedshiftAtTimeArr
    procedure :: BAO_D_v => CAMBdata_BAO_D_v
    procedure :: CosmomcTheta => CAMBdata_CosmomcTheta
    procedure :: get_lmax_lensed => CAMBdata_get_lmax_lensed
    procedure :: get_zstar => CAMBdata_get_zstar
    procedure :: DarkEnergyStressEnergy => CAMBdata_DarkEnergyStressEnergy
    procedure :: SetParams => CAMBdata_SetParams
    procedure :: Free => CAMBdata_Free
    procedure :: grho_no_de
    procedure :: GetReionizationOptDepth
    procedure :: rofChi
    procedure :: cosfunc
    procedure :: tanfunc
    procedure :: invsinfunc
    procedure :: GetComputedPKRedshifts
    procedure :: binary_search
    procedure, nopass :: PythonClass => CAMBdata_PythonClass
    procedure, nopass :: SelfPointer => CAMBdata_SelfPointer
    end type CAMBdata

    interface
    FUNCTION state_function(obj, a)
    use precision
    import
    class(CAMBdata) :: obj
    real(dl), intent(in) :: a
    real(dl) :: state_function
    END FUNCTION  state_function
    end interface

    !procedure(obj_function), private :: dtauda

    ! and
    !interface
    !subroutine Thermo_values(this, tau, a, cs2, opacity)
    !    use precision
    !    type(ThermoData), intent(in) :: this
    !    real(dl), intent(in) :: tau
    !    real(dl), intent(out) :: a
    !    real(dl), intent(out) :: cs2
    !    real(dl), intent(out) :: opacity
    !    integer i
    !    real(dl) d
    !end subroutine Thermo_values
    !end interface
    ! and

    contains

    function CAMBdata_PythonClass()
    character(LEN=:), allocatable :: CAMBdata_PythonClass
    CAMBdata_PythonClass = 'CAMBdata'
    end function CAMBdata_PythonClass

    !function CAMBdata_PythonClass() result(res)
    !    character(len=:), allocatable :: res
    !    res = 'CAMBdata'
    !end function CAMBdata_PythonClass

    subroutine CAMBdata_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (CAMBdata), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine CAMBdata_SelfPointer

    subroutine CAMBdata_SetParams(this, P, error, DoReion, call_again, background_only)
    !Initialize background variables; does not yet calculate thermal history
    use constants
    class(CAMBdata), target :: this
    type(CAMBparams), intent(in) :: P
    real(dl) fractional_number, conv
    integer, optional :: error !Zero if OK
    logical, optional :: DoReion
    logical, optional :: call_again, background_only
    logical WantReion, calling_again
    integer nu_i,actual_massless
    real(dl) nu_massless_degeneracy, neff_i, eta_k, h2
    real(dl) zpeak, sigma_z, zpeakstart, zpeakend
    Type(TRedWin), pointer :: Win
    logical back_only
    !Constants in SI units

    global_error_flag = 0

    if ((P%WantTensors .or. P%WantVectors).and. P%WantTransfer .and. .not. P%WantScalars) then
        call GlobalError( 'Cannot generate tensor C_l and transfer without scalar C_l',error_unsupported_params)
    end if

    if (.not. allocated(P%DarkEnergy)) then
        call GlobalError('DarkEnergy not set', error_darkenergy)
    end if

    if (present(error)) error = global_error_flag
    if (global_error_flag/=0) return

    WantReion = DefaultTrue(DoReion)
    calling_again= DefaultFalse(call_again)
    back_only = DefaultFalse(background_only)

    if (calling_again) then
        this%CP%WantDerivedParameters = .false.
        this%CP%Accuracy = P%Accuracy
        this%CP%Transfer%high_precision = P%Transfer%high_precision
        this%CP%Reion%Reionization = P%Reion%Reionization
        this%CP%WantTransfer =P%WantTransfer
        this%CP%WantScalars =P%WantScalars
        this%CP%WantTensors =P%WantTensors
        this%CP%WantVectors =P%WantVectors
        this%CP%WantCls = P%WantCls
    else
        this%CP=P
        this%CP%Max_eta_k = max(this%CP%Max_eta_k,this%CP%Max_eta_k_tensor)
    end if

    if (P%WantTransfer .and. .not. back_only) then
        this%CP%WantScalars=.true.
        if (.not. P%WantCls) then
            this%CP%Accuracy%AccuratePolarization = .false.
            !Sources
            this%CP%Reion%Reionization = this%CP%transfer_21cm_cl
        end if
        call this%GetComputedPKRedshifts(this%CP)
    end if
    if (this%CP%WantTransfer.and. this%CP%MassiveNuMethod==Nu_approx) then
        this%CP%MassiveNuMethod = Nu_trunc
    end if

    if (.not. calling_again) then
        this%ThermoData%HasThermoData = .false.
        if (this%CP%Num_Nu_Massive /= sum(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates))) then
            if (sum(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates))/=0) &
                call GlobalError('Num_Nu_Massive is not sum of Nu_mass_numbers', error_unsupported_params)
        end if
10      if (this%CP%Omnuh2 < 1.e-7_dl) this%CP%Omnuh2 = 0
        if (this%CP%Omnuh2==0 .and. this%CP%Num_Nu_Massive /=0) then
            if (this%CP%share_delta_neff) then
                this%CP%Num_Nu_Massless = this%CP%Num_Nu_Massless + this%CP%Num_Nu_Massive
            else
                this%CP%Num_Nu_Massless = this%CP%Num_Nu_Massless + sum(this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates))
            end if
            this%CP%Num_Nu_Massive  = 0
            this%CP%Nu_mass_numbers = 0
        end if

        nu_massless_degeneracy = this%CP%Num_Nu_massless !N_eff for massless neutrinos
        if (this%CP%Num_nu_massive > 0) then
            if (this%CP%Nu_mass_eigenstates==0) &
                call GlobalError('Have Num_nu_massive>0 but no nu_mass_eigenstates', error_unsupported_params)
            if (this%CP%Nu_mass_eigenstates==1 .and. this%CP%Nu_mass_numbers(1)==0) &
                this%CP%Nu_mass_numbers(1) = this%CP%Num_Nu_Massive
            if (all(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates)==0)) this%CP%Nu_mass_numbers=1 !just assume one for all
            if (this%CP%share_delta_neff) then
                !default case of equal heating of all neutrinos
                fractional_number = this%CP%Num_Nu_massless + this%CP%Num_Nu_massive
                actual_massless = int(this%CP%Num_Nu_massless + 1e-6_dl)
                neff_i = fractional_number/(actual_massless + this%CP%Num_Nu_massive)
                nu_massless_degeneracy = neff_i*actual_massless
                this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates) = &
                    this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates)*neff_i
            end if
            if (abs(sum(this%CP%Nu_mass_fractions(1:this%CP%Nu_mass_eigenstates))-1) > 1e-4) &
                call GlobalError('Nu_mass_fractions do not add up to 1', error_unsupported_params)
        else
            this%CP%Nu_mass_eigenstates = 0
        end if

        if (global_error_flag/=0) then
            if (present(error)) error = global_error_flag
            return
        end if


        call This%ThermoData%ScaleFactorAtTime%Clear()

        this%flat = (abs(this%CP%omk) <= OmegaKFlat)
        this%closed = this%CP%omk < -OmegaKFlat

        if (this%flat) then
            this%curv=0
            this%Ksign=0
            this%curvature_radius=1._dl !so we can use tau/curvature_radius, etc, where r's cancel
        else
            this%curv=-this%CP%omk/((c/1000)/this%CP%h0)**2
            this%Ksign =sign(1._dl,this%curv)
            this%curvature_radius=1._dl/sqrt(abs(this%curv))
        end if
        !  grho gives the contribution to the expansion rate from: (g) photons,
        !  (r) one flavor of relativistic neutrino (2 degrees of freedom),
        !  (m) nonrelativistic matter (for Omega=1).  grho is actually
        !  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
        !  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
        !  (Used only to set the initial conformal time.)

        !H0 is in km/s/Mpc

        this%grhocrit = 3*this%CP%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)

        this%grhog = kappa/c**2*4*sigma_boltz/c**3*this%CP%tcmb**4*Mpc**2 !8*pi*G/c^2*4*sigma_B/c^3 T^4
        ! grhog=1.4952d-13*tcmb**4
        this%grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*this%grhog !7/8*(4/11)^(4/3)*grhog (per neutrino species)
        !grhor=3.3957d-14*tcmb**4

        !correction for fractional number of neutrinos, e.g. 3.04 to give slightly higher T_nu hence rhor
        !for massive Nu_mass_degeneracies parameters account for heating from grhor

        this%grhornomass=this%grhor*nu_massless_degeneracy
        this%grhormass=0
        do nu_i = 1, this%CP%Nu_mass_eigenstates
            this%grhormass(nu_i)=this%grhor*this%CP%Nu_mass_degeneracies(nu_i)
        end do
        h2 = (this%CP%H0/100)**2
        this%grhoc=this%grhocrit*this%CP%omch2/h2
        this%grhob=this%grhocrit*this%CP%ombh2/h2
        this%grhok=this%grhocrit*this%CP%omk
        this%Omega_de = 1 -(this%CP%omch2 + this%CP%ombh2 + this%CP%omnuh2)/h2 - this%CP%omk  &
            - (this%grhornomass + this%grhog)/this%grhocrit
        this%grhov=this%grhocrit*this%Omega_de

        !  adotrad gives da/dtau in the asymptotic radiation-dominated era:
        this%adotrad = sqrt((this%grhog+this%grhornomass+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates)))/3)

        this%Nnow = this%CP%ombh2/h2*(1-this%CP%yhe)*this%grhocrit*c**2/kappa/m_H/Mpc**2

        this%akthom = sigma_thomson*this%Nnow*Mpc
        !sigma_T * (number density of protons now)

        this%fHe = this%CP%YHe/(mass_ratio_He_H*(1.d0-this%CP%YHe))  !n_He_tot / n_H_tot

        this%z_eq = (this%grhob+this%grhoc)/(this%grhog+this%grhornomass+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates))) -1

        if (this%CP%omnuh2/=0) then
            !Initialize things for massive neutrinos
            call ThermalNuBackground%Init()
            call this%NuPerturbations%Init(P%Accuracy%AccuracyBoost*P%Accuracy%neutrino_q_boost)
            !  nu_masses=m_nu(i)*c**2/(k_B*T_nu0)
            do nu_i=1, this%CP%Nu_mass_eigenstates
                this%nu_masses(nu_i)= ThermalNuBackground%find_nu_mass_for_rho(this%CP%omnuh2/h2*this%CP%Nu_mass_fractions(nu_i)&
                    *this%grhocrit/this%grhormass(nu_i))
            end do
            if (all(this%nu_masses(1:this%CP%Nu_mass_eigenstates)==0)) then
                !All density accounted for by massless, so just use massless
                this%CP%Omnuh2 = 0
                goto 10
            end if
            !Just prevent divide by zero
            this%nu_masses(1:this%CP%Nu_mass_eigenstates) = max(this%nu_masses(1:this%CP%Nu_mass_eigenstates),1e-3_dl)
        else
            this%nu_masses = 0
        end if
        call this%CP%DarkEnergy%Init(this)
        if (global_error_flag==0) this%tau0=this%TimeOfz(0._dl)
        if (global_error_flag==0) then
            this%chi0=this%rofChi(this%tau0/this%curvature_radius)
            this%scale= this%chi0*this%curvature_radius/this%tau0  !e.g. change l sampling depending on approx peak spacing
            if (this%closed .and. this%tau0/this%curvature_radius >3.14) then
                call GlobalError('chi >= pi in closed model not supported',error_unsupported_params)
            end if
            if (WantReion) call this%CP%Reion%Init(this)
            if (this%CP%NonLinear/=NonLinear_None .and. .not. back_only) &
                call this%CP%NonLinearModel%Init(this)
        end if
    end if
    if (allocated(this%CP%SourceWindows) .and. .not. back_only) then
        if (.not. this%CP%WantScalars) then
            this%num_redshiftwindows=0
        else
            this%num_redshiftwindows = size(this%CP%SourceWindows)
        end if
    else
        this%num_redshiftwindows = 0
        this%CP%SourceTerms%limber_windows = .false.
    endif

    if (this%CP%WantScalars .and. this%CP%WantCls .and. this%num_redshiftwindows>0) then
        eta_k = this%CP%Max_eta_k
        if (allocated(this%Redshift_W)) deallocate(this%Redshift_W)
        allocate(this%Redshift_W(this%num_redshiftwindows))
        this%num_extra_redshiftwindows = 0
        !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(zpeak, sigma_z, zpeakstart, zpeakend, nu_i, Win)
        do nu_i = 1, this%num_redshiftwindows
            Win => this%Redshift_w(nu_i)
            Win%Window => this%CP%SourceWindows(nu_i)%Window
            Win%kind = Win%Window%source_type
            call Win%Window%GetScales(zpeak, sigma_z, zpeakstart, zpeakend)
            if (FeedbackLevel > 1) then
                write(*,*) FormatString('Window scales: %d peak: %f, sigma: %f, start:%f, end %f', &
                    nu_i, zpeak, sigma_z, zpeakstart, zpeakend)
            end if
            Win%Redshift = zpeak
            Win%tau = this%TimeOfz(zpeak, tol=1e-4_dl)
            Win%sigma_tau = sigma_z*dtauda(this,1/(1+zpeak))/(1+zpeak)**2
            Win%tau_peakstart=this%TimeOfZ(zpeakstart, tol=1e-4_dl)
            Win%tau_peakend = this%TimeOfZ(max(0._dl,zpeakend), tol=1e-4_dl)
            Win%chi0 = this%tau0-Win%tau
            Win%chimin = min(Win%chi0,this%tau0 - this%TimeOfz(max(0.05_dl,zpeakend), tol=1e-4_dl))
            !$OMP CRITICAL
            this%CP%Max_eta_k = max(this%CP%Max_eta_k, this%tau0*WindowKmaxForL(Win,this%CP,this%CP%max_l))
            if (Win%Window%source_type==window_21cm) this%CP%Do21cm = .true.
            if (Win%Window%source_type==window_counts .and. P%SourceTerms%counts_lensing) then
                this%num_extra_redshiftwindows = this%num_extra_redshiftwindows + 1
                Win%mag_index = this%num_extra_redshiftwindows
            end if
            !$OMP END CRITICAL
        end do
        if (eta_k /= this%CP%Max_eta_k .and. FeedbackLevel>0) &
            write (*,*) 'source max_eta_k: ', this%CP%Max_eta_k,'kmax = ', this%CP%Max_eta_k/this%tau0
    end if

    if ((this%CP%NonLinear==NonLinear_Lens .or. this%CP%NonLinear==NonLinear_both) .and. &
        this%CP%Max_eta_k/this%tau0 > this%CP%Transfer%kmax) then
        this%CP%Transfer%kmax =this%CP%Max_eta_k/this%tau0
        if (FeedbackLevel > 0) write (*,*) 'kmax changed to ', this%CP%Transfer%kmax
    end if

    if (global_error_flag/=0) then
        if (present(error)) error = global_error_flag
        return
    end if

    if (present(error)) then
        error = 0
    else if (FeedbackLevel > 0 .and. .not. calling_again) then
        write(*,'("Om_b h^2             = ",f9.6)') P%ombh2
        write(*,'("Om_c h^2             = ",f9.6)') P%omch2
        write(*,'("Om_nu h^2            = ",f9.6)') P%omnuh2
        write(*,'("Om_darkenergy        = ",f9.6)') this%Omega_de
        write(*,'("Om_K                 = ",f9.6)') P%omk
        write(*,'("Om_m (inc Om_u)      = ",f9.6)') (P%ombh2+P%omch2+P%omnuh2)/h2
        write(*,'("100 theta (CosmoMC)  = ",f9.6)') 100*this%CosmomcTheta()
        if (this%CP%Num_Nu_Massive > 0) then
            write(*,'("N_eff (total)        = ",f9.6)') nu_massless_degeneracy + &
                sum(this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates))
            do nu_i=1, this%CP%Nu_mass_eigenstates
                conv = k_B*(8*this%grhor/this%grhog/7)**0.25*this%CP%tcmb/eV1 * &
                    (this%CP%nu_mass_degeneracies(nu_i)/this%CP%nu_mass_numbers(nu_i))**0.25 !approx 1.68e-4
                write(*,'(I2, " nu, g=",f7.4," m_nu*c^2/k_B/T_nu0= ",f9.2," (m_nu= ",f6.3," eV)")') &
                    this%CP%nu_mass_numbers(nu_i), this%CP%nu_mass_degeneracies(nu_i), &
                    this%nu_masses(nu_i),conv*this%nu_masses(nu_i)
            end do
        end if
    end if

    end subroutine CAMBdata_SetParams

    subroutine CAMBdata_Free(this)
    class(CAMBdata) :: this

    call Free_ClTransfer(this%CLdata%CTransScal)
    call Free_ClTransfer(this%ClData%CTransVec)
    call Free_ClTransfer(this%ClData%CTransTens)
    call this%MT%Free()
    if (allocated(this%CAMB_Pk)) deallocate(this%CAMB_PK)

    end subroutine CAMBdata_Free

    function CAMBdata_DeltaTime(this, a1,a2, in_tol)
    class(CAMBdata) :: this
    real(dl) CAMBdata_DeltaTime, atol
    real(dl), intent(IN) :: a1,a2
    real(dl), optional, intent(in) :: in_tol

    atol = PresentDefault(tol/1000/exp(this%CP%Accuracy%AccuracyBoost*this%CP%Accuracy%IntTolBoost-1), in_tol)
    CAMBdata_DeltaTime = Integrate_Romberg(this, dtauda,a1,a2,atol)

    end function CAMBdata_DeltaTime

    ! and
    ! andrea
    function Recombination_rayleigh_eff(this,a)
    class(TRecfast) :: this
    real(dl), intent(in) :: a
    real(dl) zst,z,az,bz,Recombination_rayleigh_eff
    integer ilo,ihi
    z=1/a-1
    associate(Calc => this%Calc)
    if (z.ge.Calc%zrec(1)) then
        !break Recombination_rayleigh_eff
        Recombination_rayleigh_eff=Calc%x_rayleigh_eff(1)
    else
        if (z.le.Calc%zrec(nz)) then
            Recombination_rayleigh_eff=Calc%x_rayleigh_eff(nz)
        else
            zst=(zinitial-z)/delta_z
            ihi= int(zst)
            ilo = ihi+1
            az=zst - int(zst)
            bz=1-az     
            Recombination_rayleigh_eff=az*Calc%x_rayleigh_eff(ilo)+bz*Calc%x_rayleigh_eff(ihi)+ &
            ((az**3-az)*Calc%dx_rayleigh_eff(ilo)+(bz**3-bz)*Calc%dx_rayleigh_eff(ihi))/6._dl
        endif
    endif
    end associate
    end function Recombination_rayleigh_eff

    function total_scattering_eff(this, State, a)
    class(TRecfast) :: this
    class(CAMBdata), intent(in) :: State
    real(dl), intent(in) :: a
    real(dl) :: a2, total_scattering_eff

    associate(Calc => this%Calc)
        if (rayleigh_back_approx) then
            a2 = a**2
            total_scattering_eff = this%x_e(a) + this%recombination_rayleigh_eff(a) * ( &
                min(1._dl, av_freq_factors(1) / a2**2 + av_freq_factors(2) / a2**3 + av_freq_factors(3) / a2**4) )
        else
            total_scattering_eff = State%CP%Recomb%x_e(a)
        end if
    end associate
    end function total_scattering_eff
    ! andrea

    subroutine CAMBdata_DeltaTimeArr(this, arr, a1, a2, n, tol)
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: a1(n), a2(n)
    real(dl), intent(in), optional :: tol
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        arr(i) = this%DeltaTime(a1(i), a2(i), tol)
    end do

    end subroutine CAMBdata_DeltaTimeArr

    function CAMBdata_TimeOfz(this, z, tol)
    class(CAMBdata) :: this
    real(dl) CAMBdata_TimeOfz
    real(dl), intent(in), optional :: tol
    real(dl), intent(IN) :: z

    CAMBdata_TimeOfz= this%DeltaTime(0._dl,1._dl/(z+1._dl), tol)
    end function CAMBdata_TimeOfz

    subroutine CAMBdata_TimeOfzArr(this, arr, z, n, tol)
    !z array must be monotonically *decreasing* so times increasing
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    real(dl), intent(in), optional :: tol
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        if (i==1) then
            arr(i) = this%DeltaTime(0._dl, 1/(1+z(1)), tol)
        else
            if (z(i) < z(i-1)) then
                arr(i) = this%DeltaTime(1/(1+z(i-1)),1/(1+z(i)),tol)
            elseif (z(i) < z(i-1) + 1e-6_dl) then
                arr(i)=0
            else
                error stop 'CAMBdata_TimeOfzArr redshifts must be decreasing'
            end if
        end if
    end do
    !$OMP END PARALLEL DO
    do i = 2, n
        arr(i) = arr(i)  + arr(i-1)
    end do

    end subroutine CAMBdata_TimeOfzArr

    function CAMBdata_DeltaPhysicalTimeGyr(this, a1,a2, in_tol)
    use constants
    class(CAMBdata) :: this
    real(dl), intent(in) :: a1, a2
    real(dl), optional, intent(in) :: in_tol
    real(dl) CAMBdata_DeltaPhysicalTimeGyr, atol

    atol = PresentDefault(1d-4/exp(this%CP%Accuracy%AccuracyBoost-1), in_tol)
    CAMBdata_DeltaPhysicalTimeGyr = Integrate_Romberg(this, dtda,a1,a2,atol)*Mpc/c/Gyr
    end function CAMBdata_DeltaPhysicalTimeGyr

    subroutine CAMBdata_DeltaPhysicalTimeGyrArr(this, arr, a1, a2, n, tol)
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: a1(n), a2(n)
    real(dl), intent(in), optional :: tol
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        arr(i) = this%DeltaPhysicalTimeGyr(a1(i), a2(i), tol)
    end do

    end subroutine CAMBdata_DeltaPhysicalTimeGyrArr


    function CAMBdata_AngularDiameterDistance(this,z)
    class(CAMBdata) :: this
    !This is the physical (non-comoving) angular diameter distance in Mpc
    real(dl) CAMBdata_AngularDiameterDistance
    real(dl), intent(in) :: z

    CAMBdata_AngularDiameterDistance = this%curvature_radius/(1+z)* &
        this%rofchi(this%ComovingRadialDistance(z) /this%curvature_radius)

    end function CAMBdata_AngularDiameterDistance

    subroutine CAMBdata_AngularDiameterDistanceArr(this, arr, z, n)
    class(CAMBdata) :: this
    !This is the physical (non-comoving) angular diameter distance in Mpc for array of z
    !z array must be monotonically increasing
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i

    call this%ComovingRadialDistanceArr(arr, z, n, 1e-4_dl)
    if (this%flat) then
        arr = arr/(1+z)
    else
        do i=1, n
            arr(i) =  this%curvature_radius/(1+z(i))*this%rofchi(arr(i)/this%curvature_radius)
        end do
    end if

    end subroutine CAMBdata_AngularDiameterDistanceArr


    function CAMBdata_AngularDiameterDistance2(this,z1, z2)
    ! z1 < z2, otherwise returns zero
    !From http://www.slac.stanford.edu/~amantz/work/fgas14/#cosmomc
    class(CAMBdata) :: this
    real(dl) CAMBdata_AngularDiameterDistance2
    real(dl), intent(in) :: z1, z2

    if (z2 < z1 + 1e-4) then
        CAMBdata_AngularDiameterDistance2=0
    else
        CAMBdata_AngularDiameterDistance2 = this%curvature_radius/(1+z2)* &
            this%rofchi( this%DeltaTime(1/(1+z2),1/(1+z1))/this%curvature_radius)
    end if

    end function CAMBdata_AngularDiameterDistance2

    subroutine CAMBdata_AngularDiameterDistance2Arr(this, arr, z1, z2, n)
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z1(n), z2(n)
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        arr(i) = this%AngularDiameterDistance2(z1(i),z2(i))
    end do
    !$OMP END PARALLEL DO

    end subroutine CAMBdata_AngularDiameterDistance2Arr


    function CAMBdata_LuminosityDistance(this,z)
    class(CAMBdata) :: this
    real(dl) CAMBdata_LuminosityDistance
    real(dl), intent(in) :: z

    CAMBdata_LuminosityDistance = this%AngularDiameterDistance(z)*(1+z)**2

    end function CAMBdata_LuminosityDistance

    function CAMBdata_ComovingRadialDistance(this, z)
    class(CAMBdata) :: this
    real(dl) CAMBdata_ComovingRadialDistance
    real(dl), intent(in) :: z

    CAMBdata_ComovingRadialDistance = this%DeltaTime(1/(1+z),1._dl)

    end function CAMBdata_ComovingRadialDistance

    subroutine CAMBdata_ComovingRadialDistanceArr(this, arr, z, n, tol)
    !z array must be monotonically increasing
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    real(dl), intent(in) :: tol
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
    do i = 1, n
        if (i==1) then
            if (z(i) < 1e-6_dl) then
                arr(i) = 0
            else
                arr(i) = this%DeltaTime(1/(1+z(i)),1._dl, tol)
            end if
        else
            if (z(i) < z(i-1)) error stop 'ComovingRadialDistanceArr redshifts out of order'
            arr(i) = this%DeltaTime(1/(1+z(i)),1/(1+z(i-1)),tol)
        end if
    end do
    !$OMP END PARALLEL DO
    do i = 2, n
        arr(i) = arr(i)  + arr(i-1)
    end do

    end subroutine CAMBdata_ComovingRadialDistanceArr

    function CAMBdata_Hofz(this,z)
    !non-comoving Hubble in MPC units, divide by MPC_in_sec to get in SI units
    !multiply by c/1e3 to get in km/s/Mpc units
    class(CAMBdata) :: this
    real(dl) CAMBdata_Hofz, a
    real(dl), intent(in) :: z

    a = 1/(1+z)
    CAMBdata_Hofz = 1/(a**2*dtauda(this,a))

    end function CAMBdata_Hofz

    subroutine CAMBdata_HofzArr(this, arr, z, n)
    !non-comoving Hubble in MPC units, divide by MPC_in_sec to get in SI units
    !multiply by c/1e3 to get in km/s/Mpc units
    class(CAMBdata) :: this
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i
    real(dl) :: a

    do i=1, n
        a = 1/(1+z(i))
        arr(i) = 1/(a**2*dtauda(this,a))
    end do

    end subroutine CAMBdata_HofzArr

    real(dl) function CAMBdata_sound_horizon(this, z)
    class(CAMBdata) :: this
    real(dl), intent(in) :: z

    CAMBdata_sound_horizon = Integrate_Romberg(this,dsound_da_exact,1d-9,1/(z+1),1e-6_dl)

    end function CAMBdata_sound_horizon

    subroutine CAMBdata_sound_horizon_zArr(this,arr, z,n)
    class(CAMBdata) :: this
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: z(n)
    integer i

    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), IF(n>4)
    do i=1,n
        arr(i) = this%sound_horizon(z(i))
    end do

    end subroutine CAMBdata_sound_horizon_zArr

    subroutine CAMBdata_RedshiftAtTimeArr(this, arr, tau, n)
    class(CAMBdata) :: this
    integer,intent(in) :: n
    real(dl), intent(out) :: arr(n)
    real(dl), intent(in) :: tau(n)
    integer i
    real(dl) om

    if (this%ThermoData%ScaleFactorAtTime%n==0) &
        call GlobalError('RedshiftAtTimeArr: background history not calculated', error_unsupported_params)
    if (global_error_flag/=0) return
    !$OMP PARALLEL DO DEFAULT(SHARED), private(om, i)
    do i=1, n
        if (tau(i) < this%ThermoData%tauminn*1.1) then
            om = (this%grhob+this%grhoc)/&
                sqrt(3*(this%grhog+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates))+this%grhornomass))
            arr(i) = 1/(this%adotrad*tau(i)*(1+om*tau(i)/4))-1
        else
            arr(i) = 1/this%ThermoData%ScaleFactorAtTime%Value(tau(i))-1
        end if
    end do

    end subroutine CAMBdata_RedshiftAtTimeArr

    real(dl) function BAO_D_v_from_DA_H(z, DA, Hz)
    real(dl), intent(in) :: z, DA, Hz
    real(dl) ADD

    ADD = DA*(1.d0+z)
    BAO_D_v_from_DA_H = ((ADD)**2.d0*z/Hz)**(1.d0/3.d0)

    end function BAO_D_v_from_DA_H

    real(dl) function CAMBdata_BAO_D_v(this,z)
    class(CAMBdata) :: this
    real(dl), intent(IN) :: z

    CAMBdata_BAO_D_v = BAO_D_v_from_DA_H(z,this%AngularDiameterDistance(z), this%Hofz(z))

    end function CAMBdata_BAO_D_v

    function CAMBdata_get_zstar(this)
    class(CAMBdata) :: this
    real(dl) CAMBdata_get_zstar
    real(dl) z_scale

    call this%CP%Recomb%Init(this)
    z_scale =  COBE_CMBTemp/this%CP%TCMB
    CAMBdata_get_zstar=this%binary_search(noreion_optdepth, 1.d0, 700.d0*z_scale, &
        2000.d0*z_scale, 1d-3*z_scale,100.d0*z_scale,5000.d0*z_scale)

    end function CAMBdata_get_zstar

    function CAMBdata_CosmomcTheta(this)
    class(CAMBdata) :: this
    real(dl) zstar, astar, atol, rs, DA
    real(dl) CAMBdata_CosmomcTheta
    real(dl) ombh2, omdmh2

    ombh2 = this%CP%ombh2
    omdmh2 = (this%CP%omch2+this%CP%omnuh2)

    !!From Hu & Sugiyama
    zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+ &
        (0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * &
        (omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))

    astar = 1/(1+zstar)
    atol = 1e-6
    rs = Integrate_Romberg(this,dsound_da_approx,1d-8,astar,atol)
    DA = this%AngularDiameterDistance(zstar)/astar
    CAMBdata_CosmomcTheta = rs/DA

    end function CAMBdata_CosmomcTheta


    subroutine CAMBdata_GetBackgroundDensities(this, n, a_arr, densities)
    ! return array of 8*pi*G*rho*a**4 for each species
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a_arr(n)
    real(dl) :: grhov_t, rhonu, grhonu, a
    real(dl), intent(out) :: densities(8,n)
    integer nu_i,i

    do i=1, n
        a = a_arr(i)
        call this%CP%DarkEnergy%BackgroundDensityAndPressure(this%grhov, a, grhov_t)
        grhonu = 0

        if (this%CP%Num_Nu_massive /= 0) then
            !Get massive neutrino density relative to massless
            do nu_i = 1, this%CP%nu_mass_eigenstates
                call ThermalNuBackground%rho(a * this%nu_masses(nu_i), rhonu)
                grhonu = grhonu + rhonu * this%grhormass(nu_i)
            end do
        end if

        densities(2,i) = this%grhok * a**2
        densities(3,i) = this%grhoc * a
        densities(4,i) = this%grhob * a
        densities(5,i) = this%grhog
        densities(6,i) = this%grhornomass
        densities(7,i) = grhonu
        densities(8,i) = grhov_t*a**2
        densities(1,i) = sum(densities(2:8,i))
    end do

    end subroutine CAMBdata_GetBackgroundDensities

    integer function CAMBdata_get_lmax_lensed(this)
    class(CAMBdata) :: this
    CAMBdata_get_lmax_lensed = this%CLdata%lmax_lensed
    end function CAMBdata_get_lmax_lensed

    !JD 08/13 New function for nonlinear lensing of CMB + MPK compatibility
    !Build master redshift array from array of desired Nonlinear lensing (NLL)
    !redshifts and an array of desired Power spectrum (PK) redshifts.
    !At the same time fill arrays for NLL and PK that indicate indices
    !of their desired redshifts in the master redshift array.
    !Finally define number of redshifts in master array. This is usually given by:
    !P%num_redshifts = P%PK_num_redshifts + NLL_num_redshifts - 1.  The -1 comes
    !from the fact that z=0 is in both arrays (when non-linear is on)
    subroutine GetComputedPKRedshifts(this, Params,eta_k_max)
    use MpiUtils, only : MpiStop
    class(CAMBdata) :: this
    Type(CAMBParams) :: Params
    integer i, iPK, iNLL
    real(dl), parameter :: tol = 1.d-5
    real(dl) maxRedshift, NL_Boost
    integer   ::  NLL_num_redshifts
    real(dl), allocatable    ::  NLL_redshifts(:), redshifts(:)
    !Sources, but unused currently
    real(dl), intent(in), optional :: eta_k_max

    NLL_num_redshifts = 0
    this%needs_good_pk_sampling = .false.
    associate(P => Params%Transfer)
        if ((Params%NonLinear==NonLinear_lens .or. Params%NonLinear==NonLinear_both) .and. &
            (Params%DoLensing .or. this%num_redshiftwindows > 0)) then
            ! Want non-linear lensing or other sources
            this%needs_good_pk_sampling = .false.
            NL_Boost = Params%Accuracy%AccuracyBoost*Params%Accuracy%NonlinSourceBoost
            if (Params%Do21cm) then
                !Sources
                if (maxval(this%Redshift_w(1:this%num_redshiftwindows)%Redshift) &
                    /= minval(this%Redshift_w(1:this%num_redshiftwindows)%Redshift))  &
                    stop 'Non-linear 21cm currently only for narrow window at one redshift'
                if (.not. present(eta_k_max)) stop 'bad call to GetComputedPKRedshifts'
                P%kmax = eta_k_max/10000.
                NLL_num_redshifts =  1
                allocate(NLL_redshifts(NLL_num_redshifts+1))
                NLL_redshifts(1) = this%Redshift_w(1)%Redshift
            else
                P%kmax = max(P%kmax,5*NL_Boost)
                maxRedshift = 10
                NLL_num_redshifts =  nint(10*5*NL_Boost)
                if (NL_Boost>=2.5) then
                    !only notionally more accuracy, more stable for RS
                    maxRedshift =15
                end if
                allocate(NLL_redshifts(NLL_num_redshifts+1)) !+1 to stop access issues below
                do i=1,NLL_num_redshifts
                    NLL_redshifts(i) = real(NLL_num_redshifts-i)/(NLL_num_redshifts/maxRedshift)
                end do
            end if
        end if
        if (allocated(this%transfer_redshifts)) deallocate(this%transfer_redshifts)
        if (NLL_num_redshifts==0) then
            this%num_transfer_redshifts=P%PK_num_redshifts
            allocate(this%transfer_redshifts(this%num_transfer_redshifts))
            this%transfer_redshifts = P%PK_redshifts(:this%num_transfer_redshifts)
            this%PK_redshifts_index(:this%num_transfer_redshifts) = (/ (i, i=1, this%num_transfer_redshifts ) /)
        else
            i=0
            iPK=1
            iNLL=1
            allocate(redshifts(NLL_num_redshifts+P%PK_num_redshifts))
            do while (iPk<=P%PK_num_redshifts .or. iNLL<=NLL_num_redshifts)
                !JD write the next line like this to account for roundoff issues with ==. Preference given to PK_Redshift
                i=i+1
                if(iNLL>NLL_num_redshifts .or. P%PK_redshifts(iPK)>NLL_redshifts(iNLL)+tol) then
                    redshifts(i)=P%PK_redshifts(iPK)
                    this%PK_redshifts_index(iPK)=i
                    iPK=iPK+1
                else if(iPK>P%PK_num_redshifts .or. NLL_redshifts(iNLL)>P%PK_redshifts(iPK)+tol) then
                    redshifts(i)=NLL_redshifts(iNLL)
                    iNLL=iNLL+1
                else
                    redshifts(i)=P%PK_redshifts(iPK)
                    this%PK_redshifts_index(iPK)=i
                    iPK=iPK+1
                    iNLL=iNLL+1
                end if
            end do
            this%num_transfer_redshifts=i
            allocate(this%transfer_redshifts(this%num_transfer_redshifts))
            this%transfer_redshifts = redshifts(:this%num_transfer_redshifts)
        end if
    end associate

    end subroutine GetComputedPKRedshifts

    subroutine CAMBdata_DarkEnergyStressEnergy(this, a, grhov_t, w, n)
    class(CAMBdata) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n)
    real(dl), intent(out) :: grhov_t(n), w(n)
    integer i

    do i=1, n
        call this%CP%DarkEnergy%BackgroundDensityAndPressure(1._dl, a(i), grhov_t(i), w(i))
    end do
    grhov_t = grhov_t/a**2

    end subroutine CAMBdata_DarkEnergyStressEnergy

    function rofChi(this,Chi) !sinh(chi) for open, sin(chi) for closed.
    class(CAMBdata) :: this
    real(dl) Chi,rofChi

    if (this%flat) then
        rofChi=chi
    else if (this%closed) then
        rofChi=sin(chi)
    else
        rofChi=sinh(chi)
    endif
    end function rofChi


    function cosfunc (this,Chi)
    class(CAMBdata) :: this
    real(dl) Chi,cosfunc

    if (this%flat) then
        cosfunc = 1._dl
    else if (this%closed) then
        cosfunc= cos(chi)
    else
        cosfunc=cosh(chi)
    endif
    end function cosfunc

    function tanfunc(this,Chi)
    class(CAMBdata) :: this
    real(dl) Chi,tanfunc
    if (this%flat) then
        tanfunc=Chi
    else if (this%closed) then
        tanfunc=tan(Chi)
    else
        tanfunc=tanh(Chi)
    end if

    end function tanfunc

    function invsinfunc(this,x)
    class(CAMBdata) :: this
    real(dl) invsinfunc,x

    if (this%flat) then
        invsinfunc = x
    else if (this%closed) then
        invsinfunc=asin(x)
    else
        invsinfunc=log((x+sqrt(1._dl+x**2)))
    endif
    end function invsinfunc

    function dsound_da_exact(this,a)
    class(CAMBdata) :: this
    real(dl) dsound_da_exact,a,R,cs

    R = 3*this%grhob*a / (4*this%grhog)
    cs=1.0d0/sqrt(3*(1+R))
    dsound_da_exact=dtauda(this,a)*cs

    end function dsound_da_exact

    function dsound_da_approx(this,a)
    !approximate form used e.g. by CosmoMC for theta
    class(CAMBdata) :: this
    real(dl) dsound_da_approx,a,R,cs

    R=3.0d4*a*this%CP%ombh2
    !          R = 3*grhob*a / (4*grhog) //above is mostly within 0.2% and used for previous consistency
    cs=1.0d0/sqrt(3*(1+R))
    dsound_da_approx=dtauda(this,a)*cs

    end function dsound_da_approx

    function dtda(this,a)
    class(CAMBdata) :: this
    real(dl) dtda,a

    dtda= dtauda(this,a)*a
    end function

    function ddamping_da(this, State, a)
    class(CAMBdata) :: this
    class(Trecfast) :: State
    real(dl) :: ddamping_da
    real(dl), intent(in) :: a
    real(dl) :: R
    R=this%ThermoData%r_drag0*a
    !ignoring reionisation, not relevant for distance measures
    ddamping_da = (R**2 + 16*(1+R)/15)/(1+R)**2*dtauda(this,a)*a**2/(State%total_scattering_eff(this, a)*this%akthom) !BB24

    end function ddamping_da

    function noreion_doptdepth_dz(this,z)
    class(CAMBdata) :: this
    real(dl) :: noreion_doptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a

    a = 1._dl/(1._dl+z)

    !ignoring reionisation, not relevant for distance measures
    noreion_doptdepth_dz = this%CP%Recomb%x_e(a)*this%akthom*dtauda(this,a)

    end function noreion_doptdepth_dz

    function noreion_optdepth(this,z)
    class(CAMBdata) :: this
    real(dl) noreion_optdepth
    real(dl),intent(in) :: z

    noreion_optdepth = Integrate_Romberg(this, noreion_doptdepth_dz, 0.d0, z, 1d-5, 20, 100)

    end function noreion_optdepth

    function ddragoptdepth_dz(this,z)
    class(CAMBdata) :: this
    real(dl) :: ddragoptdepth_dz
    real(dl), intent(in) :: z
    real(dl) :: a

    a = 1._dl/(1._dl+z)
    ddragoptdepth_dz = noreion_doptdepth_dz(this,z)/this%ThermoData%r_drag0/a

    end function ddragoptdepth_dz

    function dragoptdepth(this,z)
    class(CAMBdata) :: this
    real(dl) dragoptdepth
    real(dl),intent(in) :: z

    dragoptdepth =  Integrate_Romberg(this,ddragoptdepth_dz, 0.d0, z, 1d-5, 20, 100)

    end function dragoptdepth

    function reion_doptdepth_dz(this,z)
    class(CAMBdata) :: this
    real(dl) :: reion_doptdepth_dz
    real(dl), intent(in) :: z

    reion_doptdepth_dz = this%CP%Reion%x_e(z)*this%akthom*dtauda(this,1._dl/(1._dl+z))

    end function reion_doptdepth_dz

    function grho_no_de(this, a) result(grhoa2)
    !  Return 8*pi*G*rho_no_de*a**4 where rho_no_de includes everything except dark energy.
    class(CAMBdata) :: this
    real(dl), intent(in) :: a
    real(dl) grhoa2, rhonu
    integer nu_i

    grhoa2 = this%grhok * a**2 + (this%grhoc + this%grhob) * a + this%grhog + this%grhornomass

    if (this%CP%Num_Nu_massive /= 0) then
        !Get massive neutrino density relative to massless
        do nu_i = 1, this%CP%nu_mass_eigenstates
            call ThermalNuBack%rho(a * this%nu_masses(nu_i), rhonu)
            grhoa2 = grhoa2 + rhonu * this%grhormass(nu_i)
        end do
    end if

    end function grho_no_de

    function GetReionizationOptDepth(this)
    class(CAMBdata) :: this
    real(dl) GetReionizationOptDepth
    integer n
    real(dl) zstart, zend

    call this%CP%Reion%get_timesteps(n, zstart, zend)
    GetReionizationOptDepth = Integrate_Romberg(this, reion_doptdepth_dz,0.d0,zstart,&
        1d-5/this%CP%Accuracy%AccuracyBoost)

    end function GetReionizationOptDepth

    real(dl) function binary_search(this,func, goal, x1, x2, tol, widex1, widex2)
    !This is about twice as inefficient as Brent
    class(CAMBdata) :: this
    procedure(state_function) :: func
    real(dl), intent(in) :: goal,x1,x2,tol
    real(dl), intent(in), optional :: widex1, widex2 !Wider range in case of failure
    real(dl) try_t, try_b, avg, D_try, last_bot, last_top, diff
    integer count
    logical wide

    try_b = x1
    try_t = x2
    diff = tol*2
    count = 0
    wide = .false.
    do while (diff > tol)
        if (count>100) then
            if (.not. wide .and. present(widex1)) then
                count=0
                wide=.true.
                try_b=widex1
                try_t=widex2
            else
                call GlobalError(FormatString('binary_search (e.g for optical depth) did not converge: ' //&
                    'Base range %f-%f.',x1,x2),error_reionization)
                binary_search = 0
                return
            end if
        end if
        avg = (try_b+try_t)/2
        D_try = func(this,avg)
        count = count+1
        if (D_try < goal) then
            try_b = avg
            last_bot = D_try
        else
            try_t = avg
            last_top = D_try
        end if
        diff = abs(D_try - goal)
    end do
    if (try_b==x1) last_bot = func(this,x1)
    if (try_t==x2) last_top = func(this,x2)
    binary_search =  (try_t*(goal-last_bot) + try_b*(last_top-goal))/(last_top-last_bot)

    end function binary_search

    function WindowKmaxForL(W,CP, ell) result(res)
    class(CAMBparams), intent(in) :: CP
    Type(TRedWin), intent(in) :: W
    real(dl) res
    integer, intent(in)::  ell

    if (W%kind == window_lensing) then
        res = CP%Accuracy%AccuracyBoost*18*ell/W%chi0
    else
        !On large scales power can be aliased from smaller, so make sure k goes up until at least the turnover
        !in the matter power spectrum
        res = CP%Accuracy%AccuracyBoost*max(0.05_dl,2.5*ell/W%chimin)
    end if

    res = res* CP%Accuracy%KmaxBoost
    end function WindowKmaxForL


    function lSamples_indexOf(lSet,l)
    class(lSamples) :: lSet
    integer, intent(in) :: l
    integer lSamples_indexOf, i

    do i=2,lSet%nl
        if (l < lSet%l(i)) then
            lSamples_indexOf = i-1
            return
        end if
    end do
    lSamples_indexOf = lSet%nl

    end function  lSamples_indexOf

    subroutine lSamples_init(this, State, lmin, max_l)
    ! This subroutines initializes lSet%l arrays. Other values will be interpolated.
    class(lSamples) :: this
    class(CAMBdata), target :: State
    integer, intent(IN) :: lmin,max_l
    integer lind, lvar, step, top, bot, lmin_log
    integer, allocatable :: ls(:)
    real(dl) AScale

    allocate(ls(max_l))
    if (allocated(this%l)) deallocate(this%l)
    this%lmin = lmin
    this%use_spline_template = State%CP%use_cl_spline_template
    lmin_log = State%CP%min_l_logl_sampling
    associate(Accuracy => State%CP%Accuracy)
        Ascale=State%scale/Accuracy%lSampleBoost

        if (Accuracy%lSampleBoost >=50) then
            !just do all of them
            lind=0
            do lvar=lmin, max_l
                lind=lind+1
                ls(lind)=lvar
            end do
            this%nl=lind
            allocate(this%l(lind))
            this%l = ls(1:lind)
            return
        end if

        lind=0
        do lvar=lmin, 10
            lind=lind+1
            ls(lind)=lvar
        end do

        if (Accuracy%AccurateReionization) then
            do lvar=11, 14
                lind=lind+1
                ls(lind)=lvar
            end do
            if (Accuracy%lSampleBoost > 1) then
                do lvar=15, 37, 1
                    lind=lind+1
                    ls(lind)=lvar
                end do
            else
                do lvar=15, 37, 2
                    lind=lind+1
                    ls(lind)=lvar
                end do
            end if

            step = max(nint(5*Ascale),2)
            bot=40
            top=bot + step*10
        else
            if (Accuracy%lSampleBoost >1) then
                do lvar=11, 15
                    lind=lind+1
                    ls(lind)=lvar
                end do
            else
                lind=lind+1
                ls(lind)=12
                lind=lind+1
                ls(lind)=15
            end if
            step = max(nint(10*Ascale),3)
            bot=15+max(step/2,2)
            top=bot + step*7
        end if

        do lvar=bot, top, step
            lind=lind+1
            ls(lind)=lvar
        end do

        if (State%CP%Log_lvalues) then
            !Useful for generating smooth things like 21cm to high l
            step=max(nint(20*Ascale),4)
            do
                lvar = lvar + step
                if (lvar > max_l) exit
                lind=lind+1
                ls(lind)=lvar
                step = nint(step*1.2) !log spacing
            end do
        else
            step=max(nint(20*Ascale),4)
            bot=ls(lind)+step
            top=bot+step*2

            do lvar = bot,top,step
                lind=lind+1
                ls(lind)=lvar
            end do

            if (ls(lind)>=max_l) then
                do lvar=lind,1,-1
                    if (ls(lvar)<=max_l) exit
                end do
                lind=lvar
                if (ls(lind)<max_l) then
                    lind=lind+1
                    ls(lind)=max_l
                end if
            else
                step=max(nint(25*Ascale),4)
                !Get EE right around l=200 by putting extra point at 175
                bot=ls(lind)+step
                top=bot+step

                do lvar = bot,top,step
                    lind=lind+1
                    ls(lind)=lvar
                end do

                if (ls(lind)>=max_l) then
                    do lvar=lind,1,-1
                        if (ls(lvar)<=max_l) exit
                    end do
                    lind=lvar
                    if (ls(lind)<max_l) then
                        lind=lind+1
                        ls(lind)=max_l
                    end if
                else
                    if (.not. this%use_spline_template) then
                        step=max(nint(42*Ascale),7)
                    else
                        step=max(nint(50*Ascale),7)
                    end if
                    bot=ls(lind)+step
                    top=min(lmin_log,max_l)

                    do lvar = bot,top,step
                        lind=lind+1
                        ls(lind)=lvar
                    end do

                    if (max_l > lmin_log) then
                        !Should be pretty smooth or tiny out here
                        step=max(nint(400*Ascale),50)
                        lvar = ls(lind)
                        do
                            lvar = lvar + step
                            if (lvar > max_l) exit
                            lind=lind+1
                            ls(lind)=lvar
                            step = nint(step*1.5) !log spacing
                        end do
                        if (ls(lind) < max_l - 100) then
                            !Try to keep lensed spectra up to specified lmax
                            lind=lind+1
                            ls(lind)=max_l - lensed_convolution_margin
                        else if (ls(lind) - ls(lind-1) > lensed_convolution_margin) then
                            ls(lind)=max_l - lensed_convolution_margin
                        end if
                    end if
                end if !log_lvalues
                if (ls(lind) /=max_l) then
                    lind=lind+1
                    ls(lind)=max_l
                end if
                if (.not. State%flat .and. max_l<=lmin_log) ls(lind-1)=int(max_l+ls(lind-2))/2
                !Not in flat case so interpolation table is the same when using lower l_max
            end if
        end if
    end associate
    this%nl=lind
    allocate(this%l, source=ls(1:lind))

    end subroutine lSamples_init

    subroutine InterpolateClArr(lSet,iCl, all_Cl, max_index)
    class(lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(1:*)
    real(dl), intent(out):: all_Cl(lSet%lmin:*)
    integer, intent(in), optional :: max_index
    integer il,llo,lhi, xi
    real(dl) ddCl(lSet%nl)
    real(dl) xl(lSet%nl)
    real(dl) a0,b0,ho
    integer max_ind

    max_ind = PresentDefault(lSet%nl, max_index)

    if (max_ind > lSet%nl) call MpiStop('Wrong max_ind in InterpolateClArr')

    xl = real(lSet%l(1:lSet%nl),dl)
    call spline_def(xl,iCL,max_ind,ddCl)

    llo=1
    do il=lSet%lmin,lSet%l(max_ind)
        xi=il
        if ((xi > lSet%l(llo+1)).and.(llo < max_ind)) then
            llo=llo+1
        end if
        lhi=llo+1
        ho=lSet%l(lhi)-lSet%l(llo)
        a0=(lSet%l(lhi)-xi)/ho
        b0=(xi-lSet%l(llo))/ho

        all_Cl(il) = a0*iCl(llo)+ b0*iCl(lhi)+((a0**3-a0)* ddCl(llo) &
            +(b0**3-b0)*ddCl(lhi))*ho**2/6
    end do

    end subroutine InterpolateClArr

    subroutine InterpolateClArrTemplated(lSet,iCl, all_Cl, max_ind, template_index)
    class(lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(*)
    real(dl), intent(out):: all_Cl(lSet%lmin:*)
    integer, intent(in) :: max_ind
    integer, intent(in), optional :: template_index
    integer maxdelta, il
    real(dl) DeltaCL(lSet%nl)
    real(dl), allocatable :: tmpall(:)

    if (max_ind > lSet%nl) call MpiStop('Wrong max_ind in InterpolateClArrTemplated')
    if (lSet%use_spline_template .and. present(template_index)) then
        if (template_index<=3) then
            !interpolate only the difference between the C_l and an accurately interpolated template.
            !Using unlensed for template, seems to be good enough
            maxdelta=max_ind
            do while (lSet%l(maxdelta) > lmax_extrap_highl)
                maxdelta=maxdelta-1
            end do
            DeltaCL(1:maxdelta)=iCL(1:maxdelta)- highL_CL_template(lSet%l(1:maxdelta), template_index)

            call lSet%InterpolateClArr(DeltaCl, all_Cl, maxdelta)

            do il=lSet%lmin,lSet%l(maxdelta)
                all_Cl(il) = all_Cl(il) +  highL_CL_template(il,template_index)
            end do

            if (maxdelta < max_ind) then
                !directly interpolate high L where no t  emplate (doesn't effect lensing spectrum much anyway)
                allocate(tmpall(lSet%lmin:lSet%l(max_ind)))
                call InterpolateClArr(lSet,iCl, tmpall, max_ind)
                !overlap to reduce interpolation artefacts
                all_cl(lSet%l(maxdelta-2):lSet%l(max_ind) ) = tmpall(lSet%l(maxdelta-2):lSet%l(max_ind))
                deallocate(tmpall)
            end if
            return
        end if
    end if

    call InterpolateClArr(lSet,iCl, all_Cl, max_ind)

    end subroutine InterpolateClArrTemplated

    subroutine Thermo_values_array(this,tau,a,cs2b,opacity, dopacity) ! BB24
        !Compute unperturbed sound speed squared,
        !and ionization fraction by interpolating pre-computed tables.
        !If requested also get time derivative of opacity
        class(TThermoData) :: this
        real(dl), intent(in) :: tau
        real(dl), intent(out) :: a, cs2b
        real(dl), intent(out) :: opacity(nscatter)
        real(dl), intent(out), optional :: dopacity(nscatter)
    
        integer i
        real(dl) d
    
        d=log(tau/this%tauminn)/this%dlntau+1._dl
        i=int(d)
        d=d-i
        if (i < 1) then
            !Linear interpolation if out of bounds (should not occur).
            write(*,*) 'tau, taumin = ', tau, this%tauminn
            call MpiStop('thermo out of bounds')
        else if (i >= this%nthermo) then
            cs2b=this%cs2(this%nthermo)
            opacity=this%dotmu(this%nthermo,:)
            a=1
            if (present(dopacity)) then
                dopacity = this%ddotmu(this%nthermo,:)/(tau*this%dlntau)
            end if
        ! andrea
        else
            cs2b=this%cs2(i)+d*(this%dcs2(i)+d*(3*(this%cs2(i+1)-this%cs2(i))  &
                -2*this%dcs2(i)-this%dcs2(i+1)+d*(this%dcs2(i)+this%dcs2(i+1)  &
                +2*(this%cs2(i)-this%cs2(i+1)))))
            opacity=this%dotmu(i,:)+d*(this%ddotmu(i,:)+d*(3*(this%dotmu(i+1,:)-this%dotmu(i,:)) &
                -2*this%ddotmu(i,:)-this%ddotmu(i+1,:)+d*(this%ddotmu(i,:)+this%ddotmu(i+1,:) &
                +2*(this%dotmu(i,:)-this%dotmu(i+1,:)))))
            a = (this%ScaleFactor(i)+d*(this%dScaleFactor(i)+d*(3*(this%ScaleFactor(i+1)-this%ScaleFactor(i)) &
                -2*this%dScaleFactor(i)-this%dScaleFactor(i+1)+d*(this%dScaleFactor(i)+this%dScaleFactor(i+1) &
                +2*(this%ScaleFactor(i)-this%ScaleFactor(i+1))))))*tau
            if (present(dopacity)) then
                dopacity=(this%ddotmu(i,:)+d*(this%dddotmu(i,1)+d*(3*(this%ddotmu(i+1,:)  &
                    -this%ddotmu(i,:))-2*this%dddotmu(i,:)-this%dddotmu(i+1,:)+d*(this%dddotmu(i,1) &
                    +this%dddotmu(i+1,:)+2*(this%ddotmu(i,:)-this%ddotmu(i+1,:))))))/(tau*this%dlntau)
        ! 
            end if
        end if
    end subroutine Thermo_values_array
    
    subroutine Thermo_expansion_values(this, tau, a, adot, opacity)
        class(TThermoData) :: this
        real(dl), intent(in) :: tau
        real(dl), intent(out) :: a, adot, opacity(nscatter)
        integer i
        real(dl) d
    
        d=log(tau/this%tauminn)/this%dlntau+1._dl
        i=int(d)
        d=d-i
        if (i < 1) then
            call MpiStop('thermo out of bounds')
        else if (i >= this%nthermo) then
            opacity=this%dotmu(this%nthermo,1)
            a=1
            adot=this%adot(this%nthermo)
        else
            a = (this%ScaleFactor(i)+d*(this%dScaleFactor(i)+d*(3*(this%ScaleFactor(i+1)-this%ScaleFactor(i)) &
                -2*this%dScaleFactor(i)-this%dScaleFactor(i+1)+d*(this%dScaleFactor(i)+this%dScaleFactor(i+1) &
                +2*(this%ScaleFactor(i)-this%ScaleFactor(i+1))))))*tau
    
            adot = (this%adot(i)+d*(this%dadot(i)+d*(3*(this%adot(i+1)-this%adot(i)) &
                -2*this%dadot(i)-this%dadot(i+1)+d*(this%dadot(i)+this%dadot(i+1) &
                +2*(this%adot(i)-this%adot(i+1))))))
        ! adnrea
            opacity=this%dotmu(i,1)+d*(this%ddotmu(i,1)+d*(3*(this%dotmu(i+1,1)-this%dotmu(i,1)) &
                -2*this%ddotmu(i,1)-this%ddotmu(i+1,1)+d*(this%ddotmu(i,1)+this%ddotmu(i+1,1) &
                +2*(this%dotmu(i,1)-this%dotmu(i+1,1)))))
        ! 
        end if
    
    end subroutine Thermo_expansion_values
    
    subroutine Thermo_expansion_values_array(this, tau, a, adot, opacity) 
        class(TThermoData) :: this
        real(dl), intent(in) :: tau
        real(dl), intent(out) :: a, adot, opacity(nscatter)
        integer i
        real(dl) d
        d=log(tau/this%tauminn)/this%dlntau+1._dl
        i=int(d)
        d=d-i
        if (i < 1) then
            call MpiStop('thermo out of bounds')
        else if (i >= this%nthermo) then
            opacity=this%dotmu(this%nthermo,:)
            a=1
            adot=this%adot(this%nthermo)
        else
            a = (this%ScaleFactor(i)+d*(this%dScaleFactor(i)+d*(3*(this%ScaleFactor(i+1)-this%ScaleFactor(i)) &
                -2*this%dScaleFactor(i)-this%dScaleFactor(i+1)+d*(this%dScaleFactor(i)+this%dScaleFactor(i+1) &
                +2*(this%ScaleFactor(i)-this%ScaleFactor(i+1))))))*tau
            adot = (this%adot(i)+d*(this%dadot(i)+d*(3*(this%adot(i+1)-this%adot(i)) &
                -2*this%dadot(i)-this%dadot(i+1)+d*(this%dadot(i)+this%dadot(i+1) &
                +2*(this%adot(i)-this%adot(i+1))))))
            ! adnrea
            opacity=this%dotmu(i,:)+d*(this%ddotmu(i,:)+d*(3*(this%dotmu(i+1,:)-this%dotmu(i,:)) &
                -2*this%ddotmu(i,:)-this%ddotmu(i+1,:)+d*(this%ddotmu(i,:)+this%ddotmu(i+1,:) &
                +2*(this%dotmu(i,:)-this%dotmu(i+1,:)))))
            ! 
        end if        
    end subroutine Thermo_expansion_values_array

    function Thermo_OpacityToTime(this,opacity)
    class(TThermoData) :: this
    real(dl), intent(in) :: opacity
    integer j
    real(dl) Thermo_OpacityToTime
    !Do this the bad slow way for now..
    !The answer is approximate
    j =1
    do while(this%dotmu(j,1)> opacity)
        j=j+1
    end do

    Thermo_OpacityToTime = exp((j-1)*this%dlntau)*this%tauminn

    end function Thermo_OpacityToTime

    subroutine Thermo_Init(this, etat, State, taumin)
    !  Compute and save unperturbed baryon temperature and ionization fraction
    !  as a function of time.  With nthermo=10000, xe(tau) has a relative
    ! accuracy (numerical integration precision) better than 1.e-5.
    use constants
    use StringUtils
    class(TThermoData) :: this
    class(TRecfast) :: etat
    class(CAMBdata), target :: State
    !class(TRecfast), allocatable :: Statee
    real(dl), intent(in) :: taumin
    integer nthermo
    real(dl) tau01,a0,barssc,dtau
    real(dl) tau,a,a2
    real(dl) adot,fe,thomc0
    real(dl) dtbdla,vfi,cf1,maxvis, z_scale, vis
    ! real(dl), dimension(:,:), allocatable :: vis
    integer ncount,i,j1,iv,ns
    real(dl), allocatable :: spline_data(:,:) !BB24
    real(dl) last_dotmu, om
    real(dl) a_verydom
    real(dl) awin_lens1p,awin_lens2p,dwing_lens, rs, DA
    real(dl) a_eq, rs_eq, tau_eq, rstar
    integer noutput, f_i
    Type(CalWins), dimension(:), allocatable, target :: RW
    real(dl) awin_lens1(State%num_redshiftwindows),awin_lens2(State%num_redshiftwindows)
    real(dl) Tspin, Trad, rho_fac, window, tau_eps
    integer transfer_ix(State%CP%Transfer%PK_num_redshifts)
    integer RW_i, j2
    real(dl) Tb21cm, winamp, z, background_boost
    character(len=:), allocatable :: outstr
    real(dl), allocatable ::  taus(:)
    ! andrea
    real(dl), allocatable :: xe_a(:), opts(:)
    real(dl), allocatable :: sdotmu(:,:)
    real(dl), allocatable :: scale_factors(:), times(:), dt(:)
    Type(TCubicSpline) :: dotmuSp
    integer ninverse, nlin
    real(dl) dlna, zstar_min, zstar_max
    real(dl) reion_z_start, reion_z_complete
    Type(CAMBParams), pointer :: CP
    ! andrea
    real(dl) dq, q, dlfdlq
    logical :: plot_scatter = .false.
    real(dl) elec_fac
    ! andrea
    real(dl), parameter :: nu_eff = 3101692._dl !3125349._dl is approx from Yu paper
    ! Allocate memory for Statee
    !allocate(Statee)
    ! andrea

    CP => State%CP
    if (num_cmb_freq<10) then
        phot_freqs(1:6) = [0, 143, 217,353, 545, 857]
!       phot_freqs(1:8) = [220,265,300,320,295,460,555,660]*1.085 !Prism
        do i=1, size(phot_freqs)
            q = phot_freqs(i)/56.8
            !this should not be used, just for code consistency
            if (i==1) then
                dq= (phot_freqs(i+1)/56.8-q)*2
            elseif (i==size(phot_freqs)) then
                dq= (q-phot_freqs(i-1)/56.8)*2
            else
                dq = (phot_freqs(i+1)-phot_freqs(i-1))/2/56.8
            end if
            ! andrea
            if (q==0._dl) then
                phot_int_kernel(i)=0
            else
            ! andrea
                dlfdlq=-q/(1._dl-exp(-q))
                phot_int_kernel(i)=dq*q**3/(exp(q)-1._dl) * (-0.25_dl*dlfdlq)
            end if
        end do
    else
        dq = 18/real(num_cmb_freq)
        do i=1,num_cmb_freq
            q=(i-0.5d0)*dq
            phot_freqs(i) = 56.8*q !phot_freqs in GHz
            dlfdlq=-q/(1._dl-exp(-q))
            phot_int_kernel(i)=dq*q**3/(exp(q)-1._dl) * (-0.25_dl*dlfdlq) !now evolve 4F_l/dlfdlq(i)
        end do
        phot_int_kernel=phot_int_kernel/sum(phot_int_kernel) !  (Pi**4/15)
    end if
    print *, 'Doing frequencies: ', phot_freqs
    freq_factors(:,1) = (phot_freqs/ nu_eff)**4
    freq_factors(:,2) = (phot_freqs/ nu_eff)**6 * 638._dl/243
    freq_factors(:,3) = (phot_freqs/ nu_eff)**8 * 1299667._dl/236196 !!Fix 1626820991._dl/136048896._dl
    !These are int q^n q^3*F *(-1/4)*(d log F/dlog q) / int q^3 F
    av_freq_factors(1) = (356.88/ nu_eff)**4
    av_freq_factors(2) = (409.22/ nu_eff)**6 * 638._dl/243
    av_freq_factors(3) = (459.8/ nu_eff)**8  * 1299667._dl/236196 !!Fix 1626820991._dl/136048896._dl 
    if (.not. rayleigh_pows(1)) freq_factors(:,1)=0
    if (.not. rayleigh_pows(2)) freq_factors(:,2)=0
    if (.not. rayleigh_pows(3)) freq_factors(:,3)=0
    ! andrea

    !Allocate memory outside parallel region to keep ifort happy
    background_boost = CP%Accuracy%BackgroundTimeStepBoost*CP%Accuracy%AccuracyBoost
    if (background_boost > 20) then
        write(*,*) 'Warning: very small time steps can give less accurate spline derivatives'
        write(*,*) 'e.g. around reionization if not matched very smoothly'
    end if
    !Higher q starts earlier; scale by log(taumin) so actual step size is not made worse by increasing k_max
    nthermo = nint(thermal_history_def_timesteps*log(1.4e4/taumin)/log(1.4e4/2e-4)*background_boost)
    this%tauminn=0.95d0*taumin
    this%dlntau=log(State%tau0/this%tauminn)/(nthermo-1)

    do RW_i = 1, State%num_redshiftwindows
        !Make sure steps small enough for any features in source window functions
        associate (Win => State%Redshift_w(RW_i))
            if ((Win%kind /= window_21cm .or. .not. CP%transfer_21cm_cl) .and. &
                Win%sigma_tau/5/background_boost < Win%tau*(exp(this%dlntau)-1)) then
                this%dlntau = log(Win%sigma_tau/5/background_boost/Win%tau+1)
                nthermo = nint(log(State%tau0/this%tauminn)/this%dlntau) + 1
                this%dlntau=log(State%tau0/this%tauminn)/(nthermo-1)
            end if
        end associate
    end do
    this%nthermo = nthermo
    allocate(spline_data(nthermo, nscatter), sdotmu(nthermo, nscatter)) ! BB24

    if (allocated(this%tb) .and. this%nthermo/=size(this%tb)) then
        deallocate(this%scaleFactor, this%cs2, this%dcs2, this%ddotmu)
        deallocate(this%dscaleFactor, this%adot, this%dadot)
        deallocate(this%tb, this%xe, this%emmu, this%dotmu)
        deallocate(this%demmu, this%dddotmu, this%ddddotmu)
        if (dowinlens .and. allocated(this%winlens)) deallocate(this%winlens, this%dwinlens)
    endif
    if (.not. allocated(this%tb)) then
        allocate(this%scaleFactor(nthermo), this%cs2(nthermo), this%dcs2(nthermo), this%ddotmu(nthermo,nscatter)) !BB24
        allocate(this%dscaleFactor(nthermo), this%adot(nthermo), this%dadot(nthermo))
        allocate(this%tb(nthermo), this%xe(nthermo), this%emmu(nthermo,nscatter),this%dotmu(nthermo,nscatter))
        allocate(this%demmu(nthermo, nscatter), this%dddotmu(nthermo, nscatter), this%ddddotmu(nthermo, nscatter))
        if (dowinlens) allocate(this%winlens(nthermo), this%dwinlens(nthermo))
    end if

    if (State%num_redshiftwindows >0) then
        allocate(this%redshift_time(nthermo),this%dredshift_time(nthermo))
        allocate(this%arhos_fac(nthermo), this%darhos_fac(nthermo), this%ddarhos_fac(nthermo))
        allocate(RW(State%num_redshiftwindows))
    end if


    do RW_i = 1, State%num_redshiftwindows
        associate (RedWin => State%Redshift_w(RW_i))
            RedWin%tau_start = 0
            RedWin%tau_end = State%tau0
            if (RedWin%kind == window_lensing .or.  RedWin%kind == window_counts .and. CP%SourceTerms%counts_lensing) then
                allocate(RW(RW_i)%awin_lens(nthermo))
                allocate(RW(RW_i)%dawin_lens(nthermo))
            end if
        end associate
    end do
    om = (State%grhob+State%grhoc)/&
        sqrt(3*(State%grhog+sum(State%grhormass(1:CP%Nu_mass_eigenstates))+State%grhornomass))
    a0=this%tauminn*State%adotrad*(1+om*this%tauminn/4)
    ninverse = nint(background_boost*log(1/a0)/log(1/2d-10)*4000)
    if (.not. CP%DarkEnergy%is_cosmological_constant) ninverse = ninverse*2

    nlin = ninverse/2
    allocate(scale_factors(ninverse+nlin))
    allocate(times(ninverse+nlin))
    allocate(dt(ninverse+nlin))
    allocate(taus(nthermo), xe_a(nthermo))

    !$OMP PARALLEL SECTIONS DEFAULT(SHARED)
    !$OMP SECTION
    !call etat%Init(State,WantTSpin=CP%Do21cm) 
    call CP%Recomb%Init(State,WantTSpin=CP%Do21cm)    !almost all the time spent here
    write(*,*) 'ok'
    if (CP%Evolve_delta_xe) this%recombination_saha_tau  = State%TimeOfZ(CP%Recomb%get_saha_z(), tol=1e-4_dl)
    if (CP%Evolve_baryon_cs .or. CP%Evolve_delta_xe .or. CP%Evolve_delta_Ts .or. CP%Do21cm) &
        this%recombination_Tgas_tau = State%TimeOfz(1/CP%Recomb%min_a_evolve_Tm-1, tol=1e-4_dl)

    !$OMP SECTION
    !Do other stuff while recombination calculating
    awin_lens1=0
    awin_lens2=0
    transfer_ix =0

    call splini(spline_data,nthermo)

    this%tight_tau = 0
    this%actual_opt_depth = 0
    ncount=0
    this%z_drag=0.d0
    thomc0= Compton_CT * CP%tcmb**4
    this%r_drag0 = 3.d0/4.d0*State%grhob/State%grhog
    last_dotmu = 0

    this%matter_verydom_tau = 0
    a_verydom = CP%Accuracy%AccuracyBoost*5*(State%grhog+State%grhornomass)/(State%grhoc+State%grhob)
    if (CP%Reion%Reionization) then
        call CP%Reion%get_timesteps(State%reion_n_steps, reion_z_start, reion_z_complete)
        State%reion_tau_start = max(0.05_dl, State%TimeOfZ(reion_z_start, 1d-3))
        !Time when a very small reionization fraction (assuming tanh fitting)
        State%reion_tau_complete = min(State%tau0, &
            State%reion_tau_start+ State%DeltaTime(1/(1+reion_z_start),1/(1.d0+reion_z_complete),1d-3))
    else
        State%reion_tau_start = State%tau0
        State%reion_tau_complete = State%tau0
    end  if
    !  Initial conditions: assume radiation-dominated universe.
    !  Assume that any entropy generation occurs before tauminn.
    !  This gives wrong temperature before pair annihilation, but
    !  the error is harmless.

    !Get scale factor as function of time by inverting tau(a)
    dlna = log(0.2_dl/a0)/(ninverse-1)
    do i=2, ninverse-1
        scale_factors(1+i) = a0*exp((i-1)*dlna)
    end do
    scale_factors(1) = a0
    scale_factors(2) = a0*exp(dlna/3)
    da = 0.8_dl/(nlin-2)
    do i=1, nlin-2
        scale_factors(ninverse+i) = 0.2_dl + (i-1)*da
    end do
    scale_factors(ninverse+nlin-1) = 0.9_dl + 0.1_dl*scale_factors(ninverse+nlin-2)
    scale_factors(ninverse+nlin) = 1
    do i=1, ninverse+nlin
        dt(i) = dtauda(State,scale_factors(i))
    end do
    call this%ScaleFactorAtTime%Init(scale_factors, dt)
    call this%ScaleFactorATTime%IntegralArray(times(2), first_index=2)
    times(1) = this%tauminn
    times(2:) = times(2:) + 2*(sqrt(1 + om*scale_factors(2)/ State%adotrad) -1)/om
    times(ninverse+nlin) = State%tau0
    call This%ScaleFactorAtTime%Init(times, scale_factors)
    taus(1) = this%tauminn
    do i=2,nthermo-1
        taus(i) = this%tauminn*exp((i-1)*this%dlntau)
    end do
    taus(nthermo) = State%tau0
    call this%ScaleFactorAtTime%Array(taus(2:), this%scaleFactor(2:))
    this%scaleFactor(1) = a0
    this%scaleFactor(nthermo) = 1
    this%adot(1) = 1/dtauda(State,a0)

    tau01=this%tauminn
    do i=2,nthermo
        !Get recombination-independent parts of background now as function of conformal time tau
        !Now at the log spaced time steps
        tau=taus(i)
        dtau = tau-tau01
        a = this%scaleFactor(i)
        adot = 1/dtauda(State,a)
        this%adot(i) = adot
        if (this%matter_verydom_tau ==0 .and. a > a_verydom) then
            this%matter_verydom_tau = tau
        end if
        z= 1._dl/a-1._dl
        if (State%num_redshiftwindows>0) then
            this%redshift_time(i) = z
            do RW_i = 1, State%num_redshiftwindows
                associate (Win => RW(RW_i), RedWin => State%Redshift_w(RW_i))
                    if (a > 1d-4) then
                        window = RedWin%Window%Window_f_a(a, winamp)

                        if  (RedWin%kind == window_lensing .or.  RedWin%kind == window_counts  &
                            .and. CP%SourceTerms%counts_lensing) then
                            if (State%tau0 - tau > 2) then
                                dwing_lens =  adot * window *dtau
                                awin_lens1(RW_i) = awin_lens1(RW_i) + dwing_lens
                                awin_lens2(RW_i) = awin_lens2(RW_i) + dwing_lens/(State%tau0-tau)
                                Win%awin_lens(i) = awin_lens1(RW_i)/(State%tau0-tau) - awin_lens2(RW_i)
                            else
                                Win%awin_lens(i) = 0
                            end if
                        end if

                        if (RedWin%tau_start ==0 .and. winamp > 1e-8) then
                            RedWin%tau_start = tau01
                        else if (RedWin%tau_start /=0 .and. RedWin%tau_end==State%tau0 .and. winamp < 1e-8) then
                            RedWin%tau_end = min(State%tau0,tau + dtau)
                            if (DebugMsgs) call WriteFormat('Window %d: tau1 = %f, tau2 = %f',&
                                RW_i, RedWin%tau_start, RedWin%tau_end)
                        end if
                    else
                        if (RedWin%kind == window_lensing .or.  RedWin%kind == window_counts &
                            .and. CP%SourceTerms%counts_lensing) then
                            Win%awin_lens(i)=0
                        end if
                    end if
                end associate
            end do
        end if
        if (CP%WantTransfer .and.  CP%do21cm .and. CP%transfer_21cm_cl) then
            do RW_i = 1, CP%Transfer%PK_num_redshifts
                if (z< CP%Transfer%PK_redshifts(RW_i) .and. transfer_ix(RW_i)==0) then
                    transfer_ix(RW_i) = i
                end if
            end do
        end if
        tau01 =tau
    end do
    do RW_i = 1, State%num_redshiftwindows
        associate(Win => RW(RW_i))
            if (State%Redshift_w(RW_i)%kind == window_lensing .or. &
                State%Redshift_w(RW_i)%kind == window_counts .and. CP%SourceTerms%counts_lensing) then
                this%has_lensing_windows = .true.
                State%Redshift_w(RW_i)%has_lensing_window = .true.
                if (FeedbackLevel>0)  write(*,'(I1," Int W              = ",f9.6)') RW_i, awin_lens1(RW_i)
                Win%awin_lens=Win%awin_lens/awin_lens1(RW_i)
            else
                State%Redshift_w(RW_i)%has_lensing_window = .false.
            end if
        end associate
    end do
    !$OMP END PARALLEL SECTIONS

    if (global_error_flag/=0) return

    call CP%Recomb%xe_tm(a0,this%xe(1), this%tb(1))
    barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(1))
    this%cs2(1)=4._dl/3._dl*barssc*this%tb(1)
    ! andrea
    this%dotmu(1,1)=this%xe(1)*State%akthom/a0**2 !! BB to check
    sdotmu(1,:)=0
    this%dotmu(1,2:)=0
    this%scaleFactor(1)=a0
    ! andrea

    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC,16)
    do i=2,nthermo
        call CP%Recomb%xe_tm(this%scaleFactor(i), xe_a(i), this%tb(i))
    end do

    do i=2,nthermo
        tau =taus(i)
        a = this%scaleFactor(i)
        a2=a*a
        adot=this%adot(i)

        if (State%num_redshiftwindows>0) then
            if (a > 1d-4) then
                if (CP%Do21cm ) then
                    Tspin = CP%Recomb%T_s(a)
                    Trad = CP%TCMB/a
                    rho_fac = line21_const*State%NNow/a**3
                    tau_eps = a*line21_const*State%NNow/a**3/(adot/a)/Tspin/1000
                    this%arhos_fac(i) = (1-exp(-tau_eps))/tau_eps*a*rho_fac*(1 - Trad/Tspin)/(adot/a)
                    !         arhos_fac(i) = a*rho_fac*(1 - Trad/Tspin)/(adot/a)
                else
                    rho_fac = State%grhoc/(a2*a)
                    this%arhos_fac(i) = a*rho_fac/(adot/a)
                end if
            end if
        end if

        ! If there is re-ionization, smoothly increase xe to the
        ! requested value.
        if (CP%Reion%Reionization .and. tau > State%reion_tau_start) then
            if(ncount == 0) then
                ncount=i-1
            end if
            this%xe(i) = CP%Reion%x_e(1/a-1, tau, this%xe(ncount))
            if (CP%Accuracy%AccurateReionization .and. CP%WantDerivedParameters) then
                ! and
                !this%dotmu(i,1)=(etat%total_scattering_eff(a) - this%xe(i))*State%akthom/a2
                this%dotmu(i,1) = (etat%total_scattering_eff(State, a) - this%xe(i)) * State%akthom / a2

                ! and
                if (last_dotmu /=0) then
                    ! andrea solo change
                    this%actual_opt_depth = this%actual_opt_depth - 2._dl*(tau-taus(i-1))/(1._dl/this%dotmu(i,1)+1._dl/last_dotmu)
                end if
                last_dotmu = this%dotmu(i,1)
                ! 
            end if
        else
            this%xe(i)=etat%total_scattering_eff(State, a)
        end if

        !  approximate Baryon sound speed squared (over c**2).
        fe=(1._dl-CP%yhe)*this%xe(i)/(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(i))
        dtbdla=-2._dl*this%tb(i)
        if (a*this%tb(i)-CP%tcmb < -1e-8) then
            dtbdla= dtbdla -thomc0*fe/adot*(a*this%tb(i)-CP%tcmb)/a**3
        end if
        barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(i))
        this%cs2(i)=barssc*this%tb(i)*(1-dtbdla/this%tb(i)/3._dl)

        ! Calculation of the visibility function

        ! Andrea
        this%dotmu(i,1)=this%xe(i)*State%akthom/a2
        if (plot_scatter) then
            elec_fac= 0
        else
            elec_fac =1
        end if
        !if (.not. allocated(statee%Calc)) allocate(statee%Calc)
        !Calc => Statee%Calc
        do f_i=1,num_cmb_freq
            this%dotmu(i,1+f_i)=this%dotmu(i,1)*elec_fac + etat%recombination_rayleigh_eff(a)*State%akthom/a2*(min(1._dl,&
            freq_factors(f_i,1)/a2**2 + freq_factors(f_i,2)/a2**3  + freq_factors(f_i,3)/a2**4 ))
        end do
        ! andrea

        if (this%tight_tau==0 .and. 1/(tau*this%dotmu(i,1)) > 0.005) this%tight_tau = tau !0.005
        !
        !Tight coupling switch time when k/opacity is smaller than 1/(tau*opacity)
    end do

    if (CP%Reion%Reionization .and. (this%xe(nthermo) < 0.999d0)) then
        write(*,*)'Warning: xe at redshift zero is < 1'
        write(*,*) 'Check input parameters an Reionization_xe'
        write(*,*) 'function in the Reionization module'
    end if

    !Integrate for optical depth
    ! andrea solo change
    call dotmuSp%Init(taus(nthermo:1:-1), this%dotmu(nthermo:1:-1,1))
    ! 
    allocate(opts(nthermo))
    call dotmuSp%IntegralArray(opts)
    ! and solo change
    sdotmu(:,1) = opts(nthermo:1:-1)
    ! and solo change
    do j1=1,nthermo
        ! andrea change
        if (sdotmu(j1,1)< -69) then
            this%emmu(j1,:)=1.d-30
        else
            this%emmu(j1,:)=exp(sdotmu(j1,:))
            if (CP%Reion%Reionization .and. .not. CP%Accuracy%AccurateReionization .and. &
                this%actual_opt_depth==0 .and. this%xe(j1) < 1e-3) then
                this%actual_opt_depth = -sdotmu(j1,1)
        ! 
            end if
        end if
    end do
    z_scale =  COBE_CMBTemp/CP%TCMB
    zstar_min = 700._dl * z_scale
    zstar_max = 2000._dl * z_scale
    if ((.not. CP%Reion%Reionization .or. CP%Accuracy%AccurateReionization) .and. CP%WantDerivedParameters) then
        do j1=nint(log(100/this%tauminn)/this%dlntau),nthermo
            ! andre change
            if (-sdotmu(j1,1) - this%actual_opt_depth < 1) then
            ! 
                !Bracket z_star
                zstar_min = 1/this%scaleFactor(j1+1)-1
                zstar_max = 1/this%scaleFactor(j1-2)-1
                exit
            end if
        end do
    end if

    if (CP%WantTransfer .and.  CP%do21cm .and. CP%transfer_21cm_cl) then
        if (allocated(State%optical_depths_for21cm)) deallocate(State%optical_depths_for21cm)
        allocate(State%optical_depths_for21cm(CP%Transfer%PK_num_redshifts))
        do RW_i = 1, CP%Transfer%PK_num_redshifts
            if (CP%Transfer%PK_Redshifts(RW_i) < 1e-3) then
                State%optical_depths_for21cm(RW_i) = 0 !zero may not be set correctly in transfer_ix
            else
                ! and
                State%optical_depths_for21cm(RW_i) =  -sdotmu(transfer_ix(RW_i),1)
                ! and
            end if
        end do
    end if

    if (CP%Reion%Reionization .and. CP%Accuracy%AccurateReionization &
        .and. FeedbackLevel > 0 .and. CP%WantDerivedParameters) then
        write(*,'("Reion opt depth      = ",f7.4)') this%actual_opt_depth
    end if

    ! ! andrea
    ! if (plot_scatter) then
    !     call CreateTxtFile('c:\tmp\planck\rayleigh\fixed\visibilities_tot.txt',1)
    !     do j1=1,nthermo !!!!
    !         tau = tauminn*exp((j1-1)*dlntau)
    !         write(1,'(9E15.5)') tau, scaleFactor(j1), this%dotmu(j1,1)*scaleFactor(j1)**2/State%akthom, real(this%emmu(j1,1)*this%dotmu(j1,1)),real(this%emmu(j1,3:7)*this%emmu(j1,1)*this%dotmu(j1,3:7))
    !     end do
    !     close(1)
    !     call CreateTxtFile('c:\tmp\planck\rayleigh\fixed\taudot_tot.txt',1)
    !     do j1=1,nthermo !!!!
    !         tau = tauminn*exp((j1-1)*dlntau)
    !         write(1,'(14E15.5)') tau, scaleFactor(j1), this%dotmu(j1,1)*scaleFactor(j1)**2/State%akthom, real(this%dotmu(j1,1)),real(this%dotmu(j1,3:7)),real(this%emmu(j1,3:7))
    !     end do
    !     close(1)
    !     stop
    ! end if
    ! andrea

    iv=0
    vfi=0._dl
    ! Getting the starting and finishing times for decoupling and time of maximum visibility
    if (ncount == 0) then
        cf1=1._dl
        ns=nthermo
    else
        ! andre
        cf1=exp(-sdotmu(ncount,1))
        ! 
        ns=ncount
    end if
    maxvis = 0
    do j1=1,ns
        ! andrea
        vis = this%emmu(j1,1)*this%dotmu(j1,1)
        ! 
        !
        tau = taus(j1)
        vfi=vfi+vis*cf1*this%dlntau*tau
        if ((iv == 0).and.(vfi > 1.0d-7/CP%Accuracy%AccuracyBoost)) then
            State%taurst=9._dl/10._dl*tau
            iv=1
        elseif (iv == 1) then
            if (vis > maxvis) then
                maxvis=vis
                State%tau_maxvis = tau
            end if
            if (vfi > 0.995) then
                State%taurend=tau
                iv=2
                exit
            end if
        end if
    end do

    if (iv /= 2) then
        call GlobalError('ThemoData Init: failed to find end of recombination',error_reionization)
        return
    end if

    if (dowinlens) then
        vfi=0
        awin_lens1p=0
        awin_lens2p=0
        this%winlens=0
        do j1=1,nthermo-1
            ! andrea change
            vis = this%emmu(j1,1)*this%dotmu(j1,1)
            ! 
            tau = this%tauminn* taus(j1)
            vfi=vfi+vis*cf1*this%dlntau*tau
            if (vfi < 0.995) then
                dwing_lens =  vis*cf1*this%dlntau*tau / 0.995

                awin_lens1p = awin_lens1p + dwing_lens
                awin_lens2p = awin_lens2p + dwing_lens/(State%tau0-tau)
            end if
            this%winlens(j1)= awin_lens1p/(State%tau0-tau) - awin_lens2p
        end do
    end if

    ! Calculating the timesteps during recombination.

    if (CP%WantTensors) then
        State%dtaurec=min(State%dtaurec,State%taurst/160)/CP%Accuracy%AccuracyBoost
    else
        State%dtaurec=min(State%dtaurec,State%taurst/40)/CP%Accuracy%AccuracyBoost
        if (do_bispectrum .and. hard_bispectrum) State%dtaurec = State%dtaurec / 4
    end if

    if (CP%Reion%Reionization) State%taurend=min(State%taurend,State%reion_tau_start)

    if (DebugMsgs) then
        write (*,*) 'taurst, taurend = ', State%taurst, State%taurend
    end if

    !$OMP PARALLEL SECTIONS DEFAULT(SHARED)
    !$OMP SECTION
    ! andrea big change !
    do f_i=1, nscatter
        call splder(this%dotmu(1,f_i),this%ddotmu(1,f_i),nthermo,spline_data)
        call splder(this%ddotmu(1,f_i),this%dddotmu(1,f_i),nthermo,spline_data)
        call splder(this%dddotmu(1,f_i),this%ddddotmu(1,f_i),nthermo,spline_data)
        call splder(this%emmu(1,f_i),this%demmu(1,f_i),nthermo,spline_data)
    end do
    ! 
    if (CP%want_zstar .or. CP%WantDerivedParameters) &
        this%z_star = State%binary_search(noreion_optdepth, 1.d0, zstar_min, zstar_max, &
        & 1d-3/background_boost, 100._dl*z_scale, 4000._dl*z_scale)
    !$OMP SECTION
    call splder(this%cs2,this%dcs2,nthermo,spline_data)
    ! i deleted call splder(this%emmu(1,f_i) here !
    call splder(this%adot,this%dadot,nthermo,spline_data)
    if (dowinlens) call splder(this%winlens,this%dwinlens,nthermo,spline_data)
    if (CP%want_zdrag .or. CP%WantDerivedParameters) &
        this%z_drag = State%binary_search(dragoptdepth, 1.d0, 800*z_scale, &
        & max(zstar_max*1.1_dl,1200._dl*z_scale), 2d-3/background_boost, 100.d0*z_scale, 4000._dl*z_scale)
    !$OMP SECTION
    this%ScaleFactor(:) = this%scaleFactor/taus !a/tau
    this%dScaleFactor(:) = (this%adot - this%ScaleFactor)*this%dlntau !derivative of a/tau
    if (State%num_redshiftwindows >0) then
        call splder(this%redshift_time,this%dredshift_time,nthermo,spline_data)
        call splder(this%arhos_fac,this%darhos_fac,nthermo,spline_data)
        call splder(this%darhos_fac,this%ddarhos_fac,nthermo,spline_data)
        do j2 = 1, State%num_redshiftwindows
            if (State%Redshift_w(j2)%has_lensing_window) then
                call splder(RW(j2)%awin_lens,RW(j2)%dawin_lens,nthermo,spline_data)
            end if
        end do
    end if
    call this%SetTimeSteps(State,State%TimeSteps)
    !$OMP END PARALLEL SECTIONS

    if (State%num_redshiftwindows>0) then
        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
        do j2=1,State%TimeSteps%npoints
            call this%DoWindowSpline(State,j2,State%TimeSteps%points(j2), RW)
        end do
        !$OMP END PARALLEL DO
        call this%SetTimeStepWindows(State,State%TimeSteps)
    end if
    ! and fix segmentation
    call etat%Init(State,WantTSpin=CP%Do21cm)
    ! and
    if (CP%WantDerivedParameters) then
        associate(ThermoDerivedParams => State%ThermoDerivedParams)
            !$OMP PARALLEL SECTIONS DEFAULT(SHARED)
            !$OMP SECTION
            ThermoDerivedParams( derived_Age ) = State%DeltaPhysicalTimeGyr(0.0_dl,1.0_dl)
            rstar =State%sound_horizon(this%z_star)
            ThermoDerivedParams( derived_rstar ) = rstar
            DA = State%AngularDiameterDistance(this%z_star)/(1/(this%z_star+1))
            ThermoDerivedParams( derived_zdrag ) = this%z_drag
            !$OMP SECTION
            rs =State%sound_horizon(this%z_drag)
            ThermoDerivedParams( derived_rdrag ) = rs
            ThermoDerivedParams( derived_kD ) =  &
                sqrt(1.d0/(Integrate_Romberg(State,ddamping_da, 1d-8, 1/(this%z_star+1), 1d-6)/6))
            !$OMP SECTION
            ThermoDerivedParams( derived_zEQ ) = State%z_eq
            a_eq = 1/(1+State%z_eq)
            ThermoDerivedParams( derived_kEQ ) = 1/(a_eq*dtauda(State,a_eq))
            rs_eq = State%sound_horizon(State%z_eq)
            tau_eq = State%timeOfz(State%z_eq)
            !$OMP SECTION
            !$OMP END PARALLEL SECTIONS

            ThermoDerivedParams( derived_zstar ) = this%z_star
            ThermoDerivedParams( derived_thetastar ) = 100*rstar/DA
            ThermoDerivedParams( derived_DAstar ) = DA/1000
            ThermoDerivedParams( derived_thetaEQ ) = 100*tau_eq/DA
            ThermoDerivedParams( derived_theta_rs_EQ ) = 100*rs_EQ/DA
            ThermoDerivedParams( derived_thetaD ) =  100*const_pi/ThermoDerivedParams( derived_kD )/DA

            if (allocated(CP%z_outputs)) then
                if (allocated(State%BackgroundOutputs%H)) &
                    deallocate(State%BackgroundOutputs%H, State%BackgroundOutputs%DA, State%BackgroundOutputs%rs_by_D_v)
                    noutput = size(CP%z_outputs)
                    allocate(State%BackgroundOutputs%H(noutput), State%BackgroundOutputs%DA(noutput), &
                    State%BackgroundOutputs%rs_by_D_v(noutput))
                !$OMP PARALLEL DO DEFAULT(shared)
                    do i=1,noutput
                        State%BackgroundOutputs%H(i) = State%HofZ(CP%z_outputs(i))
                        State%BackgroundOutputs%DA(i) = State%AngularDiameterDistance(CP%z_outputs(i))
                        State%BackgroundOutputs%rs_by_D_v(i) = rs/BAO_D_v_from_DA_H(CP%z_outputs(i), &
                        State%BackgroundOutputs%DA(i),State%BackgroundOutputs%H(i))
                    end do
                end if

                if (FeedbackLevel > 0) then
                    write(*,'("Age of universe/GYr  = ",f7.3)') ThermoDerivedParams( derived_Age )
                    write(*,'("zstar                = ",f8.2)') ThermoDerivedParams( derived_zstar )
                    write(*,'("r_s(zstar)/Mpc       = ",f7.2)') ThermoDerivedParams( derived_rstar )
                    write(*,'("100*theta            = ",f9.6)') ThermoDerivedParams( derived_thetastar )
                    write(*,'("DA(zstar)/Gpc        = ",f9.5)') ThermoDerivedParams( derived_DAstar )

                    write(*,'("zdrag                = ",f8.2)') ThermoDerivedParams( derived_zdrag )
                    write(*,'("r_s(zdrag)/Mpc       = ",f7.2)') ThermoDerivedParams( derived_rdrag )

                    write(*,'("k_D(zstar) Mpc       = ",f7.4)') ThermoDerivedParams( derived_kD )
                    write(*,'("100*theta_D          = ",f9.6)') ThermoDerivedParams( derived_thetaD )

                    write(*,'("z_EQ (if v_nu=1)     = ",f8.2)') ThermoDerivedParams( derived_zEQ )
                    write(*,'("k_EQ Mpc (if v_nu=1) = ",f9.6)') ThermoDerivedParams( derived_kEQ )
                    write(*,'("100*theta_EQ         = ",f9.6)') ThermoDerivedParams( derived_thetaEQ )
                    write(*,'("100*theta_rs_EQ      = ",f9.6)') ThermoDerivedParams( derived_theta_rs_EQ )
                end if
        end associate
    end if

    !Sources
    if(State%num_redshiftwindows>0) then
        deallocate(RW,this%arhos_fac, this%darhos_fac, this%ddarhos_fac)
        deallocate(this%redshift_time, this%dredshift_time)
        do RW_i = 1, State%num_redshiftwindows
            associate (RedWin => State%Redshift_W(RW_i))
                if (RedWin%kind == window_21cm) then
                    outstr = 'z= '//trim(RealToStr(real(RedWin%Redshift),4))//': T_b = '//trim(RealToStr(real(RedWin%Fq),6))// &
                        'mK; tau21 = '//trim(RealToStr(real(RedWin%optical_depth_21),5))
                    write (*,*) RW_i,trim(outstr)
                end if
            end associate
        end do
    end if

    if (CP%Do21cm .and. CP%transfer_21cm_cl) then
        do RW_i=1,CP%Transfer%PK_num_redshifts
            a=1._dl/(1+CP%Transfer%PK_redshifts(RW_i))
            Tspin = CP%Recomb%T_s(a)
            Trad = CP%TCMB/a
            adot = 1/dtauda(State,a)
            tau_eps = a**2*line21_const*State%NNow/a**3/adot/Tspin/1000
            Tb21cm = 1000*(1-exp(-tau_eps))*a*(Tspin-Trad)
            if (FeedbackLevel>0) then
                outstr = 'z= '//trim(RealToStr(real(CP%Transfer%PK_redshifts(RW_i))))// &
                    ': tau_21cm = '//trim(RealToStr(real(tau_eps),5))//'; T_b = '//trim(RealToStr(real(Tb21cm),6))//'mK'
                write (*,*) trim(outstr)
            end if
        end do
    end if

    this%HasThermoData = .true.

    ! Deallocate memory for Statee
    ! deallocate(Statee)

    end subroutine Thermo_Init

    ! subroutine thermoarr(tau,cs2b,opacity, dopacity)
    ! !Compute unperturbed sound speed squared,
    ! !and ionization fraction by interpolating pre-computed tables.
    ! !If requested also get time derivative of opacity
    ! implicit none
    ! real(dl) tau,cs2b,opacity(nscatter)
    ! real(dl), intent(out), optional :: dopacity(nscatter)

    ! integer i
    ! real(dl) d

    ! d=log(tau/tauminn)/dlntau+1._dl
    ! i=int(d)
    ! d=d-i
    ! if (i < 1) then
    !     !Linear interpolation if out of bounds (should not occur).
    !     cs2b=cs2(1)+(d+i-1)*dcs2(1)
    !     opacity=this%dotmu(1,:)+(d-1)*ddotmu(1,:)
    !     stop 'thermo out of bounds'
    ! else if (i >= nthermo) then
    !     cs2b=cs2(nthermo)+(d+i-nthermo)*dcs2(nthermo)
    !     opacity=this%dotmu(nthermo,:)+(d-nthermo)*ddotmu(nthermo,:)
    !     if (present(dopacity)) then
    !         dopacity = 0
    !         stop 'thermo: shouldn''t happen'
    !     end if
    ! else
    !     !Cubic spline interpolation.
    !     cs2b=cs2(i)+d*(dcs2(i)+d*(3*(cs2(i+1)-cs2(i))  &
    !     -2*dcs2(i)-dcs2(i+1)+d*(dcs2(i)+dcs2(i+1)  &
    !     +2*(cs2(i)-cs2(i+1)))))
    !     opacity=dthis%otmu(i,:)+d*(ddotmu(i,:)+d*(3*(this%dotmu(i+1,:)-this%dotmu(i,:)) &
    !     -2*ddotmu(i,:)-ddotmu(i+1,:)+d*(ddotmu(i,:)+ddotmu(i+1,:) &
    !     +2*(this%dotmu(i,:)-this%dotmu(i+1,:)))))

    !     if (present(dopacity)) then
    !         dopacity=(ddotmu(i,:)+d*(dddotmu(i,:)+d*(3*(ddotmu(i+1,:)  &
    !         -ddotmu(i,:))-2*dddotmu(i,:)-dddotmu(i+1,:)+d*(dddotmu(i,:) &
    !         +dddotmu(i+1,:)+2*(ddotmu(i,:)-ddotmu(i+1,:))))))/(tau*dlntau)
    !     end if
    ! end if
    ! end subroutine thermoarr

    subroutine SetTimeSteps(this,State,TimeSteps)
    !Set time steps to use for sampling the source functions for the CMB power spectra
    class(TThermoData) :: this
    Type(TRanges) :: TimeSteps
    class(CAMBdata) State
    real(dl) dtau0
    integer nri0, nstep
    !Sources
    integer ix,i,nwindow, L_limb
    real(dl) keff, win_end, TimeSampleBoost, delta

    TimeSampleBoost = State%CP%Accuracy%AccuracyBoost*State%CP%Accuracy%TimeStepBoost
    call TimeSteps%Init()

    call TimeSteps%Add_delta(State%taurst, State%taurend, State%dtaurec)

    ! Calculating the timesteps after recombination
    if (State%CP%WantTensors) then
        dtau0=max(State%taurst/40,State%tau0/2000._dl/TimeSampleBoost)
    else
        dtau0=State%tau0/500._dl/TimeSampleBoost
        if (do_bispectrum) dtau0 = dtau0/3
        !Don't need this since adding in Limber on small scales
        !  if (CP%DoLensing) dtau0=dtau0/2
        !  if (CP%AccurateBB) dtau0=dtau0/3 !Need to get C_Phi accurate on small scales
    end if

    call TimeSteps%Add_delta(State%taurend, State%tau0, dtau0)

    !Sources

    this%tau_start_redshiftwindows = State%tau0
    this%tau_end_redshiftwindows = 0
    if (State%CP%WantScalars .or. State%CP%SourceTerms%line_reionization) then
        do ix=1, State%num_redshiftwindows
            associate (Win => State%Redshift_W(ix))

                Win%tau_start = max(Win%tau_start,state%taurst)

                if (Win%kind /= window_lensing) then
                    !Have to be careful to integrate dwinV as the window tails off
                    this%tau_end_redshiftwindows = max(Win%tau_end,this%tau_end_redshiftwindows)
                    nwindow = nint(150*TimeSampleBoost)
                    win_end = Win%tau_end
                else !lensing
                    nwindow = nint(TimeSampleBoost*Win%chi0/100)
                    win_end = State%tau0
                end if

                if (Win%kind == window_21cm .and. (State%CP%SourceTerms%line_phot_dipole .or. &
                    State%CP%SourceTerms%line_phot_quadrupole)) nwindow = nwindow *3

                L_limb = Win_limber_ell(Win,State%CP,State%CP%max_l)
                keff = WindowKmaxForL(Win,State%CP,L_limb)

                !Keep sampling in x better than Nyquist
                nwindow = max(nwindow, nint(TimeSampleBoost *(win_end- Win%tau_start)* keff/3))
                if (Win%kind /= window_lensing .and. Win%sigma_tau < (win_end - Win%tau_start)/nwindow*8) then
                    nwindow = nint(TimeSampleBoost*(win_end - Win%tau_start)/Win%sigma_tau*8)
                end if
                if (Feedbacklevel > 1 .or. DebugMsgs) call WriteFormat('nwindow %d: %d', ix, nwindow)
                delta = (win_end-Win%tau_start)/nwindow
                Win%tau_start = Win%tau_start - delta*3
                win_end = min(State%tau0, win_end+delta*2)

                this%tau_start_redshiftwindows = min(Win%tau_start,this%tau_start_redshiftwindows)

                call TimeSteps%Add(Win%tau_start, win_end, nwindow)
                !This should cover whole range where not tiny

                if (Win%kind /= window_lensing  .and. &
                    Win%tau_peakend-Win%tau_peakstart < nint(60*TimeSampleBoost) * delta) then
                    call TimeSteps%Add(Win%tau_peakstart,Win%tau_peakend, nint(60*TimeSampleBoost))
                    !This should be over peak
                end if
            end associate
        end do
    end if


    if (State%CP%Reion%Reionization) then
        nri0=int(State%reion_n_steps*State%CP%Accuracy%AccuracyBoost)
        !Steps while reionization going from zero to maximum
        call TimeSteps%Add(State%reion_tau_start,State%reion_tau_complete,nri0)
    end if

    !Sources
    if (.not. State%CP%Want_CMB .and. State%CP%WantCls) then
        if (State%num_redshiftwindows==0) then
            call GlobalError('Want_CMB=false and WantCls=true, but no redshift windows either', &
                error_unsupported_params)
        else
            call TimeSteps%Add_delta(this%tau_start_redshiftwindows, State%tau0, dtau0)
        endif
    end if

    if (global_error_flag/=0) then
        return
    end if



    !Create arrays out of the region information.
    call TimeSteps%GetArray()
    nstep = TimeSteps%npoints

    if (State%num_redshiftwindows > 0) then
        if (allocated(this%step_redshift)) deallocate(this%step_redshift, this%rhos_fac, this%drhos_fac)
        allocate(this%step_redshift(nstep), this%rhos_fac(nstep), this%drhos_fac(nstep))
        do i=1,State%num_redshiftwindows
            associate (Win => State%Redshift_W(i))
                allocate(Win%winF(nstep),Win%wing(nstep),Win%dwing(nstep),Win%ddwing(nstep), &
                    Win%winV(nstep),Win%dwinV(nstep),Win%ddwinV(nstep))
                allocate(Win%win_lens(nstep),Win%wing2(nstep),Win%dwing2(nstep),Win%ddwing2(nstep))
                allocate(Win%wingtau(nstep),Win%dwingtau(nstep),Win%ddwingtau(nstep))
                if (Win%kind == window_counts) then
                    allocate(Win%comoving_density_ev(nstep))
                end if
            end associate
        end do
    end if

    if (DebugMsgs .and. FeedbackLevel > 0) call WriteFormat('Set %d time steps', nstep)

    end subroutine SetTimeSteps

    subroutine SetTimeStepWindows(this,State,TimeSteps)
    use constants
    class(TThermoData) :: this
    class(CAMBdata) :: State
    Type(TRanges) :: TimeSteps
    integer i, j, jstart, ix
    real(dl) tau, a, a2
    real(dl) Tspin, Trad, rho_fac, tau_eps
    real(dl) window, winamp
    real(dl) z,rhos, adot, exp_fac
    real(dl) tmp(TimeSteps%npoints), tmp2(TimeSteps%npoints), hubble_tmp(TimeSteps%npoints)
    real(dl), allocatable , dimension(:,:) :: int_tmp, back_count_tmp
    integer ninterp

    ! Prevent false positive warnings for uninitialized
    Tspin = 0._dl
    Trad = 0._dl
    tau_eps = 0._dl
    exp_fac = 0._dl

    jstart = TimeSteps%IndexOf(this%tau_start_redshiftwindows)
    ninterp = TimeSteps%npoints - jstart + 1

    do i = 1, State%num_redshiftwindows
        associate (RedWin => State%Redshift_W(i))
            RedWin%wing=0
            RedWin%winV=0
            RedWin%winF=0
            RedWin%wing2=0
            RedWin%dwing=0
            RedWin%dwinV=0
            RedWin%dwing2=0
            RedWin%ddwing=0
            RedWin%ddwinV=0
            RedWin%ddwing2=0
            RedWin%wingtau=0
            RedWin%dwingtau=0
            RedWin%ddwingtau=0
            RedWin%Fq = 0
            if (RedWin%kind == window_counts) then
                RedWin%comoving_density_ev  = 0
            end if
        end associate
    end do

    allocate(int_tmp(jstart:TimeSteps%npoints,State%num_redshiftwindows))
    int_tmp = 0
    allocate(back_count_tmp(jstart:TimeSteps%npoints,State%num_redshiftwindows))
    back_count_tmp = 0

    do j=jstart, TimeSteps%npoints
        tau = TimeSteps%points(j)
        z = this%step_redshift(j)
        a = 1._dl/(1._dl+z)
        a2=a**2
        adot=1._dl/dtauda(State,a)


        if (State%CP%Do21cm) then
            Tspin = State%CP%Recomb%T_s(a)
            Trad = State%CP%TCMB/a
            rho_fac = line21_const*State%NNow/a**3 !neglect ionization fraction
            tau_eps = a*rho_fac/(adot/a)/Tspin/1000
            exp_fac =  (1-exp(-tau_eps))/tau_eps
        else
            rho_fac = State%grhoc/a**3
        end if

        hubble_tmp(j) = adot/a

        do i = 1, State%num_redshiftwindows
            associate (RedWin => State%Redshift_W(i))
                if (tau < RedWin%tau_start) cycle

                window = RedWin%Window%Window_f_a(a, winamp)

                if (RedWin%kind == window_21cm) then
                    rhos = rho_fac*(1 - Trad/Tspin)

                    !Want to integrate this...
                    int_tmp(j,i) = this%drhos_fac(j)*a*window

                    RedWin%WinV(j) = -exp(-tau_eps)*a2*rhos*window/(adot/a)

                    RedWin%wing(j) = exp_fac*a2*rhos*window

                    !The window that doesn't go to zero at T_s = T_gamma
                    RedWin%wing2(j) = exp_fac*a2*rho_fac*Trad/Tspin*window

                    !Window for tau_s for the self-absoption term
                    RedWin%wingtau(j) =  RedWin%wing(j)*(1 - exp(-tau_eps)/exp_fac)
                elseif (RedWin%kind == window_counts) then

                    !window is n(a) where n is TOTAL not fractional number
                    !delta = int wing(eta) delta(eta) deta
                    RedWin%wing(j) = adot *window

                    !Window with 1/H in
                    RedWin%wing2(j) = RedWin%wing(j)/(adot/a)

                    !winv is g/chi for the ISW and time delay terms
                    RedWin%WinV(j) = 0
                    if (tau < State%tau0 -0.1) then
                        int_tmp(j,i) = RedWin%wing(j)/(State%tau0 - tau)
                    else
                        int_tmp(j,i)=0
                    end if

                    if (State%CP%SourceTerms%counts_evolve) then
                        back_count_tmp(j,i) =  RedWin%Window%counts_background_z(1/a-1)/a
                        if (tau < State%tau0 -0.1) then
                            RedWin%comoving_density_ev(j) = back_count_tmp(j,i)*(adot/a)/(State%tau0 - tau)**2
                        else
                            RedWin%comoving_density_ev(j) = 0
                        end if
                    end if
                end if
            end associate
        end do
    end do

    do i = 1, State%num_redshiftwindows
        associate (RedWin => State%Redshift_W(i))

            ! int (a*rho_s/H)' a W_f(a) d\eta, or for counts int g/chi deta
            call spline_def(TimeSteps%points(jstart:),int_tmp(jstart:,i),ninterp,tmp)
            call spline_integrate(TimeSteps%points(jstart:),int_tmp(jstart:,i),tmp, tmp2(jstart:),ninterp)
            RedWin%WinV(jstart:TimeSteps%npoints) =  &
                RedWin%WinV(jstart:TimeSteps%npoints) + tmp2(jstart:TimeSteps%npoints)

            call spline_def(TimeSteps%points(jstart:),RedWin%WinV(jstart:),ninterp,RedWin%ddWinV(jstart:))
            call spline_deriv(TimeSteps%points(jstart:),RedWin%WinV(jstart:),RedWin%ddWinV(jstart:), RedWin%dWinV(jstart:), ninterp)

            call spline_def(TimeSteps%points(jstart:),RedWin%Wing(jstart:),ninterp,RedWin%ddWing(jstart:))
            call spline_deriv(TimeSteps%points(jstart:),RedWin%Wing(jstart:),RedWin%ddWing(jstart:), RedWin%dWing(jstart:), ninterp)

            call spline_def(TimeSteps%points(jstart:),RedWin%Wing2(jstart:),ninterp,RedWin%ddWing2(jstart:))
            call spline_deriv(TimeSteps%points(jstart:),RedWin%Wing2(jstart:),RedWin%ddWing2(jstart:), &
                RedWin%dWing2(jstart:), ninterp)

            call spline_integrate(TimeSteps%points(jstart:),RedWin%Wing(jstart:), &
                RedWin%ddWing(jstart:), RedWin%WinF(jstart:),ninterp)
            RedWin%Fq = RedWin%WinF(TimeSteps%npoints)

            if (RedWin%kind == window_21cm) then
                call spline_integrate(TimeSteps%points(jstart:),RedWin%Wing2(jstart:),&
                    RedWin%ddWing2(jstart:), tmp(jstart:),ninterp)
                RedWin%optical_depth_21 = tmp(TimeSteps%npoints) / (State%CP%TCMB*1000)
                !WinF not used.. replaced below

                call spline_def(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),ninterp,RedWin%ddWingtau(jstart:))
                call spline_deriv(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),RedWin%ddWingtau(jstart:), &
                    RedWin%dWingtau(jstart:), ninterp)
            elseif (RedWin%kind == window_counts) then

                if (State%CP%SourceTerms%counts_evolve) then
                    call spline_def(TimeSteps%points(jstart:),back_count_tmp(jstart:,i),ninterp,tmp)
                    call spline_deriv(TimeSteps%points(jstart:),back_count_tmp(jstart:,i),tmp,tmp2(jstart:),ninterp)
                    do ix = jstart, TimeSteps%npoints
                        if (RedWin%Wing(ix)==0._dl) then
                            RedWin%Wingtau(ix) = 0
                        else
                            !evo bias is computed with total derivative
                            RedWin%Wingtau(ix) =  -tmp2(ix) * RedWin%Wing(ix) / (back_count_tmp(ix,i)*hubble_tmp(ix)) &
                                !+ 5*RedWin%dlog10Ndm * ( RedWin%Wing(ix)- int_tmp(ix,i)/hubble_tmp(ix))
                                !The correction from total to partial derivative takes 1/adot(tau0-tau) cancels
                                + 10*RedWin%Window%dlog10Ndm * RedWin%Wing(ix)
                        end if
                    end do

                    !comoving_density_ev is d log(a^3 n_s)/d eta * window
                    call spline_def(TimeSteps%points(jstart:),RedWin%comoving_density_ev(jstart:),ninterp,tmp)
                    call spline_deriv(TimeSteps%points(jstart:),RedWin%comoving_density_ev(jstart:),tmp,tmp2(jstart:),ninterp)
                    do ix = jstart, TimeSteps%npoints
                        if (RedWin%Wing(ix)==0._dl) then
                            RedWin%comoving_density_ev(ix) = 0
                        elseif (RedWin%comoving_density_ev(ix)/=0._dl) then
                            !correction needs to be introduced from total derivative to partial derivative
                            RedWin%comoving_density_ev(ix) =   tmp2(ix) / RedWin%comoving_density_ev(ix) &
                                -5*RedWin%Window%dlog10Ndm * ( hubble_tmp(ix) + int_tmp(ix,i)/RedWin%Wing(ix))
                        end if
                    end do
                else
                    RedWin%comoving_density_ev=0
                    call spline_def(TimeSteps%points(jstart:),hubble_tmp(jstart:),ninterp,tmp)
                    call spline_deriv(TimeSteps%points(jstart:),hubble_tmp(jstart:),tmp, tmp2(jstart:), ninterp)

                    !assume d( a^3 n_s) of background population is zero, so remaining terms are
                    !wingtau =  g*(2/H\chi + Hdot/H^2)  when s=0; int_tmp = window/chi
                    RedWin%Wingtau(jstart:TimeSteps%npoints) = &
                        2*(1-2.5*RedWin%Window%dlog10Ndm)*int_tmp(jstart:TimeSteps%npoints,i)/&
                        hubble_tmp(jstart:TimeSteps%npoints)&
                        + 5*RedWin%Window%dlog10Ndm*RedWin%Wing(jstart:TimeSteps%npoints) &
                        + tmp2(jstart:TimeSteps%npoints)/hubble_tmp(jstart:TimeSteps%npoints)**2 &
                        *RedWin%Wing(jstart:TimeSteps%npoints)
                endif

                call spline_def(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),ninterp, &
                    RedWin%ddWingtau(jstart:))
                call spline_deriv(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),RedWin%ddWingtau(jstart:), &
                    RedWin%dWingtau(jstart:), ninterp)

                !WinF is int[ g*(...)]
                call spline_integrate(TimeSteps%points(jstart:),RedWin%Wingtau(jstart:),&
                    RedWin%ddWingtau(jstart:), RedWin%WinF(jstart:),ninterp)
            end if
        end associate
    end do

    end subroutine SetTimeStepWindows


    subroutine interp_window(TimeSteps,RedWin,tau,wing_t, wing2_t, winv_t)
    !for evolving sources for reionization we neglect wingtau self-absorption
    Type(TRanges) :: TimeSteps
    Type(TRedWin)  :: RedWin
    integer i
    real(dl) :: tau, wing_t, wing2_t,winv_t
    real(dl) a0,b0,ho

    i = TimeSteps%IndexOf(tau)
    if (i< TimeSteps%npoints) then
        ho=TimeSteps%points(i+1)-TimeSteps%points(i)
        a0=(TimeSteps%points(i+1)-tau)/ho
        b0=1-a0
        wing_t = a0*RedWin%wing(i)+ b0*RedWin%wing(i+1)+((a0**3-a0)* RedWin%ddwing(i) &
            +(b0**3-b0)*RedWin%ddwing(i+1))*ho**2/6
        wing2_t = a0*RedWin%wing2(i)+ b0*RedWin%wing2(i+1)+((a0**3-a0)* RedWin%ddwing2(i) &
            +(b0**3-b0)*RedWin%ddwing2(i+1))*ho**2/6
        winv_t = a0*RedWin%winv(i)+ b0*RedWin%winv(i+1)+((a0**3-a0)* RedWin%ddwinv(i) &
            +(b0**3-b0)*RedWin%ddwinv(i+1))*ho**2/6
    else
        wing_t = 0
        wing2_t = 0
        winv_t = 0
    end if
    end subroutine interp_window

    subroutine DoWindowSpline(this,State,j2,tau, RW)
    class(TThermoData) :: this
    class(CAMBdata) :: State
    integer j2, i, RW_i
    real(dl) d, tau
    ! andrea
    integer scat
    ! 
    Type(CalWins) :: RW(:)

    !     Cubic-spline interpolation.
    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)
    d=d-i
    if (i < this%nthermo) then
        this%step_redshift(j2) = this%redshift_time(i)+d*(this%dredshift_time(i)+ &
            d*(3._dl*(this%redshift_time(i+1)-this%redshift_time(i)) &
            -2._dl*this%dredshift_time(i)-this%dredshift_time(i+1)+d*(this%dredshift_time(i)+this%dredshift_time(i+1) &
            +2._dl*(this%redshift_time(i)-this%redshift_time(i+1)))))

        this%rhos_fac(j2) = this%arhos_fac(i)+d*(this%darhos_fac(i)+d*(3._dl*(this%arhos_fac(i+1)-this%arhos_fac(i)) &
            -2._dl*this%darhos_fac(i)-this%darhos_fac(i+1)+d*(this%darhos_fac(i)+this%darhos_fac(i+1) &
            +2._dl*(this%arhos_fac(i)-this%arhos_fac(i+1)))))
        this%drhos_fac(j2) = (this%darhos_fac(i)+d*(this%ddarhos_fac(i)+d*(3._dl*(this%darhos_fac(i+1)  &
            -this%darhos_fac(i))-2._dl*this%ddarhos_fac(i)-this%ddarhos_fac(i+1)+d*(this%ddarhos_fac(i) &
            +this%ddarhos_fac(i+1)+2._dl*(this%darhos_fac(i)-this%darhos_fac(i+1))))))/(tau &
            *this%dlntau)

        do RW_i=1, State%num_redshiftwindows
            if (State%Redshift_w(RW_i)%has_lensing_window) then
                associate(W => State%Redshift_W(RW_i), C=> RW(RW_i))

                    W%win_lens(j2) = C%awin_lens(i)+d*(C%dawin_lens(i)+d*(3._dl*(C%awin_lens(i+1)-C%awin_lens(i)) &
                        -2._dl*C%dawin_lens(i)-C%dawin_lens(i+1)+d*(C%dawin_lens(i)+C%dawin_lens(i+1) &
                        +2._dl*(C%awin_lens(i)-C%awin_lens(i+1)))))
                end associate
            end if
        end do

    else
        this%step_redshift(j2) = 0
        this%rhos_fac(j2)=0
        this%drhos_fac(j2)=0
        do RW_i=1, State%num_redshiftwindows
            associate (W => State%Redshift_W(RW_i))
                W%win_lens(j2)=0
            end associate
        end do
    end if

    end subroutine DoWindowSpline

    subroutine IonizationFunctionsAtTime(this,tau, a, opacity, dopacity, ddopacity, &
        vis, dvis, ddvis, exptau, lenswin)
    class(TThermoData) :: this
    real(dl), intent(in) :: tau
    real(dl), intent(out) :: a
    real(dl), intent(out) :: vis(nscatter), dvis(nscatter), ddvis(nscatter), &
                             exptau(nscatter), opacity(nscatter), dopacity(nscatter), ddopacity(nscatter)
    real(dl), intent(out) :: lenswin 
    real(dl) d, cs2
    integer :: i, scat

    call this%values_array(tau,a,cs2,opacity,dopacity)

    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)
    d=d-i
    ! andrea big change
    do scat = 1, nscatter
    ! 
        if (i < this%nthermo) then
            ! andrea change
            ddopacity(scat) = (this%dddotmu(i, scat) + d * (this%ddddotmu(i, scat) + d * (3.0_dl * (this%dddotmu(i + 1, scat) &
                - this%dddotmu(i, scat)) - 2.0_dl * this%ddddotmu(i, scat) - this%ddddotmu(i + 1, scat) &
                + d * (this%ddddotmu(i, scat) + this%ddddotmu(i + 1, scat) + 2.0_dl * (this%dddotmu(i, scat) &
                - this%dddotmu(i + 1, scat))))) - (this%dlntau**2) * tau * dopacity(scat)) / &
                ((tau * this%dlntau)**2)
            exptau(scat)=this%emmu(i,scat)+d*(this%demmu(i,scat)+d*(3._dl*(this%emmu(i+1,scat)-this%emmu(i,scat)) &
                -2._dl*this%demmu(i,scat)-this%demmu(i+1,scat)+d*(this%demmu(i,scat)+this%demmu(i+1,scat) &
            +   2._dl*(this%emmu(i,scat)-this%emmu(i+1,scat)))))
            ! 
            if (dowinlens) then
                lenswin=this%winlens(i)+d*(this%dwinlens(i)+d*(3._dl*(this%winlens(i+1)-this%winlens(i)) &
                    -2._dl*this%dwinlens(i)-this%dwinlens(i+1)+d*(this%dwinlens(i)+this%dwinlens(i+1) &
                    +2._dl*(this%winlens(i)-this%winlens(i+1)))))
            end if
            vis(scat)=opacity(scat)*exptau(scat)
            dvis(scat)=exptau(scat)*(opacity(scat)**2+dopacity(scat))
            ddvis(scat)=exptau(scat)*(opacity(scat)**3+3*opacity(scat)*dopacity(scat)+ddopacity(scat))
        else
            ! andrea solo change
            ddopacity(scat)=this%dddotmu(this%nthermo,scat)
            exptau(scat)=this%emmu(this%nthermo,scat)
            vis(scat)=opacity(scat)*exptau(scat)
            dvis(scat)=exptau(scat)*(opacity(scat)**2+dopacity(scat))
            ddvis(scat)=exptau(scat)*(opacity(scat)**3+3._dl*opacity(scat)*dopacity(scat)+ddopacity(scat))
            ! 
        end if
    ! 
    end do
    ! 
    end subroutine IonizationFunctionsAtTime

    subroutine Init_ClTransfer(CTrans)
    !Need to set the Ranges array q before calling this
    Type(ClTransferData) :: CTrans
    integer st

    deallocate(CTrans%Delta_p_l_k, STAT = st)
    call CTrans%q%getArray(.true.)

    allocate(CTrans%Delta_p_l_k(CTrans%NumSources,&
        min(CTrans%max_index_nonlimber,CTrans%ls%nl), CTrans%q%npoints), STAT = st)
    if (st /= 0) call MpiStop('Init_ClTransfer: Error allocating memory for transfer functions')
    CTrans%Delta_p_l_k = 0

    end subroutine Init_ClTransfer

    subroutine Init_Limber(CTrans,State)
    Type(ClTransferData) :: CTrans
    class(CAMBdata) :: State

    call Free_Limber(Ctrans)
    allocate(CTrans%Limber_l_min(CTrans%NumSources))
    CTrans%Limber_l_min = 0
    if (State%num_redshiftwindows>0 .or. State%CP%SourceTerms%limber_phi_lmin>0) then
        allocate(CTrans%Limber_windows(CTrans%NumSources,CTrans%ls%nl))
    end if

    end subroutine Init_Limber

    subroutine Free_ClTransfer(CTrans)
    Type(ClTransferData) :: CTrans

    if (allocated(CTrans%Delta_p_l_k)) deallocate(CTrans%Delta_p_l_k)
    call CTrans%q%Free()
    call Free_Limber(CTrans)

    end subroutine Free_ClTransfer

    subroutine Free_Limber(CTrans)
    Type(ClTransferData) :: CTrans

    if (allocated(CTrans%Limber_l_min)) deallocate(CTrans%Limber_l_min)
    if (allocated(CTrans%Limber_windows)) deallocate(CTrans%Limber_windows)

    end subroutine Free_Limber


    function Win_Limber_ell(W,CP,lmax) result(ell_limb)
    Type(TRedWin) :: W
    Type(CAMBParams) :: CP
    integer, intent(in) :: lmax
    integer ell_limb
    real(dl) LimBoost

    if (CP%SourceTerms%limber_windows) then
        LimBoost = CP%Accuracy%AccuracyBoost*CP%Accuracy%SourceLimberBoost
        !Turn on limber when k is a scale smaller than window width
        if (W%kind==window_lensing) then
            ell_limb = max(CP%SourceTerms%limber_phi_lmin,nint(50*LimBoost))
        else
            ell_limb = max(CP%SourceTerms%limber_phi_lmin, nint(LimBoost*6*W%chi0/W%sigma_tau))
        end if
    else
        ell_limb = lmax
    end if
    end function Win_Limber_ell


    subroutine TCLdata_InitCls(this, State)
    class(TCLData) :: this
    class(CAMBdata) :: State

    associate(CP=>State%CP)
        call CheckLoadedHighLTemplate
        if (CP%WantScalars) then
            if (allocated(this%Cl_scalar)) deallocate(this%Cl_scalar)
            allocate(this%Cl_scalar(CP%Min_l:CP%Max_l, C_Temp:State%Scalar_C_last), source=0._dl)
            if (CP%want_cl_2D_array) then
                if (allocated(this%Cl_scalar_array)) deallocate(this%Cl_scalar_array)
                ! andrea
                allocate(this%Cl_scalar_Array(CP%Min_l:CP%Max_l, &
                    3+State%num_redshiftwindows+num_cmb_freq*2+CP%CustomSources%num_custom_sources, &
                    3+State%num_redshiftwindows+num_cmb_freq*2+CP%CustomSources%num_custom_sources))
                this%Cl_scalar_array = 0
                ! andrea
            end if
        end if

        if (CP%WantVectors) then
            if (allocated(this%Cl_vector)) deallocate(this%Cl_vector)
            allocate(this%Cl_vector(CP%Min_l:CP%Max_l, CT_Temp:CT_Cross), source=0._dl)
        end if

        if (CP%WantTensors) then
            if (allocated(this%Cl_tensor)) deallocate(this%Cl_tensor)
            allocate(this%Cl_tensor(CP%Max_l_tensor, CT_Temp:CT_Cross), source=0._dl)
        end if
    end associate

    end subroutine TCLdata_InitCls

    function open_file_header(filename, Col1, Columns, n) result(unit)
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: col1
    character(LEN=name_tag_len), intent(in) :: Columns(:)
    integer, intent(in), optional :: n
    integer :: unit, nn

    nn = PresentDefault(6, n)

    open(newunit=unit,file=filename,form='formatted',status='replace')
    if (output_file_headers) then
        write(unit,'("#",1A'//Trim(IntToStr(nn-1))//'," ",*(A15))') Col1, Columns
    end if

    end function open_file_header


    function scalar_fieldname(i)
    integer, intent(in) :: i
    character(LEN=5) :: scalar_fieldname
    character(LEN=3), parameter :: scalar_fieldnames = 'TEP'

    if (i<=3) then
        scalar_fieldname = scalar_fieldnames(i:i)
    else
        scalar_fieldname = 'W'//trim(IntToStr(i-3))
    end if

    end function scalar_fieldname

    subroutine TCLdata_output_cl_files(this, State, ScalFile,ScalCovFile,TensFile, TotFile, LensFile, LensTotFile, factor)
    class(TCLData) :: this
    class(CAMBdata), target :: State
    character(LEN=*) ScalFile, TensFile, TotFile, LensFile, LensTotFile,ScalCovfile
    real(dl), intent(in), optional :: factor
    real(dl) :: fact
    integer :: last_C, il, i, j, unit
    real(dl), allocatable :: outarr(:,:)
    character(LEN=name_tag_len) :: cov_names((3+State%num_redshiftwindows)**2)
    Type(CAMBParams), pointer :: CP
    integer lmin

    ! andrea
    integer nsource_out
    integer f_i_1,f_i_2
    ! andrea

    CP=> State%CP
    lmin= CP%Min_l

    fact = PresentDefault(1._dl, factor)

    if (CP%WantScalars .and. ScalFile /= '') then
        last_C=min(C_PhiTemp,State%Scalar_C_last)
        unit = open_file_header(ScalFile, 'L', C_name_tags(:last_C))
        do il=lmin,min(10000,CP%Max_l)
            write(unit,trim(numcat('(1I6,',last_C))//'E15.6)')il ,fact*this%Cl_scalar(il,C_Temp:last_C)
        end do
        do il=10100,CP%Max_l, 100
            write(unit,trim(numcat('(1E15.6,',last_C))//'E15.6)') real(il),&
                fact*this%Cl_scalar(il,C_Temp:last_C)
        end do
        close(unit)
    end if

    if (CP%WantScalars .and. CP%want_cl_2D_array .and. ScalCovFile /= '' .and. this%CTransScal%NumSources>2) then
        ! andrea
        nsource_out = 3+State%num_redshiftwindows+num_cmb_freq*2
        allocate(outarr(1:nsource_out,1:nsource_out))
        ! andrea
        do i=1, 3+State%num_redshiftwindows
            do j=1, 3+State%num_redshiftwindows
                cov_names(j + (i-1)*(3+State%num_redshiftwindows)) = trim(scalar_fieldname(i))//'x'//trim(scalar_fieldname(j))
            end do
        end do
        unit = open_file_header(ScalCovFile, 'L', cov_names)

        do il=lmin,min(10000,CP%Max_l)
            ! andrea
            outarr=this%Cl_scalar_array(il,1:3*nsource_out,1:nsource_out)
            ! andrea
            outarr(1:2,:)=sqrt(fact)*outarr(1:2,:)
            outarr(:,1:2)=sqrt(fact)*outarr(:,1:2)
            ! andrea
            outarr(:,4:3+num_cmb_freq*2)=sqrt(fact)*outarr(:,4:3+num_cmb_freq*2)
            outarr(4:3+num_cmb_freq*2,:)=sqrt(fact)*outarr(4:3+num_cmb_freq*2,:)
            write(unit,trim(numcat('(1I6,',(nsource_out)**2))//'E15.6)') il, real(outarr)
            ! andrea
        end do
        do il=10100,CP%Max_l, 100
            outarr=this%Cl_scalar_array(il,1:nsource_out,1:nsource_out)
            outarr(1:2,:)=sqrt(fact)*outarr(1:2,:)
            outarr(:,1:2)=sqrt(fact)*outarr(:,1:2)
            ! andrea
            outarr(:,4:3+num_cmb_freq*2)=sqrt(fact)*outarr(:,4:3+num_cmb_freq*2)
            outarr(4:3+num_cmb_freq*2,:)=sqrt(fact)*outarr(4:3+num_cmb_freq*2,:)
            ! andrea
            write(unit,trim(numcat('(1E15.5,',(nsource_out)**2))//'E15.6)') real(il), real(outarr)
        end do
        close(unit)
        deallocate(outarr)
    end if

    if (CP%WantTensors .and. TensFile /= '') then
        unit = open_file_header(TensFile, 'L', CT_name_tags)
        do il=lmin,CP%Max_l_tensor
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_tensor(il, CT_Temp:CT_Cross)
        end do
        close(unit)
        ! ! andrea
        ! if (num_cmb_freq>0) then
        !     open(unit=fileio_unit,file=trim(TensFile)//'_freqs',form='formatted',status='replace')
        !     write(fileio_unit,'('//trim(IntToStr(num_cmb_freq))//'E15.5)') phot_freqs
        !     close(fileio_unit)
        !     do f_i_1=1,num_cmb_freq
        !         do f_i_2=1,num_cmb_freq
        !             open(unit=fileio_unit,file=trim(TensFile)//trim(concat('_',f_i_1,'_',f_i_2)),form='formatted',status='replace')
        !             do in=1,CP%InitPower%nn
        !                 do il=lmin, CP%Max_l_tensor
        !                     write(fileio_unit,'(1I6,4E15.5)')il, fact*Cl_tensor_freqs(il, in, CT_Temp:CT_Cross,f_i_1,f_i_2)
        !                 end do
        !             end do
        !             close(fileio_unit)
        !         end do
        !     end do
        ! end if
        ! ! andrea
    end if

    if (CP%WantTensors .and. CP%WantScalars .and. TotFile /= '') then
        unit = open_file_header(TotFile, 'L', CT_name_tags)
        do il=lmin,CP%Max_l_tensor
            write(unit,'(1I6,4E15.6)')il, fact*(this%Cl_scalar(il, C_Temp:C_E)+ this%Cl_tensor(il,C_Temp:C_E)), &
                fact*this%Cl_tensor(il, CT_B), fact*(this%Cl_scalar(il, C_Cross) + this%Cl_tensor(il, CT_Cross))
        end do
        do il=CP%Max_l_tensor+1,CP%Max_l
            write(unit,'(1I6,4E15.6)')il ,fact*this%Cl_scalar(il,C_Temp:C_E), 0._dl, fact*this%Cl_scalar(il,C_Cross)
        end do
        close(unit)
    end if

    if (CP%WantScalars .and. CP%DoLensing .and. LensFile /= '') then
        unit = open_file_header(LensFile, 'L', CT_name_tags)
        do il=lmin, this%lmax_lensed
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_lensed(il, CT_Temp:CT_Cross)
        end do
        close(unit)
        ! ! andrea
        ! if (num_cmb_freq>0) then
        !     open(unit=fileio_unit,file=trim(LensFile)//'_freqs',form='formatted',status='replace')
        !     write(fileio_unit,'('//trim(IntToStr(num_cmb_freq))//'E15.5)') phot_freqs
        !     close(fileio_unit)
        !     do f_i_1=1,num_cmb_freq
        !         do f_i_2=1,num_cmb_freq
        !             open(unit=fileio_unit,file=trim(LensFile)//trim(concat('_',f_i_1,'_',f_i_2)),form='formatted',status='replace')
        !             do in=1,CP%InitPower%nn
        !                 do il=lmin, lmax_lensed
        !                     write(unit,'(1I6,4E15.6)')il, fact*this%Cl_lensed_freqs(il, in, CT_Temp:CT_Cross,f_i_1,f_i_2)
        !                 end do
        !             end do
        !             close(fileio_unit)
        !         end do
        !     end do
        ! end if
        ! ! andrea
    end if

    if (CP%WantScalars .and. CP%WantTensors .and. CP%DoLensing .and. LensTotFile /= '') then
        unit = open_file_header(LensTotFile, 'L', CT_name_tags)
        do il=lmin,min(CP%Max_l_tensor,this%lmax_lensed)
            write(unit,'(1I6,4E15.6)')il, fact*(this%Cl_lensed(il,CT_Temp:CT_Cross)+ &
                this%Cl_tensor(il, CT_Temp:CT_Cross))
        end do
        do il=min(CP%Max_l_tensor,this%lmax_lensed)+1,this%lmax_lensed
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_lensed(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if
    end subroutine TCLdata_output_cl_files

    subroutine TCLdata_output_lens_pot_files(this,CP, LensPotFile, factor)
    !Write out L TT EE BB TE PP PT PE where P is the lensing potential, all unlensed
    !This input supported by LensPix from 2010
    class(TCLdata) :: this
    Type(CAMBParams), intent(in) :: CP
    integer :: il, unit
    real(dl), intent(in), optional :: factor
    real(dl) fact, scale, BB, TT, TE, EE
    character(LEN=*), intent(in) :: LensPotFile
    !output file of dimensionless [l(l+1)]^2 C_phi_phi/2pi and [l(l+1)]^(3/2) C_phi_T/2pi
    !This is the format used by Planck_like but original LensPix uses scalar_output_file.

    !(Cl_scalar and scalar_output_file numbers are instead l^4 C_phi and l^3 C_phi
    ! - for historical reasons)

    fact = PresentDefault(1._dl, factor)

    if (CP%WantScalars .and. CP%DoLensing .and. LensPotFile/='') then
        unit = open_file_header(LensPotFile, 'L', lens_pot_name_tags)
        do il=CP%Min_l, min(10000,CP%Max_l)
            TT = this%Cl_scalar(il, C_Temp)
            EE = this%Cl_scalar(il, C_E)
            TE = this%Cl_scalar(il, C_Cross)
            if (CP%WantTensors .and. il <= CP%Max_l_tensor) then
                TT= TT+this%Cl_tensor(il, CT_Temp)
                EE= EE+this%Cl_tensor(il, CT_E)
                TE= TE+this%Cl_tensor(il, CT_Cross)
                BB= this%Cl_tensor(il, CT_B)
            else
                BB=0
            end if
            scale = (real(il+1)/il)**2/OutputDenominator !Factor to go from old l^4 factor to new

            write(unit,'(1I6,7E15.6)') il , fact*TT, fact*EE, fact*BB, fact*TE, scale*this%Cl_scalar(il,C_Phi),&
                (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*this%Cl_scalar(il,C_PhiTemp:C_PhiE)
        end do
        do il=10100,CP%Max_l, 100
            scale = (real(il+1)/il)**2/OutputDenominator
            write(unit,'(1E15.5,7E15.6)') real(il), fact*this%Cl_scalar(il,C_Temp:C_E),0.,&
                fact*this%Cl_scalar(il,C_Cross), scale*this%Cl_scalar(il,C_Phi),&
                (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*this%Cl_scalar(il,C_PhiTemp:C_PhiE)
        end do
        close(unit)
    end if
    end subroutine TCLdata_output_lens_pot_files


    subroutine TCLdata_output_veccl_files(this, CP,VecFile, factor)
    class(TCLData) :: this
    integer :: il, unit
    Type(CAMBParams), intent(in) :: CP
    character(LEN=*), intent(in) :: VecFile
    real(dl), intent(in), optional :: factor
    real(dl) :: fact

    fact = PresentDefault(1._dl, factor)

    if (CP%WantVectors .and. VecFile /= '') then
        unit = open_file_header(VecFile, 'L', CT_name_tags)
        do il=CP%Min_l,CP%Max_l
            write(unit,'(1I6,4E15.6)')il, fact*this%Cl_vector(il, CT_Temp:CT_Cross)
        end do
        close(unit)
    end if

    end subroutine TCLdata_output_veccl_files

    subroutine TCLdata_NormalizeClsAtL(this,CP,lnorm)
    class(TCLData) :: this
    Type(CAMBParams), intent(in) :: CP
    integer, intent(IN) :: lnorm
    real(dl) Norm

    if (CP%WantScalars) then
        Norm=1/this%Cl_scalar(lnorm, C_Temp)
        this%Cl_scalar(CP%Min_l:CP%Max_l, C_Temp:C_Cross) = this%Cl_scalar(CP%Min_l:CP%Max_l, C_Temp:C_Cross) * Norm
    end if

    if (CP%WantTensors) then
        if (.not.CP%WantScalars) Norm = 1/this%Cl_tensor(lnorm, C_Temp)
        !Otherwise Norm already set correctly
        this%Cl_tensor(CP%Min_l:CP%Max_l_tensor, CT_Temp:CT_Cross) =  &
            this%Cl_tensor(CP%Min_l:CP%Max_l_tensor, CT_Temp:CT_Cross) * Norm
    end if

    end subroutine TCLdata_NormalizeClsAtL


    subroutine Thermo_values(this, tau, a, cs2b, opacity, dopacity)
        !Compute unperturbed sound speed squared,
        !and ionization fraction by interpolating pre-computed tables.
        !If requested also get time derivative of opacity
        class(TThermoData), intent(in) :: this
        real(dl), intent(in) :: tau
        real(dl), intent(out) :: a, cs2b
        real(dl), intent(out) :: opacity(nscatter)
        real(dl), intent(out), optional :: dopacity(nscatter)
        integer i
        real(dl) d
    
        d=log(tau/this%tauminn)/this%dlntau+1._dl
        i=int(d)
        d=d-i
        if (i < 1) then
            !Linear interpolation if out of bounds (should not occur).
            write(*,*) 'tau, taumin = ', tau, this%tauminn
            call MpiStop('thermo out of bounds')
        else if (i >= this%nthermo) then
            cs2b=this%cs2(this%nthermo)
            ! and guess
            opacity(1)=this%dotmu(this%nthermo,1)
            a=1
            if (present(dopacity)) then
                dopacity(1) = this%ddotmu(this%nthermo,1)/(tau*this%dlntau)
            end if
            ! and
        ! andrea
        else
            cs2b=this%cs2(i)+d*(this%dcs2(i)+d*(3*(this%cs2(i+1)-this%cs2(i))  &
                -2*this%dcs2(i)-this%dcs2(i+1)+d*(this%dcs2(i)+this%dcs2(i+1)  &
                +2*(this%cs2(i)-this%cs2(i+1)))))
            opacity(1)=this%dotmu(i,1)+d*(this%ddotmu(i,1)+d*(3*(this%dotmu(i+1,1)-this%dotmu(i,1)) &
                -2*this%ddotmu(i,1)-this%ddotmu(i+1,1)+d*(this%ddotmu(i,1)+this%ddotmu(i+1,1) &
                +2*(this%dotmu(i,1)-this%dotmu(i+1,1)))))
            a = (this%ScaleFactor(i)+d*(this%dScaleFactor(i)+d*(3*(this%ScaleFactor(i+1)-this%ScaleFactor(i)) &
                -2*this%dScaleFactor(i)-this%dScaleFactor(i+1)+d*(this%dScaleFactor(i)+this%dScaleFactor(i+1) &
                +2*(this%ScaleFactor(i)-this%ScaleFactor(i+1))))))*tau
            if (present(dopacity)) then
                dopacity(1)=(this%ddotmu(i,1)+d*(this%dddotmu(i,1)+d*(3*(this%ddotmu(i+1,1)  &
                    -this%ddotmu(i,1))-2*this%dddotmu(i,1)-this%dddotmu(i+1,1)+d*(this%dddotmu(i,1) &
                    +this%dddotmu(i+1,1)+2*(this%ddotmu(i,1)-this%ddotmu(i+1,1))))))/(tau*this%dlntau)
        ! 
            end if
        end if
    end subroutine Thermo_values

    subroutine TRecfast_ReadParams(this, Ini)
    use IniObjects
    class(TRecfast) :: this
    class(TIniFile), intent(in) :: Ini

    this%RECFAST_fudge_He = Ini%Read_Double('RECFAST_fudge_He', RECFAST_fudge_He_default)
    this%RECFAST_Heswitch = Ini%Read_Int('RECFAST_Heswitch', RECFAST_Heswitch_default)
    this%RECFAST_Hswitch = Ini%Read_Logical('RECFAST_Hswitch', RECFAST_Hswitch_default)
    this%RECFAST_fudge = Ini%Read_Double('RECFAST_fudge', RECFAST_fudge_default)
    call Ini%Read('AGauss1',this%AGauss1)
    call Ini%Read('AGauss2',this%AGauss2)
    call Ini%Read('zGauss1',this%zGauss1)
    call Ini%Read('zGauss2',this%zGauss2)
    call Ini%Read('wGauss1',this%wGauss1)
    call Ini%Read('wGauss2',this%wGauss2)
    if (this%RECFAST_Hswitch) then
        this%RECFAST_fudge = this%RECFAST_fudge - (RECFAST_fudge_default - RECFAST_fudge_default2)
    end if
    end subroutine TRecfast_ReadParams

    subroutine TRecfast_Validate(this, OK)
    class(TRecfast),intent(in) :: this
    logical, intent(inout) :: OK

    if (this%RECFAST_Heswitch<0 .or. this%RECFAST_Heswitch > 6) then
        OK = .false.
        write(*,*) 'RECFAST_Heswitch unknown'
    end if

    end subroutine TRecfast_Validate

    function TRecfast_tm(this,a)
    class(TRecfast) :: this
    real(dl), intent(in) :: a
    real(dl) zst,z,az,bz,TRecfast_tm
    integer ilo,ihi

    z=1/a-1
    associate( Calc => this%Calc)
        if (z >= Calc%zrec(1)) then
            TRecfast_tm=Calc%Tnow/a
        else
            if (z <=Calc%zrec(nz)) then
                TRecfast_tm=Calc%Tmrec(nz)
            else
                zst=(zinitial-z)/delta_z
                ihi= int(zst)
                ilo = ihi+1
                az=zst - ihi
                bz=1-az
                TRecfast_tm=az*Calc%Tmrec(ilo)+bz*Calc%Tmrec(ihi)+ &
                    ((az**3-az)*Calc%dTmrec(ilo)+(bz**3-bz)*Calc%dTmrec(ihi))/6._dl
            endif
        endif
    end associate

    end function TRecfast_tm

    function TRecfast_ts(this,a)
    class(TRecfast) :: this
    !zrec(1) is zinitial-delta_z
    real(dl), intent(in) :: a
    real(dl) zst,z,az,bz,TRecfast_ts
    integer ilo,ihi

    z=1/a-1
    associate(Calc => this%Calc)
        if (z.ge.Calc%zrec(1)) then
            TRecfast_ts=Calc%tsrec(1)
        else
            if (z.le.Calc%zrec(nz)) then
                TRecfast_ts=Calc%tsrec(nz)
            else
                zst=(zinitial-z)/delta_z
                ihi= int(zst)
                ilo = ihi+1
                az=zst - ihi
                bz=1-az

                TRecfast_ts=az*Calc%tsrec(ilo)+bz*Calc%tsrec(ihi)+ &
                    ((az**3-az)*Calc%dtsrec(ilo)+(bz**3-bz)*Calc%dtsrec(ihi))/6._dl
            endif
        endif
    end associate
    end function TRecfast_ts

    function TRecfast_xe(this,a)
    class(TRecfast) :: this
    real(dl), intent(in) :: a
    real(dl) zst,z,az,bz,TRecfast_xe
    integer ilo,ihi

    z=1/a-1
    associate(Calc => this%Calc)
        if (z.ge.Calc%zrec(1)) then
            TRecfast_xe=Calc%xrec(1)
        else
            if (z.le.Calc%zrec(nz)) then
                TRecfast_xe=Calc%xrec(nz)
            else
                zst=(zinitial-z)/delta_z
                ihi= int(zst)
                ilo = ihi+1
                az=zst - ihi
                bz=1-az
                TRecfast_xe=az*Calc%xrec(ilo)+bz*Calc%xrec(ihi)+ &
                    ((az**3-az)*Calc%dxrec(ilo)+(bz**3-bz)*Calc%dxrec(ihi))/6._dl
            endif
        endif
    end associate
    end function TRecfast_xe

    ! andrea
    ! andrea

    subroutine TRecfast_xe_Tm(this,a, xe, Tm)
    class(TRecfast) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: xe, Tm
    real(dl) zst,z,az,bz
    integer ilo,ihi

    z=1/a-1
    associate(Calc => this%Calc)
        if (z.ge.Calc%zrec(1)) then
            xe=Calc%xrec(1)
            Tm = Calc%Tnow/a
        else
            if (z.le.Calc%zrec(nz)) then
                xe=Calc%xrec(nz)
                TM =  Calc%Tmrec(nz)
            else
                zst=(zinitial-z)/delta_z
                ihi= int(zst)
                ilo = ihi+1
                az=zst - ihi
                bz=1-az
                xe=az*Calc%xrec(ilo)+bz*Calc%xrec(ihi)+ &
                    ((az**3-az)*Calc%dxrec(ilo)+(bz**3-bz)*Calc%dxrec(ihi))/6._dl
                Tm=az*Calc%Tmrec(ilo)+bz*Calc%Tmrec(ihi)+ &
                    ((az**3-az)*Calc%dTmrec(ilo)+(bz**3-bz)*Calc%dTmrec(ihi))/6._dl

            endif
        endif
    end associate
    end subroutine TRecfast_xe_Tm

    function TRecfast_version(this) result(version)
    class(TRecfast) :: this
    character(LEN=:), allocatable :: version

    version = Recfast_Version

    end function TRecfast_version

    subroutine TRecfast_init(this, State, WantTSpin)
    use MiscUtils
    implicit none
    class(TRecfast), target :: this
    class(TCAMBdata), target :: State
    real(dl) :: Trad,Tmat,Tspin
    integer :: I
    Type(RecombinationData), pointer :: Calc
    logical, intent(in), optional :: WantTSpin
    real(dl) :: z,n,x,x0,rhs,x_H,x_He,x_H0,x_He0,H, Yp
    real(dl) :: zstart,zend,z_scale
    real(dl) :: cw(24)
    real(dl), dimension(:,:), allocatable :: w
    real(dl) :: y(4)
    real(dl) :: C10, tau_21Ts
    integer :: ind, nw
    real(dl), parameter :: tol=1.D-5                !Tolerance for R-K
    procedure(TClassDverk) :: dverk


    if (.not. allocated(this%Calc)) allocate(this%Calc)
    Calc => this%Calc

    select type(State)
    class is (CAMBdata)
        Calc%State => State
        Calc%doTspin = DefaultFalse(WantTSpin)


        !       write(*,*)'recfast version 1.0'
        !       write(*,*)'Using Hummer''s case B recombination rates for H'
        !       write(*,*)' with fudge factor = 1.14'
        !       write(*,*)'and tabulated HeII singlet recombination rates'
        !       write(*,*)

        Calc%n_eq = 3
        if (Evolve_Ts) Calc%n_eq=4
        allocate(w(Calc%n_eq,9))

        Calc%Recombination_saha_z=0.d0

        Calc%Tnow = State%CP%tcmb
        !       These are easy to inquire as input, but let's use simple values
        z = zinitial
        !       will output every 1 in z, but this is easily changed also

        H = State%CP%H0/100._dl

        !Not general, but only for approx
        Calc%OmegaT=(State%CP%omch2+State%CP%ombh2)/H**2        !total dark matter + baryons
        Calc%OmegaK=State%CP%omk       !curvature


        !       convert the Hubble constant units
        Calc%HO = H*bigH
        Yp = State%CP%Yhe

        !       sort out the helium abundance parameters
        Calc%mu_H = 1.d0/(1.d0-Yp)           !Mass per H atom
        Calc%mu_T = not4/(not4-(not4-1.d0)*Yp)   !Mass per atom
        Calc%fHe = Yp/(not4*(1.d0-Yp))       !n_He_tot / n_H_tot


        Calc%Nnow = 3._dl*bigH**2*State%CP%ombh2/(const_eightpi*G*Calc%mu_H*m_H)

        n = Calc%Nnow * (1._dl+z)**3
        Calc%z_eq = State%z_eq

        !       Fudge factor to approximate for low z out of equilibrium effect
        Calc%fu=this%RECFAST_fudge

        !       Set initial matter temperature
        y(3) = Calc%Tnow*(1._dl+z)            !Initial rad. & mat. temperature
        Tmat = y(3)
        y(4) = Tmat
        Tspin = Tmat

        call get_init(Calc,z,x_H0,x_He0,x0)

        y(1) = x_H0
        y(2) = x_He0

        !       OK that's the initial conditions, now start writing output file

        !       Set up work-space stuff for DVERK
        ind  = 1
        nw   = Calc%n_eq
        do i = 1,24
            cw(i) = 0._dl
        end do

        do i = 1,Nz
            !       calculate the start and end redshift for the interval at each z
            !       or just at each z
            zstart = zinitial  - real(i-1,dl)*delta_z
            zend   = zinitial  - real(i,dl)*delta_z

            ! Use Saha to get x_e, using the equation for x_e for ionized helium
            ! and for neutral helium.
            ! Everything ionized above z=8000.  First ionization over by z=5000.
            ! Assume He all singly ionized down to z=3500, then use He Saha until
            ! He is 99% singly ionized, and *then* switch to joint H/He recombination.

            z = zend
            z_scale = Calc%Tnow/COBE_CMBTemp * (1+z) -1

            if (z_scale > 8000._dl) then

                x_H0 = 1._dl
                x_He0 = 1._dl
                x0 = 1._dl+2._dl*Calc%fHe
                y(1) = x_H0
                y(2) = x_He0
                y(3) = Calc%Tnow*(1._dl+z)
                y(4) = y(3)

            else if(z_scale > 5000._dl)then

                x_H0 = 1._dl
                x_He0 = 1._dl
                rhs = exp( 1.5d0 * log(CR*Calc%Tnow/(1._dl+z)) &
                    - CB1_He2/(Calc%Tnow*(1._dl+z)) ) / Calc%Nnow
                rhs = rhs*1._dl            !ratio of g's is 1 for He++ <-> He+
                x0 = 0.5d0 * ( sqrt( (rhs-1._dl-Calc%fHe)**2 &
                    + 4._dl*(1._dl+2._dl*Calc%fHe)*rhs) - (rhs-1._dl-Calc%fHe) )
                y(1) = x_H0
                y(2) = x_He0
                y(3) = Calc%Tnow*(1._dl+z)
                y(4) = y(3)

            else if(z_scale > 3500._dl)then

                x_H0 = 1._dl
                x_He0 = 1._dl
                x0 = x_H0 + Calc%fHe*x_He0
                y(1) = x_H0
                y(2) = x_He0
                y(3) = Calc%Tnow*(1._dl+z)
                y(4) = y(3)

            else if(y(2) > 0.99)then

                x_H0 = 1._dl
                rhs = exp( 1.5d0 * log(CR*Calc%Tnow/(1._dl+z)) &
                    - CB1_He1/(Calc%Tnow*(1._dl+z)) ) / Calc%Nnow
                rhs = rhs*4._dl            !ratio of g's is 4 for He+ <-> He0
                x_He0 = 0.5d0 * ( sqrt( (rhs-1._dl)**2 &
                    + 4._dl*(1._dl+Calc%fHe)*rhs )- (rhs-1._dl))
                x0 = x_He0
                x_He0 = (x0 - 1._dl)/Calc%fHe
                y(1) = x_H0
                y(2) = x_He0
                y(3) = Calc%Tnow*(1._dl+z)
                y(4) = y(3)

            else if (y(1) > 0.99d0) then

                rhs = exp( 1.5d0 * log(CR*Calc%Tnow/(1._dl+z)) &
                    - CB1/(Calc%Tnow*(1._dl+z)) ) / Calc%Nnow
                x_H0 = 0.5d0 * (sqrt( rhs**2+4._dl*rhs ) - rhs )

                call DVERK(this,3,ION,zstart,y,zend,tol,ind,cw,nw,w)
                y(1) = x_H0
                x0 = y(1) + Calc%fHe*y(2)
                y(4)=y(3)
            else

                call DVERK(this,nw,ION,zstart,y,zend,tol,ind,cw,nw,w)

                x0 = y(1) + Calc%fHe*y(2)

            end if

            Trad = Calc%Tnow * (1._dl+zend)
            Tmat = y(3)
            x_H = y(1)
            x_He = y(2)
            x = x0

            Calc%zrec(i)=zend
            Calc%xrec(i)=x
            ! andrea
            Calc%x_rayleigh_eff(i)=  1-x_H + HeRayleighFac*(1 - x_He)*Calc%fHe
            ! andrea
            Calc%tmrec(i) = Tmat


            if (Calc%doTspin) then
                if (Evolve_Ts .and. zend< 1/Do21cm_minev-1 ) then
                    Tspin = y(4)
                else
                    C10 = Calc%Nnow * (1._dl+zend)**3*(kappa_HH_21cm(Tmat,.false.)*(1-x_H) &
                        + kappa_eH_21cm(Tmat,.false.)*x)
                    tau_21Ts = line21_const*Calc%NNow*(1+zend)*dtauda(State,1/(1+zend))/1000

                    Tspin = Trad*( C10/Trad + A10/T_21cm)/(C10/Tmat + A10/T_21cm) + &
                        tau_21Ts/2*A10*( 1/(C10*T_21cm/Tmat+A10) -  1/(C10*T_21cm/Trad+A10) )

                    y(4) = Tspin
                end if

                Calc%tsrec(i) = Tspin

            end if

            !          write (*,'(5E15.5)') zend, Trad, Tmat, Tspin, x
        end do

        ! and
        call spline_def(Calc%zrec,Calc%x_rayleigh_eff,nz,Calc%dx_rayleigh_eff)
        call spline_def(Calc%zrec,Calc%tmrec,nz,Calc%dtmrec)
        ! and
        if (Calc%doTspin) then
            call spline_def(Calc%zrec,Calc%tsrec,nz,Calc%dtsrec)
        end if
    class default
        call MpiStop('Wrong state type')
    end select

    end subroutine TRecfast_init

    !       ===============================================================
    subroutine GET_INIT(Calc,z,x_H0,x_He0,x0)

    !       Set up the initial conditions so it will work for general,
    !       but not pathological choices of zstart
    !       Initial ionization fraction using Saha for relevant species
    Type(RecombinationData) :: Calc
    real(dl) z,x0,rhs,x_H0,x_He0, z_scale

    z_scale = Calc%Tnow/COBE_CMBTemp*(z+1)-1

    if(z_scale > 8000._dl)then
        x_H0 = 1._dl
        x_He0 = 1._dl
        x0 = 1._dl+2._dl*Calc%fHe

    else if(z_scale > 3500._dl)then

        x_H0 = 1._dl
        x_He0 = 1._dl
        rhs = exp( 1.5d0 * log(CR*Calc%Tnow/(1._dl+z)) &
            - CB1_He2/(Calc%Tnow*(1._dl+z)) ) / Calc%Nnow
        rhs = rhs*1._dl    !ratio of g's is 1 for He++ <-> He+
        x0 = 0.5d0 * ( sqrt( (rhs-1._dl-Calc%fHe)**2 &
            + 4._dl*(1._dl+2._dl*Calc%fHe)*rhs) - (rhs-1._dl-Calc%fHe) )

    else if(z_scale > 2000._dl)then

        x_H0 = 1._dl
        rhs = exp( 1.5d0 * log(CR*Calc%Tnow/(1._dl+z)) &
            - CB1_He1/(Calc%Tnow*(1._dl+z)) ) / Calc%Nnow
        rhs = rhs*4._dl    !ratio of g's is 4 for He+ <-> He0
        x_He0 = 0.5d0  * ( sqrt( (rhs-1._dl)**2 + 4._dl*(1._dl+Calc%fHe)*rhs )- (rhs-1._dl))
        x0 = x_He0
        x_He0 = (x0 - 1._dl)/Calc%fHe

    else

        rhs = exp( 1.5d0 * log(CR*Calc%Tnow/(1._dl+z)) &
            - CB1/(Calc%Tnow*(1._dl+z)) ) / Calc%Nnow
        x_H0 = 0.5d0 * (sqrt( rhs**2+4._dl*rhs ) - rhs )
        x_He0 = 0._dl
        x0 = x_H0
    end if

    end subroutine GET_INIT

    subroutine ION(this,Ndim,z,Y,f)
    class(TRecfast), target :: this
    integer Ndim

    real(dl) z,x,n,n_He,Trad,Tmat,Tspin,x_H,x_He, Hz
    real(dl) y(Ndim),f(Ndim)
    real(dl) Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz
    real(dl) timeTh,timeH
    real(dl) a_VF,b_VF,T_0,T_1,sq_0,sq_1,a_PPB,b_PPB,c_PPB,d_PPB
    real(dl) tauHe_s,pHe_s
    real(dl) a_trip,b_trip,Rdown_trip,Rup_trip
    real(dl) Doppler,gamma_2Ps,pb,qb,AHcon
    real(dl) tauHe_t,pHe_t,CL_PSt,CfHe_t,gamma_2Pt
    real(dl) epsilon
    integer Heflag
    real(dl) C10, dHdz, z_scale
    type(RecombinationData), pointer :: Recomb

    Recomb => this%Calc

    !       the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen
    a_PPB = 4.309d0
    b_PPB = -0.6166d0
    c_PPB = 0.6703d0
    d_PPB = 0.5300d0
    !       the Verner and Ferland type fitting parameters for Helium
    !       fixed to match those in the SSS papers, and now correct
    a_VF = 10.d0**(-16.744d0)
    b_VF = 0.711d0
    T_0 = 10.d0**(0.477121d0)   !3K
    T_1 = 10.d0**(5.114d0)
    !      fitting parameters for HeI triplets
    !      (matches Hummer's table with <1% error for 10^2.8 < T/K < 10^4)

    a_trip = 10.d0**(-16.306d0)
    b_trip = 0.761D0


    x_H = y(1)
    x_He = y(2)
    x = x_H + Recomb%fHe * x_He
    Tmat = y(3)
    !        Tspin = y(4)

    n = Recomb%Nnow * (1._dl+z)**3
    n_He = Recomb%fHe * Recomb%Nnow * (1._dl+z)**3
    Trad = Recomb%Tnow * (1._dl+z)

    Hz = 1/dtauda(Recomb%State,1/(1._dl+z))*(1._dl+z)**2/MPC_in_sec


    !       Get the radiative rates using PPQ fit, identical to Hummer's table

    Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB &
        /(1._dl+c_PPB*(Tmat/1.d4)**d_PPB)
    Rup = Rdown * (CR*Tmat)**(1.5d0)*exp(-CDB/Tmat)

    !       calculate He using a fit to a Verner & Ferland type formula
    sq_0 = sqrt(Tmat/T_0)
    sq_1 = sqrt(Tmat/T_1)
    !       typo here corrected by Wayne Hu and Savita Gahlaut
    Rdown_He = a_VF/(sq_0*(1.d0+sq_0)**(1.d0-b_VF))
    Rdown_He = Rdown_He/(1.d0+sq_1)**(1.d0+b_VF)
    Rup_He = Rdown_He*(CR*Tmat)**(1.5d0)*exp(-CDB_He/Tmat)
    Rup_He = 4.d0*Rup_He    !statistical weights factor for HeI
    !       Avoid overflow (pointed out by Jacques Roland)
    if((Bfact/Tmat) > 680.d0)then
        He_Boltz = exp(680.d0)
    else
        He_Boltz = exp(Bfact/Tmat)
    end if
    !   now deal with H and its fudges
    if (.not. this%RECFAST_Hswitch) then
        K = CK/Hz !Peebles coefficient K=lambda_a^3/8piH
    else
        !c  fit a double Gaussian correction function
        z_scale = this%Calc%Tnow/COBE_CMBTemp*(1+z)-1
        K = CK/Hz*(1.0d0 &
            +this%AGauss1*exp(-((log(1.0d0+z_scale)-this%zGauss1)/this%wGauss1)**2.d0) &
            +this%AGauss2*exp(-((log(1.0d0+z_scale)-this%zGauss2)/this%wGauss2)**2.d0))
    end if


    !  add the HeI part, using same T_0 and T_1 values
    Rdown_trip = a_trip/(sq_0*(1.d0+sq_0)**(1.0-b_trip))
    Rdown_trip = Rdown_trip/((1.d0+sq_1)**(1.d0+b_trip))
    Rup_trip = Rdown_trip*dexp(-h_P*C*L_He2St_ion/(k_B*Tmat))
    Rup_trip = Rup_trip*((CR*Tmat)**(1.5d0))*(4.d0/3.d0)
    !   last factor here is the statistical weight

    !       try to avoid "NaN" when x_He gets too small
    if ((x_He.lt.5.d-9) .or. (x_He.gt.0.98d0)) then
        Heflag = 0
    else
        Heflag = this%RECFAST_Heswitch
    end if
    if (Heflag.eq.0)then        !use Peebles coeff. for He
        K_He = CK_He/Hz
    else    !for Heflag>0       !use Sobolev escape probability
        tauHe_s = A2P_s*CK_He*3.d0*n_He*(1.d0-x_He)/Hz
        pHe_s = (1.d0 - dexp(-tauHe_s))/tauHe_s
        K_He = 1.d0/(A2P_s*pHe_s*3.d0*n_He*(1.d0-x_He))
        !      if (((Heflag.eq.2) .or. (Heflag.ge.5)) .and. x_H < 0.99999d0) then
        if (((Heflag.eq.2) .or. (Heflag.ge.5)) .and. x_H < 0.9999999d0) then
            !AL changed July 08 to get smoother Helium

            !   use fitting formula for continuum opacity of H
            !   first get the Doppler width parameter
            Doppler = 2.D0*k_B*Tmat/(m_H*not4*C*C)
            Doppler = C*L_He_2p*dsqrt(Doppler)
            gamma_2Ps = 3.d0*A2P_s*Recomb%fHe*(1.d0-x_He)*C*C &
                /(dsqrt(const_pi)*sigma_He_2Ps*const_eightpi*Doppler*(1.d0-x_H)) &
                /((C*L_He_2p)**2.d0)
            pb = 0.36d0  !value from KIV (2007)
            qb = this%RECFAST_fudge_He
            !   calculate AHcon, the value of A*p_(con,H) for H continuum opacity
            AHcon = A2P_s/(1.d0+pb*(gamma_2Ps**qb))
            K_He=1.d0/((A2P_s*pHe_s+AHcon)*3.d0*n_He*(1.d0-x_He))
        end if
        if (Heflag.ge.3) then     !include triplet effects
            tauHe_t = A2P_t*n_He*(1.d0-x_He)*3.d0
            tauHe_t = tauHe_t /(const_eightpi*Hz*L_He_2Pt**(3.d0))
            pHe_t = (1.d0 - dexp(-tauHe_t))/tauHe_t
            CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
            if ((Heflag.eq.3) .or. (Heflag.eq.5).or.(x_H.gt.0.99999d0)) then !Recfast 1.4.2 (?)
                !        if ((Heflag.eq.3) .or. (Heflag.eq.5) .or. x_H >= 0.9999999d0) then    !no H cont. effect
                CfHe_t = A2P_t*pHe_t*dexp(-CL_PSt/Tmat)
                CfHe_t = CfHe_t/(Rup_trip+CfHe_t)   !"C" factor for triplets
            else                  !include H cont. effect
                Doppler = 2.d0*k_B*Tmat/(m_H*not4*C*C)
                Doppler = C*L_He_2Pt*dsqrt(Doppler)
                gamma_2Pt = 3.d0*A2P_t*Recomb%fHe*(1.d0-x_He)*C*C &
                    /(dsqrt(const_pi)*sigma_He_2Pt*const_eightpi*Doppler*(1.d0-x_H)) &
                    /((C*L_He_2Pt)**2.d0)
                !   use the fitting parameters from KIV (2007) in this case
                pb = 0.66d0
                qb = 0.9d0
                AHcon = A2P_t/(1.d0+pb*gamma_2Pt**qb)/3.d0
                CfHe_t = (A2P_t*pHe_t+AHcon)*dexp(-CL_PSt/Tmat)
                CfHe_t = CfHe_t/(Rup_trip+CfHe_t)   !"C" factor for triplets
            end if
        end if
    end if


    !       Estimates of Thomson scattering time and Hubble time
    timeTh=(1._dl/(CT*Trad**4))*(1._dl+x+Recomb%fHe)/x       !Thomson time
    timeH=2./(3.*Recomb%HO*(1._dl+z)**1.5)      !Hubble time

    !       calculate the derivatives
    !       turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
    !       (clunky, but seems to work)
    if (x_H > 0.99) then   !don't change at all
        f(1) = 0._dl
        !!        else if (x_H > 0.98_dl) then
    else if (x_H.gt.0.985d0) then     !use Saha rate for Hydrogen
        f(1) = (x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-CL/Tmat)) /(Hz*(1.d0+z))
        Recomb%Recombination_saha_z = z
        !AL: following commented as not used
        !   for interest, calculate the correction factor compared to Saha
        !   (without the fudge)
        !       factor=(1.d0 + K*Lambda*n*(1.d0-x_H))
        !       /(Hz*(1.d0+z)*(1.d0+K*Lambda*n*(1.d0-x)
        !       +K*Rup*n*(1.d0-x)))
    else !use full rate for H

        f(1) = ((x*x_H*n*Rdown - Rup*(1.d0-x_H)*exp(-CL/Tmat)) &
            *(1.d0 + K*Lambda*n*(1.d0-x_H))) &
            /(Hz*(1.d0+z)*(1.d0/Recomb%fu+K*Lambda*n*(1.d0-x_H)/Recomb%fu &
            +K*Rup*n*(1.d0-x_H)))

    end if

    !       turn off the He once it is small
    if (x_He < 1.e-15) then
        f(2)=0.d0
    else

        f(2) = ((x*x_He*n*Rdown_He &
            - Rup_He*(1-x_He)*exp(-CL_He/Tmat)) &
            *(1 + K_He*Lambda_He*n_He*(1.d0-x_He)*He_Boltz)) &
            /(Hz*(1+z) &
            * (1 + K_He*(Lambda_He+Rup_He)*n_He*(1.d0-x_He)*He_Boltz))

        !   Modification to HeI recombination including channel via triplets
        if (Heflag.ge.3) then
            f(2) = f(2)+ (x*x_He*n*Rdown_trip &
                - (1.d0-x_He)*3.d0*Rup_trip*dexp(-h_P*C*L_He_2st/(k_B*Tmat))) &
                *CfHe_t/(Hz*(1.d0+z))
        end if

    end if

    if (timeTh < H_frac*timeH) then
        !                f(3)=Tmat/(1._dl+z)      !Tmat follows Trad
        !   additional term to smooth transition to Tmat evolution,
        !   (suggested by Adam Moss)
        dHdz = (Recomb%HO**2/2.d0/Hz)*(4.d0*(1.d0+z)**3/(1.d0+Recomb%z_eq)*Recomb%OmegaT &
            + 3.d0*Recomb%OmegaT*(1.d0+z)**2 + 2.d0*Recomb%OmegaK*(1.d0+z) )

        epsilon = Hz*(1.d0+x+Recomb%fHe)/(CT*Trad**3*x)
        f(3) = Recomb%Tnow &
            + epsilon*((1.d0+Recomb%fHe)/(1.d0+Recomb%fHe+x))*((f(1)+Recomb%fHe*f(2))/x) &
            - epsilon* dHdz/Hz + 3.0d0*epsilon/(1.d0+z)

    else
        f(3)= CT * (Trad**4) * x / (1._dl+x+Recomb%fHe) &
            * (Tmat-Trad) / (Hz*(1._dl+z)) + 2._dl*Tmat/(1._dl+z)
    end if

    if (evolve_Ts) then

        !       follow the matter temperature once it has a chance of diverging
        if (timeTh < H_frac*timeH) then
            f(4) = Recomb%Tnow !spin follows Trad and Tmat
        else
            if (z< 1/Do21cm_minev-1) then

                Tspin = y(4)
                C10 = n*(kappa_HH_21cm(Tmat,.false.)*(1-x_H) + kappa_eH_21cm(Tmat,.false.)*x)

                f(4) = 4*Tspin/Hz/(1+z)*( (Tspin/Tmat-1._dl)*C10 + Trad/T_21cm*(Tspin/Trad-1._dl)*A10) - f(1)*Tspin/(1-x_H)
            else
                f(4)=f(3)
            end if
        end if

    end if

    end subroutine ION


    function TRecfast_dDeltaxe_dtau(this,a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb,adotoa)
    !d x_e/d tau assuming Helium all neutral and temperature perturbations negligible
    !it is not accurate for x_e of order 1
    class(TRecfast) :: this
    real(dl) TRecfast_dDeltaxe_dtau
    real(dl), intent(in):: a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb,adotoa
    real(dl) Delta_Tg
    real(dl) xedot,z,x,n,n_He,Trad,Tmat,x_H,Hz, C_r, dlnC_r
    real(dl) Rup,Rdown,K
    real(dl) a_PPB,b_PPB,c_PPB,d_PPB
    real(dl) delta_alpha, delta_beta, delta_K, clh
    real(dl) xe

    associate(Calc=>this%Calc)
        Delta_tg =Delta_Tm
        call this%xe_Tm(a, xe, Tmat)
        x_H = min(1._dl,xe)

        !       the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen
        a_PPB = 4.309d0
        b_PPB = -0.6166d0
        c_PPB = 0.6703d0
        d_PPB = 0.5300d0

        z=1/a-1

        x = x_H

        n = Calc%Nnow /a**3
        n_He = Calc%fHe * n
        Trad = Calc%Tnow /a
        clh = adotoa !conformal time Hubble
        Hz = clh/a/MPC_in_sec !normal time in seconds

        !       Get the radiative rates using PPQ fit, identical to Hummer's table

        Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB &
            /(1._dl+c_PPB*(Tmat/1.d4)**d_PPB)   !alpha
        Rup = Rdown * (CR*Tmat)**(1.5d0)*exp(-CDB/Tmat)

        K = CK/Hz              !Peebles coefficient K=lambda_a^3/8piH


        Rdown = Rdown*Calc%fu
        Rup = Rup*Calc%fu
        C_r =  a*(1.d0 + K*Lambda*n*(1.d0-x_H)) /( 1.d0+K*(Lambda+Rup)*n*(1.d0-x_H) )*MPC_in_sec

        xedot = -(x*x_H*n*Rdown - Rup*(1.d0-x_H)*exp(-CL/Tmat))*C_r

        delta_alpha = (b_PPB + c_PPB*(Tmat/1d4)**d_PPB*(b_PPB-d_PPB))/(1+c_PPB*(Tmat/1d4)**d_PPB)*Delta_Tg
        delta_beta = delta_alpha + (3./2 + CDB/Tmat)*delta_Tg !(Rup = beta)
        delta_K = - hdot/clh - kvb/clh/3


        dlnC_r = -Rup*K*n*( (Delta_nH+Delta_K + Delta_beta*(1+K*Lambda*n*(1-x_H)))*(1-x_H) - x_H*Delta_xe) &
            / ( 1.d0+K*(Lambda+Rup)*n*(1.d0-x_H) ) /(1.d0 + K*Lambda*n*(1.d0-x_H))

        TRecfast_dDeltaxe_dtau= xedot/x_H*(dlnC_r +Delta_alpha - Delta_xe) &
            - C_r*( (2*Delta_xe + Delta_nH)*x_H*n*Rdown + (Delta_xe - (3./2+ CB1/Tmat)*(1/x_H-1)*Delta_Tg)*Rup*exp(-CL/Tmat))
    end associate

    !Approximate form valid at late times
    !        dDeltaxe_dtau= xedot/x_H*(Delta_alpha + Delta_xe + Delta_nH)


    end function TRecfast_dDeltaxe_dtau

    real(dl) function TRecfast_Get_Saha_z(this)
    class(TRecfast) :: this
    TRecfast_Get_Saha_z =  this%Calc%recombination_saha_z
    end function


    subroutine TRecfast_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TRecfast), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TRecfast_SelfPointer

    end module Recombination
