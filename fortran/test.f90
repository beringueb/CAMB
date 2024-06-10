 subroutine TRecfast_init(this, State, WantTSpin)
    use MiscUtils
    implicit none
    class(TRecfast), target :: this
    class(TCAMBdata), target :: State

subroutine Thermo_Init(this, etat, State, taumin)
    !  Compute and save unperturbed baryon temperature and ionization fraction
    !  as a function of time.  With nthermo=10000, xe(tau) has a relative
    ! accuracy (numerical integration precision) better than 1.e-5.
    use constants
    use StringUtils
    class(TThermoData) :: this
    class(TRecfast) :: etat
    class(CAMBdata), target :: State

call etat%Init(State,WantTSpin=CP%Do21cm)
call etat%Init(State,WantTSpin=CP%Do21cm)

subroutine InitVars(this,state,etat)
    class(TTHermodata) :: this
    class(TRecfast) :: etat
    type(CAMBdata) :: state
