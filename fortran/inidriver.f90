    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
    !     This is a sample driver routine that reads
    !     in one set of parameters and produdes the corresponding output.

    program driver
    use CAMB
    use Recombination
    implicit none
    class(TThermoData), allocatable :: this
    class(TRecfast), allocatable :: etat
    character(len=:), allocatable :: InputFile

    allocate(this)
    allocate(etat)

    InputFile = ''
    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') error stop 'No parameter input file'

    call CAMB_CommandLineRun(this, etat, InputFile)
    deallocate(InputFile) ! Just so no memory leaks in valgrind

    deallocate(this)
    deallocate(etat)

    end program driver

