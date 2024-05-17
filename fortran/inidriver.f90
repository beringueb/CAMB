    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
    !     This is a sample driver routine that reads
    !     in one set of parameters and produdes the corresponding output.

    program driver
    use CAMB
    implicit none
    character(len=:), allocatable :: InputFile
    
    ! andrea
    character(LEN=Ini_max_string_len) RayleighTerms
    ! andrea

    InputFile = ''
    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') error stop 'No parameter input file'
    ! andrea
    RayleighTerms = Ini_read_String('rayleigh_pows')
    rayleigh_diff = Ini_read_logical('rayleigh_diff',rayleigh_diff)
    if (RayleighTerms/='') read(RayleighTerms,*) rayleigh_pows
    ! andrea
    call CAMB_CommandLineRun(InputFile)
    deallocate(InputFile) ! Just so no memory leaks in valgrind

    end program driver

