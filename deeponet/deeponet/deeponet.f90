program deeponet
    use torch_ftn
    use iso_fortran_env

    implicit none

    integer :: n
    type(torch_module) :: torch_mod
    type(torch_tensor_wrap) :: input_tensors
    type(torch_tensor) :: out_tensor
    integer :: arglen, stat

    integer :: NX, NY, NT, NXY, NXYT
    
    real(real32), allocatable :: data_trunk(:)
    
    real(real32), allocatable :: data_t(:, :)

    real(real32), allocatable :: input_data(:) ! 248 x 36
    
    real(real32), allocatable :: input_f(:, :, :, :)

    real(real32), pointer :: output(:, :, :, :)

    character(:), allocatable :: filename
    
    NX = 6
    NY = 6
    NT = 248
    NXY  = NX*NY
    NXYT = NX * NY * NT 

    allocate (data_trunk(NT))
    allocate (data_t(1, NT))
    allocate (input_data(NXYT))
    allocate (input_f(NXY, NT, 1, 1))

    if (command_argument_count() /= 1) then
        print *, "Need to pass a single argument: Pytorch model file name"
        stop
    end if

    call get_command_argument(number=1, length=arglen)
    allocate(character(arglen) :: filename)
    call get_command_argument(number=1, value=filename, status=stat)


    open(51, file="x_trunk.out")
    read(51, *) data_trunk
    close(51)


    open(1001, file="U_branch_rev2.out")
    read(1001, *) input_data
    close(1001)


    input_f = reshape(input_data, (/NXY, NT, 1, 1/))

    data_t = reshape(data_trunk, (/1, NT/))

    call input_tensors%create
    call input_tensors%add_array(input_f)
    call input_tensors%add_array(data_t)
    call torch_mod%load(filename)
    call torch_mod%forward(input_tensors, out_tensor, 1)
    call out_tensor%to_array(output)
    print *, "Output is", shape(output)
    open(18, file = "dnet.dat", action="write", status = "replace")
    write(18, *) output
    close(18)
end program deeponet
