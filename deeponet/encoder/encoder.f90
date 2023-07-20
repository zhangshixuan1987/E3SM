program encoder
    use torch_ftn
    use iso_fortran_env
    implicit none
    integer :: n
    type(torch_module) :: torch_mod
    type(torch_tensor_wrap) :: input_tensors
    type(torch_tensor) :: out_tensor
     integer :: arglen, stat
    integer ::      NX, NY, NXY
    real(real32), allocatable :: data_in(:)
    real(real32), allocatable :: input(:, :)
    real(real32), allocatable ::  input_f(:, :, :, :)
    character(:), allocatable :: filename
    real(real32), pointer :: output(:, :)

   
    NX = 70 ! No. of Longitude
    NY = 70 ! No of Latitude 
    NXY = NX*NY

    allocate (data_in(NXY))
    allocate (input(NX, NY))
    allocate (input_f(NX, NY, 1, 1)) ! (FD(NX) , FD(NY), NCCHANNEL, BATCHSIZE)
   
    if (command_argument_count() /= 1) then
        print *, "Need to pass a single argument: Pytorch model file name"
        stop
    end if

    call get_command_argument(number=1, length=arglen)
    allocate(character(arglen) :: filename)
    call get_command_argument(number=1, value=filename, status=stat)


    open(12, file="U_bf_revised_checked.out")
    read(12, *) data_in
    close(12)

    print *, "Shape is: ", shape(data_in)

    input = reshape(data_in, (/NX, NY/))

    input_f = reshape(data_in, (/NX, NY, 1, 1/))
    
    call input_tensors%create
    call input_tensors%add_array(input_f)
    call torch_mod%load(filename)
    call torch_mod%forward(input_tensors, out_tensor)
    call out_tensor%to_array(output)
    
    open(unit=13, file = "enc.dat", action="write", status="new")
    write(13, *) output
    close(13)
end program encoder
