module m_tests
    use, intrinsic :: iso_c_binding
    use m_fftw3
    implicit none
    
    contains
    subroutine test_1d_guru()
        integer, parameter :: NX = 8, NY = 16, NZ = 4
        real(c_double), parameter :: eps = epsilon(1.0_C_DOUBLE)
        real(c_double), dimension(NX) :: in_x, out_x
        real(c_double), dimension(NY) :: in_y, out_y
        type(c_ptr) :: plan_x_1d, plan_y_1d
        integer :: iunit, i, j, k
        
        real(c_double), pointer, dimension(:, :, :) :: out
        type(c_ptr) :: out_ptr
        
        ! guru information
        integer(C_INT) :: rank, howmany_rank
        type(fftw_iodim), allocatable, dimension(:) :: dims, howmany_dims
        integer(C_FFTW_R2R_KIND), allocatable, dimension(:) :: kind_x, kind_y
        type(c_ptr) :: plan_x, plan_y
        
        ! check equivalence
        logical :: is_same
        
        
        print '(A)', "test_1d_guru"
        
        ! Perform 1D transform to the 1D data (baseline)
        plan_x_1d = fftw_plan_r2r_1d(NX, in_x, out_x, FFTW_R2HC, FFTW_ESTIMATE)
        plan_y_1d = fftw_plan_r2r_1d(NY, in_y, out_y, FFTW_R2HC, FFTW_ESTIMATE)
        
        ! I put some number in out_x and out_y
        in_x = [1.0_C_DOUBLE, 2.3_C_DOUBLE, 1.4_C_DOUBLE, 4.0_C_DOUBLE,  &
        1.32_C_DOUBLE, 3.0_C_DOUBLE, 1.0_C_DOUBLE,  4.2_C_DOUBLE]
        in_y = [1.0_C_DOUBLE, 2.3_C_DOUBLE, 1.4_C_DOUBLE, 4.0_C_DOUBLE,  &
        1.32_C_DOUBLE, 3.0_C_DOUBLE, 1.0_C_DOUBLE,  4.2_C_DOUBLE, &
        2.45_C_DOUBLE, 0.32_C_DOUBLE, 3.4_C_DOUBLE, 0.22_C_DOUBLE, &
        1.45_C_DOUBLE, 0.98_C_DOUBLE, 2.23_C_DOUBLE, 1.02_C_DOUBLE]
        
        call fftw_execute_r2r(plan_x_1d, in_x, out_x)
        call fftw_execute_r2r(plan_y_1d, in_y, out_y)
        
        call fftw_destroy_plan(plan_x_1d)
        call fftw_destroy_plan(plan_y_1d)
        
        ! memory for out
        out_ptr = fftw_alloc_real(int(NX*NY*NZ, c_size_t))
        call c_f_pointer(out_ptr, out, [NX, NY, NZ])
        
        rank = 1
        howmany_rank = 2
        allocate(dims(rank), howmany_dims(howmany_rank))
        allocate(kind_x(rank), kind_y(rank))
        
        ! We first do the 1D transform along x
        ! The logical is as
        ! do k = 1, Nz
        !   do j = 1, Ny
        !     work tranform on u(:, j, k)
        !   end do
        ! end do
        dims(1)%n = NX
        dims(1)%is = 1 ! element of u(:, j, k) is continuous, a.k.a, stride = 1
        dims(1)%os = 1
        
        ! do j = 1, Ny
        howmany_dims(1)%n  = NY
        ! increment j by 1, the stride is NX
        howmany_dims(1)%is = NX 
        howmany_dims(1)%os = NX
        
        ! do k = 1, Nz
        howmany_dims(2)%n = NZ
        ! increment k by 1, the stride is NX*NY
        howmany_dims(2)%is = NX*NY
        howmany_dims(2)%os = NX*NY
        
        kind_x(1) = FFTW_R2HC
        plan_x = fftw_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, out, out, kind_x, FFTW_ESTIMATE)
        do k = 1, NZ
            do j = 1, NY
                out(:, j, k) = in_x(:)
            end do
        end do
        call fftw_execute_r2r(plan_x, out, out)
        open(newunit=iunit, file="test_1d_guru_x.txt")
        write(iunit, '(8E20.12)') out_x
        is_same = .true.
        do k = 1, NZ
            do j = 1, NY
                write(iunit, '(8E20.12)') out(:, j, k)
                if (norm2(out_x - out(:, j, k)) > eps) then
                    is_same = .false.
                    exit
                end if
            end do
        end do
        if (is_same) then
            print *, "1D in x-direction: PASS"
        else
            print *, "1D in x-direction: FAIL"
        end if
        close(iunit)
        call fftw_destroy_plan(plan_x)
        
        ! Now, we do the 1D transform along y
        ! The logical is as
        ! do k = 1, Nz
        !   do i = 1, Nx
        !     work transform on u(i, :, k)
        !   end do
        ! end do
        ! for u(i, :, k)
        ! u(i, j, k) and u(i, j + 1, k) is separated by NX
        dims(1)%n  = NY
        dims(1)%is = NX
        dims(1)%os = NX
        
        ! do i = 1, Nx
        howmany_dims(1)%n = NX
        ! increment i by 1, the stride is 1
        howmany_dims(1)%is = 1
        howmany_dims(1)%os = 1
        
        ! do k = 1, Nz
        howmany_dims(2)%n = NZ
        ! increment k by 1, the stride is NX*NY
        howmany_dims(2)%is = NX*NY
        howmany_dims(2)%os = NX*NY
        
        kind_y(1) = FFTW_R2HC
        plan_y = fftw_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, out, out, kind_y, FFTW_ESTIMATE)
        do k = 1, NZ
            do i = 1, NX
                out(i, :, k) = in_y(:)
            end do
        end do
        call fftw_execute_r2r(plan_y, out, out)
        open(newunit=iunit, file="test_1d_guru_y.txt")
        write(iunit, '(16E20.12)') out_y
        is_same = .true.
        do k = 1, NZ
            do i = 1, NX
                write(iunit, '(8E20.12)') out(i, :, k)
                if (norm2(out_y - out(i, :, k)) > eps) then
                    is_same = .false.
                    exit
                end if
            end do
        end do
        if (is_same) then
            print *, "1D in y-direction: PASS"
        else
            print *, "1D in y-direction: FAIL"
        end if
        close(iunit)
        call fftw_destroy_plan(plan_y)
        call fftw_free(out_ptr)
        deallocate(dims, howmany_dims)
        deallocate(kind_x, kind_y)
    end subroutine test_1d_guru
    
    subroutine test_1d_r2r()
        integer, parameter :: N = 8
        type(c_ptr)    :: plan
        real(c_double), dimension(N) :: in
        real(c_double), dimension(N) :: out
        integer :: i, iunit
        
        print '(A)', "test_1d_r2r"
        plan = fftw_plan_r2r_1d(N, in, out, FFTW_R2HC, FFTW_ESTIMATE)
        
        in = [1.0_c_double,  &
        2.3_c_double,  &
        1.4_c_double,  &
        4.0_c_double,  &
        1.32_c_double, &
        3.0_c_double,  &
        1.0_c_double,  &
        4.2_c_double]
        call fftw_execute_r2r(plan, in, out)
        
        
        
        ! Print the output
        
        open(newunit=iunit, file="test_1d_r2r.txt")
        write(iunit, '(A)')  "# Output of FFT:"
        write(iunit, '(A)')  "#   i         out(i)"
        do i = 1, N
            write(iunit, '(i5, es23.15)') i - 1, out(i)
        end do
        close(iunit)
        
        call fftw_destroy_plan(plan)
    end subroutine test_1d_r2r
    
    ! Solve 3D Poisson equation with Pseudo Spectral Method
    subroutine test_3d_r2c()
        real(c_double), parameter :: PI = 4.0_c_double*atan(1.0_c_double)
        real(c_double), parameter :: PI2 = 2.0_c_double*PI
        integer, parameter :: N = 128
        real(c_double), parameter :: dx = PI2/real(N, c_double)
        real(c_double), parameter :: factor = 1.0_c_double/real(N*N*N, c_double)
        type(c_ptr) :: fwd_plan, bwd_plan
        integer :: iunit, i, j, k
        
        ! reference data
        real(c_double), allocatable, dimension(:) :: x, y, z
        real(c_double), allocatable, dimension(:, :, :) :: p_ref
        real(c_double), allocatable, dimension(:) :: akx, aky, akz
        real(c_double), allocatable, dimension(:, :, :) :: laplace
        
        

        ! fftw stuff
        real(c_double), pointer :: p(:, :, :)
        complex(c_double_complex), pointer :: p_hat(:, :, :)
        type(c_ptr) :: work_r, work_c
        print '(A)', "test_3d_r2c"
        ! allocate meomery
        allocate(x(N), y(N), z(N))
        allocate(p_ref(N, N, N))
        allocate(akx(N), aky(N), akz(N))
        allocate(laplace(N/2 + 1, N, N))

        work_r = fftw_alloc_real(int(N*N*N, c_size_t))
        call c_f_pointer(work_r, p, [N, N, N])
        work_c = fftw_alloc_complex(int((N/2 + 1)*N*N, c_size_t))
        call c_f_pointer(work_c, p_hat, [N/2 + 1, N, N])
        
        ! plan DFT
        fwd_plan = fftw_plan_dft_r2c_3d(N, N, N, p, p_hat, FFTW_ESTIMATE)
        bwd_plan = fftw_plan_dft_c2r_3d(N, N, N, p_hat, p, FFTW_ESTIMATE)
        
        
        ! geometry
        do i = 1, N
            x(i) = real(i - 1, c_double)*dx; y(i) = x(i); z(i) = x(i)
        end do
        do k = 1, N
            do j = 1, N
                do i = 1, N
                    p_ref(i, j, k) = (cos(x(i)) + cos(y(j)))*cos(z(k))
                    ! RHS of nabla^2 p = f
                    p(i, j, k) = -2.0_c_double*p_ref(i, j, k)
                end do
            end do
        end do
        
        call fftw_execute_dft_r2c(fwd_plan, p, p_hat)
        
        call fftfreq(n, akx); call fftfreq(n, aky); call fftfreq(n, akz)
        do k = 1, N
            do j = 1, N
                do i = 1, N/2 + 1
                    if (i == 1 .and. j == 1 .and. k == 1) then
                        laplace(i, j, k) = 0.0_c_double
                    else
                        laplace(i, j, k) = -1.0_c_double/(akx(i)**2 + aky(j)**2 + akz(k)**2)
                    end if
                end do
            end do
        end do
        
        do k = 1, N
            do j = 1, N
                do i = 1, N/2 + 1
                    p_hat(i, j, k) = laplace(i, j, k)*p_hat(i, j, k)
                end do
            end do
        end do
        
        call fftw_execute_dft_c2r(bwd_plan, p_hat, p)
        p(:, :, :) = factor*p(:, :, :)
        
        open(newunit=iunit, file="test_3d_r2c.txt")
        write(iunit, *) norm2(p - p_ref)/norm2(p_ref)
        close(iunit)
        
        call fftw_destroy_plan(fwd_plan)
        call fftw_destroy_plan(bwd_plan)
        call fftw_free(work_r)
        call fftw_free(work_c)

        deallocate(x, y, z, p_ref, akx, aky, akz, laplace)
    end subroutine test_3d_r2c

    subroutine fftfreq(n, ak)
        integer, intent(in) :: n
        real(c_double), intent(out), dimension(:) :: ak

        integer :: l
        do l = 1, n/2
            ak(l) = real(l - 1, c_double)
        end do

        do l = n/2 + 1, n
            ak(l) = real(-(n - l + 1), c_double)
        end do
    end subroutine fftfreq
    
end module m_tests