module input_grid_and_boundary

	contains
 	subroutine generate_uniform_grid_1D(dz,z1,z2,gridpts_z, nz)
 	implicit none
	 INTEGER 					:: i
	 integer, intent(out)				:: nz
	 REAL*8,intent(in) 				:: z1, z2,dz
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:)
	 print*,dz,(z2 - z1)/dz
	 nz =nint((z2 - z1)/dz)
	 allocate(gridpts_z(0:nz))
	 print*,nz
	 do i=0,nz 
	   gridpts_z(i)= z1+ dz*(i)  
	 end do

	 print*, "Grid_generation done"
	end subroutine generate_uniform_grid_1D

	subroutine apply_dir_bc(b1, b2,rhs)
		implicit none
		REAL*8,intent(in) 	:: b1, b2
		real*8,allocatable,intent(inout)	:: rhs(:)
		rhs(1) = rhs(1) - b1
		rhs(size(rhs)) = rhs(size(rhs)) - b2
	end subroutine apply_dir_bc

	subroutine read_input(L,dz,k1,k2,k3,bnd1,bnd2, source,x_i1,x_i2)
	  implicit none
	  real*8 	:: L,dz,k1,k2,k3,bnd1,bnd2,source,x_i1,x_i2
	 
	  namelist/input_parameters/L,dz,k1,k2,k3,bnd1,bnd2,source,x_i1,x_i2
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input

	subroutine create_linear_system(vec_b, vec_d, vec_a, rhs, n, source,dz,k1,k2,k3,x_i1,x_i2,gridpts_z)
		implicit none
		integer, intent(in)	::	n
		integer			::	m,j
		real*8, intent(in)	::	source,dz,k1,k2,k3,x_i1,x_i2
		REAL(kind=8),allocatable, intent(in) 		:: gridpts_z(:)
		real*8, allocatable,intent(out)	::	rhs(:)
		real*8, allocatable,intent(out)	::	vec_b(:)
		real*8, allocatable,intent(out)	::	vec_d(:)
		real*8, allocatable,intent(out)	::	vec_a(:)
		
		m = n-2
		allocate(rhs(1:m),vec_b(1:m),vec_d(1:m),vec_a(1:m))
		rhs = 0.0d0!-source*dz*dz/k
		vec_b = 1.0d0
		vec_a = 1.0d0
		vec_d = -2.0d0
		vec_b(1) = 0.0d0
		vec_a(m) = 0.0d0

		do j=1,m
			if(abs(gridpts_z(j)-x_i1) < 1.0d-8) then
				vec_b(j) = k1
				vec_d(j) = -(k1+k2)
				vec_a(j) = k2
				print*,vec_d(j)," x_i1 done",j
			else if(abs(gridpts_z(j)-x_i2) < 1.0d-8) then
				vec_b(j) = k2
				vec_d(j) = -(k3+k2)
				vec_a(j) = k3
				print*,vec_d(j)," x_i2 done",j
			end if
		end do
	end subroutine create_linear_system
	
end module input_grid_and_boundary


module thomas_solver
	contains
	subroutine vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
		implicit none
		real*8, allocatable,intent(in)	::	vec_b(:)
		real*8, allocatable,intent(in)	::	vec_d(:)
		real*8, allocatable,intent(in)	::	vec_a(:)
		real*8, allocatable,intent(out)	::	bd_inv(:)
		real*8, allocatable,intent(out)	::	d_inv(:)
		integer			::	i,N
		
		N = size(vec_d)
		allocate(bd_inv(1:N),d_inv(1:N))
		d_inv(1) = 1.0d0/vec_d(1)
		bd_inv(1) = 0.0d0
		do i =2,N
			bd_inv(i) = vec_b(i)*d_inv(i-1)
			d_inv(i) = 1.0d0/(vec_d(i)-bd_inv(i)*vec_a(i-1))
		end do
	end subroutine vectors_out_for_solver
	
	subroutine back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
		implicit none
		real*8, allocatable,intent(inout)	::	rhs(:)
		real*8, allocatable,intent(in)	::	bd_inv(:)
		real*8, allocatable,intent(in)	::	d_inv(:)
		real*8, allocatable,intent(in)	::	vec_a(:)
		real*8, allocatable,intent(out)	::	solution(:)
		integer			::	i,N
		
		N=size(rhs)
		allocate(solution(1:N))
		do i =2,N
			rhs(i) = rhs(i)-bd_inv(i)*rhs(i-1)
		end do
		
		solution(N) = rhs(N)*d_inv(N)
		do i=N-1,1,-1
			solution(i) = (rhs(i) - vec_a(i)*solution(i+1))*d_inv(i)
		end do
	end subroutine back_substitution
			
end module thomas_solver

program Steady_temp_dist
	use input_grid_and_boundary
	use thomas_solver
	implicit none
	real*8			::	L,dz,k1,k2,k3,bnd1,bnd2,source,z1,z2,x_i1,x_i2
	real*8,allocatable	::	gridpts_z(:)
	real*8,allocatable	::	vec_b(:),vec_d(:),vec_a(:),rhs(:),solution(:),bd_inv(:),d_inv(:)
	integer			::	nz,m
	
	call read_input(L,dz,k1,k2,k3,bnd1,bnd2, source,x_i1,x_i2)
	z1=0.0d0; z2=0.0d0+L
	call generate_uniform_grid_1D(dz,z1,z2,gridpts_z, nz)
	call create_linear_system(vec_b, vec_d, vec_a, rhs, nz+1, source,dz,k1,k2,k3,x_i1,x_i2,gridpts_z)
	m=0
	call apply_dir_bc(bnd1, bnd2, rhs)
	call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
	call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
	
	print*,"Solution Done"
	
	open(unit =2, file='Temerature_dist.dat')
	write(2,*),'index',"	",'cordinate',"	",'temperature'
	write(2,*),m, gridpts_z(m),bnd1
	do m=1,nz-1
		write(2,*),m,gridpts_z(m),solution(m)
	end do
	write(2,*),nz,gridpts_z(nz),bnd2
	close(2)
	
	Print*,"write_to_file Done"
	deallocate(gridpts_z,vec_b,vec_d,vec_a,rhs,solution,bd_inv,d_inv)
	

end program Steady_temp_dist
