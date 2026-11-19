module input_grid_and_boundary

	contains
 	subroutine generate_uniform_grid_1D(dz,z1,z2,gridpts_z, nz)
 	implicit none
	 INTEGER 					:: i
	 integer, intent(out)				:: nz
	 REAL*8,intent(in) 				:: z1, z2,dz
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:)
	 
	 nz = nint((z2 - z1)/dz)
	 allocate(gridpts_z(0:nz))
	 
	 do i=0,nz 
	   gridpts_z(i)= z1+ dz*(i)  
	 end do

	 print*, "Grid_generation done"
	end subroutine generate_uniform_grid_1D

	subroutine apply_dir_bc(b1, b2, dz,k, rhs)
		implicit none
		REAL*8,intent(in) 	:: b1, b2,dz,k
		real*8,allocatable,intent(inout)	:: rhs(:)
		rhs(1) = rhs(1) - b1
		rhs(size(rhs)) = rhs(size(rhs)) - b2
	end subroutine apply_dir_bc

	subroutine read_input(L,dz,k,bnd1,bnd2, source)
	  implicit none
	  real*8 	:: L,dz,k,bnd1,bnd2,source
	 
	  namelist/input_parameters/L,dz,k,bnd1,bnd2,source
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input

	subroutine create_linear_system(vec_b, vec_d, vec_a, rhs, n,k,dz,source)
		implicit none
		integer, intent(in)	::	n
		real*8, intent(in)	::	dz, source,k
		integer			::	m
		real*8, allocatable,intent(out)	::	rhs(:)
		real*8, allocatable,intent(out)	::	vec_b(:)
		real*8, allocatable,intent(out)	::	vec_d(:)
		real*8, allocatable,intent(out)	::	vec_a(:)
		
		m = n-2
		allocate(rhs(1:m),vec_b(1:m),vec_d(1:m),vec_a(1:m))
		rhs = -source*dz*dz/k
		vec_b = 1.0d0
		vec_a = 1.0d0
		vec_d = -2.0d0
		vec_b(1) = 0.0d0
		vec_a(m) = 0.0d0
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
	real*8			::	L,dz,k,bnd1,bnd2,source,z1,z2
	real*8,allocatable	::	gridpts_z(:)
	real*8,allocatable	::	vec_b(:),vec_d(:),vec_a(:),rhs(:),solution(:),bd_inv(:),d_inv(:)
	integer			::	nz,m
	
	call read_input(L,dz,k,bnd1,bnd2, source)
	z1=0.0d0; z2=0.0d0+L
	call generate_uniform_grid_1D(dz,z1,z2,gridpts_z, nz)
	call create_linear_system(vec_b, vec_d, vec_a, rhs, nz+1,k,dz,source)
	m=0
	call apply_dir_bc(bnd1, bnd2, dz,k, rhs)
	call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
	call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
	
	print*,"Solution Done"
	
	open(unit =2, file='Temerature_dist2.dat')
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
