module input_grid_and_boundary

	contains
 	subroutine generate_uniform_grid_1D(dz,z1,z2,gridpts_z, nz)
 	implicit none
	 INTEGER 					:: i
	 integer, intent(out)				:: nz
	 REAL*8,intent(in) 				:: z1, z2,dz
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:)
	 
	 nz = (z2 - z1)/dz
	 allocate(gridpts_z(0:nz))
	 
	 do i=0,nz 
	   gridpts_z(i)= z1+ dz*(i)  
	 end do

	 print*, "Grid_generation done"
	end subroutine generate_uniform_grid_1D

	subroutine apply_dir_bc(b1, b2,rhs,vec_b,vec_a)
		implicit none
		REAL*8,intent(in) 	:: b1, b2
		real*8,allocatable,intent(inout)	:: rhs(:),vec_b(:),vec_a(:)
		rhs(1) = rhs(1) - vec_b(1)*b1
		rhs(size(rhs)) = rhs(size(rhs)) - vec_a(size(rhs))*b2
	end subroutine apply_dir_bc

	subroutine read_input(L,dz,T_initial,bnd1,bnd2, source,beta,eps_tol,betaVal)
	  implicit none
	  real*8 	:: L,dz,bnd1,bnd2,source,T_initial,beta,eps_tol
	  character	:: betaVal*10
	 
	  namelist/input_parameters/L,dz,bnd1,bnd2,source,T_initial,beta,eps_tol,betaVal
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input

	subroutine create_linear_system(vec_b, vec_d, vec_a, rhs, n,k_arr)
		implicit none
		integer, intent(in)	::	n
		integer			::	m
		real*8, allocatable,intent(in)	::	k_arr(:)
		real*8, allocatable,intent(out)	::	rhs(:)
		real*8, allocatable,intent(out)	::	vec_b(:)
		real*8, allocatable,intent(out)	::	vec_d(:)
		real*8, allocatable,intent(out)	::	vec_a(:)
		
		m = n-2
		allocate(rhs(1:m),vec_b(1:m),vec_d(1:m),vec_a(1:m))
		rhs = 0.0d0
		vec_b(1:m) = 2.0d0/(1.0d0/k_arr(0:m-1)+1.0d0/k_arr(1:m))
		vec_a(1:m) = 2.0d0/(1.0d0/k_arr(1:m)+1.0d0/k_arr(2:m+1))
		vec_d(1:m) = -(vec_b+vec_a)
	end subroutine create_linear_system
	
	subroutine create_k_arr(beta,bnd1,bnd2,solution,k_arr)
		implicit none
		real*8, intent(in)	::	beta,bnd1,bnd2
		real*8, allocatable,intent(in)	::solution(:)
		real*8, allocatable,intent(out)	::k_arr(:)
		allocate(k_arr(0:size(solution)+1))
		k_arr(0) = 0.1d0*(1.0d0+beta*bnd1)
		k_arr(1:size(solution))=0.1d0*(1.0d0+beta*solution)
		k_arr(size(solution)+1)=0.1d0*(1.0d0+beta*bnd2)
		!print*,size(solution),size(k_arr),k_arr(0),"k_arr created"
	end subroutine create_k_arr
	
end module input_grid_and_boundary


module thomas_solver_and_matops
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
	
	subroutine product_lhs_routine(vec_b,vec_d,vec_a,solution,product_lhs)
		implicit none
		real*8, allocatable,intent(in)	::	vec_b(:)
		real*8, allocatable,intent(in)	::	vec_d(:)
		real*8, allocatable,intent(in)	::	vec_a(:)
		real*8, allocatable,intent(in)	::	solution(:)
		real*8, allocatable,intent(out)		::	product_lhs(:)
		integer				::	j,m
		
		m =size(solution)
		allocate(product_lhs(1:m))
		product_lhs(1)=vec_d(1)*solution(1)+vec_a(1)*solution(2)
		product_lhs(m)=vec_d(m)*solution(m)+vec_b(m)*solution(m-1)
		do j = 2,m-1
			product_lhs(j)=vec_b(j)*solution(j-1)+vec_d(j)*solution(j)+vec_a(j)*solution(j+1)
		end do
		
	end subroutine product_lhs_routine
	
end module thomas_solver_and_matops

program Steady_temp_dist

	use input_grid_and_boundary
	use thomas_solver_and_matops
	implicit none
	real*8			::	L,dz,bnd1,bnd2,source,z1,z2,T_initial,beta,eps_tol,stop_value,old_solution_value,new_solution_value
	real*8,allocatable	::	gridpts_z(:)
	real*8,allocatable	::	vec_b(:),vec_d(:),vec_a(:),rhs(:),solution(:),bd_inv(:),d_inv(:),k_arr(:),product_lhs(:),rhs_val(:)
	integer			::	nz,m,iter=0
	character		::	betaVal*10
	call read_input(L,dz,T_initial,bnd1,bnd2, source,beta,eps_tol,betaVal)
	z1=0.0d0; z2=0.0d0+L
	call generate_uniform_grid_1D(dz,z1,z2,gridpts_z, nz)
	allocate(solution(1:nz-1))
	allocate(rhs_val(1:nz-1))
	solution(1:nz-1)=T_initial
	old_solution_value=norm2(solution)
	
	do m=1,40
	
		call create_k_arr(beta,bnd1,bnd2,solution,k_arr)
		deallocate(solution)
		old_solution_value=new_solution_value
		call create_linear_system(vec_b, vec_d, vec_a, rhs, nz+1,k_arr)
		call apply_dir_bc(bnd1, bnd2,rhs,vec_b,vec_a)
		rhs_val=rhs
		call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
		call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
		call product_lhs_routine(vec_b,vec_d,vec_a,solution,product_lhs)
		deallocate(vec_b,vec_d,vec_a,rhs,bd_inv,d_inv,k_arr)
		stop_value=norm2(product_lhs-rhs_val)/norm2(rhs_val)
		new_solution_value=norm2(solution)
		deallocate(product_lhs)
		iter = iter+1
		if(abs(new_solution_value-old_solution_value)<eps_tol) exit
		
	end do
	
	print*,"Solution Done","	iter=",iter,abs(new_solution_value-old_solution_value),stop_value,eps_tol
	m=0
	open(unit =2, file='Temerature_dist_'//trim(betaVal)//'.dat')
	write(2,*),'index',"	",'cordinate',"	",'temperature'
	write(2,*),m, gridpts_z(m),bnd1
	do m=1,nz-1
		write(2,*),m,gridpts_z(m),solution(m)
	end do
	write(2,*),nz,gridpts_z(nz),bnd2
	close(2)
	
	Print*,"write_to_file Done"
	deallocate(gridpts_z,solution)
	

end program Steady_temp_dist
