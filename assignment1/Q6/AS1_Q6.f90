module input_and_grid

	contains
 	subroutine generate_uniform_grid(L,dz,t_span,dt,gridpts_z,gridpts_t, nz,nt)
 	implicit none
	 INTEGER 					:: i
	 integer, intent(out)				:: nz,nt
	 REAL*8,intent(in) 				:: L,dz,t_span,dt
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:),gridpts_t(:)
	 
	 nz = nint(L/dz)
	 nt = nint(t_span/dt)
	 allocate(gridpts_z(0:nz),gridpts_t(0:nt))
	 
	 do i=0,nz 
	   gridpts_z(i)= dz*(i)  
	 end do
	 
	 do i=0,nt
	 	gridpts_t(i)= dt*(i)
	end do

	 print*, "Grid_generation done"," nz",nz," nt",nt
	end subroutine generate_uniform_grid

	subroutine read_input(L,dz,time_span,time_step, diffusivity,T_initial,T_surface)
	  implicit none
	  real*8 	:: L,dz,time_span,time_step, diffusivity,T_initial,T_surface
	 
	  namelist/input_parameters/L,dz,time_span,time_step, diffusivity,T_initial,T_surface
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input

end module input_and_grid


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

module Explicit_solver
	contains
	subroutine explicit_scheme(T_explicit,T_initial,T_surface,nz,nt,alpha,dz,dt)
		implicit none
		integer, intent(in)	::	nz,nt
		real*8, intent(in)	::	T_initial,T_surface,alpha,dz,dt
		integer			::	m,p
		real*8			::	c
		real*8, allocatable,intent(out)	::	T_explicit(:,:)
		
		c=alpha*dt/(dz*dz)
		allocate(T_explicit(0:nz,0:nt))
		T_explicit(1:nz-1,0)=T_initial; T_explicit(0,0:nt)=T_surface; T_explicit(nz,0:nt)=T_surface
		do p=1,nt
			do m=1,nz-1
				T_explicit(m,p)=c*T_explicit(m-1,p-1)+(1.0d0-2*c)*T_explicit(m,p-1)+c*T_explicit(m+1,p-1)
			end do
		end do
		print*,"explicit solution done"
	end subroutine explicit_scheme
	
end module Explicit_solver

module Implicit_solver_module
	contains
	subroutine apply_dir_bc(T_surface,dz,dt,alpha, rhs)
		implicit none
		REAL*8,intent(in) 	:: T_surface,dz,dt,alpha
		real*8,allocatable,intent(inout)	:: rhs(:)
		real*8			:: c
		c=alpha*dt/(dz*dz)
		rhs(1) = rhs(1) - c*T_surface
		rhs(size(rhs)) = rhs(size(rhs)) - c*T_surface
	end subroutine apply_dir_bc
	
	subroutine build_implicit_system(vec_b, vec_d, vec_a,alpha,dz,dt, n)
		implicit none
		integer, intent(in)	::	n
		real*8, intent(in)	::	alpha,dz,dt
		real*8			::	c
		real*8, allocatable,intent(out)	::	vec_b(:)
		real*8, allocatable,intent(out)	::	vec_d(:)
		real*8, allocatable,intent(out)	::	vec_a(:)
		
		c=dt*alpha/(dz*dz)
		allocate(vec_b(1:n),vec_d(1:n),vec_a(1:n))
		vec_b = c
		vec_a = c
		vec_d = -(1.0d0+2*c)
	end subroutine build_implicit_system
	
	subroutine implicit_solver(T_implicit,T_initial,T_surface,nz,nt,alpha,dz,dt)
		use thomas_solver
		implicit none
		integer, intent(in)	::	nz,nt
		real*8, intent(in)	::	T_initial,T_surface,alpha,dz,dt
		integer			::	m,p
		real*8			::	c
		real*8, allocatable,intent(out)	::	T_implicit(:,:)
		real*8, allocatable	::	vec_b(:),vec_d(:),vec_a(:),rhs(:),bd_inv(:),d_inv(:),solution(:)
		
		c=dt*alpha/(dz*dz)
		allocate(rhs(1:nz-1))
		allocate(T_implicit(0:nz,0:nt))
		call build_implicit_system(vec_b,vec_d,vec_a,alpha,dz,dt,nz-1)
		call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
		T_implicit(1:nz-1,0)=T_initial; T_implicit(0,0:nt)=T_surface; T_implicit(nz,0:nt)=T_surface
		do p=1,nt
			rhs(1:nz-1)=-T_implicit(1:nz-1,p-1)
			call apply_dir_bc(T_surface,dz,dt,alpha, rhs)
			call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
			T_implicit(1:nz-1,p)=solution(1:nz-1)
		end do
		print*,"implicit solution done"
	end subroutine implicit_solver

end module Implicit_solver_module

module CN_explicit_implicit_solver
	contains
	subroutine apply_dir_bc_CN(T_surface,dz,dt,alpha, rhs)
		implicit none
		REAL*8,intent(in) 	:: T_surface,dz,dt,alpha
		real*8,allocatable,intent(inout)	:: rhs(:)
		real*8			:: c
		c=alpha*dt/(dz*dz)
		rhs(1) = rhs(1) + c*T_surface/2.0d0
		rhs(size(rhs)) = rhs(size(rhs)) + c*T_surface/2.0d0
	end subroutine apply_dir_bc_CN
	
	subroutine build_imex_system(vec_b, vec_d, vec_a,alpha,dz,dt, n)
		implicit none
		integer, intent(in)	::	n
		real*8, intent(in)	::	alpha,dz,dt
		real*8			::	c
		real*8, allocatable,intent(out)	::	vec_b(:)
		real*8, allocatable,intent(out)	::	vec_d(:)
		real*8, allocatable,intent(out)	::	vec_a(:)
		
		c=dt*alpha/(dz*dz)
		allocate(vec_b(1:n),vec_d(1:n),vec_a(1:n))
		vec_b = -c/2.0d0
		vec_a = -c/2.0d0
		vec_d = (1.0d0+c)
	end subroutine build_imex_system
	
	subroutine imex_solver(T_imex,T_initial,T_surface,nz,nt,alpha,dz,dt)
		use thomas_solver
		implicit none
		integer, intent(in)	::	nz,nt
		real*8, intent(in)	::	T_initial,T_surface,alpha,dz,dt
		integer			::	m,p
		real*8			::	c
		real*8, allocatable,intent(out)	::	T_imex(:,:)
		real*8, allocatable	::	vec_b(:),vec_d(:),vec_a(:),rhs(:),bd_inv(:),d_inv(:),solution(:)
		
		c=dt*alpha/(dz*dz)
		allocate(rhs(1:nz-1))
		allocate(T_imex(0:nz,0:nt))
		call build_imex_system(vec_b,vec_d,vec_a,alpha,dz,dt,nz-1)
		call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
		T_imex(1:nz-1,0)=T_initial; T_imex(0,0:nt)=T_surface; T_imex(nz,0:nt)=T_surface
		do p=1,nt
			rhs(1:nz-1)=c*T_imex(0:nz-2,p-1)/2.0d0 + (1.0d0-c)*T_imex(1:nz-1,p-1) + c*T_imex(2:nz,p-1)/2.0d0
			call apply_dir_bc_CN(T_surface,dz,dt,alpha, rhs)
			call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
			T_imex(1:nz-1,p)=solution(1:nz-1)
		end do
		print*,"imex solution done"
	end subroutine imex_solver

end module CN_explicit_implicit_solver

program OneD_temp_dist
	use input_and_grid
	use thomas_solver
	use Explicit_solver
	use Implicit_solver_module
	use CN_explicit_implicit_solver
	implicit none
	real*8			::	L,dz,time_span,time_step, diffusivity,T_initial,T_surface
	real*8,allocatable	::	gridpts_z(:),gridpts_t(:),T_explicit(:,:),T_imex(:,:),T_implicit(:,:)
	integer			::	nz,nt,m,p
	character(len=10)	::	dt_str
	integer			::	t_array(1:5)
	
	call read_input(L,dz,time_span,time_step, diffusivity,T_initial,T_surface)
	call generate_uniform_grid(L,dz,time_span,time_step,gridpts_z,gridpts_t, nz,nt)
	call explicit_scheme(T_explicit,T_initial,T_surface,nz,nt,diffusivity,dz,time_step)
	call imex_solver(T_imex,T_initial,T_surface,nz,nt,diffusivity,dz,time_step)
	call implicit_solver(T_implicit,T_initial,T_surface,nz,nt,diffusivity,dz,time_step)
	
	write(dt_str,'(F5.4)') time_step
	
	do p=0,nt
		if(abs(gridpts_t(p)-0.1d0)<1.0d-6) t_array(1)=p
		if(abs(gridpts_t(p)-0.2d0)<1.0d-6) t_array(2)=p
		if(abs(gridpts_t(p)-0.3d0)<1.0d-6) t_array(3)=p
		if(abs(gridpts_t(p)-0.4d0)<1.0d-6) t_array(4)=p
		if(abs(gridpts_t(p)-0.5d0)<1.0d-6) t_array(5)=p
	end do
	
	
	open(unit =2, file='Temerature_dist_'//trim(adjustl(dt_str))//'.dat')
	write(2,*),'index',"	",'cordinate',"	",'ex_sol0.1--->0.5',"	",'imex_sol0.1--->0.5',"	",'im_sol0.1--->0.5'
	do m=0,nz
		write(2, '(I0, 1X, F6.4, 15(1X, G0))') m,gridpts_z(m),T_explicit(m,t_array(1:5)),T_imex(m,t_array(1:5)),T_implicit(m,t_array(1:5))
	end do
	close(2)
	
	Print*,"write_to_file Done"

end program OneD_temp_dist
