module input_grid_and_boundary
	contains
 	subroutine generate_uniform_grid(dy,dz,dt,y1,y2,z1,z2,t_ss,gridpts_y,gridpts_z,gridpts_t,ny,nz,nt)
 	implicit none
	 INTEGER 					:: i
	 integer, intent(out)				:: nz,ny,nt
	 REAL*8,intent(in) 				:: z1, z2,y1,y2,t_ss
	 REAL*8, intent(inout)				:: dz,dy,dt
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:), gridpts_y(:), gridpts_t(:)
	 
	 nz = nint((z2 - z1)/dz)
	 ny = nint((y2 - y1)/dy)
	 nt = nint(t_ss/dt)
	 allocate(gridpts_z(0:nz),gridpts_y(0:ny),gridpts_t(0:nt))
	 dz = (z2-z1)/nz
	 dy = (y2-y1)/ny
	 
	 do i=0,nz
	   gridpts_z(i)= z1 + dz*i  
	 end do
	 
	 do i=0,ny 
	   gridpts_y(i)= y1+ dy*i
	 end do
	 
	 do i=0,nt
	 	gridpts_t(i)= dt*i
	 end do
	 
	 print*, "Grid_generation done"
	end subroutine generate_uniform_grid


	subroutine read_input(L,H,k,alpha,Tb,Tl,Tt,Tr,dx,dy,dt,t_ss,w0,wf,w_step,w_fixed,res_tol,T_ini)
	  implicit none
	  real*8 	:: L,H,k,alpha,Tb,Tl,Tt,Tr,dx,dy,dt,t_ss,w0,wf,w_step,w_fixed,res_tol,T_ini
	 
	  namelist/input_parameters/L,H,k,alpha,Tb,Tl,Tt,Tr,dx,dy,w0,wf,w_step,w_fixed,res_tol,T_ini,t_ss,dt
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input
	
	
	subroutine Temp_ini_bc(nx,ny,Tb,Tl,Tt,Tr,T_ini,Temp)
		implicit none
		real*8,intent(in)	::	Tb,Tl,Tt,Tr,T_ini
		integer, intent(in)	::	nx,ny
		real*8,allocatable, intent(out)	::	Temp(:,:)
		
		allocate(Temp(0:nx,0:ny))
		Temp = T_ini
		Temp(0:nx,0)=Tb; Temp(0:nx,ny)=Tt
		Temp(0,1:ny-1)=Tl; Temp(nx,1:ny-1)=Tr
	end subroutine Temp_ini_bc
		
		
end module input_grid_and_boundary


module linear_solvers
	contains
	subroutine create_tridiag_system(vec_b, vec_d, vec_a,rhs, rhs_vec, bnd1, bnd2, c1)
		implicit none
		real*8, intent(in)			::	c1, bnd1, bnd2
		real*8, intent(in)			::	rhs_vec(:)
		real*8, allocatable,intent(out)		::	vec_b(:)
		real*8, allocatable,intent(out)		::	vec_d(:)
		real*8, allocatable,intent(out)		::	vec_a(:)
		real*8, allocatable,intent(out)		::	rhs(:)
		
		allocate(rhs(size(rhs_vec)))
		rhs=rhs_vec
		allocate(vec_b(size(rhs)),vec_d(size(rhs)),vec_a(size(rhs)))
		rhs(1) = rhs(1) - c1*bnd1
		rhs(size(rhs)) = rhs(size(rhs)) - c1*bnd2
		vec_b = c1
		vec_a = c1
		vec_d = +1.0d0
	end subroutine create_tridiag_system
	
	
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
	
	
	subroutine linear_ADI1(nx,ny,dx,dy,dt,res_tol,Temp_in,Temp_out,w,alpha)
		implicit none
		real*8,intent(in)			::	w,dx,dy,res_tol,alpha,dt
		integer, intent(in)			::	nx,ny
		real*8, intent(in)			::	Temp_in(0:nx,0:ny)
		real*8, intent(out)			::	Temp_out(0:nx,0:ny)
		real*8					::	temp_res_norm
		integer					::	j,k,i,n,iter
		real*8					::	Temp_prerelax,c1,c2,r1,r2
		real*8, allocatable			::	vec_b(:), vec_d(:), vec_a(:),rhs(:), bd_inv(:), d_inv(:), solution(:)
		
		r1=alpha*dt/(dx*dx); r2=alpha*dt/(dy*dy)
		c1 = -r1/(1.0d0+2*(r1+r2)); c2 = -r2/(1.0d0+2*(r1+r2))
		iter =0
		n=3000
		Temp_out = Temp_in

		do i=1,n
			temp_res_norm=0.0d0
			do k=1,ny-1
				call create_tridiag_system(vec_b, vec_d, vec_a,rhs, -c2*(Temp_out(1:nx-1,k-1)+Temp_out(1:nx-1,k+1))+Temp_in(1:nx-1,k)/(1.0d0+2*(r1+r2)), Temp_out(0,k), Temp_out(nx,k), c1)
				call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
				call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
				Temp_out(1:nx-1,k) = Temp_out(1:nx-1,k) + w*(solution - Temp_out(1:nx-1,k))
			end do
			do j=1,nx-1
				call create_tridiag_system(vec_b, vec_d, vec_a,rhs, -c1*(Temp_out(j-1,1:ny-1)+Temp_out(j+1,1:ny-1))+Temp_in(j,1:ny-1)/(1.0d0+2*(r1+r2)), Temp_out(j,0), Temp_out(j,ny), c2)
				call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
				call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
				Temp_out(j,1:ny-1) = Temp_out(j,1:ny-1) + w*(solution - Temp_out(j,1:ny-1))
			end do
			
			iter =iter + 1
			do k=1,ny-1
				do j=1,nx-1
					temp_res_norm = temp_res_norm + (c1*(Temp_out(j-1,k)+Temp_out(j+1,k))+c2*(Temp_out(j,k-1)+Temp_out(j,k+1))+Temp_out(j,k)-Temp_in(j,k)/(1.0d0+2*(r1+r2)))**2
				end do
			end do
			temp_res_norm = sqrt(temp_res_norm)
			if(temp_res_norm<res_tol) exit
		end do
		
		print*, "iter: ",iter,"res_norm: ",temp_res_norm
		
	end subroutine linear_ADI1
	
	
	subroutine linear_ADI2(nx,ny,dx,dy,dt,res_tol,Temp_in,Temp_out,w,alpha)
		implicit none
		real*8,intent(in)			::	w,dx,dy,res_tol,alpha,dt
		integer, intent(in)			::	nx,ny
		real*8, intent(in)			::	Temp_in(0:nx,0:ny)
		real*8, intent(out)			::	Temp_out(0:nx,0:ny)
		real*8					::	temp_res_norm
		integer					::	j,k,i,n,iter
		real*8					::	Temp_prerelax,c1,c2,r1,r2,d
		real*8, allocatable			::	vec_b(:), vec_d(:), vec_a(:),rhs(:), bd_inv(:), d_inv(:), solution(:)
		
		r1=alpha*dt/(2*dx*dx); r2=alpha*dt/(2*dy*dy)
		c1 = -r1/(1.0d0+2*(r1+r2)); c2 = -r2/(1.0d0+2*(r1+r2)); d = (1.0d0-2*(r1+r2))/(1.0d0+2*(r1+r2))
		iter =0
		n=3000
		Temp_out = Temp_in

		do i=1,n
			temp_res_norm=0.0d0
			do k=1,ny-1
				call create_tridiag_system(vec_b, vec_d, vec_a,rhs, -c1*(Temp_in(0:nx-2,k)+Temp_in(2:nx,k))-c2*(Temp_out(1:nx-1,k-1)+Temp_out(1:nx-1,k+1)+Temp_in(1:nx-1,k-1)+Temp_in(1:nx-1,k+1))+d*Temp_in(1:nx-1,k), Temp_out(0,k), Temp_out(nx,k), c1)
				call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
				call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
				Temp_out(1:nx-1,k) = Temp_out(1:nx-1,k) + w*(solution - Temp_out(1:nx-1,k))
			end do
			do j=1,nx-1
				call create_tridiag_system(vec_b, vec_d, vec_a,rhs, -c2*(Temp_in(j,0:ny-2)+Temp_in(j,2:ny))-c1*(Temp_out(j-1,1:ny-1)+Temp_out(j+1,1:ny-1)+Temp_in(j-1,1:ny-1)+Temp_in(j+1,1:ny-1))+d*Temp_in(j,1:ny-1), Temp_out(j,0), Temp_out(j,ny), c2)
				call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
				call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
				Temp_out(j,1:ny-1) = Temp_out(j,1:ny-1) + w*(solution - Temp_out(j,1:ny-1))
			end do
			
			iter =iter + 1
			do k=1,ny-1
				do j=1,nx-1
					temp_res_norm = temp_res_norm + (c1*(Temp_out(j-1,k)+Temp_out(j+1,k)+Temp_in(j-1,k)+Temp_in(j+1,k))+c2*(Temp_out(j,k-1)+Temp_out(j,k+1)+Temp_in(j,k-1)+Temp_in(j,k+1))+Temp_out(j,k)-d*Temp_in(j,k))**2
				end do
			end do
			temp_res_norm = sqrt(temp_res_norm)
			if(temp_res_norm<res_tol) exit
		end do
		
		print*, "iter: ",iter,"res_norm: ",temp_res_norm
		
	end subroutine linear_ADI2
	
	
end module linear_solvers


module transient_solvers
	contains
	subroutine explicit_scheme(T_explicit,T_initial,nx,ny,nt,dx,dy,dt,alpha)
		implicit none
		integer, intent(in)		::	nx,ny,nt
		real*8, intent(in)		::	alpha,dx,dy,dt
		real*8, allocatable, intent(in)	::	T_initial(:,:)
		integer				::	j,k,n
		real*8				::	r1,r2
		real*8, allocatable,intent(out)	::	T_explicit(:,:,:)
		
		allocate(T_explicit(0:nx,0:ny,0:nt))
		do n=0,nt
			T_explicit(:,:,n) = T_initial
		end do
		r1 = alpha*dt/(dx*dx); r2 = alpha*dt/(dy*dy)
		do n=1,nt
			do k=1,ny-1
				do j=1,nx-1
					T_explicit(j,k,n) = r1*(T_explicit(j-1,k,n-1)+T_explicit(j+1,k,n-1)) + (1.0d0-2*(r1+r2))*T_explicit(j,k,n-1) + r2*(T_explicit(j,k-1,n-1)+T_explicit(j,k+1,n-1))
				end do
			end do
		end do
		print*,"explicit solution done"
	end subroutine explicit_scheme
	
	
	subroutine implicit_scheme(T_implicit,T_initial,nx,ny,nt,dx,dy,dt,alpha,res_tol,w)
		use linear_solvers
		implicit none
		integer, intent(in)		::	nx,ny,nt
		real*8, intent(in)		::	alpha,dx,dy,dt,res_tol,w
		real*8, allocatable, intent(in)	::	T_initial(:,:)
		integer				::	j,k,n
		real*8				::	r1,r2
		real*8, allocatable,intent(out)	::	T_implicit(:,:,:)
		
		allocate(T_implicit(0:nx,0:ny,0:nt))
		do n=0,nt
			T_implicit(:,:,n) = T_initial
		end do
		
		do n=1,nt
			call linear_ADI1(nx,ny,dx,dy,dt,res_tol,T_implicit(:,:,n-1),T_implicit(:,:,n),w,alpha)
		end do
		
		print*,"implicit solution done"
	end subroutine implicit_scheme
	
	
	subroutine CN_imex_scheme(T_CN_imex,T_initial,nx,ny,nt,dx,dy,dt,alpha,res_tol,w)
		use linear_solvers
		implicit none
		integer, intent(in)		::	nx,ny,nt
		real*8, intent(in)		::	alpha,dx,dy,dt,res_tol,w
		real*8, allocatable, intent(in)	::	T_initial(:,:)
		integer				::	j,k,n
		real*8				::	r1,r2
		real*8, allocatable,intent(out)	::	T_CN_imex(:,:,:)
		
		allocate(T_CN_imex(0:nx,0:ny,0:nt))
		do n=0,nt
			T_CN_imex(:,:,n) = T_initial
		end do
		
		do n=1,nt
			call linear_ADI2(nx,ny,dx,dy,dt,res_tol,T_CN_imex(:,:,n-1),T_CN_imex(:,:,n),w,alpha)
		end do
		print*,"CN_imex solution done"
	end subroutine CN_imex_scheme
	
	
	subroutine ADI_scheme(T_ADI,T_initial,nx,ny,nt,dx,dy,dt,alpha)
		use linear_solvers
		implicit none
		integer, intent(in)		::	nx,ny,nt
		real*8, intent(in)		::	alpha,dx,dy,dt
		real*8, allocatable, intent(in)	::	T_initial(:,:)
		integer				::	j,k,n
		real*8				::	rx,ry,cxx,cxy,cyy
		real*8, allocatable,intent(out)	::	T_ADI(:,:,:)
		real*8, allocatable		::	vec_b(:), vec_d(:), vec_a(:),rhs(:), bd_inv(:), d_inv(:), solution(:), du(:,:)
		
		allocate(T_ADI(0:nx,0:ny,0:nt), du(0:nx,0:ny))
		do n=0,nt
			T_ADI(:,:,n) = T_initial
		end do
		du=0.0d0
		rx = alpha*dt/(dx*dx); ry = alpha*dt/(dy*dy)
		cxx = -rx/(1+rx)/2.0d0; cxy = -ry/(1+rx)/2.0d0; cyy = -ry/(1+ry)/2.0d0
		
		do n=1,nt
			do k=1,ny-1
				call create_tridiag_system(vec_b, vec_d, vec_a,rhs, -2*cxx*(T_ADI(0:nx-2,k,n-1)+T_ADI(2:nx,k,n-1))+4*(cxx+cxy)*T_ADI(1:nx-1,k,n-1)-2*cxy*(T_ADI(1:nx-1,k-1,n-1)+T_ADI(1:nx-1,k+1,n-1)), 0.0d0, 0.0d0, cxx)
				call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
				call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
				du(1:nx-1,k) = solution
			end do
			do j=1,nx-1
				call create_tridiag_system(vec_b, vec_d, vec_a,rhs, du(j,1:ny-1)/(1.0d0+ry), 0.0d0, 0.0d0, cyy)
				call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
				call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
				du(j,1:ny-1) = solution
			end do
			T_ADI(1:nx-1,1:ny-1,n) = T_ADI(1:nx-1,1:ny-1,n-1) + du(1:nx-1,1:ny-1)
		end do
		
		print*,"ADI solution done"
	end subroutine ADI_scheme
	
end module transient_solvers


program diffusion
	use input_grid_and_boundary
	use linear_solvers
	use transient_solvers
	implicit none
	real*8			::	L,H,cond_k,alpha,Tb,Tl,Tt,Tr,dx,dy,dt,t_ss,w0,wf,w_step,w_fixed,res_tol,T_ini,x1,x2,y1,y2,w
	real*8,allocatable	::	gridpts_x(:), gridpts_y(:),gridpts_t(:),Temp_ini(:,:)
	real*8,allocatable	::	T_explicit(:,:,:),T_implicit(:,:,:),T_CN_imex(:,:,:), T_ADI(:,:,:)
	integer			::	nx,ny,nt,j,k
	character(len=10)	::	dt_str
	
	
	call read_input(L,H,cond_k,alpha,Tb,Tl,Tt,Tr,dx,dy,dt,t_ss,w0,wf,w_step,w_fixed,res_tol,T_ini)
	
	x1=0.0d0; x2=x1+L; y1=0.0d0; y2=y1+H
	call generate_uniform_grid(dx,dy,dt,x1,x2,y1,y2,t_ss,gridpts_x,gridpts_y,gridpts_t,nx,ny,nt)
	call Temp_ini_bc(nx,ny,Tb,Tl,Tt,Tr,T_ini,Temp_ini)

	w = w_fixed
	call explicit_scheme(T_explicit,Temp_ini,nx,ny,nt,dx,dy,dt,alpha)
	call implicit_scheme(T_implicit,Temp_ini,nx,ny,nt,dx,dy,dt,alpha,res_tol,w)
	call CN_imex_scheme(T_CN_imex,Temp_ini,nx,ny,nt,dx,dy,dt,alpha,res_tol,w)
	call ADI_scheme(T_ADI,Temp_ini,nx,ny,nt,dx,dy,dt,alpha)
	
	write(dt_str,'(F3.1)') dt
	open(unit=2,file='Temp_dist_explicit_dt_'//trim(adjustl(dt_str))//'.dat')
	open(unit=3,file='Temp_dist_implicit_dt_'//trim(adjustl(dt_str))//'.dat')
	open(unit=4,file='Temp_dist_CN_imex_dt_'//trim(adjustl(dt_str))//'.dat')
	open(unit=5,file='Temp_dist_ADI_dt_'//trim(adjustl(dt_str))//'.dat')
	
	write(2,'(*(6X,A18))'),"grid_pt_x","grid_pt_y","Temp_dist-->t_ss"
	write(3,'(*(6X,A18))'),"grid_pt_x","grid_pt_y","Temp_dist-->t_ss"
	write(4,'(*(6X,A18))'),"grid_pt_x","grid_pt_y","Temp_dist-->t_ss"
	write(5,'(*(6X,A18))'),"grid_pt_x","grid_pt_y","Temp_dist-->t_ss"
	do k=0,ny
		do j=0,nx
			write(2,'(1X, F6.4, 2X, F6.4, *(2X, F10.6))'),gridpts_x(j),gridpts_y(k),T_explicit(j,k,:)
			write(3,'(1X, F6.4, 2X, F6.4, *(2X, F10.6))'),gridpts_x(j),gridpts_y(k),T_implicit(j,k,:)
			write(4,'(1X, F6.4, 2X, F6.4, *(2X, F10.6))'),gridpts_x(j),gridpts_y(k),T_CN_imex(j,k,:)
			write(5,'(1X, F6.4, 2X, F6.4, *(2X, F10.6))'),gridpts_x(j),gridpts_y(k),T_ADI(j,k,:)
		end do
		write(2,*), "	"
		write(3,*), "	"
		write(4,*), "	"
		write(5,*), "	"
	end do
	
	close(2)
	close(3)
	close(4)
	close(5)
		
end program diffusion
