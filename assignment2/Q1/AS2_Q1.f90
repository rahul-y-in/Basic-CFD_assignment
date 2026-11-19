module input_grid_and_boundary
	contains
 	subroutine generate_uniform_grid(dy,dz,y1,y2,z1,z2,gridpts_y,gridpts_z,ny,nz)
 	implicit none
	 INTEGER 					:: i
	 integer, intent(out)				:: nz,ny
	 REAL*8,intent(in) 				:: z1, z2,y1,y2
	 REAL*8, intent(inout)				:: dz,dy
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:), gridpts_y(:)
	 
	 nz = nint((z2 - z1)/dz)
	 ny = nint((y2 - y1)/dy)
	 allocate(gridpts_z(0:nz),gridpts_y(0:ny))
	 dz = (z2-z1)/nz
	 dy = (y2-y1)/ny
	 
	 do i=0,nz
	   gridpts_z(i)= z1 + dz*i  
	 end do
	 
	 do i=0,ny 
	   gridpts_y(i)= y1+ dy*i
	 end do
	 
	 print*, "Grid_generation done"
	end subroutine generate_uniform_grid


	subroutine read_input(L,H,k,alpha,Tb,Tl,Tt,Tr,dx,dy,w0,wf,w_step,w_fixed,res_tol,T_ini)
	  implicit none
	  real*8 	:: L,H,k,alpha,Tb,Tl,Tt,Tr,dx,dy,w0,wf,w_step,w_fixed,res_tol,T_ini
	 
	  namelist/input_parameters/L,H,k,alpha,Tb,Tl,Tt,Tr,dx,dy,w0,wf,w_step,w_fixed,res_tol,T_ini
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


module solvers
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
		vec_d = -1.0d0
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
	
	
	subroutine point_GS(nx,ny,dx,dy,res_tol,Temp_in,Temp_out,w,res_norm,iter)
		implicit none
		real*8,intent(in)			::	w,dx,dy,res_tol
		integer, intent(in)			::	nx,ny
		real*8, allocatable, intent(in)		::	Temp_in(:,:)
		real*8, allocatable, intent(out)	::	Temp_out(:,:),res_norm(:)
		real*8, allocatable			::	temp_res_norm(:)
		integer					::	j,k,i,n
		integer,intent(out)			::	iter
		real*8					::	Temp_prerelax,c1,c2
		
		c1=(dy**2)/(2*(dx**2+dy**2)); c2=(dx**2)/(2*(dx**2+dy**2))
		iter =0
		n=3000
		allocate(Temp_out(0:nx,0:ny), temp_res_norm(n))
		Temp_out = Temp_in
		temp_res_norm=0.0d0

		do i=1,n
			do k=1,ny-1
				do j=1,nx-1
					Temp_prerelax=c1*(Temp_out(j-1,k)+Temp_out(j+1,k))+c2*(Temp_out(j,k-1)+Temp_out(j,k+1))
					Temp_out(j,k)=Temp_out(j,k)+w*(Temp_prerelax-Temp_out(j,k))
				end do
			end do
			iter =iter + 1
			do k=1,ny-1
				do j=1,nx-1
					temp_res_norm(iter) = temp_res_norm(iter) + (c1*(Temp_out(j-1,k)+Temp_out(j+1,k))+c2*(Temp_out(j,k-1)+Temp_out(j,k+1))-Temp_out(j,k))**2
				end do
			end do
			temp_res_norm(iter) = sqrt(temp_res_norm(iter))
			if(temp_res_norm(iter)<res_tol) exit
		end do
		allocate(res_norm(iter))
		res_norm(1:iter)=temp_res_norm(1:iter)
		deallocate(temp_res_norm)
		print*, "iter: ",iter,"res_norm: ",res_norm(iter)
		
	end subroutine point_GS
	
	
	subroutine line_GS(nx,ny,dx,dy,res_tol,Temp_in,Temp_out,w,res_norm,iter)
		implicit none
		real*8,intent(in)			::	w,dx,dy,res_tol
		integer, intent(in)			::	nx,ny
		real*8, allocatable, intent(in)		::	Temp_in(:,:)
		real*8, allocatable, intent(out)	::	Temp_out(:,:),res_norm(:)
		real*8, allocatable			::	temp_res_norm(:)
		integer					::	j,k,i,n
		integer,intent(out)			::	iter
		real*8					::	Temp_prerelax,c1,c2
		real*8, allocatable			::	vec_b(:), vec_d(:), vec_a(:),rhs(:), bd_inv(:), d_inv(:), solution(:)
		
		c1=(dy**2)/(2*(dx**2+dy**2)); c2=(dx**2)/(2*(dx**2+dy**2))
		iter =0
		n=3000
		
		allocate(Temp_out(0:nx,0:ny), temp_res_norm(n))
		Temp_out = Temp_in
		temp_res_norm=0.0d0

		do i=1,n
			do k=1,ny-1
				call create_tridiag_system(vec_b, vec_d, vec_a,rhs, -c2*(Temp_out(1:nx-1,k-1)+Temp_out(1:nx-1,k+1)), Temp_out(0,k), Temp_out(nx,k), c1)
				call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
				call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
				Temp_out(1:nx-1,k) = Temp_out(1:nx-1,k) + w*(solution - Temp_out(1:nx-1,k))
			end do
			iter =iter + 1
			do k=1,ny-1
				do j=1,nx-1
					temp_res_norm(iter) = temp_res_norm(iter) + (c1*(Temp_out(j-1,k)+Temp_out(j+1,k))+c2*(Temp_out(j,k-1)+Temp_out(j,k+1))-Temp_out(j,k))**2
				end do
			end do
			temp_res_norm(iter) = sqrt(temp_res_norm(iter))
			if(temp_res_norm(iter)<res_tol) exit
		end do
		allocate(res_norm(iter))
		res_norm(1:iter)=temp_res_norm(1:iter)
		deallocate(temp_res_norm)
		print*, "iter: ",iter,"res_norm: ",res_norm(iter)
		
	end subroutine line_GS
	
	
	subroutine ADI(nx,ny,dx,dy,res_tol,Temp_in,Temp_out,w,res_norm,iter)
		implicit none
		real*8,intent(in)			::	w,dx,dy,res_tol
		integer, intent(in)			::	nx,ny
		real*8, allocatable, intent(in)		::	Temp_in(:,:)
		real*8, allocatable, intent(out)	::	Temp_out(:,:),res_norm(:)
		real*8, allocatable			::	temp_res_norm(:)
		integer					::	j,k,i,n
		integer,intent(out)			::	iter
		real*8					::	Temp_prerelax,c1,c2
		real*8, allocatable			::	vec_b(:), vec_d(:), vec_a(:),rhs(:), bd_inv(:), d_inv(:), solution(:)
		
		c1=(dy**2)/(2*(dx**2+dy**2)); c2=(dx**2)/(2*(dx**2+dy**2))
		iter =0
		n=3000
		
		allocate(Temp_out(0:nx,0:ny), temp_res_norm(n))
		Temp_out = Temp_in
		temp_res_norm=0.0d0

		do i=1,n
			do k=1,ny-1
				call create_tridiag_system(vec_b, vec_d, vec_a,rhs, -c2*(Temp_out(1:nx-1,k-1)+Temp_out(1:nx-1,k+1)), Temp_out(0,k), Temp_out(nx,k), c1)
				call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
				call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
				Temp_out(1:nx-1,k) = Temp_out(1:nx-1,k) + w*(solution - Temp_out(1:nx-1,k))
			end do
			do j=1,nx-1
				call create_tridiag_system(vec_b, vec_d, vec_a,rhs, -c1*(Temp_out(j-1,1:ny-1)+Temp_out(j+1,1:ny-1)), Temp_out(j,0), Temp_out(j,ny), c2)
				call vectors_out_for_solver(vec_b, vec_d, vec_a, bd_inv, d_inv)
				call back_substitution(rhs, bd_inv, d_inv, vec_a, solution)
				Temp_out(j,1:ny-1) = Temp_out(j,1:ny-1) + w*(solution - Temp_out(j,1:ny-1))
			end do
			
			iter =iter + 1
			do k=1,ny-1
				do j=1,nx-1
					temp_res_norm(iter) = temp_res_norm(iter) + (c1*(Temp_out(j-1,k)+Temp_out(j+1,k))+c2*(Temp_out(j,k-1)+Temp_out(j,k+1))-Temp_out(j,k))**2
				end do
			end do
			temp_res_norm(iter) = sqrt(temp_res_norm(iter))
			if(temp_res_norm(iter)<res_tol) exit
		end do
		allocate(res_norm(iter))
		res_norm(1:iter)=temp_res_norm(1:iter)
		deallocate(temp_res_norm)
		print*, "iter: ",iter,"res_norm: ",res_norm(iter)
		
	end subroutine ADI			
end module solvers


program diffusion
	use input_grid_and_boundary
	use solvers
	implicit none
	real*8			::	L,H,cond_k,alpha,Tb,Tl,Tt,Tr,dx,dy,w0,wf,w_step,w_fixed,res_tol,T_ini,x1,x2,y1,y2,w
	real*8,allocatable	::	gridpts_x(:), gridpts_y(:),Temp_ini(:,:)
	real*8,allocatable	::	final_Temp1(:,:),res_norm1(:),final_Temp2(:,:),res_norm2(:),final_Temp3(:,:),res_norm3(:)
	integer			::	nx,ny,j,k,iter1,iter2,iter3
	character(len=10)	::	w_str
	
	
	call read_input(L,H,cond_k,alpha,Tb,Tl,Tt,Tr,dx,dy,w0,wf,w_step,w_fixed,res_tol,T_ini)
	
	x1=0.0d0; x2=x1+L; y1=0.0d0; y2=y1+H
	call generate_uniform_grid(dx,dy,x1,x2,y1,y2,gridpts_x,gridpts_y,nx,ny)
	call Temp_ini_bc(nx,ny,Tb,Tl,Tt,Tr,T_ini,Temp_ini)
	
	w = w_fixed
	call point_GS(nx,ny,dx,dy,res_tol,Temp_ini,final_Temp1,w,res_norm1,iter1)
	call line_GS(nx,ny,dx,dy,res_tol,Temp_ini,final_Temp2,w,res_norm2,iter2)
	call ADI(nx,ny,dx,dy,res_tol,Temp_ini,final_Temp3,w,res_norm3,iter3)
	
	print*, "Iter 1,2,3 ",iter1,iter2,iter3, " res_norm 1,2,3", res_norm1(iter1),res_norm2(iter2),res_norm3(iter3)
	
	open(unit=2,file='Temp_dist.dat')
	open(unit=3,file='iter_vs_res_norm_pointGS.dat')
	open(unit=4,file='iter_vs_res_norm_lineGS.dat')
	open(unit=5,file='iter_vs_res_norm_ADI.dat')
	
	write(2,'(*(6X,A10))'),"grid_pt_X","grid_pt_y","point_GS","line_GS","ADI"
	do k=0,ny
		do j=0,nx
			write(2,'(1X, F6.4, 3X, F6.4, *(3X, F10.6))'),gridpts_x(j),gridpts_y(k),final_Temp1(j,k),final_Temp2(j,k),final_Temp3(j,k)
		end do
		write(2,*), "	"
	end do
	
	write(3,*),"iter	","res_norm"
	write(4,*),"iter	","res_norm"
	write(5,*),"iter	","res_norm"
	do j=1,max(iter1,iter2,iter3)
		if(j<=iter1) then
			write(3,"(1X, I5,3X, F10.6)"),j,res_norm1(j)
		end if
		if(j<=iter2) then
			write(4,"(1X, I5,3X, F10.6)"),j,res_norm2(j)
		end if
		if(j<=iter3) then
			write(5,"(1X, I5,3X, F10.6)"),j,res_norm3(j)
		end if
	end do
	close(2)
	close(3)
	close(4)
	close(5)
	
	open(unit=7, file="iter_and_residuals.dat")
	write(7,*),"iters	","residual	","PointGS-->LineGS-->ADI"
	w=w0
	do while (w0<=w .and. w<=wf)
		write(w_str,'(F3.1)') w
		open(unit=6,file='Temp_dist_w_'//trim(adjustl(w_str))//'.dat')
		call point_GS(nx,ny,dx,dy,res_tol,Temp_ini,final_Temp1,w,res_norm1,iter1)
		call line_GS(nx,ny,dx,dy,res_tol,Temp_ini,final_Temp2,w,res_norm2,iter2)
		call ADI(nx,ny,dx,dy,res_tol,Temp_ini,final_Temp3,w,res_norm3,iter3)
		
		write(6,'(*(6X,A10))'),"grid_pt_X","grid_pt_y","point_GS","line_GS","ADI"
		do k=0,ny
			do j=0,nx
				write(6,'(1X, F6.4, 3X, F6.4, *(3X, F10.6))'),gridpts_x(j),gridpts_y(k),final_Temp1(j,k),final_Temp2(j,k),final_Temp3(j,k)
			end do
			write(6,*), "	"
		end do
		close(6)
		
		write(7,*),"**************	","w = ",trim(adjustl(w_str)),"	*************"
		write(7,"(1X, I5,3X, F10.6)"),iter1,res_norm1(iter1)
		write(7,"(1X, I5,3X, F10.6)"),iter2,res_norm2(iter2)
		write(7,"(1X, I5,3X, F10.6)"),iter3,res_norm3(iter3)
		
		w=w+w_step
	end do
	close(7)
		
end program diffusion
