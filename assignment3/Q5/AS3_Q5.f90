module input_grid_and_boundary
	contains
 	subroutine generate_uniform_grid(dz,z1,z2,gridpts_z,nz)
 	implicit none
	 INTEGER 					:: i
	 integer, intent(out)				:: nz
	 REAL*8,intent(in) 				:: z1, z2
	 REAL*8, intent(in)				:: dz
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:)
	 
	 nz = nint((z2 - z1)/dz)
	 allocate(gridpts_z(0:nz))

	 do i=0,nz
	   gridpts_z(i)= z1 + dz*i  
	 end do
	
	 print*, "Grid_generation done: z(nz)= ",gridpts_z(nz)
	end subroutine generate_uniform_grid


	subroutine read_input(z1,z2,f1,f2,f2_inf,f31,f32,dz,res_tol,max_iter,search_val)
	  implicit none
	  real*8 		:: z1,z2,f1,f2,f2_inf,f31,f32,dz,res_tol,search_val
	  integer		:: max_iter
	  
	  namelist/input_parameters/z1,z2,f1,f2,f2_inf,f31,f32,dz,res_tol,max_iter,search_val
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input
	
	
	subroutine initial_cond(u,f1,f2,f3)
		implicit none
		real*8,intent(inout)	:: u(1:3)
		real*8, intent(in)	:: f1,f2,f3
		u(1)=f1;u(2)=f2;u(3)=f3
	end subroutine initial_cond
	

end module input_grid_and_boundary


module bnd_layer_solvers
	contains
	
	subroutine eval_F(u,F)			!Routine to calculate derivatives of u
		implicit none
		real*8, intent(in)		:: u(:)
		real*8, allocatable,intent(out)	:: F(:)
		allocate(F(1:3))
		F(1)=u(2);F(2)=u(3);F(3)=-u(1)*u(3)/2.0d0
	end subroutine eval_F
	
	
	subroutine RK4(U,nz,dz,f1,f2,f3)
		use input_grid_and_boundary
		implicit none
		integer, intent(in)			:: nz
		real*8, intent(in)			:: dz,f1,f2,f3
		real*8, allocatable, intent(out)	:: U(:,:)
		real*8, allocatable			:: K1(:),K2(:),K3(:),K4(:),temp_u(:)
		integer					:: j
		
		allocate(U(1:3,0:nz),K1(1:3),K2(1:3),K3(1:3),K4(1:3),temp_u(1:3))
		call initial_cond(U(1:3,0),f1,f2,f3)
		
		do j=0,nz-1
			call eval_F(U(1:3,j),K1)
			temp_u(1:3)=U(1:3,j)+(dz/2.0d0)*K1(1:3)
			call eval_F(temp_u,K2)
			temp_u(1:3)=U(1:3,j)+(dz/2.0d0)*K2(1:3)
			call eval_F(temp_u,K3)
			temp_u(1:3)=U(1:3,j)+dz*K3(1:3)
			call eval_F(temp_u,K4)
			U(1:3,j+1)=U(1:3,j)+(dz/6.0d0)*(K1(1:3)+2.0d0*K2(1:3)+2.0d0*K3(1:3)+K4(1:3))
		end do
	end subroutine RK4
	
	
	subroutine bisection_and_solver(f31,f32,f2_inf,res_tol,max_iter,nz,dz,U,f1,f2)
		use input_grid_and_boundary
		implicit none
		integer, intent(in)			:: nz,max_iter
		real*8, intent(in)			:: f31,f32,f2_inf,res_tol,dz,f1,f2
		real*8, allocatable, intent(out)	:: U(:,:)
		real*8, allocatable			:: V(:,:),W(:,:)
		integer					:: n
		real*8					:: l_lim,r_lim,t_lim,iter_count
		
		call RK4(V,nz,dz,f1,f2,f31)
		call RK4(W,nz,dz,f1,f2,f32)
		if((V(2,nz)-f2_inf)*(W(2,nz)-f2_inf)>0.0d0) then
			print*,"Both the initial conditions are on the same side. Aborted!"
			return
		elseif((V(2,nz)-f2_inf)==0.0d0) then
			U=V
			print*,"f31 condition exact. solved!"
			return
		elseif((W(2,nz)-f2_inf)==0.0d0) then
			U=W
			print*,"f32 condition exact. solved!"
			return
		end if
		
		l_lim=f31;r_lim=f32
		do n=1,max_iter
			t_lim=(l_lim+r_lim)/2.0d0
			call RK4(U,nz,dz,f1,f2,t_lim)
			iter_count=n
			if(abs(U(2,nz)-f2_inf)<res_tol) exit
			if((V(2,nz)-f2_inf)*(U(2,nz)-f2_inf)>0.0d0) then
				l_lim=t_lim
				V=U
			else
				r_lim=t_lim
				!W=U
			end if
		end do
	print*,"solution done in ",iter_count
	end subroutine bisection_and_solver
	
	
	subroutine search(gridpts_z,z_val,U,search_val)
		real*8, intent(in)		:: search_val
		real*8, intent(out)		:: z_val
		integer				:: j
		real*8, allocatable,intent(in)	:: gridpts_z(:),U(:,:)
		do j=1,size(gridpts_z)-1
			if(U(2,j)>=search_val) then
				z_val=gridpts_z(j)
				print*,"u/U = ",search_val,"is at ",z_val
				return
			end if
		end do
	end subroutine search
	
	
end module bnd_layer_solvers


program bnd_layer
	use input_grid_and_boundary
	use bnd_layer_solvers
	implicit none
	real*8 			:: z1,z2,f1,f2,f2_inf,f31,f32,dz,res_tol,z_val,search_val
	integer			:: nz,j,k,max_iter
	real*8, allocatable	:: U(:,:),gridpts_z(:)
	
	call read_input(z1,z2,f1,f2,f2_inf,f31,f32,dz,res_tol,max_iter,search_val)
	call generate_uniform_grid(dz,z1,z2,gridpts_z,nz)
	call bisection_and_solver(f31,f32,f2_inf,res_tol,max_iter,nz,dz,U,f1,f2)
	call search(gridpts_z,z_val,U,search_val)
	
		open(unit=2,file='BL_Solution.dat')
		
		write(2,*),"coordinate	","f	","f'	","f''	"
		
		
		do k=0,nz
			write(2,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(k),U(1:3,k)
		end do
		close(2)

end program bnd_layer
