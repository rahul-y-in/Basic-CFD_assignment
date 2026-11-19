module input_grid_and_boundary
	contains
 	subroutine generate_uniform_grid(dz,dt,z1,z2,t_final,gridpts_z,gridpts_t,nz,nt)
 	implicit none
	 INTEGER 					:: i
	 integer, intent(out)				:: nz,nt
	 REAL*8,intent(in) 				:: z1, z2,t_final
	 REAL*8, intent(in)				:: dz,dt
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:), gridpts_t(:)
	 
	 nz = nint((z2 - z1)/dz)
	 nt = nint(t_final/dt)
	 allocate(gridpts_z(0:nz),gridpts_t(0:nt))

	 do i=0,nz
	   gridpts_z(i)= z1 + dz*i  
	 end do
	
	 do i=0,nt
	 	gridpts_t(i)= dt*i
	 end do
	 
	 print*, "Grid_generation done: z(nz)= ",gridpts_z(nz)," t(nt)= ",gridpts_t(nt)
	end subroutine generate_uniform_grid


	subroutine read_input(z1,z2,t_final,dz,dt,nu,u_inf,t)
	  implicit none
	  real*8 		:: z1,z2,t_final,dz,dt,nu,u_inf,t
	  
	  namelist/input_parameters/z1,z2,t_final,dz,dt,nu,u_inf,t
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input
	
	
	subroutine initial_cond(u,v,nz,u_inf)
		implicit none
		integer, intent(in)	:: nz
		real*8,intent(inout)	:: u(0:nz),v(0:nz)
		real*8, intent(in)	:: u_inf
		u(0)=0.0d0
		u(1:nz) =u_inf 	
		v(0:nz)=0.0d0
	end subroutine initial_cond
	
	
	subroutine b_cond(u,v,nz,u_inf)
		implicit none
		integer, intent(in)	:: nz
		real*8,intent(inout)	:: u(0:nz),v(0:nz)
		real*8, intent(in)	:: u_inf
		u(0)=0.0d0; v(0)=0.0d0
		u(nz) =u_inf 			
	end subroutine b_cond
		
		
end module input_grid_and_boundary


module bnd_layer_solvers
	contains
	
	subroutine FTCS(U,V,nu,dz,dt,nz,n)
		implicit none
		integer, intent(in)			:: nz,n
		real*8, intent(in)			:: dz,dt,nu
		real*8, allocatable, intent(inout)	:: U(:,:),V(:,:)
		integer					:: j
		real*8					:: c,r

		c=dt/2.0d0/dz; r=nu*dt/dz/dz
		
		do j=1,nz-1
			U(j,n)=U(j,n-1) - c*(V(j,n-1)/U(j,n-1))*(U(j+1,n-1)-U(j-1,n-1)) + r*(U(j+1,n-1)/U(j,n-1) - 2.0d0 + U(j-1,n-1)/U(j,n-1))
		end do
	end subroutine FTCS
	
	
	subroutine solver(U,V,nz,nt,dz,dt,u_inf,nu)
		use input_grid_and_boundary
		implicit none
		integer, intent(in)			:: nz,nt
		real*8, intent(in)			:: dz,dt,u_inf,nu
		real*8, allocatable, intent(out)	:: U(:,:),V(:,:)
		integer					:: j,n
		real*8					:: c,r
		allocate(U(0:nz,0:nt),V(0:nz,0:nt))
		call initial_cond(U(0:nz,0),V(0:nz,0),nz,u_inf)
		do n=1,nt
			call b_cond(U(0:nz,n),V(0:nz,n),nz,u_inf)
			call FTCS(U,V,nu,dz,dt,nz,n)
			V(1:nz,n)=V(0:nz-1,n) + (dz/dt)*(U(1:nz,n)-U(1:nz,n-1))
		end do
	print*,"solution done"
	end subroutine solver
	
	
	real*8 function delta(t,dt,gridpts_z,u_inf,U)				!Boundary layer thickness
		real*8, intent(in)			:: t,dt,u_inf
		real*8, allocatable, intent(in)		:: gridpts_z(:),U(:,:)
		integer					:: j,n
		n=nint(t/dt)
		do j=1,size(gridpts_z)-1
			if(U(j,n)>=0.99d0*u_inf) then
				delta = gridpts_z(j)
				exit
			elseif(j==size(gridpts_z)-1) then
				print*,"not found"
			end if
		end do
	end function delta
	
	
end module bnd_layer_solvers


program bnd_layer
	use input_grid_and_boundary
	use bnd_layer_solvers
	implicit none
	real*8 			:: dz,dt,z1,z2,t_final,u_inf,nu,t
	integer			:: nz,nt,np,j,k
	real*8, allocatable	:: U(:,:),V(:,:),gridpts_z(:),gridpts_t(:)
	
	call read_input(z1,z2,t_final,dz,dt,nu,u_inf,t)
	call generate_uniform_grid(dz,dt,z1,z2,t_final,gridpts_z,gridpts_t,nz,nt)
	call solver(U,V,nz,nt,dz,dt,u_inf,nu)
	
		open(unit=2,file='BL_Solution_U.dat')
		open(unit=3,file='BL_Solution_U_non-dimensional.dat')
		open(unit=4,file='BL_Solution_V.dat')
		
		write(2,*),"y-coordinate	","solution, x=0 --> 1"
		write(4,*),"y-coordinate	","solution, x=0 --> 1"
		write(3,*),"y-coordinate	","solution, x=0 --> 1"
		
		do k=0,nz
			write(2,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(k),U(k,0:nt)
			write(4,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(k),V(k,0:nt)
			write(3,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(k)/delta(t,dt,gridpts_z,u_inf,U),U(k,0:nt)/u_inf
		end do
		close(2);close(3);close(4)

end program bnd_layer
