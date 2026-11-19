module input_grid_and_boundary
	contains
 	subroutine generate_uniform_grid(x0,x1,np,gridpts_z)
 	implicit none
	 INTEGER, intent(in)				:: np
	 REAL*8,intent(in) 				:: x0,x1
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:)
	 integer					:: nd,j
	 real*8						:: dz
	 nd = np-1
	 dz = (x1 - x0)/nd
	 allocate(gridpts_z(0:nd))
	 
	 do j=0,nd
	   gridpts_z(j)= x0 + dz*j  
	 end do
	 
	 print*, "Grid_generation done, z(0)= ", gridpts_z(0), "z(nd)= ", gridpts_z(nd)
	end subroutine generate_uniform_grid


	subroutine read_input(x0,x1,ux0,u0t,u1t,r,Co,np,res_tol,max_t)
	  implicit none
	  real*8 	:: x0,x1,ux0,u0t,u1t,r,Co,res_tol
	  integer	:: np,max_t
	  
	  namelist/input_parameters/x0,x1,ux0,u0t,u1t,r,Co,np,res_tol,max_t
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input
	
	
	subroutine initial_cond(u,ux0)
		implicit none
		real*8,intent(inout)	::u(:)
		real*8, intent(in)	::ux0
		u=ux0
	end subroutine initial_cond
		
		
end module input_grid_and_boundary


module burgers_eqn_schemes
	contains
	subroutine steady_state(u,u_prev,res_tol,ss)
		implicit none
		real*8, intent(in)	:: u(:),u_prev(:)
		real*8, intent(in)	:: res_tol
		logical, intent(out)	:: ss
		ss = .false.
		if(norm2(u-u_prev)/norm2(u_prev)<res_tol) ss = .true.
	end subroutine steady_state
	
	
	subroutine FTCS(U,Co,r,nd,max_t,u0t,u1t,res_tol,ux0,ssc)
		use input_grid_and_boundary
		implicit none
		integer, intent(in)			:: nd,max_t
		real*8, intent(in)			:: Co,r,u0t,u1t,res_tol,ux0
		real*8, allocatable, intent(out)	:: U(:,:)
		integer					:: j,n
		logical					:: ss		!steady state
		integer, intent(out)			:: ssc		!steady state count i.e. time-steps taken till ss
		character(len=6)			:: ssc_str
		
		allocate(U(0:nd,0:max_t))
		U(0,0)=u0t; U(nd,0)=u1t
		call initial_cond(U(1:nd-1,0),ux0)
		do n=1,max_t
			U(0,n)=u0t; U(nd,n)=u1t
			U(1:nd-1,n)=(r+Co/2.0d0)*U(0:nd-2,n-1)+(1.0d0-2.0d0*r)*U(1:nd-1,n-1)+(r-Co/2.0d0)*U(2:nd,n-1)
			call steady_state(U(1:nd-1,n),U(1:nd-1,n-1),res_tol,ss)
			ssc = n
			if(ss) exit
		end do
		write(ssc_str,'(I4)') ssc
		print*,"FTCS solution done: "//trim(adjustl(ssc_str))//" time-steps"
	end subroutine FTCS
	
	
	subroutine MacCormack(U,Co,r,nd,max_t,u0t,u1t,res_tol,ux0,ssc)
		use input_grid_and_boundary
		implicit none
		integer, intent(in)			:: nd,max_t
		real*8, intent(in)			:: Co,r,u0t,u1t,res_tol,ux0
		real*8, allocatable, intent(out)	:: U(:,:)
		real*8, allocatable			:: temp_U(:)
		integer					:: j,n
		logical					:: ss		!steady state
		integer, intent(out)			:: ssc		!steady state count i.e. time-steps taken till ss
		character(len=6)			:: ssc_str
		
		allocate(U(0:nd,0:max_t),temp_U(0:nd))
		U(0,0)=u0t; U(nd,0)=u1t
		call initial_cond(U(1:nd-1,0),ux0)
		do n=1,max_t
			temp_U(0)=u0t; temp_U(nd)=u1t
			temp_U(1:nd-1)=r*U(0:nd-2,n-1)+(1.0d0-2.0d0*r+Co)*U(1:nd-1,n-1)+(r-Co)*U(2:nd,n-1)					!predictor
			U(0,n)=u0t; U(nd,n)=u1t
			U(1:nd-1,n)=0.5d0*(U(1:nd-1,n-1) + (r+Co)*temp_U(0:nd-2) + (1.0d0-2.0d0*r-Co)*temp_U(1:nd-1) + r*temp_U(2:nd))		!corrector
			call steady_state(U(1:nd-1,n),U(1:nd-1,n-1),res_tol,ss)
			ssc = n
			if(ss) exit
		end do
		write(ssc_str,'(I4)') ssc
		print*,"MacCormack solution done: "//trim(adjustl(ssc_str))//" time-steps"
	end subroutine MacCormack
	
	
end module burgers_eqn_schemes


program Burgers
	use input_grid_and_boundary
	use burgers_eqn_schemes
	implicit none
	real*8 			:: x0,x1,ux0,u0t,u1t,r,Co,res_tol
	integer			:: np,max_t,ssc1,ssc2,j
	real*8, allocatable	:: FTCS_sol(:,:),McK_sol(:,:),gridpts_z(:)
	character(len=6)	:: ssc_str
	
	call read_input(x0,x1,ux0,u0t,u1t,r,Co,np,res_tol,max_t)
	call generate_uniform_grid(x0,x1,np,gridpts_z)
	call FTCS(FTCS_sol,Co,r,np-1,max_t,u0t,u1t,res_tol,ux0,ssc1)
	call MacCormack(McK_sol,Co,r,np-1,max_t,u0t,u1t,res_tol,ux0,ssc2)
	
	write(ssc_str,'(I4)') ssc1
	open(unit=2,file='FTCS_Solution_'//trim(adjustl(ssc_str))//'.dat')
	write(ssc_str,'(I4)') ssc2
	open(unit=3,file='MacCormack_Solution_'//trim(adjustl(ssc_str))//'.dat')
	
	write(2,*),"coordinate	","solution, t=0 --> ssc"
	write(3,*),"coordinate	","solution, t=0 --> ssc"

	do j=0,np-1
		write(2,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(j),FTCS_sol(j,0:ssc1)
		write(3,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(j),McK_sol(j,0:ssc2)
	end do
	close(2)
	close(3)
end program Burgers
