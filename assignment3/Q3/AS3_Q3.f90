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


	subroutine read_input(z1,z2,t_final,dz,np,dt,mu,res_tol)
	  implicit none
	  real*8 		:: z1,z2,t_final,res_tol
	  integer		:: np
	  real*8,allocatable	:: dt(:),dz(:),mu(:)
	  
	  namelist/input_parameters/z1,z2,t_final,np,res_tol
	  namelist/list/dt,dz,mu
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  allocate(dt(1:np),dz(1:np),mu(1:np))
	  read(1, nml=list)
	  close(1)
	end subroutine read_input
	
	
	subroutine initial_cond(u,nz,gridpts_z)
		implicit none
		integer, intent(in)	:: nz
		real*8,intent(inout)	:: u(0:nz)
		real*8, intent(in)	:: gridpts_z(0:nz)
		real*8, parameter	:: pi = 4.0d0*atan(1.0d0)
		
		u(0:nz) =0.0d0 			!zero initial condition
		
	end subroutine initial_cond
	
	
	subroutine steady_state(u,u_prev,res_tol,ss)
		implicit none
		real*8, intent(in)	:: u(:),u_prev(:)
		real*8, intent(in)	:: res_tol
		logical, intent(out)	:: ss
		ss = .false.
		if(norm2(u-u_prev)/norm2(u_prev)<res_tol) ss = .true.
	end subroutine steady_state
		
		
end module input_grid_and_boundary


module math
	
	abstract interface
		function advection(x) result(res)
		    real*8, intent(in)	:: x
		    real*8		:: res
		end function advection
    	end interface
	
	contains
	
	
	real*8 function F(u)			!Function F, conservative form of advection term
		implicit none
		real*8, intent(in) :: u
		F = u*u/2.0d0
	end function F
	
	real*8 function source(z)
		implicit none
		real*8, intent(in) 	:: z
		real*8, parameter	:: pi = 4.0d0*atan(1.0d0)
		source = 4.0d0*sin(2.0d0*pi*z)	!𝑓 = 4sin(2𝜋𝑥)
	end function source
		
	integer function indx(i,nb)
		implicit none
		integer,intent(in) :: nb
		integer,intent(in) :: i

		if(i < 0) then
			indx = nb + i
		elseif(i > nb) then 
			indx = i - nb            
		else  				! This handles 0 <= indx <= nb
			indx = i
		endif
  	end function indx  
end module math	


module burgers_eqn_schemes
	contains
	
	subroutine FTCS(U,F,mu,dz,dt,nz,nt,gridpts_z,source,res_tol,ssc)
		use input_grid_and_boundary
		use math, only: indx, advection
		implicit none
		integer, intent(in)			:: nz,nt
		integer, intent(out)			:: ssc
		real*8, intent(in)			:: dz,dt,mu,res_tol
		real*8, allocatable, intent(in)		:: gridpts_z(:)
		real*8, allocatable, intent(out)	:: U(:,:)
		procedure(advection)			:: F,source
		integer					:: j,n
		real*8					:: c,r
		logical					:: ss
		c=dt/2.0d0/dz; r=mu*dt/dz/dz
		
		allocate(U(0:nz,0:nt))
		call initial_cond(U(0:nz,0),nz,gridpts_z)
		do n=1,nt
			do j=0,nz
				U(j,n)=U(j,n-1) - c*(F(U(indx(j+1,nz),n-1)) - F(U(indx(j-1,nz),n-1))) + r*(U(indx(j+1,nz),n-1) - 2.0d0*U(indx(j,nz),n-1) + U(indx(j-1,nz),n-1)) + source(gridpts_z(j))*dt
			end do
			call steady_state(U(0:nz,n),U(0:nz,n-1),res_tol,ss)
			ssc = n
			if(ss) exit
		end do
		print*,"FTCS solution done"
	end subroutine FTCS
	
	
	subroutine MacCormack(U,F,mu,dz,dt,nz,nt,gridpts_z,source,res_tol,ssc)
		use input_grid_and_boundary
		use math, only: indx, advection
		implicit none
		integer, intent(in)			:: nz,nt
		integer, intent(out)			:: ssc
		real*8, intent(in)			:: dz,dt,mu,res_tol
		real*8, allocatable, intent(in)		:: gridpts_z(:)
		real*8, allocatable, intent(out)	:: U(:,:)
		real*8, allocatable			:: pred_U(:)
		procedure(advection)			:: F,source
		integer					:: j,n
		real*8					:: c,r
		logical					:: ss
		c=dt/dz; r=mu*dt/dz/dz
		
		allocate(U(0:nz,0:nt),pred_U(0:nz))
		call initial_cond(U(0:nz,0),nz,gridpts_z)
		do n=1,nt
			do j=0,nz
				pred_U(j)=U(j,n-1) - c*(F(U(indx(j+1,nz),n-1)) - F(U(indx(j,nz),n-1))) + r*(U(indx(j+1,nz),n-1) - 2.0d0*U(indx(j,nz),n-1) + U(indx(j-1,nz),n-1)) + source(gridpts_z(j))*dt
			end do
			do j=0,nz
				U(j,n)=(U(j,n-1) + pred_U(j) - c*(F(pred_U(indx(j,nz))) - F(pred_U(indx(j-1,nz)))) + r*(pred_U(indx(j+1,nz)) - 2.0d0*pred_U(indx(j,nz)) + pred_U(indx(j-1,nz))) + source(gridpts_z(j))*dt)*0.5d0
			end do
			call steady_state(U(0:nz,n),U(0:nz,n-1),res_tol,ss)
			ssc = n
			if(ss) exit
		end do
		print*,"MacCormack solution done"
	end subroutine MacCormack
	
	
end module burgers_eqn_schemes


program Burgers
	use input_grid_and_boundary
	use burgers_eqn_schemes
	use math
	implicit none
	real*8 			:: z1,z2,t_final,res_tol
	integer			:: nz,nt,np,j,k,ssc1,ssc2
	real*8, allocatable	:: FTCS_sol(:,:),McK_sol(:,:),gridpts_z(:),gridpts_t(:),dt(:),dz(:),mu(:)
	character(len=10)	:: ssc_str,pr
	
	call read_input(z1,z2,t_final,dz,np,dt,mu,res_tol)
	do j=1,np
		call generate_uniform_grid(dz(j),dt(j),z1,z2,t_final,gridpts_z,gridpts_t,nz,nt)
		call FTCS(FTCS_sol,F,mu(j),dz(j),dt(j),nz,nt,gridpts_z,source,res_tol,ssc1)
		call MacCormack(McK_sol,F,mu(j),dz(j),dt(j),nz,nt,gridpts_z,source,res_tol,ssc2)
		write(pr,'(I3)') j
		write(ssc_str,'(I6)') ssc1
		open(unit=2,file='Pr_'//trim(adjustl(pr))//'_FTCS_Solution_ssc_'//trim(adjustl(ssc_str))//'.dat')
		write(ssc_str,'(I6)') ssc2
		open(unit=3,file='Pr_'//trim(adjustl(pr))//'_MacCormack_Solution_ssc_'//trim(adjustl(ssc_str))//'.dat')
		
		print*, "steady-state-time:	",gridpts_t(ssc1),gridpts_t(ssc2)
		write(2,*),"coordinate	","solution, t=0 --> ssc"
		write(3,*),"coordinate	","solution, t=0 --> ssc"
		
		do k=0,nz
			write(2,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(k),FTCS_sol(k,0:ssc1)
			write(3,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(k),McK_sol(k,0:ssc2)
		end do
		close(2);close(3)
		deallocate(FTCS_sol,McK_sol)
	end do

end program Burgers
