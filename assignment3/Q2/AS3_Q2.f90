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


	subroutine read_input(z1,z2,t_final,dz,ndt,dt,mu)
	  implicit none
	  real*8 		:: z1,z2,t_final,dz,mu
	  integer		:: ndt
	  real*8,allocatable	:: dt(:)
	  
	  namelist/input_parameters/z1,z2,t_final,dz,ndt,mu
	  namelist/list/dt
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  allocate(dt(1:ndt))
	  read(1, nml=list)
	  close(1)
	end subroutine read_input
	
	
	subroutine initial_cond(u,nz,gridpts_z)
		implicit none
		integer, intent(in)	:: nz
		real*8,intent(inout)	:: u(0:nz)
		real*8, intent(in)	:: gridpts_z(0:nz)
		real*8, parameter	:: pi = 4.0d0*atan(1.0d0)
		
		u(0:nz) = 1.0d0 + sin(4.0d0*pi*gridpts_z(0:nz))  		!𝑢 = 1 + sin(4𝜋𝑥)
		
	end subroutine initial_cond
		
		
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
	
	subroutine FTCS(U,F,mu,dz,dt,nz,nt,gridpts_z)
		use input_grid_and_boundary
		use math, only: indx, advection
		implicit none
		integer, intent(in)			:: nz,nt
		real*8, intent(in)			:: dz,dt,mu
		real*8, allocatable, intent(in)		:: gridpts_z(:)
		real*8, allocatable, intent(out)	:: U(:,:)
		procedure(advection)			:: F
		integer					:: j,n
		real*8					:: c,r
		
		c=dt/2.0d0/dz; r=mu*dt/dz/dz
		
		allocate(U(0:nz,0:nt))
		call initial_cond(U(0:nz,0),nz,gridpts_z)
		do n=1,nt
			do j=0,nz
				U(j,n)=U(j,n-1) - c*(F(U(indx(j+1,nz),n-1)) - F(U(indx(j-1,nz),n-1))) + r*(U(indx(j+1,nz),n-1) - 2.0d0*U(indx(j,nz),n-1) + U(indx(j-1,nz),n-1))
			end do
		end do
		print*,"FTCS solution done"
	end subroutine FTCS
	
	
	subroutine MacCormack(U,F,mu,dz,dt,nz,nt,gridpts_z)
		use input_grid_and_boundary
		use math, only: indx, advection
		implicit none
		integer, intent(in)			:: nz,nt
		real*8, intent(in)			:: dz,dt,mu
		real*8, allocatable, intent(in)		:: gridpts_z(:)
		real*8, allocatable, intent(out)	:: U(:,:)
		real*8, allocatable			:: pred_U(:)
		procedure(advection)			:: F
		integer					:: j,n
		real*8					:: c,r
		
		c=dt/dz; r=mu*dt/dz/dz
		
		allocate(U(0:nz,0:nt),pred_U(0:nz))
		call initial_cond(U(0:nz,0),nz,gridpts_z)
		do n=1,nt
			do j=0,nz
				pred_U(j)=U(j,n-1) - c*(F(U(indx(j+1,nz),n-1)) - F(U(indx(j,nz),n-1))) + r*(U(indx(j+1,nz),n-1) - 2.0d0*U(indx(j,nz),n-1) + U(indx(j-1,nz),n-1))
			end do
			do j=0,nz
				U(j,n)=(U(j,n-1) + pred_U(j) - c*(F(pred_U(indx(j,nz))) - F(pred_U(indx(j-1,nz)))) + r*(pred_U(indx(j+1,nz)) - 2.0d0*pred_U(indx(j,nz)) + pred_U(indx(j-1,nz))))*0.5d0
			end do
		end do
		print*,"MacCormack solution done"
	end subroutine MacCormack
	
	
	subroutine FTUWCS(U,F,mu,dz,dt,nz,nt,gridpts_z)
		use input_grid_and_boundary
		use math, only: indx, advection
		implicit none
		integer, intent(in)			:: nz,nt
		real*8, intent(in)			:: dz,dt,mu
		real*8, allocatable, intent(in)		:: gridpts_z(:)
		real*8, allocatable, intent(out)	:: U(:,:)
		procedure(advection)			:: F
		integer					:: j,n
		real*8					:: c,r
		
		c=dt/dz/2.0d0; r=mu*dt/dz/dz
		
		allocate(U(0:nz,0:nt))
		call initial_cond(U(0:nz,0),nz,gridpts_z)
		do n=1,nt
			do j=0,nz
				U(j,n)=U(j,n-1) - c*((U(j,n-1)+abs(U(j,n-1)))*(U(indx(j,nz),n-1) - U(indx(j-1,nz),n-1)) + (U(j,n-1)-abs(U(j,n-1)))*(U(indx(j+1,nz),n-1) - U(indx(j,nz),n-1))) + r*(U(indx(j+1,nz),n-1) - 2.0d0*U(indx(j,nz),n-1) + U(indx(j-1,nz),n-1))
			end do
		end do
		print*,"FTUWCS solution done"
	end subroutine FTUWCS
	
	
end module burgers_eqn_schemes


program Burgers
	use input_grid_and_boundary
	use burgers_eqn_schemes
	use math
	implicit none
	real*8 			:: dz,z1,z2,t_final,mu
	integer			:: nz,nt,ndt,j,k
	real*8, allocatable	:: FTCS_sol(:,:),McK_sol(:,:),FTUWCS_sol(:,:),gridpts_z(:),gridpts_t(:),dt(:)
	character(len=10)	:: ssc_str
	
	call read_input(z1,z2,t_final,dz,ndt,dt,mu)
	do j=1,size(dt)
		call generate_uniform_grid(dz,dt(j),z1,z2,t_final,gridpts_z,gridpts_t,nz,nt)
		call FTCS(FTCS_sol,F,mu,dz,dt(j),nz,nt,gridpts_z)
		call MacCormack(McK_sol,F,mu,dz,dt(j),nz,nt,gridpts_z)
		call FTUWCS(FTUWCS_sol,F,mu,dz,dt(j),nz,nt,gridpts_z)

		write(ssc_str,'(E8.1)') dt(j)
		open(unit=2,file='FTCS_Solution_dt_'//trim(adjustl(ssc_str))//'.dat')
		open(unit=3,file='MacCormack_Solution_dt_'//trim(adjustl(ssc_str))//'.dat')
		open(unit=4,file='FTUWCS_Solution_dt_'//trim(adjustl(ssc_str))//'.dat')
		
		write(2,*),"coordinate	","solution, t=0 --> nt"
		write(3,*),"coordinate	","solution, t=0 --> nt"
		write(4,*),"coordinate	","solution, t=0 --> nt"

		do k=0,nz
			write(2,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(k),FTCS_sol(k,0:nt)
			write(3,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(k),McK_sol(k,0:nt)
			write(4,'(1X, F8.4, *(2X, E15.6))'),gridpts_z(k),FTUWCS_sol(k,0:nt)
		end do
		close(2);close(3);close(4)
	end do

end program Burgers
