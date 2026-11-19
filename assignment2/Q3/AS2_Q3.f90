module input_grid_and_boundary
	contains
 	subroutine generate_uniform_grid(dz,dt,z1,z2,t_ss,gridpts_z,gridpts_t,nz,nt)
 	implicit none
	 INTEGER 					:: i
	 integer, intent(out)				:: nz,nt
	 REAL*8,intent(in) 				:: z1, z2,t_ss
	 REAL*8, intent(inout)				:: dz,dt
	 REAL(kind=8),allocatable, intent(out) 		:: gridpts_z(:), gridpts_t(:)
	 
	 nz = nint((z2 - z1)/dz)
	 nt = nint(t_ss/dt)
	 allocate(gridpts_z(0:nz),gridpts_t(0:nt))
	 !dz = (z2-z1)/nz
	 !dt = t_ss/nt
	 do i=0,nz
	   gridpts_z(i)= z1 + dz*i  
	 end do
	
	 do i=0,nt
	 	gridpts_t(i)= dt*i
	 end do
	 
	 print*, "Grid_generation done"
	end subroutine generate_uniform_grid


	subroutine read_input(L,t_final,n_size,dz_size,Co_size,n,dz,Co)
	  implicit none
	  real*8 	:: L,t_final
	  integer	:: n_size,dz_size,Co_size
	  integer, allocatable	:: n(:)
	  real*8, allocatable	:: dz(:),Co(:)
	 
	  namelist/input_parameters/L,t_final,n_size,dz_size,Co_size
	  namelist/arrays/n,dz,Co
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  allocate(n(n_size),dz(dz_size),Co(Co_size))
	  read(1, nml=arrays)
	  close(1)
	end subroutine read_input
	
	
	subroutine initial_cond(nz,L,gridpts_z,n,sinus,square)
		implicit none
		real*8,allocatable,intent(in)	::	gridpts_z(:)
		integer, intent(in)		::	nz,n
		real*8,allocatable, intent(out)	::	sinus(:),square(:)
		real*8, parameter		::	pi=4*atan(1.0d0)
		real*8, intent(in)		::	L
		integer				::	j
		
		allocate(sinus(0:nz),square(0:nz))
		do j=0,nz
			sinus(j) = sin(2.0d0*n*pi*gridpts_z(j)/L)
			if(gridpts_z(j)<=10.0d0) then
				square(j) = 1.0d0
			else
				square(j) = 0.0d0
			end if
		end do
		square(nz) = square(0)
		
	end subroutine initial_cond
		
		
end module input_grid_and_boundary


module wave_eqn_schemes
	contains
	subroutine fo_upwind(nz,nt,Co,sinus,square,sinus_sol,square_sol)
		implicit none
		integer, intent(in)		::	nz,nt
		real*8, intent(in)		::	Co
		real*8, allocatable, intent(in)	::	sinus(:),square(:)
		integer				::	j,n
		real*8, allocatable,intent(out)	::	sinus_sol(:,:),square_sol(:,:)
		
		allocate(sinus_sol(0:nz,0:nt),square_sol(0:nz,0:nt))
		sinus_sol(:,0) = sinus; square_sol(:,0) = square
		
		do n=1,nt
			do j=1,nz
				sinus_sol(j,n) = sinus_sol(j,n-1) - Co*(sinus_sol(j,n-1)-sinus_sol(j-1,n-1))
				square_sol(j,n) = square_sol(j,n-1) - Co*(square_sol(j,n-1)-square_sol(j-1,n-1))
			end do
			sinus_sol(0,n) = sinus_sol(nz,n)
			square_sol(0,n) = square_sol(nz,n)
		end do
		
		print*,"fo_upwind solution done"
	end subroutine fo_upwind
	
	
	subroutine lax(nz,nt,Co,sinus,square,sinus_sol,square_sol)
		implicit none
		integer, intent(in)		::	nz,nt
		real*8, intent(in)		::	Co
		real*8, allocatable, intent(in)	::	sinus(:),square(:)
		integer				::	j,n
		real*8, allocatable,intent(out)	::	sinus_sol(:,:),square_sol(:,:)
		
		allocate(sinus_sol(0:nz,0:nt),square_sol(0:nz,0:nt))
		sinus_sol(:,0) = sinus; square_sol(:,0) = square
		
		do n=1,nt
			do j=1,nz-1
				sinus_sol(j,n) = (1.0d0+Co)*sinus_sol(j-1,n-1)/2.0d0 + (1.0d0-Co)*sinus_sol(j+1,n-1)/2.0d0
				square_sol(j,n) = (1.0d0+Co)*square_sol(j-1,n-1)/2.0d0 + (1.0d0-Co)*square_sol(j+1,n-1)/2.0d0
			end do
			sinus_sol(nz,n) = (1.0d0+Co)*sinus_sol(nz-1,n-1)/2.0d0 + (1.0d0-Co)*sinus_sol(1,n-1)/2.0d0
			square_sol(nz,n) = (1.0d0+Co)*square_sol(nz-1,n-1)/2.0d0 + (1.0d0-Co)*square_sol(1,n-1)/2.0d0
			sinus_sol(0,n) = sinus_sol(nz,n)
			square_sol(0,n) = square_sol(nz,n)
		end do
		
		print*,"lax solution done"
	end subroutine lax
	
	
	subroutine MacCormack(nz,nt,Co,sinus,square,sinus_sol,square_sol)
		implicit none
		integer, intent(in)		::	nz,nt
		real*8, intent(in)		::	Co
		real*8, allocatable, intent(in)	::	sinus(:),square(:)
		integer				::	j,n
		real*8, allocatable,intent(out)	::	sinus_sol(:,:),square_sol(:,:)
		real*8, allocatable		::	sinus_pred(:),square_pred(:)
		
		allocate(sinus_sol(0:nz,0:nt),square_sol(0:nz,0:nt),sinus_pred(0:nz),square_pred(0:nz))
		sinus_sol(:,0) = sinus; square_sol(:,0) = square
		
		do n=1,nt
			do j=0,nz-1	!predictor
				sinus_pred(j) = sinus_sol(j,n-1) - Co*(sinus_sol(j+1,n-1) - sinus_sol(j,n-1))
				square_pred(j) = square_sol(j,n-1) - Co*(square_sol(j+1,n-1) - square_sol(j,n-1))
			end do
			sinus_pred(nz) = sinus_sol(nz,n-1) - Co*(sinus_sol(1,n-1) - sinus_sol(nz,n-1))
			square_pred(nz) = square_sol(nz,n-1) - Co*(square_sol(1,n-1) - square_sol(nz,n-1))
			
			do j=1,nz	!corrector
				sinus_sol(j,n) = 0.5d0*(sinus_sol(j,n-1) + sinus_pred(j) - Co*(sinus_pred(j) - sinus_pred(j-1)))
				square_sol(j,n) = 0.5d0*(square_sol(j,n-1) + square_pred(j) - Co*(square_pred(j) - square_pred(j-1)))
			end do
			sinus_sol(0,n) = sinus_sol(nz,n)
			square_sol(0,n) = square_sol(nz,n)
		end do
		
		print*,"MacCormack solution done"
	end subroutine MacCormack
	
	
	subroutine so_upwind(nz,nt,Co,sinus,square,sinus_sol,square_sol)
		implicit none
		integer, intent(in)		::	nz,nt
		real*8, intent(in)		::	Co
		real*8, allocatable, intent(in)	::	sinus(:),square(:)
		integer				::	j,n
		real*8, allocatable,intent(out)	::	sinus_sol(:,:),square_sol(:,:)
		
		allocate(sinus_sol(0:nz,0:nt),square_sol(0:nz,0:nt))
		sinus_sol(:,0) = sinus; square_sol(:,0) = square
		
		do n=1,nt
			do j=2,nz
				sinus_sol(j,n) = sinus_sol(j,n-1) - Co*(sinus_sol(j,n-1) - sinus_sol(j-1,n-1)) + CO*(Co-1.0d0)*(sinus_sol(j,n-1)-2*sinus_sol(j-1,n-1)+sinus_sol(j-2,n-1))/2.0d0
				square_sol(j,n) = square_sol(j,n-1) - Co*(square_sol(j,n-1) - square_sol(j-1,n-1)) + CO*(Co-1.0d0)*(square_sol(j,n-1)-2*square_sol(j-1,n-1)+square_sol(j-2,n-1))/2.0d0
			end do
			sinus_sol(0,n) = sinus_sol(nz,n)
			square_sol(0,n) = square_sol(nz,n)
			sinus_sol(1,n) = sinus_sol(1,n-1) - Co*(sinus_sol(1,n-1) - sinus_sol(0,n-1)) + CO*(Co-1.0d0)*(sinus_sol(1,n-1)-2*sinus_sol(0,n-1)+sinus_sol(nz-1,n-1))/2.0d0
			square_sol(1,n) = square_sol(1,n-1) - Co*(square_sol(1,n-1) - square_sol(0,n-1)) + CO*(Co-1.0d0)*(square_sol(1,n-1)-2*square_sol(0,n-1)+square_sol(nz-1,n-1))/2.0d0
		end do
		
		print*,"so_upwind solution done"
	end subroutine so_upwind
	
end module wave_eqn_schemes


program wave
	use input_grid_and_boundary
	use wave_eqn_schemes
	implicit none
	real*8 			:: L,t_final,dt,z1,z2
	integer			:: n_size,dz_size,Co_size,j,k,i,nz,nt,p,q
	integer, allocatable	:: n(:)
	real*8, allocatable	:: dz(:),Co(:),sinus(:),square(:),gridpts_z(:),gridpts_t(:)
	real*8, allocatable	:: sinus_sol1(:,:),square_sol1(:,:),sinus_sol2(:,:),square_sol2(:,:),sinus_sol3(:,:),square_sol3(:,:),sinus_sol4(:,:),square_sol4(:,:)
	character(len=10)	:: n_str,dz_str,Co_str
	
	
	call read_input(L,t_final,n_size,dz_size,Co_size,n,dz,Co)
	z1 = 0.0d0
	z2 = Z1+L
	
	do j=1,dz_size
		write(dz_str,'(F3.1)') dz(j)
		do k=1,Co_size
			write(Co_str,'(F3.1)') Co(k)
			dt=Co(k)*dz(j)
			call generate_uniform_grid(dz(j),dt,z1,z2,t_final,gridpts_z,gridpts_t,nz,nt)
			do i=1,n_size
				call initial_cond(nz,L,gridpts_z,n(i),sinus,square)
				call fo_upwind(nz,nt,Co(k),sinus,square,sinus_sol1,square_sol1)
				call lax(nz,nt,Co(k),sinus,square,sinus_sol2,square_sol2)
				call MacCormack(nz,nt,Co(k),sinus,square,sinus_sol3,square_sol3)
				call so_upwind(nz,nt,Co(k),sinus,square,sinus_sol4,square_sol4)
				print*,n(i)
				write(n_str,'(I2)') n(i)
				open(unit=2,file='sol_dz_'//trim(adjustl(dz_str))//'_C0_'//trim(adjustl(Co_str))//'_n_'//trim(adjustl(n_str))//'.dat')
				write(2,'(1X,A18,*(6X,A18))'),"time_t","position_z","fo_upwind","lax","MacCormack","so_upwind"
				do p = 0,nt
					do q=0,nz
						write(2,'(1X, F8.4, 2X, F8.4, *(2X, F15.6))'),gridpts_t(p),gridpts_z(q),sinus_sol1(q,p),square_sol1(q,p),sinus_sol2(q,p),square_sol2(q,p),sinus_sol3(q,p),square_sol3(q,p),sinus_sol4(q,p),square_sol4(q,p)
					end do
				end do
				close(2)
			end do
		end do
	end do

end program wave
