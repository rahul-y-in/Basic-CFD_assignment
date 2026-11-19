module input_ibcs_and_fns
	contains
 	 subroutine generate_grid(ny_l,ny_r,nz_u,nz_l,Ly_l,Ly_r,Lz_u,Lz_l,gridpts_y,gridpts_z)
		 implicit none
		 INTEGER,intent(in) 				:: ny_l,ny_r,nz_u,nz_l
		 INTEGER 					:: i,j,k
		 REAL*8,intent(in) 				:: Ly_l,Ly_r,Lz_u,Lz_l
		 REAL(kind=8) 					:: dy_l,dy_r,dz_l,dz_u
		 REAL(kind=8),allocatable, intent(out) 		:: gridpts_y(:),gridpts_z(:)
		 allocate(gridpts_y(-ny_l:ny_r),gridpts_z(-nz_l:nz_u))
		 dy_l=Ly_l/ny_l
		 dy_r=Ly_r/ny_r
		 dz_l=Lz_l/nz_l
		 dz_u=Lz_u/nz_u
		 
		 do i=0,-ny_l,-1
		   gridpts_y(i)= dy_l*(i)
		 end do
		 do i=0,ny_r
		   gridpts_y(i)= dy_r*(i)
		 end do
		 do i=0,-nz_l,-1
		   gridpts_z(i)= dz_l*(i)  
		 end do
		 do i=0,nz_u
		   gridpts_z(i)= dz_u*(i)  
		 end do
		 print*, "Grid_generation done"
	end subroutine generate_grid


	integer*4 function jb(indx,ny_l,ny_r,nz_l,nz_u)
	      implicit none
	      integer*4,intent(in)::indx,ny_l,ny_r,nz_l,nz_u
	      if(indx<nz_u+1) then
		 jb=-ny_l 
	      elseif(indx<(nz_u+ny_l+ny_r+1)) then
		 jb=indx-ny_l-nz_u-1          
	      elseif(indx<(nz_l+2*nz_u+ny_l+ny_r+1)) then
		 jb=ny_r          
	      elseif(indx<(nz_l+2*nz_u+ny_l+2*ny_r+1)) then
		 jb=nz_l+2*nz_u+ny_l+2*ny_r+1-indx
	      elseif(indx<(2*nz_l+2*nz_u+ny_l+2*ny_r+1)) then
		 jb=0
	      else
		 jb=2*nz_l+2*nz_u+ny_l+2*ny_r+1-indx
	      endif
	    return
    	end function jb


	integer*4 function kb(indx,ny_l,ny_r,nz_l,nz_u)
	      implicit none
	      integer*4,intent(in)::indx,ny_l,ny_r,nz_l,nz_u
	      if(indx<nz_u+1) then
		 kb=indx
	      elseif(indx<(nz_u+ny_l+ny_r+1)) then
		 kb=nz_u
	      elseif(indx<(nz_l+2*nz_u+ny_l+ny_r+1)) then
		 kb=2*nz_u+ny_l+ny_r+1 - indx         
	      elseif(indx<(nz_l+2*nz_u+ny_l+2*ny_r+1)) then
		 kb=-nz_l
	      elseif(indx<(2*nz_l+2*nz_u+ny_l+2*ny_r+1)) then
		 kb=indx-(2*nz_l+2*nz_u+ny_l+2*ny_r+1)
	      else
		 kb=0
	      endif
	    return
    	end function kb


	real*8 function h(j1,k1,j2,k2,gridpts_y,gridpts_z)
		implicit none
	  real*8  ,allocatable, intent(in)              ::gridpts_y(:),gridpts_z(:)
	  integer, intent(in)				::j1,k1,j2,k2
		h = sqrt((gridpts_y(j1)-gridpts_y(j2))**2 + (gridpts_z(k1)-gridpts_z(k2))**2)
	end function h


	Subroutine apply_bcs(psi,omega,ny_l,ny_r,nz_l,nz_u,gridpts_y,gridpts_z,U0,u,v)
		implicit none
		integer, intent(in) 			::	ny_l,ny_r,nz_l,nz_u
		real*8, intent(in)			::	U0
		integer					::	i,j,k,nb
		real*8, allocatable, intent(inout)	::	psi(:,:),omega(:,:),u(:,:),v(:,:)
		real*8, allocatable, intent(in)		::	gridpts_y(:),gridpts_z(:)
		
		nb=2*(ny_l+ny_r+nz_l+nz_u)
		do i=1,nb
			j=jb(i,ny_l,ny_r,nz_l,nz_u)
			k=kb(i,ny_l,ny_r,nz_l,nz_u)
			psi(j,k)=0.0d0
			u(j,k)=0.0d0; v(j,k)=0.0d0
			
			if(i<=nz_u+1) then
				u(j,k)=U0
				psi(j,k)=gridpts_z(k)*U0
	      		elseif(i<=(nz_u+ny_l+ny_r+1)) then
		 		psi(j,k)=psi(-ny_l,k)
	      		elseif(i<(nz_l+2*nz_u+ny_l+ny_r+1)) then
	      			psi(j,k)=((gridpts_y(j)-gridpts_y(j-1))*psi(j-2,k) - (gridpts_y(j)-gridpts_y(j-2))*psi(j-1,k))/(gridpts_y(j-2)-gridpts_y(j-1))
		 		u(j,k)=u(j-1,k)
		 		v(j,k)=v(j-1,k)
	      		end if
	      		if(i==1 .or. i==nz_u+1) u(j,k)=0.0d0
		end do
		do i=1,nb
			j=jb(i,ny_l,ny_r,nz_l,nz_u)
			k=kb(i,ny_l,ny_r,nz_l,nz_u)
	      		if(i<=nz_u+1) then
		 		omega(j,k)=2.0d0*(psi(j,k)-psi(j+1,k))/(h(j+1,k,j,k,gridpts_y,gridpts_z)**2)
	      		elseif(i<=(nz_u+ny_l+ny_r+1)) then
		 		omega(j,k)=2.0d0*(psi(j,k)-psi(j,k-1))/(h(j,k-1,j,k,gridpts_y,gridpts_z)**2)
	      		elseif(i<(nz_l+2*nz_u+ny_l+ny_r+1)) then
		 		omega(j,k)=2.0d0*(psi(j,k)-psi(j-1,k)+v(j,k)*h(j-1,k,j,k,gridpts_y,gridpts_z))/(h(j-1,k,j,k,gridpts_y,gridpts_z)**2)         
	      		elseif(i<=(nz_l+2*nz_u+ny_l+2*ny_r+1)) then
		 		omega(j,k)=2.0d0*(psi(j,k)-psi(j,k+1))/(h(j,k+1,j,k,gridpts_y,gridpts_z)**2)
	      		elseif(i<(2*nz_l+2*nz_u+ny_l+2*ny_r+1)) then
		 		omega(j,k)=2.0d0*(psi(j,k)-psi(j+1,k))/(h(j+1,k,j,k,gridpts_y,gridpts_z)**2)
	      		else
				omega(j,k)=2.0d0*(psi(j,k)-psi(j,k+1))/(h(j,k+1,j,k,gridpts_y,gridpts_z)**2)
	      		endif
		end do
	end subroutine apply_bcs
	

	Subroutine apply_ics(psi,omega,ny_l,ny_r,nz_l,nz_u,gridpts_y,gridpts_z,U0,u,v)
		implicit none
		integer, intent(in) 			::	ny_l,ny_r,nz_l,nz_u
		real*8, intent(in)			::	U0
		integer					::	i,j,k
		real*8, allocatable, intent(inout)	::	psi(:,:),omega(:,:),u(:,:),v(:,:)
		real*8, allocatable, intent(in)		::	gridpts_y(:),gridpts_z(:)
		
		psi=0.0d0
		omega=0.0d0
		u=0.0d0;v=0.0d0
		call apply_bcs(psi,omega,ny_l,ny_r,nz_l,nz_u,gridpts_y,gridpts_z,U0,u,v)
	end subroutine apply_ics


	subroutine read_input(ny_l,ny_r,nz_l,nz_u,Ly_l,Ly_r,Lz_u,Lz_l,U0,Re,iter1,iter2,dt,dt2,res_tol,D)
	  implicit none
	  real*8 	:: Ly_l,Ly_r,Lz_u,Lz_l,D,U0,dt,dt2,res_tol,Re
	  integer	:: ny_l,ny_r,nz_l,nz_u,iter1,iter2
	  
	  namelist/input_parameters/ny_l,ny_r,nz_l,nz_u,Ly_l,Ly_r,Lz_u,Lz_l,U0,Re,iter1,iter2,dt,dt2,res_tol,D
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input

end module input_ibcs_and_fns


module solvers
	contains
	subroutine time_marching(omega,u,v,Re,ny_l,ny_r,nz_l,nz_u,dt,gridpts_y,gridpts_z)	!essentially FTCS
		implicit none
		integer, intent(in)			:: ny_l,ny_r,nz_l,nz_u
		real*8, intent(in)			:: dt,Re
		real*8, allocatable			:: omega_out(:,:)
		real*8, allocatable, intent(inout)	:: omega(:,:)
		real*8, allocatable, intent(in)		:: u(:,:),v(:,:),gridpts_y(:),gridpts_z(:)
		integer					:: j,k
		real*8					:: diff,adv,ce,we,ea,no,so
		!allocate(omega_out(0:ny,0:nz))
		do k=-nz_l+1,nz_u-1
			do j=-ny_l+1,ny_r-1
				if(k<=0 .and. j<=0) then
				
				else
				ce=2.0d0/(gridpts_y(j+1)-gridpts_y(j))/(gridpts_y(j)-gridpts_y(j-1)) + 2.0d0/(gridpts_z(k+1)-gridpts_z(k))/(gridpts_z(k)-gridpts_z(k-1))
				we=2.0d0/(gridpts_y(j+1)-gridpts_y(j-1))/(gridpts_y(j)-gridpts_y(j-1))
				ea=2.0d0/(gridpts_y(j+1)-gridpts_y(j))/(gridpts_y(j+1)-gridpts_y(j-1))
				no=2.0d0/(gridpts_z(k+1)-gridpts_z(k))/(gridpts_z(k+1)-gridpts_z(k-1))
				so=2.0d0/(gridpts_z(k+1)-gridpts_z(k-1))/(gridpts_z(k)-gridpts_z(k-1))
				diff=(omega(j+1,k)*ea + omega(j-1,k)*we + omega(j,k+1)*no + omega(j,k-1)*so - omega(j,k)*ce)/Re
				adv=u(j,k)*(omega(j+1,k)-omega(j-1,k))/(gridpts_y(j+1)-gridpts_y(j-1)) + v(j,k)*(omega(j,k+1)-omega(j,k-1))/(gridpts_z(k+1)-gridpts_z(k-1))
				omega(j,k)=omega(j,k)+(diff-adv)*dt
				end if
			end do
		end do
		!omega=omega_out
		!deallocate(omega_out)
	end subroutine time_marching
	
	
	subroutine poisson_solver(omega,ny_l,ny_r,nz_l,nz_u,dt,gridpts_y,gridpts_z,psi,iter,res_tol,U0,Ly)		!based on pseudo-transient approach
		implicit none
		integer, intent(in)			:: ny_l,ny_r,nz_l,nz_u,iter
		real*8, intent(in)			:: dt,res_tol,U0,Ly
		real*8, allocatable			:: psi_out(:,:)!,diff_arr(:,:)
		real*8, allocatable, intent(inout)	:: psi(:,:)
		real*8, allocatable, intent(in)		:: omega(:,:),gridpts_y(:),gridpts_z(:)
		integer					:: j,k,i,iter_count
		real*8					:: diff,ce,we,ea,no,so
		allocate(psi_out(-ny_l:ny_r,-nz_l:nz_u))!,diff_arr(0:ny,0:nz))
		psi_out=0.0d0
		!diff_arr=0.0d0
		iter_count=0
		do i=1,iter
			psi_out=psi
			iter_count=iter_count+1
			do k=-nz_l+1,nz_u-1
				do j=-ny_l+1,ny_r-1
					if(k<=0 .and. j<=0) then
				
					else
					ce=2.0d0/(gridpts_y(j+1)-gridpts_y(j))/(gridpts_y(j)-gridpts_y(j-1)) + 2.0d0/(gridpts_z(k+1)-gridpts_z(k))/(gridpts_z(k)-gridpts_z(k-1))
					we=2.0d0/(gridpts_y(j+1)-gridpts_y(j-1))/(gridpts_y(j)-gridpts_y(j-1))
					ea=2.0d0/(gridpts_y(j+1)-gridpts_y(j))/(gridpts_y(j+1)-gridpts_y(j-1))
					no=2.0d0/(gridpts_z(k+1)-gridpts_z(k))/(gridpts_z(k+1)-gridpts_z(k-1))
					so=2.0d0/(gridpts_z(k+1)-gridpts_z(k-1))/(gridpts_z(k)-gridpts_z(k-1))
					diff=psi(j+1,k)*ea + psi(j-1,k)*we + psi(j,k+1)*no + psi(j,k-1)*so - psi(j,k)*ce + omega(j,k)
					psi(j,k)=psi(j,k)+diff*dt/(U0*Ly)
					end if
				end do
			end do
			if(norm2(psi_out-psi)<norm2(psi)*res_tol) exit
		end do
		print*, iter_count, "iter in poisson"
	end subroutine poisson_solver
	

end module solvers


module linear_systems
	contains
	subroutine u_and_v(u,v,psi,ny_l,ny_r,nz_l,nz_u,gridpts_y,gridpts_z)
		implicit none
		integer, intent(in)			:: ny_l,ny_r,nz_l,nz_u
		real*8, allocatable, intent(inout)	:: u(:,:),v(:,:)
		real*8, allocatable, intent(in)		:: psi(:,:),gridpts_y(:),gridpts_z(:)
		integer					:: j,k
		
		do k=-nz_l+1,nz_u-1
			do j=-ny_l+1,ny_r-1
				if(k<=0 .and. j<=0) then
				
				else
				u(j,k)=(psi(j,k+1)-psi(j,k-1))/(gridpts_z(k+1)-gridpts_z(k-1))
				v(j,k)=-(psi(j+1,k)-psi(j-1,k))/(gridpts_y(j+1)-gridpts_y(j-1))
				end if
			end do
		end do
	end subroutine u_and_v
	
	
	subroutine lin_sys_and_solve(omega,psi,u,v,U0,Re,ny_l,ny_r,nz_l,nz_u,dt,dt2,gridpts_y,gridpts_z,iter1,iter2,res_tol,D)
		use input_ibcs_and_fns
		use solvers
		implicit none
		integer, intent(in)			:: ny_l,ny_r,nz_l,nz_u,iter1,iter2
		real*8, intent(in)			:: dt,dt2,res_tol,U0,Re,D
		real*8, allocatable, intent(out)	:: psi(:,:),omega(:,:),u(:,:),v(:,:)
		real*8, allocatable, intent(in)		:: gridpts_y(:),gridpts_z(:)
		real*8, allocatable			:: temp_omega(:,:)
		integer					:: j,k,i,iter_count=0
		real*8					:: U_nd,dt_nd,dt2_nd
		
		allocate(psi(-ny_l:ny_r,-nz_l:nz_u),omega(-ny_l:ny_r,-nz_l:nz_u),u(-ny_l:ny_r,-nz_l:nz_u),v(-ny_l:ny_r,-nz_l:nz_u), temp_omega(-ny_l:ny_r,-nz_l:nz_u))
		U_nd=U0/U0
		dt_nd=dt/(D/U0)
		dt2_nd=dt2/(D/U0)
		psi=0.0d0;omega=0.0d0;u=0.0d0;v=0.0d0;temp_omega=0.0d0
		call apply_ics(psi,omega,ny_l,ny_r,nz_l,nz_u,gridpts_y,gridpts_z,U_nd,u,v)
		do i=1,iter1
			iter_count=iter_count+1
			temp_omega=omega
			call time_marching(omega,u,v,Re,ny_l,ny_r,nz_l,nz_u,dt_nd,gridpts_y,gridpts_z)
			call poisson_solver(omega,ny_l,ny_r,nz_l,nz_u,dt2_nd,gridpts_y,gridpts_z,psi,iter2,res_tol,U0,D)
			call apply_bcs(psi,omega,ny_l,ny_r,nz_l,nz_u,gridpts_y,gridpts_z,U_nd,u,v)
			call u_and_v(u,v,psi,ny_l,ny_r,nz_l,nz_u,gridpts_y,gridpts_z)
			if(norm2(omega-temp_omega)<res_tol) exit
		end do
		print*, iter_count, "iter in time marching"
		deallocate(temp_omega)
	end subroutine lin_sys_and_solve
end module linear_systems


program vorticity_stream_function
	use input_ibcs_and_fns
	use solvers
	use linear_systems
	implicit none
	real*8 			:: Ly_l,Ly_r,Lz_u,Lz_l,U0,Re,dt,dt2,res_tol,D
	integer			:: ny_l,ny_r,nz_l,nz_u,iter1,iter2,n,j,k
	real*8, allocatable	:: omega(:,:),psi(:,:),u(:,:),v(:,:),gridpts_y(:),gridpts_z(:)
	character(len=10)	:: Re_str
	
	call read_input(ny_l,ny_r,nz_l,nz_u,Ly_l,Ly_r,Lz_u,Lz_l,U0,Re,iter1,iter2,dt,dt2,res_tol,D)
	call generate_grid(ny_l,ny_r,nz_u,nz_l,Ly_l,Ly_r,Lz_u,Lz_l,gridpts_y,gridpts_z)
	call lin_sys_and_solve(omega,psi,u,v,U0,Re,ny_l,ny_r,nz_l,nz_u,dt,dt2,gridpts_y,gridpts_z,iter1,iter2,res_tol,D)
	
	write(Re_str,'(F6.2)') Re
	open(unit=2,file='backward_step_Re_'//trim(adjustl(Re_str))//'.dat')
	open(unit=4,file='separation_Re_'//trim(adjustl(Re_str))//'.dat')
	
		
	write(2,*),"y","	z","	psi","	omega","	u","	v"
	write(4,*),"index","	y","	z","	psi","	omega","	u","	v"
	do k=-nz_l,nz_u
	do j=-ny_l,ny_r
		write(2,'(1X, F8.4, 2X, F8.4, *(2X, F15.6))'),gridpts_y(j),gridpts_z(k),psi(j,k), &
		omega(j,k),u(j,k),v(j,k)
	end do
		write(2,*) " "
	end do

	do k=-nz_l,nz_u
		write(4,'(1X,I6,2X, F8.4, 2X, F8.4, *(2X, F15.6))'),k,gridpts_y(0),gridpts_z(k),psi(0,k), &
		omega(0,k),u(0,k),v(0,k)
	end do
	
	do j=0,ny_r
		if(u(j-1,-nz_l+1)<0 .and. u(j,-nz_l+1)>=0) then
			print*,"reattachment length = ",gridpts_y(j),"	U = ", u(j,-nz_l+1)
		exit
		end if
	end do
	
	close(2)
	close(4)
end program vorticity_stream_function
