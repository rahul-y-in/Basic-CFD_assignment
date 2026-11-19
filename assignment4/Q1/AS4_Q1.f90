module input_ibcs_and_fns
	contains
 	 subroutine generate_grid(ny,nz,Ly,Lz,A,gridpts_y,gridpts_z)
		 implicit none
		 INTEGER,intent(in) 				:: ny,nz
		 INTEGER 					:: i,j,k
		 REAL*8,intent(in) 				:: Ly,Lz,A
		 REAL(kind=8) 					:: dy,dz
		 REAL(kind=8),allocatable, intent(out) 		:: gridpts_y(:),gridpts_z(:)
		 allocate(gridpts_y(0:ny),gridpts_z(0:nz))
		 dy=Ly/(ny)
		 dz=Lz/(nz)
		 
		 do i=0,ny
		   gridpts_y(i)= -(Ly/2.0d0) + dy*(i)
		 end do
		 do i=0,nz 
		   gridpts_z(i)= -(Lz/2.0d0) + dz*(i)  
		 end do
		 !This loop stretches the grid point
		  gridpts_y= tanh(A*gridpts_y)/tanh(A)/Ly
		  gridpts_z= tanh(A*gridpts_z)/tanh(A)/Ly 
		 print*, "Grid_generation done"
	end subroutine generate_grid


	integer*4 function corner(indx,ny,nz)
	      implicit none
	      integer*4,intent(in)::indx,ny,nz
	      if((indx==1) .OR. (indx==nz+1) .OR. (indx==nz+ny+1) .OR. (indx==2*nz+ny+1)) then
		corner = 1
	      else
		corner = 0
	      endif
	    return
	end function corner


	integer*4 function jb(indx,ny,nz)
	      implicit none
	      integer*4,intent(in)::indx,ny,nz
	      if(indx<nz+1) then
		 jb=0 
	      elseif(indx<(nz+ny+1)) then
		 jb=indx-nz-1          
	      elseif(indx<(2*nz+ny+1)) then
		 jb=ny          
	      else
		 jb=2*ny+2*nz-indx+1     
	      endif
	    return
    	end function jb


	integer*4 function kb(indx,ny,nz)
	      implicit none
	      integer*4,intent(in)::indx,ny,nz
	      if(indx<nz+1) then
		 kb=indx-1
	      elseif(indx<(nz+ny+1)) then 
		 kb=nz          
	      elseif(indx<(2*nz+ny+1)) then
		 kb=2*nz+ny-indx+1          
	      else
		 kb=0         
	      endif
	    return
	end function kb


      	integer*4 function ibound(j,k,ny,nz)
	      implicit none
	      integer*4 :: j,k,ny,nz
		 if(j==0) then 
		   ibound = k+1
		 elseif(k==nz) then
		   ibound = nz+1+j
		 elseif(j==ny) then
		   ibound = nz+ny+1 + (nz-k)
		 else
		   ibound = 2*nz+ny+1+ (ny-j)
		 endif
	      return
	end function ibound


	integer*4 function indx_j(indx,ny,nz)
	      implicit none
	      integer*4,intent(in)::indx,ny,nz
	      indx_j=mod(indx-1,ny+1)
	    return
    	end function indx_j
    	
    	
    	integer*4 function indx_k(indx,ny,nz)
	      implicit none
	      integer*4,intent(in)::indx,ny,nz
	      indx_k=(indx-1)/(ny+1)
	    return
    	end function indx_k


	real*8 function h(j1,k1,j2,k2,gridpts_y,gridpts_z)
		implicit none
	  real*8  ,allocatable, intent(in)              ::gridpts_y(:),gridpts_z(:)
	  integer, intent(in)				::j1,k1,j2,k2
		h = sqrt((gridpts_y(j1)-gridpts_y(j2))**2 + (gridpts_z(k1)-gridpts_z(k2))**2)
	end function h


	Subroutine apply_bcs(psi,omega,ny,nz,gridpts_y,gridpts_z,U0,u,v)
		implicit none
		integer, intent(in) 			::	ny,nz
		real*8, intent(in)			::	U0
		integer					::	i,j,k,nb
		real*8, allocatable, intent(inout)	::	psi(:,:),omega(:,:),u(:,:),v(:,:)
		real*8, allocatable, intent(in)		::	gridpts_y(:),gridpts_z(:)
		
		nb=2*(ny+nz)
		do i=1,nb
			j=jb(i,ny,nz)
			k=kb(i,ny,nz)
			psi(j,k)=0.0d0
			u(j,k)=0.0d0; v(j,k)=0.0d0
			if(i<nz+1) then
				omega(j,k)=2.0d0*(psi(j,k)-psi(j+1,k))/(h(j+1,k,j,k,gridpts_y,gridpts_z)**2)
			elseif(i<=ny+nz+1) then
				omega(j,k)=2.0d0*(psi(j,k)-psi(j,k-1)-U0*h(j,k-1,j,k,gridpts_y,gridpts_z))/(h(j,k-1,j,k,gridpts_y,gridpts_z)**2)
				u(j,k)=U0
			elseif(i<ny+2*nz+1) then
				omega(j,k)=2.0d0*(psi(j,k)-psi(j-1,k))/(h(j-1,k,j,k,gridpts_y,gridpts_z)**2)
			else
				omega(j,k)=2.0d0*(psi(j,k)-psi(j,k+1))/(h(j,k+1,j,k,gridpts_y,gridpts_z)**2)
			end if
		end do
	end subroutine apply_bcs
	

	Subroutine apply_ics(psi,omega,ny,nz,gridpts_y,gridpts_z,U0,u,v)
		implicit none
		integer, intent(in) 			::	ny,nz
		real*8, intent(in)			::	U0
		integer					::	i,j,k
		real*8, allocatable, intent(inout)	::	psi(:,:),omega(:,:),u(:,:),v(:,:)
		real*8, allocatable, intent(in)		::	gridpts_y(:),gridpts_z(:)
		
		psi=0.0d0
		omega=0.0d0
		u=0.0d0;v=0.0d0
		call apply_bcs(psi,omega,ny,nz,gridpts_y,gridpts_z,U0,u,v)
	end subroutine apply_ics


	subroutine read_input(ny,nz,Ly,Lz,A,U0,Re,iter1,iter2,dt,dt2,res_tol)
	  implicit none
	  real*8 	:: Ly,Lz,A,U0,dt,dt2,res_tol,Re
	  integer	:: ny,nz,iter1,iter2
	  
	  namelist/input_parameters/ny,nz,Ly,Lz,A,U0,iter1,iter2,dt,dt2,res_tol,Re
	  open(unit=1, file='input.txt')
	  read(1, nml=input_parameters)
	  close(1)
	end subroutine read_input

end module input_ibcs_and_fns


module solvers
	contains
	subroutine time_marching(omega,u,v,Re,ny,nz,dt,gridpts_y,gridpts_z)	!essentially FTCS
		implicit none
		integer, intent(in)			:: ny,nz
		real*8, intent(in)			:: dt,Re
		real*8, allocatable			:: omega_out(:,:)
		real*8, allocatable, intent(inout)	:: omega(:,:)
		real*8, allocatable, intent(in)		:: u(:,:),v(:,:),gridpts_y(:),gridpts_z(:)
		integer					:: j,k
		real*8					:: diff,adv,ce,we,ea,no,so
		!allocate(omega_out(0:ny,0:nz))
		do k=1,nz-1
			do j=1,ny-1
				ce=2.0d0/(gridpts_y(j+1)-gridpts_y(j))/(gridpts_y(j)-gridpts_y(j-1)) + 2.0d0/(gridpts_z(k+1)-gridpts_z(k))/(gridpts_z(k)-gridpts_z(k-1))
				we=2.0d0/(gridpts_y(j+1)-gridpts_y(j-1))/(gridpts_y(j)-gridpts_y(j-1))
				ea=2.0d0/(gridpts_y(j+1)-gridpts_y(j))/(gridpts_y(j+1)-gridpts_y(j-1))
				no=2.0d0/(gridpts_z(k+1)-gridpts_z(k))/(gridpts_z(k+1)-gridpts_z(k-1))
				so=2.0d0/(gridpts_z(k+1)-gridpts_z(k-1))/(gridpts_z(k)-gridpts_z(k-1))
				diff=(omega(j+1,k)*ea + omega(j-1,k)*we + omega(j,k+1)*no + omega(j,k-1)*so - omega(j,k)*ce)/Re
				adv=u(j,k)*(omega(j+1,k)-omega(j-1,k))/(gridpts_y(j+1)-gridpts_y(j-1)) + v(j,k)*(omega(j,k+1)-omega(j,k-1))/(gridpts_z(k+1)-gridpts_z(k-1))
				omega(j,k)=omega(j,k)+(diff-adv)*dt
			end do
		end do
		!omega=omega_out
		!deallocate(omega_out)
	end subroutine time_marching
	
	
	subroutine poisson_solver(omega,ny,nz,dt,gridpts_y,gridpts_z,psi,iter,res_tol,U0,Ly)		!based on pseudo-transient approach
		implicit none
		integer, intent(in)			:: ny,nz,iter
		real*8, intent(in)			:: dt,res_tol,U0,Ly
		real*8, allocatable			:: psi_out(:,:)!,diff_arr(:,:)
		real*8, allocatable, intent(inout)	:: psi(:,:)
		real*8, allocatable, intent(in)		:: omega(:,:),gridpts_y(:),gridpts_z(:)
		integer					:: j,k,i,iter_count
		real*8					:: diff,ce,we,ea,no,so
		allocate(psi_out(0:ny,0:nz))!,diff_arr(0:ny,0:nz))
		!diff_arr=0.0d0
		iter_count=0
		do i=1,iter
			psi_out=psi
			iter_count=iter_count+1
			do k=1,nz-1
				do j=1,ny-1
					ce=2.0d0/(gridpts_y(j+1)-gridpts_y(j))/(gridpts_y(j)-gridpts_y(j-1)) + 2.0d0/(gridpts_z(k+1)-gridpts_z(k))/(gridpts_z(k)-gridpts_z(k-1))
					we=2.0d0/(gridpts_y(j+1)-gridpts_y(j-1))/(gridpts_y(j)-gridpts_y(j-1))
					ea=2.0d0/(gridpts_y(j+1)-gridpts_y(j))/(gridpts_y(j+1)-gridpts_y(j-1))
					no=2.0d0/(gridpts_z(k+1)-gridpts_z(k))/(gridpts_z(k+1)-gridpts_z(k-1))
					so=2.0d0/(gridpts_z(k+1)-gridpts_z(k-1))/(gridpts_z(k)-gridpts_z(k-1))
					diff=psi(j+1,k)*ea + psi(j-1,k)*we + psi(j,k+1)*no + psi(j,k-1)*so - psi(j,k)*ce + omega(j,k)
					psi(j,k)=psi(j,k)+diff*dt/(U0*Ly)
				end do
			end do
			if(norm2(psi_out-psi)<norm2(psi)*res_tol) exit
		end do
		print*, iter_count, "iter in poisson"
	end subroutine poisson_solver
	

end module solvers


module linear_systems
	contains
	subroutine u_and_v(u,v,psi,ny,nz,gridpts_y,gridpts_z)
		implicit none
		integer, intent(in)			:: ny,nz
		real*8, allocatable, intent(inout)	:: u(:,:),v(:,:)
		real*8, allocatable, intent(in)		:: psi(:,:),gridpts_y(:),gridpts_z(:)
		integer					:: j,k
		
		do k=1,nz-1
			do j=1,ny-1
				u(j,k)=(psi(j,k+1)-psi(j,k-1))/(gridpts_z(k+1)-gridpts_z(k-1))
				v(j,k)=-(psi(j+1,k)-psi(j-1,k))/(gridpts_y(j+1)-gridpts_y(j-1))
			end do
		end do
	end subroutine u_and_v
	
	
	subroutine lin_sys_and_solve(omega,psi,u,v,U0,Re,ny,nz,dt,dt2,gridpts_y,gridpts_z,iter1,iter2,res_tol,Ly)
		use input_ibcs_and_fns
		use solvers
		implicit none
		integer, intent(in)			:: ny,nz,iter1,iter2
		real*8, intent(in)			:: dt,dt2,res_tol,U0,Re,Ly
		real*8, allocatable, intent(out)	:: psi(:,:),omega(:,:),u(:,:),v(:,:)
		real*8, allocatable, intent(in)		:: gridpts_y(:),gridpts_z(:)
		real*8, allocatable			:: temp_omega(:,:)
		integer					:: j,k,i,iter_count=0
		real*8					:: U_nd,dt_nd,dt2_nd
		
		allocate(psi(0:ny,0:nz),omega(0:ny,0:nz),u(0:ny,0:nz),v(0:ny,0:nz), temp_omega(0:ny,0:nz))
		U_nd=U0/U0
		dt_nd=dt/(Ly/U0)
		dt2_nd=dt2/(Ly/U0)
		call apply_ics(psi,omega,ny,nz,gridpts_y,gridpts_z,U_nd,u,v)
		do i=1,iter1
			iter_count=iter_count+1
			temp_omega=omega
			call time_marching(omega,u,v,Re,ny,nz,dt_nd,gridpts_y,gridpts_z)
			call poisson_solver(omega,ny,nz,dt2_nd,gridpts_y,gridpts_z,psi,iter2,res_tol,U0,Ly)
			call apply_bcs(psi,omega,ny,nz,gridpts_y,gridpts_z,U_nd,u,v)
			call u_and_v(u,v,psi,ny,nz,gridpts_y,gridpts_z)
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
	real*8 			:: Ly,Lz,A,U0,Re,dt,dt2,res_tol
	integer			:: ny,nz,iter1,iter2,n,j,k
	real*8, allocatable	:: omega(:,:),psi(:,:),u(:,:),v(:,:),gridpts_y(:),gridpts_z(:)
	character(len=10)	:: Re_str
	
	call read_input(ny,nz,Ly,Lz,A,U0,Re,iter1,iter2,dt,dt2,res_tol)
	call generate_grid(ny,nz,Ly,Lz,A,gridpts_y,gridpts_z)
	call lin_sys_and_solve(omega,psi,u,v,U0,Re,ny,nz,dt,dt2,gridpts_y,gridpts_z,iter1,iter2,res_tol,Ly)
	
	n=(ny+1)*(nz+1)
	write(Re_str,'(F6.2)') Re
	open(unit=2,file='driven_cavity_Re_'//trim(adjustl(Re_str))//'.dat')
	open(unit=4,file='vertical_centreline_Re_'//trim(adjustl(Re_str))//'.dat')
	open(unit=3,file='horizontal_centreline_Re_'//trim(adjustl(Re_str))//'.dat')
		
	write(2,*),"y","	z","	psi","	omega","	u","	v"
	write(3,*),"index","	y","	z","	psi","	omega","	u","	v"
	write(4,*),"index","	y","	z","	psi","	omega","	u","	v"
	do j=1,n
		write(2,'(1X, F8.4, 2X, F8.4, *(2X, F15.6))'),gridpts_y(indx_j(j,ny,nz)),gridpts_z(indx_k(j,ny,nz)),psi(indx_j(j,ny,nz),indx_k(j,ny,nz)), &
		omega(indx_j(j,ny,nz),indx_k(j,ny,nz)),u(indx_j(j,ny,nz),indx_k(j,ny,nz)),v(indx_j(j,ny,nz),indx_k(j,ny,nz))
		 	if(mod(j, ny+1) == 0) then
				write(2,*), "	"
			end if
	end do
	do j=0,ny
		write(3,'(1X,I6,2X, F8.4, 2X, F8.4, *(2X, F15.6))'),j,gridpts_y(j),gridpts_z(nz/2),psi(j,nz/2), &
		omega(j,nz/2),u(j,nz/2),v(j,nz/2)
	end do
	do k=0,nz
		write(4,'(1X,I6,2X, F8.4, 2X, F8.4, *(2X, F15.6))'),k,gridpts_y(ny/2),gridpts_z(k),psi(ny/2,k), &
		omega(ny/2,k),u(ny/2,k),v(ny/2,k)
	end do
	
	close(2)
	close(3)
	close(4)
end program vorticity_stream_function
