      SUBROUTINE intp3d_many ( x, y, z, f, kt, ft, nx, ny, nz, nvars, xt, yt, zt)
!
      implicit none
!                                                          
!---------------------------------------------------------------------
!
!     purpose: interpolation of a function of three variables in an
!              equidistant(!!!) table.
!
!     method:  8-point Lagrange linear interpolation formula          
!
!     x        input vector of first  variable
!     y        input vector of second variable
!     z        input vector of third  variable
!
!     f        output vector of interpolated function values
!
!     kt       vector length of input and output vectors
!
!     ft       3d array of tabulated function values
!     nx       x-dimension of table
!     ny       y-dimension of table
!     nz       z-dimension of table
!     xt       vector of x-coordinates of table
!     yt       vector of y-coordinates of table
!     zt       vector of z-coordinates of table
!
!---------------------------------------------------------------------

      integer kt,nx,ny,nz,iv,nvars
      real*8 :: ft(nx,ny,nz,nvars)

      real*8 x(kt),y(kt),z(kt),f(kt,nvars)
      real*8 xt(nx),yt(ny),zt(nz)
      real*8 d1,d2,d3
!
!
      integer,parameter :: ktx = 1
      real*8  fh(ktx,8,nvars), delx(ktx), dely(ktx), delz(ktx), &
           a1(ktx,nvars), a2(ktx,nvars), a3(ktx,nvars), a4(ktx,nvars), &
           a5(ktx,nvars), a6(ktx,nvars), a7(ktx,nvars), a8(ktx,nvars)

      real*8 dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer n,ix,iy,iz

      IF (kt .GT. ktx)  STOP '***KTX**'
!
!
!------  determine spacing parameters of (equidistant!!!) table
!
      dx    = (xt(nx) - xt(1)) / FLOAT(nx-1)
      dy    = (yt(ny) - yt(1)) / FLOAT(ny-1)
      dz    = (zt(nz) - zt(1)) / FLOAT(nz-1)
!
      dxi   = 1. / dx
      dyi   = 1. / dy
      dzi   = 1. / dz
!
      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi
!
      dxyzi = dxi * dyi * dzi
!
!
!------- loop over all points to be interpolated
!
      dO  n = 1, kt                                            
!
!------- determine location in (equidistant!!!) table 
!                                                                  
         ix = 2 + INT( (x(n) - xt(1) - 1.e-10) * dxi )
         iy = 2 + INT( (y(n) - yt(1) - 1.e-10) * dyi )
         iz = 2 + INT( (z(n) - zt(1) - 1.e-10) * dzi )
!~          ix = 1 + NINT( (x(n) - xt(1) - 1.e-10) * dxi )
!~          iy = 1 + NINT( (y(n) - yt(1) - 1.e-10) * dyi )
!~          iz = 1 + NINT( (z(n) - zt(1) - 1.e-10) * dzi )
!                                                     
         ix = MAX( 2, MIN( ix, nx ) )
         iy = MAX( 2, MIN( iy, ny ) )
         iz = MAX( 2, MIN( iz, nz ) )
!
!~          write(*,*) "3d",ix,iy,iz , xt(ix),yt(iy),zt(iz)
!
!------- set-up auxiliary arrays for Lagrange interpolation
!                                                                 
         delx(n) = xt(ix) - x(n)
         dely(n) = yt(iy) - y(n)
         delz(n) = zt(iz) - z(n)
!      
         do iv = 1, nvars
            fh(n,1,iv) = ft(ix  , iy  , iz, iv  )                             
            fh(n,2,iv) = ft(ix-1, iy  , iz, iv  )                             
            fh(n,3,iv) = ft(ix  , iy-1, iz, iv  )                             
            fh(n,4,iv) = ft(ix  , iy  , iz-1, iv)                             
            fh(n,5,iv) = ft(ix-1, iy-1, iz, iv  )                             
            fh(n,6,iv) = ft(ix-1, iy  , iz-1, iv)                             
            fh(n,7,iv) = ft(ix  , iy-1, iz-1, iv)                             
            fh(n,8,iv) = ft(ix-1, iy-1, iz-1, iv)                             
!              
!------ set up coefficients of the interpolation polynomial and 
!       evaluate function values 
            !                                                    
            a1(n,iv) = fh(n,1,iv)                             
            a2(n,iv) = dxi   * ( fh(n,2,iv) - fh(n,1,iv) )       
            a3(n,iv) = dyi   * ( fh(n,3,iv) - fh(n,1,iv) )       
            a4(n,iv) = dzi   * ( fh(n,4,iv) - fh(n,1,iv) )       
            a5(n,iv) = dxyi  * ( fh(n,5,iv) - fh(n,2,iv) - fh(n,3,iv) + fh(n,1,iv) )
            a6(n,iv) = dxzi  * ( fh(n,6,iv) - fh(n,2,iv) - fh(n,4,iv) + fh(n,1,iv) )
            a7(n,iv) = dyzi  * ( fh(n,7,iv) - fh(n,3,iv) - fh(n,4,iv) + fh(n,1,iv) )
            a8(n,iv) = dxyzi * ( fh(n,8,iv) - fh(n,1,iv) + fh(n,2,iv) + fh(n,3,iv) + &
                 fh(n,4,iv) - fh(n,5,iv) - fh(n,6,iv) - fh(n,7,iv) )
!
            f(n,iv)  = a1(n,iv) +  a2(n,iv) * delx(n)                         &
                 +  a3(n,iv) * dely(n)                         &
                 +  a4(n,iv) * delz(n)                         &
                 +  a5(n,iv) * delx(n) * dely(n)               &
                 +  a6(n,iv) * delx(n) * delz(n)               &
                 +  a7(n,iv) * dely(n) * delz(n)               &
                 +  a8(n,iv) * delx(n) * dely(n) * delz(n)     
!
         enddo

      enddo
!~       write(*,*) ix,iy,iz,xt(ix),yt(iy),zt(iz),delz
!~       write(*,*) zt
!~       write(*,*) ft(:,iy,iz,2) 
!~       write(*,*) f(1,2) 
!~       write(*,*) f
!~       stop                           
	  
    end SUBROUTINE intp3d_many
 
	SUBROUTINE intp4d_many ( x, y, z, w, f, ft, nx, ny, nz, nw, xt, yt, zt,wt)
!
      implicit none
!                                                          
!---------------------------------------------------------------------
!
!     purpose: interpolation of a function of three variables in an
!              equidistant(!!!) table.
!
!     method:  8-point Lagrange linear interpolation formula          
!
!     x        input vector of first  variable
!     y        input vector of second variable
!     z        input vector of third  variable
!     w        input vector of fourth  variable
!
!     f        interpolated function value
!
!     ft       4d array of tabulated function values
!     nx       x-dimension of table
!     ny       y-dimension of table
!     nz       z-dimension of table
!     nw       w-dimension of table
!     xt       vector of x-coordinates of table
!     yt       vector of y-coordinates of table
!     zt       vector of z-coordinates of table
!     wt       vector of w-coordinates of table
!
!---------------------------------------------------------------------

      integer,intent(in) :: nx,ny,nz,nw
      real*8,intent(in) :: ft(nx,ny,nz,nw)
	  
      real*8,intent(in) :: x,y,z,w
      real*8,intent(in) ::  xt(nx),yt(ny),zt(nz),wt(nw)
	  real*8,intent(out) :: f

      real*8 :: d1,d2,d3
!
!
      integer,parameter :: ktx = 1
      real*8  fh(16), delx, dely, delz,delw, &
           a1, a2, a3, a4, &
           a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16

      real*8 dx,dy,dz,dw,dxi,dyi,dzi,dwi,dxyi,dxzi&
			,dyzi,dxwi,dywi,dzwi,dxyzi,dxywi,dxwzi,dwyzi,dxyzwi
      integer n,ix,iy,iz,iw

      
!
!
!------  determine spacing parameters of (equidistant!!!) table
!
      dx    = (xt(nx) - xt(1)) / FLOAT(nx-1)
      dy    = (yt(ny) - yt(1)) / FLOAT(ny-1)
      dz    = (zt(nz) - zt(1)) / FLOAT(nz-1)
      dw    = (wt(nw) - wt(1)) / FLOAT(nw-1)
!
      dxi   = 1. / dx
      dyi   = 1. / dy
      dzi   = 1. / dz
      dwi   = 1. / dw
!
      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi
      dxwi  = dxi * dwi
      dywi  = dyi * dwi
      dzwi  = dzi * dwi
!
      dxyzi = dxi * dyi * dzi
      dxywi = dxi * dyi * dwi
      dxwzi = dxi * dwi * dzi
      dwyzi = dwi * dyi * dzi
!
	  dxyzwi = dxi *dyi * dzi * dwi
!
!------- loop over all points to be interpolated
!
!
!------- determine location in (equidistant!!!) table 
!                                                                  
         ix = 2 + INT( (x - xt(1) - 1.e-10) * dxi )
         iy = 2 + INT( (y - yt(1) - 1.e-10) * dyi )
         iz = 2 + INT( (z - zt(1) - 1.e-10) * dzi )
         iw = 2 + INT( (w - wt(1) - 1.e-10) * dwi )
!                                                     
         ix = MAX( 2, MIN( ix, nx ) )
         iy = MAX( 2, MIN( iy, ny ) )
         iz = MAX( 2, MIN( iz, nz ) )
         iw = MAX( 2, MIN( iw, nw ) )
!
!         write(*,*) iy-1,iy,iy+1

!~ 		write(*,*) "4d",ix,iy,iz,iw

!------- set-up auxiliary arrays for Lagrange interpolation
!~ !                                                                 
         delx = xt(ix) - x
         dely = yt(iy) - y
         delz = zt(iz) - z
         delw = wt(iw) - w
!      
		fh(1) = ft(ix  , iy  , iz, iw  )                             
		fh(2) = ft(ix-1, iy  , iz, iw  )                             
		fh(3) = ft(ix  , iy-1, iz, iw  )                             
		fh(4) = ft(ix  , iy  , iz-1, iw)                             
		fh(5) = ft(ix  , iy  , iz, iw-1)                             
		fh(6) = ft(ix-1, iy-1, iz, iw )                             
		fh(7) = ft(ix-1, iy  , iz-1, iw)                             
		fh(8) = ft(ix  , iy-1, iz-1, iw)                             
		fh(9) = ft(ix  , iy, iz-1, iw-1)                             
		fh(10) = ft(ix  , iy-1, iz, iw-1)                             
		fh(11) = ft(ix-1, iy, iz, iw-1)                             
		fh(12) = ft(ix-1, iy-1, iz-1, iw)                             
		fh(13) = ft(ix, iy-1, iz-1, iw-1)                             
		fh(14) = ft(ix-1, iy, iz-1, iw-1)                                                          
		fh(15) = ft(ix-1, iy-1, iz, iw-1)                             
		fh(16) = ft(ix-1, iy-1, iz-1, iw-1)                             
!              
!------ set up coefficients of the interpolation polynomial and 
!       evaluate function values 
            !                                                    
		a1 = fh(1)                             
		a2 = dxi   * ( fh(2) - fh(1) )       
		a3 = dyi   * ( fh(3) - fh(1) )       
		a4 = dzi   * ( fh(4) - fh(1) )       
		a5 = dwi   * ( fh(5) - fh(1) )
		
		
		
		a6 = dxyi  * ( fh(6) - fh(2) - fh(3) + fh(1) )
		a7 = dxzi  * ( fh(7) - fh(2) - fh(4) + fh(1) )
		a8 = dyzi  * ( fh(8) - fh(3) - fh(4) + fh(1) )
		a9 = dxwi  * ( fh(11) - fh(5) - fh(2) + fh(1) )
		a10 = dywi  * ( fh(10) - fh(5) - fh(3) + fh(1) )
		a11 = dzwi  * ( fh(9) - fh(5) - fh(4) + fh(1) )
		
		a12 = dxyzi * ( fh(12) - fh(1) + fh(2) + fh(3) + &
			 fh(4) - fh(6) - fh(7) - fh(8) )
			 
		a13 = dxywi * ( fh(15) - fh(1) + fh(2) + fh(3) + &
			 fh(5) - fh(10) - fh(11) - fh(6) ) 
			 
		a14 = dxwzi * ( fh(14) - fh(1) + fh(2) + fh(4) + &
			 fh(5) - fh(7) - fh(9) - fh(11) )
		
		a15 = dwyzi * ( fh(13) - fh(1) + fh(3) + fh(4) + &
			 fh(5) - fh(8) - fh(9) - fh(10) )
		
		a16 = dxyzwi * ( fh(16) &
						-  fh(12) -  fh(13) -  fh(14) -  fh(15) &
						+ fh(1) & 
						-  fh(2) -  fh(3) -  fh(4) -  fh(5) &
						+  fh(11) -  fh(6)  &
						+  fh(10) -  fh(7)  &
						+  fh(9)  -  fh(8) )
						
!
	   
		f  = a1 +  a2 * delx  	&
			 +  a3 * dely               &
			 +  a4 * delz               &
			 +  a5 * delw               &
			 +  a6 * delx * dely        &
			 +  a7 * delx * delz        &
			 +  a8 * dely * delz	    &
			 +  a9 * delx * delw 	    &
			 +  a10 * dely * delw        &
			 +  a11 * delz * delw        &
			 +  a12 * delx * dely * delz &
			 +  a13 * delx * dely * delw &
			 +  a14 * delx * delz * delw &
			 +  a15 * dely * delz * delw &
			 +  a16 * delx * dely * delz * delw		
!~          f = ft(ix,iy,iz,iw)
		if ( f .LT. 0.0d0 ) then 
			write(*,*) (f-ft(ix,iy,iz,iw))/ft(ix,iy,iz,iw),ix,iy,iz,iw
			write(*,*) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16
			write(*,*) fh(1),fh(2),fh(3),fh(4),fh(5),fh(6),fh(7),fh(8),fh(9) &
						,fh(10),fh(11),fh(12),fh(13),fh(14),fh(15),fh(16)
			write(*,*) dxi,dyi,dzi,dwi,dxyi,dxzi,dyzi,dxwi,dywi,dzwi,dxyzi &
						, dxywi,dxwzi,dwyzi,dxyzwi	
			write(*,*) delx,dely,delz,delw
			write(*,*)		
		endif                                       
      
    end SUBROUTINE intp4d_many
 
