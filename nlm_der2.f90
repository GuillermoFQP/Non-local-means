! Compile with Intel Fortran:
! ifort -O3 nlm_der2.f90 -I/usr/local/src/Healpix_3.80/includeifort -cm -w -sox -qopt-report=0 -qopenmp -fPIC -o nlm_der2 -L/usr/local/src/Healpix_3.80/libifort -L/usr/lib -L/usr/local/src/Healpix_3.80/lib -lhealpix -lhpxgif -lsharp -lcfitsio -Wl,-R/usr/lib -Wl,-R/usr/local/src/Healpix_3.80/lib -Wl,-R/usr/local/src/Healpix_3.80/libifort -lcurl

program nlm_der

use healpix_modules

implicit none

real(DP), allocatable     :: oldmap(:,:), map(:,:), newmap(:,:), d1(:,:), d1inv(:,:), d2(:,:), d2inv(:,:), lmv(:,:)
complex(DPC), allocatable :: alm(:,:,:) 
integer                   :: lmax, nside, npix, nmaps, mapsn, ord, n
real(DP)                  :: sigma
character(len=80)         :: fin, fout, sigma_arg, mapsn_arg, header(43)

!==================================================
! Input parameters
!--------------------------------------------------
if (nArguments() /= 4) then
	write (*,*) "Usage: nlm_der input_file output_file sigma nmaps"
	call fatal_error('Not enough arguments')
end if

call getArgument(1, fin)       ! Input file name
call getArgument(2, fout)      ! Output file name
call getArgument(3, sigma_arg) ! Sigma parameter
call getArgument(4, mapsn_arg) ! Analysis type (1=T and 3=TQU)

read(sigma_arg,*) sigma
read(mapsn_arg,*) mapsn

npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
lmax = 3*nside
n    = nside2npix(nside) - 1

if (nmaps < mapsn) call fatal_error("Invalid argument for nmaps.")
if (mapsn /= 1 .and. mapsn /= 3) call fatal_error("Invalid argument for nmaps.")

allocate(oldmap(0:n,nmaps))
call input_map(fin, oldmap, npix, nmaps)

allocate(map(0:n,mapsn))
map = oldmap

write (*,*) "Map read successfully."
!==================================================

!==================================================
! Computing gradient magnitudes
!--------------------------------------------------
! alm coefficients
allocate(alm(mapsn,0:lmax,0:lmax))
if (mapsn==1) call map2alm(nside, lmax, lmax, map(:,1), alm)
if (mapsn==3) call map2alm(nside, lmax, lmax, map, alm)
write (*,*) "alm coefficients generated successfully."

! Map and derivatives from alm coefficients
allocate(d1(0:n,2*mapsn), d2(0:n,3*mapsn))
if (mapsn==1) call alm2map_der(nside, lmax, lmax, alm, map(:,1), d1, d2)
if (mapsn==3) call alm2map_der(nside, lmax, lmax, alm, map, d1, d2)
write (*,*) "1st and 2nd derivatives generated successfully."

map = oldmap

! Computing the gradient magnitudes
allocate(d1inv(0:n,mapsn), d2inv(0:n,mapsn))
call der_invariants(nside, mapsn, d1, d2, d1inv, d2inv)
!==================================================

!==================================================
! Non-local means denoising
!--------------------------------------------------
! Converting order to nested
select case (ord)
	case (1)
		call convert_ring2nest(nside, map)
		call convert_ring2nest(nside, d1inv)
		call convert_ring2nest(nside, d2inv)
		ord = 2
	case (2)
		continue
	case default
		call fatal_error("Unknown ordering of input map.")
end select

! Computing local mean values
allocate(lmv(0:n,mapsn))
call local_mean_values(nside, mapsn, map, lmv)

! Non-local means denoising procedure
allocate(newmap(0:n,mapsn))
call nlm_denoising(nside, mapsn, sigma, map, d1inv, d2inv, lmv, newmap)

! Converting order to ring
call convert_nest2ring(nside, newmap)
call convert_nest2ring(nside, d1inv)
ord = 1

write (*,*) "Non-local means procedure applied successfully."
!==================================================

!==================================================
! Generating output file
!--------------------------------------------------
call write_minimal_header(header, 'map', nside=nside, order=ord)
call output_map(newmap, header, fout)
write (*,*) "Output file generated successfully."
!==================================================

contains

!subroutine linspace(xi, xf, array)
!	real(DP), intent(in)  :: xi, xf
!	real(DP), intent(out) :: array(:)
!	real(DP)              :: ran
!	integer               :: n, i
!	
!	n   = size(array)
!	ran = xf - xi
!	
!	if (n == 0) return
!	
!	if (n == 1) then
!		array(1) = xi
!		return
!	end if
!	
!	do i = 1, n
!		array(i) = xi + ran*(i-1)/(n-1)
!	end do
!end subroutine

subroutine der_invariants(nside, nmaps, d1, d2, d1inv, d2inv)
	integer, intent(in)   :: nside, nmaps
	real(DP), intent(in)  :: d1(:,:), d2(:,:)
	real(DP), intent(out) :: d1inv(:,:), d2inv(:,:)
	real(DP), allocatable :: cot(:), c2(:,:)
	real(DP)              :: theta, phi
	integer               :: n, m, p
	
	n = nside2npix(nside) - 1
	
	allocate(cot(0:n), c2(0:n,3*nmaps))
	
	do p = 0, n
		call pix2ang_ring(nside, p, theta, phi)
		cot(p) = cotan(theta)
	end do
	
	do m = 1, nmaps
		! Second covariant derivatives
		c2(:,3*m-2) = d2(:,3*m-2)
		c2(:,3*m-1) = d2(:,3*m-1) - cot(:)*d1(:,3*m-1)
		c2(:,3*m  ) = d2(:,3*m  ) + cot(:)*d1(:,3*m  )
		! Invariants constructed from 1st and 2nd derivatives
		do p = 0, n
			d1inv(p,m) = sqrt(d1(p,2*m-1)**2 + d1(p,2*m)**2)
			d2inv(p,m) = -(c2(p,3*m-2)*d1(p,2*m)**2 - 2.0*c2(p,3*m-1)*d1(p,2*m-1)*d1(p,2*m) + c2(p,3*m)*d1(p,2*m-1)**2)/d1inv(p,m)
		end do
		
	end do
	
	deallocate(cot, c2)
end subroutine

subroutine local_mean_values(nside, nmaps, map, lmv)
	integer, intent(in)   :: nside, nmaps
	real(DP), intent(in)  :: map(:,:)
	real(DP), intent(out) :: lmv(:,:)
	real(DP)              :: sumi
	integer               :: list(8), nneigh, p, i, m
	
	n = nside2npix(nside) - 1
	
	do m = 1, nmaps
		do p = 0, n
			call neighbours_nest(nside, p, list, nneigh)
			sumi = map(p,m)
			do i = 1, nneigh
				sumi = sumi + map(list(i),m)
			end do
			lmv(p,m) = sumi / nneigh
		end do
	end do
end subroutine

subroutine nlm_denoising(nside, nmaps, sig, map, d1inv, d2inv, lmv, newmap)
	integer, intent(in)   :: nmaps, nside
	real(DP), intent(in)  :: sig, map(:,:), d1inv(:,:), d2inv(:,:), lmv(:,:)
	real(DP), intent(out) :: newmap(:,:)
	real(DP), allocatable :: nf(:), ws(:)
	integer               :: p, q, m
	
	n = nside2npix(nside) - 1
	
	allocate(nf(0:n), ws(0:n))
	
	do m = 1, nmaps
		do p = 0, n
			ws(p) = 0.0
			nf(p) = 0.0
			do q = 0, n
				! Normalization factor
				nf(p) = nf(p) +            W(lmv(p,m),lmv(q,m),sig) * W(d1inv(p,m),d1inv(q,m),sig) * W(d2inv(p,m),d2inv(q,m),sig)
				! Weighted sum
				ws(p) = ws(p) + map(q,m) * W(lmv(p,m),lmv(q,m),sig) * W(d1inv(p,m),d1inv(q,m),sig) * W(d2inv(p,m),d2inv(q,m),sig)
			end do
			newmap(p,m) = ws(p) / nf(p)
		end do
	end do
	
	deallocate(nf, ws)
end subroutine

! Weight function
real function W(x, y, h)
	real(DP) :: x, y, h, arg
	arg = - abs(x-y)**2 / h**2
	W   = exp(arg)
end function W

end program
