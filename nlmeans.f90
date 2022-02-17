! Non local means denoising for HEALPix maps

program nlmeans

use healpix_modules

implicit none

!===============================================================
real(DP), allocatable :: map_in(:,:), map_out(:,:), F(:,:)
integer               :: nside, npix, nmaps, ord, lmax, bmax, n, i
real(DP)              :: FP, FWHM, FWHM_rad, w(3), radius
character(len=80)     :: fin, fout, fres, FP_arg, FWHM_arg, header(43)
!===============================================================

! Input parameters
if (nArguments() /= 4 .and. nArguments() /= 5) then
	write (*,*) "Usage: nlmeans input_file output_file FWHM filtering_parameter [residual_file]"
	call fatal_error('Not enough arguments')
end if

call getArgument(1, fin); call getArgument(2, fout)        ! Input and output file names
call getArgument(3, FWHM_arg); read(FWHM_arg,*) FWHM       ! FWHM size of the gaussian beam in arcminutes
call getArgument(4, FP_arg); read(FP_arg,*) FP             ! Filtering parameter
if (nArguments() == 5) call getArgument(5, fres)

FWHM_rad = (FWHM/60)*(pi/180)

npix   = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
n      = nside2npix(nside) - 1
bmax   = int(4*sqrt(-log(2.0)*log(1.0d-7))/FWHM_rad + 0.5)
lmax   = min(3*nside-1, bmax, 1000)
radius = min(max(64.0/nside, FWHM_rad), 0.5)

! Allocating arrays
allocate(map_in(0:n,nmaps), map_out(0:n,nmaps), F(0:n,3), source=0.0)

! Reading input map
call input_map(fin, map_in, npix, nmaps)

write (*,*) "Map read successfully."

! Computing Minkowski functionals
call invariant_features(map_in(:,1), nside, ord, lmax, FWHM, F, w)

write (*,*) "Minkowski functionals generated successfully."

! Non-local means denoising procedure
call nlmeans_denoising(map_in, nside, nmaps, ord, radius, F, FP*w, map_out)

write (*,*) "Non-local means procedure applied successfully."

! Generating output file
call write_minimal_header(header, 'map', nside=nside, order=ord)
call output_map(map_out, header, fout)
write (*,*) "Output file generated successfully."

! Generating residual map file
if (nArguments() == 5) then
	call output_map(map_in - map_out, header, fres)
	write (*,*) "Residual map file generated successfully."
end if

deallocate(map_in, map_out, F)

contains

subroutine invariant_features(map_in, nside, ord, lmax, FWHM, F, w)
	integer, intent(in)       :: nside, ord, lmax
	real(DP), intent(in)      :: map_in(0:12*nside**2-1), FWHM
	real(DP), intent(out)     :: F(0:12*nside**2-1,3), w(3)
	real(DP), allocatable     :: I(:), D1(:,:), D2(:,:), cot(:), CG2(:)
	complex(DPC), allocatable :: alm(:,:,:), rho_array(:)
	real(DP)                  :: theta, phi, delta, rho, sigma
	integer                   :: n, p
	
	n = nside2npix(nside) - 1
	
	allocate(I(0:n), alm(1,0:lmax,0:lmax), D1(0:n,2), D2(0:n,3), cot(0:n), CG2(0:n), rho_array(0:n))
	
	do p = 0, n
		if (ord==1) call pix2ang_ring(nside, p, theta, phi)
		if (ord==2) call pix2ang_nest(nside, p, theta, phi)
		cot(p) = cotan(theta)
	end do
	
	call map2alm(nside, lmax, lmax, map_in, alm)
	call alter_alm(nside, lmax, lmax, FWHM, alm)
	call alm2map_der(nside, lmax, lmax, alm, I, D1, D2)
	
	! Covariant gradiend squared
	CG2 = D1(:,1)**2 + D1(:,2)**2
	
	! Second covariant derivatives
	! D2(:,1) does not need any change
	D2(:,2) = D2(:,2) - cot * D1(:,2)
	D2(:,3) = D2(:,3) + cot * D1(:,1)	
	
	! First minkowski functional
	F(:,1) = I
	
	! Second Minkowski functional
	F(:,2) = sqrt(CG2)
	
	! Skeleton invariant
	F(:,3) = (D1(:,1)*D1(:,2)*(D2(:,3) - D2(:,1)) + D2(:,2)*(D1(:,1)**2 - D1(:,2)**2)) / CG2
	
	! Quantities related to the eatures metric (for 1st Mink functional + 2nd Mink functional + skeleton invariant)
	delta     = ((FWHM / 60) * (pi / 180)) / sqrt(8.0 * log(2.0))
	rho_array = ((D2(:,1) - D2(:,3))*(D1(:,1)**2 - D1(:,2)**2) + 4.0*D2(:,2)*D1(:,1)*D2(:,2))**2 / CG2**3
	rho       = sum(rho_array) / (n+1)
	
	! Features metric
	sigma = norm2(map_in - F(:,1)) / (n + 1)
	w = [1.0, 2.0 * delta**2, 12.0 * delta**4 / (3.0 + (6.0 * rho - 1.0) * delta**2)] / sigma**2
	
	deallocate(I, alm, D1, D2, cot, CG2, rho_array)
end subroutine

subroutine nlmeans_denoising(map_in, nside, nmaps, ord, radius, F, w, map_out)
	integer, intent(in)   :: nside, nmaps, ord
	real(DP), intent(in)  :: map_in(0:12*nside**2-1,nmaps), F(0:12*nside**2-1,3), w(3), radius
	real(DP), intent(out) :: map_out(0:12*nside**2-1,nmaps)
	integer               :: p, n
	
	n     = nside2npix(nside) - 1
	
	!$OMP PARALLEL DO
	do p = 0, n; map_out(p,:) = nlmean(p, nside, nmaps, ord, radius, w, map_in, F(p,:), F); end do
	!$OMP END PARALLEL DO
end subroutine

! Non-local mean at a given pixel
function nlmean(pixel, nside, nmaps, ord, radius, w, map_in, V, F)
	integer, intent(in)  :: nside, nmaps, ord, pixel
	real(DP), intent(in) :: map_in(0:12*nside**2-1,nmaps), V(3), F(0:12*nside**2-1,3), w(3), radius
	real(DP)             :: nlmean(nmaps), WS(0:nmaps), vec(3)
	integer, allocatable :: list(:)
	integer              :: p, n, nlist, l
	
	n = nside2npix(nside) - 1
	l = nside2npix(nside) * sin(radius/2.0)**2
	
	allocate(list(0:3*l/2), source=0)
	
	if (ord == 1) then
		call pix2vec_ring(nside, pixel, vec)
		call query_disc(nside, vec, radius, list, nlist)
	end if
	
	if (ord == 2) then
		call pix2vec_nest(nside, pixel, vec)
		call query_disc(nside, vec, radius, list, nlist, nest=1)
	end if
	
	WS = 0.0
	
	do p = 0, nlist-1
		WS = WS + weight(w, V-F(list(p),:)) * [1.0, map_in(list(p),:)]
	end do
	
	nlmean = WS(1:nmaps) / WS(0)
	
	deallocate(list)
end function	

! Weighting function
function weight(w, x)
	real(DP), intent(in) :: w(3), x(3)
	real(DP)             :: weight
	
	weight = exp(-sum(w*x*x)/2.0)
end function

end program
