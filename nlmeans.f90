! Non-local means denoising method for HEALPix maps 

program nlmeans

use healpix_modules

implicit none

!===============================================================
real(DP), allocatable :: map_in(:,:), map_out(:,:), F(:,:)
integer               :: nside, npix, nmaps, ord, n, lcut, lmax, i
real(DP)              :: alpha, fwhm, w(3), radius, fwhm_rad, bl_min
character(len=80)     :: fin, fout, fres, arg, header(43)
!===============================================================

! Input parameters
if (nArguments() /= 4 .and. nArguments() /= 5) then
	write (*,'(/,X,A)') "Usage: nlmeans input_file output_file FWHM filtering_parameter [residual_file]"
	call fatal_error('Not enough arguments')
end if

call getArgument(1, fin)                         ! Input file name
call getArgument(2, fout)                        ! Output file name
call getArgument(3, arg)
read(arg,*) fwhm                                 ! FWHM size of the gaussian beam in arcminutes
call getArgument(4, arg)
read(arg,*) alpha                                ! Filtering parameter
if (nArguments() == 5) call getArgument(5, fres) ! Residual file name

! Map parameters
npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
n = nside2npix(nside) - 1                                            ! Total number of pixels minus one
fwhm_rad = (fwhm/60.0) * (pi/180.0)                                  ! FWHM parameter in radians
bl_min = 1.0D-8                                                      ! Cutoff value for Gaussian beam
lcut = int(sqrt(0.25 - 16.0*log(bl_min)*log(2.0)/fwhm_rad**2) - 0.5) ! Cutoff value for "l" due to Gaussian beam
lmax = min(3*nside - 1, lcut, 4000)                                  ! Maximum "l" for interpolation
radius = min(max(64.0/nside, fwhm_rad), pi/6.0)                      ! Radius for filtering

! Allocating arrays
allocate(map_in(0:n,nmaps), map_out(0:n,nmaps), F(0:n,3), source=0.0)

! Reading input map
call input_map(fin, map_in, npix, nmaps)

write (*,'(/,X,A,/,X,A)') "Map read successfully.", "Computing Minkowski functionals."

! Computing Minkowski functionals
call invariant_features(map_in(:,1), nside, ord, lmax, fwhm, F, w)

! All the following subroutines are designed for maps in NESTED ordering
if (ord == 1) call convert_ring2nest(nside, map_in)

write (*,*) "Map denoising in progress."

! Non-local means denoising subroutine
call non_local_means(map_in(:,1), nside, radius, F, w/(alpha**2), map_out(:,1))

! Go back to RING ordering if necessary
if (ord == 1) then; call convert_nest2ring(nside, map_in); call convert_nest2ring(nside, map_out); end if

! Generating output file
call write_minimal_header(header, 'map', nside=nside, order=ord); call output_map(map_out, header, fout)
write (*,*) "Output file generated successfully."

! Generating residual map file
if (nArguments() == 5) then
	call output_map(map_in - map_out, header, fres)
	write (*,*) "Residual map file generated successfully."
end if

deallocate(map_in, map_out, F)

contains

subroutine invariant_features(map_in, nside, ord, lmax, fwhm, F, w)
	integer, intent(in)       :: nside, ord, lmax
	real(DP), intent(inout)   :: map_in(0:12*nside**2-1)
	real(DP), intent(in)      :: fwhm
	real(DP), intent(out)     :: F(0:12*nside**2-1,3), w(3)
	real(DP), allocatable     :: I(:), D1(:,:), D2(:,:), cot(:), CG2(:)
	complex(DPC), allocatable :: alm(:,:,:), rho_map(:)
	real(DP)                  :: theta, phi, delta, rho, sigma
	integer                   :: j, p
	
	n = nside2npix(nside) - 1
	
	allocate(I(0:n), alm(1,0:lmax,0:lmax), D1(0:n,2), D2(0:n,3), cot(0:n), CG2(0:n), rho_map(0:n))
	
	! The subroutines used here are designed for maps in RING ordering
	if (ord == 2) call convert_nest2ring(nside, map_in)
	
	! Cotangents of the polar angles
	do p = 0, n; call pix2ang_ring(nside, p, theta, phi); cot(p) = cotan(theta); end do
	
	! Computing Gaussian-smoothed map and derivatives
	call map2alm(nside, lmax, lmax, map_in, alm); call alter_alm(nside, lmax, lmax, fwhm, alm); call alm2map_der(nside, lmax, lmax, alm, I, D1, D2)
	
	! Covariant gradiend squared and second covariant derivatives
	CG2 = D1(:,1)**2 + D1(:,2)**2; D2(:,2) = D2(:,2) - cot * D1(:,2); D2(:,3) = D2(:,3) + cot * D1(:,1)
	
	! Feature space (Gaussian-smoothed map, 1st Minkowski functional and Skeleton maps)
	F(:,1) = I; F(:,2) = sqrt(CG2); F(:,3) = (D1(:,1)*D1(:,2)*(D2(:,3) - D2(:,1)) + D2(:,2)*(D1(:,1)**2 - D1(:,2)**2)) / CG2
	
	! Noise variance parameters
	sigma   = norm2(map_in-F(:,1)) / (n+1)
	delta   = (fwhm/60.0) * (pi/180.0) / sqrt(8.0*log(2.0))
	rho_map = ((D2(:,1)-D2(:,3))*(D1(:,1)**2-D1(:,2)**2)+4.0*D2(:,2)*D1(:,1)*D2(:,2))**2 / CG2**3
	rho     = sum(rho_map) / (n+1)
	
	! Features metric
	w = [1.0, 2.0 * delta**2, 12.0*delta**4 / (3.0+(6.0*rho-1.0)*delta**2)] / sigma**2
    
	! Go back to NESTED ordering if necessary
	if (ord == 2) then; call convert_ring2nest(nside, map_in); do j = 1, 3; call convert_ring2nest(nside, F(:,j)); end do; end if
	
	deallocate(I, alm, D1, D2, cot, CG2, rho_map)

end subroutine invariant_features

! Non-local means denoising subroutine for an input map in NESTED ordering
subroutine non_local_means(map_in, nside, radius, F, w, map_out)
	integer, intent(in)   :: nside
	real(DP), intent(in)  :: map_in(0:12*nside**2-1), F(0:12*nside**2-1,3), w(3), radius
	real(DP), intent(out) :: map_out(0:12*nside**2-1)
	real(DP)              :: vj(3), WS(0:1)
	integer, allocatable  :: list(:)
	integer               :: i, j, n, l, nlist
	
	n = nside2npix(nside) - 1; l = nside2npix(nside) * sin(radius/2.0)**2
	
	allocate(list(0:3*l/2))
	
	!$OMP PARALLEL PRIVATE(i, vj, WS, list, nlist, j) SHARED(n, nside, radius, w, F, map_out)
	!$OMP DO SCHEDULE(DYNAMIC)
	do i = 0, n
		! Disc around pixel "i"
		call pix2vec_nest(nside, i, vj); call query_disc(nside, vj, radius, list, nlist, nest=1)
		
		! Weighted sum WS(1) and normalization constant WS(0)
		WS = 0.0; do j = 0, nlist-1; WS = WS + weight(w, F(i,:)-F(list(j),:)) * [1.0, map_in(list(j))]; end do
		
		! Filtered value
		map_out(i) = WS(1) / WS(0)
	end do
	!$OMP END DO
	!$OMP END PARALLEL
	
	deallocate(list)
    
end subroutine non_local_means

! Weight function
pure function weight(w, x)
	real(DP), intent(in) :: w(3), x(3)
	real(DP)             :: weight
	
	weight = exp(-sum(w*x*x)/2.0)
	
end function weight

end program
