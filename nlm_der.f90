! Compile with Intel Fortran:
! ifort -O3 nlm_der.f90 -I/usr/local/src/Healpix_3.80/includeifort -cm -w -sox -qopt-report=0 -qopenmp -fPIC -o nlm_der -L/usr/local/src/Healpix_3.80/libifort -L/usr/lib -L/usr/local/src/Healpix_3.80/lib -lhealpix -lhpxgif -lsharp -lcfitsio -Wl,-R/usr/lib -Wl,-R/usr/local/src/Healpix_3.80/lib -Wl,-R/usr/local/src/Healpix_3.80/libifort -lcurl

program nlm_der

use healpix_modules

implicit none

real(DP), allocatable     :: oldmap(:,:), map(:,:), newmap(:,:), der1(:,:), der1mag(:,:), lmv(:,:)
complex(DPC), allocatable :: alm(:,:,:) 
integer                   :: lmax, nside, npix, pixn, nmaps, mapsn, ord
real(DP)                  :: sigma
character(len=80)         :: fin, fout, sigma_arg, mapsn_arg, header(1:43)

!==================================================
! Input parameters
!--------------------------------------------------
if (nArguments() /= 4) then
	write (*,*) "Usage: nlm_der input_file output_file sigma nmaps"
	call fatal_error('Not enough arguments')
end if

call getArgument(1, fin)   ! Input file name
call getArgument(2, fout)  ! Output file name
call getArgument(3, sigma_arg) ! Sigma parameter
call getArgument(4, mapsn_arg) ! Analysis type (1=T and 3=TQU)

read(sigma_arg,*) sigma
read(mapsn_arg,*) mapsn

pixn = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
npix = nside2npix(nside)
lmax = 3*nside

if (nmaps < mapsn) call abort("ERROR: Invalid number.")
if (mapsn /= 1 .and. mapsn /= 3) call abort("ERROR: Invalid number.")
if (npix /= pixn) call abort("ERROR: Invalid NSIDE.")

allocate(oldmap(0:npix-1,1:nmaps))
call input_map(fin, oldmap, npix, nmaps)

allocate(map(0:npix-1,1:nmaps))
map = oldmap

write (*,*) "Map read successfully."
call sleep(1)
!==================================================

!==================================================
! Computing gradient magnitudes
!--------------------------------------------------
! alm coefficients
allocate(alm(1:mapsn,0:lmax,0:lmax))
if (mapsn==1) call map2alm(nside, lmax, lmax, map(:,1), alm)
if (mapsn==3) call map2alm(nside, lmax, lmax, map, alm)
write (*,*) "alm coefficients generated successfully."
call sleep(1)

! Map and derivatives from alm coefficients
allocate(der1(0:npix-1,2*mapsn))
if (mapsn==1) call alm2map_der(nside, lmax, lmax, alm, map(:,1), der1)
if (mapsn==3) call alm2map_der(nside, lmax, lmax, alm, map, der1)
write (*,*) "1st derivatives generated successfully."
call sleep(1)

map = oldmap

! Computing the gradient magnitudes
allocate(der1mag(0:npix-1,1:mapsn))
call der1_magnitudes(mapsn, der1, der1mag)
!==================================================

!==================================================
! Non-local means procedure
!--------------------------------------------------
! Converting order to nested
select case (ord)
	case (1)
		call convert_ring2nest(nside, map)
		call convert_ring2nest(nside, der1mag)
		ord = 2
	case (2)
		continue
	case default
		call abort("ERROR: Unknown ordering of input map.")
end select

! Computing local mean values
allocate(lmv(0:npix-1,mapsn))
call local_mean_values(nside, mapsn, map, lmv)

! Computing the filtered values
allocate(newmap(0:npix-1,mapsn))
call filtered_values(nside, mapsn, sigma, map, der1mag, lmv, newmap)

call convert_nest2ring(nside, newmap)
call convert_nest2ring(nside, der1mag)
ord = 1

write (*,*) "Non-local means procedure applied successfully."
call sleep(1)
!==================================================

!==================================================
! Generating output file
!--------------------------------------------------
call write_minimal_header(header, 'map', nside=nside, order=ord)
call output_map(newmap, header, fout)
write (*,*) "Output file generated successfully."
call sleep(1)
!==================================================

contains

subroutine linspace(xi, xf, array)
	real(DP), intent(in)  :: xi, xf
	real(DP), intent(out) :: array(:)
	real(DP)              :: ran
	integer               :: n, i

	n   = size(array)
	ran = xf - xi
	
	if (n == 0) return

	if (n == 1) then
		array(1) = xi
		return
	end if

	do i = 1, n
		array(i) = xi + ran*(i-1)/(n-1)
	end do
end subroutine

subroutine der1_magnitudes(nmaps, der1, der1mag)
	integer, intent(in)   :: nmaps
	real(DP), intent(in)  :: der1(:,:)
	real(DP), intent(out) :: der1mag(:,:)
	integer               :: i
	
	do i = 1, nmaps
		der1mag(:,i) = sqrt(der1(:,2*i-1)**2 + (der1(:,2*i)**2))
	end do
end subroutine

subroutine local_mean_values(nside, nmaps, map, lmv)
	integer, intent(in)   :: nside, nmaps
	real(DP), intent(in)  :: map(:,:)
	real(DP), intent(out) :: lmv(:,:)
	real(DP)              :: sumi
	integer               :: list(8), npix, nneigh, p, i, j
	
	npix = nside2npix(nside)
	
	do j = 1, nmaps
		do p = 0, npix-1
			call neighbours_nest(nside, p, list, nneigh)
			sumi = map(p,j)
			do i = 1, nneigh
				sumi = sumi + map(list(i),j)
			end do
			lmv(p,j) = sumi/nneigh
		end do
	end do
end subroutine

subroutine filtered_values(nside, nmaps, sigma, map, der1mag, lmv, newmap)
	integer, intent(in)   :: nmaps, nside
	real(DP), intent(in)  :: sigma, map(:,:), der1mag(:,:), lmv(:,:)
	real(DP), intent(out) :: newmap(:,:)
	real(DP), allocatable :: nf(:), wsum(:)
	integer               :: p, q, j
	
	npix = nside2npix(nside)
	
	allocate(nf(0:npix-1), wsum(0:npix-1))
	
	do j = 1, nmaps
		do p = 0, npix-1
			wsum(p) = 0.0
			nf(p)   = 0.0
			do q = 0, npix-1
				nf(p)   = nf(p)   +          weight(lmv(p,j), lmv(q,j), sigma) * weight(der1mag(p,j), der1mag(q,j), sigma)
				wsum(p) = wsum(p) + map(q,j)*weight(lmv(p,j), lmv(q,j), sigma) * weight(der1mag(p,j), der1mag(q,j), sigma)
			end do
			newmap(p,j) = wsum(p)/nf(p)
		end do
	end do
end subroutine

real function weight(x, y, h)
	real(DP) :: x, y, h, arg
	arg = -abs(x-y)**2/h**2
	weight = exp(arg)
end function weight

end program
