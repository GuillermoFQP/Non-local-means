! Compile with Intel Fortran:
! ifort -O3 nlmeans.f90 -I/usr/local/src/Healpix_3.80/includeifort -cm -w -sox -qopt-report=0 -qopenmp -fPIC -o nlmeans -L/usr/local/src/Healpix_3.80/libifort -L/usr/lib -L/usr/local/src/Healpix_3.80/lib -lhealpix -lhpxgif -lsharp -lcfitsio -Wl,-R/usr/lib -Wl,-R/usr/local/src/Healpix_3.80/lib -Wl,-R/usr/local/src/Healpix_3.80/libifort -lcurl
!
! To activate Intel Fortran:
! source /opt/intel/oneapi/setvars.sh

program nlmeans

use healpix_modules

implicit none

character(len=80)     :: fin, fout, header(1:60)
integer               :: nmaps, nside, npix, ord, npixtot, nneigh, p, q, n, i, list(1:8)
integer, parameter    :: ch=1       ! CHANNEL TO BE FILTERED
real(DP), allocatable :: map(:,:), newmap(:,:), lmv(:), nf(:), wsum(:)
real(DP)              :: sumi
real(DP), parameter   :: sigma=1.d0 ! FILTERING PARAMETER

! Parse arguments
if (nArguments() /= 2) call abort("ERROR: Cannot parse arguments.")

! Input map parameters
call getArgument(1, fin)
call getArgument(2, fout)
npixtot = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
npix    = nside2npix(nside)
n       = npix-1

! Allocating the maps
allocate(map(0:n, 1:nmaps), newmap(0:n, 1:nmaps))

! Reading the input map
call input_map(fin, map, npix, nmaps, header=header)
newmap = map

! Converting order to nested
select case (ord)
	case (1); call convert_ring2nest(nside, map); ord = 2
	case (2); continue
	case default; call abort("ERROR: Unknown ordering of input map.")
end select

! Allocating local mean value, normalizing factor and weigthed sum arrays
allocate(lmv(0:n), nf(0:n), wsum(0:n))

! Computing local mean values
do p = 0, n
	call neighbours_nest(nside, p, list, nneigh)
	sumi = map(p,ch)
	do i = 1, nneigh
		sumi = sumi+map(list(i),ch)
	end do
	lmv(p) = sumi/nneigh
end do

! Computing the filtered values
do p = 0, n
	wsum(p) = 0.0
	nf(p)   = 0.0
	do q = 0, n
		nf(p)   = nf(p)   + weight(lmv(p), lmv(q), sigma)
		wsum(p) = wsum(p) + map(q, ch)*weight(lmv(p), lmv(q), sigma)
	end do
	newmap(p,1) = wsum(p)/nf(p)
end do

! Deallocating the previous map and converting the order of the new map to ring
deallocate(map)
call convert_nest2ring(nside, newmap); ord = 1

! Generating the output file
call output_map(newmap, header, fout) 

contains

real function weight(x, y, h)
	real(DP) x, y, h, arg
	arg = -abs(x-y)**2/h**2
	weight = exp(arg)
end function weight

end program nlmeans
