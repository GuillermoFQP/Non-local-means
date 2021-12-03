! Non local means denoising for HEALPix maps

program nlmeans

use healpix_modules

implicit none

!===============================================================
real(DP), allocatable :: map_in(:,:), map_out(:,:), F(:,:)
integer, allocatable  :: idx(:,:), rank(:,:)
integer               :: nside, npix, nmaps, ord, lmax, bmax, n, i
real(DP)              :: amount, fwhm, sigma(3)
character(len=80)     :: fin, fout, amount_arg, fwhm_arg, header(43)
real(DP), parameter   :: eps = 1.0d-7
!===============================================================
type grid
	integer               :: nside, n                            ! Map parameters
	real(DP), allocatable :: DM(:,:), F(:,:), UF(:,:), MM(:,:,:) ! Downgraded map, features, upgraded features and moment matrix
end type
!===============================================================

! Input parameters
if (nArguments() /= 4) then
	write (*,*) "Usage: nlmeans input_file output_file amount FWHM"
	call fatal_error('Not enough arguments')
end if

call getArgument(1, fin); call getArgument(2, fout)        ! Input and output file names
call getArgument(3, amount_arg); read(amount_arg,*) amount ! Sigma parameter scaling
call getArgument(4, fwhm_arg); read(fwhm_arg,*) fwhm       ! FWHM size of the gaussian beam in arcminutes

npix = getsize_fits(fin, nmaps=nmaps, nside=nside, ordering=ord)
n    = nside2npix(nside) - 1
bmax = int(180*240*sqrt(-log(2.0)*log(eps))/(pi*fwhm) + 0.5)
lmax = min(3*nside-1, bmax, 1000)

! Reading input map
allocate(map_in(0:n,nmaps)); call input_map(fin, map_in, npix, nmaps)

write (*,*) "Map read successfully."

! Computing Minkowski functionals
allocate(F(0:n,3)); call minkowski_functionals(map_in(:,1), nside, ord, lmax, fwhm, F(:,1), F(:,2), F(:,3))

write (*,*) "Minkowski functionals generated successfully."

! Estimating feature variances using L2-moment
allocate(idx(n+1,3), rank(n+1,3))

do i = 1, 3; call indexx(n+1, F(:,i), idx(:,i)); call rankx(n+1, idx(:,i), rank(:,i)); end do
 
sigma = sum(F*(2*rank-n-1), 1)/(n-1.0)**2 * sqrt(pi)
write (*,*) "Variances:", sigma

! Non-local means denoising procedure
allocate(map_out, mold=map_in)
! Normal procedure
!call nlmeans_denoising(map_in, nside, nmaps, 1.0/(amount*sigma)**2, F, map_out)
! Multigrid procedure
call nlmeans_multigrid(nside, nmaps, ord, 1.0/(amount*sigma)**2, F, map_in, map_out)

write (*,*) "Non-local means procedure applied successfully."

! Generating output file
call write_minimal_header(header, 'map', nside=nside, order=ord)
call output_map(map_out, header, fout)
write (*,*) "Output file generated successfully."

contains

subroutine minkowski_functionals(map_in, nside, ord, lmax, fwhm, M1, M2, M3)
	integer, intent(in)       :: nside, ord, lmax
	real(DP), intent(in)      :: map_in(0:12*nside**2-1), fwhm
	real(DP), intent(out)     :: M1(0:12*nside**2-1), M2(0:12*nside**2-1), M3(0:12*nside**2-1)
	real(DP), allocatable     :: I(:), D1(:,:), D2(:,:), cot(:), C(:)
	complex(DPC), allocatable :: alm(:,:,:)
	real(DP)                  :: theta, phi
	integer                   :: n, p
	
	n = nside2npix(nside) - 1
	
	allocate(I(0:n), alm(1,0:lmax,0:lmax), D1(0:n,2), D2(0:n,3), cot(0:n), C(0:n))
	
	do p = 0, n
		if (ord==1) call pix2ang_ring(nside, p, theta, phi)
		if (ord==2) call pix2ang_nest(nside, p, theta, phi)
		cot(p) = cotan(theta)
	end do
	
	call map2alm(nside, lmax, lmax, map_in, alm)
	call alter_alm(nside, lmax, lmax, fwhm, alm)
	call alm2map_der(nside, lmax, lmax, alm, I, D1, D2)
	
	! Covariant gradiend squared
	C = D1(:,1)**2 + D1(:,2)**2
	
	! Second covariant derivatives
	! D2(:,1) does not need any change
	D2(:,2) = D2(:,2) - cot * D1(:,2)
	D2(:,3) = D2(:,3) + cot * D1(:,1)	
	
	! Minkowski functional maps
	M1 = I
	M2 = sqrt(C)
	M3 = - (D2(:,1)*D1(:,2)**2 - 2.0*D2(:,2)*D1(:,1)*D1(:,2) + D2(:,3)*D1(:,1)**2)/C
	
	deallocate(I, alm, D1, D2, cot, C)
end subroutine

function nlmean1(nside, nmaps, sigma, map_in, V, F)
	integer, intent(in)  :: nside, nmaps
	real(DP), intent(in) :: sigma(3), map_in(0:12*nside**2-1,nmaps), V(3), F(0:12*nside**2-1,3)
	real(DP)             :: nlmean1(nmaps), WS(0:nmaps)
	integer              :: p, n
	
	n = nside2npix(nside) - 1; WS = 0.0
	
	do p = 0, n
		WS = WS + weight(sigma, V-F(p,:)) * [1.0, map_in(p,:)]
	end do
	
	nlmean1 = WS(1:nmaps) / WS(0)
end function

function nlmean2(nside, nmaps, w, map_in, V, F)
	integer, intent(in)  :: nside, nmaps
	real(DP), intent(in) :: w(3), map_in(0:12*nside**2-1,nmaps), V(3), F(0:12*nside**2-1,3)
	real(DP)             :: nlmean2(0:nmaps,10), u, WS(0:nmaps,10), q(3), d(0:nmaps), g(10)
	integer              :: p, n, k, j
	
	n = nside2npix(nside) - 1; WS = 0.0
	
	do p = 0, n
		q = w * (V - F(p,:))
		u = weight(w, V-F(p,:))
		g = [1.0, -q, q(1)*q(1), q(1)*q(2), q(1)*q(3), q(2)*q(2), q(2)*q(3), q(3)*q(3)]
		d = [1.0, map_in(p,:)]
		
		forall (j=0:nmaps, k=1:10) WS(j,k) = WS(j,k) + u*g(k)*d(j)
	end do
	
	nlmean2 = WS
end function	

subroutine nlmeans_denoising(map_in, nside, nmaps, w, F, map_out)
	integer, intent(in)   :: nside, nmaps
	real(DP), intent(in)  :: w(3), map_in(0:12*nside**2-1,nmaps), F(0:12*nside**2-1,3)
	real(DP), intent(out) :: map_out(0:12*nside**2-1,nmaps)
	integer               :: p, n
	
	n = nside2npix(nside) - 1
	
	!$OMP PARALLEL DO
	do p = 0, n; map_out(p,:) = nlmean1(nside, nmaps, w, map_in, F(p,:), F); end do
	!$OMP END PARALLEL DO
end subroutine

subroutine nlmeans_multigrid(fnside, nmaps, ord, w, F, map_in, map_out)
	integer, intent(in)     :: fnside, nmaps, ord
	real(DP), intent(in)    :: w(3), map_in(0:12*fnside**2-1,nmaps), F(0:12*fnside**2-1,3)
	real(DP), intent(out)   :: map_out(0:12*fnside**2-1,nmaps)
	integer                 :: fn, m, p, l, lvls
	type(grid), allocatable :: mg(:)
	
	fn = nside2npix(fnside) - 1
	
	! Total multigrid levels
	lvls = log(fnside/4.0) / log(2.0); if (lvls < 1) lvls = 1
	
	allocate(mg(lvls))
	
	! Initalizing grid structure
	do l = 1, lvls
		mg(l)%nside = ishft(fnside,1-l); mg(l)%n = nside2npix(mg(l)%nside) - 1
		
		! Downgraded map, upgraded features, features and accumulation matrix
		allocate(mg(l)%DM(0:mg(l)%n,nmaps), mg(l)%UF(0:mg(l)%n,3), mg(l)%F(0:mg(l)%n,3), mg(l)%MM(0:mg(l)%n,0:nmaps,10))
		
		if (l==1) then
			mg(l)%DM = map_in
			mg(l)%F   = F
			if (ord==1) then
				call convert_ring2nest(mg(l)%nside, mg(l)%DM)
				call convert_ring2nest(mg(l)%nside, mg(l)%F)
			end if
		end if
		if (l>1) then
			call udgrade_nest(mg(1)%DM, mg(1)%nside, mg(l)%DM, mg(l)%nside)
			call udgrade_nest(mg(1)%F, mg(1)%nside, mg(l)%F, mg(l)%nside)
		end if
	end do
	
	! Evaluating accumulation matrix and derivatives at the coarsest level
	!$OMP PARALLEL DO
	do p = 0, mg(lvls)%n
		mg(lvls)%MM(p,:,:) = nlmean2(mg(lvls)%nside, nmaps, w, mg(lvls)%DM, mg(lvls)%F(p,:), mg(lvls)%F)
	end do
	!$OMP END PARALLEL DO
	
	! Upgrading filtered map and correcting for differences in features
	do l = lvls - 1, 1, -1
		call udgrade_nest(mg(l+1)%F, mg(l+1)%nside, mg(l)%UF, mg(l)%nside)
		
		do i = 1, 10
			call udgrade_nest(mg(l+1)%MM(:,:,i), mg(l+1)%nside, mg(l)%MM(:,:,i), mg(l)%nside)
		end do
		
		!$OMP PARALLEL DO
		do p = 0, mg(l)%n
			! if (some criteria)
			call taylor_expansion(mg(l)%MM(p,:,:), mg(l)%F(p,:), mg(l)%UF(p,:), nmaps, mg(l)%DM(p,:))
			! else
			! mg(l)%MM(p,:,:) = nlmean2(mg(l)%nside, nmaps, w, mg(l)%DM, mg(l)%F(p,:), mg(l)%F)
			! end if
		end do
		!$OMP END PARALLEL DO
	end do
	
	! Output map
	!$OMP PARALLEL DO
	do p = 0, fn
		map_out(p,:) = mg(1)%MM(p,1:nmaps,1) / mg(1)%MM(p,0,1)
	end do
	!$OMP END PARALLEL DO
	
	if (ord==1) call convert_nest2ring(fnside, map_out)
	
	! Cleaning up
	do l = 1, lvls
		deallocate(mg(l)%DM, mg(l)%F, mg(l)%UF, mg(l)%MM)
	end do
	
	deallocate(mg)
end subroutine

subroutine taylor_expansion(moments, features, upgraded_features, nmaps, map_out)
	integer, intent(in)   :: nmaps
	real(DP), intent(in)  :: moments(0:nmaps,10), features(3), upgraded_features(3)
	real(DP), intent(out) :: map_out(nmaps)
	real(DP)              :: W
	real(DP), allocatable :: Wa(:), Wab(:), Q(:), Qa(:,:), Qab(:,:), S(:), grad(:,:), hess(:,:), df(:)
	integer               :: m
	
	allocate(Wa(3), Wab(6), Q(nmaps), Qa(nmaps,3), Qab(nmaps,6), S(nmaps), grad(nmaps,3), hess(nmaps,6), df(3))
	
	W        = moments(0,1)
	Wa(:)    = moments(0,2:4)
	Wab(:)   = moments(0,5:10)
	Q(:)     = moments(1:nmaps,1)
	Qa(:,:)  = moments(1:nmaps,2:4)
	Qab(:,:) = moments(1:nmaps,5:10)
	S(:)     = Q(:) / W
	df(:)    = features(:) - upgraded_features(:)
	
	forall (m=1:nmaps)
		! Computing the gradient and Hessian
		! Gradient components G1, G2 and G3
		grad(m,:) = (Qa(m,:) - S(m) * Wa(:)) / W
		! Hessian component H11
		hess(m,1) = (Qab(m,1) - S(m) * Wab(1) - grad(m,1) * Wa(1) - grad(m,1) * Wa(1)) / W
		! Hessian component H12
		hess(m,2) = (Qab(m,2) - S(m) * Wab(2) - grad(m,1) * Wa(2) - grad(m,2) * Wa(1)) / W
		! Hessian component H13
		hess(m,3) = (Qab(m,3) - S(m) * Wab(3) - grad(m,1) * Wa(3) - grad(m,3) * Wa(1)) / W
		! Hessian component H22
		hess(m,4) = (Qab(m,4) - S(m) * Wab(4) - grad(m,2) * Wa(2) - grad(m,2) * Wa(2)) / W
		! Hessian component H23
		hess(m,5) = (Qab(m,5) - S(m) * Wab(5) - grad(m,2) * Wa(3) - grad(m,3) * Wa(2)) / W
		! Hessian component H33
		hess(m,6) = (Qab(m,6) - S(m) * Wab(6) - grad(m,3) * Wa(3) - grad(m,3) * Wa(3)) / W
		! Correcting the gradient components and map values (Hessian components do not change)
		! Gradient component G1
		grad(m,1) = grad(m,1) + hess(m,1) * df(1) + hess(m,2) * df(2) + hess(m,3) * df(3)
		! Gradient component G2
		grad(m,2) = grad(m,2) + hess(m,2) * df(1) + hess(m,4) * df(2) + hess(m,5) * df(3)
		! Gradient component G3
		grad(m,3) = grad(m,3) + hess(m,3) * df(1) + hess(m,5) * df(2) + hess(m,6) * df(3)
		! Map values
		S(m) = S(m) + sum(grad(m,:)*df(:)) - 0.5*(hess(m,1)*df(1)*df(1) + hess(m,4)*df(2)*df(2) + hess(m,6)*df(3)*df(3)) &
		                                   -     (hess(m,2)*df(1)*df(2) + hess(m,3)*df(1)*df(3) + hess(m,5)*df(2)*df(3))
	end forall
	
	map_out = S
	
	deallocate(Wa, Wab, Q, Qa, Qab, S, grad, hess, df)
end subroutine

! Weighting function
function weight(w, x)
	real(DP), intent(in) :: w(3), x(3)
	real(DP)             :: weight
	
	weight = exp(-sum(w*x*x)/2.0)
end function

end program

