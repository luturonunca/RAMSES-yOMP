module imf_commons
  use amr_parameters, only: dp
  implicit none

  ! Lookup table: eta_sn(M_cl) and M_SNII(M_cl) via Chabrier + Weidner-Kroupa
  ! Second index: high-mass slope column, 1 = Chabrier (chab_alpha),
  ! 2 = top-heavy (sf_cluster_imf_s_th), assigned at birth from alpha_vir
  integer, parameter :: nimf = 200
  real(dp) :: imf_mcl  (nimf)     ! cluster mass grid [Msun], log-spaced
  real(dp) :: imf_eta  (nimf,2)   ! N_SN per Msun of stars formed
  real(dp) :: imf_msnii(nimf,2)   ! average SN progenitor mass [Msun]
  real(dp) :: imf_s    (2)        ! high-mass slope s (dN/dm ~ m^-s) of each column

  ! Maximum clusters drawn from CMF per SF event (safety ceiling)
  integer, parameter :: max_clusters_per_event = 10000

  ! Per-cell characteristic cluster mass [Msun], lognormal width, and
  ! fragmentation coherence radius [code units], for sf_cluster_kernel=
  ! 'lognormal'. Filled in starform2 (sf_model=1 branch), consumed by the
  ! CbC pre-pass in star_formation. Indexed like flag2.
  real(dp), allocatable :: cbc_mclchar_buf(:), cbc_sigmalnm_buf(:), cbc_rchar_buf(:)

  ! Per-cell continuous SF-law gas mass [code units], captured in starform2
  ! before the Poisson quantization (mgas=nstar*mstar) overwrites it. Used
  ! as the CbC mass budget M_sf instead of the quantized n*mstar. Filled
  ! whenever sf_cluster_sampling=.true. (both kernels). Indexed like flag2.
  real(dp), allocatable :: cbc_mgas_buf(:)

  ! Per-cell virial parameter alpha0 (2Ekin/|Egrav| at cell scale), filled
  ! in starform2 (sf_model=1 branch). Sets the IMF slope of each CbC event
  ! at birth. Indexed like flag2.
  real(dp), allocatable :: cbc_avir_buf(:)

  ! Chabrier (2003) single-object IMF constants
  real(dp), parameter :: mstar_low  = 0.08d0   ! IMF lower stellar mass limit [Msun]
  real(dp), parameter :: mstar_high = 150.0d0  ! IMF upper stellar mass limit [Msun]
  real(dp), parameter :: msn_thresh = 8.0d0    ! minimum SN progenitor mass [Msun]
  real(dp), parameter :: chab_mc    = 0.079d0  ! lognormal characteristic mass [Msun]
  real(dp), parameter :: chab_sig   = 0.69d0   ! lognormal sigma (in log10)
  real(dp), parameter :: chab_alpha = 2.3d0    ! power-law slope for m > 1 Msun

contains

  !------------------------------------------------------------------
  ! Chabrier (2003) IMF: dN/dm (unnormalized)
  ! Lognormal for m <= 1 Msun, power-law m^(-s) for m > 1 Msun,
  ! continuous at m = 1 Msun. Only s changes between IMF columns.
  !------------------------------------------------------------------
  pure function xi_chabrier(m, s) result(xi)
    real(dp), intent(in) :: m, s
    real(dp) :: xi
    real(dp), parameter :: cont = &
         exp(-(log10(chab_mc))**2 / (2.0d0*chab_sig**2))   ! continuity at m=1
    if(m <= 1.0d0) then
      xi = (1.0d0/m) * exp(-(log10(m) - log10(chab_mc))**2 &
           & / (2.0d0*chab_sig**2))
    else
      xi = cont * m**(-s)
    endif
  end function xi_chabrier

  !------------------------------------------------------------------
  ! Number integral: int_{m1}^{m2} xi(m) dm  (midpoint, log-spaced)
  !------------------------------------------------------------------
  function integ_n(m1, m2, s) result(res)
    real(dp), intent(in) :: m1, m2, s
    real(dp) :: res
    integer, parameter :: np = 4000
    integer :: k
    real(dp) :: lm1, lm2, dlm, mk
    res = 0.0d0
    if(m2 <= m1) return
    lm1 = log(m1);  lm2 = log(m2)
    dlm = (lm2 - lm1) / np
    do k = 1, np
      mk  = exp(lm1 + (k - 0.5d0)*dlm)
      res = res + xi_chabrier(mk, s) * mk * dlm   ! dm = m d(ln m)
    enddo
  end function integ_n

  !------------------------------------------------------------------
  ! Mass integral: int_{m1}^{m2} m * xi(m) dm  (midpoint, log-spaced)
  !------------------------------------------------------------------
  function integ_m(m1, m2, s) result(res)
    real(dp), intent(in) :: m1, m2, s
    real(dp) :: res
    integer, parameter :: np = 4000
    integer :: k
    real(dp) :: lm1, lm2, dlm, mk
    res = 0.0d0
    if(m2 <= m1) return
    lm1 = log(m1);  lm2 = log(m2)
    dlm = (lm2 - lm1) / np
    do k = 1, np
      mk  = exp(lm1 + (k - 0.5d0)*dlm)
      res = res + mk * xi_chabrier(mk, s) * mk * dlm
    enddo
  end function integ_m

  !------------------------------------------------------------------
  ! Weidner-Kroupa (2004): m_max(M_cl) via bisection in log-space.
  ! Solves: integ_n(mmax, mstar_high) = integ_m(mstar_low, mstar_high) / mcl
  ! i.e., one expected star more massive than mmax in a cluster of mass mcl.
  !------------------------------------------------------------------
  function wk_mmax(mcl, s) result(mmax)
    real(dp), intent(in) :: mcl, s
    real(dp) :: mmax
    real(dp) :: itot, target_n, nmid, mlo, mhi, mmid
    integer :: iter

    itot     = integ_m(mstar_low, mstar_high, s)
    target_n = itot / mcl

    ! Cluster large enough: mmax saturates at stellar cap
    if(target_n <= 0.0d0) then
      mmax = mstar_high; return
    endif
    ! Cluster too small: no star can be expected above mstar_low
    if(target_n >= integ_n(mstar_low, mstar_high, s)) then
      mmax = mstar_low; return
    endif

    ! Bisect in log-space
    mlo = mstar_low
    mhi = mstar_high
    do iter = 1, 80
      mmid = sqrt(mlo * mhi)
      nmid = integ_n(mmid, mstar_high, s)
      if(nmid > target_n) then
        mlo = mmid
      else
        mhi = mmid
      endif
      if((mhi - mlo) / mhi < 1.0d-7) exit
    enddo
    mmax = min(sqrt(mlo * mhi), mstar_high)
  end function wk_mmax

  !------------------------------------------------------------------
  ! Build the IMF lookup table.
  ! Called once at startup when sf_cluster_sampling=.true.
  !------------------------------------------------------------------
  subroutine init_imf_table(mcl_min, mcl_max, s_th)
    use amr_commons, only: myid
    real(dp), intent(in) :: mcl_min, mcl_max   ! cluster mass range [Msun]
    real(dp), intent(in) :: s_th               ! top-heavy high-mass slope
    integer :: i, itab
    real(dp) :: mcl, mmax, nsn, mtot, msnii_tot, dlm, eta_lo, eta_hi

    dlm = log(mcl_max / mcl_min) / dble(nimf - 1)
    imf_s(1) = chab_alpha
    imf_s(2) = s_th

    do itab = 1, 2
    do i = 1, nimf
      mcl  = mcl_min * exp(dble(i-1) * dlm)
      mmax = wk_mmax(mcl, imf_s(itab))

      imf_mcl(i) = mcl

      ! Same xi for the SN count and the total mass: each column is
      ! renormalised to the cluster birth mass
      if(mmax > msn_thresh) then
        nsn       = integ_n(msn_thresh, mmax, imf_s(itab))
        mtot      = integ_m(mstar_low,  mmax, imf_s(itab))
        msnii_tot = integ_m(msn_thresh, mmax, imf_s(itab))
        imf_eta  (i,itab) = nsn       / max(mtot,  1.0d-30)
        imf_msnii(i,itab) = msnii_tot / max(nsn,   1.0d-30)
      else
        imf_eta  (i,itab) = 0.0d0
        imf_msnii(i,itab) = 0.0d0
      endif
    enddo
    enddo

    if(myid == 1) then
      write(*,*) '>>> CbC: IMF table initialized (Chabrier + Weidner-Kroupa)'
      write(*,'(A,2ES12.4)') '    M_cl range [Msun]        : ', mcl_min, mcl_max
      do itab = 1, 2
        call eta_sn_cluster(1.0d4, itab, eta_lo, mmax)
        call eta_sn_cluster(1.0d6, itab, eta_hi, mmax)
        write(*,'(A,F6.3)')    '    high-mass slope s        : ', imf_s(itab)
        write(*,'(A,ES12.4)')  '    eta_sn at M_cl=1e4 Msun  : ', eta_lo
        write(*,'(A,ES12.4)')  '    eta_sn at M_cl=1e6 Msun  : ', eta_hi
      enddo
    endif
  end subroutine init_imf_table

  !------------------------------------------------------------------
  ! Interpolate eta_sn and M_SNII for a given cluster mass [Msun]
  ! in IMF column itab (1 = Chabrier, 2 = top-heavy).
  !------------------------------------------------------------------
  subroutine eta_sn_cluster(mcl, itab, eta, msnii)
    real(dp), intent(in)  :: mcl
    integer,  intent(in)  :: itab
    real(dp), intent(out) :: eta, msnii
    integer :: i
    real(dp) :: t

    if(mcl <= imf_mcl(1)) then
      eta = imf_eta(1,itab);  msnii = imf_msnii(1,itab);  return
    endif
    if(mcl >= imf_mcl(nimf)) then
      eta = imf_eta(nimf,itab);  msnii = imf_msnii(nimf,itab);  return
    endif
    i = 1
    do while(imf_mcl(i+1) < mcl .and. i < nimf-1)
      i = i + 1
    enddo
    t     = log(mcl / imf_mcl(i)) / log(imf_mcl(i+1) / imf_mcl(i))
    eta   = imf_eta  (i,itab) + t*(imf_eta  (i+1,itab) - imf_eta  (i,itab))
    msnii = imf_msnii(i,itab) + t*(imf_msnii(i+1,itab) - imf_msnii(i,itab))
  end subroutine eta_sn_cluster

  !------------------------------------------------------------------
  ! Decompose a SF event into clusters from CMF (dN/dM ~ M^{-2}).
  ! Environment-dependent caps privilege massive clusters at high Mj:
  !   mmin_eff = max(mmin0, mmin0*(Mj/Mj_ref)^delta)
  !   mmax_eff = min(mmax,  f_cap*Mj)
  ! Mass is exactly conserved: remainder added to last cluster.
  ! If M_sf < mmin_eff, returns one cluster with the full mass.
  ! All masses in Msun.
  !------------------------------------------------------------------
  subroutine sample_cmf_clusters(M_sf, Mj, mmin0, mmax, f_cap, Mj_ref, delta, &
                                  seed, cluster_masses, n_cl)
    use random
    real(dp), intent(in)    :: M_sf, Mj, mmin0, mmax, f_cap, Mj_ref, delta
    integer,  intent(inout) :: seed(IRandNumSize)
    real(dp), intent(out)   :: cluster_masses(max_clusters_per_event)
    integer,  intent(out)   :: n_cl

    real(dp)     :: M_rem, mmin_eff, mmax_eff, mmax_draw, inv_mmin, inv_mmax, mcl
    real(kind=8) :: RandNum

    ! Lower cap: raises floor at high Mj, privileges massive clusters
    mmin_eff = mmin0
    if(Mj > 0.0d0 .and. Mj_ref > 0.0d0) &
      mmin_eff = mmin0 * (Mj / Mj_ref)**delta
    mmin_eff = max(mmin_eff, mmin0)

    ! Upper cap
    mmax_eff = mmax
    if(f_cap * Mj < mmax_eff) mmax_eff = f_cap * Mj

    n_cl     = 0
    M_rem    = M_sf
    inv_mmin = 1.0d0 / mmin_eff

    do while(M_rem >= mmin_eff .and. n_cl < max_clusters_per_event)
      mmax_draw = min(mmax_eff, M_rem)
      if(mmax_draw <= mmin_eff) exit

      ! Inverse CDF for dN/dM ~ M^{-2}:
      ! m = 1 / (1/mmin - u*(1/mmin - 1/mmax))
      inv_mmax = 1.0d0 / mmax_draw
      call ranf(seed, RandNum)
      mcl = 1.0d0 / (inv_mmin - dble(RandNum) * (inv_mmin - inv_mmax))

      n_cl = n_cl + 1
      cluster_masses(n_cl) = mcl
      M_rem = M_rem - mcl
    end do

    ! Exact mass conservation: remainder into last cluster
    if(n_cl > 0) then
      cluster_masses(n_cl) = cluster_masses(n_cl) + M_rem
    else
      ! M_sf below mmin_eff: single cluster with full mass
      n_cl = 1
      cluster_masses(1) = M_sf
    end if

  end subroutine sample_cmf_clusters

  !------------------------------------------------------------------
  ! Decompose a SF event into clusters drawn from a lognormal kernel
  ! tied to the unresolved turbulent density PDF (subgrid fragmentation
  ! model): ln M_cl ~ N(ln Mcl_char, sigma_lnM^2).
  ! Mcl_char and sigma_lnM are precomputed per SF event (see starform2,
  ! sf_model=1 branch) from the same multi-ff sigs/scrit already used
  ! for the star formation rate.
  ! Mass is exactly conserved: remainder added to last cluster.
  ! mmin/mmax bound both the sampling center Mcl_char and every draw,
  ! matching the caps already enforced by sample_cmf_clusters (and the
  ! range the eta_sn_cluster IMF table is built over).
  ! If M_event < mmin, returns one cluster with the full mass.
  ! All masses in Msun.
  !------------------------------------------------------------------
  subroutine sample_lognormal_clusters(M_event, Mcl_char, sigma_lnM, mmin, mmax, &
                                        seed, cluster_masses, n_cl)
    use random
    real(dp), intent(in)    :: M_event, Mcl_char, sigma_lnM, mmin, mmax
    integer,  intent(inout) :: seed(IRandNumSize)
    real(dp), intent(out)   :: cluster_masses(max_clusters_per_event)
    integer,  intent(out)   :: n_cl

    real(dp)     :: M_rem, mcl, lnMcl_char
    real(kind=8) :: RandNum1, RandNum2, v1, v2, rsq, fac, GaussNum

    lnMcl_char = log(min(max(Mcl_char, mmin), mmax))

    n_cl  = 0
    M_rem = M_event

    do while(M_rem >= mmin .and. n_cl < max_clusters_per_event)
      ! Box-Muller, inlined from random::gaussdev but without its module-level
      ! cache (IGauss/GaussBak): that cache lives outside seed(:), so it is not
      ! restored by the CbC pre-pass's phase-1/phase-2 seed replay and desyncs
      ! the two passes. Drawing (and discarding) a fresh pair every time keeps
      ! the entire RNG state in seed(:), as sample_cmf_clusters already assumes.
      rsq = 0.0d0
      do while(rsq.ge.1.0d0 .or. rsq.le.0.0d0)
        call ranf(seed, RandNum1)
        call ranf(seed, RandNum2)
        v1  = 2.0d0*RandNum1 - 1.0d0
        v2  = 2.0d0*RandNum2 - 1.0d0
        rsq = v1**2 + v2**2
      end do
      fac      = sqrt(-2.0d0*log(rsq)/rsq)
      GaussNum = v2*fac

      mcl = exp(lnMcl_char + sigma_lnM*dble(GaussNum))
      mcl = min(max(mcl, mmin), mmax)
      mcl = min(mcl, M_rem)

      n_cl = n_cl + 1
      cluster_masses(n_cl) = mcl
      M_rem = M_rem - mcl
    end do

    ! Exact mass conservation: remainder into last cluster
    if(n_cl > 0) then
      cluster_masses(n_cl) = cluster_masses(n_cl) + M_rem
    else
      ! M_event below mmin: single cluster with full mass
      n_cl = 1
      cluster_masses(1) = M_event
    end if

  end subroutine sample_lognormal_clusters

end module imf_commons
