module myradex_wrapper

implicit none


logical :: flag_good

integer, parameter :: n_item_column = 20
character(len=64) :: &
    About = 'Author: Fujun Du (fjdu@pmo.ac.cn, fujun.du@gmail.com)'

character(len=128) :: column_names = &
    'iup' //' '// 'ilow' //' '// 'Eup' //' '// 'freq' //' '// 'lam' //' ' &
    // 'Tex' //' '// 'tau' //' '// 'Tr' //' '// &
    'fup' //' '// 'flow' //' '// 'flux_K' //' '// 'flux_int' //' '// 'flux_Jy' //' '// 'beta' //' ' &
    // 'Jnu' //' '// 'gup' //' '// 'glow' //' '// 'Aul' //' '// 'Bul' //' '// 'Blu'

character(len=32) :: molecule_name = ''
contains


subroutine run_one_params( Tkin, &
    Ncol_X_CGS, H2_density_CGS, HI_density_CGS, Hp_density_CGS, E_density_CGS, Tbg, &
    mol_name, data_shape, n_transition, partner_names, colli_shape, &
      level_data, rad_data, colli_T, colli_data, ini_occ,&
      Tb, f_occupation)
  !
  use my_radex
  use statistic_equilibrium
  !
  double precision, intent(in) :: Tkin, Ncol_X_CGS, H2_density_CGS, HI_density_CGS, Hp_density_CGS, E_density_CGS, Tbg
  integer, intent(in) :: n_transition
  integer, dimension(3), intent(in) :: data_shape    ! level number, transition number, partner number
  integer, dimension(:,:), intent(in) :: colli_shape     ! partner number * (colli_transitions , temperature number)
  character(len=8), intent(in) :: mol_name
  character(len=24), intent(in) :: partner_names
  double precision, dimension(:,:), intent(in) :: level_data    !level number * (energy, weight)
  double precision, dimension(:,:), intent(in) :: rad_data      !transition number * (iup, ilow, Aul, freq, Eup, Elow,
                                                                            !lambda, Bul, Blu)
  double precision, dimension(:,:), intent(in) :: colli_T
  double precision, dimension(:,:,:), intent(in) :: colli_data            !partner number * transition number * temperature number
  double precision, dimension(:), intent(in) :: ini_occ
  double precision, dimension(n_transition), intent(out) :: Tb
  double precision, dimension(n_transition+1), intent(out) :: f_occupation

  !
  integer, parameter :: nstr_split = 32
  character(len=32), dimension(nstr_split) :: name_split

  type(type_rad_transition) r
  double precision fup, flow, gup, glow, Tex, Tr, flux_CGS, flux_K_km_s, beam_area, flux_Jy
  double precision Inu_t, tau, t1, t2
  integer i, nout
  !
  rdxx_cfg%geotype = 'lvg'
  rdxx_cfg%nTkin   = 1
  rdxx_cfg%ndv     = 1
  rdxx_cfg%nn_x    = 1
  rdxx_cfg%nNcol_x = 1
  rdxx_cfg%ndens   = 1
  !
  rdxx_cfg%iTkin   = 1
  rdxx_cfg%idv     = 1
  rdxx_cfg%in_x    = 1
  rdxx_cfg%iNcol_x = 1
  rdxx_cfg%idens   = 1
  !
  rdxx_cfg%Tkin(1) = Tkin
  rdxx_cfg%dv(1) = 1D5
  rdxx_cfg%n_x(1) = 1D6
  rdxx_cfg%Ncol_x(1) = Ncol_X_CGS
  !
  rdxx_cfg%n_H2(1) = H2_density_CGS
  rdxx_cfg%n_HI(1) = HI_density_CGS
  rdxx_cfg%n_oH2(1) = 0.0
  rdxx_cfg%n_pH2(1) = 0.0
  rdxx_cfg%n_Hplus(1) = Hp_density_CGS
  rdxx_cfg%n_E(1) = E_density_CGS
  rdxx_cfg%beam_FWHM_in_arcsec=0.4

  rdxx_cfg%nTbg = 1
  rdxx_cfg%Tbg(1) = Tbg
  rdxx_cfg%verbose = .true.
  !
!======load molecule file========
  allocate(a_mol_using)
  a_mol_using%Tkin=Tkin
  a_mol_using%name_molecule=mol_name
  a_mol_using%n_level=data_shape(1)
  allocate(a_mol_using%level_list(a_mol_using%n_level),a_mol_using%f_occupation(a_mol_using%n_level))
  a_mol_using%level_list%energy=level_data(:,1)
  a_mol_using%level_list%weight=level_data(:,2)
  allocate(a_mol_using%rad_data)
  a_mol_using%rad_data%n_transition=data_shape(2)
  allocate(a_mol_using%rad_data%list(a_mol_using%rad_data%n_transition))
  a_mol_using%rad_data%list%iup=int(rad_data(:,1))
  a_mol_using%rad_data%list%ilow=int(rad_data(:,2))
  a_mol_using%rad_data%list%Aul=rad_data(:,3)
  a_mol_using%rad_data%list%freq=rad_data(:,4)
  a_mol_using%rad_data%list%lambda=rad_data(:,5)
  a_mol_using%rad_data%list%Eup=rad_data(:,6)
  a_mol_using%rad_data%list%Elow=rad_data(:,7)
  a_mol_using%rad_data%list%Bul=rad_data(:,8)
  a_mol_using%rad_data%list%Blu=rad_data(:,9)
  allocate(a_mol_using%colli_data)
  a_mol_using%colli_data%n_partner=data_shape(3)
  allocate(a_mol_using%colli_data%list(a_mol_using%colli_data%n_partner))
  call split_str_by_space(partner_names, name_split, nstr_split, nout)
  do i=1, a_mol_using%colli_data%n_partner
    a_mol_using%colli_data%list(i)%name_partner=name_split(i)
    a_mol_using%colli_data%list(i)%n_transition=colli_shape(i,1)
    a_mol_using%colli_data%list(i)%n_T=colli_shape(i,2)
    allocate(a_mol_using%colli_data%list(i)%iup(colli_shape(i,1)), &
             a_mol_using%colli_data%list(i)%ilow(colli_shape(i,1)), &
             a_mol_using%colli_data%list(i)%T_coll(colli_shape(i,2)), &
             a_mol_using%colli_data%list(i)%Cul(colli_shape(i,2),colli_shape(i,1)))
    a_mol_using%colli_data%list(i)%T_coll=colli_T(i,:colli_shape(i,2))
    a_mol_using%colli_data%list(i)%iup=int(colli_data(i,1,:colli_shape(i,1)))
    a_mol_using%colli_data%list(i)%ilow=int(colli_data(i,2,:colli_shape(i,1)))
    a_mol_using%colli_data%list(i)%Cul=colli_data(i,3:colli_shape(i,2)+3,:colli_shape(i,1))
  end do
  a_mol_using%f_occupation=ini_occ

  !======end load molecule file========
  !
  call my_radex_prepare
  call my_radex_prepare_molecule
  call statistic_equil_solve
  call calc_cooling_rate
  !
  flag_good = statistic_equil_params%is_good
  !
  !
  rdxx_cfg%freqmin = 1D-99
  rdxx_cfg%freqmax = 1D99
  !
  do i=1, a_mol_using%rad_data%n_transition
    !associate(r => a_mol_using%rad_data%list(i))
      r = a_mol_using%rad_data%list(i)
      if ((r%freq .lt. rdxx_cfg%freqmin) .or. &
          (r%freq .gt. rdxx_cfg%freqmax)) then
        cycle
      end if
      !
      fup  = a_mol_using%f_occupation(r%iup)
      flow = a_mol_using%f_occupation(r%ilow)
      gup  = a_mol_using%level_list(r%iup)%weight
      glow = a_mol_using%level_list(r%ilow)%weight
      Tex = -(r%Eup - r%Elow) / log(fup*glow / (flow*gup))
      !
      tau = r%tau
      t1 = exp(-tau)
      if (abs(tau) .lt. 1D-6) then
        t2 = tau
      else
        t2 = 1D0 - t1
      end if
      Inu_t = planck_B_nu(Tex, r%freq) * t2 + t1 * r%J_cont_bg
      !
      Tr = (Inu_t - r%J_cont_bg) * phy_SpeedOfLight_CGS**2 / &
        (2D0 * r%freq**2 * phy_kBoltzmann_CGS)
      flux_K_km_s = Tr * a_mol_using%dv / 1D5 * phy_GaussFWHM_c
      flux_CGS = (Inu_t - r%J_cont_bg) * &
        a_mol_using%dv * r%freq / phy_SpeedOfLight_CGS
      beam_area = FWHM_to_area(rdxx_cfg%beam_FWHM_in_arcsec)
      flux_Jy = (Inu_t - r%J_cont_bg) * beam_area / phy_jansky2CGS
      Tb(i)=Tr
      !
    !end associate
  end do
  f_occupation=a_mol_using%f_occupation
  deallocate (a_mol_using)
end subroutine run_one_params


subroutine calc_critical_density(tau, n_partner, n_transitions, &
    critical_densities, iup, ilow)
  use statistic_equilibrium
  double precision, intent(in) :: tau
  integer, intent(in) :: n_partner, n_transitions
  double precision, dimension(n_partner,n_transitions), intent(out) :: critical_densities
  integer, dimension(n_transitions), intent(out) :: iup, ilow
  integer ipt, itr
  call calc_critical_density_f(tau)
  do ipt=1, n_partner
    do itr=1, n_transitions
      critical_densities(ipt, itr) = a_mol_using%rad_data%list(itr)%critical_densities(ipt)
    end do
  end do
  do itr=1, n_transitions
    iup(itr) = a_mol_using%rad_data%list(itr)%iup
    ilow(itr) = a_mol_using%rad_data%list(itr)%ilow
  end do
end subroutine calc_critical_density



subroutine calc_critical_density_old_def(tau, n_partner, n_transitions, &
    critical_densities, iup, ilow)
  use statistic_equilibrium
  double precision, intent(in) :: tau
  integer, intent(in) :: n_partner, n_transitions
  double precision, dimension(n_partner,n_transitions), intent(out) :: critical_densities
  integer, dimension(n_transitions), intent(out) :: iup, ilow
  integer ipt, itr
  call calc_critical_density_old_def_f(tau)
  do ipt=1, n_partner
    do itr=1, n_transitions
      critical_densities(ipt, itr) = a_mol_using%rad_data%list(itr)%critical_densities(ipt)
    end do
  end do
  do itr=1, n_transitions
    iup(itr) = a_mol_using%rad_data%list(itr)%iup
    ilow(itr) = a_mol_using%rad_data%list(itr)%ilow
  end do
end subroutine calc_critical_density_old_def

end module myradex_wrapper
