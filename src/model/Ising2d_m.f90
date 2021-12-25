! skew boundary condition Square lattice Ising model.
module Ising2d_m
  use, intrinsic :: iso_fortran_env
  use utility_m
  use random_base_m
  implicit none

  integer, parameter, private :: num_neighbors = 4                          ! 最近接格子数.
  integer, parameter, private :: J_interaction = 1.0
  integer, parameter, private :: array_discrete_energy_diff(-8:8) = &
       (/-8, -7, -6, -5, -4, -3, -2, -1, 0&
       ,  1,  2,  3,  4,  5,  6,  7,  8/)

  type Ising2d
     private
     integer(ikind)              :: nx_s, ny_s
     real(rkind)                 :: kbt_s, beta_s
     integer(ikind)              :: offset_s, particles_s
     integer(ikind)              :: begin_s, end_s, all_s
     integer(ikind)              :: neighbors_indices_array_s(num_neighbors)
     integer(ikind), allocatable :: spin_s(:)
     real(rkind)                 :: ising_exp_s(-8:8)                        ! 更新確率の配列(使うのは5パターン-4,-2,0,+2,+4).
     procedure(updater), pointer :: updater_s => null()
   contains
     ! getter
     procedure,pass :: x                 => get_x
     procedure,pass :: y                 => get_y
     procedure,pass :: offset            => get_offset                       ! のりしろ用のoffset.
     procedure,pass :: particles         => get_particles                    ! 全体の格子数.
     procedure,pass :: kbt               => get_kbt                          ! J/T.
     procedure,pass :: beta              => get_beta                         ! 逆温度.
     procedure,pass :: indices_neighbors => get_indices_neighbors            ! 周りの格子点のindexの配列を返す.
     procedure,pass :: spin              => get_spin
     procedure,pass :: spin_row          => get_spin_row
     ! setter
     procedure,pass :: set_kbt           => set_kbt_s
     procedure,pass :: set_beta          => set_beta_s
     procedure,pass :: set_order_spin    => set_order_spin_ising
     procedure,pass :: set_random_spin   => set_random_spin_ising
     procedure,pass :: norishiro         => norishiro_ising
     ! updater
     procedure,pass :: set_updater                    => set_updater_ising
     procedure,pass :: update_one_mcs                 => update_one_mcs_ising
     procedure,pass :: update_with_Metropolis_one_mcs => update_Metropolis_one_mcs_ising
  end type Ising2d

  private :: init_Ising2d
  ! デフォルトコンストラクタ.
  interface Ising2d
     module procedure init_Ising2d
  end interface Ising2d

  interface
     !! updater: update method.
     subroutine updater(this, rng)
       import Ising2d
       import random_base_t
       class(Ising2d)      , intent(inout) :: this
       class(random_base_t), intent(inout) :: rng
     end subroutine updater
  end interface

contains
  ! デフォルトコンストラクタ.
  !! init_Ising2d: kbtとx, yを引数で取る.
  !! 逆温度や総粒子数, 始まりと終わりのインデックス, 最近接格子までのインデックスの距離と更新確率を初期化する.
  impure type(Ising2d) function init_Ising2d(kbt,x,y) result(res_sp)
    real(rkind), intent(in) :: kbt
    integer    , intent(in) :: x, y
    call res_sp%set_updater("Metropolis")
    res_sp%kbt_s  = kbt
    res_sp%beta_s = 1.0_rkind/kbt
    res_sp%nx_s   = x
    res_sp%ny_s   = y

    res_sp%particles_s = res_sp%nx_s * res_sp%ny_s                                             ! 総粒子数.
    res_sp%offset_s    = res_sp%nx_s                                                           ! 上と下ののりしろ分.
    res_sp%begin_s     = res_sp%offset_s   + 1                                                 ! 始まりのインデックス.
    res_sp%end_s       = res_sp%offset_s   + res_sp%particles_s                                ! 終わりのインデックス.
    res_sp%all_s       = res_sp%offset_s*2 + res_sp%particles_s                                ! 配列の大きさ.
    allocate(res_sp%spin_s( (res_sp%begin_s-res_sp%offset_s):(res_sp%end_s+res_sp%offset_s) )) ! 総粒子数とのりしろの領域をとる.
    if (.not. allocated(res_sp%spin_s)) then
       call util_error_stop("res_sp%spin_s is not allocated."&
            , __LINE__, __FILE__)
    end if
    ! 遷移確率の配列.
    res_sp%ising_exp_s = exp( -res_sp%beta_s * array_discrete_energy_diff)
    ! 最近接格子までのインデックスの距離の初期化.
    res_sp%neighbors_indices_array_s(1) = +1
    res_sp%neighbors_indices_array_s(2) = -1
    res_sp%neighbors_indices_array_s(3) = +res_sp%nx_s
    res_sp%neighbors_indices_array_s(4) = -res_sp%nx_s
  end function init_Ising2d

  ! getters
  pure elemental integer(ikind) function get_x(this) result(res_i)
    class(Ising2d), intent(in) :: this
    res_i = this%nx_s
  end function get_x
  pure elemental integer(ikind) function get_y(this) result(res_i)
    class(Ising2d), intent(in) :: this
    res_i = this%ny_s
  end function get_y
  pure elemental real(rkind) function get_kbt(this) result(res_r)
    class(Ising2d), intent(in) :: this
    res_r = this%kbt_s
  end function get_kbt
  pure elemental real(rkind) function get_beta(this) result(res_r)
    class(Ising2d), intent(in) :: this
    res_r = this%beta_s
  end function get_beta
  pure elemental integer(ikind) function get_offset(this) result(res_i)
    class(Ising2d), intent(in) :: this
    res_i = this%offset_s
  end function get_offset
  pure elemental integer(ikind) function get_particles(this) result(res_i)
    class(Ising2d), intent(in) :: this
    res_i = this%particles_s
  end function get_particles
  ! 最近接格子の配列を返す関数. 要素数はnum_neighbors... (配列返さなくても良さそう...)
  impure function get_indices_neighbors(this,index)
    class(Ising2d), intent(in) :: this
    integer       , intent(in) :: index
    integer(ikind)             :: get_indices_neighbors(num_neighbors) ! 最近接格子の配列
    get_indices_neighbors(1:num_neighbors) = index + this%neighbors_indices_array_s(1:num_neighbors)
  end function get_indices_neighbors
  !! get_spin: begin()からのindexのスピンを返す. [begin_s, end_s].
  pure integer(ikind) function get_spin(this,index) result(res_p)
    class(Ising2d), intent(in) :: this
    integer       , intent(in) :: index
    res_p = this%spin_s(this%begin_s -1 + index)
  end function get_spin
  !! get_spin_row: 生のspinを返す. [begin_s-offset_s, end_s+offset_s].
  integer(ikind) function get_spin_row(this,index) result(res_p)
    class(Ising2d), intent(in) :: this
    integer       , intent(in) :: index
    res_p = this%spin_s(index)
  end function get_spin_row

  ! setter
  subroutine set_kbt_s(this,kbt)
    class(Ising2d), intent(inout) :: this
    real(rkind)   , intent(in)    :: kbt
    this%kbt_s  = kbt
    this%beta_s = 1.0_rkind/kbt
    this%ising_exp_s = exp( -this%beta_s * array_discrete_energy_diff)
  end subroutine set_kbt_s
  subroutine set_beta_s(this,beta)
    class(Ising2d), intent(inout) :: this
    real(rkind)   , intent(in)    :: beta
    this%kbt_s  = 1.0_rkind/beta
    this%beta_s = beta
    this%ising_exp_s = exp( -this%beta_s * array_discrete_energy_diff)
  end subroutine set_beta_s
  ! initializer of system.
  !! set_order_spin_ising: 秩序状態になるようにスピンを配置する.
  subroutine set_order_spin_ising(this)
    class(Ising2d), intent(inout) :: this
    this%spin_s(:) = 1
  end subroutine set_order_spin_ising
  !! set_random_spin_ising: ランダムにスピンを配置する.
  subroutine set_random_spin_ising(this, rand_gen)
    class(Ising2d)      , intent(inout) :: this
    class(random_base_t), intent(inout) :: rand_gen
    integer                             :: i
    real(rkind)                         :: rnd(this%begin_s:this%end_s)
    call rand_gen%random_arr(rnd)
    do i = this%begin_s, this%end_s
       if (rnd(i) > 0.5_rkind) then
          this%spin_s(i) = 1
       else
          this%spin_s(i) = -1
       end if
    end do
    call this%norishiro()
  end subroutine set_random_spin_ising
  !! norishiro_ising: のりしろの更新.
  subroutine norishiro_ising(this)
    class(Ising2d), intent(inout) :: this
    integer :: nori_up_start, most_down, nori_down_start, most_up
    nori_up_start   = this%end_s   + 1
    most_down       = this%begin_s
    nori_down_start = this%begin_s - this%offset_s
    most_up         = this%end_s   - this%offset_s + 1
    this%spin_s       (nori_up_start   : nori_up_start   + this%offset_s-1)&
         = this%spin_s(most_down       : most_down       + this%offset_s-1)
    this%spin_s       (nori_down_start : nori_down_start + this%offset_s-1)&
         = this%spin_s(most_up         : most_up         + this%offset_s-1)
  end subroutine norishiro_ising
  ! updater
  !! set_updater_ising: set `system%updater_s` by argument `method`.
  subroutine set_updater_ising(system, method)
    class(Ising2d)  , intent(inout) :: system
    character(len=*), intent(in)    :: method
    select case (method)
    case("Metropolis"); system%updater_s => update_Metropolis_one_mcs_ising
    case default
       call util_error_stop("unknown method: "//method&
            , __LINE__, __FILE__)
    end select
  end subroutine set_updater_ising
  !! update_one_mcs_ising: 系をupdater_sで一回更新する.
  subroutine update_one_mcs_ising(system, rand_gen)
    class(Ising2d)      , intent(inout) :: system
    class(random_base_t), intent(inout) :: rand_gen
    call system%updater_s(rand_gen)
  end subroutine update_one_mcs_ising

  !! update_Metropolis_one_mcs_ising: 系をMetropolis法で 1MCS だけ更新する.
  subroutine update_Metropolis_one_mcs_ising(system, rand_gen)
    class(Ising2d)      , intent(inout) :: system
    class(random_base_t), intent(inout) :: rand_gen
    integer(ikind)                      :: energy_diff
    real(real64)                        :: rnd(system%begin_s:system%end_s)
    integer                             :: i, j
    call rand_gen%random_arr(rnd)
    do j = 0, 1
       do i = j+system%begin_s  , system%end_s, 2
          energy_diff = energy_onespin(i, system%neighbors_indices_array_s(:), system%spin_s(:))
          !! Metropolis法の遷移確率に従ってスピンを反転させる.
          if ( rnd(i) < system%ising_exp_s(energy_diff) ) then
             system%spin_s(i) = - system%spin_s(i)
          end if
       end do
       call system%norishiro()
    end do
  end subroutine update_Metropolis_one_mcs_ising
  !! energy_onespin: 局所的な格子のエネルギー
  pure integer(ikind) function energy_onespin(index, neighbor_indices, spins)
    integer       , intent(in) :: index, neighbor_indices(:)
    integer(ikind), intent(in) :: spins(:)
    energy_onespin = 2*J_interaction * spins(index) *&
         ( spins(index + neighbor_indices(1))&
         + spins(index + neighbor_indices(2))&
         + spins(index + neighbor_indices(3))&
         + spins(index + neighbor_indices(4)))
  end function energy_onespin
  !! calc_magne: 磁化の計算, sum()でOK.
  pure real(rkind) function calc_magne(system) result(magne)
    type(Ising2d), intent(in) :: system
    magne = real(sum(system%spin_s(system%begin_s:system%end_s)), rkind) / system%particles()
  end function calc_magne
  !! calc_energy: エネルギーの計算, ループ回す.
  pure real(rkind) function calc_energy(system) result(energy)
    type(Ising2d), intent(in) :: system
    real(rkind)               :: energy_tmp
    integer                   :: i
    energy_tmp = 0.0_rkind
    do i = 1, system%particles()
       energy_tmp   = energy_tmp &
            - J_interaction*system%spin(i) * &
            ( system%spin(i + system%neighbors_indices_array_s(1))&
            + system%spin(i + system%neighbors_indices_array_s(3)))
    end do
    energy = real(energy_tmp, rkind) / system%particles_s
    energy = real(energy_tmp, rkind) / system%particles()
  end function calc_energy
  !! calc_magne_and_energy: 両方計算する, 1回のループで両方計算する.
  pure subroutine calc_magne_and_energy(system, magne, energy)
    type(Ising2d), intent(in)  :: system
    real(rkind)  , intent(out) :: magne
    real(rkind)  , intent(out) :: energy
    integer(ikind)             :: magne_tmp
    real(rkind)                :: energy_tmp
    integer                    :: i
    magne_tmp  = 0.0_rkind
    energy_tmp = 0.0_rkind
    do i = 1, system%particles()
       magne_tmp  = magne_tmp + system%spin(i)

       energy_tmp = energy_tmp &
            - J_interaction*system%spin(i) * &
            ( system%spin(i + system%neighbors_indices_array_s(1))&
            + system%spin(i + system%neighbors_indices_array_s(3)))
    end do
    magne  = real(magne_tmp, rkind)  / system%particles()
    energy = real(energy_tmp, rkind) / system%particles()
    return
  end subroutine calc_magne_and_energy
end module Ising2d_m
