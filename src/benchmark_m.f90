module benchmark_m
  use, intrinsic :: iso_fortran_env, only : int64, error_unit
  use utility_m, only : util_debug_print, util_error_stop
  implicit none
  private

  type :: stamp_time_list_t
     private
     character(len=128)               :: stamp
     real                             :: time
     type(stamp_time_list_t), pointer :: next => null()
   contains
     procedure, pass :: head_stmp => head_stmp_st_lst
     procedure, pass :: head_time => head_time_st_lst
     procedure, pass :: tail      => tail_st_lst
     procedure, pass :: add_front => add_front_st, add_front_lst
     procedure, pass :: dump      => dump_st_lst
     final           :: destroy_st_lst
  end type stamp_time_list_t

  ! default constructor.
  interface stamp_time_list_t
     module procedure init_stamp_time_list_t
  end interface Stamp_Time_List_T

  public :: benchmark_t, destroy_benchmark_t
  type :: benchmark_t
     private
     integer(int64)                   :: time_beg_cnt, time_end_cnt, cnt_per_sec, cnt_max
     type(stamp_time_list_t), pointer :: stamp_lst => null()
   contains
     procedure, pass :: stamp => stamp_counter_bench
     procedure, pass :: dump  => dump_counter_bench
     final           :: destroy_benchmark_t
  end type benchmark_t

  ! default constructor.
  interface benchmark_t
     module procedure init_benchmark_t
  end interface Benchmark_T

contains

  impure function init_stamp_time_list_t(stamp, time) result(st_lst)
    !! `init_stamp_time_list_t`: initialize with `stamp` and `time`.
    character(len=*), intent(in)     :: stamp
    real            , intent(in)     :: time
    type(stamp_time_list_t), pointer :: st_lst
    allocate(st_lst)
    if (.not. associated(st_lst)) then
       call util_error_stop("failed: init of stamp_time_list_t"&
            , __LINE__, __FILE__)
    end if
    st_lst%stamp = stamp
    st_lst%time  = time
    nullify(st_lst%next)
  end function init_stamp_time_list_t

  pure function head_stmp_st_lst(this) result(stmp)
    !! `head_stmp_st_lst`: return `stamp` in head of list.
    class(stamp_time_list_t), intent(in) :: this
    character(len=128)                   :: stmp
    stmp = this%stamp
  end function head_stmp_st_lst
  pure real function head_time_st_lst(this) result(time)
    !! `head_time_st_lst`: return `time` in head of list.
    class(stamp_time_list_t), intent(in) :: this
    time = this%time
  end function head_time_st_lst
  function tail_st_lst(this) result(lst)
    !! `tail_st_lst`: return pointer to `this%next`.
    class(stamp_time_list_t), intent(in) :: this
    class(stamp_time_list_t), pointer    :: lst
    lst => this%next
  end function tail_st_lst

  subroutine add_front_st(this, stamp, time)
    !! `add_front_st`: use `stamp` and `time` to add.
    class(stamp_time_list_t), target , intent(inout) :: this
    character(len=*)                 , intent(in)    :: stamp
    real                             , intent(in)    :: time
    class(stamp_time_list_t), pointer                :: lst_tmp
    lst_tmp      => stamp_time_list_t(this%stamp, this%time)
    lst_tmp%next => this%next
    this%next    => lst_tmp
    this%stamp   =  stamp
    this%time    =  time
    nullify(lst_tmp)
  end subroutine add_front_st
  subroutine add_front_lst(this, lst)
    !! `add_front_lst`: use list of `stamp` and `time` to add.
    class(stamp_time_list_t), intent(inout) :: this
    class(stamp_time_list_t), intent(in)    :: lst
    class(stamp_time_list_t), pointer       :: lst_tmp
    lst_tmp       => stamp_time_list_t(this%stamp, this%time)
    lst_tmp%next  => this%next
    this%next     => lst_tmp
    this%stamp    =  lst%stamp
    this%time     =  lst%time
    nullify(lst_tmp)
  end subroutine add_front_lst

  recursive subroutine dump_st_lst(this)
    !! `dump_st_lst`: print all elements of list in reverse (this is in order in case of only use of `add_front`).
    class(stamp_time_list_t), target, intent(in) :: this
    class(stamp_time_list_t), pointer            :: lst
    lst => this
    if (.not. associated(lst)) then
       return
    end if
    call lst%next%dump()
    write(error_unit, '(a, f8.3, a)') trim(lst%stamp), lst%time, " sec"
  end subroutine dump_st_lst

  subroutine destroy_st_lst(lst)
    !! `destroy_st_lst`: deallocate all elements of list.
    type(stamp_time_list_t), intent(inout) :: lst
    type(stamp_time_list_t), pointer       :: lst_tmp
    lst_tmp => lst%next
    call util_debug_print("finalize stamp_time_list_t"&
         , __LINE__, __FILE__)
    if (associated(lst_tmp)) then
       deallocate(lst_tmp)
       nullify(lst_tmp)
    end if
  end subroutine destroy_st_lst

  impure type(benchmark_t) function init_benchmark_t() result(bench)
    !! `init_conuter_bench`: initialize benchmark_t.
    call system_clock(bench%time_beg_cnt, bench%cnt_per_sec, bench%cnt_max)
    bench%time_end_cnt =  bench%time_beg_cnt
    bench%stamp_lst    => stamp_time_list_t("init", 0.0)
  end function init_benchmark_t

  subroutine stamp_counter_bench(this, message)
    !! `stamp_conuter_bench`: output elapsed seconds from last stamp to current stamp.
    !! call bench%stamp("start")
    !! (some sequence...)
    !! call bench%stamp("stop")
    class(benchmark_t), intent(inout) :: this
    character(len=*)  , intent(in)    :: message
    real                              :: elapsed_sec
    call system_clock(this%time_end_cnt)
    elapsed_sec = real(this%time_end_cnt-this%time_beg_cnt)/real(this%cnt_per_sec)
    write(error_unit, '(a, f8.3, a)')  message&
         , elapsed_sec, " sec"
    this%time_beg_cnt =  this%time_end_cnt
    call this%stamp_lst%add_front(message, elapsed_sec)
  end subroutine stamp_counter_bench

  subroutine dump_counter_bench(this)
    class(benchmark_t), intent(in) :: this
    call this%stamp_lst%dump()
  end subroutine dump_counter_bench

  subroutine destroy_benchmark_t(bm)
    !! `destroy_benchmark_t`: deallocate `stamp_lst`.
    type(benchmark_t), intent(inout) :: bm
    if (associated(bm%stamp_lst)) then
       deallocate(bm%stamp_lst)
    end if
    call util_debug_print("finalize benchmark_t"&
         , __LINE__, __FILE__)
  end subroutine destroy_benchmark_t

end module benchmark_m
