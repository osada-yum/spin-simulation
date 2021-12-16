      program Ising2d_equilibrium_f77_withfunc
      use utility_m
      use benchmark_m
      implicit real(8) (a-h,o-z)
      parameter(nx=501, ny=500, N=nx*ny, noff=nx, nall=N+2*noff)
      parameter(ilb=-noff+1, irb=ilb+nall)
      parameter(nkbt=100, dkbt_beg=1.7d0, dkbt_end=2.4d0)
      parameter(mcs_relx=1000, mcs_smpl=1000)
      dimension Ising(ilb:irb)
      dimension exparr(-8:8)
      dimension dkbts(nkbt), dmagne(nkbt), energy(nkbt)
      dimension rnd(N)
      type(benchmark_t) :: bm
      write(6, '(a,i0)'     ) "# N = ", N
      write(6, '(a,i0)'     ) "# MCS(relaxation) = ", mcs_relx
      write(6, '(a,i0)'     ) "# MCS(sampling)   = ", mcs_smpl
      write(0, '(3(a, i0))' ) "nx", nx, " ny", ny
      write(0, '(a, f30.18)') "method: METROPOLIS"
      bm = benchmark_t()
c initialize dkbts, Ising, energy, dmagne.
      call linspace(nkbt, dkbt_beg, dkbt_end, dkbts)

      call random_number(rnd)
      call init_Ising(ilb, N, rnd, Ising)

      energy = 0.0d0
      dmagne = 0.0d0
c
c update norishiro.
      do i = 1, noff
         Ising(i-noff) = Ising(i+N-noff)
         Ising(i+N)    = Ising(i)
      end do
c update Ising in each temperatures.
      exparr = 1.0d0
      energy = 0.0d0
      magne  = 0.0d0
      call bm%stamp("start update")
      do k = 1, nkbt
         dbeta = 1/dkbts(k)
         write(0, '(a, i7, es23.15)') "iterate: ", k, dkbts(k)
c        Initialize exparr.
         call init_exparr(dbeta, exparr)
c        relax Ising with Metropolis.
         do j = 1, mcs_relx
c           update Ising with checkerboard pattern.
            call random_number(rnd)
            call Metropolis(ilb, N, nx, rnd, exparr, Ising)
         end do ! end j (1:mcs_relx)
c        relax Ising with Metropolis and calculate parameters.
         dm = 0.0d0
         e  = 0.0d0
         do j = 1, mcs_smpl
c           update Ising with checkerboard pattern.
            call random_number(rnd)
            call Metropolis(ilb, N, nx, rnd, exparr, Ising)
            magne   = 0
            ienergy = 0
            do i = 1, N
               magne   = magne   + Ising(i)
               ienergy = ienergy - Ising(i)*(Ising(i+1)+Ising(i+nx))
            end do
            dm = dm + 1.0d0*magne   / N
            e  = e  + 1.0d0*ienergy / N
         end do ! end j (1:mcs_smpl)
         dmagne(k) = dm / mcs_smpl
         energy(k) = e  / mcs_smpl
      end do ! end k (1:nkbt)
      call bm%stamp("end update")
c print parameters.
      write(6, '(a)') "# temperature, magne, energy"
      do k = 1, nkbt
         write(6, '(i0,a,i0," ",i0," ",*(es20.12))')
     &        nx, "x", ny,
     &        k, dkbts(k), dmagne(k), energy(k)
      end do
      call bm%dump()
      call destroy_benchmark_t(bm)
c print all spins including norishiro
c$$$      call print_Ising(ilb, nx, ny, Ising)
      end program Ising2d_equilibrium_f77_withfunc
c
      subroutine linspace(n, beg, end, arr)
      implicit real(8) (a-h,o-z)
      dimenstion arr(n)
c internal division point.
      do i = 1, n
         arr(i) =
     &        ( (i-1)       *beg
     &        + (n-1-i+1)*end ) / (n-1)
      end do
      end subroutine
c
      subroutine init_Ising(ilb, n, rnd, I)
      dimension rnd(n), Ising(ilb:)
      do i = 1, N
         if (rnd(i) < 0.5d0) then
            I(i) = +1
         else
            I(i) = -1
         end if
      end do
c$$$      call print_Ising(ilb, nx, ny, Ising)
      end subroutine
c
      subroutine init_exparr(dbeta, exparr)
      dimension exparr(-8:8)
      do i = 1, 8
         exparr(i) = exp(-dbeta*i)
      end do
      end subroutine
c
      subroutine Metropolis(ilb, n, nx, rnd, exparr, Ising)
      dimension rnd(n), exparr(-8:8), Ising(ilb:)
      do is = 1, 2
         do i = is, n, 2
            ide = 2*Ising(i)*
     &           (Ising(i+1)+Ising(i-1)+Ising(i+nx)+Ising(i-nx))
            if (rnd(i) .lt. exparr(ide)) Ising(i) = -Ising(i)
         end do
c        update norishiro.
         call norishiro(ilb, n, noff, Ising)
      end do
      end subroutine
c
      subroutine norishiro(ilb, n, noff, Ising)
      dimension Ising(ilb:)
      do i = 1, noff
         Ising(i-noff) = Ising(i+n-noff)
         Ising(i+n)    = Ising(i)
      end do
      end subroutine
c
      subroutine print_Ising(ilb, nx, ny, Ising)
      dimension Ising(ilb:)
      do j = 0, ny+1
         write(0, '(i4, a)', advance="NO") j, ": "
         do i = 1, nx
            write(0, '(i3)', advance="NO") Ising((j-1)*nx+i)
         end do
         write(0, *)
      end do
      end subroutine
