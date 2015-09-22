

      
      
      module geom_module
      implicit none
      !
      !integer, parameter :: nx = 64, ny = 64, nz = 64
      !integer, parameter :: nx = 50, ny = 50, nz = 50
      integer :: nx = 0, ny = 0, nz = 0 
      
      end module geom_module


      
      
      program main
      !
      use const_module
      use lbm3d_module
      use geom_module
      !
      implicit none
      ! 
      double precision, parameter :: G = -1.2d0 ! molecular interaction 
      double precision, parameter :: omega1 = 1.0d0, omega2 = 1.0d0
      !
      double precision, parameter :: drho = 1.0d-3
      !
      double precision, allocatable :: fIn(:,:,:,:), gIn(:,:,:,:)
      double precision, allocatable :: fOut(:,:,:,:), gOut(:,:,:,:)
      double precision, allocatable :: rho1(:,:,:), rho2(:,:,:)
      double precision, allocatable :: jx1(:,:,:), jy1(:,:,:), jz1(:,:,:)
      double precision, allocatable :: jx2(:,:,:), jy2(:,:,:), jz2(:,:,:)
      double precision :: fEq, gEq

      !
      double precision :: rhotot_omega
      double precision :: utotx, utoty, utotz
      double precision :: rhoContrib1x, rhoContrib1y, rhoContrib1z
      double precision :: rhoContrib2x, rhoContrib2y, rhoContrib2z
      double precision :: utotx1, utoty1, utotz1
      double precision :: utotx2, utoty2, utotz2
      double precision :: cu1, cu2
      !
      double precision, allocatable :: xh(:), yh(:), zh(:)
      
      !
      integer :: nstep = 100000
      integer :: istep = 0
      integer :: plot_int = 100
      integer :: iplot = 0
      double precision :: time = 0, dt = 1.0d0
      
      !
      integer :: domlo(nd), domhi(nd)
      integer :: i,j,k
      integer :: l
      integer :: iup,jup,kup
      double precision :: tmp
      double precision :: tmpvec(nq)
      
      ! for initialization with random number
      integer :: seed_size
      
      !
      integer :: untin
      namelist /fortin/ &
      nx, ny, nz, &
      nstep, plot_int
      
      ! read namelists
      untin = 9
      open(unit=untin, file='prob.fortin', form='formatted', status='old')
      read(unit=untin, nml=fortin)
      close(unit=untin)
      
      
      if (nx.le.0 .or. ny.le.0 .or. nz.le.0) then
        print *, 'Invalid cell number'
        stop
      endif
      
      !
      allocate(fIn(nq,nx+2,ny+2,nz+2), gIn(nq,nx+2,ny+2,nz+2))
      allocate(fOut(nq,nx+2,ny+2,nz+2), gOut(nq,nx+2,ny+2,nz+2))
      
      allocate(rho1(nx+2,ny+2,nz+2), rho2(nx+2,ny+2,nz+2))
      allocate(jx1(nx+2,ny+2,nz+2), jy1(nx+2,ny+2,nz+2), jz1(nx+2,ny+2,nz+2))
      allocate(jx2(nx+2,ny+2,nz+2), jy2(nx+2,ny+2,nz+2), jz2(nx+2,ny+2,nz+2))
      
      allocate(xh(nx+1), yh(ny+1), zh(nz+1))
      
      domlo(1) = 2
      domhi(1) = nx+1
      domlo(2) = 2
      domhi(2) = ny+1
      domlo(3) = 2
      domhi(3) = nz+1
      
      
      !
      do i = 1, nx+1
        xh(i) = 1.0d0/nx * (i - 1)
      enddo
      do j = 1, ny+1
        yh(j) = 1.0d0/ny * (j - 1)
      enddo
      do k = 1, nz+1
        zh(k) = 1.0d0/nz * (k - 1)
      enddo
      
      ! initial condition for F and G distribution functions
      if (.true.) then 
        call random_seed(size=seed_size)
        print *, 'seed_size=', seed_size
        call random_seed()
      endif
      do k = domlo(3), domhi(3)
      do j = domlo(2), domhi(2)
      do i = domlo(1), domhi(1)
        call random_number(tmp)
        tmp = drho * (2.0d0*tmp - 1.0d0)
        fIn(:,i,j,k) = qw(:) * (1.0d0 + tmp)
        gIn(:,i,j,k) = qw(:) * (1.0d0 - tmp)
      enddo
      enddo
      enddo
      ! fill BC
      call fill_pbc(fIn, nx,ny,nz)
      call fill_pbc(gIn, nx,ny,nz)
      ! derive macro variables
      call derive_macroscopic(fIn, rho1,jx1,jy1,jz1, nx,ny,nz)
      call derive_macroscopic(gIn, rho2,jx2,jy2,jz2, nx,ny,nz)
      
      ! initial output
      time = 0.0d0
      istep = 0
      iplot = 0
      call output_vtk3d(iplot,time,&
      rho1, &
      xh,yh,zh, &
      nx,ny,nz)
      iplot = iplot + 1
            
      ! main loop
      do istep = 1, nstep
        print *, 'step=',istep, 'time=',time
        
        !
        do k = domlo(3), domhi(3)
        do j = domlo(2), domhi(2)
        do i = domlo(1), domhi(1)
          !
          rhotot_omega = rho1(i,j,k)*omega1 + rho2(i,j,k)*omega2
          utotx = (jx1(i,j,k)*omega1 + jx2(i,j,k)*omega2) / rhotot_omega
          utoty = (jy1(i,j,k)*omega1 + jy2(i,j,k)*omega2) / rhotot_omega
          utotz = (jz1(i,j,k)*omega1 + jz2(i,j,k)*omega2) / rhotot_omega
          
          !
          rhoContrib1x = 0.0d0
          rhoContrib1y = 0.0d0
          rhoContrib1z = 0.0d0
          rhoContrib2x = 0.0d0
          rhoContrib2y = 0.0d0
          rhoContrib2z = 0.0d0
          do l = 2, nq
            iup = i - qex(l)
            jup = j - qey(l)
            kup = k - qez(l)
            rhoContrib1x = rhoContrib1x + rho1(iup,jup,kup)*qw(l)*qex(l)
            rhoContrib1y = rhoContrib1y + rho1(iup,jup,kup)*qw(l)*qey(l)
            rhoContrib1z = rhoContrib1z + rho1(iup,jup,kup)*qw(l)*qez(l)
            rhoContrib2x = rhoContrib2x + rho2(iup,jup,kup)*qw(l)*qex(l)
            rhoContrib2y = rhoContrib2y + rho2(iup,jup,kup)*qw(l)*qey(l)
            rhoContrib2z = rhoContrib2z + rho2(iup,jup,kup)*qw(l)*qez(l)
          enddo
          
          ! potential contribution
          utotx1 = utotx - G/omega1*rhoContrib2x
          utoty1 = utoty - G/omega1*rhoContrib2y
          utotz1 = utotz - G/omega1*rhoContrib2z
          !
          utotx2 = utotx - G/omega2*rhoContrib1x
          utoty2 = utoty - G/omega2*rhoContrib1y
          utotz2 = utotz - G/omega2*rhoContrib1z
          
          ! collision
          do l = 1,nq
            cu1 = 3.0d0 * (qex(l)*utotx1 + qey(l)*utoty1 + qez(l)*utotz1)
            cu2 = 3.0d0 * (qex(l)*utotx2 + qey(l)*utoty2 + qez(l)*utotz2)
            
            fEq = rho1(i,j,k) * qw(l) * &
            (1.0d0 + cu1 + 0.5d0*cu1*cu1 - 1.5d0*(utotx1**2+utoty1**2+utotz1**2)) 
            gEq = rho2(i,j,k) * qw(l) * &
            (1.0d0 + cu2 + 0.5d0*cu2*cu2 - 1.5d0*(utotx2**2+utoty2**2+utotz2**2))
            
            fOut(l,i,j,k) = fIn(l,i,j,k) - omega1 * (fIn(l,i,j,k)-fEq)
            gOut(l,i,j,k) = gIn(l,i,j,k) - omega2 * (gIn(l,i,j,k)-gEq)
          enddo
        enddo
        enddo
        enddo
        
        !
        call fill_pbc(fOut, nx,ny,nz)
        call fill_pbc(gOut, nx,ny,nz)
        do k = domlo(3), domhi(3)
        do j = domlo(2), domhi(2)
        do i = domlo(1), domhi(1)
          do l = 1,nq
            iup = i - qex(l)
            jup = j - qey(l)
            kup = k - qez(l)
            
            fIn(l,i,j,k) = fOut(l,iup,jup,kup)
            gIn(l,i,j,k) = gOut(l,iup,jup,kup)
          enddo
        enddo
        enddo
        enddo
        
        !
        call fill_pbc(fIn, nx,ny,nz)
        call fill_pbc(gIn, nx,ny,nz)
        call derive_macroscopic(fIn, rho1,jx1,jy1,jz1, nx,ny,nz)
        call derive_macroscopic(gIn, rho2,jx2,jy2,jz2, nx,ny,nz)
        
        time = time + dt
        
        if (mod(istep,plot_int)==0 .or. istep==nstep) then
          call output_vtk3d(iplot,time,&
          rho1, &
          xh,yh,zh, &
          nx,ny,nz)
          iplot = iplot + 1
        endif
      enddo ! main loop
      
      
      end program main
      
      subroutine fill_pbc(phi, nx,ny,nz)
      use const_module
      use lbm3d_module
      implicit none
      integer, intent(in) :: nx, ny, nz
      double precision, intent(inout) :: phi(nq,nx+2,ny+2,nz+2)
      !
      integer :: i,j,k
      
      ! x
      do k = 1, nz+2
      do j = 1, ny+2
        phi(:,1,j,k) = phi(:,nx+1,j,k)
        phi(:,nx+2,j,k) = phi(:,2,j,k)
      enddo
      enddo
      ! y
      do k = 1, nz+2
      do i = 1, nx+2
        phi(:,i,1,k) = phi(:,i,ny+1,k)
        phi(:,i,ny+2,k) = phi(:,i,2,k)
      enddo
      enddo
      !z
      do j = 1, ny+2
      do i = 1, nx+2
        phi(:,i,j,1) = phi(:,i,j,nz+1)
        phi(:,i,j,nz+2) = phi(:,i,j,2)
      enddo
      enddo
      
      return
      end subroutine fill_pbc
      
      subroutine derive_macroscopic(f, rho,jx,jy,jz, nx,ny,nz)
      use const_module
      use lbm3d_module
      implicit none
      integer, intent(in) :: nx, ny, nz
      double precision, intent(in) :: f(nq,nx+2,ny+2,nz+2)
      double precision, intent(out) :: rho(nx+2,ny+2,nz+2)
      double precision, intent(out) :: jx(nx+2,ny+2,nz+2), jy(nx+2,ny+2,nz+2), jz(nx+2,ny+2,nz+2)
      !
      integer :: i,j,k
      
      do k = 1, nz+2
      do j = 1, ny+2
      do i = 1, nx+2
        rho(i,j,k) = sum(f(:,i,j,k))
        jx(i,j,k) = sum(qex * f(:,i,j,k))
        jy(i,j,k) = sum(qey * f(:,i,j,k))
        jz(i,j,k) = sum(qez * f(:,i,j,k))
      enddo
      enddo
      enddo
      
      return
      endsubroutine derive_macroscopic
      
      
