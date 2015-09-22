

      
      
      module geom_module
      implicit none
      !
      integer, parameter :: nx = 200, ny = 200
      end module geom_module


      subroutine fill_pbc(phi, nx,ny)
      use const_module
      use lbm2d_module
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(inout) :: phi(nq,nx+2,ny+2)
      !
      integer :: i,j
      
      ! x
      do j = 1, ny+2
        phi(:,1,j) = phi(:,nx+1,j)
        phi(:,nx+2,j) = phi(:,2,j)
      enddo
      ! y
      do i = 1, nx+2
        phi(:,i,1) = phi(:,i,ny+1)
        phi(:,i,ny+2) = phi(:,i,2)
      enddo
      
      return
      end subroutine fill_pbc
      
      subroutine derive_macroscopic(f, rho, jx, jy, nx,ny)
      use const_module
      use lbm2d_module
      implicit none
      integer, intent(in) :: nx, ny
      double precision, intent(in) :: f(nq,nx+2,ny+2)
      double precision, intent(out) :: rho(nx+2,ny+2)
      double precision, intent(out) :: jx(nx+2,ny+2), jy(nx+2,ny+2)
      !
      integer :: i,j
      !integer :: l
      
      do j = 1, ny+2
      do i = 1, nx+2
        rho(i,j) = sum(f(:,i,j))
        jx(i,j) = sum(qex * f(:,i,j))
        jy(i,j) = sum(qey * f(:,i,j))
      enddo
      enddo
      
      return
      endsubroutine derive_macroscopic
      
      program main2d
      !
      use const_module
      use lbm2d_module
      use geom_module
      !
      implicit none
      ! 
      double precision, parameter :: G = -1.2d0 ! molecular interaction 
      double precision, parameter :: omega1 = 1.0d0, omega2 = 1.0d0
      !
      double precision, parameter :: drho = 1.0d-3
      !
      double precision, allocatable :: fIn(:,:,:), gIn(:,:,:)
      double precision, allocatable :: fOut(:,:,:), gOut(:,:,:)
      double precision, allocatable :: rho1(:,:), rho2(:,:)
      double precision, allocatable :: jx1(:,:), jy1(:,:)
      double precision, allocatable :: jx2(:,:), jy2(:,:)
      double precision :: fEq, gEq

      !
      double precision :: rhotot_omega
      double precision :: utotx, utoty
      double precision :: rhoContrib1x, rhoContrib1y
      double precision :: rhoContrib2x, rhoContrib2y
      double precision :: utotx1, utoty1
      double precision :: utotx2, utoty2
      double precision :: cu1, cu2
      !
      double precision, allocatable :: xh(:), yh(:)
      
      !
      integer :: nstep = 100000
      integer :: istep = 0
      integer :: plot_int = 100
      integer :: iplot = 0
      double precision :: time = 0, dt = 1.0d0
      
      !
      integer :: domlo(nd), domhi(nd)
      integer :: i,j
      integer :: l
      integer :: iup,jup
      double precision :: tmp
      double precision :: tmpvec(nq)
      
      domlo(1) = 2
      domhi(1) = nx+1
      domlo(2) = 2
      domhi(2) = ny+1
      
      !
      allocate(fIn(nq,nx+2,ny+2), gIn(nq,nx+2,ny+2))
      allocate(fOut(nq,nx+2,ny+2), gOut(nq,nx+2,ny+2))
      
      allocate(rho1(nx+2,ny+2), rho2(nx+2,ny+2))
      allocate(jx1(nx+2,ny+2), jy1(nx+2,ny+2))
      allocate(jx2(nx+2,ny+2), jy2(nx+2,ny+2))
      
      allocate(xh(nx+1), yh(ny+1))
      
      !
      do i = 1, nx+1
        xh(i) = i - 1
      enddo
      do j = 1, ny+1
        yh(j) = j - 1
      enddo
      
      ! initial condition for F and G distribution functions
      if (.true.) then 
        call random_seed()
      endif
      do j = domlo(2), domhi(2)
      do i = domlo(1), domhi(1)
        call random_number(tmp)
        tmp = drho * (2.0d0*tmp - 1.0d0)
        fIn(:,i,j) = qw * (1.0d0 + tmp)
        gIn(:,i,j) = qw * (1.0d0 - tmp)
      enddo
      enddo
      ! fill BC
      call fill_pbc(fIn, nx,ny)
      call fill_pbc(gIn, nx,ny)
      ! derive macro variables
      call derive_macroscopic(fIn, rho1,jx1,jy1, nx,ny)
      call derive_macroscopic(gIn, rho2,jx2,jy2, nx,ny)
      
      ! initial output
      time = 0.0d0
      istep = 0
      iplot = 0
      call output_csv(iplot,time,&
      rho1,rho1,rho1,rho1, &
      rho1,rho1, &
      xh,yh,xh,yh, &
      nx,ny)
      iplot = iplot + 1
            
      ! main loop
      do istep = 1, nstep
        print *, 'step=',istep, 'time=',time
        
        ! macroscopic
        !derive_macroscopic(fIn, rho1,jx1,jy1, nx,ny)
        !derive_macroscopic(gIn, rho2,jx2,jy2, nx,ny)
        
        !
        do j = domlo(2), domhi(2)
        do i = domlo(1), domhi(1)
          !
          rhotot_omega = rho1(i,j)*omega1 + rho2(i,j)*omega2
          utotx = (jx1(i,j)*omega1 + jx2(i,j)*omega2) / rhotot_omega
          utoty = (jy1(i,j)*omega1 + jy2(i,j)*omega2) / rhotot_omega
          
          !
          rhoContrib1x = 0.0d0
          rhoContrib1y = 0.0d0
          rhoContrib2x = 0.0d0
          rhoContrib2y = 0.0d0
          do l = 2, nq
            iup = i - qex(l)
            jup = j - qey(l)
            rhoContrib1x = rhoContrib1x + rho1(iup,jup)*qw(l)*qex(l)
            rhoContrib1y = rhoContrib1y + rho1(iup,jup)*qw(l)*qey(l)
            rhoContrib2x = rhoContrib2x + rho2(iup,jup)*qw(l)*qex(l)
            rhoContrib2y = rhoContrib2y + rho2(iup,jup)*qw(l)*qey(l)
          enddo
          
          ! potential contribution
          utotx1 = utotx - G/omega1*rhoContrib2x
          utoty1 = utoty - G/omega1*rhoContrib2y
          !
          utotx2 = utotx - G/omega2*rhoContrib1x
          utoty2 = utoty - G/omega2*rhoContrib1y
          
          ! collision
          do l = 1,nq
            cu1 = 3.0d0 * (qex(l)*utotx1 + qey(l)*utoty1)
            cu2 = 3.0d0 * (qex(l)*utotx2 + qey(l)*utoty2)
            
            fEq = rho1(i,j) * qw(l) * &
            (1.0d0 + cu1 + 0.5d0*cu1*cu1 - 1.5d0*(utotx1**2+utoty1**2)) 
            gEq = rho2(i,j) * qw(l) * &
            (1.0d0 + cu2 + 0.5d0*cu2*cu2 - 1.5d0*(utotx2**2+utoty2**2))
            
            fOut(l,i,j) = fIn(l,i,j) - omega1 * (fIn(l,i,j)-fEq)
            gOut(l,i,j) = gIn(l,i,j) - omega2 * (gIn(l,i,j)-gEq)
          enddo
        enddo
        enddo
        
        !
        call fill_pbc(fOut, nx,ny)
        call fill_pbc(gOut, nx,ny)
        do j = domlo(2), domhi(2)
        do i = domlo(1), domhi(1)
          do l = 1,nq
            iup = i - qex(l)
            jup = j - qey(l)
            
            fIn(l,i,j) = fOut(l,iup,jup)
            gIn(l,i,j) = gOut(l,iup,jup)
          enddo
        enddo
        enddo
        
        !
        call fill_pbc(fIn, nx,ny)
        call fill_pbc(gIn, nx,ny)
        call derive_macroscopic(fIn, rho1,jx1,jy1, nx,ny)
        call derive_macroscopic(gIn, rho2,jx2,jy2, nx,ny)
        
        time = time + dt
        
        if (mod(istep,plot_int)==0 .or. istep==nstep) then
          call output_csv(iplot,time,&
      rho1,rho1,rho1,rho1, &
      rho1,rho1, &
      xh,yh,xh,yh, &
      nx,ny)
      iplot = iplot + 1
        endif
      enddo ! main loop
      
      
      end program main2d
