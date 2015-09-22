
      module const_module
      implicit none
      !
      double precision, parameter :: PI = 4.0d0 * atan(1.0d0)
      double precision, parameter :: zero = 0.0d0
      double precision, parameter :: half = 0.5d0
      double precision, parameter :: one = 1.0d0, two = 2.0d0
      double precision, parameter :: three = 3.0d0, four = 4.0d0
      !double precision, parameter :: five = 5.0d0, six = 6.0d0
      !double precision, parameter :: seven = 7.0d0, eight = 8.0d0
      double precision, parameter :: nine = 9.0d0
      
      
      !
      integer, parameter :: COORDSYS_CART = 0
      !integer, parameter :: COORDSYS_RZ = 1
      !
      integer, parameter :: BCTYPE_PER = 0
      !integer, parameter :: BCTYPE_SYM = 1
      
      end module const_module
      
      
      
      
      module lbm2d_module
      implicit none
      ! D2Q9
      integer, parameter :: nd = 2, nq = 9
      double precision, parameter :: qw(9) = (/ &
      4.0d0/9.0d0, &
      1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,&
      1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0 /)
      integer, parameter :: qex(9) = (/ 0, 1, 0, -1, 0, 1, -1, -1, 1 /)
      integer, parameter :: qey(9) = (/ 0, 0, 1, 0, -1, 1, 1, -1, -1 /)
      integer, parameter :: qopp(9) = (/ 1, 4, 5, 2, 3, 8, 9, 6, 7 /)
      
      
      
      end module lbm2d_module
      
      module lbm3d_module
      implicit none
      ! D3Q19
      integer, parameter :: nd = 3, nq = 19
      double precision, parameter :: qw(nq) = (/ 1.0d0/3.0d0, &
      1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, &
      1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, &
      1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)
      integer, parameter :: qex(nq) = (/ 0, -1, 1,  0, 0,  0, 0, &
      -1,  1, -1, 1, -1,  1, -1, 1,  0,  0,  0, 0 /)
      integer, parameter :: qey(nq) = (/ 0,  0, 0, -1, 1,  0, 0, &
      -1, -1,  1, 1,  0,  0,  0, 0, -1,  1, -1, 1 /)
      integer, parameter :: qez(nq) = (/ 0,  0, 0,  0, 0, -1, 1, &
       0,  0,  0, 0, -1, -1,  1, 1, -1, -1,  1, 1 /)
      
      integer, parameter :: qe(nd,nq) = reshape((/ &
       0,  0,  0, &
      -1,  0,  0, &
       1,  0,  0, &
       0, -1,  0, &
       0,  1,  0, &
       0,  0, -1, &
       0,  0,  1, &
      -1, -1,  0, &
       1, -1,  0, &
      -1,  1,  0, &
       1,  1,  0, &
      -1,  0, -1, &
       1,  0, -1, &
      -1,  0,  1, &
       1,  0,  1, &
       0, -1, -1, &
       0,  1, -1, &
       0, -1,  1, &
       0,  1,  1  /), shape(qe))
      
      end module lbm3d_module
