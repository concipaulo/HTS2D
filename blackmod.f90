module solver
implicit none
contains
!
!This module contains all the subroutines used to solve the linear system. This
!module, is portable, and work by itself, for example, you wrote a FEM program
!you could use it to solve the linear equation arrived from that program. The
!subroutines will be explained so you can call it and tweak it as you pleased.
subroutine bound_cond(msh, nodes)
!
! This subroutine sets the values of temperature on boundaries.
implicit none
    integer(kind=4), intent(in) :: nodes
    real(kind=8),dimension(nodes, 4), intent(inout) :: msh
!
!
    integer(kind=4) :: i
    real(kind=8) :: t, mxval, myval, crdx
    real(kind=8), dimension(nodes) :: temp
!
    mxval = maxval(msh(1:nodes,2))
    myval = maxval(msh(1:nodes,3))
!
! creating a dummy array for workspace
    temp(1:nodes) = 0.0d0 ! initializing
!
    do i=1,nodes
!x = xmax
        if(msh(i,2) .eq. mxval)then
        temp(i) = 0.0d0
! y = ymax and x =/ xmax and x=/ 0
        elseif(msh(i,3) .eq. myval .and. msh(i,2) .ne. mxval &
        .and. msh(i,2) .ne. 0.0)then
        temp(i) = 0.0d0
! care with steady state is needed
        crdx = msh(i,2)
            call temp_dist(crdx, t)
        temp(i) = t
!y = ymax
        elseif(msh(i,3) .eq. myval)then
       temp(i) = 0.0d0
! x = 0
        elseif(msh(i,2) .eq. 0.0d0)then
        temp(i) = 0.0d0
! y = 0
        elseif(msh(i,3) .eq. 0.0d0)then
        temp(i) = 0.0d0
!        crdx = msh(i,2)
!            call temp_dist(crdx, t)
!        temp(i) = t
        else
        temp(i) = 0.0d0
        endif
    enddo
!
!copying into the returning array
    msh(1:nodes,4) = temp(1:nodes)
return
end subroutine bound_cond
!********************************************************************************
subroutine temp_dist(crdx, temp)
!
! This subroutine is capable of set a temperature on certain nodes. You need to 
!give the coordinate and it'll return the temperature on that node.
implicit none
!
    real(kind=8), intent(in) :: crdx
    real(kind=8), intent(out) :: temp
    real(kind=8), parameter :: pi = 4.0d0*datan(1.0_8)
!
    ! temp = dsin(pi*crdx)
    temp = 200.d0*dsin((pi*crdx)/30.0d0)
return
end subroutine temp_dist
!********************************************************************************
subroutine temp_exact(msh, nodes)
!
! On steady or simple problems were is known the function that represents the 
!temperature distribution, you give and array of nodes and it'll set the exact
!solution or analytical.
implicit none
!
    integer(kind=4), intent(in) :: nodes
    real(kind=8), dimension(nodes,4), intent(inout) :: msh
!
!local variables
    integer(kind=4) :: i
    real(kind=8), dimension(nodes) :: temp
    real(kind=8) :: tm = 200.0d0
    real(kind=8) :: L = 30.0d0
    real(kind=8) :: b = 20.0d0
    real(kind=8), parameter :: pi = 4.0d0*datan(1.0_8)
!
    do i = 1,nodes
        temp(i) = tm*(dsinh(pi*msh(i,3)/L)/dsinh(b*pi/L))*&
                    dsin(pi*msh(i,2)/L)
    end do
!
    msh(1:nodes,4) = temp(1:nodes)
!
return
end subroutine temp_exact
!********************************************************************************
subroutine exactreddy(crdx, crdy, t)
  use physics
  implicit none
!
!External variables
  real(kind=8), intent(in) :: crdx
  real(kind=8), intent(in) :: crdy
!
! Internal variables
  real(kind=8), parameter :: pi = 4.0d0*datan(1.0_8)
  integer(kind=4) :: n
  real(kind=8) :: an, s, summ, t
!
  summ = 0.d0
  do n = 1, 50
    an = 0.5d0*(2.d0*dfloat(n) - 1.d0)*pi
    s = ((-1**n)*cos(an*crdy)*cosh(an*crdx))/((an**3)*cosh(an))
    summ = summ + s
    ! write(*,*) "sum is =", summ
  enddo
!
  t = (e/2.d0*k)*((1-crdy**2) + 4.d0*summ)
return
stop
end subroutine exactreddy
!********************************************************************************
subroutine reddy(msh, n)
  implicit none
!
!External variables
integer(kind=4), intent(in) :: n
real(kind=8), intent(inout), dimension(n,4) :: msh
!
!Internal variables
integer(kind=4) :: i
real(kind=8) :: t
!
  do i = 1, n
    call exactreddy(msh(i,2), msh(i,3), t)
    msh(i,4) = t
  enddo
!
return
stop
end subroutine reddy
!********************************************************************************
subroutine logerror(x, x_exact, nodes, dxdy, dirname)
!
! This subroutine is responsable to mensurate the error between the analytical 
!and the approximate solution.
use parmesh
implicit none
!
    character(*) :: dirname
    integer(kind=4), intent(in) :: nodes
    real(kind=8) :: dxdy
    real(kind=8), dimension(nodes), intent(in) :: x
    real(kind=8), dimension(nodes), intent(in) :: x_exact
!
!local variables
    real(kind=8), dimension(nodes) :: dummy
    integer(kind=4) :: i
    real(kind=8) :: loge2
    real(kind=8) :: loge1
    real(kind=8) :: loge_inf
    real(kind=8) :: summ, l_two, l_one, l_inf, logh
!
    do i=1,nodes
        dummy(i) = x(i) - x_exact(i)
    enddo
!
    summ = 0.0d0
    do i=1, nodes
        summ = summ + dummy(i)*dummy(i)
    enddo
!
    l_two = dsqrt(summ/dfloat(nodes))
!
    l_one = abs(sum(dummy(1:nodes))/dfloat(nodes))
    l_inf = maxval(dummy(1:nodes))
!
    loge1 = dlog10(l_one)
    loge2 = dlog10(l_two)
    loge_inf = dlog10(l_inf)
    logh = abs(dlog10(dxdy))
!
    open(unit=12, file= trim(adjustl(trim(dirname)//'logfile.dat')), status='unknown')
    write(12,*) 'L1 norm'
    write(12,10) logh, loge1
    write(12,10) dxdy, l_one
    write(12,*) 'L2 norm'
    write(12,10) logh, loge2
    write(12,10) dxdy, l_two
    write(12,*) 'L inf'
    write(12,10) logh, loge_inf
    write(12,10) dxdy, l_inf
10 format(2(F30.16))
!
   close(unit=12)
!
return
end subroutine logerror
!********************************************************************************
subroutine gs(a, rhs, nodes, x_exact, fout, dirname)
use parsolvers
!
!The main is to implement the gauss siedel interation method to solve
! linear systems Ax = b
!       Programmer: Paulo Conci                 Date: 20/06/2016
!                                        Last update: 21/06/2016
implicit none
!
    character(*) :: dirname
    integer*4, intent(in) :: nodes
    real*8, dimension(nodes,nodes), intent(in) :: a
    real*8, dimension(nodes), intent(in) :: rhs
    real*8, dimension(nodes) :: x
    real*8, dimension(nodes) :: x_init
    real*8, dimension(nodes), intent(inout) :: x_exact
!
    integer*4 :: i, j, cont
    real*8, dimension(nodes) :: dummy
    real*8 :: x_error, summ
!
   character(len=*) :: fout
!
        open(unit=11, file=dirname//fout, status='unknown')
!
        x(1:nodes) = x_exact(1:nodes)
!
        cont = 0
!
10 continue
        x_init(1:nodes) = x(1:nodes)
!
!
        do i= 1, nodes
                call mat_mul(a, x, nodes, i, summ)
                x(i) = (1.0d0 / a(i,i))*(rhs(i) - summ)
        end do
!
!
        do j = 1, nodes
        dummy(j) = abs(x(j) - x_init(j))
        end do
!
        x_error = maxval(dummy(1:nodes))
!
        if(x_error .le. tol_abs)then
        go to 15
        else if(x_error .gt. tol_abs)then
        cont = cont + 1
        write(11,*)  cont, x_error
        go to 10
        end if
!
!
15 continue
        write(11,*)  cont, x_error
        write(*,*) 'Final error of this iteration is =', x_error
!
        x_exact = x
!
return
end subroutine gs
subroutine mat_mul(a, x, nodes, i, summ)
!
! Subroutine that multiplies the a and x
implicit none
!
    integer*4 :: nodes, i, j
    real*8, dimension(nodes,nodes), intent(in) :: a
    real*8, dimension(nodes), intent(in) :: x
    real*8, intent(out):: summ
!
!
        summ = 0.0d0
        do j = 1, nodes
                if (j .ne. i)then
                summ = summ + (a(i,j) * x(j))
                end if
        end do
!
!
return
end subroutine mat_mul
!********************************************************************************
subroutine cvgc(x, x0, t, nodes, normt)
!
!This subroutine will test if the heat transfer alread reached the
!stead state
implicit none
!
    integer*4 :: nodes
    real*8, dimension(nodes), intent(in) :: x
    real*8, dimension(nodes), intent(in) :: x0
    real(kind=8), dimension(nodes), intent(in) :: t
    real*8, intent(out) :: normt
!
!
    real*8 :: nom, denom, tnow, told
!    real(kind=8), dimension(nodes) :: dummy
!
!
    tnow = dsqrt(dot_product(x,x))
!    write(*,*) 'tnow = ', tnow
    told = dsqrt(dot_product(t,t))
!    write(*,*) 'told = ', told
    denom = sqrt(dot_product(x0,x0))
!
    nom = abs(tnow - told)
!
    normt = nom/tnow
!
    write(*,*) 'Temperature Norm = ', normt
!
return
end subroutine cvgc
!********************************************************************************
subroutine temp_time_plot(msh, nodes, dxdy)
!
! This is a subrotine for ploting.
use parmesh
implicit none
    integer(kind=4), intent(in) :: nodes
    real(kind=8) :: dxdy
    real(kind=8), dimension(nodes,4) :: msh
!
! internal variables
    integer(kind=4) :: i, r
    integer(kind=4) :: ioerr
    real(kind=8) :: xmax, ymax, eps, crdx
    integer(kind=4) :: row, col
!
        eps = dxdy/2.0d0
!        msh(1:nodes, 4) = t(1:nodes)
        xmax = maxval(msh(1:nodes,2))
        ymax = maxval(msh(1:nodes,3))
!
!        write(*,*)'xm', xmax, 'ym', ymax
        row = (int(ymax/dxdy)) + 1
        col = (int(xmax/dxdy)) + 1
!
    open(unit = 16, file = 'graph1.plt', status ='unknown', iostat = ioerr)

    write(16,*) 'VARIABLES = "X", "T"'
    write(16,*) 'ZONE    I =', col, '    F=POINT'
            do i = 1, col
                crdx = (dfloat(i-1))*dxdy
                    do r = 1, nodes
                        if ((abs(msh(r,2) - crdx)) .lt. eps .and. (msh(r,3) .eq. 0.0)) then
                            write(16,320) msh(r,2), abs(msh(r,4))
                        end if
                    end do
            end do
320 format (3(F30.16))
!    do i = 1, nodes
!        if(msh(i,3) .eq. 0.d0)then
!            write(16, 100) msh(i,2), abs(msh(i,4))
!        endif
!    enddo
!100 format(2(F30.16))
return
!
end subroutine temp_time_plot
!********************************************************************************
subroutine temp_origin(msh, nodes, delta, dxdy)
!
! Ploting values of temperature.
use physics
use parmesh
implicit none
    integer(kind=4), intent(in) :: nodes
    real(kind=8) :: dxdy
    real(kind=8), dimension(nodes,4) :: msh
!
! internal variables
    integer(kind=4) :: i
    integer(kind=4) :: ioerr, delta
    real(kind=8) :: xmax, ymax, time
!
    open(unit = 17, file = 'graph2.plt', status ='unknown', iostat = ioerr)
!
    time = (float(delta))*deltat
    do i = 1, nodes
        if((msh(i,3) .eq. 0.d0) .and. (msh(i,2) .eq. 0.d0))then
            write(17, 100)  time, abs(msh(i,4))
        endif
    enddo

!
100 format(2(F30.16))
!
return
end subroutine temp_origin
!********************************************************************************
! subroutine connect(msh, conn, nodes, dxdy)
! implicit none
! !
! !External variables
! integer(kind=4) :: nodes
! real(kind=8), intent(in) :: dxdy
! real(kind=8), intent(in), dimension(nodes,4) :: msh
! integer(kind=4), intent(inout), dimension(nodes,5) :: conn
! !
! ! Internal variables
! real(kind=8) :: crdx, crdy, plus_x, plus_y, minus_x, minus_y, eps, xmax, ymax
! integer(kind=4) :: i,j
! !
!     eps = dxdy/2
!     xmax = maxval(msh(1:nodes,2))
!     ymax = maxval(msh(1:nodes,3))
! !
!     do i=1,nodes
!         conn(i,1) = msh(i,1)
!         crdx = msh(i,2)
!         crdy = msh(i,3)
!         plus_x = msh(i,2) + dxdy
!         minus_x = msh(i,2) - dxdy
!         plus_y = msh(i,3) + dxdy
!         minus_y = msh(i,3) - dxdy
! !
!         do j = 1, nodes
!             if (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
!             ((abs(msh(j,3) - crdy)) .lt. eps)) then
!                 conn(i,2) = nint(msh(j,1))
!             elseif(((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
!             ((abs(msh(j,3) - crdy)) .lt. eps))then
!                 conn(i,4) = nint(msh(j,1))
!             elseif (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!             ((abs(msh(j,3) - plus_y)) .lt. eps)) then
!                 conn(i,3) = nint(msh(j,1))
!             elseif (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!             ((abs(msh(j,3) - minus_y)) .lt. eps)) then
!                  conn(i,5) = nint(msh(j,1))
!             elseif((abs(msh(j,3) - crdy) .lt. eps) .and. (plus_x .gt. xmax))then
!                 conn(i,2) = 0
!             elseif ((abs(msh(j,2) - crdx) .lt. eps) .and. (plus_y .gt. ymax))then
!                 conn(i,3) = 0
!             elseif((abs(msh(j,3) - crdy) .lt. eps) .and. minus_x .lt. 0)then
!                 conn(i,4) = 0
!             elseif((abs(msh(j,2) - crdx) .lt. eps) .and. minus_y .lt. 0)then
!                 conn(i,5) = 0
!             end if
!         enddo
!     enddo
! return
! end subroutine connect
!********************************************************************************
subroutine connect(msh, conn, nodes, dxdy)
implicit none
!
!External variables
integer(kind=4) :: nodes
real(kind=8), intent(in) :: dxdy
real(kind=8), intent(in), dimension(nodes,4) :: msh
integer(kind=4), intent(inout), dimension(nodes,9) :: conn
!
! Internal variables
real(kind=8) :: crdx, crdy, plus_x, plus_y, minus_x, minus_y, eps, xmax, ymax
integer(kind=4) :: i,j
!
    eps = dxdy/2
    xmax = maxval(msh(1:nodes,2))
    ymax = maxval(msh(1:nodes,3))
!
  conn = 0
    do i=1,nodes
    ! write(*,*) 'dxdy=', dxdy
      conn(i,1) = msh(i,1)
      crdx = msh(i,2)
      crdy = msh(i,3)
      plus_x = msh(i,2) + dxdy
      minus_x = msh(i,2) - dxdy
      plus_y = msh(i,3) + dxdy
      minus_y = msh(i,3) - dxdy
!
! write(*,*) 'crdx=', crdx, 'crdy=', crdy, 'plusx=', plus_x, 'plusy', plus_y
      do j = 1, nodes
        if (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
        ((abs(msh(j,3) - crdy)) .lt. eps)) then
            conn(i,2) = nint(msh(j,1))
        elseif(((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
        ((abs(msh(j,3) - crdy)) .lt. eps))then
            conn(i,4) = nint(msh(j,1))
        elseif (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
        ((abs(msh(j,3) - plus_y)) .lt. eps)) then
            conn(i,3) = nint(msh(j,1))
        elseif (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
        ((abs(msh(j,3) - minus_y)) .lt. eps)) then
             conn(i,5) = nint(msh(j,1))
        elseif (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
        ((abs(msh(j,3) - minus_y)) .lt. eps)) then
             conn(i,6) = nint(msh(j,1))
        elseif (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
        ((abs(msh(j,3) - plus_y)) .lt. eps)) then
             conn(i,7) = nint(msh(j,1))
        elseif (((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
        ((abs(msh(j,3) - plus_y)) .lt. eps)) then
             conn(i,8) = nint(msh(j,1))
        elseif (((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
        ((abs(msh(j,3) - minus_y)) .lt. eps)) then
             conn(i,9) = nint(msh(j,1))
        elseif((abs(msh(j,3) - crdy) .lt. eps) .and. (plus_x .gt. xmax))then
            conn(i,2) = 0
            conn(i,6) = 0
            conn(i,7) = 0
        elseif ((abs(msh(j,2) - crdx) .lt. eps) .and. (plus_y .gt. ymax))then
            conn(i,3) = 0
            conn(i,7) = 0
            conn(i,8) = 0
        elseif((abs(msh(j,3) - crdy) .lt. eps) .and. minus_x .lt. 0)then
            conn(i,4) = 0
            conn(i,8) = 0
            conn(i,9) = 0
        elseif((abs(msh(j,2) - crdx) .lt. eps) .and. minus_y .lt. 0)then
            conn(i,5) = 0
            conn(i,6) = 0
            conn(i,9) = 0
        end if
      enddo
  enddo
return
end subroutine connect
!********************************************************************************
subroutine ccall(c, n, a, rhs, t)
use constants
implicit none
!
! External Variables
integer(kind=4), intent(in) :: n
integer(kind=4), intent(in), dimension(n,9) :: c
real(kind=8), intent(inout), dimension(n,n) :: a
real(kind=8), intent(inout), dimension(n) :: rhs
real(kind=8), intent(inout), dimension(n) :: t
!
!Internal variables
integer(kind=4) :: i, ii
integer(kind=4), dimension(9) :: v
real(kind=8) :: n1, n2, n3, n4, n5, n6, n7, n8, n9
!
  do i=1, n
    if(c(i,2) .eq. 0) then
      if(c(i,3) .eq. 0) then 
        n1 = xmax_ymax_n
        n2 = 0
        n3 = 0
        n4 = xmax_ymax_tl
        n5 = xmax_ymax_td
        n6 = 0
        n7 = 0
        n8 = 0
        n9 = xmax_ymax_tld
        ! write(*,*) "xmax_ymax"
        v(1:9) = c(i,1:9)
        ! write(*,*) (v(ii),ii=1,5)
        call create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
      elseif(c(i,5) .eq. 0) then
        ! write(*,*) "oy_xmax"
        n1 = oy_xmax_n
        n2 = 0
        n3 = oy_xmax_tt
        n4 = oy_xmax_tl
        n5 = 0
        n6 = 0
        n7 = 0
        n8 = oy_xmax_tlt
        n9 = 0
        v(1:9) = c(i,1:9)
        ! write(*,*) (v(ii),ii=1,5)
        call create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
      else
        n1 = wall_one_n
        n2 = 0
        n3 = wall_one_tt
        n4 = wall_one_tl
        n5 = wall_one_td
        n6 = 0
        n7 = wall_one_trt
        n8 = wall_one_tlt
        n9 = 0
        ! write(*,*) "wall_1"
        v(1:9) = c(i,1:9)
        ! write(*,*) (v(ii),ii=1,5)
        call create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
      endif
    elseif(c(i,3) .eq. 0) then
      if(c(i,4) .eq. 0) then
        n1 = ox_ymax_n
        n2 = ox_ymax_tr
        n3 = 0
        n4 = 0
        n5 = ox_ymax_td
        n6 = ox_ymax_trd
        n7 = 0
        n8 = 0
        n9 = 0
        ! write(*,*) "ox_ymax"
        v(1:9) = c(i,1:9)
        ! write(*,*) (v(ii),ii=1,5)
        call create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
      elseif(c(i,2) .eq. 0) then
        cycle
      else
        n1 = wall_two_n
        n2 = wall_two_tr
        n3 = 0
        n4 = wall_two_tl
        n5 = wall_two_td
        n6 = 0
        n7 = 0
        n8 = wall_two_tlt
        n9 = wall_two_tld
        ! write(*,*) "wall_2"
        v(1:9) = c(i,1:9)
        ! write(*,*) (v(ii),ii=1,5)
        call create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
      endif
    elseif(c(i,4) .eq. 0) then
      if(c(i,5) .eq. 0) then
        n1 = oxoy_n
        n2 = oxoy_tr
        n3 = oxoy_tt
        n4 = 0
        n5 = 0
        n6 = 0
        n7 = oxoy_trt
        n8 = 0
        n9 = 0
        ! write(*,*) "ox_oy"
        v(1:9) = c(i,1:9)
        ! write(*,*) (v(ii),ii=1,5)
        call create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
      elseif(c(i,3) .eq. 0) then
        cycle
      else
        n1 = wall_three_n
        n2 = wall_three_tr
        n3 = wall_three_tt
        n4 = 0
        n5 = wall_three_td
        n6 = wall_three_trd
        n7 = 0
        n8 = 0
        n9 = wall_three_tld
        ! write(*,*) "wall_3"
        v(1:9) = c(i,1:9)
        ! write(*,*) (v(ii),ii=1,5)
        call create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
      endif
    elseif(c(i,5) .eq. 0)then
      if(c(i,4) .eq. 0) then
        cycle
      elseif(c(i,2) .eq. 0) then
        cycle
      else
        n1 = wall_four_n
        n2 = wall_four_tr
        n3 = wall_four_tt
        n4 = wall_four_tl
        n5 = 0
        n6 = wall_four_trd
        n7 = wall_four_trt
        n8 = 0
        n9 = 0
        ! write(*,*) "wall_4"
        v(1:9) = c(i,1:9)
        ! write(*,*) (v(ii),ii=1,5)
        call create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
      endif
    else
        n1 = interior_n
        n2 = interior_tr
        n3 = interior_tt
        n4 = interior_tl
        n5 = interior_td
        n6 = interior_trd
        n7 = interior_trt
        n8 = interior_tlt
        n9 = interior_tld
        ! write(*,*) "interior"
        v(1:9) = c(i,1:9)
        ! write(*,*) (v(ii),ii=1,5)
        call create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
    endif
  enddo
!
return
end subroutine ccall
!********************************************************************************
subroutine create (v, a, n, rhs, t, n1, n2, n3, n4, n5, n6, n7, n8, n9)
implicit none
!
! External Variable
integer(kind=4), intent(in) :: n
integer(kind=4), intent(in), dimension(9) :: v
real(kind=8), intent(out), dimension(n,n) :: a
real(kind=8), intent(out), dimension(n) :: rhs
real(kind=8), intent(out), dimension(n) :: t
real(kind=8), intent(in) :: n1, n2, n3, n4, n5, n6, n7, n8, n9
!
! Internal variables
!
  a(v(1),v(1)) = n1
  rhs(v(1)) = t(v(1))
!
  if(v(2) .ne. 0) then
    a(v(1),v(2)) = n2
  ! elseif(v(2) .eq. 0) then
  !   a(v(1),v(2)) = 0.d0
  endif
!
  if(v(3) .ne. 0)then
    a(v(1),v(3)) = n3
  ! elseif(v(3) .eq. 0) then
  !   a(v(1),v(3)) = 0
  endif
!
  if(v(4) .ne. 0) then
    a(v(1),v(4)) = n4
  ! elseif(v(4) .eq. 0) then
  !   a(v(1),v(4)) = 0
  endif 
!
  if(v(5) .ne. 0) then
    a(v(1),v(5)) = n5
  ! elseif(v(5) .eq. 0) then
  !   a(v(1),v(5)) = 0
  endif
!
  if(v(6) .ne. 0) then
    a(v(1),v(6)) = n6
  ! elseif(v(6) .eq. 0) then
  !   a(v(1),v(6)) = 0
  endif
!
  if(v(7) .ne. 0) then
    a(v(1),v(7)) = n7
  ! elseif(v(7) .eq. 0) then
  !   a(v(1),v(7)) = 0
  endif
!
  if(v(8) .ne. 0) then
    a(v(1),v(8)) = n8
  ! elseif(v(8) .eq. 0) then
  !   a(v(1),v(8)) = 0
  endif
!
  if(v(9) .ne. 0) then
    a(v(1),v(9)) = n9
  ! elseif(v(9) .eq. 0) then
  !   a(v(1),v(9)) = 0
  endif
!
return
end subroutine create
!********************************************************************************
subroutine create_a_steady(msh, t, nodes, alin, rhs, dxdy)
use parmesh
!
!This is a subroutine to generate the matrix of equations and the rhs
!for the numerical solver, this require the mesh array!
!
!       Programmer: Paulo Conci         Last Modify: August 29
!
implicit none
! Variables o main prog
    integer*4, intent(in) :: nodes
    real(kind=8) :: dxdy
    real*8, dimension(nodes,4), intent(in) :: msh
    real*8, dimension(nodes), intent(in) :: t
    real*8, dimension(nodes,nodes), intent(out) :: alin
    real*8, dimension(nodes), intent(out) :: rhs
!
!
!Variables of subroutine
    integer*4 :: i, j, c, d, kt
    integer*4 :: td, tl, tr, tt, tempxy
    real*8 :: ymax, xmax, crdx, crdy, eps, temp
    real*8 :: plus_x, minus_x, plus_y, minus_y
    real(kind=8) :: cdx, cdy, temprhs
!
!
!finding max values of matrix of mesh
        xmax = maxval(msh(1:nodes,2))
        ymax = maxval(msh(1:nodes,3))
!
!
        eps = dxdy/2.0
!
        alin = 0.0d0
!
!
!
 do i = 1, nodes
!=====================================1========================================
! x = xmax and y=/ 0 and y =/ ymax
    if ((msh(i,2) .eq. xmax) .and. msh(i,3) .ne. 0.0 .and. &
        msh(i,3) .ne. ymax) then
            tempxy = nint(msh(i,1))
            temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    endif
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) =  temp
                endif
        enddo
!
!
!
! isolated nodes implemetation in y
! *******************************************************************
!        minus_x = msh(i,2) - dxdy
!        minus_y = msh(i,3) - dxdy
!        plus_y = msh(i,3) + dxdy
!        crdx = msh(i,2)
!        crdy = msh(i,3)
!        tempxy = nint(msh(i,1))
!
!        tt = 0
!        tl = 0
!        td = 0
!
!    do j = 1, nodes
!        if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!            ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                tt = nint(msh(j,1))
!        else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!            ((abs(msh(j,2) - minus_x)) .lt. eps))then
!                tl = nint(msh(j,1))
!        else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_y)) .lt. eps))then
!                td = nint(msh(j,1))
!        end if
!    end do
!
!
!    do c = 1, nodes
!        do d = 1, nodes
!            if((c .eq. tempxy) .and. (d .eq. tt)) then
!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tl))then
!               alin(c,d) = -1.0d0/(1.0d0*l**2)
!            else if ((c .eq. tempxy) .and. (d .eq. td)) then
!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 2.0d0/l**2
!            end if
!        enddo
!    enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) = 0.0d0
!                end if
!        end do
!=====================================2========================================
!x = 0  and y =/0 and y =/ ymax
    else if ((msh(i,2) .eq. 0.0) .and. (msh(i,3) .ne. ymax) &
        .and. msh(i,3) .ne. 0.0) then
!
!        plus_x = msh(i,2) + dxdy
!        minus_y = msh(i,3) - dxdy
!        plus_y = msh(i,3) + dxdy
!        crdx = msh(i,2)
!        crdy = msh(i,3)
!        tempxy = nint(msh(i,1))
!        temp = t(i)
!Implementation that do not change the temperature
!********************************************************************
        tempxy = nint(msh(i,1))
        temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                if((c .eq. tempxy).and.(d .eq. tempxy))then
                    alin(c,d) = 1.0d0
                end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!********************************************************************
!        tt = 0
!        tr = 0
!        td = 0
!!
!    do j = 1, nodes
!        if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!            ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                tt = nint(msh(j,1))
!!
!        else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!            ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                tr = nint(msh(j,1))
!!
!        else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_y)) .lt. eps))then
!                td = nint(msh(j,1))
!!
!        end if
!    end do
!!
!!
!    do c = 1, nodes
!        do d = 1, nodes
!            if((c .eq. tempxy) .and. (d .eq. tt))then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. td)) then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 2.0d0/dxdy**2
!            end if
!        enddo
!    enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                    rhs(kt) =  - temp + 1.0
!                end if
!        end do
!
!
!=====================================3========================================
! y = ymax and x=/0 and x=/ max
    else if ((msh(i,3) .eq. ymax) .and. (msh(i,2) .ne. xmax) .and. &
            msh(i,2) .ne. 0.0) then
! na parte superior
!
!
    tempxy = nint (msh(i,1))
    temp = t(i)
        do c = 1, nodes
            do d = 1, nodes
                if((c .eq. tempxy).and.(d .eq. tempxy))then
                    alin(c,d) = 1.0d0
                end if
            enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!
!
!                plus_x = msh(i,2) + dxdy
!                minus_x = msh(i,2) - dxdy
!                minus_y = msh(i,3) - dxdy
!                crdx = msh(i,2)
!                crdy = msh(i,3)
!                tempxy = nint(msh(i,1))
!                temp = t(i)
!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy).and.(d .eq. tempxy))then
!                alin(c,d) = 1.0d0
!                end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                rhs(j) = temp
!                end if
!        enddo
!
!        write(*,*) 'no',tempxy
!
!                tt = 0
!                tr = 0
!                tl = 0

!        do j = 1, nodes
!                if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_y)) .lt. eps))then
!                td = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                tr = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - minus_x)) .lt. eps))then
!                tl = nint(msh(j,1))
!                end if
!        end do



!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy) .and. (d .eq. td)) then
!                alin(c,d) = 2.0d0*fo
!                else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                alin(c,d) = fo
!                else if ((c .eq. tempxy) .and. (d .eq. tl)) then
!                alin(c,d) = fo
!                else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = a2
!                end if
!                enddo
!        enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) = b2 - t(kt)
!                end if
!        end do

!=====================================4========================================
! x =/ 0 and x =/ xmax and y = 0
else if ((msh(i,2) .ne. 0.0) .and. (msh(i,2) .ne. xmax) .and.&
        (msh(i,3) .eq. 0.0)) then ! nós isolados na parte inferior
!
!                plus_x = msh(i,2) + dxdy
!                minus_x = msh(i,2) - dxdy
!                plus_y = msh(i,3) + dxdy
!                crdx = msh(i,2)
!                crdy = msh(i,3)
                tempxy = nint(msh(i,1))
                temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                if((c .eq. tempxy).and.(d .eq. tempxy))then
                    alin(c,d) = 1.0d0
                end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo

!        write(*,*) 'no',tempxy
!
!                tt = 0
!                tr = 0
!                tl = 0
!!
!        do j = 1, nodes
!                if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                    tt = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                    tr = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - minus_x)) .lt. eps))then
!                    tl = nint(msh(j,1))
!                end if
!        end do
!!
!!
!!
!        do c = 1, nodes
!        do d = 1, nodes
!            if((c .eq. tempxy) .and. (d .eq. tt))then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tl)) then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 2.0d0/dxdy**2
!            end if
!        enddo
!        enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) =  - temp + 1.0
!                end if
!        end do
!=====================================5========================================
! x= 0 and y= 0
elseif ((msh(i,2) .eq. 0.0) .and. msh(i,3) .eq. 0.0)then
        tempxy = nint(msh(i,1))
        temp = t(i)
!        crdx = msh(i,2)
!        crdy = msh(i,3)
!        plus_x = msh(i,2) + dxdy
!        plus_y = msh(i,3) + dxdy
!
!
!        tr = 0
!        tt = 0
!!
!!
!        do j = 1, nodes
!                if (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
!                ((abs(msh(j,3) - crdy)) .lt. eps)) then
!                    tr = nint(msh(j,1))
!                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - plus_y)) .lt. eps)) then
!                    tt = nint(msh(j,1))
!                end if
!        enddo
!!
!!        write(*,*) 'coordenada do nó central', tempxy
!!        write(*,*) 'coordenadas do nós adjacentes', tl, td, tr, tt
!!
!        do c = 1, nodes
!                do d = 1, nodes
!                    if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                        alin(c,d) = 1.0d0/dxdy**2
!                    else if ((c .eq. tempxy) .and. (d .eq. tr)) then
!                        alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!                    else if ((c .eq. tempxy) .and. (d .eq. tt)) then
!                        alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!                    end if
!                enddo
!        enddo
!        do kt = 1, nodes
!                if(kt .eq. tempxy)then
!                    rhs(kt) = - temp + 1.0
!                end if
!        end do
!
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!=====================================6========================================
!x = xmax and y=0
elseif((msh(i,2) .eq. xmax) .and. (msh(i,3) .eq. 0.0))then
    tempxy = nint(msh(i,1))
    temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) =  temp
                end if
        enddo
!=====================================7========================================
! x = xmax and y = ymax
elseif((msh(i,2) .eq. xmax) .and. msh(i,3) .eq. ymax)then
    tempxy = nint(msh(i,1))
    temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!=====================================8========================================
!x = 0 and y = ymax
elseif ((msh(i,2) .eq. 0.0) .and. (msh(i,3) .eq. ymax))then
    tempxy = nint(msh(i,1))
    temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!=====================================i========================================
! interior nodes !!!!
else if ((msh(i,2) .ne. 0.0) .and. (msh(i,2) .ne. xmax) .and.&
        (msh(i,3) .ne. 0.0) .and. (msh(i,3) .ne. ymax))then
                crdx = msh(i,2)
                crdy = msh(i,3)
                plus_x = msh(i,2) + dxdy
                minus_x = msh(i,2) - dxdy
                plus_y = msh(i,3) + dxdy
                minus_y = msh(i,3) - dxdy
                tempxy = nint(msh(i,1))
!
!
        tr = 0
        tl = 0
        tt = 0
        td = 0
!
!
        do j = 1, nodes
                if (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - crdy)) .lt. eps)) then
                    tr = nint(msh(j,1))
                else if(((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - crdy)) .lt. eps))then
                    tl = nint(msh(j,1))
                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
                ((abs(msh(j,3) - plus_y)) .lt. eps)) then
                    tt = nint(msh(j,1))
                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
                ((abs(msh(j,3) - minus_y)) .lt. eps)) then
                    td = nint(msh(j,1))
                end if
        enddo
!
!        write(*,*) 'coordenada do nó central', tempxy
!        write(*,*) 'coordenadas do nós adjacentes', tl, td, tr, tt
!
       do c = 1, nodes
                do d = 1, nodes
                if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
                alin(c,d) = 4.0d0/(1.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. td)) then
                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. tr)) then
                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. tt)) then
                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. tl)) then
                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
                end if
                enddo
        enddo
        do kt = 1, nodes
                if(kt .eq. tempxy)then
                    rhs(kt) = 0.0d0
                end if
        end do
!
end if
!
!
!
!do i = 1, nodes
!!==============================================================================
!! x = xmax and y=/ 0 and y =/ ymax
!if ((msh(i,2) .eq. xmax) .and. msh(i,3) .ne. 0.0 .and. &
!     msh(i,3) .ne. ymax) then
!        tempxy = nint(msh(i,1))
!        temp = t(i)
!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy).and.(d .eq. tempxy))then
!                alin(c,d) = 1.0d0
!                end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                rhs(j) = temp
!                end if
!        enddo
!!
!!
!!
!! isolated nodes implemetation in y
!! *******************************************************************
!!        minus_x = msh(i,2) - dxdy
!!        minus_y = msh(i,3) - dxdy
!!        plus_y = msh(i,3) + dxdy
!!        crdx = msh(i,2)
!!        crdy = msh(i,3)
!!        tempxy = nint(msh(i,1))
!!
!!        tt = 0
!!        tl = 0
!!        td = 0
!!
!!    do j = 1, nodes
!!        if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!!            ((abs(msh(j,3) - plus_y)) .lt. eps))then
!!                tt = nint(msh(j,1))
!!        else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!!            ((abs(msh(j,2) - minus_x)) .lt. eps))then
!!                tl = nint(msh(j,1))
!!        else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!!                ((abs(msh(j,3) - minus_y)) .lt. eps))then
!!                td = nint(msh(j,1))
!!        end if
!!    end do
!!
!!
!!    do c = 1, nodes
!!        do d = 1, nodes
!!            if((c .eq. tempxy) .and. (d .eq. tt)) then
!!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!!            else if ((c .eq. tempxy) .and. (d .eq. tl))then
!!               alin(c,d) = -1.0d0/(1.0d0*l**2)
!!            else if ((c .eq. tempxy) .and. (d .eq. td)) then
!!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!!                alin(c,d) = 2.0d0/l**2
!!            end if
!!        enddo
!!    enddo
!!        do kt = 1, nodes
!!                if (kt .eq. tempxy)then
!!                rhs(kt) = 0.0d0
!!                end if
!!        end do
!!==============================================================================
!!x = 0  and y =/0 and y =/ ymax
!else if ((msh(i,2) .eq. 0.0) .and. (msh(i,3) .ne. ymax) &
!       .and. msh(i,3) .ne. 0.0) then
!!
!        plus_x = msh(i,2) + dxdy
!        minus_x = msh(i,2) - dxdy
!        plus_y = msh(i,3) + dxdy
!        crdx = msh(i,2)
!        crdy = msh(i,3)
!        tempxy = nint(msh(i,1))
!        temp = t(i)
!
!!Implementation that do not change the temperature
!!******************************************************************************
!!        tempxy = nint(msh(i,1))
!!        temp = t(i)
!!        do c = 1, nodes
!!                do d = 1, nodes
!!                if((c .eq. tempxy).and.(d .eq. tempxy))then
!!                alin(c,d) = 1.0d0
!!                end if
!!                enddo
!!        enddo
!!
!!        do j= 1, nodes
!!                if(j .eq. tempxy)then
!!                rhs(j) = temp
!!                end if
!!        enddo
!!********************************************************************
!        tt = 0
!        tr = 0
!        tl = 0

!    do j = 1, nodes
!        if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!            ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                tt = nint(msh(j,1))
!!
!        else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!            ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                tr = nint(msh(j,1))
!!
!        else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_x)) .lt. eps))then
!                tl = nint(msh(j,1))
!!
!        end if
!    end do
!!
!    do c = 1, nodes
!        do d = 1, nodes
!            if((c .eq. tempxy) .and. (d .eq. tt))then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tl)) then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 2.0d0/dxdy**2
!            end if
!        enddo
!    enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) = temp
!                end if
!        end do
!!
!!else if ((msh(i,2) .eq. 0.0) .and. (msh(i,3) .eq. 0.0))then
!!
!!        tempxy = nint(msh(i,1))
!!        plus_x = msh(i,2) + l
!!        plus_y = msh(i,3) + l
!!        crdx = msh(i,2)
!!        crdy = msh(i,3)
!!
!!    tt = 0
!!    tr = 0
!!
!!    do j = 1, nodes
!!        if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!!            ((abs(msh(j,3) - plus_y)) .lt. eps))then
!!                tt = nint(msh(j,1))
!!        else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!!            ((abs(msh(j,2) - plus_x)) .lt. eps))then
!!                tr = nint(msh(j,1))
!!        endif
!!    enddo
!!
!!    do c = 1, nodes
!!        do d = 1, nodes
!!            if((c .eq. tempxy) .and. (d .eq. tt)) then
!!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!!            else if ((c .eq. tempxy) .and. (d .eq. tr))then
!!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!!                alin(c,d) = 1.0d0/(1.0d0*l**2)
!!            end if
!!        enddo
!!    enddo
!!        do kt = 1, nodes
!!                if (kt .eq. tempxy)then
!!                rhs(kt) = 0.0d0
!!                end if
!!        end do
!!
!!==============================================================================
!! y = ymax and x=/0 and x=/ max
!else if ((msh(i,3) .eq. ymax) .and. (msh(i,2) .ne. xmax) .and. &
!        msh(i,2) .ne. 0.0) then
!! na parte superior
!    tempxy = nint (msh(i,1))
!    temp = t(i)
!        do c = 1, nodes
!            do d = 1, nodes
!                if((c .eq. tempxy).and.(d .eq. tempxy))then
!                alin(c,d) = 1.0d0
!                end if
!            enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                rhs(j) = temp
!                end if
!        enddo
!!
!!==============================================================================
!! x =/ 0 and x =/ xmax and y = 0
!else if ((msh(i,2) .ne. 0.0) .and. (msh(i,2) .ne. xmax) .and.&
!        (msh(i,3) .eq. 0.0)) then ! nós isolados na parte inferior
!                plus_x = msh(i,2) + dxdy
!                minus_x = msh(i,2) - dxdy
!                plus_y = msh(i,3) + dxdy
!                crdx = msh(i,2)
!                crdy = msh(i,3)
!                tempxy = nint(msh(i,1))
!                temp = t(i)
!!******************************************************************************
!!        do c = 1, nodes
!!                do d = 1, nodes
!!                if((c .eq. tempxy).and.(d .eq. tempxy))then
!!                alin(c,d) = 1.0d0
!!                end if
!!                enddo
!!        enddo
!!!
!!        do j= 1, nodes
!!                if(j .eq. tempxy)then
!!                rhs(j) = temp
!!                end if
!!        enddo
!!
!!        write(*,*) 'no',tempxy
!!
!                tt = 0
!                tr = 0
!                tl = 0
!!
!        do j = 1, nodes
!                if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                    tt = nint(msh(j,1))
!!
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                    tr = nint(msh(j,1))
!!
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - minus_x)) .lt. eps))then
!                    tl = nint(msh(j,1))
!!
!                end if
!        end do
!!
!!
!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy) .and. (d .eq. tt)) then
!                    alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!                else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                    alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!                else if ((c .eq. tempxy) .and. (d .eq. tl)) then
!                    alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!                else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                    alin(c,d) = 2.0d0/dxdy**2
!                end if
!                enddo
!        enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                    rhs(kt) = 0.0d0
!                end if
!        end do
!!==============================================================================
!! x= 0 and y= 0
!    elseif ((msh(i,2) .eq. 0.0) .and. msh(i,3) .eq. 0.0)then
!        tempxy = nint(msh(i,1))
!        temp = t(i)
!        plus_x = msh(i,2) + dxdy
!        plus_y = msh(i,3) + dxdy
!        crdx = msh(i,2)
!        crdy = msh(i,3)
!!
!        tr = 0
!        tt = 0
!!
!        do j = 1, nodes
!                if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                tt = nint(msh(j,1))
!!
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                tr = nint(msh(j,1))
!!
!                end if
!        end do
!!
!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy) .and. (d .eq. tt)) then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!                else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!                else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 1.0d0/dxdy**2
!                end if
!                enddo
!        enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) = 0.0d0
!                end if
!        end do
!!        do c = 1, nodes
!!                do d = 1, nodes
!!                    if((c .eq. tempxy).and.(d .eq. tempxy))then
!!                        alin(c,d) = 1.0d0
!!                    end if
!!                enddo
!!        enddo
!!!
!!        do j= 1, nodes
!!                if(j .eq. tempxy)then
!!                    rhs(j) = temp
!!                end if
!!        enddo
!!====================================================================
!!x = xmax and y=0
!elseif((msh(i,2) .eq. xmax) .and. (msh(i,3) .eq. 0.0))then
!    tempxy = nint(msh(i,1))
!    temp = t(i)
!        do c = 1, nodes
!                do d = 1, nodes
!                    if((c .eq. tempxy).and.(d .eq. tempxy))then
!                        alin(c,d) = 1.0d0
!                    end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                    rhs(j) = temp
!                end if
!        enddo
!!====================================================================
!! x = xmax and y = ymax
!elseif((msh(i,2) .eq. xmax) .and. msh(i,3) .eq. ymax)then
!    tempxy = nint(msh(i,1))
!    temp = t(i)
!        do c = 1, nodes
!                do d = 1, nodes
!                    if((c .eq. tempxy).and.(d .eq. tempxy))then
!                        alin(c,d) = 1.0d0
!                    end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                    rhs(j) = temp
!                end if
!        enddo
!!====================================================================
!!x = 0 and y = ymax
!elseif ((msh(i,2) .eq. 0.0) .and. (msh(i,3) .eq. ymax))then
!    tempxy = nint(msh(i,1))
!    temp = t(i)
!        do c = 1, nodes
!                do d = 1, nodes
!                    if((c .eq. tempxy).and.(d .eq. tempxy))then
!                        alin(c,d) = 1.0d0
!                    end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                    rhs(j) = temp
!                end if
!        enddo
!====================================================================
!! interior nodes !!!!
!else if ((msh(i,2) .ne. 0.0) .and. (msh(i,2) .ne. xmax) .and.&
!        (msh(i,3) .ne. 0.0) .and. (msh(i,3) .ne. ymax))then
!                crdx = msh(i,2)
!                crdy = msh(i,3)
!                plus_x = msh(i,2) + dxdy
!                minus_x = msh(i,2) - dxdy
!                plus_y = msh(i,3) + dxdy
!                minus_y = msh(i,3) - dxdy
!                tempxy = nint(msh(i,1))
!!
!!
!        tr = 0
!        tl = 0
!        tt = 0
!        td = 0
!!
!!
!        do j = 1, nodes
!                if (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
!                ((abs(msh(j,3) - crdy)) .lt. eps)) then
!                tr = nint(msh(j,1))
!                else if(((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
!                ((abs(msh(j,3) - crdy)) .lt. eps))then
!                tl = nint(msh(j,1))
!                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - plus_y)) .lt. eps)) then
!                tt = nint(msh(j,1))
!                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_y)) .lt. eps)) then
!                td = nint(msh(j,1))
!                end if
!        enddo
!!
!!        write(*,*) 'coordenada do nó central', tempxy
!!        write(*,*) 'coordenadas do nós adjacentes', tl, td, tr, tt
!!
!        do c = 1, nodes
!                do d = 1, nodes
!                if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 4.0d0/(1.0d0*dxdy**2)
!                else if ((c .eq. tempxy) .and. (d .eq. td)) then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!                else if ((c .eq. tempxy) .and. (d .eq. tr)) then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!                else if ((c .eq. tempxy) .and. (d .eq. tt)) then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!                else if ((c .eq. tempxy) .and. (d .eq. tl)) then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!                end if
!                enddo
!        enddo
!        do kt = 1, nodes
!                if(kt .eq. tempxy)then
!                rhs(kt) = 0.0d0
!                end if
!        end do
!!
!end if
end do
!
return
end subroutine create_a_steady
!********************************************************************************
subroutine create_4_steady(msh, t, nodes, alin, rhs, dxdy)
!
!This is a subroutine to generate the matrix of equations and the rhs
!for the numerical solver, this require the mesh array!
!
!       Programmer: Paulo Conci         Last Modify: August 29
!
implicit none
! Variables o main prog
    integer*4, intent(in) :: nodes
    real(kind=8), intent(in) :: dxdy
    real*8, dimension(nodes,4), intent(in) :: msh
    real*8, dimension(nodes), intent(in) :: t
    real*8, dimension(nodes,nodes), intent(out) :: alin
    real*8, dimension(nodes), intent(out) :: rhs
!
!
!Variables of subroutine
    integer*4 :: i, j, c, d, kt
    integer*4 :: td, tl, tr, tt, tempxy, trt, trd, tlt, tld
    real*8 :: ymax, xmax, crdx, crdy, eps, temp
    real*8 :: plus_x, minus_x, plus_y, minus_y
    real(kind=8) :: cdx, cdy, temprhs
!
!
!finding max values of matrix of mesh
        xmax = maxval(msh(1:nodes,2))
        ymax = maxval(msh(1:nodes,3))
!
!
        eps = dxdy/3.0
!
        alin = 0.0d0
!
!
!
 do i = 1, nodes
!=====================================1========================================
! x = xmax and y=/ 0 and y =/ ymax
    if ((msh(i,2) .eq. xmax) .and. msh(i,3) .ne. 0.0 .and. &
        msh(i,3) .ne. ymax) then
            tempxy = nint(msh(i,1))
            temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    endif
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) =  temp
                endif
        enddo
!
!
!
! isolated nodes implemetation in y
! *******************************************************************
!        minus_x = msh(i,2) - dxdy
!        minus_y = msh(i,3) - dxdy
!        plus_y = msh(i,3) + dxdy
!        crdx = msh(i,2)
!        crdy = msh(i,3)
!        tempxy = nint(msh(i,1))
!
!        tt = 0
!        tl = 0
!        td = 0
!
!    do j = 1, nodes
!        if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!            ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                tt = nint(msh(j,1))
!        else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!            ((abs(msh(j,2) - minus_x)) .lt. eps))then
!                tl = nint(msh(j,1))
!        else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_y)) .lt. eps))then
!                td = nint(msh(j,1))
!        end if
!    end do
!
!
!    do c = 1, nodes
!        do d = 1, nodes
!            if((c .eq. tempxy) .and. (d .eq. tt)) then
!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tl))then
!               alin(c,d) = -1.0d0/(1.0d0*l**2)
!            else if ((c .eq. tempxy) .and. (d .eq. td)) then
!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 2.0d0/l**2
!            end if
!        enddo
!    enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) = 0.0d0
!                end if
!        end do
!=====================================2========================================
!x = 0  and y =/0 and y =/ ymax
    else if ((msh(i,2) .eq. 0.0) .and. (msh(i,3) .ne. ymax) &
        .and. msh(i,3) .ne. 0.0) then
!
!        plus_x = msh(i,2) + dxdy
!        minus_y = msh(i,3) - dxdy
!        plus_y = msh(i,3) + dxdy
!        crdx = msh(i,2)
!        crdy = msh(i,3)
!        tempxy = nint(msh(i,1))
!        temp = t(i)
!Implementation that do not change the temperature
!********************************************************************
        tempxy = nint(msh(i,1))
        temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                if((c .eq. tempxy).and.(d .eq. tempxy))then
                    alin(c,d) = 1.0d0
                end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!********************************************************************
!        tt = 0
!        tr = 0
!        td = 0
!!
!    do j = 1, nodes
!        if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!            ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                tt = nint(msh(j,1))
!!
!        else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!            ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                tr = nint(msh(j,1))
!!
!        else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_y)) .lt. eps))then
!                td = nint(msh(j,1))
!!
!        end if
!    end do
!!
!!
!    do c = 1, nodes
!        do d = 1, nodes
!            if((c .eq. tempxy) .and. (d .eq. tt))then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. td)) then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 2.0d0/dxdy**2
!            end if
!        enddo
!    enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                    rhs(kt) =  - temp + 1.0
!                end if
!        end do
!
!
!=====================================3========================================
! y = ymax and x=/0 and x=/ max
    else if ((msh(i,3) .eq. ymax) .and. (msh(i,2) .ne. xmax) .and. &
            msh(i,2) .ne. 0.0) then
! na parte superior
!
!
    tempxy = nint (msh(i,1))
    temp = t(i)
        do c = 1, nodes
            do d = 1, nodes
                if((c .eq. tempxy).and.(d .eq. tempxy))then
                    alin(c,d) = 1.0d0
                end if
            enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!
!
!                plus_x = msh(i,2) + dxdy
!                minus_x = msh(i,2) - dxdy
!                minus_y = msh(i,3) - dxdy
!                crdx = msh(i,2)
!                crdy = msh(i,3)
!                tempxy = nint(msh(i,1))
!                temp = t(i)
!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy).and.(d .eq. tempxy))then
!                alin(c,d) = 1.0d0
!                end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                rhs(j) = temp
!                end if
!        enddo
!
!        write(*,*) 'no',tempxy
!
!                tt = 0
!                tr = 0
!                tl = 0

!        do j = 1, nodes
!                if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_y)) .lt. eps))then
!                td = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                tr = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - minus_x)) .lt. eps))then
!                tl = nint(msh(j,1))
!                end if
!        end do



!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy) .and. (d .eq. td)) then
!                alin(c,d) = 2.0d0*fo
!                else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                alin(c,d) = fo
!                else if ((c .eq. tempxy) .and. (d .eq. tl)) then
!                alin(c,d) = fo
!                else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = a2
!                end if
!                enddo
!        enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) = b2 - t(kt)
!                end if
!        end do

!=====================================4========================================
! x =/ 0 and x =/ xmax and y = 0
else if ((msh(i,2) .ne. 0.0) .and. (msh(i,2) .ne. xmax) .and.&
        (msh(i,3) .eq. 0.0)) then ! nós isolados na parte inferior
!
!                plus_x = msh(i,2) + dxdy
!                minus_x = msh(i,2) - dxdy
!                plus_y = msh(i,3) + dxdy
!                crdx = msh(i,2)
!                crdy = msh(i,3)
                tempxy = nint(msh(i,1))
                temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                if((c .eq. tempxy).and.(d .eq. tempxy))then
                    alin(c,d) = 1.0d0
                end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo

!        write(*,*) 'no',tempxy
!
!                tt = 0
!                tr = 0
!                tl = 0
!!
!        do j = 1, nodes
!                if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                    tt = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                    tr = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - minus_x)) .lt. eps))then
!                    tl = nint(msh(j,1))
!                end if
!        end do
!!
!!
!!
!        do c = 1, nodes
!        do d = 1, nodes
!            if((c .eq. tempxy) .and. (d .eq. tt))then
!                alin(c,d) = -1.0d0/(1.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tl)) then
!                alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 2.0d0/dxdy**2
!            end if
!        enddo
!        enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) =  - temp + 1.0
!                end if
!        end do
!=====================================5========================================
! x= 0 and y= 0
elseif ((msh(i,2) .eq. 0.0) .and. msh(i,3) .eq. 0.0)then
        tempxy = nint(msh(i,1))
        temp = t(i)
!        crdx = msh(i,2)
!        crdy = msh(i,3)
!        plus_x = msh(i,2) + dxdy
!        plus_y = msh(i,3) + dxdy
!
!
!        tr = 0
!        tt = 0
!!
!!
!        do j = 1, nodes
!                if (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
!                ((abs(msh(j,3) - crdy)) .lt. eps)) then
!                    tr = nint(msh(j,1))
!                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - plus_y)) .lt. eps)) then
!                    tt = nint(msh(j,1))
!                end if
!        enddo
!!
!!        write(*,*) 'coordenada do nó central', tempxy
!!        write(*,*) 'coordenadas do nós adjacentes', tl, td, tr, tt
!!
!        do c = 1, nodes
!                do d = 1, nodes
!                    if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                        alin(c,d) = 1.0d0/dxdy**2
!                    else if ((c .eq. tempxy) .and. (d .eq. tr)) then
!                        alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!                    else if ((c .eq. tempxy) .and. (d .eq. tt)) then
!                        alin(c,d) = -1.0d0/(2.0d0*dxdy**2)
!                    end if
!                enddo
!        enddo
!        do kt = 1, nodes
!                if(kt .eq. tempxy)then
!                    rhs(kt) = - temp + 1.0
!                end if
!        end do
!
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!=====================================6========================================
!x = xmax and y=0
elseif((msh(i,2) .eq. xmax) .and. (msh(i,3) .eq. 0.0))then
    tempxy = nint(msh(i,1))
    temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) =  temp
                end if
        enddo
!=====================================7========================================
! x = xmax and y = ymax
elseif((msh(i,2) .eq. xmax) .and. msh(i,3) .eq. ymax)then
    tempxy = nint(msh(i,1))
    temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!=====================================8========================================
!x = 0 and y = ymax
elseif ((msh(i,2) .eq. 0.0) .and. (msh(i,3) .eq. ymax))then
    tempxy = nint(msh(i,1))
    temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) = temp
                end if
        enddo
!=====================================i========================================
! interior nodes !!!!
else if ((msh(i,2) .ne. 0.0) .and. (msh(i,2) .ne. xmax) .and.&
        (msh(i,3) .ne. 0.0) .and. (msh(i,3) .ne. ymax))then
                crdx = msh(i,2)
                crdy = msh(i,3)
                plus_x = msh(i,2) + dxdy
                minus_x = msh(i,2) - dxdy
                plus_y = msh(i,3) + dxdy
                minus_y = msh(i,3) - dxdy
                tempxy = nint(msh(i,1))
!
!
        tr = 0
        tl = 0
        tt = 0
        td = 0
        trt = 0
        trd = 0
        tlt = 0
        tld = 0
!
!
        do j = 1, nodes
                if (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - crdy)) .lt. eps)) then
                    tr = nint(msh(j,1))
                else if(((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - crdy)) .lt. eps))then
                    tl = nint(msh(j,1))
                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
                ((abs(msh(j,3) - plus_y)) .lt. eps)) then
                    tt = nint(msh(j,1))
                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
                ((abs(msh(j,3) - minus_y)) .lt. eps)) then
                    td = nint(msh(j,1))
                elseif(((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - plus_y)) .lt. eps)) then
                    trt = msh(j,1)
                elseif(((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - minus_y)) .lt. eps)) then
                    trd = msh(j,1)
                elseif (((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - minus_y)) .lt. eps)) then
                    tld = msh(j,1)
                elseif(((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - plus_y)) .lt. eps)) then
                    tlt = msh(j,1)
                end if
        enddo
!
!        write(*,*) 'coordenada do nó central', tempxy
!        write(*,*) 'coordenadas do nós adjacentes', tl, td, tr, tt, trt, trd, tld, tlt
!
       do c = 1, nodes
                do d = 1, nodes
                if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
                alin(c,d) = 20.0d0/(6.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. td)) then
                alin(c,d) = -4.0d0/(6.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. tr)) then
                alin(c,d) = -4.0d0/(6.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. tt)) then
                alin(c,d) = -4.0d0/(6.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. tl)) then
                alin(c,d) = -4.0d0/(6.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. trt)) then
                alin(c,d) = -1.0d0/(6.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. trd)) then
                alin(c,d) = -1.0d0/(6.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. tlt)) then
                alin(c,d) = -1.0d0/(6.0d0*dxdy**2)
                else if ((c .eq. tempxy) .and. (d .eq. tld)) then
                alin(c,d) = -1.0d0/(6.0d0*dxdy**2)
                end if
                enddo
        enddo
        do kt = 1, nodes
                if(kt .eq. tempxy)then
                    rhs(kt) = 0.0d0
                end if
        end do
!
end if
end do
!
return
end subroutine create_4_steady
!********************************************************************************
subroutine cte_gen(fo, a1, b1, a2, b2, b3, b4, b5, b6, dxdy)
use physics
use parmesh
!
!All constants should be declared here, this subroutine will store
!and calculate all constantes used in another subroutines
!
implicit none
!
! Constants
!
    real(kind=8) :: dxdy
    real*8,intent(out) :: fo, a1, a2, b1, b2, b3, b4, b5, b6
!
    real*8 :: rho_c
!
! Defining some constants
!
        Fo = (alpha*deltat)/dxdy**2
        rho_c = (k/alpha)
        a1 = (-1.0d0 - 4.0d0*Fo)
        b1 = -((alpha*deltat/k)*e)
        a2 = -(1.0 + 4.0d0*fo + 2.0d0*deltat*alpha*h/k*dxdy)
        b2 = -(2.0*deltat*alpha*h*temp_inf/k*dxdy + alpha*deltat*e/k)
        b3 = -(deltat*e/rho_c)
        b4 = 1.d0/(deltat*k)
        b5 = e/k
        b6 = 1.d0/(deltat*k)
!
!
!
!
return
end subroutine cte_gen
!********************************************************************************
subroutine create_a_unsteady(msh, t, nodes, alin, rhs, dxdy)
use parmesh
!
!This is a subroutine to generate the matrix of equations and the rhs
!for the numerical solver, this require the mesh array!
!
implicit none
! Variables o main prog
    integer*4, intent(in) :: nodes
    real(kind=8) :: dxdy
    real*8, dimension(nodes,4), intent(in) :: msh
    real*8, dimension(nodes), intent(in) :: t
    real*8, dimension(nodes,nodes), intent(out) :: alin
    real*8, dimension(nodes), intent(out) :: rhs
!
!
!Variables of subroutine
    integer*4 :: i, j, c, d, kt
    integer*4 :: td, tl, tr, tt, tempxy
    real*8 :: ymax, xmax, crdx, crdy, eps, temp
    real*8 :: plus_x, minus_x, plus_y, minus_y
    real(kind=8) :: cdx, cdy, temprhs
    real*8 :: fo, a1, b1, a2, b2, b3, b4, b5, b6
!
!    square domain and boundary conditions
!       8  ________3___________ 7
!         |                   |
!         |                   |
!       2 |        i          | 1
!         |                   |
!       5 |___________________| 6
!                  4
!
!
!
!finding max values of matrix of mesh
        xmax = maxval(msh(1:nodes,2))
        ymax = maxval(msh(1:nodes,3))
!
      call cte_gen(fo, a1, b1, a2, b2, b3, b4, b5, b6, dxdy)
!
!
        eps = dxdy/2.0
!
        alin = 0.0d0
!
    do i = 1, nodes
!=====================================1========================================
! x = xmax and y=/ 0 and y =/ ymax
    if ((msh(i,2) .eq. xmax) .and. msh(i,3) .ne. 0.0 .and. &
        msh(i,3) .ne. ymax) then
            tempxy = nint(msh(i,1))
            temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    endif
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) =  temp
                endif
        enddo
!
!
!
! isolated nodes implemetation in y
! *******************************************************************
!        minus_x = msh(i,2) - dxdy
!        minus_y = msh(i,3) - dxdy
!        plus_y = msh(i,3) + dxdy
!        crdx = msh(i,2)
!        crdy = msh(i,3)
!        tempxy = nint(msh(i,1))
!
!        tt = 0
!        tl = 0
!        td = 0
!
!    do j = 1, nodes
!        if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!            ((abs(msh(j,3) - plus_y)) .lt. eps))then
!                tt = nint(msh(j,1))
!        else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!            ((abs(msh(j,2) - minus_x)) .lt. eps))then
!                tl = nint(msh(j,1))
!        else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_y)) .lt. eps))then
!                td = nint(msh(j,1))
!        end if
!    end do
!
!
!    do c = 1, nodes
!        do d = 1, nodes
!            if((c .eq. tempxy) .and. (d .eq. tt)) then
!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tl))then
!               alin(c,d) = -1.0d0/(1.0d0*l**2)
!            else if ((c .eq. tempxy) .and. (d .eq. td)) then
!                alin(c,d) = -1.0d0/(2.0d0*l**2)
!            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = 2.0d0/l**2
!            end if
!        enddo
!    enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) = 0.0d0
!                end if
!        end do
!=====================================2========================================
!x = 0  and y =/0 and y =/ ymax
    else if ((msh(i,2) .eq. 0.0) .and. (msh(i,3) .ne. ymax) &
        .and. msh(i,3) .ne. 0.0) then
!
        plus_x = msh(i,2) + dxdy
        minus_y = msh(i,3) - dxdy
        plus_y = msh(i,3) + dxdy
        crdx = msh(i,2)
        crdy = msh(i,3)
        tempxy = nint(msh(i,1))
        temp = t(i)
!Implementation that do not change the temperature
!********************************************************************
!        tempxy = nint(msh(i,1))
!        temp = t(i)
!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy).and.(d .eq. tempxy))then
!                alin(c,d) = 1.0d0
!                end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                rhs(j) = temp
!                end if
!        enddo
!********************************************************************
        tt = 0
        tr = 0
        td = 0
!
    do j = 1, nodes
        if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
            ((abs(msh(j,3) - plus_y)) .lt. eps))then
                tt = nint(msh(j,1))
!
        else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
            ((abs(msh(j,2) - plus_x)) .lt. eps))then
                tr = nint(msh(j,1))
!
        else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
                ((abs(msh(j,3) - minus_y)) .lt. eps))then
                td = nint(msh(j,1))
!
        end if
    end do
!
!
    do c = 1, nodes
        do d = 1, nodes
            if((c .eq. tempxy) .and. (d .eq. tt)) then
                alin(c,d) = fo
            else if ((c .eq. tempxy) .and. (d .eq. tr))then
                alin(c,d) = 2.0*fo
            else if ((c .eq. tempxy) .and. (d .eq. td)) then
                alin(c,d) = fo
            else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
                alin(c,d) = a1
            end if
        enddo
    enddo
    do kt = 1, nodes
            if (kt .eq. tempxy)then
                rhs(kt) = - t(kt) + b1
            end if
    end do
!
!
!=====================================3========================================
! y = ymax and x=/0 and x=/ max
    else if ((msh(i,3) .eq. ymax) .and. (msh(i,2) .ne. xmax) .and. &
            msh(i,2) .ne. 0.0) then
! na parte superior
!
!
    tempxy = nint (msh(i,1))
    temp = t(i)
        do c = 1, nodes
            do d = 1, nodes
                if((c .eq. tempxy).and.(d .eq. tempxy))then
                    alin(c,d) = 1.0d0
                end if
            enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                rhs(j) =  temp
                end if
        enddo
!
!
!                plus_x = msh(i,2) + dxdy
!                minus_x = msh(i,2) - dxdy
!                minus_y = msh(i,3) - dxdy
!                crdx = msh(i,2)
!                crdy = msh(i,3)
!                tempxy = nint(msh(i,1))
!                temp = t(i)
!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy).and.(d .eq. tempxy))then
!                alin(c,d) = 1.0d0
!                end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                rhs(j) = temp
!                end if
!        enddo
!
!        write(*,*) 'no',tempxy
!
!                tt = 0
!                tr = 0
!                tl = 0

!        do j = 1, nodes
!                if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
!                ((abs(msh(j,3) - minus_y)) .lt. eps))then
!                td = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - plus_x)) .lt. eps))then
!                tr = nint(msh(j,1))
!                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
!                ((abs(msh(j,2) - minus_x)) .lt. eps))then
!                tl = nint(msh(j,1))
!                end if
!        end do



!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy) .and. (d .eq. td)) then
!                alin(c,d) = 2.0d0*fo
!                else if ((c .eq. tempxy) .and. (d .eq. tr))then
!                alin(c,d) = fo
!                else if ((c .eq. tempxy) .and. (d .eq. tl)) then
!                alin(c,d) = fo
!                else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
!                alin(c,d) = a2
!                end if
!                enddo
!        enddo
!        do kt = 1, nodes
!                if (kt .eq. tempxy)then
!                rhs(kt) = b2 - t(kt)
!                end if
!        end do

!=====================================4========================================
! x =/ 0 and x =/ xmax and y = 0
else if ((msh(i,2) .ne. 0.0) .and. (msh(i,2) .ne. xmax) .and.&
        (msh(i,3) .eq. 0.0)) then ! nós isolados na parte inferior
!
                plus_x = msh(i,2) + dxdy
                minus_x = msh(i,2) - dxdy
                plus_y = msh(i,3) + dxdy
                crdx = msh(i,2)
                crdy = msh(i,3)
                tempxy = nint(msh(i,1))
                temp = t(i)
!        do c = 1, nodes
!                do d = 1, nodes
!                if((c .eq. tempxy).and.(d .eq. tempxy))then
!                alin(c,d) = 1.0d0
!                end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                rhs(j) = temp
!                end if
!        enddo
!
!        write(*,*) 'no',tempxy
!
                tt = 0
                tr = 0
                tl = 0
!
        do j = 1, nodes
                if(((abs(msh(j,2) - crdx)) .lt. eps) .and. &
                ((abs(msh(j,3) - plus_y)) .lt. eps))then
                    tt = nint(msh(j,1))
                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
                ((abs(msh(j,2) - plus_x)) .lt. eps))then
                    tr = nint(msh(j,1))
                else if (((abs(msh(j,3) - crdy)) .lt. eps) .and. &
                ((abs(msh(j,2) - minus_x)) .lt. eps))then
                    tl = nint(msh(j,1))
                end if
        end do
!
!
!
        do c = 1, nodes
                do d = 1, nodes
                if((c .eq. tempxy) .and. (d .eq. tt)) then
                    alin(c,d) = 2.0d0*fo
                else if ((c .eq. tempxy) .and. (d .eq. tr))then
                    alin(c,d) = fo
                else if ((c .eq. tempxy) .and. (d .eq. tl)) then
                    alin(c,d) = fo
                else if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
                    alin(c,d) = a1
                end if
                enddo
        enddo
        do kt = 1, nodes
                if (kt .eq. tempxy)then
                    rhs(kt) = b1 - t(kt)
                end if
        end do
!=====================================5========================================
! x= 0 and y= 0
elseif ((msh(i,2) .eq. 0.0) .and. msh(i,3) .eq. 0.0)then
        tempxy = nint(msh(i,1))
        temp = t(i)
        crdx = msh(i,2)
        crdy = msh(i,3)
        plus_x = msh(i,2) + dxdy
        plus_y = msh(i,3) + dxdy
!
!
        tr = 0
        tt = 0
!
!
        do j = 1, nodes
                if (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - crdy)) .lt. eps)) then
                    tr = nint(msh(j,1))
                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
                ((abs(msh(j,3) - plus_y)) .lt. eps)) then
                    tt = nint(msh(j,1))
                end if
        enddo
!
!        write(*,*) 'coordenada do nó central', tempxy
!        write(*,*) 'coordenadas do nós adjacentes', tl, td, tr, tt
!
        do c = 1, nodes
                do d = 1, nodes
                    if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
                        alin(c,d) = a1
                    else if ((c .eq. tempxy) .and. (d .eq. tr)) then
                        alin(c,d) = 2.d0*fo
                    else if ((c .eq. tempxy) .and. (d .eq. tt)) then
                        alin(c,d) = 2.d0*fo
                    end if
                enddo
        enddo
        do kt = 1, nodes
                if(kt .eq. tempxy)then
                    rhs(kt) = b3 - t(kt)
                end if
        end do
!
!        do c = 1, nodes
!                do d = 1, nodes
!                    if((c .eq. tempxy).and.(d .eq. tempxy))then
!                        alin(c,d) = 1.0d0
!                    end if
!                enddo
!        enddo
!!
!        do j= 1, nodes
!                if(j .eq. tempxy)then
!                    rhs(j) = temp
!                end if
!        enddo
!=====================================6========================================
!x = xmax and y=0
elseif((msh(i,2) .eq. xmax) .and. (msh(i,3) .eq. 0.0))then
    tempxy = nint(msh(i,1))
    temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) =  temp
                end if
        enddo
!=====================================7========================================
! x = xmax and y = ymax
elseif((msh(i,2) .eq. xmax) .and. msh(i,3) .eq. ymax)then
    tempxy = nint(msh(i,1))
    temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) =  temp
                end if
        enddo
!=====================================8========================================
!x = 0 and y = ymax
elseif ((msh(i,2) .eq. 0.0) .and. (msh(i,3) .eq. ymax))then
    tempxy = nint(msh(i,1))
    temp = t(i)
        do c = 1, nodes
                do d = 1, nodes
                    if((c .eq. tempxy).and.(d .eq. tempxy))then
                        alin(c,d) = 1.0d0
                    end if
                enddo
        enddo
!
        do j= 1, nodes
                if(j .eq. tempxy)then
                    rhs(j) =  temp
                end if
        enddo
!=====================================i========================================
! interior nodes !!!!
else if ((msh(i,2) .ne. 0.0) .and. (msh(i,2) .ne. xmax) .and.&
        (msh(i,3) .ne. 0.0) .and. (msh(i,3) .ne. ymax))then
                crdx = msh(i,2)
                crdy = msh(i,3)
                plus_x = msh(i,2) + dxdy
                minus_x = msh(i,2) - dxdy
                plus_y = msh(i,3) + dxdy
                minus_y = msh(i,3) - dxdy
                tempxy = nint(msh(i,1))
!
!
        tr = 0
        tl = 0
        tt = 0
        td = 0
!
!
        do j = 1, nodes
                if (((abs(msh(j,2) - plus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - crdy)) .lt. eps)) then
                    tr = nint(msh(j,1))
                else if(((abs(msh(j,2) - minus_x)) .lt. eps) .and. &
                ((abs(msh(j,3) - crdy)) .lt. eps))then
                    tl = nint(msh(j,1))
                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
                ((abs(msh(j,3) - plus_y)) .lt. eps)) then
                    tt = nint(msh(j,1))
                else if (((abs(msh(j,2) - crdx)) .lt. eps) .and. &
                ((abs(msh(j,3) - minus_y)) .lt. eps)) then
                    td = nint(msh(j,1))
                end if
        enddo
!
!        write(*,*) 'coordenada do nó central', tempxy
!        write(*,*) 'coordenadas do nós adjacentes', tl, td, tr, tt
!
        do c = 1, nodes
                do d = 1, nodes
                if ((c .eq. tempxy) .and. (d .eq. tempxy)) then
                alin(c,d) = a1
                else if ((c .eq. tempxy) .and. (d .eq. td)) then
                alin(c,d) = fo
                else if ((c .eq. tempxy) .and. (d .eq. tr)) then
                alin(c,d) = fo
                else if ((c .eq. tempxy) .and. (d .eq. tt)) then
                alin(c,d) = fo
                else if ((c .eq. tempxy) .and. (d .eq. tl)) then
                alin(c,d) = fo
                end if
                enddo
        enddo
        do kt = 1, nodes
                if(kt .eq. tempxy)then
                rhs(kt) = b3 - t(kt)
                end if
        end do
!
end if
end do
!
return
end subroutine create_a_unsteady
!********************************************************************************
subroutine tcplt(msh, nodes, iounit, dxdy, dirname)
use parmesh
!this is a subroutine do ordering de mesh nodes for tecplot reading
!
implicit none
!
!
    character(*) :: dirname
    integer*4, intent(in) :: nodes
    real(kind=8) :: dxdy
    real*8, dimension(nodes,4), intent(in) :: msh
    character(*), intent(in) :: iounit
!
!
    integer*4 :: row, col, i, j, r
    real*8 ::  crdy, crdx
    real*8 :: eps, xmax, ymax
!
        eps = dxdy/2.0d0
!        msh(1:nodes, 4) = t(1:nodes)
        xmax = maxval(msh(1:nodes,2))
        ymax = maxval(msh(1:nodes,3))
!
!        write(*,*)'xm', xmax, 'ym', ymax
        row = (int(ymax/dxdy)) + 1
        col = (int(xmax/dxdy)) + 1
!
        open(unit=7, file=trim(adjustl(trim(dirname)//trim(iounit))), status='unknown')
    write(7,*) 'VARIABLES = "X", "Y", "T"'
    write(7,*) 'ZONE    I =', col, '    J =',  row,'    F=POINT'
        do j = 1, row
            do i = 1, col
            crdy = (dfloat(j-1))*dxdy
            crdx = (dfloat(i-1))*dxdy
               do r = 1, nodes
               if (((abs(msh(r,3) - crdy)) .lt. eps) .and. &
               ((abs(msh(r,2) - crdx)) .lt. eps)) then
               write(7,320) msh(r,2), msh(r,3), abs(msh(r,4))
               end if
               end do
            end do
        end do
320 format (3(F30.16))
!
        close(unit=7)
!
return
end subroutine tcplt
!********************************************************************************
subroutine pre_gmres(alin, nodes, a, ia, ja, nz_num)
!
!This is a subroutine to create the parameters to call gmres solver
!
implicit none
!
!Intent in/out variables
!
    integer(kind=4) :: nodes
    real(kind=8), dimension(nodes, nodes), intent(in) :: alin
    real(kind=8), allocatable, dimension(:), intent(out) :: a
    integer(kind=4), allocatable, dimension(:), intent(out) :: ia
    integer(kind=4), allocatable, dimension(:), intent(out) :: ja
    integer(kind=4), intent(out) :: nz_num
!
! Local variables
    integer(kind=4) :: i, j, k
!
        k = 0
        do i = 1, nodes
                do j = 1, nodes
                        if(alin(i,j) .ne. 0.0)then
                            k = k + 1
                        end if
                enddo
        enddo
        nz_num = k
!
!
        allocate(a(nz_num))
        allocate(ia(nz_num))
        allocate(ja(nz_num))
!
        k = 1
        do i = 1, nodes
                do j = 1, nodes
                        if((alin(i,j)) .ne. 0.0 )then
                            ia(k) = i
                            ja(k) = j
                            a(k) = alin(i,j)
                            k = k + 1
                        end if
                enddo
        enddo
!
return
end subroutine pre_gmres
!********************************************************************************
subroutine kicks(nodes, x)
!
!this subroutine will set the initial kick for the interative methods
implicit none
    integer(kind=4), intent(in) :: nodes
    real(kind=8), dimension(nodes), intent(out) :: x
!
!        if(opened(unit=80))then


        x(1:nodes) = 1.0d0
!
return
end subroutine kicks
!********************************************************************************
subroutine iluk(n, a, ja, ia, alu, jlu, ju, levs, iwk, w, jw, ierr)
use parsolvers
implicit none
integer*4 :: n
real*8 :: a(*),alu(*),w(*)
integer*4 :: ja(*),ia(*),jlu(*),ju(*),levs(*),jw(*),iwk,ierr
!c----------------------------------------------------------------------*
!c     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) *
!c----------------------------------------------------------------------*
!c
!c on entry:
!c==========
!c n       = integer. The row dimension of the matrix A. The matrix
!c
!c a,ja,ia = matrix stored in Compressed Sparse Row format.
!c
!c lfil    = integer. The fill-in parameter. Each element whose
!c           leve-of-fill exceeds lfil during the ILU process is dropped.
!c           lfil must be .ge. 0
!c
!c tol     = real*8. Sets the threshold for dropping small terms in the
!c           factorization. See below for details on dropping strategy.
!c
!c iwk     = integer. The minimum length of arrays alu, jlu, and levs.
!c
!c On return:
!c===========
!c
!c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!c           the L and U factors together. The diagonal (stored in
!c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!c           contains the i-th row of L (excluding the diagonal entry=1)
!c           followed by the i-th row of U.
!c
!c ju      = integer array of length n containing the pointers to
!c           the beginning of each row of U in the matrix alu,jlu.
!c
!c levs    = integer (work) array of size iwk -- which contains the
!c           levels of each element in alu, jlu.
!c
!c ierr    = integer. Error message with the following meaning.
!c           ierr  = 0    --> successful return.
!c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!c           ierr  = -1   --> Error. input matrix may be wrong.
!c                            (The elimination process has generated a
!c                            row in L or U whose length is .gt.  n.)
!c           ierr  = -2   --> The matrix L overflows the array al.
!c           ierr  = -3   --> The matrix U overflows the array alu.
!c           ierr  = -4   --> Illegal value for lfil.
!c           ierr  = -5   --> zero row encountered in A or U.
!c
!c work arrays:
!c=============
!c jw      = integer work array of length 3*n.
!c w       = real work array of length n
!c
!c Notes/known bugs: This is not implemented efficiently storage-wise.
!c       For example: Only the part of the array levs(*) associated with
!c       the U-matrix is needed in the routine.. So some storage can
!c       be saved if needed. The levels of fills in the LU matrix are
!c       output for information only -- they are not needed by LU-solve.
!c
!c----------------------------------------------------------------------
!c w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
!c jw(n+1:2n)  stores the nonzero indicator.
!c
!c Notes:
!c ------
!c All the diagonal elements of the input matrix must be  nonzero.
!c
!c----------------------------------------------------------------------*
!c     locals
      integer*4 :: ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,n2, &
                jlev, min
      real*8 :: t, s, fact
      if (lfil .lt. 0) goto 998
!c-----------------------------------------------------------------------
!c     initialize ju0 (points to next element to be added to alu,jlu)
!c     and pointer array.
!c-----------------------------------------------------------------------
      n2 = n+n
      ju0 = n+2
      jlu(1) = ju0
!c
!c     initialize nonzero indicator array + levs array --
!c
      do 1 j=1,2*n
         jw(j)  = 0
 1    continue
!c-----------------------------------------------------------------------
!c     beginning of main loop.
!c-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
!c
!c     unpack L-part and U-part of row of A in arrays w
!c
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
!c
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (t .eq. 0.0) goto 170
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n2+lenl) = 0
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
               jw(n2+ii) = 0
            else
               lenu = lenu+1
               jpos = ii+lenu-1
               jw(jpos) = k
               w(jpos) = t
               jw(n2+jpos) = 0
               jw(n+k) = jpos
            endif
 170     continue
!c
         jj = 0
!c
!c     eliminate previous rows
!c
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
!c-----------------------------------------------------------------------
!c     in order to do the elimination in the correct order we must select
!c     the smallest column index among jw(k), k=jj+1, ..., lenl.
!c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
!c
!c     determine smallest column index
!c
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
!c
         if (k .ne. jj) then
!c     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
!c     exchange in jw(n+  (pointers/ nonzero indicator).
            jw(n+jrow) = jj
            jw(n+j) = k
!c     exchange in jw(n2+  (levels)
            j = jw(n2+jj)
            jw(n2+jj)  = jw(n2+k)
            jw(n2+k) = j
!c     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
!c
!c     zero out element in row by resetting jw(n+jrow) to zero.
!c
         jw(n+jrow) = 0
!c
!c     get the multiplier for row to be eliminated (jrow) + its level
!c
         fact = w(jj)*alu(jrow)
         jlev = jw(n2+jj)
         if (jlev .gt. lfil) goto 150
!c
!c     combine current row and row jrow
!c
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
!c
!c     dealing with upper part.
!c
               if (jpos .eq. 0) then
!c
!c     this is a fill-in element
!c
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
                  jw(n2+i) = jlev+levs(k)+1
               else
!c
!c     this is not a fill-in element
!c
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            else
!c
!c     dealing with lower part.
!c
               if (jpos .eq. 0) then
!c
!c     this is a fill-in element
!c
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
                  jw(n2+lenl) = jlev+levs(k)+1
               else
!c
!c     this is not a fill-in element
!c
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            endif
 203     continue
         w(jj) = fact
         jw(jj)  = jrow
         goto 150
 160     continue
!c
!c     reset double-pointer to zero (U-part)
!c
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
!c
!c     update l-matrix
!c
         do 204 k=1, lenl
            if (ju0 .gt. iwk) goto 996
            if (jw(n2+k) .le. lfil) then
               alu(ju0) =  w(k)
               jlu(ju0) =  jw(k)
               ju0 = ju0+1
            endif
 204     continue
!c
!c     save pointer to beginning of row ii of U
!c
         ju(ii) = ju0
!c
!c     update u-matrix
!c
         do 302 k=ii+1,ii+lenu-1
            if (jw(n2+k) .le. lfil) then
               jlu(ju0) = jw(k)
               alu(ju0) = w(k)
               levs(ju0) = jw(n2+k)
               ju0 = ju0+1
            endif
 302     continue

         if (w(ii) .eq. 0.0) goto 999
!c
         alu(ii) = 1.0d0/ w(ii)
!c
!c     update pointer to beginning of next row of U.
!c
         jlu(ii+1) = ju0
!c-----------------------------------------------------------------------
!c     end main loop
!c-----------------------------------------------------------------------
 500  continue
      ierr = 0
      return
!c
!c     incomprehensible error. Matrix must be wrong.
!c
 995  ierr = -1
      return
!c
!c     insufficient storage in L.
!c
 996  ierr = -2
      return
!c
!c     insufficient storage in U.
!c
 997  ierr = -3
      return
!c
!c     illegal lfil entered.
!c
 998  ierr = -4
      return
!c
!c     zero row encountered in A or U.
!c
 999  ierr = -5
      return
!c----------------end-of-iluk--------------------------------------------
!c-----------------------------------------------------------------------
end subroutine iluk
!********************************************************************************
subroutine ilu0 ( n, a, ja, ia, alu, jlu, ju, iw, ierr )
!********************************************************************************
!
!! ILU0 is an ILU(0) preconditioner.
!
!  Discussion:
!
!    Note that this has been coded in such a way that it can be used
!    with PGMRES.  Normally, since the data structure of a, ja, ia is
!    the same as that of a, ja, ia, savings can be made. In fact with
!    some definitions (not correct for general sparse matrices) all we
!    need in addition to a, ja, ia is an additional diagonal.
!    Ilu0 is not recommended for serious problems. It is only provided
!    here for comparison purposes.
!
!    It is assumed that the the elements in the input matrix are stored
!    in such a way that in each row the lower part comes first and
!    then the upper part. To get the correct ILU factorization, it is
!    also necessary to have the elements of L sorted by increasing
!    column number. It may therefore be necessary to sort the
!    elements of a, ja, ia prior to calling ilu0. This can be
!    achieved by transposing the matrix twice using csrcsc.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju        = pointer to the diagonal elements in alu, jlu.
!
! ierr        = integer ( kind = 4 ) indicating error code on return
!           ierr = 0 --> normal return
!           ierr = k --> code encountered a zero pivot at step k.
! work arrays:
!
! iw          = integer ( kind = 4 ) work array of length n.
!
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) alu(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iw(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jf
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jlu(*)
  integer ( kind = 4 ) jm
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) js
  integer ( kind = 4 ) ju(*)
  integer ( kind = 4 ) ju0
  integer ( kind = 4 ) jw
  real ( kind = 8 ) tl

  ju0 = n + 2
  jlu(1) = ju0
!
!  Initialize the work vector.
!
  iw(1:n) = 0
!
!  The main loop.
!
  do ii = 1, n

    js = ju0
!
!  Generating row II of L and U.
!
    do j = ia(ii), ia(ii+1)-1
!
!  Copy row II of A, JA, IA into row II of ALU, JLU (L/U) matrix.
!
      jcol = ja(j)

      if ( jcol == ii ) then
        alu(ii) = a(j)
        iw(jcol) = ii
        ju(ii) = ju0
      else
        alu(ju0) = a(j)
        jlu(ju0) = ja(j)
        iw(jcol) = ju0
        ju0 = ju0 + 1
      end if

    end do

    jlu(ii+1) = ju0
    jf = ju0 - 1
    jm = ju(ii) - 1
!
!  Exit if the diagonal element is reached.
!
    do j = js, jm

      jrow = jlu(j)
      tl = alu(j) * alu(jrow)
      alu(j) = tl
!
!  Perform linear combination.
!
      do jj = ju(jrow), jlu(jrow+1)-1
        jw = iw(jlu(jj))
        if ( jw /= 0 ) then
          alu(jw) = alu(jw) - tl * alu(jj)
        end if
      end do

    end do
!
!  Invert and store the diagonal element.
!
    if ( alu(ii) == 0.0D+00 ) then
      ierr = ii
      return
    end if

    alu(ii) = 1.0D+00 / alu(ii)
!
!  Reset pointer IW to zero.
!
    iw(ii) = 0
    do i = js, jf
      iw(jlu(i)) = 0
    end do

  end do

  ierr = 0
  return
end subroutine ilu0
!********************************************************************************
subroutine sttocsr ( nrow, nnz, a, ir, jc, ao, jao, iao )
!******************************************************************************
!
!! COOCSR converts COO to CSR.
!
!  Discussion:
!
!    This routine converts a matrix that is stored in COO coordinate format
!    a, ir, jc into a CSR row general sparse ao, jao, iao format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NNZ, the number of nonzero elements.
!
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! on return:
!
! ir       is destroyed
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iad
  integer ( kind = 4 ) iao(nrow+1)
  integer ( kind = 4 ) ir(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jao(*)
  integer ( kind = 4 ) jc(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  integer ( kind = 4 ) nnz
  real ( kind = 8 ) x

  iao(1:nrow+1) = 0
!
!  Determine the row lengths.
!
  do k = 1, nnz
    iao(ir(k)) = iao(ir(k)) + 1
  end do
!
!  The starting position of each row.
!
  k = 1
  do j = 1, nrow+1
     k0 = iao(j)
     iao(j) = k
     k = k + k0
  end do
!
!  Go through the structure once more.  Fill in output matrix.
!
  do k = 1, nnz
     i = ir(k)
     j = jc(k)
     x = a(k)
     iad = iao(i)
     ao(iad) = x
     jao(iad) = j
     iao(i) = iad + 1
  end do
!
!  Shift back IAO.
!
  do j = nrow, 1, -1
    iao(j+1) = iao(j)
  end do
  iao(1) = 1
!
  return
end subroutine sttocsr
!*****************************************************************************80
subroutine csrmsr ( n, a, ja, ia, ao, jao, ju, wk, iwk )
!*****************************************************************************80
!
!! CSRMSR converts Compressed Sparse Row to Modified Sparse Row.
!
!  Discussion:
!
!    This routine converts a general sparse matrix a, ja, ia into
!    a compressed matrix using a separated diagonal (referred to as
!    the bell-labs format as it is used by bell labs semi conductor
!    group. We refer to it here as the modified sparse row format.
!
!    This has been coded in such a way that one can overwrite
!    the output matrix onto the input matrix if desired by a call of
!    the form
!
!     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!
!    In case ao, jao, are different from a, ja, then one can
!    use ao, jao as the work arrays in the calling sequence:
!
!     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!
!    Algorithm is in place.  i.e. both:
!
!          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!          (in which  ao, jao, are different from a, ja)
!           and
!          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!          (in which  wk, jwk, are different from a, ja)
!        are OK.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return :
!
! ao, jao  = sparse matrix in modified sparse row storage format:
!         +  ao(1:n) contains the diagonal of the matrix.
!         +  ao(n+2:nnz) contains the nondiagonal elements of the
!             matrix, stored rowwise.
!         +  jao(n+2:nnz) : their column indices
!         +  jao(1:n+1) contains the pointer array for the nondiagonal
!             elements in ao(n+1:nnz) and jao(n+2:nnz).
!             i.e., for i <= n+1 jao(i) points to beginning of row i
!            in arrays ao, jao.
!             here nnz = number of nonzero elements+1
!
!    Work array, real WK(N).
!
!    Work array, integer ( kind = 4 ) IWK(N+1).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iptr
  integer ( kind = 4 ) iwk(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jao(*)
  integer ( kind = 4 ) k
  real ( kind = 8 ) wk(n)
!
!*****************************************************************************80
! first of all creating the ju matrix
   integer ( kind = 4 ) ju(*)
   integer ( kind = 4 ) ju0
   integer ( kind = 4 ) ik
   integer ( kind = 4 ) jcol
!
   ju0 = n + 2
      do ik = 1, n
         do j= ia(ik), ia(ik+1)-1
            jcol = ja(j)
!
            if(jcol .eq. ii)then
            ju(ik) = ju0
            else
            ju0 = ju0 + 1
            end if
         enddo
      enddo
      j = 0
!*****************************************************************************80
!
!
!
  icount = 0
!
!  Store away diagonal elements and count nonzero diagonal elements.
!
  do i = 1, n
    wk(i) = 0.0D+00
    iwk(i+1) = ia(i+1) - ia(i)
    do k = ia(i), ia(i+1)-1
      if ( ja(k) == i ) then
        wk(i) = a(k)
        icount = icount + 1
        iwk(i+1) = iwk(i+1) - 1
      end if
    end do
  end do
!
!  Compute total length.
!
  iptr = n + ia(n+1) - icount
!
!  Copy backwards, to avoid collisions.
!
  do ii = n, 1, -1
    do k = ia(ii+1)-1, ia(ii), -1
      j = ja(k)
      if ( j /= ii ) then
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr - 1
      end if
    end do
  end do
!
!  Compute the pointer values and copy WK.
!
  jao(1) = n + 2
  do i = 1, n
    ao(i) = wk(i)
    jao(i+1) = jao(i) + iwk(i+1)
  end do

  return
end subroutine csrmsr
!*****************************************************************************80
subroutine bsort2 ( w, ind, n, ncut )
!*****************************************************************************80
!
!! BSORT2 returns the NCUT largest elements of an array, using bubble sort.
!
!  Discussion:
!
!    This routine carries out a simple bubble sort for getting the NCUT largest
!    elements in modulus, in array W.  IND is sorted accordingly.
!    (Ought to be replaced by a more efficient sort especially
!    if NCUT is not that small).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(*)
  integer ( kind = 4 ) iswp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ncut
  logical test
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) wswp

  i = 1

  do

    test = .false.

    do j = n-1, i, -1

      if ( abs ( w(j) ) < abs ( w(j+1) ) ) then
!
!  Swap.
!
        wswp = w(j)
        w(j) = w(j+1)
        w(j+1) = wswp
!
!  Reorder the original ind array accordingly.
!
        iswp = ind(j)
        ind(j) = ind(j+1)
        ind(j+1) = iswp
!
!  Set indicator that sequence is still unsorted.
!
        test = .true.

      end if

    end do

    i = i + 1

    if ( .not. test .or. ncut < i ) then
      exit
    end if

  end do

  return
end subroutine bsort2
!*****************************************************************************80
subroutine daxpy ( n, da, dx, incx, dy, incy )
!
!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    Uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Dongarra, Moler, Bunch, Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979.
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in DX and DY.
!
!    Input, real ( kind = 8 ) DA, the multiplier of DX.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries of DX.
!
!    Input/output, real ( kind = 8 ) DY(*), the second vector.
!    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive entries of DY.
!
  implicit none

  real ( kind = 8 ) da
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    return
  end if

  if ( da  == 0.0D+00 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dy(iy) = dy(iy) + da * dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 4 )

    do i = 1, m
      dy(i) = dy(i) + da * dx(i)
    end do

    do i = m+1, n, 4
      dy(i  ) = dy(i  ) + da * dx(i  )
      dy(i+1) = dy(i+1) + da * dx(i+1)
      dy(i+2) = dy(i+2) + da * dx(i+2)
      dy(i+3) = dy(i+3) + da * dx(i+3)
    end do

  end if

  return
end subroutine daxpy
!*****************************************************************************80
subroutine lusol0 ( n, y, x, alu, jlu, ju )
!*****************************************************************************
!
!! LUSOL0 performs a forward followed by a backward solve
! for LU matrix as produced by  ILUT
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real Y(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution.
!
!    ALU, JLU, JU, ...
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alu(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jlu(*)
  integer ( kind = 4 ) ju(*)
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
!
!  Forward solve
!
  do i = 1, n
    x(i) = y(i)
    do k = jlu(i), ju(i)-1
      x(i) = x(i) - alu(k) * x(jlu(k))
    end do
  end do
!
!  Backward solve.
!
  do i = n, 1, -1
    do k = ju(i), jlu(i+1)-1
      x(i) = x(i) - alu(k) * x(jlu(k))
    end do
    x(i) = alu(i) * x(i)
  end do

  return
end subroutine lusol0
!*****************************************************************************80
subroutine ope ( n, x, y, a, ja, ia )
!*****************************************************************************80
!
!! OPE sparse matrix * vector multiplication
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real Y(N), the product A * X.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) -1
    y(i) = 0.0D+00
    do k = k1, k2
      y(i) = y(i) + a(k) * x(ja(k))
    end do
  end do

  return
end subroutine ope
!*****************************************************************************80
subroutine pgmres( n, im, rhs, sol, vv, aa, ja, ia, alu, jlu, ju, fout, ierr, dirname)
use parsolvers
!******************************************************************************
!
!! PGMRES is an ILUT - Preconditioned GMRES solver.
!
!  Discussion:
!
!    This is a simple version of the ILUT preconditioned GMRES algorithm.
!    The ILUT preconditioner uses a dual strategy for dropping elements
!    instead  of the usual level of-fill-in approach. See details in ILUT
!    subroutine documentation. PGMRES uses the L and U matrices generated
!    from the subroutine ILUT to precondition the GMRES algorithm.
!    The preconditioning is applied to the right. The stopping criterion
!    utilized is based simply on reducing the residual norm by epsilon.
!    This preconditioning is more reliable than ilu0 but requires more
!    storage. It seems to be much less prone to difficulties related to
!    strong nonsymmetries in the matrix. We recommend using a nonzero tol
!    (tol=.005 or .001 usually give good results) in ILUT. Use a large
!    lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the
!    more reliable the code is. Efficiency may also be much improved.
!    Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as
!    Gaussian elimination without pivoting.
!
!    ILU(0) and MILU(0) are also provided for comparison purposes
!    USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and
!    then call pgmres.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IM, the size of the Krylov subspace.  IM
!    should not exceed 50 in this version.  This restriction can be reset by
!    changing the parameter command for KMAX below.
!
!    Input/output, real RHS(N), on input, the right hand side vector.
!    On output, the information in this vector has been destroyed.
!
! sol   == real vector of length n containing an initial guess to the
!          solution on input. approximate solution on output
!
! eps   == tolerance for stopping criterion. process is stopped
!          as soon as ( ||.|| is the euclidean norm):
!          || current residual||/||initial residual|| <= eps
!
! maxits== maximum number of iterations allowed
!
! iout  == output unit number number for printing intermediate results
!          if (iout <= 0) nothing is printed out.
!
!    Input, real AA(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!
! alu,jlu== A matrix stored in Modified Sparse Row format containing
!           the L and U factors, as computed by routine ilut.
!
! ju     == integer ( kind = 4 ) array of length n containing the pointers to
!           the beginning of each row of U in alu, jlu as computed
!           by routine ILUT.
!
! on return:
!
! sol   == contains an approximate solution (upon successful return).
! ierr  == integer ( kind = 4 ). Error message with the following meaning.
!          ierr = 0 --> successful return.
!          ierr = 1 --> convergence not achieved in itmax iterations.
!          ierr =-1 --> the initial guess seems to be the exact
!                       solution (initial residual computed was zero)
!
! work arrays:
!
! vv    == work array of length  n x (im+1) (used to store the Arnoli
!          basis)
  implicit none
  character(*) :: dirname
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: kmax = 1000
  real ( kind = 8 ) aa(*)
  real ( kind = 8 ) alu(*)
  real ( kind = 8 ) c(kmax)
  real ( kind = 8 ) eps1
  real ( kind = 8 ), parameter :: epsmac = 1.0D-16
  real ( kind = 8 ) gam
  real ( kind = 8 ) hh(kmax+1,kmax)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) im
  integer ( kind = 4 ) its
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jlu(*)
  integer ( kind = 4 ) ju(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) n1
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) ro
  real ( kind = 8 ) rs(kmax+1)
  real ( kind = 8 ) s(kmax)
  real ( kind = 8 ) sol(n)
  real ( kind = 8 ) t
  real ( kind = 8 ) vv(n,*)
!
   character(len=*) :: fout
!
!  Arnoldi size should not exceed KMAX=100 in this version.
!  To reset modify parameter KMAX accordingly.
!   kmax = n
!
    if (iout .ne. 0)then
        open(unit=iout, file=trim(adjustl(trim(dirname)//fout)), status='unknown')
    endif
  n1 = n + 1
  its = 0
!
!  Outer loop starts here.
!  Compute initial residual vector.
!
  call ope ( n, sol, vv, aa, ja, ia )

  vv(1:n,1) = rhs(1:n) - vv(1:n,1)

  do

    ro = sqrt ( ddot ( n, vv, 1, vv, 1 ) )

    if ( 0 < iout .and. its == 0 ) then
      write(iout, *) its, ro
      write(*,*) 'Iteration= ',its,'Residual=', ro
    end if

    if ( ro == 0.0D+00 ) then
      ierr = -1
      exit
    end if

    t = 1.0D+00 / ro
    vv(1:n,1) = vv(1:n,1) * t

    if ( its == 0 ) then
      eps1 = tol_abs * ro
    end if
!
!  Initialize first term of RHS of Hessenberg system.
!
     rs(1) = ro
     i = 0

 4   continue

     i = i + 1
     its = its + 1
     i1 = i + 1
     call lusol0 ( n, vv(1,i), rhs, alu, jlu, ju )
     call ope ( n, rhs, vv(1,i1), aa, ja, ia )
!
!  Modified Gram - Schmidt.
!
     do j = 1, i
       t = ddot ( n, vv(1,j), 1, vv(1,i1), 1 )
       hh(j,i) = t
       call daxpy ( n, -t, vv(1,j), 1, vv(1,i1), 1 )
     end do

     t = sqrt ( ddot ( n, vv(1,i1), 1, vv(1,i1), 1 ) )
     hh(i1,i) = t

     if ( t /= 0.0D+00 ) then
       t = 1.0D+00 / t
       vv(1:n,i1) = vv(1:n,i1) * t
     end if
!
!  Update factorization of HH.
!
    if ( i == 1 ) then
      go to 121
    end if
!
!  Perform previous transformations on I-th column of H.
!
    do k = 2, i
       k1 = k-1
       t = hh(k1,i)
       hh(k1,i) = c(k1) * t + s(k1) * hh(k,i)
       hh(k,i) = -s(k1) * t + c(k1) * hh(k,i)
    end do

121 continue

    gam = sqrt ( hh(i,i)**2 + hh(i1,i)**2 )
!
!  If GAMMA is zero then any small value will do.
!  It will affect only residual estimate.
!
    if ( gam == 0.0D+00 ) then
      gam = epsmac
    end if
!
!  Get the next plane rotation.
!
    c(i) = hh(i,i) / gam
    s(i) = hh(i1,i) / gam
    rs(i1) = -s(i) * rs(i)
    rs(i) = c(i) * rs(i)
!
!  Determine residual norm and test for convergence.
!
    hh(i,i) = c(i) * hh(i,i) + s(i) * hh(i1,i)
    ro = abs ( rs(i1) )
131 format(1h ,2e14.4)

    if ( 0 < iout ) then
        write(iout, *) its, ro
        write(*,*) 'Iteration= ',its,'Residual=', ro
    end if

    if ( i < im .and. eps1 < ro ) then
      go to 4
    end if
!
!  Now compute solution.  First solve upper triangular system.
!
    rs(i) = rs(i) / hh(i,i)

    do ii = 2, i
      k = i - ii + 1
      k1 = k + 1
      t = rs(k)
      do j = k1, i
        t = t - hh(k,j) * rs(j)
      end do
      rs(k) = t / hh(k,k)
    end do
!
!  Form linear combination of V(*,i)'s to get solution.
!
    t = rs(1)
    rhs(1:n) = vv(1:n,1) * t

    do j = 2, i
      t = rs(j)
      rhs(1:n) = rhs(1:n) + t * vv(1:n,j)
    end do
!
!  Call preconditioner.
!
    call lusol0 ( n, rhs, rhs, alu, jlu, ju )

    sol(1:n) = sol(1:n) + rhs(1:n)
!
!  Restart outer loop when necessary.
!
    if ( ro <= eps1 ) then
      ierr = 0
      exit
    end if

    if ( maxits < its ) then
      ierr = 1
      exit
    end if
!
!  Else compute residual vector and continue.
!
    do j = 1, i
      jj = i1 - j + 1
      rs(jj-1) = -s(jj-1) * rs(jj)
      rs(jj) = c(jj-1) * rs(jj)
    end do

    do j = 1, i1
      t = rs(j)
      if ( j == 1 ) then
        t = t - 1.0D+00
      end if
      call daxpy ( n, t, vv(1,j), 1,  vv, 1 )
    end do

199 format(' its =', i4, ' res. norm =', G14.6)

  end do
!
    if ( 0 < iout ) then
      write(iout, *) its, ro
      write ( *, '(a,i4,a,g14.6)' ) '  K = ', its, '  Residual = ', ro
    end if
!
   close(unit=iout)
  return
end subroutine pgmres
!*****************************************************************************80
subroutine ilut ( n, a, ja, ia, alu, jlu, ju, iwk, wu, wl, jr, jwl, jwu, ierr )
use parsolvers
!*****************************************************************************80
!
!! ILUT is an ILUT preconditioner.
!
!  Discussion:
!
!    This routine carries out incomplete LU factorization with dual
!    truncation mechanism.  Sorting is done for both L and U.
!
!    The dual drop-off strategy works as follows:
!
!    1) Theresholding in L and U as set by TOL.  Any element whose size
!       is less than some tolerance (relative to the norm of current
!       row in u) is dropped.
!
!    2) Keeping only the largest lenl0+lfil elements in L and the
!       largest lenu0+lfil elements in U, where lenl0=initial number
!       of nonzero elements in a given row of lower part of A
!       and lenlu0 is similarly defined.
!
!    Flexibility: one can use tol=0 to get a strategy based on keeping the
!    largest elements in each row of L and U. Taking tol /= 0 but lfil=n
!    will give the usual threshold strategy (however, fill-in is then
!    unpredictible).
!
!    A must have all nonzero diagonal elements.
!
!  Modified:
!
!    21 January 2015
!
!  Author:
!
!    Youcef Saad
!
!  Reference:
!
!    Youcef Saad,
!    Sparsekit: a basic tool kit for sparse matrix computations,
!    Technical Report, Computer Science Department,
!    University of Minnesota, June 1994
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(*), integer ( kind = 4 ) JA(*), IA(N+1),
!    the matrix in Compressed Sparse Row (CSR) format.
!
!    Input, integer ( kind = 4 ) LFIL, the fill-in parameter.  Each row of L
!    and each row of U will have a maximum of LFIL elements in addition to
!    the original number of nonzero elements.  Thus storage can be
!    determined beforehand.
!    0 <= LFIL.
!
!    Input, real ( kind = 8 ) TOL, the tolerance.
!
!    Output, real ( kind = 8 ) ALU(*), integer ( kind = 4 ) JLU(*),
!    the matrix stored in Modified Sparse Row (MSR) format, containing
!    the L and U factors together. The diagonal (stored in alu(1:n) ) is
!    inverted. Each I-th row of the ALU, JLU matrix contains the I-th row
!    of L (excluding the diagonal entry=1) followed by the I-th row of U.
!
!    Output, integer ( kind = 4 ) JU(N), pointers to the beginning of each
!    row of U in the matrix ALU, JLU.
!
!    Input, integer ( kind = 4 ) IWK, the minimum length of arrays ALU and
!    JLU.
!
!    Workspace, real ( kind = 8 ) WU(N+1), WL(N).
!
!    Workspace, integer ( kind = 4 ) JR(N), JWL(N), JWU(N).
!
!    Output, integer IERR. Error message with the following meaning.
!    ierr  = 0    --> successful return.
!    ierr > 0  --> zero pivot encountered at step number ierr.
!    ierr  = -1   --> Error. input matrix may be wrong.
!                     (The elimination process has generated a
!                     row in L or U whose length is >  n.)
!    ierr  = -2   --> The matrix L overflows the array alu.
!    ierr  = -3   --> The matrix U overflows the array alu.
!    ierr  = -4   --> Illegal value for lfil.
!    ierr  = -5   --> zero pivot encountered.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) alu(*)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) idiag
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iwk
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jlu(*)
  integer ( kind = 4 ) jpos
  integer ( kind = 4 ) jr(*)
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) ju(*)
  integer ( kind = 4 ) ju0
  integer ( kind = 4 ) jwl(n)
  integer ( kind = 4 ) jwu(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) len
  integer ( kind = 4 ) lenl
  integer ( kind = 4 ) lenl0
  integer ( kind = 4 ) lenu
  integer ( kind = 4 ) lenu0
!  integer ( kind = 4 ) lfil
  integer ( kind = 4 ) nl
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) tnorm
!  real ( kind = 8 ) tol
  real ( kind = 8 ) wl(n)
  real ( kind = 8 ) wu(n+1)

  if ( lfil < 0 ) then
    ierr = -4
    return
  end if
!
!  Initialize JU0 (points to next element to be added to ALU, JLU)
!  and pointer.
!
  ju0 = n + 2
  jlu(1) = ju0
!
!  integer ( kind = 4 ) double pointer array.
!
  jr(1:n) = 0
!
!  The main loop.
!
  do ii = 1, n

    j1 = ia(ii)
    j2 = ia(ii+1) - 1
    lenu = 0
    lenl = 0

    tnorm = 0.0D+00
    do k = j1, j2
      tnorm = tnorm + abs ( a(k) )
    end do
    tnorm = tnorm / real ( j2-j1+1, kind = 8 )
!
!  Unpack L-part and U-part of row of A in arrays WL, WU.
!
    do j = j1, j2

      k = ja(j)
      t = a(j)

      if ( tol * tnorm <= abs ( t ) ) then

        if ( k < ii ) then
          lenl = lenl + 1
          jwl(lenl) = k
          wl(lenl) = t
          jr(k) = lenl
        else
          lenu = lenu+1
          jwu(lenu) = k
          wu(lenu) = t
          jr(k) = lenu
        end if

      end if

    end do

    lenl0 = lenl
    lenu0 = lenu
    jj = 0
    nl = 0
!
!  Eliminate previous rows.
!
150 continue

    jj = jj + 1

    if ( lenl < jj ) then
      go to 160
    end if
!
!  In order to do the elimination in the correct order we need to
!  exchange the current row number with the one that has
!  smallest column number, among JJ, JJ+1, ..., LENL.
!
    jrow = jwl(jj)
    k = jj
!
!  Determine the smallest column index.
!
    do j = jj+1, lenl
       if ( jwl(j) < jrow ) then
          jrow = jwl(j)
          k = j
       end if
    end do
!
!  Exchange in JWL.
!
    j = jwl(jj)
    jwl(jj) = jrow
    jwl(k) = j
!
!  Exchange in JR.
!
    jr(jrow) = jj
    jr(j) = k
!
!  Exchange in WL.
!
    s = wl(k)
    wl(k) = wl(jj)
    wl(jj) = s

    if ( ii <= jrow ) then
      go to 160
    end if
!
!  Get the multiplier for row to be eliminated: JROW.
!
    fact = wl(jj) * alu(jrow)
    jr(jrow) = 0

    if ( abs ( fact ) * wu(n+2-jrow) <= tol * tnorm ) then
      go to 150
    end if
!
!  Combine current row and row JROW.
!
    do k = ju(jrow), jlu(jrow+1)-1
       s = fact * alu(k)
       j = jlu(k)
       jpos = jr(j)
!
!  If fill-in element and small disregard.
!
       if ( abs ( s ) < tol * tnorm .and. jpos == 0 ) then
         cycle
       end if

       if ( ii <= j ) then
!
!  Dealing with upper part.
!
          if ( jpos == 0 ) then
!
!  This is a fill-in element.
!
             lenu = lenu + 1

             if ( n < lenu ) then
               go to 995
             end if

             jwu(lenu) = j
             jr(j) = lenu
             wu(lenu) = - s
          else
!
!  No fill-in element.
!
             wu(jpos) = wu(jpos) - s
          end if
       else
!
!  Dealing with lower part.
!
          if ( jpos == 0 ) then
!
!  This is a fill-in element.
!
             lenl = lenl + 1

             if ( n < lenl ) then
               go to 995
             end if

             jwl(lenl) = j
             jr(j) = lenl
             wl(lenl) = -s
          else
!
!  No fill-in element.
!
             wl(jpos) = wl(jpos) - s
          end if
       end if

  end do

    nl = nl + 1
    wl(nl) = fact
    jwl(nl) = jrow
  go to 150
!
!  Update the L matrix.
!
 160 continue

    len = min ( nl, lenl0 + lfil )

    call bsort2 ( wl, jwl, nl, len )

    do k = 1, len

       if ( iwk < ju0 ) then
         ierr = -2
         return
       end if

       alu(ju0) =  wl(k)
       jlu(ju0) =  jwl(k)
       ju0 = ju0 + 1

    end do
!
!  Save pointer to beginning of row II of U.
!
    ju(ii) = ju0
!
!  Reset double pointer JR to zero (L-part - except first
!  JJ-1 elements which have already been reset).
!
  do k = jj, lenl
    jr(jwl(k)) = 0
  end do
!
!  Be sure that the diagonal element is first in W and JW.
!
    idiag = jr(ii)

    if ( idiag == 0 ) then
      go to 900
    end if

    if ( idiag /= 1 ) then

       s = wu(1)
       wu(j) = wu(idiag)
       wu(idiag) = s

       j = jwu(1)
       jwu(1) = jwu(idiag)
       jwu(idiag) = j

    end if

    len = min ( lenu, lenu0 + lfil )

    call bsort2 ( wu(2), jwu(2), lenu-1, len )
!
! Update the U-matrix.
!
    t = 0.0D+00

    do k = 2, len

       if ( iwk < ju0 ) then
         ierr = -3
         return
       end if

       jlu(ju0) = jwu(k)
       alu(ju0) = wu(k)
       t = t + abs ( wu(k) )
       ju0 = ju0 + 1

    end do
!
!  Save norm in WU (backwards). Norm is in fact average absolute value.
!
    wu(n+2-ii) = t / real ( len + 1, kind = 8 )
!
!  Store inverse of diagonal element of U.
!
    if ( wu(1) == 0.0D+00 ) then
      ierr = -5
      return
    end if

    alu(ii) = 1.0D+00 / wu(1)
!
!  Update pointer to beginning of next row of U.
!
  jlu(ii+1) = ju0
!
!  Reset double pointer JR to zero (U-part).
!
  do k = 1, lenu
    jr(jwu(k)) = 0
  end do

  end do

  ierr = 0

  return
!
!  Zero pivot :
!
 900    ierr = ii
    return
!
!  Incomprehensible error. Matrix must be wrong.
!
 995    ierr = -1
    return
end subroutine ilut
!*****************************************************************************80
subroutine milu0 ( n, a, ja, ia, alu, jlu, ju, iw, ierr )
!*****************************************************************************80
!
!! MILU0 is a simple milu(0) preconditioner.
!
!  Discussion:
!
!    Note that this has been coded in such a way that it can be used
!    with pgmres. Normally, since the data structure of a, ja, ia is
!    the same as that of a, ja, ia, savings can be made. In fact with
!    some definitions (not correct for general sparse matrices) all we
!    need in addition to a, ja, ia is an additional diagonal.
!    Ilu0 is not recommended for serious problems. It is only provided
!    here for comparison purposes.
!
!    It is assumed that the the elements in the input matrix are ordered
!    in such a way that in each row the lower part comes first and
!    then the upper part. To get the correct ILU factorization, it is
!    also necessary to have the elements of L ordered by increasing
!    column number. It may therefore be necessary to sort the
!    elements of a, ja, ia prior to calling milu0. This can be
!    achieved by transposing the matrix twice using csrcsc.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju        = pointer to the diagonal elements in alu, jlu.
!
!    Workspace, integer ( kind = 4 ) IW(N).
!
! ierr        = integer ( kind = 4 ) indicating error code on return
!           ierr = 0 --> normal return
!           ierr = k --> code encountered a zero pivot at step k.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) alu(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iw(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jf
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jlu(*)
  integer ( kind = 4 ) jm
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) js
  integer ( kind = 4 ) ju(*)
  integer ( kind = 4 ) ju0
  integer ( kind = 4 ) jw
  real ( kind = 8 ) s
  real ( kind = 8 ) tl

  ju0 = n + 2
  jlu(1) = ju0
!
!  Initialize work vector.
!
  iw(1:n) = 0
!
!  The main loop.
!
  do ii = 1, n

    js = ju0
!
!  Generating row II or L and U.
!
    do j = ia(ii), ia(ii+1)-1
!
!  Copy row II of A, JA, IA into row II of ALU, JLU (L/U) matrix.
!
      jcol = ja(j)

      if ( jcol == ii ) then
        alu(ii) = a(j)
        iw(jcol) = ii
        ju(ii) = ju0
      else
        alu(ju0) = a(j)
        jlu(ju0) = ja(j)
        iw(jcol) = ju0
        ju0 = ju0 + 1
      end if

    end do

    jlu(ii+1) = ju0
    jf = ju0 - 1
    jm = ju(ii) - 1
!
!  S accumulates fill-in values.
!
    s = 0.0D+00

    do j = js, jm

      jrow = jlu(j)
      tl = alu(j) * alu(jrow)
      alu(j) = tl
!
!  Perform linear combination.
!
      do jj = ju(jrow), jlu(jrow+1)-1
        jw = iw(jlu(jj))
        if ( jw /= 0 ) then
          alu(jw) = alu(jw) - tl * alu(jj)
        else
          s = s + tl * alu(jj)
        end if
      end do

    end do
!
!  Invert and store diagonal element.
!
    alu(ii) = alu(ii) - s

    if ( alu(ii) == 0.0D+00 ) then
      ierr = ii
      return
    end if

    alu(ii) = 1.0D+00 / alu(ii)
!
!  Reset pointer IW to zero.
!
    iw(ii) = 0
    do i = js, jf
      iw(jlu(i)) = 0
    end do

  end do

  ierr = 0

  return
end subroutine milu0
!*****************************************************************************80
subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )
!******************************************************************80
!
!! ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) W(N), the value of A'*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(ja(k1:k2)) = w(ja(k1:k2)) + a(k1:k2) * x(i)
  end do

  return
end subroutine atx_cr
!*****************************************************************************80
subroutine atx_st ( n, nz_num, ia, ja, a, x, w )
!*******************************************************************************
!
!! ATX_ST computes A'*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) W(N), the value of A'*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(j) = w(j) + a(k) * x(i)
  end do

  return
end subroutine atx_st
!*****************************************************************************80
subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )
!*******************************************************************************
!
!! AX_CR computes A*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) W(N), the value of A*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

  return
end subroutine ax_cr
!*****************************************************************************80
subroutine ax_st ( n, nz_num, ia, ja, a, x, w )
!*******************************************************************************
!
!! AX_ST computes A*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) W(N), the value of A*X.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(i) = w(i) + a(k) * x(j)
  end do

  return
end subroutine ax_st
!*****************************************************************************80
subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
!*******************************************************************************
!
!! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!    The array UA can be used to locate the diagonal elements of the matrix.
!
!    It is assumed that every row of the matrix includes a diagonal element,
!    and that the elements of each row have been ascending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!    On output, the order of the entries of JA may have changed because of
!    the sorting.
!
!    Output, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) ua(n)

  ua(1:n) = -1

  do i = 1, n
    do k = ia(i), ia(i+1) - 1
      if ( ja(k) == i ) then
        ua(i) = k
      end if
    end do
  end do

  return
end subroutine diagonal_pointer_cr
!*****************************************************************************80
subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )
!*******************************************************************************
!
!! ILU_CR computes the incomplete LU factorization of a matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iw(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) jw
  integer ( kind = 4 ) k
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) tl
  integer ( kind = 4 ) ua(n)
!
!  Copy A.
!
  l(1:nz_num) = a(1:nz_num)

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i+1) - 1
      iw(ja(k)) = k
    end do

    do j = ia(i), ia(i+1) - 1
      jrow = ja(j)
      if ( i <= jrow ) then
        exit
      end if
      tl = l(j) * l(ua(jrow))
      l(j) = tl
      do jj = ua(jrow) + 1, ia(jrow+1) - 1
        jw = iw(ja(jj))
        if ( jw /= -1 ) then
          l(jw) = l(jw) - tl * l(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow /= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( l(j) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    l(j) = 1.0D+00 / l(j)

  end do

  l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

  return
end subroutine ilu_cr
!*****************************************************************************80
subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )
!*******************************************************************************
!
!! LUS_CR applies the incomplete LU preconditioner.
!
!  Discussion:
!
!    The linear system M * Z = R is solved for Z.  M is the incomplete
!    LU preconditioner matrix, and R is a vector supplied by the user.
!    So essentially, we're solving L * U * Z = R.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) L(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Input, real ( kind = 8 ) R(N), the right hand side.
!
!    Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) z(n)
!
!  Copy R in.
!
  w(1:n) = r(1:n)
!
!  Solve L * w = w where L is unit lower triangular.
!
  do i = 2, n
    do j = ia(i), ua(i) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
  end do
!
!  Solve U * w = w, where U is upper triangular.
!
  do i = n, 1, -1
    do j = ua(i) + 1, ia(i+1) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
    w(i) = w(i) / l(ua(i))
  end do
!
!  Copy Z out.
!
  z(1:n) = w(1:n)

  return
end subroutine lus_cr
!*****************************************************************************80
subroutine mgmres_st ( n, nz_num, ia, ja, a, x, rhs, fout, dirname)
use parsolvers
!*****************************************************************************80
!
!! MGMRES_ST applies restarted GMRES to a sparse triplet matrix.
!
!  Discussion:
!
!    The linear system A*X=B is solved iteratively.
!
!    The matrix A is assumed to be stored in sparse triplet form.  Only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
!
!    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer)
!    iterations to take.
!
!    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations
!    to take.  0 < MR <= N.
!
!    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none
!
    character(*) :: dirname
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) av
  real ( kind = 8 ) c(1:mr)
  real ( kind = 8 ), parameter :: delta = 1.0D-03
  real ( kind = 8 ) g(1:mr+1)
  real ( kind = 8 ) h(1:mr+1,1:mr)
  real ( kind = 8 ) htmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) itr
  integer ( kind = 4 ) itr_used
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_copy
  real ( kind = 8 ) mu
  real ( kind = 8 ) r(1:n)
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_tol
  real ( kind = 8 ) rhs(1:n)
  real ( kind = 8 ) s(1:mr)
  real ( kind = 8 ) v(1:n,1:mr+1)
  logical, parameter :: verbose = .true.
  real ( kind = 8 ) x(1:n)
  real ( kind = 8 ) y(1:mr+1)
!
   character(len=*) :: fout
!
!
! printing the residual
    if(iout .ne. 0)then
        open(unit=iout, file=trim(adjustl(trim(dirname)//fout)), status='unknown')
    endif
  itr_used = 0

  if ( n < mr ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
    write ( *, '(a)' ) '  N < MR.'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  MR = ', mr
    stop
  end if

  do itr = 1, maxits

    call ax_st ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    rho = sqrt ( dot_product ( r(1:n), r(1:n) ) )

    if ( verbose ) then
    write ( *, '(a,i8,a,g14.6)' )' ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_st ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - h(j,k) * v(1:n,j)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( av + delta * h(k+1,k) == av ) then

        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do

        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then

        y(1:k+1) = h(1:k+1,k)

        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y(1:k+1) )
        end do

        h(1:k+1,k) = y(1:k+1)

      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )
      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g(1:k+1) )
      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
      write(iout,'(i8,g14.6)') itr_used, rho
      write ( *, '(a,i8,a,g14.6)' ) '  K =   ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)'       ) ' '
    write ( *, '(a)'       ) 'MGMRES_ST:'
    write ( *, '(a,i8)'    ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end subroutine mgmres_st
!*****************************************************************************80
subroutine mult_givens ( c, s, k, g )
!*****************************************************************************80
!
!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
!
!  Discussion:
!
!    In order to make it easier to compare this code with the Original C,
!    the vector indexing is 0-based.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens
!    rotation.
!
!    Input, integer ( kind = 4 ) K, indicates the location of the first
!    vector entry.
!
!    Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
!
  implicit none

  integer ( kind = 4 ) k

  real ( kind = 8 ) c
  real ( kind = 8 ) g(1:k+1)
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) s

  g1 = c * g(k) - s * g(k+1)
  g2 = s * g(k) + c * g(k+1)

  g(k)   = g1
  g(k+1) = g2

  return
end subroutine mult_givens
!*****************************************************************************80
subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs)
use parsolvers
!***********************************************************************
!
!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
!
!    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer)
!    iterations to take.
!
!    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations
!    to take.  MR must be less than N.
!
!    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none
!
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) av
  real ( kind = 8 ) c(mr+1)
  real ( kind = 8 ), parameter :: delta = 1.0D-03
  real ( kind = 8 ) g(mr+1)
  real ( kind = 8 ) h(mr+1,mr)
  real ( kind = 8 ) htmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) itr
  integer ( kind = 4 ) itr_used
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_copy
  real ( kind = 8 ) l(ia(n+1)+1)
  real ( kind = 8 ) mu
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_tol
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) s(mr+1)
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) v(n,mr+1);
  logical, parameter :: verbose = .true.
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(mr+1)
!
    if(iout .ne. 0)then
        open(unit=iout, file='residual.dat', status='unknown')
    endif
  itr_used = 0

  call rearrange_cr ( n, nz_num, ia, ja, a )

  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

  call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

  do itr = 1, maxits

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
    write(iout, *) itr_used, rho
!      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) )

      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
      write(iout,*) itr_used, rho
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do
!
    if ( verbose ) then
        write(iout,*) itr_used, rho
    endif
!
  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end subroutine pmgmres_ilu_cr
!*****************************************************************************80
subroutine rearrange_cr ( n, nz_num, ia, ja, a )
!*****************************************************************************80
!
!! REARRANGE_CR sorts a sparse compressed row matrix.
!
!  Discussion:
!
!    This routine guarantees that the entries in the CR matrix
!    are properly sorted.
!
!    After the sorting, the entries of the matrix are rearranged in such
!    a way that the entries of each column are listed in ascending order
!    of their column values.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
!
!    Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
!    On output, these may have been rearranged by the sorting.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) i4temp
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) r8temp

  do i = 1, n

    do k = ia(i), ia(i+1) - 2
      do l = k + 1, ia(i+1) - 1

        if ( ja(l) < ja(k) ) then
          i4temp = ja(l)
          ja(l)  = ja(k)
          ja(k)  = i4temp

          r8temp = a(l)
          a(l)   = a(k)
          a(k)   = r8temp
        end if

      end do
    end do

  end do

  return
end subroutine rearrange_cr
!*****************************************************************************80
function ddot ( n, dx, incx, dy, incy )
!*****************************************************************************80
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Dongarra, Moler, Bunch, Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries in X.
!
!    Input, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive entries in Y.
!
!    Output, real DDOT, the sum of the product of the corresponding
!    entries of X and Y.
!
  implicit none

  real ( kind = 8 ) ddot
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  ddot = 0.0D+00
  dtemp = 0.0D+00

  if ( n <= 0 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dtemp = dtemp + dx(ix) * dy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 5 )

    do i = 1, m
      dtemp = dtemp + dx(i) * dy(i)
    end do

    do i = m+1, n, 5

      dtemp = dtemp + dx(i  ) * dy(i  ) &
                    + dx(i+1) * dy(i+1) &
                    + dx(i+2) * dy(i+2) &
                    + dx(i+3) * dy(i+3) &
                    + dx(i+4) * dy(i+4)
    end do

  end if

  ddot = dtemp

  return
end function ddot
end module solver
