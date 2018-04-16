program hts2d
! Adding modules 
use solver
use parmesh
use chars
use physics
use parsolvers
use constants
implicit none
!
! HTS2D is a Fortran 95/03 implementation, of GMRES method to solve a linear
!system that comes from the discretization of a physical domain, under heat 
!transfer problems. This program is modular so feel free to mess around and 
!create another subroutines. I'll look forward for it.
!
!       Programmer: Paulo Conci                Date: 23/03/16
!                                       Last modify: 29/03/18
!
! Declaring all variables used in this file(program), and some used in subrou-
!tines.
!
! Counting variables
   integer*4 :: stat, lin, cls, i, j, nz_num, delt
   integer*4 :: nodes
!
! Time variables
    real(kind=8) :: time_phsc, tm
    integer(kind=8) :: cnt1, cnt2, rate_cnt, cnt_max, nb_ticks
    real(kind=8), allocatable, dimension(:) :: time
!
! Mesh variables
   real*8, allocatable, dimension(:,:) :: msh
   real*8, allocatable, dimension(:,:) :: msh_copy
   real*8 :: n, dxdy
!
! Iterative Method Variables
   character(len=30) :: fout
   character(len=60) :: dirname, dname, dir
   integer(kind=4) :: iwk, ierr, im, getcwd
   real*8 ::  normt
   real*8, allocatable, dimension(:,:) :: alin
   integer(kind=4), allocatable, dimension(:,:) :: conn
   integer(kind=4), allocatable, dimension(:,:) :: connfth
   real*8, allocatable, dimension(:) :: rhs
   real*8, allocatable, dimension(:) :: x
   real*8, allocatable, dimension(:) :: x0
   real*8, allocatable, dimension(:) :: x_exact
   integer*4, allocatable, dimension(:) :: ia
   integer*4, allocatable, dimension(:) :: ja
   real*8, allocatable, dimension(:) :: a
   real*8, allocatable, dimension(:) :: t
   real(kind=8), allocatable, dimension(:) :: ao
   integer(kind=4), allocatable, dimension(:) :: iao
   integer(kind=4), allocatable, dimension(:) :: jao
   integer(kind=4), allocatable, dimension(:) :: ju
   real(kind=8), allocatable, dimension(:) :: wu
   real(kind=8), allocatable, dimension(:) :: wl
   real(kind=8), allocatable, dimension(:) :: alu
   real(kind=8), allocatable, dimension(:,:) :: vv
   real(kind=8), allocatable, dimension(:) ::w
   integer(kind=4), allocatable, dimension(:) :: jr
   integer(kind=4), allocatable, dimension(:) :: jwl
   integer(kind=4), allocatable, dimension(:) :: jwu
   integer(kind=4), allocatable, dimension(:) :: jlu
   integer(kind=4), allocatable, dimension(:) :: iw
   integer(kind=4), allocatable, dimension(:) :: levs
   integer(kind=4), allocatable, dimension(:) :: jw
!
! Defining some name, for lazyness at the time, this was the better solution,
! it will be appreciated ideas to automated this.
    if(choice .eq. 1)then
      write(*,*) 'Solving w/ Gauss Siedel'
      fout = 'residual_gauss_siedel.dat'
      dname = 'gauss'
    elseif(choice .eq. 2)then
      write(*,*) 'Solving w/ GMRES'
      fout = 'residual_gmres.dat'
      dname = 'gmres'
    elseif(choice .eq. 3)then
      write(*,*) 'Solving w/ GMRES p. ILU0'
      fout = 'residual_ilu0.dat'
      dname = 'ilu0'
    elseif(choice .eq. 4)then
      write(*,*) 'Solving w/ GMRES p. ILUT'
      fout = 'residual_ilut.dat'
      dname = 'ilut'
    elseif(choice .eq. 5)then
      write(*,*) 'Solving w/ GMRES p. MILU0'
      fout = 'residual_milu0.dat'
      dname = 'milu'
    elseif(choice .eq. 6)then
      write(*,3) 'Solving w/ GMRES p. ILU', lfil
        if(lfil .eq. 0)then
        dname = 'ilu0'
        fout = 'residual_ilu0.dat'
        elseif(lfil .eq. 1)then
        dname = 'ilu1'
        fout = 'residual_ilu1.dat'
        elseif(lfil .eq. 2)then
        dname = 'ilu2'
        fout = 'residual_ilu2.dat'
        elseif(lfil .eq. 3)then
        dname = 'ilu3'
        fout = 'residual_ilu3.dat'
        elseif(lfil .eq. 4)then
        dname = 'ilu4'
        fout = 'residual_ilu4.dat'
        else
        dname = 'iluk'
        fout = 'residual_iluk.dat'
        endif
    endif
3 format(a, I1)
    if(choice2 .eq. 1)then
      write(*,*) 'Unsteady State'
    elseif((choice2 .eq. 2) .or. (choice2 .eq. 3))then
      write(*,*) 'Steady State'
    endif
!
! Creating directories and open files to read the mesh file and write out the 
!solution. This is literally where everything starts.
    call execute_command_line('mkdir -p ' // adjustl(trim(dname)))
!    call execute_command_line('pwd ')
    stat = getcwd(dir)
    if(stat .ne. 0 )stop 'getcwd: error'
!
    dirname = adjustl(adjustr(trim(trim(dir)//'/'//trim(dname)//'/')))
    dirname = trim(adjustl(trim(dirname)))
!
    open(unit=2, file=mshdir//meshfile, status= 'old', iostat = stat)
    open(unit=43, file= trim(adjustl(trim(dirname))//'residual_timemachine.plt'), status='unknown')
    open(unit=44, file= trim(adjustl(trim(dirname))//'residual_iteration.plt'), status='unknown')
    open(unit=45, file= trim(adjustl(trim(dirname))//'residual_timephysics.plt'), status='unknown')
    open(unit=80, file= trim(adjustl(trim(dirname))//savename), status='unknown')
!
! Counting how many number are on the mesh file
!for allocate after (this will use only the necessary memory)
  if(stat .eq. 0) then
    nodes=0
    do
      read (2, *, iostat = stat) n
      if(stat .ne. 0)exit
      nodes = nodes + 1
    end do
  endif
!
!Allocating the memory for mesh, linear system and solution.
   allocate (msh(nodes,4))
   allocate (msh_copy(nodes,4))
   allocate (alin(nodes,nodes))
   allocate (rhs(nodes))
   allocate (t(nodes))
   allocate (x(nodes))
   allocate (x0(nodes))
   allocate (time(itr_max))
   allocate (conn(nodes,9))
   ! allocate (connfth(nodes,9))
!
!This will reopen unit two for avoiding reading errors.
  rewind (unit=2)
!
!Reading matrix of nodes
  read(2,*, iostat=stat) ((msh(lin,cls), cls=1,4), lin=1,nodes)
  if (stat .eq. 0)then
! This print the mesh file on screen, uncomment for debug
!	write(*,300) ((msh(lin,cls), cls=1,4), lin=1,nodes)
!300 format(4(F20.6))
! This is VERY important for the program to work, dxdy should always be the 
! distance between two mesh points, if your mesh file is not ordered in the 
! correct way this value will be wrong and crash the solution.
    dxdy = msh(2,2)
  else
    write(*,*)'Error reading matrix', stat
  end if
  do i=1, nodes
!	    msh(i,1) = i
    msh(i,4) = 1.0d0
  enddo
!
  call connect(msh, conn, nodes, dxdy)
  ! call connect4th(msh,connfth,nodes,dxdy)
  write(*,13) ((conn(i,j),j=1,9),i=1,nodes)
!
!
13 format(9(I4))
! proportion bb
!    do i = 1, nodes
!        msh(i,2) = 1.d3 * msh(i,2)
!        msh(i,3) = 1.d3 * msh(i,3)
!    enddo
!
!Defining initial parameter for denominator of norm temperature
    x0(1:nodes) = msh(1:nodes,4)
!Setting the boundary conditions, msh is an intent(inout)
    call bound_cond(msh, nodes)
!Iniciating the allocatable matrix for avoiding errors
    msh_copy = msh
    x(1:nodes) = msh(1:nodes,4)
!
! write(*,*) "HEEEEELLLLLOOOOO MATTY D"
!Here is the critical part, it's where we read all nodes of mesh and
!them generate one equation for each node, that is a finite difference
!equation, we are allocating in alin(i,j) only the the constants terms
!that multiplicate the temperature value, this alin matrix have the
!dimension of (nodes,nodes) and must be square matrix, here we have a
!little problem because we must identifying the neibourhood nodes,
!and for that it is necessary to identifying them.
! Initializind the time variables
    delt = 0
    cnt1 = 0
    cnt2 = 0
    cnt_max = 1000000000
    rate_cnt = 1000
    call system_clock(count_rate=rate_cnt, count_max=cnt_max)
!
!    call cpu_time(start_time)
    call system_clock(count=cnt1)
!    call cpu_time(init_time)
! LOOOOOOOOOOOOOOOOOOOP
5 continue
!
!getting out of the loop on steady state
    if((delt .ne. 0) .and. (choice2 .eq. 2))then
        goto 25
    endif
!
! Saving a copy of the initial values of solution to mensurate how the method
!is behaving
    t(1:nodes) = x(1:nodes)
!
! Starting matrix of constants, at each time step for security
   alin(1:nodes,1:nodes) = 0.0d0
   rhs(1:nodes) = 0.0d0
!
! Subroutine responsible to generate constants, this is still an odd way to do
!it, take and example, if you want to simulate diferent boundary conditions you
!must change this subroutine. At this time is working, but it'll be implemented
!in another more user fliendly way.
    if (choice2 .eq. 1)then
      call create_a_unsteady(msh, t, nodes, alin, rhs, dxdy)
    elseif(choice2 .eq. 2)then
      call ccall(conn, nodes, alin, rhs, t)
    ! write(80,12) ((alin(i,j),j=1,nodes),i=1,nodes)
    ! 12 format(25(F10.5))
     ! call create_a_steady(msh, t, nodes, alin, rhs, dxdy)
    elseif(choice2 .eq. 3)then
      call create_4_steady(msh, t, nodes, alin, rhs, dxdy)
    end if
!Iterative methods
!==============================================================================
!solve the linear system with gauss-siedel iterative method
   if (choice .eq. 1)then
!
      call kicks(nodes, x)
!
      call gs(alin, rhs, nodes, x, fout, dirname)
      if(allocated(alin)) deallocate(alin)
!******************************************************************************
! GMRES without preconditioning with no restart, is not working at the time
   elseif (choice .eq. 2)then
!
!      call pre_gmres(alin, nodes, a, ia, ja, nz_num)
!      if(allocated(alin)) deallocate(alin)
!!
!   ind = (nodes+2) - (nz_num+1)
!!
!       if (delt .eq. 0)then
!        allocate(ao(nz_num))
!        allocate(jao(nz_num))
!        allocate(iao(nodes+1))
!      endif
!!
!      call sttocsr ( nodes, nz_num, a, ia, ja, ao, jao, iao )
!      call kicks(nodes, x)
!!
!       if(delt .eq. 0)then
!        iwk = nodes*nodes
!        allocate(alu(nodes*nodes))
!        allocate(jlu(nodes*nodes))
!        allocate(ju(nodes))
!        allocate(jr(nodes+1))
!        allocate(wl(nodes))
!        allocate(jw(nodes*nodes))
!    endif
!
!      call csrmsr ( nodes, ao, jao, iao, alu, jlu, ju, wl, jr )
!
!      if (ierr .eq. 0)then
!            im = 100
!!
!            if(delt .eq. 0)then
!                allocate(vv(nodes,im+1))
!            endif
!!
!!
!      call pgmres(nodes, im, rhs, x, vv, ao, jao, iao, alu, jlu, ju, ierr )
!!
!      if(ierr .ne. 0)then
!         write(*,*) 'Error on GMRES = ', ierr
!      elseif(ierr .eq. 0)then
!         write(*,*) 'Sucessful Solving!'
!      endif
!      endif
!
! Instead we are using a restarted version of GMRES, is no preconditioned, so
!work as test case for the preconditioned version, although, no all orthogonal
!vector are stored. This saves memory, but not improve convergence.
      call pre_gmres(alin, nodes, a, ia, ja, nz_num)
      if(allocated(alin)) deallocate(alin)
!
      call kicks(nodes, x)
!
      call mgmres_st ( nodes, nz_num, ia, ja, a, x, rhs, fout, dirname)
!******************************************************************************
! GMRES preconditioned with ILU0
   elseif(choice .eq. 3) then
!
      call pre_gmres(alin, nodes, a, ia, ja, nz_num)
      if(allocated(alin)) deallocate(alin)
!
    if (delt .eq. 0)then
        allocate(ao(nz_num))
        allocate(jao(nz_num))
        allocate(iao(nodes+1))
    endif
      call sttocsr ( nodes, nz_num, a, ia, ja, ao, jao, iao )
!
   if(delt .eq. 0)then
      allocate(iw(nodes))
      allocate(alu(nodes*nodes))
      allocate(jlu(nodes*nodes))
      allocate(ju(nodes))
   endif
      call ilu0 ( nodes, ao, jao, iao, alu, jlu, ju, iw, ierr )
   if(ierr .ne. 0)then
      write(*,*) 'Error after ILU0 = ', ierr
   endif
!
      call kicks(nodes, x)
!
!
      if (ierr .eq. 0)then
            im = 500
!
            if(delt .eq. 0)then
                allocate(vv(nodes,im+1))
            endif
!
      call pgmres(nodes, im, rhs, x, vv, ao, jao, iao, alu, jlu, ju, fout, ierr, dirname)
!
      if(ierr .ne. 0)then
         write(*,*) 'Error on GMRES ILU0 = ', ierr
      elseif(ierr .eq. 0)then
         write(*,*) 'Sucessful Solving!'
      endif
      endif
!******************************************************************************
!GMRES with ILUT
    elseif(choice .eq. 4)then
!
!storing alin in sparse triplet form
        call pre_gmres(alin, nodes, a, ia, ja, nz_num)
        if(allocated(alin)) deallocate(alin)
!
!storing a(nz_num), ja(nz_num) and ia(nz_num) ia CSR
    if (delt .eq. 0)then
        allocate(ao(nz_num))
        allocate(jao(nz_num))
        allocate(iao(nodes+1))
    endif
        call sttocsr ( nodes, nz_num, a, ia, ja, ao, jao, iao )
!
!setting some parameters of the preconditioner
    if(delt .eq. 0)then
        iwk = nodes*nodes
        allocate(alu(nodes*nodes))
        allocate(jlu(nodes*nodes))
        allocate(ju(nodes))
        allocate(wu(nodes+1))
        allocate(wl(nodes))
        allocate(jr(nodes))
        allocate(jwl(nodes))
        allocate(jwu(nodes))
    endif
!calling the preconditioner ILUT
        call ilut ( nodes, ao, jao, iao, alu, jlu, ju, &
                    iwk, wu, wl, jr, jwl, jwu, ierr )
!
    if(ierr .ne. 0)then
      write(*,*) 'Error after ILUT = ', ierr
    endif
        call kicks(nodes, x)
!
        if (ierr .eq. 0)then
            im = 500
!
            if(delt .eq. 0)then
                allocate(vv(nodes,im+1))
            endif
!
      call pgmres(nodes, im, rhs, x, vv, ao, jao, iao, alu, jlu, ju, fout, ierr, dirname)
!
      if(ierr .ne. 0)then
         write(*,*) 'Error on GMRES w/ ILUT = ', ierr
      elseif(ierr .eq. 0)then
         write(*,*) 'Sucessful Solving!'
      endif
   endif
!******************************************************************************
! GMRES with MILU0
    else if(choice .eq. 5)then
        call pre_gmres(alin, nodes, a, ia, ja, nz_num)
        if(allocated(alin)) deallocate(alin)
!
        if (delt .eq. 0)then
            allocate(ao(nz_num))
            allocate(jao(nz_num))
            allocate(iao(nodes+1))
        endif
            call sttocsr ( nodes, nz_num, a, ia, ja, ao, jao, iao )
!
        if(delt .eq. 0)then
            allocate(alu(nz_num+2))
            allocate(jlu(nz_num+2))
            allocate(ju(nodes))
            allocate(iw(nodes))
        endif
            call milu0( nodes, ao, jao, iao, alu, jlu, ju, iw, ierr )

!
   if (ierr .eq. 0)then
                im = 500
!
        call kicks(nodes, x)
!
                if(delt .eq. 0)then
                    allocate(vv(nodes,im+1))
                endif
!
      call pgmres(nodes, im, rhs, x, vv, ao, jao, iao, alu, jlu, ju, fout, ierr, dirname)
!
      if(ierr .ne. 0)then
         write(*,*) 'Error on GMRES w/ ILU0 = ', ierr
      elseif(ierr .eq. 0)then
         write(*,*) 'Sucessful Solving!'
      endif
   endif
!******************************************************************************
!Gmres w/ iluk
        elseif(choice .eq. 6)then
!
!storing alin in sparse triplet form
        call pre_gmres(alin, nodes, a, ia, ja, nz_num)
        ! if(allocated(alin)) deallocate(alin)
!
!storing a(nz_num), ja(nz_num) and ia(nz_num) ia CSR
    if (delt .eq. 0)then
        allocate(ao(nz_num))
        allocate(jao(nz_num))
        allocate(iao(nodes+1))
    endif
        call sttocsr ( nodes, nz_num, a, ia, ja, ao, jao, iao )
!
!setting some parameters of the preconditioner
    if(delt .eq. 0)then
        iwk = nodes*nodes
        allocate(alu(iwk))
        allocate(jlu(iwk))
        allocate(ju(nodes))
        allocate(levs(iwk))
        allocate(w(nodes))
        allocate(jw(3*nodes))
    endif
!calling the preconditioner ILUk
!
      call iluk(nodes, ao, jao, iao, alu, jlu, ju, levs, iwk, w, jw, ierr)
!
   if(ierr .ne. 0)then
      write(*,*) 'Error on ILU',lfil,' = ', ierr
   elseif(ierr .eq. 0)then
      write(*,*) 'Sucessful precondiotioning'
   endif
!
   call kicks(nodes, x)
!
   if (ierr .eq. 0)then
            im = 500
!
      if(delt .eq. 0)then
          allocate(vv(nodes,im+1))
      endif
!
   call pgmres(nodes, im, rhs, x, vv, ao, jao, iao, alu, jlu, ju, fout, ierr, dirname)
      if(ierr .ne. 0)then
         write(*,*) 'Error on GMRES w/ ILUK = ', ierr
      elseif(ierr .eq. 0)then
         write(*,*) 'Sucessful Solving'
      endif
   endif
   endif
!******************************************************************************
!Checking the convergence, this is made with a norm L2, with both initial and
!final solution. If the norm is less the *epsolon* (yes I know this is not the
!real word) the program exits the loop and is take the last vector as solution
!for the problem.
   call cvgc(x, x0, t, nodes, normt)
!
   if((normt .lt. epsolon) .and. (delt .ge. itr_max))then
      goto 25
   else if((normt .gt. epsolon) .and. (delt .lt. itr_max))then
!
! Counting how many interations has taken to obtain the solution.
  delt = delt + 1
!
! Mensuring the time your CPU took to obtain each iteration. This is still in
!an early state. But the main goal is to compare how the method takes less time
!to obtain a solution once it get close to the analytical solution.
        call system_clock(count=cnt2)
        nb_ticks = cnt2 - cnt1
        if(cnt2 .lt. cnt1)then
            nb_ticks = nb_ticks + cnt_max
        endif
            time(delt) = dfloat(nb_ticks)/rate_cnt
            time_phsc = delt*deltat
!
        if(delt .eq. 1)then
            tm = time(delt)
        else
            tm = time(delt) - time(delt - 1)
        endif
!
! Writing to files the number of the iteration, and time taken to obtain, the 
!norm value and the physical CPU time.
    write(43,*) delt, tm
    write(44,*) delt, normt
    write(45,*) time_phsc, normt
!writing a backup file for power loss or program interruption. Early stages, 
!this needs a subroutine to read the file and a parameter to tell if it's a 
!new simulation or a old one.
! do i=1, nodes
!     write(80,*) i, x(i)
! enddo
!
! I use TecPlot to plot the graphs, so the value file needs to be formated 
!before ploting, to use another ploting program probably you will need 
!another format, so feel free to implemente your own subroutine and call it.
        msh(1:nodes,4) = x(1:nodes)
        call temp_time_plot(msh, nodes, dxdy)
        call temp_origin(msh, nodes, delt, dxdy)
!
! Tell the user that and iteration has occur.
        write (*,*) 'Outer Iteration =', delt
!
! Once we deallocate alin inside the iteration process to save memory ( meshs
!with 30k nodes, storing real*8 variables, will need up to 6GB of RAM), we need
!to reallocate that space.
! write(*,*) "heellloooooo"
        ! allocate(alin(nodes, nodes))
!
    goto 5
   end if
!
!END OF LOOOOOOOOOOOOOOOP
25 continue
!
! Verifying how much time has taken to solve the problem.
    call system_clock(cnt2, rate_cnt, cnt_max)
    write(*,*) 'System Clock = ', real(cnt2 - cnt1) / real(rate_cnt)
!    call cpu_time(finish_time)
!    write(*,*) 'Total time is = ', finish_time - init_time
    msh(1:nodes,4) = x(1:nodes)
!
!    delt = delt + 1
!    write(*,*) 'Outer Iteration =', delt
!
!Printing the values of solution in a file for opening with tecplot
   call tcplt(msh, nodes, filename, dxdy, dirname)
   if(allocated(msh)) deallocate(msh)
!Only for steady cases
!==============================================================================
! Once the problem is simple and has a analytical solution, is possible to 
!verified how close are the approximate solution by the numerical method and 
!the exact solution to that problem.
  if (choice3 .eq. 1)then
!Calculating the exact solution for steady state case
    allocate(x_exact(nodes))
    if(choice2 .eq. 2 .or. choice2 .eq. 3) then
        call temp_exact(msh_copy, nodes)
        x_exact(1:nodes) = msh_copy(1:nodes,4)
!printing the exact solution
        call tcplt(msh_copy, nodes, filename_exact, dxdy, dirname)
        if(allocated(msh_copy)) deallocate(msh_copy)
!Calculate the logarithm error for the steady state
        call logerror(x, x_exact, nodes, dxdy, dirname)
!    write(*,*) 'type the next mesh file name'
!    read (*,*) meshfile
!    goto3
    elseif(choice2 .eq. 1) then
        call reddy(msh_copy, nodes)
        x_exact(1:nodes) = msh_copy(1:nodes,4)
        call tcplt(msh_copy, nodes, filename_exact, dxdy, dirname)
        if(allocated(msh_copy)) deallocate(msh_copy)
!Calculate the logarithm error for the steady state
        call logerror(x, x_exact, nodes, dxdy, dirname)
    endif
!
  endif
stop
!I guess this is THE END.
end program hts2d
