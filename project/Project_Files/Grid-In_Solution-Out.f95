
!==========================SAMPLE GRID INPUT==========================
  ! Note: this Fortran 95 code will read in the grid as follows:
  ! Read all x coordinates, looping over i first, then j
  !     => it will have to loop k = 1 to 2 if kmax_grid = 2
  ! Read all y coordinates, looping over i first, then j
  !     => it will have to loop k = 1 to 2 if kmax_grid = 2
  ! The z coordinates are read into the zztemp variable, but not used

  open(unit=12,file='grid.dat')
  read(12,*) nzones
  read(12,*) imax, jmax, kmax 
  read(12,*) (((x(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
             (((y(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
             (((zztemp,i=1,imax),j=1,jmax),k=1,kmax)

  !======ALTERNATE FORM (for reading 2D Grid from 3D File======
  open(unit=12,file='grid.dat')
  read(12,*) nzones
  read(12,*) imax, jmax, kmax 
  ! Read in x-coordinate
  do k = 1, kmax
  do j = 1, jmax
  do i = 1, imax
    read(12,*) x(i,j)
  enddo
  enddo
  enddo
  ! Read in y-coordinate
  do k = 1, kmax
  do j = 1, jmax
  do i = 1, imax
    read(12,*) y(i,j)
  enddo
  enddo
  enddo


!==========================SAMPLE SOLUTION OUTPUT==========================
!======== This writes out a Tecplot output file in "point" format =========

  ! Note: write this once at the beginning
  open(40,file='2dEuler-curvilinear.dat',status='unknown')
  write(40,*) 'TITLE = "2D Curvilinear Euler Equations Field Data"'
  write(40,*) 'variables="x(m)""y(m)""rho(kg/m^3)""u(m/s)""v(m/s)""press(N/m^2)"'

  ! Write this every time you want to output the solution
  !    => assumes n has been passed as the iteration #
  write(40,*) 'zone T="',n,'" '
  write(40,*) 'I=',imax,' J=',jmax
  write(40,*) 'DATAPACKING=POINT'  ! Tecplot "point" format
  do j = 1, jmax
  do i = 1, imax
    ! Note: V is assumed to the vector of primitive variables V  
    !       = [rho, u, v, p]^T interpolated to the grid nodes
    write(40,*)x(i,j),y(i,j),V(i,j,1),V(i,j,2),V(i,j,3),V(i,j,4)
  enddo
  enddo
