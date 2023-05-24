! This is a more detailed example of how to write a output file for use with Tecplot.
! This will use Tecplot's ability to plot using cell centered data instead of having to do the interpolation yourself.
! If you want to make a more speciallized version for Tecplot you can look at the data format manual
! at http://download.tecplot.com/360/current/360_data_format_guide.pdf in the ASCII section.

! First we need to set up the header for the file.
! This tells tecplot the title of the data set and what variables you will be giving it.

open(10,file="output.dat")
write(10,*) 'TITLE = "Dataset Title"'
write(10,*) 'variables = "x","y","var1","var2","var3"'

! Next we will output the data. In tecplot this is done in zones
! If you want solution data as you converge you would write one zone for each time step.
! (Please do not output every time step this will take far too long)
! Or if you had say Exact Soln, Numerical Soln, and DE those could be 3 seperate zones.

write(10,*) 'ZONE'
write(10,*) 'T = "Zone Title"'                             ! Maybe this is the iteration number
write(10,*) 'I =', imax, 'J =', jmax                       ! This is the number of nodes in the I and J directions
write(10,*) 'DT = DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE'  ! This is the data type for each variable
                                                           ! (Must tell tecplot it is double presision default is single)
write(10,*) 'DATAPACKING = BLOCK'                          ! This is either Block or Point. Block must be used for cell centered
                                                           ! BLOCK means you give it all of var1 then all of var2
                                                           ! POINT means you give it var1 and var2 at one point at a time
write(10,*) 'VARLOCATION = ([3-4]=CELLCENTERED)            ! This tells tecplot which variables are nodal and which are cell centered
                                                           ! In this example var1 and var2 are cell centered and x, y, and var3 are nodal
write(10,*) 'SOLUTIONTIME =', time                         ! If you want to use time animation in tecplot you need to give it a time for the zone
                                                           ! For simplicity this can just be the iteration number if you want
                                                           ! You do not have to use this option if you dont want to

! Now that our zone is setup we output the data
! write out x
do i=1,imax
  do j=1,jmax
    write(10,*) x(i,j)
  end do
end do

! write out y
do i=1,imax
  do j=1,jmax
    write(10,*) y(i,j)
  end do
end do

! write out var1
do i=1,imax-1             ! note imax-1 since these are cell centered values
  do j=1,jmax-1
    write(10,*) var1(i,j)
  end do
end do

! write out var2
do i=1,imax-1
  do j=1,jmax-1
    write(10,*) var2(i,j)
  end do
end do

! write out var3
do i=1,imax
  do j=1,jmax
    write(10,*) var3(i,j)
  end do
end do


! You can then write out another zone if you want otherwise you are done
close(10)


! Now in Tecplot you can read this file in directly and it will make pretty pictures.
! To read in File > Load Data File > Tecplot Data Loader and choose output.dat
! Then choose 2D Cartesian as plot type.
! Now you should see a blank plot with only the borders shown.
! You can check that you read in and output the grid correctly by checking mesh in the layer list
! You can plot contour plots by checking the contour check box.
! To change what you are plotting you can push the button next to the checkbox
! and choose what variable to plot from the drop down menu in the upper right.
! Finally you can play your animation (if you set this up in the output file) by
! using the play controller below the checkboxes. Make sure you are on your final solution before
! you make your plots!

! Tecplot can also make some of the checking of your code easier with several tools.
! First you can use the primary cell value to see what the actual solution is in the cell
! (as opposed to the interpolated contours).
! To do this press the zone style button, then go to the contour tab and change the contour type
! to primary value flood.
! Next you can use the streamline functionality to check your boundary conditions, especially your
! wall boundaries and your periodic boundaries. Streamlines should go smoothly through the periodic
! boundaries and they should not go through walls (of course). Thus you can plot streamlines to make
! sure this is the case.
! To do this press the button with a picture next to streamlines checkbox. Then click where you
! want to start a stream line or you can click and drag to start several streamlines.

! If you want more information or you can google it, email Chip, or go his office hours.