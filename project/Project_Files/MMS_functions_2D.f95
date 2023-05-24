!=============================================================================80
module set_precision

  implicit none

  save

  integer, parameter :: sngl = selected_real_kind( 6,  37)
  integer, parameter :: dbl  = selected_real_kind(15, 307)
  integer, parameter :: dp = dbl

end module
!=============================================================================80
module constants

  use set_precision, only : dp

  implicit none

  real(dp), parameter :: one   = 1.0_dp
  real(dp), parameter :: two   = 2.0_dp
  real(dp), parameter :: three = 3.0_dp
  real(dp), parameter :: four  = 4.0_dp
  real(dp), parameter :: five  = 5.0_dp
  real(dp), parameter :: six   = 6.0_dp
  
end module
!=============================================================================80
module set_inputs

  use set_precision, only : dp

  private

  public :: pi
  public :: initialize_constants
  
  real(dp) :: pi

  contains
  
  subroutine initialize_constants

    use constants,  only : one

    implicit none

    pi = acos(-one)
  
  end subroutine initialize_constants 

end module
!=============================================================================80
module mms_constants

  use set_precision, only : dp

  implicit none

 ! NOTE: These are currently set up to run the supersonic manufactured solution
  real(dp), parameter :: rho0   = 1.0_dp
  real(dp), parameter :: rhox   = 0.15_dp
  real(dp), parameter :: rhoy   = -0.1_dp
  real(dp), parameter :: uvel0  = 800.0_dp
  real(dp), parameter :: uvelx  = 50.0_dp
  real(dp), parameter :: uvely  = -30.0_dp
  real(dp), parameter :: vvel0  = 800.0_dp
  real(dp), parameter :: vvelx  = -75.0_dp
  real(dp), parameter :: vvely  = 40.0_dp
  real(dp), parameter :: wvel0  = 0.0_dp
  real(dp), parameter :: wvelx  = 0.0_dp
  real(dp), parameter :: wvely  = 0.0_dp
  real(dp), parameter :: press0 = 100000.0_dp
  real(dp), parameter :: pressx = 20000.0_dp
  real(dp), parameter :: pressy = 50000.0_dp
  
end module mms_constants
!=============================================================================80
!!!!MMS functions
pure function rho_mms(length, x, y)

  use set_precision, only : dp
  use constants,     only : two
  use mms_constants, only : rho0, rhox, rhoy
  use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(in) :: y
  real(dp), intent(in) :: length
  real(dp) :: rho_mms

  rho_mms = rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length)

end function rho_mms
!=============================================================================80
pure function uvel_mms(length, x, y)

  use set_precision, only : dp
  use constants,     only : two, three, five
  use mms_constants, only : uvel0, uvelx, uvely
  use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(in) :: y
  real(dp), intent(in) :: length
  real(dp) :: uvel_mms
  
  uvel_mms = uvel0 + uvely*cos((three*pi*y)/(five*length)) +                   &
                 uvelx*sin((three*pi*x)/(two*length))

end function uvel_mms
!=============================================================================80
pure function vvel_mms(length,x,y)

  use set_precision, only : dp
  use constants,     only : two, three
  use mms_constants, only : vvel0, vvelx, vvely
  use set_inputs,    only : pi

  implicit none

  real(dp), intent(in)  :: x
  real(dp), intent(in)  :: y
  real(dp), intent(in)  :: length
  real(dp) :: vvel_mms

  vvel_mms = vvel0 + vvelx*cos((pi*x)/(two*length)) +                          &
                 vvely*sin((two*pi*y)/(three*length))
end function vvel_mms
!=============================================================================80
pure function press_mms(length,x,y)

  use set_precision, only : dp
  use constants,     only : two
  use mms_constants, only : press0, pressx, pressy
  use set_inputs,    only : pi
  
  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(in) :: y
  real(dp), intent(in) :: length
  real(dp) :: press_mms

  press_mms = press0 + pressx*cos((two*pi*x)/length) + pressy*sin((pi*y)/length)

end function press_mms
!=============================================================================80
pure function rmassconv(length,x,y)

  use set_precision, only : dp
  use constants,     only : two, three, five
  use mms_constants, only : rho0, rhox, rhoy, uvel0, uvelx, uvely,             &
                            vvel0, vvelx, vvely, press0, pressx, pressy
  use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(in) :: y
  real(dp), intent(in) :: length
  real(dp) :: rmassconv

  rmassconv = (three*pi*uvelx*cos((three*pi*x)/(two*length)) *                 &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))) /        &
    (two*length) + (two*pi*vvely*cos((two*pi*y)/(three*length)) *              &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))) /        &
    (three*length) + (pi*rhox*cos((pi*x)/length) *                             &
    (uvel0 + uvely*cos((three*pi*y)/(five*length)) + uvelx*sin((three*pi*x)/   &
    (two*length))))/length - (pi*rhoy*sin((pi*y)/(two*length)) *               &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) + vvely*sin((two*pi*y) /           &
    (three*length))))/(two*length)

end function rmassconv
!=============================================================================80
pure function xmtmconv(length,x,y)

  use set_precision, only : dp
  use constants,     only : two, three, five
  use mms_constants, only : rho0, rhox, rhoy, uvel0, uvelx, uvely,             &
                            vvel0, vvelx, vvely, press0, pressx, pressy
  use set_inputs,    only : pi
  
  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(in) :: y
  real(dp), intent(in) :: length
  real(dp) :: xmtmconv

  xmtmconv = (three*pi*uvelx*cos((three*pi*x)/(two*length)) *                  &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length)) *         &
    (uvel0 + uvely*cos((three*pi*y)/(five*length)) +                           &
    uvelx*sin((three*pi*x)/(two*length))))/length +                            &
    (two*pi*vvely*cos((two*pi*y) /                                             &
    (three*length))*(rho0 + rhoy*cos((pi*y)/(two*length)) +                    &
    rhox*sin((pi*x)/length))*(uvel0 + uvely*cos((three*pi*y) /                 &
    (five*length)) + uvelx*sin((three*pi*x)/(two*length))))/(three*length) +   &
    (pi*rhox*cos((pi*x)/length)*(uvel0 + uvely*cos((three*pi*y) /              &
    (five*length)) + uvelx*sin((three*pi*x)/(two*length)))**2)/length -        &
    (two*pi*pressx*sin((two*pi*x)/length))/length -                            &
    (pi*rhoy*(uvel0 + uvely*cos((three*pi*y)/(five*length)) +                  &
    uvelx*sin((three*pi*x)/(two*length)))*sin((pi*y)/(two*length))*            &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length))))/(two*length) -                      &
    (three*pi*uvely*(rho0 + rhoy*cos((pi*y)/(two*length)) +                    &
    rhox*sin((pi*x)/length))*sin((three*pi*y)/(five*length))*(vvel0 + vvelx *  &
    cos((pi*x)/(two*length)) + vvely*sin((two*pi*y)/(three*length)))) /        &
    (five*length)

end function xmtmconv
!=============================================================================80
pure function ymtmconv(length,x,y)

  use set_precision, only : dp
  use constants,     only : two, three, four, five
  use mms_constants, only : rho0, rhox, rhoy, uvel0, uvelx, uvely,             &
                            vvel0, vvelx, vvely, press0, pressx, pressy
  use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(in) :: y
  real(dp), intent(in) :: length
  real(dp) :: ymtmconv

  ymtmconv = (pi*pressy*cos((pi*y)/length))/length -                           &
    (pi*vvelx*sin((pi*x)/(two*length))*(rho0 + rhoy*cos((pi*y)/(two*length)) + &
    rhox*sin((pi*x)/length))*(uvel0 + uvely*cos((three*pi*y)/(five*length)) +  &
    uvelx*sin((three*pi*x)/(two*length))))/(two*length) +                      &
    (three*pi*uvelx*cos((three*pi*x)/(two*length)) *                           &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length)) *         &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length))))/(two*length) +                      &
    (four*pi*vvely*cos((two*pi*y) /                                            &
    (three*length))*(rho0 + rhoy*cos((pi*y)/(two*length)) +                    &
    rhox*sin((pi*x)/length))*(vvel0 + vvelx*cos((pi*x)/(two*length)) +         &
    vvely*sin((two*pi*y)/(three*length))))/(three*length) +                    &
    (pi*rhox*cos((pi*x)/length) *                                              &
    (uvel0 + uvely*cos((three*pi*y)/(five*length)) +                           &
    uvelx*sin((three*pi*x)/(two*length))) *                                    &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length))))/length -                            &
    (pi*rhoy*sin((pi*y)/(two*length)) *                                        &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length)))**2)/(two*length)

end function ymtmconv
!=============================================================================80
pure function energyconv(gamma, length,x,y)

  use set_precision, only : dp
  use constants,     only : one, two, three, four, five, six
  use mms_constants, only : rho0, rhox, rhoy, uvel0, uvelx, uvely,             &
                            vvel0, vvelx, vvely, wvel0, wvelx, wvely,          &
                            press0, pressx, pressy
  use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: gamma
  real(dp), intent(in) :: length
  real(dp), intent(in) :: x
  real(dp), intent(in) :: y
  real(dp) :: energyconv

  energyconv = (uvel0 + uvely*cos((three*pi*y)/(five*length)) +                &
    uvelx*sin((three*pi*x)/(two*length)))*((-two*pi*pressx*sin((two*pi*x) /    &
    length))/length + (rho0 + rhoy*cos((pi*y)/(two*length)) +                  &
    rhox*sin((pi*x)/length))*((-two*pi*pressx*sin((two*pi*x)/length))/         &
    ((-one + gamma)*length*(rho0 + rhoy*cos((pi*y)/(two*length)) +             &
    rhox*sin((pi*x)/length))) + ((three*pi*uvelx*cos((three*pi*x) /            &
    (two*length))*(uvel0 + uvely*cos((three*pi*y)/(five*length)) +             &
    uvelx*sin((three*pi*x)/(two*length))))/length - (pi*vvelx*sin((pi*x) /     &
    (two*length))*(vvel0 + vvelx*cos((pi*x)/(two*length)) +                    &
    vvely*sin((two*pi*y)/(three*length))))/length)/two - (pi*rhox*cos((pi*x) / &
    length)*(press0 + pressx*cos((two*pi*x)/length) +                          & 
    pressy*sin((pi*y)/length)))/((-one + gamma)*length*(rho0 + rhoy*cos((pi*y)/&
    (two*length)) + rhox*sin((pi*x)/length))**2)) +                            &
    (pi*rhox*cos((pi*x)/length)*((wvel0**2 + (uvel0 + uvely*cos((three*pi*y) / &
    (five*length)) + uvelx*sin((three*pi*x)/(two*length)))**2 +                &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) + vvely*sin((two*pi*y) /           &
    (three*length)))**2)/two + (press0 + pressx*cos((two*pi*x)/length) +       &
    pressy*sin((pi*y)/length))/((-one + gamma) *                               &
    (rho0 + rhoy*cos((pi*y)/(two*length)) +                                    &
    rhox*sin((pi*x)/length)))))/length) +                                      &
    (three*pi*uvelx*cos((three*pi*x)/(two*length)) *                           &
    (press0 + pressx*cos((two*pi*x)/length) + pressy*sin((pi*y)/length) +      &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))*          &
    ((wvel0**2 + (uvel0 + uvely*cos((three*pi*y)/(five*length)) +              &
    uvelx*sin((three*pi*x)/(two*length)))**2 +                                 &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length)))**2)/two +                            &
    (press0 + pressx*cos((two*pi*x)/length) +                                  &
    pressy*sin((pi*y)/length))/((-one + gamma) *                               &
    (rho0 + rhoy*cos((pi*y)/(two*length)) +                                    &
    rhox*sin((pi*x)/length))))))/(two*length) +                                &
    (two*pi*vvely*cos((two*pi*y)/(three*length)) *                             &
    (press0 + pressx*cos((two*pi*x)/length) +                                  &
    pressy*sin((pi*y)/length) + (rho0 + rhoy*cos((pi*y)/(two*length)) +        &
    rhox*sin((pi*x)/length))*((wvel0**2 +                                      &
    (uvel0 + uvely*cos((three*pi*y)/(five*length)) +                           &
    uvelx*sin((three*pi*x)/(two*length)))**2 +                                 &
    (vvel0 + vvelx*cos((pi*x)/(two*length)) +                                  &
    vvely*sin((two*pi*y)/(three*length)))**2)/two +                            &
    (press0 + pressx*cos((two*pi*x)/length) + pressy*sin((pi*y)/length)) /     &
	((-one + gamma)*(rho0 + rhoy*cos((pi*y)/(two*length)) +                    &
    rhox*sin((pi*x)/length))))))/(three*length) + (vvel0 + vvelx*cos((pi*x) /  &
    (two*length)) + vvely*sin((two*pi*y)/(three*length))) *                    &
    ((pi*pressy*cos((pi*y)/length))/length - (pi*rhoy*sin((pi*y)/(two*length))*&
    ((wvel0**2 + (uvel0 + uvely*cos((three*pi*y)/(five*length)) +              &
    uvelx*sin((three*pi*x)/(two*length)))**2 + (vvel0 + vvelx *                &
    cos((pi*x)/(two*length)) + vvely*sin((two*pi*y)/(three*length)))**2)/two + &
    (press0 + pressx*cos((two*pi*x)/length) +                                  &
    pressy*sin((pi*y)/length))/((-one + gamma) *                               &
    (rho0 + rhoy*cos((pi*y)/(two*length)) +                                    &
    rhox*sin((pi*x)/length)))))/(two*length) +                                 &
    (rho0 + rhoy*cos((pi*y)/(two*length)) +                                    &
    rhox*sin((pi*x)/length))*((pi*pressy*cos((pi*y)/length)) /                 &
    ((-one + gamma)*length*(rho0 + rhoy*cos((pi*y)/(two*length)) +             &
    rhox*sin((pi*x)/length))) +                                                &
    ((-six*pi*uvely*(uvel0 + uvely*cos((three*pi*y) /                          &
    (five*length)) + uvelx*sin((three*pi*x)/(two*length))) *                   &
    sin((three*pi*y)/(five*length)))/(five*length) +                           &
    (four*pi*vvely*cos((two*pi*y) /                                            &
    (three*length))*(vvel0 + vvelx*cos((pi*x)/(two*length)) +                  &
    vvely*sin((two*pi*y)/(three*length))))/(three*length))/two +               &
    (pi*rhoy*sin((pi*y)/(two*length))*(press0 + pressx*cos((two*pi*x)/length) +&
    pressy*sin((pi*y)/length)))/(two*(-one + gamma)*length*                    &
    (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))**2)))

end function energyconv