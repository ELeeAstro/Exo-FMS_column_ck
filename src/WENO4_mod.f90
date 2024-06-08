module WENO4_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private
  public :: interpolate_weno4

contains

  function interpolate_weno4(xs, xp, fp, extrapolate) result(result)
    implicit none
    real(dp), intent(in) :: xs(:)
    real(dp), intent(in) :: xp(:)
    real(dp), intent(in) :: fp(:)
    logical, intent(in), optional :: extrapolate
    real(dp) :: result(size(xs))
    integer :: Ngrid, i, idx, prevB
    real(dp) :: eps, B2, B3
    logical :: use_extrapolate

    Ngrid = size(xp)
    eps = 1.0e-6_dp
    prevB = -1
    B2 = 0.0_dp
    B3 = 0.0_dp
    use_extrapolate = .false.
    if (present(extrapolate)) use_extrapolate = extrapolate

    do idx = 1, size(xs)
      if (xs(idx) < xp(1)) then
        if (.not. use_extrapolate) then
          result(idx) = fp(1)
          cycle
        end if
        i = 1  ! quadratic extrapolation from first three points
      elseif (xs(idx) > xp(Ngrid)) then
        if (.not. use_extrapolate) then
          result(idx) = fp(Ngrid)
          cycle
        end if
        i = Ngrid - 2  ! quadratic extrapolation from last three points
      else
        i = binary_search(xp, xs(idx))
      end if

      if (i >= Ngrid - 1) i = Ngrid - 2

      ! set stencil, pad with zeros when not relevant
      result(idx) = compute_weno4(xs(idx), xp, fp, i, Ngrid, eps, B2, B3, prevB)
    end do

  end function interpolate_weno4

  function binary_search(a, x) result(index)
    implicit none
    real(dp), intent(in) :: a(:)
    real(dp), intent(in) :: x
    integer :: index
    integer :: ileft, i, mid, n

    n = size(a)
    ileft = 1
    i = n

    do while (ileft <= i)
      mid = (ileft + i) / 2
      if (x < a(mid)) then
        i = mid - 1
      else
        ileft = mid + 1
      end if
    end do

    index = max(1, i)

  end function binary_search

  function compute_weno4(x, xp, fp, i, Ngrid, eps, B2, B3, prevB) result(y)
    implicit none
    real(dp), intent(in) :: x, xp(:), fp(:), eps
    integer, intent(in) :: i, Ngrid
    real(dp), intent(inout) :: B2, B3
    integer, intent(inout) :: prevB
    real(dp) :: y
    real(dp) :: xim, xi, xip, xipp, yim, yi, yip, yipp
    real(dp) :: q2, q3, gam2, gam3, al2, al3, om2, om3

    xi = xp(i)
    xip = xp(i + 1)
    yi = fp(i)
    yip = fp(i + 1)

    if (i == 1) then
      xim = 0.0_dp
      xipp = xp(i + 2)
      yim = 0.0_dp
      yipp = fp(i + 2)
    elseif (i == Ngrid - 1) then
      xim = xp(i - 1)
      xipp = 0.0_dp
      yim = fp(i - 1)
      yipp = 0.0_dp
    else
      xim = xp(i - 1)
      xipp = xp(i + 2)
      yim = fp(i - 1)
      yipp = fp(i + 2)
    end if

    call weno4_q(x, xim, xi, xip, xipp, yim, yi, yip, yipp, q2, q3)

    if (i == 1) then
      y = q3
    elseif (i == Ngrid - 1) then
      y = q2
    else
      if (i /= prevB) then
        call weno4_B(xim, xi, xip, xipp, yim, yi, yip, yipp, B2, B3)
        prevB = i
      end if
      gam2 = -(x - xipp) / (xipp - xim)
      gam3 = (x - xim) / (xipp - xim)
      al2 = gam2 / (eps + B2)
      al3 = gam3 / (eps + B3)
      om2 = al2 / (al2 + al3)
      om3 = al3 / (al2 + al3)
      y = om2 * q2 + om3 * q3
    end if

  end function compute_weno4

  subroutine weno4_B(xim, xi, xip, xipp, yim, yi, yip, yipp, B2, B3)
    implicit none
    real(dp), intent(in) :: xim, xi, xip, xipp, yim, yi, yip, yipp
    real(dp), intent(out) :: B2, B3
    real(dp) :: him, hi, hip, H
    real(dp) :: yyim, yyi, yyip, yyipp

    him = xi - xim
    hi = xip - xi
    hip = xipp - xip
    H = him + hi + hip

    yyim = -((2.0_dp * him + hi) * H + him * (him + hi)) / (him * (him + hi) * H) * yim
    yyim = yyim + ((him + hi) * H) / (him * hi * (hi + hip)) * yi
    yyim = yyim - (him * H) / ((him + hi) * hi * hip) * yip
    yyim = yyim + (him * (him + hi)) / ((hi + hip) * hip * H) * yipp

    yyi = -(hi * (hi + hip)) / (him * (him + hi) * H) * yim
    yyi = yyi + (hi * (hi + hip) - him * (2.0_dp * hi + hip)) / (him * hi * (hi + hip)) * yi
    yyi = yyi + (him * (hi + hip)) / ((him + hi) * hi * hip) * yip
    yyi = yyi - (him * hi) / ((hi + hip) * hip * H) * yipp

    yyip = (hi * hip) / (him * (him + hi) * H) * yim
    yyip = yyip - (hip * (him + hi)) / (him * hi * (hi + hip)) * yi
    yyip = yyip + ((him + 2.0_dp * hi) * hip - (him + hi) * hi) / ((him + hi) * hi * hip) * yip
    yyip = yyip + ((him + hi) * hi) / ((hi + hip) * hip * H) * yipp

    yyipp = -((hi + hip) * hip) / (him * (him + hi) * H) * yim
    yyipp = yyipp + (hip * H) / (him * hi * (hi + hip)) * yi
    yyipp = yyipp - ((hi + hip) * H) / ((him + hi) * hi * hip) * yip
    yyipp = yyipp + ((2.0_dp * hip + hi) * H + hip * (hi + hip)) / ((hi + hip) * hip * H) * yipp

    B2 = (hi + hip)**2 * (abs(yyip - yyi) / hi - abs(yyi - yyim) / him)** 2
    B3 = (him + hi)**2 * (abs(yyipp - yyip) / hip - abs(yyip - yyi) / hi)** 2

  end subroutine weno4_B

  subroutine weno4_q(x, xim, xi, xip, xipp, yim, yi, yip, yipp, q2, q3)
    implicit none
    real(dp), intent(in) :: x, xim, xi, xip, xipp, yim, yi, yip, yipp
    real(dp), intent(out) :: q2, q3
    real(dp) :: him, hi, hip

    him = xi - xim
    hi = xip - xi
    hip = xipp - xip

    q2 = yim * ((x - xi) * (x - xip)) / (him * (him + hi))
    q2 = q2 - yi * ((x - xim) * (x - xip)) / (him * hi)
    q2 = q2 + yip * ((x - xim) * (x - xi)) / ((him + hi) * hi)

    q3 = yi * ((x - xip) * (x - xipp)) / (hi * (hi + hip))
    q3 = q3 - yip * ((x - xi) * (x - xipp)) / (hi * hip)
    q3 = q3 + yipp * ((x - xi) * (x - xip)) / ((hi + hip) * hip)

  end subroutine weno4_q

end module WENO4_mod
