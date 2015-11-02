Module shiny_mod
  Use kind_mod, Only: rknd, iknd
  Implicit None
  Contains

    Function getPsi(x,y,nx,ny)
      Use phys_const, Only: pi
      Implicit None
      Real(rknd) :: getPsi(nx,ny) 
      Integer, Intent(In) :: nx, ny
      Real(rknd), Intent(In) :: x(nx,ny), y(nx,ny)
      getPsi = cos(pi*x) * cos(pi*y)
    End Function getPsi

    Function getQ(x,y,nx,ny)
      Use phys_const, Only: pi
      Implicit None
      Real(rknd) :: getQ(nx,ny)
      Real(rknd), Intent(In) :: x(nx,ny),y(nx,ny)
      Integer(iknd), Intent(In) :: nx,ny
      Real(rknd) :: lap(nx,ny)

      lap = -2.d0*pi**2*getPsi(x,y,nx,ny) ! analytic laplacian (psi)
      getQ = -lap
    End Function getQ

    Function tfac(kper,t)
      Use phys_const, Only: pi
      Implicit None
      Real(rknd), Intent(In) :: kper, t
      Real(rknd) :: tfac
      tfac = (1.d0 - Exp( -2.d0*kper*pi**2*t) ) / kper
    End Function tfac

    Function get_b(x,y)
      Use phys_const, Only: pi
      Implicit None
      Real(rknd) :: get_b(3), b(3), modb
      Real(rknd), Intent(In) :: x,y
      b(1) =  pi*cos(pi*x)*sin(pi*y) !bx
      b(2) = -pi*cos(pi*y)*sin(pi*x) !by
      b(3) = 0.d0                    !bz
      modb = mag(b)
      get_b = b/modb
      If (any(isnan(get_b))) get_b = b
    End Function get_b

    Function mag(v)
      Implicit None
      Real(rknd) :: mag
      Real(rknd), Intent(In) :: v(3)
      mag = Sqrt(Sum(v**2))
    End Function mag

    Function get_D(kper,kpar,bu)
      Implicit None
      Real(rknd), Intent(In) :: kper, kpar, bu(3)
      Real(rknd) :: get_D(2,2)
      get_D(1,1) = bu(1)**2*(kpar-kper) + kper
      get_D(1,2) = bu(1)*bu(2)*(kpar-kper)
      get_D(2,1) = bu(1)*bu(2)*(kpar-kper)
      get_D(2,2) = bu(2)**2*(kpar-kper) + kper
    End Function Get_D

    Function ddotgradt(d,dtdx,dtdy)      
      Implicit None
      Real(rknd) :: ddotgradt(2)
      Real(rknd), Intent(In) :: d(2,2), dtdx,dtdy
      ddotgradt(1) = d(1,1) * dtdx + d(2,1) * dtdy !qx
      ddotgradt(2) = d(1,2) * dtdx + d(2,2) * dtdy !qy
    End Function Ddotgradt

    Function getT00(nx,ny,T)
      Implicit None
      Integer(iknd), Intent(In) :: nx, ny
      Real(rknd), Intent(In) :: T(nx,ny)
      Real(rknd) :: getT00
      If (Mod(nx,2) == 1 .AND. Mod(ny,2) == 1) Then ! odd nx, odd ny
        getT00 = T((nx-1)/2,(ny-1)/2)
      ElseIf (Mod(nx,2) == 1 .AND. Mod(ny,2) == 0) Then ! odd nx, even ny
        getT00 = (T((nx-1)/2,ny/2)+T((nx-1)/2,ny/2-1))/2.d0
      ElseIf (Mod(nx,2) == 0 .AND. Mod(ny,2) == 1) Then ! even nx, odd ny
        getT00 = (T(nx/2,(ny-1)/2)+T(nx/2-1,(ny-1)/2))/2.d0
      ElseIf (Mod(nx,2) == 0 .AND. Mod(ny,2) == 0) Then ! even nx, even ny
        getT00 = (T(nx/2,ny/2)+T(nx/2-1,ny/2)+T(nx/2,ny/2-1)+T(nx/2-1,ny/2-1))/4.d0
    Endif
  End Function getT00

!  Function divq_asym
    
End Module shiny_mod
