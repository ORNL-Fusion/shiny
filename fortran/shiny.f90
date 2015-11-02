Program shiny
  Use kind_mod, Only: iknd, rknd, iknd15
  Use math_geo_module, Only: rlinspace
  Use shiny_mod, Only: getPsi, getQ, tfac, get_b, get_d, ddotgradt, getT00
  Use phys_const, Only: pi
  Implicit None
  Integer(iknd15) :: nt, kt
  Integer(iknd) :: nx, ny, i, j
  Real(rknd) :: CFL, y0, x0, xmin, xmax, dx, ymin, ymax, dy, endtime, T00, T00a
  Real(rknd) :: kper, kpar, D, dt, xp, xm, yp, ym, divq, x0_test, y0_test
  Real(rknd) :: b(3), dkper, l2
  Real(rknd), Dimension(2,2) :: D_ipj, D_ijp, D_imj, D_ijm, D_imjp, D_ipjm, D_imjm, D_ipjp
  Real(rknd), Dimension(2) :: q_ipj, q_ijm, q_ijp, q_imj, q_ipjp, q_ipjm, q_imjp, q_imjm
  Real(rknd) :: dTdx_ipj, dTdy_ipj, dTdx_ijp, dTdy_ijp, dTdx_imj, dTdy_imj, dTdx_ijm, dTdy_ijm
  Real(rknd) :: dTdx_imjm, dTdx_imjp, dTdx_ipjm, dTdx_ipjp, dTdy_ipjp, dTdy_ipjm, dTdy_imjp, dTdy_imjm
  Real(rknd), Allocatable :: x(:), x2d(:,:), y(:), y2d(:,:), psi(:,:)
  Real(rknd), Allocatable, Dimension(:,:) :: bx, by, bz, Q, T, T2, Tsolution, Tupdate
  Logical :: eqdsk, closedcase, dogunter, symmetric, do1Dset
  Character(Len=120) :: outfile
  Real(kind=4) :: tarray(2), tres, tres0

  Call Etime(tarray,tres0)
  
  nx = 32
  ny = 32

  dogunter = .true.
  symmetric = .true.
  do1Dset = .true.
  eqdsk = .false.

  kper = 1.d0
  kpar = 1.d2

  Write(*,*) 'dogunter: ',dogunter
  Write(*,*) 'symmetric', symmetric

  If (eqdsk) Then
  Else
    
    closedcase = .true.
    If (closedcase) Then 
      x0 = 0.d0
      y0 = 0.d0
    Else
      x0 = -0.5d0
      y0 = -0.5d0
    Endif
  
    xmin = -0.5d0
    xmax =  0.5d0
    Allocate(x(0:nx-1),x2D(0:nx-1,0:ny-1))
    x = rlinspace(xmin,xmax,nx)
    Do i = 0,ny-1
      x2d(:,i) = x
    Enddo
    dx = x(1) - x(0)
    Open(99,file='x.out')
    Write(99,*) x
    Close(99)
    
    ymin = -0.5d0
    ymax =  0.5d0
    Allocate(y(0:ny-1),y2d(0:nx-1,0:ny-1))
    y = rlinspace(ymin,ymax,ny)
    Do i = 0,nx-1
      y2d(i,:) = y
    Enddo
    dY = y(1) - y(0)
    Open(99,file='y.out')
    Write(99,*) y
    Close(99)
    
    Allocate(psi(0:nx-1,0:ny-1))
    psi = getPsi(x2d,y2d,nx,ny)
    
    ! Use analytic expression for b instead
    
    Allocate(bx(0:nx-1,0:ny-1),by(0:nx-1,0:ny-1),bz(0:nx-1,0:ny-1))
    bx =  pi*cos(pi*x2d)*sin(pi*y2d)
    by = -pi*cos(pi*y2d)*sin(pi*x2d)
    bz =  bx*0.d0
    
    
    Allocate(Q(0:nx-1,0:ny-1))
    Q = getQ(x2d,y2d,nx,ny)
    Allocate(T(0:nx-1,0:ny-1),T2(0:nx,0:ny-1))
    Allocate(tsolution(0:nx-1,0:ny-1),Tupdate(0:nx-1,0:ny-1))    
    T = 0.d0
    T2 = T
    Tsolution = tfac(kper,0.d0) * psi

  Endif
    
  CFL = 0.9d0 ! must be < 1.0
  D = Max(Abs(kper),Abs(kpar))
  dt = CFL * 0.125d0 * (dx**2 + dy**2) / D
  Write(*,*) 'dt = ',dt

  endtime = 1.d0
  !nt = Ceiling(endtime/dt,iknd15)
  nt = 10000_iknd15
  Write(*,*) 'nt = ',nt

  Open(99,file='run_info.out')
  Write(99,*) nx,ny,1.d0/Real(nx,rknd),1.d0/Real(ny,rknd)
  Write(99,*) kper, kpar
  Write(99,*) nt
  Write(99,*) dt, endtime
  Close(99)  

  Call Flush
  
  If (dogunter) Then    
    Do kt = 0, nt -1
      !!Tsolution = tfac(kper,dt*kt)*psi  ! dont need this unless comparing as function of time

      Tupdate = 0.d0

      Do i = 1,nx-2
        xp = x(i)+dx/2.d0
        xm = x(i)-dx/2.d0

        Do j = 1,ny-2          
          yp = y(j)+dy/2.d0
          ym = y(j)-dy/2.d0

          If (.not. symmetric) Then
            ! ASYMMETRIC SCHEME

            ! Temperature gradient terms

            dTdx_ipj = (T(i+1,j) - T(i,j))/dx
            dTdy_ipj = (T(i+1,j+1) + T(i,j+1) - T(i,j-1) - T(i+1,j-1)) / (4.d0*dy)
            dTdx_ijp = (T(i+1,j+1) + T(i+1,j) - T(i-1,j+1) - T(i-1,j)) / (4.d0*dx)
            dTdy_ijp = (T(i,j+1) - T(i,j))/dy
            
            dTdx_imj = (T(i,j) - T(i-1,j))/dx
            dTdy_imj = (T(i,j+1) + T(i-1,j+1) - T(i-1,j-1) - T(i,j-1)) / (4.d0*dy)
            dTdx_ijm = (T(i+1,j) + T(i+1,j-1) - T(i-1,j) - T(i-1,j-1)) / (4.d0*dx)
            dTdy_ijm = (T(i,j) - T(i,j-1))/dy
            
            ! Heat conduction term
            ! q = -D.\/(T) where D is a 2x2 tensor
            
            ! x plus / minus terms
            b = get_b(xp,y(j))

            D_ipj = get_D(kPer,kPar,b)
            q_ipj = -DdotGradT( D_ipj, dTdx_ipj, dTdy_ipj )

            b = get_b(xm,y(j))

            D_imj = get_D(kPer,kPar,b)
            q_imj = -DdotGradT( D_imj, dTdx_imj, dTdy_imj )

            ! y plus / minus terms
            b = get_b(x(i),yp)

            D_ijp = get_D(kPer,kPar,b)
            q_ijp = -DdotGradT( D_ijp, dTdx_ijp, dTdy_ijp )

            b = get_b(x(i),ym)            

            D_ijm = get_D(kPer,kPar,b)
            q_ijm = -DdotGradT( D_ijm, dTdx_ijm, dTdy_ijm )

            ! Diffusion term
            ! -\/.q

            divq = (q_ipj(1) - q_imj(1))/dx + (q_ijp(2)-q_ijm(2))/dy

            outfile = 'Tasym.out'
          Else
              
            ! SYMMETRIC
            dTdx_ipjp = (T(i+1,j+1) + T(i+1,j) - T(i,j+1) - T(i,j)) / (2.d0*dx)    
            dTdy_ipjp = (T(i,j+1) + T(i+1,j+1) - T(i+1,j) -T(i,j)) / (2.d0*dy)
            
            dTdx_ipjm = (T(i+1,j) + T(i+1,j-1) - T(i,j) - T(i,j-1)) / (2.d0*dx)    
            dTdy_ipjm = (T(i,j) + T(i+1,j) - T(i+1,j-1) -T(i,j-1)) / (2.d0*dy)
            
            dTdx_imjp = (T(i,j+1) + T(i,j) - T(i-1,j+1) - T(i-1,j)) / (2.d0*dx)    
            dTdy_imjp = (T(i-1,j+1) + T(i,j+1) - T(i,j) -T(i-1,j)) / (2.d0*dy)
            
            dTdx_imjm = (T(i,j) + T(i,j-1) - T(i-1,j) - T(i-1,j-1)) / (2.d0*dx)    
            dTdy_imjm = (T(i-1,j) + T(i,j) - T(i,j-1) -T(i-1,j-1)) / (2.d0*dy)
                        
            b = get_b(xp,yp)            
            D_ipjp = get_D(kPer,kPar,b)
            q_ipjp = -DdotGradT( D_ipjp, dTdx_ipjp, dTdy_ipjp )
            
            b = get_b(xp,ym)
            D_ipjm = get_D(kPer,kPar,b)
            q_ipjm = -DdotGradT( D_ipjm, dTdx_ipjm, dTdy_ipjm )
            
            b = get_b(xm,yp)            
            D_imjp = get_D(kPer,kPar,b)
            q_imjp = -DdotGradT( D_imjp, dTdx_imjp, dTdy_imjp )
            
            b = get_b(xm,ym)            
            D_imjm = get_D(kPer,kPar,b)
            q_imjm = -DdotGradT( D_imjm, dTdx_imjm, dTdy_imjm )
            
            
            divq = (q_ipjp(1) + q_ipjm(1) - q_imjp(1) - q_imjm(1)) / (2.d0*dx) &
                 + (q_ipjp(2) + q_imjp(2) - q_imjm(2) - q_ipjm(2)) / (2.d0*dy)
            
            outfile = 'Tsym.out'
          Endif

          Tupdate(i,j) = dt * (-divq + Q(i,j))
          
        Enddo
      Enddo

      T = T + Tupdate
      If (Mod(kt,Int(nt*0.01,iknd15)) .eq. 0) Then
        Write(*,'(i0,x,i0,x,f5.1)') kt,nt, Real(kt,rknd)/Real(nt,rknd)*100.d0
        Call flush
      Endif
      
    Enddo !time

    ! Analytic value
    
    Tsolution = tfac(kper,nt*dt)*psi

    ! Evaluate T00 -- depends if grid is odd or even
    T00 = getT00(nx,ny,T)
    T00a = getT00(nx,ny,Tsolution)
    dkPer = 1.d0/T00-kper          

    Write(*,*) 'T00  ',T00
    Write(*,*) 'T00a ',T00a
    
    Write(*,*) 'endtime =',nt*dt  
      
    l2 = Sqrt(Sum((Tsolution-T)**2))
    Write(*,*) 'dkper ',dkper
    Write(*,*) 'L2 ',l2
    
    Open(99,file=outfile)
    Write(99,*) T
    Close(99)

    Open(99,file='Tanalytic.out')
    Write(99,*) Tsolution
    Close(99)
    
  Endif ! dogunter
  
!  If (do1dset) Then
!    Do i = 1,nx-2
!      Do j =1,ny-2
!        
!
!  Endif

  
  Deallocate(bx,by,bz,Q,T,T2,Tsolution,Tupdate)
  Deallocate(x,y,x2d,y2d,psi)
  Call Flush
  Call Etime(tarray,tres)
  Write(*,*) ' This took ',tres-tres0,' seconds'

End Program shiny
