function testheat_1d,T0,S,nt,nx

dene=5d19 ;Assume contant density
gamma = 7 ;(total) sheath heat transmission coeff
L = 20.0d ;Total length of 1D system (typical d3d connection legnth)
k0=2000.0d0 ;factor for parallel e-heat conductivity
cfl = 0.5 ;Mine seems to be meaningul, don't know what DLG's deal is...

x =L*dindgen(nx)/(nx-1)
dx = L/(nx-1)

T = dblarr(nx,nt)
T(*,0)=T0
qp12 = dblarr(nx-1)
kap = qp12

T(0,0) = 2.0d
T(nx-1,0) = 2.0d
time = 0.0d

for n=0L,nt-2 do begin
   for i = 0,nx-2 do begin
      Tp12 = 0.5d0*(T(i,n)+T(i+1,n))
      kap(i) = k0*Tp12^2.5 ;~Spitzer conductivity
      qp12(i) = -kap(i)*(T(i+1,n)-T(i,n))/dx 
   endfor
   dt = min((cfl*dx^2)/kap)
   time = [time,time(n)+dt]
; Get particle fluxes based on assumed density, current temperature,
; +Bohm condition
   pflux_left = dene*sqrt((2*T(0,n)*1.6d-19)/(2*1.67d-27)) 
   pflux_right = dene*sqrt((2*T(nx-1,n)*1.6d-19)/(2*1.67d-27))
; Get sheath heat flux based on particle flux, current temperature,
; sheath heat transmission factor
   qleft = -gamma*T(0,n)*pflux_left*1.6d-19
   qright = gamma*T(nx-1,n)*pflux_right*1.6d-19

   for i = 1,nx-2 do begin
      T(i,n+1) = T(i,n)+dt*(S-(qp12(i)-qp12(i-1))/dx)
   endfor
;Apply sheath heat transmission BC's
   T(0,n+1) = (qleft*dx*2.0d)/(kap(0)*3.0d) + 4.0*T(1,n+1)/3.0d - T(2,n+1)/3.0d
   T(nx-1,n+1) = (-qright*dx*2.0d)/(kap(nx-2)*3.0d) + 4.0*T(nx-2,n+1)/3.0d - T(nx-3,n+1)/3.0d
endfor

return,{T:T,qp12:qp12,time:time,x:x,pflux_left:pflux_left,pflux_right:pflux_right}
end



function parallel_heat,T0,S,tfinal,x,bctype=bctype,dene=dene
; Now do it again but more flexible for open/closed field lines
; Still using uniform length spacing, but L is now an input

if not keyword_set(bctype) then bctype=1; (1=field line on surface; 2=closed; dTds=0 BC for now)

if not keyword_set(dene) then dene=5d19 ;Assume contant density
gamma = 7 ;(total) sheath heat transmission coeff
;L = 20.0d ;Total length of 1D system (typical d3d connection legnth)
k0=2000.0d0 ;factor for parallel e-heat conductivity
cfl = 0.5 ;Mine seems to be meaningul, don't know what DLG's deal is...
itmax=10000
L = max(x)-min(x)
nx = n_elements(x)
;x =L*dindgen(nx)/(nx-1)
dx = L/(nx-1)

;T = dblarr(nx,nt)
;T(*,0)=T0
Told = T0
Tnew = T0
;T = T0
qp12 = dblarr(nx-1)
kap = qp12

;T(0) = 2.0d
;T(nx-1) = 2.0d
time = 0.0d
iter = 0

while time lt tfinal and (iter lt itmax) do begin
   iter = iter + 1
   for i = 0,nx-2 do begin
      Tp12 = 0.5d0*(Told(i)+Told(i+1))
      kap(i) = k0*Tp12^2.5 ;~Spitzer conductivity
      qp12(i) = -kap(i)*(Told(i+1)-Told(i))/dx 
   endfor
   dt = min((cfl*dx^2)/kap)
   time = time+dt
; Get particle fluxes based on assumed density, current temperature,
; +Bohm condition
   pflux_left = dene*sqrt((2*Told(0)*1.6d-19)/(2*1.67d-27)) 
   pflux_right = dene*sqrt((2*Told(nx-1)*1.6d-19)/(2*1.67d-27))
; Get sheath heat flux based on particle flux, current temperature,
; sheath heat transmission factor
   qleft = -gamma*Told(0)*pflux_left*1.6d-19
   qright = gamma*Told(nx-1)*pflux_right*1.6d-19
   if bctype eq 2 then begin
      qleft=0.0d
      qright=0.0d
   endif

   for i = 1,nx-2 do begin
      Tnew(i) = Tnew(i)+dt*(S(i)-(qp12(i)-qp12(i-1))/dx)
   endfor
;Apply sheath heat transmission BC's
   Tnew(0) = (qleft*dx*2.0d)/(kap(0)*3.0d) + 4.0*Tnew(1)/3.0d - Tnew(2)/3.0d
   Tnew(nx-1) = (-qright*dx*2.0d)/(kap(nx-2)*3.0d) + 4.0*Tnew(nx-2)/3.0d - Tnew(nx-3)/3.0d
   Told = Tnew
endwhile
T = Tnew

return,{T:T,qp12:qp12,time:time,x:x,pflux_left:pflux_left,pflux_right:pflux_right,niter:iter}
end

function build_fieldlines,g,r0=r0,r1=r1,nr=nr,z0=z0,z1=z1,nz=nz,nx=nx

  if not keyword_set(r0) then r0 = 1.0
  if not keyword_set(r1) then r1 = 2.4
  if not keyword_set(nr) then nr = 10
  if not keyword_set(z0) then z0 = -1.5
  if not keyword_set(z1) then z1 = 1.5
  if not keyword_set(nz) then nz = 10
  if not keyword_set(nx) then nx = 20

  r = r0 + (r1-r0)*dindgen(nr)/(nr-1)
  z = z0 + (z1-z0)*dindgen(nz)/(nz-1)
  
  rmp = build_nstx_rwmcoils(0d3)

  rfl = dblarr(nr,nz,nx)
  zfl = dblarr(nr,nz,nx)
  xfl = dblarr(nr,nz,nx)

  il = where(g.lim(0,*) gt 1e-4)
  rlim = g.lim(0,il)
  zlim = g.lim(1,il)
  vesobj = obj_new('IDLanROI',rlim,zlim)
  
  for i = 0,nz-1 do begin
     f=follow_fieldlines_phi(g,rmp,0.5,0.5,1,1,0.0d,!dpi/180,100*max(g.qpsi),rstart=r,zstart=z(i),/quiet)
     f1=follow_fieldlines_phi(g,rmp,0.5,0.5,1,1,0.0d,-!dpi/180,100*max(g.qpsi),rstart=r,zstart=z(i),/quiet)
     nfl = n_elements(f.r(0,*))
     x = f.r*cos(f.phi)
     y = f.r*sin(f.phi)
     z = f.z
     x1 = f1.r*cos(f1.phi)
     y1 = f1.r*sin(f1.phi)
     z1 = f1.z
;stop
     for j=0,nr-1 do begin
        ds = [0,reform(sqrt((x(j,1:nfl-1)-x(j,0:nfl-2))^2 + (y(j,1:nfl-1)-y(j,0:nfl-2))^2 + (z(j,1:nfl-1)-z(j,0:nfl-2))^2))]
        ds1 = [0,reform(sqrt((x1(j,1:nfl-1)-x1(j,0:nfl-2))^2 + (y1(j,1:nfl-1)-y1(j,0:nfl-2))^2 + (z1(j,1:nfl-1)-z1(j,0:nfl-2))^2))]
        stop
     endfor
        
  endfor


  stop
end
