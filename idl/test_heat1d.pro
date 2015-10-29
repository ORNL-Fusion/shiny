function Ta, x,t,L,D

	return, sin(!pi*x/L) * exp(-D * !pi^2 * t / L^2 )		
	
end

pro test_heat1d

	; Using http://www.nada.kth.se/~jjalap/numme/FDheat.pdf as the test case.

	nX = 100
	minX = 0.0
	maxX = 1.0
	L = maxX-minX
	dx = L/(nX-1)
	x = fIndGen(nX)*dX+minX

	kx = 0.1 
	k = fltArr(nX)+kx
	T = Ta(x,0,L,kx)

	Q = fltArr(nX)

	; CFL = D*dt/dx^2
	
	CFL = 0.4
	dt_est = CFL * dx^2 / max(k)

	time = 0.5
	nT = ceil(time/dt_est)
	dt = time/nT 

	print, 'Euler time advance parameters ...'
	print, 'nT: ', nT
	print, 'dt: ', dt
	print, 'CFL: ', kx*dt/dx^2

	Ti = T
	Tex = heat1d(x,Ti,Q,k,dT,nT) 

	nT_CN = 10

	Ti = T
	Tim = heat1d(x,Ti,Q,k,time/nT_CN,nT_CN,CN=1)

	Ti = T
	Tbt = heat1d(x,Ti,Q,k,time/nT_CN,nT_CN,BT=1)

	Tana = Ta(x,time,L,kx)

	range = max(Tana)
	p=plot(x,T)
	p=plot(x,Tana,/over,color='y',thick=4,yRange=[0,1]*range)
	p=plot(x,Tex,/over,color='r',lineStyle='--',thick=2)
	p=plot(x,Tim,/over,color='b')
	p=plot(x,Tbt,/over,color='g')

	print, 'Time: ', time
	print, 'r: ', kx*dt/(dx^2)
	print, 'L2 Norm Euler: ',norm(Tex-Tana,lNorm=2)
	print, 'L2 Norm Implicit Euler: ',norm(Tbt-Tana,lNorm=2)
	print, 'L2 Norm C-N: ',norm(Tim-Tana,lNorm=2)

	print, 'dt_CN / dt_Euler: ', (time/nT_CN) / dt
	stop

end
