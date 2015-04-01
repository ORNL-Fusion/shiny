pro heat2d

	nX = 100
	nY = 50

	xMin = -1 
	xMax = +1
	x = fIndGen(nX)/(nX-1)*(xMax-xMin)+xMin
	x2D = rebin(x,nX,nY)
	dX = x[1]-x[0]

	yMin = -1 
	yMax = +1
	y = fIndGen(nY)/(nY-1)*(yMax-yMin)+yMin
	y2D = transpose(rebin(y,nY,nX))
	dY = y[1]-y[0]
	

	T = cos(sqrt(x2D^2+y2D^2)*3)>0 ; initial condition
	F = T*0
	kx = fltArr(nX,nY) + 0.06 ; diffusion coefficent
	ky = fltArr(nX,nY) + 0.22 ; diffusion coefficent

	cfl = 4 ; must be < 0.5 for this shitty explicit forward Euler time differencing

	; seems to work fine for large CFL numbers here, not sure why exactly

	dt = min( [cfl * dX^2 / kx, cfl * dY^2 / ky] )

	nT = 500L

	oVid = IDLffVideoWrite('heat2d.webm')
	width = 400
	height = 400
	fps = 20 
	vidStream = oVid.AddVideoStream(width, height, fps)

	nLevs = 11
	scale = max(T)
	levels = fIndGen(nLevs)/(nLevs-1)*scale
	colors = reverse(bytScl(levels, top=253)+1)

	for _t = 0, nT - 1 do begin

		; plot time evolving solution at a subset of times
		if _t mod 5 eq 0 then begin
			if _t gt 0 then c.erase

			c=contour(T, x, y, aspect_ratio=1.0, /buffer, $
					dimensions=[width,height], $
					/fill, c_value=levels, rgb_indices=colors,$
					rgb_table=3)
			frame = c.CopyWindow()
			!null = oVid.put(vidStream, frame)
		endif

		T[1:-2,1:-2] = $
				T[1:-2,1:-2] $
				+ dt * kx[1:-2] * dt / dX^2  * ( T[0:-3,1:-2] - 2*T[1:-2,1:-2] + T[2:-1,1:-2] ) $
				+ dt * ky[1:-2] * dt / dY^2  * ( T[1:-2,0:-3] - 2*T[1:-2,1:-2] + T[1:-2,2:-1] ) $
				+ dt * F[1:-2]

		;; Apply BCs

		;; left side : dT/dX = 0 second order accurate forward difference
		;T[0,*] = (-2*T[1,*] + 0.5*T[2,*])/(-1.5) ;  

		;; right side : dT/dX = 0 second order accurate backward difference
		;T[-1,*] = (+2*T[-2,*] - 0.5*T[-3,*])/(+1.5) 

		;; top : dT/dY = 0 second order accurate forward difference
		;T[*,0] = (-2*T[*,1] + 0.5*T[*,2])/(-1.5)  

		;; bottom : dT/dY = 0 second order accurate backward difference
		;T[*,-1] = (+2*T[*,-2] - 0.5*T[*,-3])/(+1.5) ; 

	endfor	

	oVid.cleanup
stop
end
