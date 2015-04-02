pro heat2d

	nX = 40
	nY = 20

	eqdskFileName = 'g122976.03021'
	g = readgeqdsk(eqdskFileName,/noTor)

	size = 0.2
	r0 = g.rmaxis
	z0 = 0.0

	xMin = r0-size 
	xMax = r0+size
	x = fIndGen(nX)/(nX-1)*(xMax-xMin)+xMin
	x2D = rebin(x,nX,nY)
	dX = x[1]-x[0]

	yMin = z0-size 
	yMax = z0+size
	y = fIndGen(nY)/(nY-1)*(yMax-yMin)+yMin
	y2D = transpose(rebin(y,nY,nX))
	dY = y[1]-y[0]
	
	_s = sqrt((x2D-r0)^2+(y2D-z0)^2)
	_sig = 0.5
	T = exp(-_s/_sig^2) ; initial condition
	F = T*0
	kx = fltArr(nX,nY) + 0.06 ; diffusion coefficent
	ky = fltArr(nX,nY) + 20.42 ; diffusion coefficent

	cfl = 40 ; must be < 0.5 for this shitty explicit forward Euler time differencing

	; seems to work fine for large CFL numbers here, not sure why exactly

	dt = min( [cfl * dX^2 / kx, cfl * dY^2 / ky] )

	nT = 500L

	oVid = IDLffVideoWrite('heat2d.webm')
	width = 400
	height = 400
	fps = 20 
	vidStream = oVid.AddVideoStream(width, height, fps)

	xRange=[xMin,xMax]
	yRange=[yMin,yMax]
	nLevs = 11
	scale = max(T)
	levels = fIndGen(nLevs)/(nLevs-1)*scale
	colors = reverse(bytScl(levels, top=253)+1)
	c=contour(T, x, y, aspect_ratio=1.0, $
			dimensions=[width,height], $
			/fill, c_value=levels, rgb_indices=colors,$
			rgb_table=3,yRange=yRange,xRange=xRange)
	psi_levels = (fIndGen(10)+1)/10
	c=contour(g.psizr, g.r, g.z,/over,color='w',n_levels = 200)

	; Solve using the 1D set

	for i=0,nX-1 do begin
		for j=0,nY-1 do begin

			; generate the field line

			_n = 100
			dS = 0.01
			ThisPoint = [x[i],0,y[j]]
			g=readgeqdsk(eqdskFileName,/noTor, $
				FieldLineIn = ThisPoint, $
		    	FieldLineTraceDir = 1, $
    			FieldLineTraceDS = dS, $
    			FieldLineTraceNSteps = _n, $
				FieldLine_CYL = _line1)

			g=readgeqdsk(eqdskFileName,/noTor, $
				FieldLineIn = ThisPoint, $
		    	FieldLineTraceDir = -1, $
    			FieldLineTraceDS = dS, $
    			FieldLineTraceNSteps = _n, $
				FieldLine_CYL = _line2)
		
			fLine = [[reverse(_line1[*,0:-2],2)],[_line2[*,1:-2]]]
	
		endfor
	endfor
stop

	; Solve the 2D problem directly on a Cartesian grid
	
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
