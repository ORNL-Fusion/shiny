function dr_ds, s, r

	common bfield, b

	_x = (r[0]-min(b.x))/(max(b.x)-min(b.x))*(n_elements(b.x)-1)
	_y = (r[1]-min(b.y))/(max(b.y)-min(b.y))*(n_elements(b.y)-1)

	bx = interpolate(b.bx,_x,_y,cubic=-0.5)	
	by = interpolate(b.by,_x,_y,cubic=-0.5)	

	return, [bx,by] 

end

pro getDomain, x0, y0, d, _b, bndry, perp, crash=crash

    bndry_x = bndry.x
    bndry_y = bndry.y
    boundary = bndry.b

	_n = (d.N+1)/2
	; generate the field line

	ThisPoint = [x0,y0]

	_line1 = dlg_fieldLineTrace_xy(_b,ThisPoint, dir = -1, dS=d.dS, nS=_n, perp=perp )
	_line2 = dlg_fieldLineTrace_xy(_b,ThisPoint, dir = +1, dS=d.dS, nS=_n, perp=perp )

	fLine_XY = [[reverse(_line1[*,0:-2],2)],[_line2[*,1:-2]]]
    fLine_XY = fLine_XY[0:1,*]
	_s = fIndGen(_n)*d.dS
	s = [reverse(-_s[1:-1]),_s[0:-1]]

	; Check intersection with boundary

	isInsideDomain = boundary->ContainsPoints(fLine_XY[0,*],fLine_XY[1,*])
	iiOutside = where(isInsideDomain eq 0, iiOutsideCnt)
	iiThisPt = _n-1 ; this is the middle point

    if keyword_set(crash) then stop

	; Extract that part of the line within the boundary
	; accounting for craziness, multiple intersections, etc

	nF = n_elements(fLine_XY[0,*])

	iin = !null
	_cont = 1
	_ii = 1
	while _cont eq 1 do begin
		if iiThisPt-_ii lt 0 then begin
				_cont = 0
		endif else begin
			if isInsideDomain[iiThisPt-_ii] then iin=[iin,iiThisPt-_ii] else _cont = 0
			++_ii
		endelse
	endwhile

	iip = !null
	_cont = 1
	_ii = 1
	while _cont eq 1 do begin
		if iiThisPt+_ii ge nF then begin
				_cont = 0
		endif else begin
			if isInsideDomain[iiThisPt+_ii] then iip=[iip,iiThisPt+_ii] else _cont = 0
			++_ii
		endelse
	endwhile

	iiUse = [reverse(iin),iiThisPt,iip]
	fLine_XY_all = fLine_XY
	fLine_XY = fLine_XY[*,iiUse]
	s_all = s
	s = s[iiUse]
	sLeft = s[0]
	sRight = s[-1]
	n_fl = d.N

	LeftEnd = [fLine_XY[0,0],fLine_XY[1,0]]	
	RightEnd = [fLine_XY[0,-1],fLine_XY[1,-1]]	


	; Now find actual intersection points with the wall(s)

	ReGrid = 0

	if min(iiUse) ne 0 then begin ; left end needs fixing 

		ReGrid = 1

		iiOut = min(iiUse)-1
		iiIn = min(iiUse)
		x1 = fLine_XY_all[0,iiOut]
		y1 = fLine_XY_all[1,iiOut]
		x2 = fLine_XY_all[0,iiIn]
		y2 = fLine_XY_all[1,iiIn]

		LeftEnd = dlg_line_polygon(x1,y1,x2,y2,bndry_x,bndry_y)

		sExtraBit = sqrt ( (x2-LeftEnd[0])^2+(y2-LeftEnd[1])^2) 
		sLeft = s_all[iiIn] - sExtraBit

	endif 

	if max(iiUse) ne n_elements(fLine_XY_all[0,*])-1 then begin ; right end needs fixing 

		ReGrid = 1

		iiOut = max(iiUse)+1
		iiIn = max(iiUse)
		x1 = fLine_XY_all[0,iiOut]
		y1 = fLine_XY_all[1,iiOut]
		x2 = fLine_XY_all[0,iiIn]
		y2 = fLine_XY_all[1,iiIn]

		RightEnd = dlg_line_polygon(x1,y1,x2,y2,bndry_x,bndry_y,crash=crash)

		sExtraBit = sqrt ( (x2-RightEnd[0])^2+(y2-RightEnd[1])^2 ) 
		sRight = s_all[iiIn] + sExtraBit

	endif 

	; Interpolate field line to new constant dS such that it ends at the boundary(s)

	if ReGrid then begin

		sNew = fIndGen(n_fl)/(n_fl-1)*(sRight-sLeft)+sLeft

		_fLine_x = interpol(fLine_XY_all[0,*],s_all,sNew)
		_fLine_y = interpol(fLine_XY_all[1,*],s_all,sNew)

		s = sNew
		fLine_XY_all[0,*] = _fLine_x
		fLine_XY_all[1,*] = _fLine_y

	endif

    d.s = s
    d.x = fLine_XY_all[0,*]
    d.y = fLine_XY_all[1,*]
    if keyword_set(crash) then stop

end


pro grad, f,x,y,gradX,gradY

    nX = n_elements(x)
    nY = n_elements(y)

    hx = x[1]-x[0]
    hy = y[1]-y[0]

    gradX = fltArr(nX,nY)
    gradY = fltArr(nX,nY)

    for i=0,nX-1 do begin
        for j=0,nY-1 do begin

            ; df_dX 

            if i gt 0 and i lt nX-1 then begin
                gradX[i,j] = (-f[i-1,j]+f[i+1,j])/(2*hx)
            endif

            if i eq 0 then begin
                gradX[i,j] = (-f[i,j]+f[i+1,j])/(hx)
            endif

            if i eq nX-1 then begin
                gradX[i,j] = (-f[i-1,j]+f[i,j])/(hx)
            endif

            ; df_dY

            if j gt 0 and j lt nY-1 then begin
                gradY[i,j] = (-f[i,j-1]+f[i,j+1])/(2*hy)
            endif

            if j eq 0 then begin
                gradY[i,j] = (-f[i,j]+f[i,j+1])/(hy)
            endif

            if j eq nY-1 then begin
                gradY[i,j] = (-f[i,j-1]+f[i,j])/(hy)
            endif


        endfor
    endfor

end

function laplacian, f, x, y

    nX = n_elements(x)
    nY = n_elements(y)

    hx = x[1]-x[0]
    hy = y[1]-y[0]

    lapX = fltArr(nX,nY)
    lapY = fltArr(nX,nY)

    for i=0,nX-1 do begin
        for j=0,nY-1 do begin

            ; d2f_dX 

            if i gt 0 and i lt nX-1 then begin
                lapX[i,j] = (+f[i-1,j] - 2*f[i,j] + f[i+1,j])/(hx^2)
            endif

            if i eq 0 then begin
                lapX[i,j] = (+f[i,j] -2*f[i+1,j] + f[i+2,j])/(hx^2)
            endif

            if i eq nX-1 then begin
                lapX[i,j] = (+f[i-2,j] -2*f[i-1,j] + f[i,j])/(hx^2)
            endif

            ; d2f_dY

            if j gt 0 and j lt nY-1 then begin
                lapY[i,j] = (+f[i,j-1] - 2*f[i,j] + f[i,j+1])/(hy^2)
            endif

            if j eq 0 then begin
                lapY[i,j] = (+f[i,j] -2*f[i,j+1] + f[i,j+2])/(hy^2)
            endif

            if j eq nY-1 then begin
                lapY[i,j] = (+f[i,j-2] -2*f[i,j-1] + f[i,j])/(hy^2)
            endif

        endfor
    endfor

    return, lapX + lapY

end


pro heat2d

    ; I'm unsure why I need to do this, bah.
    resolve_routine, 'interpb', /either, /compile_full_file, /no_recompile

	nX = 40
	nY = 41

    eqdsk = 0
    if eqdsk then begin
	    eqdskFileName = 'g122976.03021'
	    g = readgeqdsk(eqdskFileName,/noTor)

        psi = -g.psizr

	    size = 0.19
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

    endif else begin

        closedCase = 1

        if closedCase then begin
            x0 = 0
            y0 = 0
        endif else begin
            x0 = -0.5 
            y0 = -0.5
        endelse

        xMin = -0.5 
	    xMax = +0.5
	    x = fIndGen(nX)/(nX-1)*(xMax-xMin)+xMin
	    x2D = rebin(x,nX,nY)
	    dX = x[1]-x[0]

	    yMin = -0.75 
	    yMax = +0.75
	    y = fIndGen(nY)/(nY-1)*(yMax-yMin)+yMin
	    y2D = transpose(rebin(y,nY,nX))
	    dY = y[1]-y[0]
	
        psi = sin(2*!pi*x2D) * cos(2*!pi*y2D)
        grad, psi, x, y, gradX, gradY
        lap = laplacian( psi, x, y)


        bx = -gradY
        by = +gradX
        bz = bx*0

        ; Generate kx / ky from kPer / kPar

        kPer = 1
        kPar = 1e9

        bMag = sqrt(bx^2+by^2+bz^2) 
        bxU = bx / bMag
        byU = by / bMag
        bzU = bz / bMag

        kx = fltArr(nX,nY)
        ky = fltArr(nX,nY)

        for i=0,nX-1 do begin
            for j=0,nY-1 do begin
        
                parU = [bxU[i,j],byU[i,j],bzU[i,j]]
                zU = [0,0,1]
                perU = crossp(parU,zU)
                perU = perU / sqrt(perU[0]^2+perU[1]^2+perU[2]^2)

                kx[i,j] = 0 
                kx[i,j] = kx[i,j] + abs(kPer * perU[0])
                kx[i,j] = kx[i,j] + abs(kPar * parU[0])

                ky[i,j] = 0
                ky[i,j] = ky[i,j] + abs(kPer * perU[1])
                ky[i,j] = ky[i,j] + abs(kPar * parU[1])

            endfor
        endfor 

        ;Q = -kPer * lap * 5e1
        Q = lap * 5e1

        ;Q[*] = 0

        r = sqrt((x2D-x0)^2+(y2D-y0)^2)
        T = (1-r^3)*0
        T2 = T 
        TSolution = sin(2*!pi*x2D) * cos(2*!pi*y2D)

        c=contour(psi,x,y,layout=[3,2,1],/fill,title='psi')
        c=contour(TSolution,x,y,layout=[3,2,2],/current,/fill,title='T')
        c=contour(kx,x,y,layout=[3,2,3],/current,/fill,title='kx',rgb_table=50)
        c=contour(ky,x,y,layout=[3,2,4],/current,/fill,title='ky',rgb_table=50)
        c=contour(Q,x,y,layout=[3,2,5],/current,/fill,title='Q')
        v=vector(bx,by,x,y,layout=[3,2,6],/current,auto_color=1,rgb_table=10,auto_subsample=1,title='B')

    endelse


	CFL = 0.9 ; must be < 1 for this shitty explicit forward Euler time differencing
    _D = max(abs([kx,ky]))
    dt = CFL * ( 1.0 / 8.0 ) * (dx^2 + dy^2) / _D

	nT = 500L

	width = 400
	height = 400
	oVid = IDLffVideoWrite('heat2d.webm')
	fps = 4 
	vidStream = oVid.AddVideoStream(width, height, fps)

	xRange=[x[0],x[-1]]
	yRange=[y[0],y[-1]]
	nLevs = 11
    c=contour(T,x,y,/fill,/buffer,dimensions=[width,height],rgb_table=50)

    frame = c.CopyWindow()
	!null = oVid.put(vidStream, frame)
	
	; Solve the 2D problem directly on a Cartesian grid

	T[0,*]  = TSolution[0,*]
    T[-1,*] = TSolution[-1,*]
    T[*,0]  = TSolution[*,0]
    T[*,-1] = TSolution[*,-1]

	for _t = 0, nT - 1 do begin


		T[1:-2,1:-2] = $
				T[1:-2,1:-2] $
				+ kx[1:-2,1:-2] * dt / dX^2  * ( T[0:-3,1:-2] - 2*T[1:-2,1:-2] + T[2:-1,1:-2] ) $
				+ ky[1:-2,1:-2] * dt / dY^2  * ( T[1:-2,0:-3] - 2*T[1:-2,1:-2] + T[1:-2,2:-1] ) $
				+ dt * Q[1:-2,1:-2]

		;; Apply BCs

		;; left side : dT/dX = 0 second order accurate forward difference
		;T[0,*] = (-2*T[1,*] + 0.5*T[2,*])/(-1.5) ;  

		;; right side : dT/dX = 0 second order accurate backward difference
		;T[-1,*] = (+2*T[-2,*] - 0.5*T[-3,*])/(+1.5) 

		;; top : dT/dY = 0 second order accurate forward difference
		;T[*,0] = (-2*T[*,1] + 0.5*T[*,2])/(-1.5)  

		;; bottom : dT/dY = 0 second order accurate backward difference
		;T[*,-1] = (+2*T[*,-2] - 0.5*T[*,-3])/(+1.5) ; 

		; plot time evolving solution at a subset of times
		if _t mod 50 eq 0 then begin
			if _t gt 0 then c.erase
            print, _t, nT    
            c=contour(T,x,y,/fill,/buffer,/current,dimensions=[width,height],rgb_table=50)
			frame = c.CopyWindow()
			!null = oVid.put(vidStream, frame)
		endif

	endfor	

	oVid.cleanup

	p=plot(TSolution[*],T[*],symbol="Circle",lineStyle='none')


	; Solve using the 1D set
    ; ----------------------

	bndry_x = [xMin,xMax,xMax,xMin,xMin]
	bndry_y = [yMin,yMin,yMax,yMax,yMin]
	boundary = Obj_New('IDLanROI',bndry_x,bndry_y)

    bndry = {x:bndry_x,y:bndry_y,b:boundary}

	common bfield, b
	b = { x: x, y: y, bx : bx, by : by, bz: bx*0, psi:psi }

	_n = 30
	dS = 0.025

	__n = 2*_n-1
	d1 = { $
		N: __n, $
		x: fltArr(__n), $
		y: fltArr(__n), $
		dS: dS, $
		s: fltArr(__n), $
		kx: fltArr(__n), $
		ky: fltArr(__n) }	

	pt = { x: 0.0, y:0.0, par : d1, per : d1 } ; grid point structure contains both per and par domains	   

	points = replicate( pt, nX, nY )

    c=contour(psi,x,y,/fill)

	for i=1,nX-2 do begin
		for j=1,nY-2 do begin
			points[i,j].x = x[i]
			points[i,j].y = y[j]

            par = points[i,j].par
			getDomain, x[i],y[j], par, b, bndry, 0
            points[i,j].par = par
            
            per = points[i,j].per
            crash=0
            ;if i eq 17 and j eq 10 then crash = 1
			getDomain, x[i],y[j], per, b, bndry, 1, crash=crash
            points[i,j].per = per

            plotMod = 4 
            if i mod plotMod eq 0 and j mod plotMod eq 0 then begin
                    p=plot(points[i,j].per.x,points[i,j].per.y,/over)
                    p=plot(points[i,j].par.x,points[i,j].par.y,/over)
            endif 
		endfor
	endfor

	T2[0,*]  = TSolution[0,*]
    	T2[-1,*] = TSolution[-1,*]
    	T2[*,0]  = TSolution[*,0]
    	T2[*,-1] = TSolution[*,-1]

    nItr = 500 

    c=contour(T2,x,y,/fill,dimensions=[width,height],rgb_table=50)

    for itr=0, nItr-1 do begin

        ; Solve parallel 

        T2_copy = T2
	    for i=1,nX-2 do begin ; don't do the boundary pts
	    	for j=1,nY-2 do begin

                d = points[i,j].par

	    		; Get T along parallel domain 

        		_T  = interpolate ( T2_copy, ( d.x - x[0] ) / (x[-1]-x[0]) * (nX-1.0), ( d.y - y[0] ) / (y[-1]-y[0]) * (nY-1.0), cubic = -0.5 )
                	_Q  = interpolate ( Q, ( d.x - x[0] ) / (x[-1]-x[0]) * (nX-1.0), ( d.y - y[0] ) / (y[-1]-y[0]) * (nY-1.0), cubic = -0.5 )

	    		k = fltArr(n_elements(d.s)) + kPar ; diffusion coefficent
	    		nT = 1
	    		_T = heat1d(d.s,_T,_Q,k,dT/2,nT,cfl=cfl,plot=0) 

	    		T2[i,j] = interpol(_T,d.s,0) ; get T at the actual point

	    	endfor
	    endfor

        c=contour(T2,x,y,/fill,/current,dimensions=[width,height],rgb_table=50)

        ; Solve perp

        T2_copy = T2
	    for i=1,nX-2 do begin ; don't do the boundary pts
	    	for j=1,nY-2 do begin

                d = points[i,j].per

	    		; Get T along parallel domain 

        		_T  = interpolate ( T2_copy, ( d.x - x[0] ) / (x[-1]-x[0]) * (nX-1.0), ( d.y - y[0] ) / (y[-1]-y[0]) * (nY-1.0), cubic = -0.5 )
                	_Q  = interpolate ( Q, ( d.x - x[0] ) / (x[-1]-x[0]) * (nX-1.0), ( d.y - y[0] ) / (y[-1]-y[0]) * (nY-1.0), cubic = -0.5 )

	    		k = fltArr(n_elements(d.s)) + kPer ; diffusion coefficent
	    		nT = 1
	    		_T = heat1d(d.s,_T,_Q,k,dt/2,nT,cfl=cfl,plot=0) 

	    		T2[i,j] = interpol(_T,d.s,0) ; get T at the actual point

	    	endfor
	    endfor

        c.erase
        c=contour(T2,x,y,/fill,/current,dimensions=[width,height],rgb_table=50)

    endfor

	p=plot(TSolution[*],T2[*],symbol="Circle",lineStyle='none')

stop
end
