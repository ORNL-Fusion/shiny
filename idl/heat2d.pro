function dr_ds, s, r

	common bfield, b

	_x = (r[0]-min(b.x))/(max(b.x)-min(b.x))*(n_elements(b.x)-1)
	_y = (r[1]-min(b.y))/(max(b.y)-min(b.y))*(n_elements(b.y)-1)

	bx = interpolate(b.bx,_x,_y,cubic=-0.5)	
	by = interpolate(b.by,_x,_y,cubic=-0.5)	

	return, [bx,by] 

end

pro getDomain, p, d, _b, bndry, perp

    bndry_x = bndry.x
    bndry_y = bndry.y
    boundary = bndry.b

	_n = (d.N+1)/2
	; generate the field line

	ThisPoint = [p.x,p.y]

	_line1 = dlg_fieldLineTrace_xy(_b,ThisPoint, dir = -1, dS=d.dS, nS=_n, perp=perp )
	_line2 = dlg_fieldLineTrace_xy(_b,ThisPoint, dir = +1, dS=d.dS, nS=_n, perp=perp )

	fLine_XY = [[reverse(_line1[*,0:-2],2)],[_line2[*,1:-2]]]
    fLine_XY = fLine_XY[0:1,*]
	_s = fIndGen(_n)*d.dS
	s = [reverse(-_s[1:-1]),_s[0:-1]]

	; Check intersection with boundary

	isInsideDomain = boundary->ContainsPoints(fLine_XY[0,*],fLine_XY[1,*])
	iiOutside = where(isInsideDomain eq 0, iiOutsideCnt)
	iiThisPt = _n-1

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
	n_fl = d.N;n_elements(fLine_XY[0,*])

	LeftEnd = [fLine_XY[0,0],fLine_XY[1,0]]	
	RightEnd = [fLine_XY[0,-1],fLine_XY[1,-1]]	


	; Now find actual intersection points with the wall(s)

	ReGrid = 0

	if min(iiUse) ne 0 then begin ; left end needs fixing 

		print, 'Fixing left end'

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

		print, 'Fixing right end'

		ReGrid = 1

		iiOut = max(iiUse)+1
		iiIn = max(iiUse)
		x1 = fLine_XY_all[0,iiOut]
		y1 = fLine_XY_all[1,iiOut]
		x2 = fLine_XY_all[0,iiIn]
		y2 = fLine_XY_all[1,iiIn]

		RightEnd = dlg_line_polygon(x1,y1,x2,y2,bndry_x,bndry_y)

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

	;print, LeftEnd, fLine_XY[0,0],fLine_XY[1,0]
	;print, RightEnd, fLine_XY[0,-1],fLine_XY[1,-1]

    d.s = s
    d.x = fLine_XY_all[0,*]
    d.y = fLine_XY_all[1,*]

	;p=plot(d.x,d.y,/over,color='b')
	;p=plot([1,1]*ThisPoint[0],[1,1]*ThisPoint[1],symbol='o',/sym_filled,/over)

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

	nX = 60
	nY = 61

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
	;scale = max(abs(T))
	;levels = fIndGen(nLevs)/(nLevs-1)*scale
	;colors = reverse(bytScl(levels, top=253)+1)
	;c=contour(T, x, y, aspect_ratio=1.0, $
	;		dimensions=[width,height], $
	;		/fill, c_value=levels, rgb_indices=colors,$
	;		rgb_table=3,yRange=yRange,xRange=xRange,/buffer)
    c=contour(T,x,y,/fill,/buffer,dimensions=[width,height],rgb_table=50)

    frame = c.CopyWindow()
	!null = oVid.put(vidStream, frame)
	
	; Solve the 2D problem directly on a Cartesian grid

	T[0,*] = TSolution[0,*]
        T[-1,*] = TSolution[-1,*]
        T[*,0] = TSolution[*,0]
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
	        ;scale = max(abs(T-mean(T)))
	        ;levels = fIndGen(nLevs)/(nLevs-1)*scale+mean(T)
	        ;colors = (bytScl(levels-mean(T), top=253)+1)
			;c=contour(T, x, y, aspect_ratio=1.0, $
			;		dimensions=[width,height], $
			;		/fill, c_value=levels, rgb_indices=colors,$
			;		rgb_table=3,/current,/buffer)
			;c=contour(-T, x, y, aspect_ratio=1.0, $
			;		dimensions=[width,height], $
			;		/fill, c_value=levels, rgb_indices=colors,$
			;		rgb_table=1,/current,/buffer,/over)
			frame = c.CopyWindow()
			!null = oVid.put(vidStream, frame)
		endif



	endfor	

	oVid.cleanup

	p=plot(TSolution[*],T[*],symbol="Circle",lineStyle='none')


	; Solve using the 1D set

	bndry_x = [xMin,xMax,xMax,xMin,xMin]
	bndry_y = [yMin,yMin,yMax,yMax,yMin]
	boundary = Obj_New('IDLanROI',bndry_x,bndry_y)

    bndry = {x:bndry_x,y:bndry_y,b:boundary}

	common bfield, b
	b = { x: x, y: y, bx : bx, by : by, bz: bx*0, psi:psi }

	_n = 30
	dS = 0.005

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
			getDomain, points[i,j], par, b, bndry, 0
            points[i,j].par = par

            per = points[i,j].per
			getDomain, points[i,j], per, b, bndry, 1
            points[i,j].per = per

            plotMod = 4 
            if i mod plotMod eq 0 and j mod plotMod eq 0 then begin
                    p=plot(points[i,j].per.x,points[i,j].per.y,/over)
                    p=plot(points[i,j].par.x,points[i,j].par.y,/over)
            endif 
		endfor
	endfor

stop
	fl_nX = 10
	fl_nY = 10

	fl_size = size

	fl_xMin = r0-fl_size 
	fl_xMax = r0+fl_size
	fl_x = fIndGen(fl_nX)/(fl_nX-1)*(fl_xMax-fl_xMin)+fl_xMin
	fl_x2D = rebin(fl_x,fl_nX,fl_nY)
	fl_dX = fl_x[1]-fl_x[0]

	fl_yMin = z0-fl_size 
	fl_yMax = z0+fl_size
	fl_y = fIndGen(fl_nY)/(fl_nY-1)*(fl_yMax-fl_yMin)+fl_yMin
	fl_y2D = transpose(rebin(fl_y,fl_nY,fl_nX))
	fl_dY = fl_y[1]-fl_y[0]
	
	__s = fIndGen(_n)*dS
	__s = [reverse(-__s[1:-1]),__s[0:-1]]

    _b = { bR : g.bR, $
		bt : g.bPhi, $
		bz : g.bz, $
		r : g.r, $
		z : g.z, $
		rsize : g.r[-1]-g.r[0], $
		zsize : g.z[-1]-g.z[0], $
		nR : n_elements(g.r), $
		nZ : n_elements(g.z)}   

	fl_T = fltArr(fl_nX,fl_nY)

	for i=1,fl_nX-2 do begin ; don't do the boundary pts
		for j=1,fl_nY-2 do begin

				; Get T along said field line

    		_T  = interpolate ( T, ( fLine_CYL[0,*] - x[0] ) / (xMax-xMin) * (nX-1.0), $
        		( fLine_CYL[2,*] - y[0] ) / (yMax-yMin) * (nY-1.0), cubic = -0.5 )

			k = fltArr(n_elements(s)) + 0.06 ; diffusion coefficent
			nT = 5000
			_T = heat1d(s,_T,k,nT,cfl=0.4,plot=0) 

			fl_T[i,j] = interpol(_T,s,0,/spline) ; get T at the actual point

		endfor
		
	endfor

	flyRange=[fl_yMin,fl_yMax]
	flxRange=[fl_xMin,fl_xMax]
	c=contour(fl_T, fl_x, fl_y, aspect_ratio=1.0, $
			dimensions=[width,height], $
			/fill, c_value=levels, rgb_indices=colors,$
			rgb_table=3,yRange=flyRange,xRange=flxRange)
	c=contour(alog(-g.psizr), g.r, g.z,/over,color='w', c_value=psi_levels)

stop

stop
end
