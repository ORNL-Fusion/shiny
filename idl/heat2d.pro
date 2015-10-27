function DdotGradT, D, gradT

	q = [0,0]
	q[0] = D[0,0] * gradT[0] + D[1,0] * gradT[1]
	q[1] = D[0,1] * gradT[0] + D[1,1] * gradT[1]

	return, q
	
end


function D, kPer, kPar, bu 
	b1 = bu[0]
	b2 = bu[1]
	return, [$
		[ kPar*b1^2 + kPer*b2^2, 	(kPar-kPer)*b1*b2 ],$
		[ (kPar-kPer)*b1*b2, 		kPer*b1^2 + kPar*b2^2 ]] 
end


function b, x, y

	bx = +!pi*cos(!pi*x)*sin(!pi*y)
	by = -!pi*cos(!pi*y)*sin(!pi*x) 
	bz = bx*0

	return, [[bx],[by],[bz]]

end

function tFac, kPer, t
    
    return, (1d0 - exp( -2d0*kPer*!dpi^2*t) ) / kPer

end

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

	_n = fix(d.L/d.dS)+1
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

	ReGrid = 1

	if min(iiUse) ne 0 then begin ; left end needs fixing 

		ReGrid = 1

		iiOut = min(iiUse)-1
		iiIn = min(iiUse)
		x1 = fLine_XY_all[0,iiOut]
		y1 = fLine_XY_all[1,iiOut]
		x2 = fLine_XY_all[0,iiIn]
		y2 = fLine_XY_all[1,iiIn]

		if crash then stop
		LeftEnd = dlg_line_polygon(x1,y1,x2,y2,bndry_x,bndry_y,crash=crash)

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

		if crash then stop
		RightEnd = dlg_line_polygon(x1,y1,x2,y2,bndry_x,bndry_y,crash=crash)

		sExtraBit = sqrt ( (x2-RightEnd[0])^2+(y2-RightEnd[1])^2 ) 
		sRight = s_all[iiIn] + sExtraBit

	endif 

	; Interpolate field line to new constant dS such that it ends at the boundary(s)

	sNew = fIndGen(n_fl)/(n_fl-1)*(sRight-sLeft)+sLeft

	_fLine_x = interpol(fLine_XY_all[0,*],s_all,sNew)
	_fLine_y = interpol(fLine_XY_all[1,*],s_all,sNew)

    d.s = sNew
    d.x = _fLine_x
    d.y = _fLine_y

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

	nX = 41 
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
	
	yMin = -0.5 
	yMax = +0.5
	y = fIndGen(nY)/(nY-1)*(yMax-yMin)+yMin
	y2D = transpose(rebin(y,nY,nX))
	dY = y[1]-y[0]
	
        psi = cos(!pi*x2D) * cos(!pi*y2D)
        grad, psi, x, y, gradX, gradY
        ;lap = laplacian( psi, x, y)
	lap = -2*!pi^2*psi ; analytic Laplacian(psi)

        ; b = zUnit x grad(psi)

        bx = -gradY
        by = +gradX
        bz = bx*0

	; Use analytic expression for b instead

	bx = +!pi*cos(!pi*x2d)*sin(!pi*y2d)
        by = -!pi*cos(!pi*y2d)*sin(!pi*x2d) 
        bz = bx*0

        ; Generate kx / ky from kPer / kPar

        kPer = 1
        kPar = 1e3

        bMag = sqrt(bx^2+by^2+bz^2) 
        bxU = bx / bMag
        byU = by / bMag
        bzU = bz / bMag

       	Q = -lap 

        T = fltArr(nX,nY)
        T2 = T 
        TSolution = tFac(kPer,0) * psi

        c=contour(psi,x,y,layout=[2,2,1],/fill,title='psi')
        c=contour(TSolution,x,y,layout=[2,2,2],/current,/fill,title='T')
        c=contour(Q,x,y,layout=[2,2,3],/current,/fill,title='Q')
        v=vector(bx,by,x,y,layout=[2,2,4],/current,auto_color=1,rgb_table=10,auto_subsample=1,title='B')

    endelse


	CFL = 0.3 ; must be < 0.5 for this shitty explicit forward Euler time differencing
    	_D = max(abs([kPer,kPar]))
    	dt = CFL * ( 1.0 / 8.0 ) * (dx^2 + dy^2) / _D

	nT = 5000L

	width = 400
	height = 400
	oVid = IDLffVideoWrite('asymmetric.webm')
	fps = 4 
	vidStream = oVid.AddVideoStream(width, height, fps)

	xRange=[x[0],x[-1]]
	yRange=[y[0],y[-1]]
	nLevs = 11
    	c=contour(T,x,y,/fill,/buffer,dimensions=[width,height],rgb_table=50)

	; Solve the 2D problem directly on a Cartesian grid

	for _t = 0, nT - 1 do begin

		TSolution = tFac(kPer,_t*dt) * psi

		; Try the scheme(s) of Gunter et al. 

		TUpdate = fltArr(nX,nY)*0

		for i=1,nX-2 do begin
		for j=1,nY-2 do begin

			; NOTE : 
			; ipj -> [i+1/2,j]
			; ijp -> [i,j+1/2]
			; imj -> [i-1/2,j]
			; ijm -> [i,j-1/2]

			xp = x[i]+dx/2
			xm = x[i]-dx/2
			yp = y[i]+dy/2
			ym = y[i]-dy/2

			;; Asymmetric scheme
			;; -----------------

			;; Temperature gradient terms

			;dTdx_ipj = (T[i+1,j] - T[i,j])/dx
			;dTdy_ipj = (T[i+1,j+1] + T[i,j+1] - T[i,j-1] - T[i+1,j-1]) / (4*dy)
			;dTdx_ijp = (T[i+1,j+1] + T[i+1,j] - T[i-1,j+1] - T[i-1,j]) / (4*dx)
			;dTdy_ijp = (T[i,j+1] - T[i,j])/dy

			;dTdx_imj = (T[i,j] - T[i-1,j])/dx
			;dTdy_imj = (T[i,j+1] + T[i-1,j+1] - T[i-1,j-1] - T[i,j-1]) / (4*dy)
			;dTdx_ijm = (T[i+1,j] + T[i+1,j-1] - T[i-1,j] - T[i-1,j-1]) / (4*dx)
			;dTdy_ijm = (T[i,j] - T[i,j-1])/dy

			;; Heat conduction term
			;; q = -D.\/(T) where D is a 2x2 tensor

			;; x plus / minus terms
			;_b = b(xp,y[j])
			;_b = _b/sqrt(_b[0]^2+_b[1]^2+_b[2]^2)
			;b1 = _b[0]
			;b2 = _b[1]

			;D_ipj = D(kPer,kPar,b1,b2)
			;q_ipj = [0,0]
			;q_ipj[0] = -( D_ipj[0,0] * dTdx_ipj + D_ipj[1,0] * dTdy_ipj )
			;q_ipj[1] = -( D_ipj[0,1] * dTdx_ipj + D_ipj[1,1] * dTdy_ipj )

			;_b = b(xm,y[j])
			;_b = _b/sqrt(_b[0]^2+_b[1]^2+_b[2]^2)
			;b1 = _b[0]
			;b2 = _b[1]

			;D_imj = D(kPer,kPar,b1,b2)
			;q_imj = [0,0]
			;q_imj[0] = -( D_imj[0,0] * dTdx_imj + D_imj[1,0] * dTdy_imj )
			;q_imj[1] = -( D_imj[0,1] * dTdx_imj + D_imj[1,1] * dTdy_imj )

			;; y plus / minus terms
			;_b = b(x[i],yp)
			;_b = _b/sqrt(_b[0]^2+_b[1]^2+_b[2]^2)
			;b1 = _b[0]
			;b2 = _b[1]

			;D_ijp = D(kPer,kPar,b1,b2)
			;q_ijp = [0,0]
			;q_ijp[0] = -( D_ijp[0,0] * dTdx_ijp + D_ijp[1,0] * dTdy_ijp )
			;q_ijp[1] = -( D_ijp[0,1] * dTdx_ijp + D_ijp[1,1] * dTdy_ijp )

			;_b = b(x[i],ym)
			;_b = _b/sqrt(_b[0]^2+_b[1]^2+_b[2]^2)
			;b1 = _b[0]
			;b2 = _b[1]

			;D_ijm = D(kPer,kPar,b1,b2)
			;q_ijm = [0,0]
			;q_ijm[0] = -( D_ijm[0,0] * dTdx_ijm + D_ijm[1,0] * dTdy_ijm )
			;q_ijm[1] = -( D_ijm[0,1] * dTdx_ijm + D_ijm[1,1] * dTdy_ijm )

			;; Diffusion term
			;; -\/.q

			;divq = (q_ipj[0] - q_imj[0])/dx + (q_ijp[1]-q_ijm[1])/dy


			; Symmetric scheme	
			; ----------------

			dTdx_ipjp = (T[i+1,j+1] + T[i+1,j] - T[i,j+1] - T[i,j]) / (2*dx)	
			dTdy_ipjp = (T[i,j+1] + T[i+1,j+1] - T[i+1,j] -T[i,j]) / (2*dy)

			dTdx_ipjm = (T[i+1,j] + T[i+1,j-1] - T[i,j] - T[i,j-1]) / (2*dx)	
			dTdy_ipjm = (T[i,j] + T[i+1,j] - T[i+1,j-1] -T[i,j-1]) / (2*dy)

			dTdx_imjp = (T[i,j+1] + T[i,j] - T[i-1,j+1] - T[i-1,j]) / (2*dx)	
			dTdy_imjp = (T[i-1,j+1] + T[i,j+1] - T[i,j] -T[i-1,j]) / (2*dy)

			dTdx_imjm = (T[i,j] + T[i,j-1] - T[i-1,j] - T[i-1,j-1]) / (2*dx)	
			dTdy_imjm = (T[i-1,j] + T[i,j] - T[i,j-1] -T[i-1,j-1]) / (2*dy)


			_b = b(xp,yp)
			_b = _b/sqrt(_b[0]^2+_b[1]^2+_b[2]^2)

			D_ipjp = D(kPer,kPar,_b)
			q_ipjp = -DdotGradT( D_ipjp, [dTdx_ipjp,dTdy_ipjp] )

			_b = b(xp,ym)
			_b = _b/sqrt(_b[0]^2+_b[1]^2+_b[2]^2)

			D_ipjm = D(kPer,kPar,_b)
			q_ipjm = -DdotGradT( D_ipjm, [dTdx_ipjm,dTdy_ipjm] )

			_b = b(xm,yp)
			_b = _b/sqrt(_b[0]^2+_b[1]^2+_b[2]^2)

			D_imjp = D(kPer,kPar,_b)
			q_imjp = -DdotGradT( D_imjp, [dTdx_imjp,dTdy_imjp] )

			_b = b(xm,ym)
			_b = _b/sqrt(_b[0]^2+_b[1]^2+_b[2]^2)

			D_imjm = D(kPer,kPar,_b)
			q_imjm = -DdotGradT( D_imjm, [dTdx_imjp,dTdy_imjm] )


			divq = (q_ipjp[0] + q_ipjm[0] - q_imjp[0] - q_imjm[0]) / (2*dx) $
				+ (q_ipjp[1] + q_imjp[1] - q_imjm[1] - q_ipjm[1]) / (2*dy)


			TUpdate[i,j] = dt * (-divq + Q[i,j])

		endfor
		endfor

		T = T + TUpdate

		; plot time evolving solution at a subset of times
		if _t mod 50 eq 0 and _t gt 0 then begin
			if _t gt 0 then c.erase
            		print, _t, nT    
            		c=contour(T,x,y,/fill,/buffer,/current,dimensions=[width,height],rgb_table=50)
			frame = c.CopyWindow()
			!null = oVid.put(vidStream, frame)
		endif

	endfor	

	oVid.cleanup

	l2 = norm(TSolution-T,lNorm=2)
	dkPer = 1/T[nX/2,nY/2]-kPer
	p=plot(TSolution[*],T[*],symbol="Circle",lineStyle='none',aspect_ratio=1.0)
	r=regress(TSolution[*],T[*],yfit=fit)
	p=plot(TSolution[*],fit,/over,$
		title='Slope: '+string(r)+'    L2Norm: '+string(l2),color='y')


stop
	; Solve using the 1D set
    	; ----------------------

	oVid2 = IDLffVideoWrite('operator-split.webm')
	fps = 4 
	vidStream2 = oVid2.AddVideoStream(width, height, fps)

	bndry_x = [xMin,xMax,xMax,xMin,xMin]
	bndry_y = [yMin,yMin,yMax,yMax,yMin]
	boundary = Obj_New('IDLanROI',bndry_x,bndry_y)

    	bndry = {x:bndry_x,y:bndry_y,b:boundary}

	common bfield, b
	b = { x: x, y: y, bx : bx, by : by, bz: bx*0, psi:psi }

	; Estimate appropriate 1-D length for a single timestep
	; i.e., information only flows at the CFL condition, so 
	; ...

	; CFL = kPar * dt / dS^2

	nCFL = 10 ; i.e., number of grid points in the 1-D domain
	
	lPar = nCFL * sqrt(kPar * dt / 0.4)	

	dS = 0.005

	__n = 2*nCFL-1
	d1 = { $
		N: __n, $
		x: fltArr(__n), $
		y: fltArr(__n), $
		dS: dS, $
		s: fltArr(__n), $
		kx: fltArr(__n), $
		ky: fltArr(__n), $
		L: lpar }	

	pt = { x: 0.0, y:0.0, par : d1, per : d1 } ; grid point structure contains both per and par domains	   

	points = replicate( pt, nX, nY )

    	c=contour(psi,x,y,/fill)

	for i=1,nX-2 do begin
		for j=1,nY-2 do begin
			points[i,j].x = x[i]
			points[i,j].y = y[j]
 			crash=0 
             		if i eq 19 and j eq 14 then crash = 1

    			par = points[i,j].par
			getDomain, x[i],y[j], par, b, bndry, 0, crash=crash
            		points[i,j].par = par
            
            		per = points[i,j].per
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

    	nItr = nT 

	c=contour(T2,x,y,/fill,/buffer,dimensions=[width,height],rgb_table=51)

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

		; plot time evolving solution at a subset of times
		if itr mod 50 eq 0 and itr gt 0 then begin
			if itr gt 0 then c.erase
            		print, itr, nItr    
            		c=contour(T2,x,y,/fill,/buffer,/current,dimensions=[width,height],rgb_table=51)
			frame = c.CopyWindow()
			!null = oVid2.put(vidStream2, frame)
		endif
    endfor

	TSolution = tFac(kPer,nItr*dt) * psi

	l2 = norm(TSolution-T2,lNorm=2)
	p=plot(TSolution[*],T2[*],symbol="Circle",lineStyle='none',aspect_ratio=1.0)
	r=regress(TSolution[*],T2[*],yfit=fit)
	p=plot(TSolution[*],fit,/over,$
		title='Slope: '+string(r)+'    L2Norm: '+string(l2),color='y')

stop
end
