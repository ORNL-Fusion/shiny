pro shiny

    resolve_routine, 'interpb', /either, /compile_full_file

    nX = 32 
    nY = 32  

	noPlot = 0		

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

        xMin = -0.5d0 
        xMax = +0.5d0
        x = dIndGen(nX)/(nX-1)*(xMax-xMin)+xMin
        x2D = rebin(x,nX,nY)
        dX = x[1]-x[0]
        
        yMin = -0.5d0 
        yMax = +0.5d0
        y = dIndGen(nY)/(nY-1)*(yMax-yMin)+yMin
        y2D = transpose(rebin(y,nY,nX))
        dY = y[1]-y[0]
    
        psi = getPsi(x2D,y2D)

        ; b = zUnit x grad(psi)

        ; Use analytic expression for b instead

        bx = +!pi*cos(!pi*x2d)*sin(!pi*y2d)
        by = -!pi*cos(!pi*y2d)*sin(!pi*x2d) 
        bz = bx*0

        kPer = 1
        kPar = 1e6

        bMag = sqrt(bx^2+by^2+bz^2) 
        bxU = bx / bMag
        byU = by / bMag
        bzU = bz / bMag

        Q = getQ(x2d,y2d) 

        T = fltArr(nX,nY)
        T2 = T 
        TSolution = tFac(kPer,0) * psi

		if not noPlot then begin
			c=contour(psi,x,y,layout=[2,2,1],/fill,title='psi')
			c=contour(TSolution,x,y,layout=[2,2,2],/current,/fill,title='T')
			c=contour(Q,x,y,layout=[2,2,3],/current,/fill,title='Q')
			v=vector(bx,by,x,y,layout=[2,2,4],/current,auto_color=1,rgb_table=10,auto_subsample=1,title='B')
		endif

    endelse


    CFL = 0.9 ; must be < 1.0 for this shitty explicit forward Euler time differencing
    _D = max(abs([kPer,kPar]))
    dt = CFL * ( 1.0 / 8.0 ) * (dx^2 + dy^2) / _D

	EndTime = 1d0
    nT = ceil(EndTime/dt);1000000L

    width = 400
    height = 400
    oVid = IDLffVideoWrite('asymmetric.webm')
    fps = 4 
    vidStream = oVid.AddVideoStream(width, height, fps)

    xRange=[x[0],x[-1]]
    yRange=[y[0],y[-1]]
    nLevs = 11
    c=contour(T,x,y,/fill,/buffer,dimensions=[width,height],rgb_table=7)

	doGunter = 0

    ; Solve the 2D problem directly on a Cartesian grid
	if doGunter then begin
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
            yp = y[j]+dy/2
            ym = y[j]-dy/2

            symmetric = 0

            if not symmetric then begin

                ; Asymmetric scheme
                ; -----------------

                ; Temperature gradient terms

                dTdx_ipj = (T[i+1,j] - T[i,j])/dx
                dTdy_ipj = (T[i+1,j+1] + T[i,j+1] - T[i,j-1] - T[i+1,j-1]) / (4*dy)
                dTdx_ijp = (T[i+1,j+1] + T[i+1,j] - T[i-1,j+1] - T[i-1,j]) / (4*dx)
                dTdy_ijp = (T[i,j+1] - T[i,j])/dy

                dTdx_imj = (T[i,j] - T[i-1,j])/dx
                dTdy_imj = (T[i,j+1] + T[i-1,j+1] - T[i-1,j-1] - T[i,j-1]) / (4*dy)
                dTdx_ijm = (T[i+1,j] + T[i+1,j-1] - T[i-1,j] - T[i-1,j-1]) / (4*dx)
                dTdy_ijm = (T[i,j] - T[i,j-1])/dy

                ; Heat conduction term
                ; q = -D.\/(T) where D is a 2x2 tensor

                ; x plus / minus terms
                _b = get_b(xp,y[j])
                _b = _b/mag(_b)

                D_ipj = get_D(kPer,kPar,_b)
                q_ipj = -DdotGradT( D_ipj, [dTdx_ipj, dTdy_ipj] )

                _b = get_b(xm,y[j])
                _b = _b/mag(_b)

                D_imj = get_D(kPer,kPar,_b)
                q_imj = -DdotGradT( D_imj, [dTdx_imj, dTdy_imj] )

                ; y plus / minus terms
                _b = get_b(x[i],yp)
                _b = _b/mag(_b)

                D_ijp = get_D(kPer,kPar,_b)
                q_ijp = -DdotGradT( D_ijp, [dTdx_ijp, dTdy_ijp] )

                _b = get_b(x[i],ym)
                _b = _b/mag(_b)

                D_ijm = get_D(kPer,kPar,_b)
                q_ijm = -DdotGradT( D_ijm, [dTdx_ijm, dTdy_ijm] )

                ; Diffusion term
                ; -\/.q

                divq = (q_ipj[0] - q_imj[0])/dx + (q_ijp[1]-q_ijm[1])/dy

            endif else begin

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


                _b = get_b(xp,yp)
                _b = _b/mag(_b)

                D_ipjp = get_D(kPer,kPar,_b)
                q_ipjp = -DdotGradT( D_ipjp, [dTdx_ipjp,dTdy_ipjp] )

                _b = get_b(xp,ym)
                _b = _b/mag(_b)

                D_ipjm = get_D(kPer,kPar,_b)
                q_ipjm = -DdotGradT( D_ipjm, [dTdx_ipjm,dTdy_ipjm] )

                _b = get_b(xm,yp)
                _b = _b/mag(_b)

                D_imjp = get_D(kPer,kPar,_b)
                q_imjp = -DdotGradT( D_imjp, [dTdx_imjp,dTdy_imjp] )

                _b = get_b(xm,ym)
                _b = _b/mag(_b)

                D_imjm = get_D(kPer,kPar,_b)
                q_imjm = -DdotGradT( D_imjm, [dTdx_imjm,dTdy_imjm] )


                divq = (q_ipjp[0] + q_ipjm[0] - q_imjp[0] - q_imjm[0]) / (2*dx) $
                    + (q_ipjp[1] + q_imjp[1] - q_imjm[1] - q_ipjm[1]) / (2*dy)

            endelse

            TUpdate[i,j] = dt * (-divq + Q[i,j])

        endfor
        endfor

        T = T + TUpdate

        ; plot time evolving solution at a subset of times
        if _t mod 50 eq 0 and _t gt 0 then begin
            if _t gt 0 then c.erase
                    print, _t, nT    
                    c=contour(T,x,y,/fill,/buffer,/current,dimensions=[width,height],rgb_table=7)
            frame = c.CopyWindow()
            !null = oVid.put(vidStream, frame)
			time = dt * nT
			save, T, time, _t, dt, nT, fileName='gunter.sav'
        endif

    endfor  

    oVid.cleanup

    ; Compare with analytic solution

    l2 = norm(TSolution-T,lNorm=2)
    dkPer = 1/T[nX/2,nY/2]-kPer
    r=regress(TSolution[*],T[*],yfit=fit)

	if not noPlot then begin
		p=plot(TSolution[*],T[*],symbol="Circle",lineStyle='none')
		p=plot(TSolution[*],fit,/over,$
			title='Slope: '+string(r)+'    L2Norm: '+string(l2),color='y')
	endif
	
	endif ; doGunter

    ; Solve using the 1D set
    ; ----------------------

    oVid2 = IDLffVideoWrite('operator-split.webm')
    fps = 4 
    vidStream2 = oVid2.AddVideoStream(width*2, height, fps)

    bndry_x = [xMin,xMax,xMax,xMin,xMin]
    bndry_y = [yMin,yMin,yMax,yMax,yMin]
    boundary = Obj_New('IDLanROI',bndry_x,bndry_y)

    bndry = {x:bndry_x,y:bndry_y,b:boundary}

    b = { x: x, y: y, bx : bx, by : by, bz: bx*0, psi:psi }

    ; Estimate appropriate 1-D length for a single timestep
    ; i.e., information only flows at the CFL condition, so 
    ; ...

    ; CFL = kPar * dt / dS^2

	nT_im = 1000L 
	dT_im = EndTime/nT_im 

	restorePoints = 0

	if restorePoints then begin

		restore, 'points.sav'

	endif else begin

		nCFL = 15 ; i.e., number of grid points in the 1-D domain
		
		lPar = 5*sqrt(kPar * dt_im / 0.4) 
		lPer = 2.0;nCFL * sqrt(kPer * dt / 0.4)

		n1DTrace = 5000
		dSPar = lPar / n1dTrace
		dSPer = lPer / n1dTrace

		__n = 2*nCFL-1
		d1 = { $
			N: __n, $
			x: fltArr(__n), $
			y: fltArr(__n), $
			dS: dSPar, $
			s: fltArr(__n), $
			kx: fltArr(__n), $
			ky: fltArr(__n), $
			L: lpar }   

		pt = { x: 0.0, y:0.0, par : d1, per : d1 } ; grid point structure contains both per and par domains    

		points = replicate( pt, nX, nY )

		if not noPlot then c=contour(psi,x,y,/fill)

		print, 'Tracing 1-D per & par domains ...'

		for i=1,nX-2 do begin
			for j=1,nY-2 do begin

				points[i,j].x = x[i]
				points[i,j].y = y[j]
				points[i,j].par.L = lPar
				points[i,j].per.L = lPer
				points[i,j].par.dS = dsPar
				points[i,j].per.dS = dsPer
		
				crash=0 
				;if i eq 19 and j eq 14 then crash = 1

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

		save, points, fileName = 'points.sav'

	endelse ; if restorePoints

	print, 'dt_im / dt : ',dT_im / dt

    c=contour(T2,x,y,/fill,/buffer,dimensions=[width*2,height],rgb_table=1,layout=[2,1,1])

	lookAt1DPar = 0
	lookAt1DPer = 0
	cubic = -0.5 ; Setting this to be non-zero causes problems. Do not do it :)
    useAnalyticBCs = 1
	plotMod = 1
	ExponentiallyIncreasingScheme = 1
	nSubPer = 10
	nSubPar = 10
	if ExponentiallyIncreasingScheme then begin
		this_dtPer = dT_im
		this_dtPar = dT_im
	endif else begin
		this_dtPer = dT_im/nSubPer
		this_dtPar = dT_im/nSubPar
	endelse

    for itr=0, nT_im-1 do begin

		tNow = (itr)*dt_im 
		tNext = tNow + dt_im

		TSolution = getTa(x2d,y2d,kPer,tNext)

        TiPer = T2
        TiPar = T2

        TPer = fltArr(nX,nY)
        TPar = fltArr(nX,nY)

        ; Solve perp

        for i=1,nX-2 do begin ; don't do the boundary pts
            for j=1,nY-2 do begin

                d = points[i,j].per

                ; Get T along parallel domain 

				_i = ( d.x - x[0] ) / (x[-1]-x[0]) * (nX-1.0)
				_j = ( d.y - y[0] ) / (y[-1]-y[0]) * (nY-1.0)

                _T  = interpolate ( TiPer, _i, _j, cubic = cubic, /double )
				_Q  = getQ(d.x,d.y)

                k = fltArr(n_elements(d.s)) + kPer ; diffusion coefficent

				if lookAt1DPer and i eq nX/3 and j eq nY/4 then begin
				
					_TaNow = getTa(d.x,d.y,kPer,tNow) 

					_Ti = _T
					_T2 = _T

					_T = heat1d(d.s,_T,_Q,k,this_dtPer,nSubPer,tNow,$
							cfl=cfl,plot=0,CN=1,BT=0,d=d,$
							useAnalyticBCs=useAnalyticBCs,_debug=0,$
							AnalyticBCTime = tNow+dt_im, $
							ExponentiallyIncreasing=ExponentiallyIncreasingScheme)
					;_T2 = heat1d(d.s,_T2,_Q,k,dt,ceil(dt_im/dt),tNow,cfl=cfl,plot=0,CN=0,BT=0,d=d,useAnalyticBCs=useAnalyticBCs) 

					_Ta = getTa(d.x,d.y,kPer,tNext) 
					__Q  = interpolate ( Q, _i, _j, cubic = cubic, /double )

					p=plot(d.s,_Ti,layout=[2,1,1],color='r',thick=3)
					p=plot(d.s,_Ta,/over,color='g',thick=3)
					p=plot(d.s,_TaNow,/over,color='orange',thick=2)

					p=plot(d.s,_T,/over)
					;p=plot(d.s,_T2,/over,color='b')

					p=plot(d.s,_Q,layout=[3,1,3],/current)
					p=plot(d.s,__Q,/over,color='b')
					stop
				endif else begin
					_T = heat1d(d.s,_T,_Q,k,this_dtPer,nSubPer,tNow,$
							cfl=cfl,plot=0,CN=1,BT=0,d=d,$
							useAnalyticBCs=useAnalyticBCs, $ 
							AnalyticBCTime = tNow+dt_im,$
							ExponentiallyIncreasing=ExponentiallyIncreasingScheme) 
				endelse

                ;TPer[i,j] = _T[n_elements(d.s)/2]; get T at the actual point
                TPer[i,j] = interpol(_T,d.s,0) ; get T at the actual point

            endfor
        endfor

        ; Solve parallel 

        for i=1,nX-2 do begin ; don't do the boundary pts
            for j=1,nY-2 do begin

                d = points[i,j].par

                ; Get T along parallel domain 

				_i = ( d.x - x[0] ) / (x[-1]-x[0]) * (nX-1.0)
				_j = ( d.y - y[0] ) / (y[-1]-y[0]) * (nY-1.0)

                _T  = interpolate ( TPer, _i, _j, cubic = cubic )

				_Q  = getQ(d.x,d.y)*0 ; Source is applied only to the perp solve

                k = fltArr(n_elements(d.s)) + kPar ; diffusion coefficent

				if lookAt1DPar and i eq nX/4 and j eq nY/4 then begin
					_Ti = _T
                    _T2 = _T
                    _TiBL  = (bilinear ( TPer, transpose(_i), transpose(_j) ))[*]

					_T = heat1d(d.s,_T,_Q,k,this_dtPar,nSubPar,tNow,$
							cfl=cfl,plot=0,CN=1,BT=0,$
							useAnalyticBCs=useAnalyticBCs,d=d,_debug=0, $
							AnalyticBCTime = tNow+dt_im, $
							ExponentiallyIncreasing=ExponentiallyIncreasingScheme) 
					;_T2 = heat1d(d.s,_T2,_Q,k,dt_im/2000d0,2000,tNow,cfl=cfl,plot=0,CN=0,BT=0,useAnalyticBCs=useAnalyticBCs,d=d,_debug=0) 

					_Ta = getTa(d.x,d.y,kPer,tNext) 
					_Tai = getTa(d.x,d.y,kPer,tNow) 

					__Q  = interpolate ( Q, _i, _j, cubic = cubic, /double ) * 0

					p=plot(d.s,_Ti,layout=[2,1,1],color='r',thick=3)
					;p=plot(d.s,_Tai,/over,color='y',thick=2)

                    p=plot(d.s,_TiBL,/over,color='orange',thick=2)
					p=plot(d.s,_Ta,/over,color='g',thick=3)
					p=plot(d.s,_T,/over)
					;p=plot(d.s,_T2,/over,color='r')
					p=plot(d.s,_Q,layout=[3,1,3],/current)
					p=plot(d.s,__Q,/over,color='b')

					stop

				endif else begin
					_T = heat1d(d.s,_T,_Q,k,this_dtPar,nSubPar,tNow,$
							cfl=cfl,plot=0,CN=1,BT=0,$
							useAnalyticBCs=useAnalyticBCs,d=d,$ 
							AnalyticBCTime = tNow+dt_im, $
							ExponentiallyIncreasing=ExponentiallyIncreasingScheme) 
				endelse

                TPar[i,j] = interpol(_T,d.s,0) ; get T at the actual point

            endfor
        endfor

        T2 = TPar

        ; plot time evolving solution at a subset of times
        if itr mod plotMod eq 0 and itr gt 0 then begin
            if itr gt 0 then c.erase
            print, itr, nT_im    
            c=contour(T2,x,y,/fill,/buffer,/current,dimensions=[width*2,height],rgb_table=1,layout=[2,1,1])
            c=contour(T2-TSolution,x,y,$
				/fill,/buffer,/current,dimensions=[width*2,height],rgb_table=1,layout=[2,1,2])

            frame = c.CopyWindow()
            !null = oVid2.put(vidStream2, frame)
			time = dt_im * nT_im 
			save, T2, time, itr, dt_im, nT_im, fileName='op-split.sav'

			_i = ( 0 - x[0] ) / (x[-1]-x[0]) * (nX-1.0)
			_j = ( 0 - y[0] ) / (y[-1]-y[0]) * (nY-1.0)

			T_00 = interpolate ( T2, _i, _j, cubic = cubic, /double )
			T_00a = getTa(0.0,0.0,kPer,tNext)
			print, 1d0/T_00 - kPer, ((1d0/T_00-kPer)-(1d0/T_00a-kPer))/(1d0/T_00a-kPer)

        endif

    endfor

    ; Compare with analytic solution

    TSolution = getTa(x2d,y2d,kPer,tNext)

    l2 = norm(TSolution-T2,lNorm=2)
    r=regress(TSolution[*],T2[*],yfit=fit)
	if not noPlot then begin
		p=plot(TSolution[*],T2[*],symbol="Circle",lineStyle='none')
		p=plot(TSolution[*],fit,/over,$
			title='Slope: '+string(r)+'    L2Norm: '+string(l2),color='y')
	endif
	print, 'nT_im : ', nT_im
	print, 'nCFL : ', nCFL
	print, 'lPar : ', lPar
	print, 'lPer : ', lPer
	print, 'n1DTrace : ', n1DTrace
	print, 'kPar : ', kPar
	print, 'L2Norm : ', l2
	print, 'OS Time : ', nT_im * dt_im
	print, 'Time for information along parallel domain : ', 1/(kPar / (points[0,0].par.l)^2)
	print, 'Time per parallel iteration : ', dt_im / 2 
stop
end
