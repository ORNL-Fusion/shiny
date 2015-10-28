pro shiny

    resolve_routine, 'interpb', /either, /compile_full_file

    nX = 40 
    nY = 40

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
    
        psi = cos(!dpi*x2D) * cos(!dpi*y2D)
        grad, psi, x, y, gradX, gradY
        ;lap = laplacian( psi, x, y)
        lap = -2*!pi^2*psi ; analytic Laplacian(psi)

        ; b = zUnit x grad(psi)

        ;bx = -gradY
        ;by = +gradX
        ;bz = bx*0

        ; Use analytic expression for b instead

        bx = +!pi*cos(!pi*x2d)*sin(!pi*y2d)
        by = -!pi*cos(!pi*y2d)*sin(!pi*x2d) 
        bz = bx*0

        kPer = 1
        kPar = 1e12

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


    CFL = 0.9 ; must be < 1.0 for this shitty explicit forward Euler time differencing
    _D = max(abs([kPer,kPar]))
    dt = CFL * ( 1.0 / 8.0 ) * (dx^2 + dy^2) / _D

    nT = 100L

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
                q_imjm = -DdotGradT( D_imjm, [dTdx_imjp,dTdy_imjm] )


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
                    c=contour(T,x,y,/fill,/buffer,/current,dimensions=[width,height],rgb_table=50)
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
    p=plot(TSolution[*],T[*],symbol="Circle",lineStyle='none')
    r=regress(TSolution[*],T[*],yfit=fit)
    p=plot(TSolution[*],fit,/over,$
        title='Slope: '+string(r)+'    L2Norm: '+string(l2),color='y')


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

    nCFL = 10 ; i.e., number of grid points in the 1-D domain
    
    lPar = nCFL * sqrt(kPar * dt / 0.4) 
    lPer = nCFL * sqrt(kPer * dt / 0.4) 

	n1DTrace = 300
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

    c=contour(psi,x,y,/fill)

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


    nItr = nT 

    c=contour(T2,x,y,/fill,/buffer,dimensions=[width*2,height],rgb_table=51,layout=[2,1,1])

    for itr=0, nItr-1 do begin

    	TSolution = tFac(kPer,nItr*dt) * psi

		TSolution[0,*]=0
		TSolution[-1,*]=0
		TSolution[*,0]=0
		TSolution[*,-1]=0

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
                _T = heat1d(d.s,_T,_Q,k,dT/2,nT,cfl=cfl,plot=0,CN=1) 

                T2[i,j] = interpol(_T,d.s,0,/spline) ; get T at the actual point

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

                T2[i,j] = interpol(_T,d.s,0,/spline) ; get T at the actual point

            endfor
        endfor

        ; plot time evolving solution at a subset of times
        if itr mod 50 eq 0 and itr gt 0 then begin
            if itr gt 0 then c.erase
                    print, itr, nItr    
                    c=contour(T2,x,y,/fill,/buffer,/current,dimensions=[width*2,height],rgb_table=51,layout=[2,1,1])
                    c=contour(T2[1:-2,1:-2]/TSolution[1:-2,1:-2],x[1:-2],y[1:-2],$
						/fill,/buffer,/current,dimensions=[width*2,height],rgb_table=51,layout=[2,1,2])

            frame = c.CopyWindow()
            !null = oVid2.put(vidStream2, frame)
			time = dt * itr
			save, T2, time, itr, dt, nItr, fileName='op-split.sav'
        endif
    endfor


    ; Compare with analytic solution

    TSolution = tFac(kPer,nItr*dt) * psi

	TSolution[0,*]=0
	TSolution[-1,*]=0
	TSolution[*,0]=0
	TSolution[*,-1]=0

    l2 = norm(TSolution-T2,lNorm=2)
    p=plot(TSolution[*],T2[*],symbol="Circle",lineStyle='none')
    r=regress(TSolution[*],T2[*],yfit=fit)
    p=plot(TSolution[*],fit,/over,$
        title='Slope: '+string(r)+'    L2Norm: '+string(l2),color='y')


stop
end
