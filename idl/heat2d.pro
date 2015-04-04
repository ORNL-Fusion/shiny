pro heat2d

	nX = 100
	nY = 140

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

	width = 400
	height = 400
	;oVid = IDLffVideoWrite('heat2d.webm')
	;fps = 20 
	;vidStream = oVid.AddVideoStream(width, height, fps)

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
	psi_min = -0.5
	psi_max = -0.49
	psi_n = 30
	psi_levels = [-1,-0.9,-0.8,-0.75,-0.72,-0.71,-0.705]
	c=contour(alog(-g.psizr), g.r, g.z,/over,color='w',c_value = psi_levels)

	; Solve using the 1D set

	fl_nX = 10
	fl_nY = 20

	fl_size = 0.15

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
	
	_n = 100
	dS = 0.1
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

	bndry_x = [xMin,xMax,xMax,xMin,xMin]
	bndry_y = [yMin,yMin,yMax,yMax,yMin]
	boundary = Obj_New('IDLanROI',bndry_x,bndry_y)

	for i=0,fl_nX-1 do begin
		for j=0,fl_nY-1 do begin

			s = __s

			print, i, j

			; generate the field line

			ThisPoint = [fl_x[i],0,fl_y[j]]

			_line1 = dlg_fieldLineTrace(_b,ThisPoint,$
					dir = -1, dS=dS, nS=_n )
			_line2 = dlg_fieldLineTrace(_b,ThisPoint,$
					dir = +1, dS=dS, nS=_n )

			fLine_CYL = [[reverse(_line1[*,0:-2],2)],[_line2[*,1:-2]]]

			; Check intersection with boundary

			isInsideDomain = boundary->ContainsPoints(fLine_CYL[0,*],fLine_CYL[2,*])
			iiOutside = where(isInsideDomain eq 0, iiOutsideCnt)
			iiThisPt = _n-1

			; Extract that part of the line within the boundary
			; accounting for craziness, multiple intersections, etc

			nF = n_elements(fLine_CYL[0,*])

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
			fLine_CYL_all = fLine_CYL
			fLine_CYL = fLine_CYL[*,iiUse]
			s_all = s
			s = s[iiUse]
			sLeft = s[0]
			sRight = s[-1]
			n_fl = n_elements(fLine_CYL[0,*])

			LeftEnd = [fLine_CYL[0,0],fLine_CYL[2,0]]	
			RightEnd = [fLine_CYL[0,-1],fLine_CYL[2,-1]]	

			p=plot(fLine_CYL_all[0,*],fLine_CYL_all[2,*],/over,thick=2)
			p=plot(fLine_CYL[0,*],fLine_CYL[2,*],/over,color='b')

			; Now find actual intersection points with the wall(s)

			ReGrid = 0

			if min(iiUse) ne 0 then begin ; left end needs fixing 

				print, 'Fixing left end'

				ReGrid = 1

				iiOut = min(iiUse)-1
				iiIn = min(iiUse)
				x1 = fLine_CYL_all[0,iiOut]
				y1 = fLine_CYL_all[2,iiOut]
				x2 = fLine_CYL_all[0,iiIn]
				y2 = fLine_CYL_all[2,iiIn]

				LeftEnd = dlg_line_polygon(x1,y1,x2,y2,bndry_x,bndry_y)

				t2 = fLine_CYL_all[1,iiIn]
				; this is a bit of a hack, but it's pretty close
				tLeft = interpol(fLine_CYL_all[1,iiOut:iiIn],[x1,x2],LeftEnd[0])

				_x1 = LeftEnd[0]*cos(tLeft)
				_y1 = LeftEnd[0]*sin(tLeft)
				_z1 = LeftEnd[1]

				_x2 = x2*cos(t2)
				_y2 = x2*sin(t2)
				_z2 = y2

				sExtraBit = sqrt ( (_x2-_x1)^2+(_y2-_y1)^2+(_z2-_z1)^2 ) 

				sLeft = s_all[iiIn] - sExtraBit

			endif 

			if max(iiUse) ne n_elements(fLine_CYL_all[0,*])-1 then begin ; right end needs fixing 

				print, 'Fixing right end'

				ReGrid = 1

				iiOut = max(iiUse)+1
				iiIn = max(iiUse)
				x1 = fLine_CYL_all[0,iiOut]
				y1 = fLine_CYL_all[2,iiOut]
				x2 = fLine_CYL_all[0,iiIn]
				y2 = fLine_CYL_all[2,iiIn]

				RightEnd = dlg_line_polygon(x1,y1,x2,y2,bndry_x,bndry_y)

				t2 = fLine_CYL_all[1,iiIn]
				tRight = interpol(fLine_CYL_all[1,iiIn:iiOut],[x2,x1],RightEnd[0])

				_x1 = RightEnd[0]*cos(tRight)
				_y1 = RightEnd[0]*sin(tRight)
				_z1 = RightEnd[1]

				_x2 = x2*cos(t2)
				_y2 = x2*sin(t2)
				_z2 = y2

				sExtraBit = sqrt ( (_x2-_x1)^2+(_y2-_y1)^2+(_z2-_z1)^2 ) 

				sRight = s_all[iiIn] + sExtraBit

			endif 

			; Interpolate field line to new constant dS such that it ends at the boundary(s)

			if ReGrid then begin

				sNew = fIndGen(n_fl)/(n_fl-1)*(sRight-sLeft)+sLeft

				_fLine_r = interpol(fLine_CYL_all[0,*],s_all,sNew,/spline)
				_fLine_z = interpol(fLine_CYL_all[2,*],s_all,sNew,/spline)

				s = sNew
				fLine_CYL[0,*] = _fLine_r
				fLine_CYL[2,*] = _fLine_z

			endif

			print, LeftEnd, fLine_CYL[0,0],fLine_CYL[2,0]
			print, RightEnd, fLine_CYL[0,-1],fLine_CYL[2,-1]

			; Get T along said field line

    		_T  = interpolate ( T, ( fLine_CYL[0,*] - x[0] ) / (xMax-xMin) * (nX-1.0), $
        		( fLine_CYL[2,*] - y[0] ) / (yMax-yMin) * (nY-1.0), cubic = -0.5 )

			k = fltArr(n_elements(s)) + 0.06 ; diffusion coefficent
			nT = 35000
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
