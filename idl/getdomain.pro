pro getDomain, x0, y0, d, _b, bndry, perp, crash=crash

    bndry_x = bndry.x
    bndry_y = bndry.y
    boundary = bndry.b

	_n = fix(d.L/d.dS)+1
	; generate the field line

	ThisPoint = [x0,y0]

	_line1 = dlg_fieldLineTrace_xy(_b,ThisPoint, dir = -1, dS=d.dS, nS=_n, perp=perp, analytic_b=1 )
	_line2 = dlg_fieldLineTrace_xy(_b,ThisPoint, dir = +1, dS=d.dS, nS=_n, perp=perp, analytic_b=1 )

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
		LeftEnd = intersect_line_polygon(x1,y1,x2,y2,bndry_x,bndry_y,crash=crash)

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
		RightEnd = intersect_line_polygon(x1,y1,x2,y2,bndry_x,bndry_y,crash=crash)

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


