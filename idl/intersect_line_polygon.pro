; Find intersection point of a line [(x1,y1),(x2,y2)] and a polygon

function intersect_line_polygon, x1, y1, x2, y2, bndry_x, bndry_y, crash=crash

	px = !null
	py = !null

	for bb=0,n_elements(bndry_x)-2 do begin ; assumes closed polygon bndry
		x3=bndry_x[bb]
		y3=bndry_y[bb]
		x4=bndry_x[bb+1]
		y4=bndry_y[bb+1]
        ; https://en.wikipedia.org/wiki/Line-line_intersection
		px = [px, ( (x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4) ) $
				/ ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)  )]
		py = [py, ( (x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4)  ) $
				/ ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)  )]
	endfor

	;pick the closest of those intersection points to the interior pt 

	_DistToIntersections = sqrt((x2-px)^2+(y2-py)^2)

    ; catch for parallel line with polyon edge
    if keyword_set(crash) then stop

    iiGood = where(_distToIntersections eq _distToIntersections, iiGoodCnt)
    if iiGoodCnt ne n_elements(_distToIntersections) then begin
        _DistToIntersections = _DistToIntersections[iiGood]        
        px = px[iiGood]
        py = py[iiGood]
    endif

	iiPt = where(_DistToIntersections eq min(_DistToIntersections), iiPtCnt)
	;if iiPtCnt gt 1 or iiPt[0] lt 0 then 

	xLeft = px[iiPt[0]]
	yLeft = py[iiPt[0]]

    if keyword_set(crash) then stop
	return, [xLeft,yLeft]

end
