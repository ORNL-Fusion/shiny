function dlg_line_polygon, x1, y1, x2, y2, bndry_x, bndry_y 

	px = !null
	py = !null

	for bb=0,n_elements(bndry_x)-2 do begin ; assumes closed polygon bndry
		x3=bndry_x[bb]
		y3=bndry_y[bb]
		x4=bndry_x[bb+1]
		y4=bndry_y[bb+1]
		px = [px, ( (x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4) ) $
				/ ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)  )]
		py = [py, ( (x1*y1-y2*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4)  ) $
				/ ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)  )]
	endfor

	;pick the closest of those intersection points to the interior pt 

	_DistToIntersections = sqrt((x2-px)^2+(y2-py)^2)
	iiPt = where(_DistToIntersections eq min(_DistToIntersections), iiPtCnt)
	if iiPtCnt gt 1 then stop ; bugger

	xLeft = px[iiPt]
	yLeft = py[iiPt]

	return, [xLeft,yLeft]
end
