function get_b, x, y

	bx = +!pi*cos(!pi*x)*sin(!pi*y)
	by = -!pi*cos(!pi*y)*sin(!pi*x) 
	bz = bx*0

	return, [[bx],[by],[bz]]

end


