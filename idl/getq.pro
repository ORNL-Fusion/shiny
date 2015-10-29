function getQ, x, y
	lap = -2*!pi^2*getPsi(x,y) ; analytic Laplacian(psi)
	return, -lap
end


