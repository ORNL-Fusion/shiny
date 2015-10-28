function heat1d, s, T, Q, k, dt, nT, cfl = _cfl, plot = _plot, CN=CN

	if keyword_set(CN) then useCN=CN else useCN = 0	

	nX = n_elements(s)

	; Sanity checking
	if n_elements(k) ne nX then stop
	if n_elements(Q) ne nX then stop
	if n_elements(T) ne nX then stop
	iiBadK = where(k lt 0, iiBadKCnt)
	if iiBadKCnt gt 0 then stop
	
	dS = s[1]-s[0]

	if keyword_set(_plot) then p=plot(s,T)

	; Apply BCs

	;T[0] = (-2*T[1] + 0.5*T[2])/(-1.5) ; dT/dS = 0 second order accurate forward difference 
	;T[-1] = (+2*T[-2] - 0.5*T[-3])/(+1.5) ; dT/dS = 0 second order accurate backward difference

    T[0] = 0.0 
	T[-1] = 0.0

	if useCN then begin
	; Implicit Crank-Nicolson temporal scheme

		idx = IndGen(nX)
		inr = idx[1:-2]

		A = fltArr(nX,nX)
		B = fltArr(nX)

		;bb = -k[inr]/(2*dS^2) ; super diagonal	
		;cc = bb ; sub diagonal
		;aa = 1/dt - (bb+cc) ; main diagonal

		alp = k[inr]*dt/(ds^2)
		bb = -alp
		cc = -alp
		aa = 2*(1+alp)

		; Fill matrix
		A[inr,inr+1] = bb
		A[inr,inr-1] = cc
		A[inr,inr] = aa

		for _t = 0, nT-1 do begin

			;dd = -cc * T[inr-1] + (1/dt + bb + cc) * T[inr] - bb*T[inr+1]
			;B[inr] = dd

			dd = alp * T[inr+1] + 2*(1-alp) * T[inr] + alp*T[inr-1]
			B[inr] = dd

			TNew = la_linear_equation(A[1:-2,1:-2],B[1:-2],/double,status=status)
			if status ne 0 then stop
	
			T[inr] = TNew

		endfor

	endif else begin
	; Explicit Euler temporal scheme

		_cfl = dt * max(k) / dS^2
		if _cfl gt 0.5 then stop

		for _t = 0, nT - 1 do begin

			T[1:-2] = $
					T[1:-2] $
					+ k[1:-2] * dt / dS^2  * ( T[0:-3] - 2*T[1:-2] + T[2:-1] ) $
					+ dt * Q[1:-2]

		endfor	

	endelse

	return, T

end
