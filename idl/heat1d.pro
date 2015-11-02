function heat1d, s, T, Q, k, dt, nT, tNow, cfl = _cfl, plot = _plot, $
CN=CN, BT=BT, BC=BC, $
d=d, useAnalyticBCs=useAnalyticBCs 

	if keyword_set(CN) then useCN=CN else useCN = 0	
	if keyword_set(BT) then useBT=BT else useBT = 0	

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

	if useCN or useBT then begin
	; Crank-Nicolson or Implicit Euler (BTCS) temporal scheme

		idx = IndGen(nX)
		inr = idx[1:-2]

		A = fltArr(nX,nX)
		B = fltArr(nX)

		if useCN then begin

			alp = k[inr]*dt/(ds^2)
			bb = -alp
			cc = -alp
			aa = 2*(1+alp)

			A[inr,inr+1] = bb 
			A[inr,inr-1] = cc 
			A[inr,inr] = aa 

			A[0,0] = 1
			A[-1,-1] = 1

			;cc = -k[inr]/(2*ds^2)
			;bb = 1/dt + k[inr]/ds^2 
			;aa = cc 

			;A[inr,inr+1] = cc
			;A[inr,inr] = bb 
			;A[inr,inr-1] = aa

		endif else begin

            ; http://www.nada.kth.se/~jjalap/numme/FDheat.pdf         

			cc = -k[inr]/ds^2 
			bb = 1/dt + 2*k[inr]/ds^2
			aa = cc 

			A[inr,inr+1] = cc
			A[inr,inr] = bb 
			A[inr,inr-1] = aa

		endelse

		; Fill matrix

		for _t = 0, nT-1 do begin

			if keyword_set(useAnalyticBCs) then begin
				TLNow = getTa(d.x[0],d.y[0],1,tNow+_t*dt);BC.T_LEnd[1]
				TRNow = getTa(d.x[-1],d.y[-1],1,tNow+_t*dt);BC.T_REnd[1]
				TLNex = getTa(d.x[0],d.y[0],1,tNow+_t*dt+dt);BC.T_LEnd[1]
				TRNex = getTa(d.x[-1],d.y[-1],1,tNow+_t*dt+dt);BC.T_REnd[1]
			endif else begin
				TLNow = T[0]
				TRNow = T[-1]
				TLNex = 0
				TRNex = 0
			endelse

			T[0] = TLNex
			T[-1] = TRNex
	
			if useCN then begin
				dd = alp * T[inr+1] + 2*(1-alp) * T[inr] + alp*T[inr-1] + 2*dt*Q[inr]
				;dd = aa * T[inr-1] + (1/dt + aa + cc) * T[inr] + cc*T[inr+1] + 2*Q[inr]
			endif else begin
				dd = 1/dt * T[inr] + Q[inr]				
			endelse

			; See http://people.sc.fsu.edu/~jpeterson/5-CrankNicolson.pdf
			B[inr] = dd
			B[1] = B[1] + alp[0]*TLNow
			B[-2] = B[-2] + alp[-1]*TRNow
		
			TNew = la_linear_equation(A[1:-2,1:-2],B[1:-2],/double,status=status)
			if status ne 0 then stop
	
			T[inr] = TNew

		endfor

	endif else begin
	; Explicit Euler temporal scheme

		_cfl = dt * max(k) / dS^2
		if _cfl gt 0.5 then stop

		for _t = 0, nT - 1 do begin

			if keyword_set(useAnalyticBCs) then begin
				T[0] = getTa(d.x[0],d.y[0],1,tNow+_t*dt+dt);BC.T_LEnd[1]
				T[-1] = getTa(d.x[-1],d.y[-1],1,tNow+_t*dt+dt);BC.T_REnd[1]
			endif

			T[1:-2] = $
					T[1:-2] $
					+ k[1:-2] * dt / dS^2  * ( T[0:-3] - 2*T[1:-2] + T[2:-1] ) $
					+ dt * Q[1:-2]

		endfor	

	endelse

	return, T

end
