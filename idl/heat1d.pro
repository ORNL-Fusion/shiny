function heat1d, s, T, Q, k, dt, nT, cfl = _cfl, plot = _plot

	dS = s[1]-s[0]

	;if keyword_set(_cfl) then cfl = _cfl else cfl = 10.0 ; must be < 0.5 for this shitty explicit forward Euler time differencing
	;dt = cfl * dS^2 / min(k) 

	_cfl = dt * max(k) / dS^2
	if _cfl gt 0.5 then stop

	if keyword_set(_plot) then p=plot(s,T)

	; Apply BCs

	;T[0] = (-2*T[1] + 0.5*T[2])/(-1.5) ; dT/dS = 0 second order accurate forward difference 
	;T[-1] = (+2*T[-2] - 0.5*T[-3])/(+1.5) ; dT/dS = 0 second order accurate backward difference

    T[0] = 0.0 
	T[-1] = 0.0

	for _t = 0, nT - 1 do begin

		if keyword_set(_plot) then begin			
			; plot time evolving solution at a subset of times
			if _t mod (nT/10) eq 0 then begin
				p=plot(s,T,/over)
			endif
		endif

		T[1:-2] = $
				T[1:-2] $
				+ k[1:-2] * dt / dS^2  * ( T[0:-3] - 2*T[1:-2] + T[2:-1] ) $
				+ dt * Q[1:-2]

	endfor	

	return, T

end
