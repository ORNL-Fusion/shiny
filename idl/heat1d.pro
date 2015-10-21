function heat1d, s, T, k, nT, cfl = _cfl, plot = _plot

	;nPts = 100
	;sMin = 0
	;sMax = 1
	;s = fIndGen(nPts)/(nPts-1)*(sMax-sMin)+sMin
	dS = s[1]-s[0]

	;T = cos(2*s)>0 ; initial condition
	F = T*0
	;k = fltArr(nPts) + 0.06 ; diffusion coefficent

	if keyword_set(_cfl) then cfl = _cfl else cfl = 10.0 ; must be < 0.5 for this shitty explicit forward Euler time differencing

	; seems to work fine for large CFL numbers here, not sure why exactly

	dt = min( cfl * dS^2 / k)

	;nT = 50000L

	if keyword_set(_plot) then p=plot(s,T)

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
				+ dt * F[1:-2]

		; Apply BCs

		T[0] = (-2*T[1] + 0.5*T[2])/(-1.5) ; dT/dS = 0 second order accurate forward difference 
		T[-1] = (+2*T[-2] - 0.5*T[-3])/(+1.5) ; dT/dS = 0 second order accurate backward difference

	endfor	

	return, T

end
