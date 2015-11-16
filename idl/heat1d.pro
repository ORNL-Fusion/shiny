function heat1d, s, T, Q, k, dt, nT, tNow, cfl = _cfl, plot = _plot, $
CN=CN, BT=BT, $
d=d, useAnalyticBCs=useAnalyticBCs, _debug=debug, $
analyticBCTime=analyticBCTime, $
ExponentiallyIncreasing=ExponentiallyIncreasing 

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

	dTArr = FltArr(nT)+dT
	tNowArr = tNow + IntArr(nT)*dT 
	tNexArr = tNowArr + dT

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

		; Fill matrix

		if keyword_set(ExponentiallyIncreasing) then begin

			; See "Five Ways of Reducing the Crank–Nicolson Oscillations"
			; by Osterby (2003) that describes how the problems of oscillations
			; stem from the initial condition, and that you can take several small
			; initial timesteps to damp them, then ramp up to your desired timestep.
			; The value of bMu below should be <0.5 for several timesteps to start.

			; https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&cad=rja&uact=8&ved=0CDkQFjADahUKEwif9pLAppXJAhVI2T4KHZTWCuw&url=http%3A%2F%2Flink.springer.com%2Farticle%2F10.1023%252FB%253ABITN.0000009942.00540.94&usg=AFQjCNEYUfFh7y1sYGLloQ3a_s9h8ELPuA&sig2=Jhg-ABcFBF5cPtTSImh57g
			
			mu = dt / dS^2
			bMu = max(k * mu) 

			common EI, EI_M, EI_b, EI_mu
			EI_betaGuess = 2d0
			EI_M = nT
			EI_b = max(k)
			EI_mu = mu
			EI_beta = newton(EI_betaGuess,'ei_beta_rhs',/double,stepmax=10,check=local)
			if local gt 0 then stop
			EI_dt = DblArr(EI_M)
			tmp = 0
			k1 = dt*(EI_beta-1d0)/(EI_beta^EI_M-1d0)
			for ee = 0,EI_M-1 do begin
				tmp = tmp + EI_beta^ee*k1
				EI_dt[ee] = tmp
			endfor

			dTArr = EI_dt

		endif

		for _t = 0, nT-1 do begin

			if useCN then begin

        	    ; http://people.sc.fsu.edu/~jpeterson/5-CrankNicolson.pdf

				alp = k[inr]*dtArr[_t]/(ds^2)
				bb = -alp
				cc = -alp
				aa = 2*(1+alp)

				A[inr,inr+1] = bb 
				A[inr,inr-1] = cc 
				A[inr,inr] = aa 

			endif else begin

        	    ; http://www.nada.kth.se/~jjalap/numme/FDheat.pdf         

				cc = -k[inr]/ds^2 
				bb = 1/dtArr[_t] + 2*k[inr]/ds^2
				aa = cc 

				A[inr,inr+1] = cc
				A[inr,inr] = bb 
				A[inr,inr-1] = aa

			endelse

			if keyword_set(useAnalyticBCs) then begin
				TLNex = getTa(d.x[0],d.y[0],1,tNexArr[_t])
				TRNex = getTa(d.x[-1],d.y[-1],1,tNexArr[_t])
				if keyword_set(AnalyticBCTime) then begin
					TLNex = getTa(d.x[0],d.y[0],1,AnalyticBCTime)
					TRNex = getTa(d.x[-1],d.y[-1],1,AnalyticBCTime)
				endif
			endif else begin
				TLNex = T[0]
				TRNex = T[-1]
			endelse

            if keyword_set(debug) then stop	
	
			if useCN then begin
				dd = alp * T[inr+1] + 2*(1-alp) * T[inr] + alp*T[inr-1] + 2*dtArr[_t]*Q[inr]
			endif else begin
				dd = 1/dtArr[_t] * T[inr] + Q[inr]				
			endelse

			; See http://people.sc.fsu.edu/~jpeterson/5-CrankNicolson.pdf
			B[inr] = dd

			B[1] = B[1] + alp[0]*TLNex
			B[-2] = B[-2] + alp[-1]*TRNex
	
			TNew = la_linear_equation(A[1:-2,1:-2],B[1:-2],/double,status=status)
			if status ne 0 then stop
            if keyword_set(debug) then stop	

			T[inr] = TNew
            T[0] = TLNex
            T[-1] = TRNex
		endfor

	endif else begin
	; Explicit Euler temporal scheme

		_cfl = dt * max(k) / dS^2
		if _cfl gt 0.5 then begin
            print, 'dt too large, CFL violated, reduce dt'
            print, 'CFL: ',_CFL
            stop
        endif

		for _t = 0, nT - 1 do begin

			if keyword_set(useAnalyticBCs) then begin
				T[0] = getTa(d.x[0],d.y[0],1,tNexArr[_t])
				T[-1] = getTa(d.x[-1],d.y[-1],1,tNexArr[_t])
				if keyword_set(AnalyticBCTime) then begin		
					T[0] = getTa(d.x[0],d.y[0],1,AnalyticBCTime)
					T[-1] = getTa(d.x[-1],d.y[-1],1,AnalyticBCTime)
				endif
			endif

			T[1:-2] = $
					T[1:-2] $
					+ k[1:-2] * dt / dS^2  * ( T[0:-3] - 2*T[1:-2] + T[2:-1] ) $
					+ dt * Q[1:-2]

		endfor	

	endelse

	return, T

end

