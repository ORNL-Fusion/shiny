function get_D, kPer, kPar, bu 
	b1 = bu[0]
	b2 = bu[1]
	return, [$
		[ kPar*b1^2 + kPer*b2^2, 	(kPar-kPer)*b1*b2 ],$
		[ (kPar-kPer)*b1*b2, 		kPer*b1^2 + kPar*b2^2 ]] 
end


