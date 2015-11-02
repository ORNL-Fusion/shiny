function getTa2, args

	x = args.x
	t = args.t
	L = args.L
	D = args.D

	return, sin(!pi*x/L) * exp(-D * !pi^2 * t / L^2 )		
	
end


