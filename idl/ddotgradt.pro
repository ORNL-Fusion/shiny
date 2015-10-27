function DdotGradT, D, gradT

	q = [0,0]
	q[0] = D[0,0] * gradT[0] + D[1,0] * gradT[1]
	q[1] = D[0,1] * gradT[0] + D[1,1] * gradT[1]

	return, q
	
end
