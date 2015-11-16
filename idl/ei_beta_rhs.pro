function EI_beta_RHS, x
	common EI
	return, EI_b * EI_mu * ( x^2-1d0) - (x^EI_M-1d0)
end


