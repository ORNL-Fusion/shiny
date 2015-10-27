function tfac, kPer, t
    compile_opt IDL2 
    return, (1d0 - exp( -2d0*kPer*!dpi^2*t) ) / kPer

end


