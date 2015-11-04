function get_bCrossDotGradPsi, x, y

    return, -!dpi^2*cos(!dpi*y)^2*sin(!dpi*x)^2 - !dpi^2*cos(!dpi*x)^2*sin(!dpi*y)^2

end
