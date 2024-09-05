
function GM = GravitationalParameter(Body)

GM = cspice_bodvrd(upper(Body),'GM',1);

end