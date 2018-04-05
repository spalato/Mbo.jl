export ev2angphz, phz2ev, ev2phz

#=                  UNIT CONVERSION                                     =#
function ev2phz(x)
    eV = 1.6021766208E-19
    h = 6.62607E-34
    femto = 1E-15
    x*eV*femto/h
end

ev2angphz(x) = ev2phz(x)*2*pi

function phz2ev(x)
    eV = 1.6021766208E-19
    h = 6.62607E-34
    femto = 1E-15
    x/eV/femto*h
end
