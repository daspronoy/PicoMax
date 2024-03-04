function [k,f] = fcc(point)
% k = [kx;ky;kz] (cartesian coordinate in 2*pi/a units)
% f = [f1;f2;f3] (fractional coordinate in b units)
switch point
    case 'G'
        k = [0;0;0];
    case 'X'
        k = [1;0;0];
    case 'L'
        k = [0.5;0.5;0.5];
    case 'K'
        k = [0.75;0.75;0];
    case 'W'
        k = [1;0.5;0];
    case 'U'
        k = [1;0.25;0.25];
    otherwise
        error('check highsymmetry point (fcc)')
end
if nargout<2; return; end

switch point
    case 'G'
        f = [0;0;0];
    case 'X'
        f = [0;0.5;0.5];
    case 'L'
        f = [0.5;0.5;0.5];
    case 'K'
        f = [3/8;3/8;3/4];
    case 'W'
        f = [0.25;0.5;0.75];
    case 'U'
        f = [0.25;5/8;5/8];
    otherwise
        error('check highsymmetry point (fcc)')
end

end