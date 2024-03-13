









function [V] = pseudopotential(G,v)

Ry2eV = 13.605699169854621; % Ry/eV

vec_tau = [0.125;0.125;0.125];
GMAG = norm(G)^2;

vs = 0.0;
va = 0.0;
if abs(GMAG-3.0)<1e-10
    vs = v(1); va = v(4);
elseif abs(GMAG-4.0)<1e-10
               va = v(5);
elseif abs(GMAG-8.0)<1e-10
    vs = v(2);
elseif abs(GMAG-11.0)<1e-10
    vs = v(3); va = v(6);
end

V = (vs*cos(2*pi*dot(G,vec_tau)) ...
    +1i*va*sin(2*pi*dot(G,vec_tau)) ...
    )*Ry2eV;

V = V*exp(-1i*2*pi*dot(G,vec_tau));


% vec_tau = [0.25;0.25;0.25];
% 
% V = (vs+va) ...
%     +(vs-va)*exp(2i*pi*dot(G,vec_tau));
% V = Ry2eV*V;



