


% a : lattice constnat [angstrom]

function [H] = HamiltonianEPM(G,K,v,a)


KCON = 3.809982111485962; % hbar^2/(2*m_e*eV*angstrom^2)
KINETIC_CONST = KCON*(4*pi*pi)/(a*a);
NPW = size(G,2);

H = zeros(NPW,NPW);
for i=1:NPW 
    for j=1:i 
        if i==j
            H(i,j) = KINETIC_CONST*norm(G(:,i)+K)^2;
        else
            dG = G(:,i)-G(:,j);
            if norm(dG)<3.32
                H(i,j) = picomax.electron.pseudopotential(dG,v);
                H(j,i) = conj(H(i,j));
            end
        end
    end
end








