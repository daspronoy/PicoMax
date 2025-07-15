%===================================================================
%  derive_lambda_SA_Te.m
%-------------------------------------------------------------------
%  Purpose
%   ▸  Take the latest sin(bx)/(bx)^n analytic fit for the local Te
%      pseudopotential, V(q), and derive the reciprocal‑space spin–orbit
%      form‑factor amplitudes λ_S and λ_A (in Rydberg).
%   ▸  λ_S multiplies  Σ cos(G·τ_i);   λ_A multiplies Σ sin(G·τ_i)  in the
%      SOC term  V_SO(G) = λ_S cos_sum + i λ_A sin_sum.
%   ▸  The script
%        1)  Evaluates λ(q) = C⋅(1/q) dV/dq  (C = 0.025, from your code)
%        2)  Builds a linear system  A·λ = b  over a shell of |G| ≤ G_max
%        3)  Solves for the constants λ_S, λ_A  (least‑squares) and prints
%           them both in eV and in Ry.
%-------------------------------------------------------------------
%  Usage
%    •  Ensure the coefficient set  a_0…a_2, b_0…b_2  in the PARAMETER
%       BLOCK reflects your most recent fit.
%    •  Run the script – it prints λ_S_eV, λ_A_eV and their Ry equivalents.
%    •  Copy the Ry numbers into your C++ routine:
%          lambda_S_eV =  ...;   // actually in Ry now
%          lambda_A_eV =  ...;
%===================================================================
%% ---------------- LATTICE & BASIS -------------------------------------
u = 0.263170;           % internal trigonal parameter (a units)
f = 1.330207;           % c/a ratio

vecp1 = [-1/2, -sqrt(3)/2, 0];
vecp2 = [ 1  ,      0     , 0];
vecp3 = [ 0  ,      0     , f];

% reciprocal primitive vectors (2π/a units)
vec1 = [ 1 , -1/sqrt(3) , 0   ];
vec2 = [ 1 ,  1/sqrt(3) , 0   ];
vec3 = [ 0 ,  0         , 1/f ];
K     = @(h,k,l) h*vec1 + k*vec2 + l*vec3;     % |G(hkl)|   (2π/a)

% atomic positions τ_i  (in units of the direct a‑lattice)
vecT1 = [ 0               ,  0.5*sqrt(3)*u , -1/3*f];
vecT2 = [ 1.5*u           ,  0             ,  0    ];
vecT3 = [ 0               , -0.5*sqrt(3)*u ,  1/3*f];
Tlist = [vecT1; vecT2; vecT3];

%% ---------------- PSEUDOPOTENTIAL PARAMETERS --------------------------
% Coefficients from your latest sin(bx)/(bx)^n fit  (energy in eV)
% a_0 =  0.108;
% a_1 = -0.5043;
% a_2 = -0.1753;
% b_0 =  0.5751;
% b_1 = -1.011 ;
% b_2 =  1.412 ;
% 
% f_sinc  = @(B,X)  sin(B.*X)./(B.*X);                                       % sinc(bx)
% df_sinc = @(B,X)  B.*((B.*X).*cos(B.*X) - sin(B.*X))./((B.*X).^2);         % d/dx
% 
% Vq       = @(x) a_0*f_sinc(b_0,x)           + a_1*(f_sinc(b_1,x)).^2 + ...
%                       a_2*(f_sinc(b_2,x)).^3;
% dVdq     = @(x) a_0*df_sinc(b_0,x)          + a_1*2*f_sinc(b_1,x).*df_sinc(b_1,x) + ...
%                       a_2*3*(f_sinc(b_2,x)).^2.*df_sinc(b_2,x);




w      = 0.6636;                 % (2π/a)-units
vshift = 0;                      % keep zero for Te three-fold axis
orderN = 8;

% parameter set [a0  a1 b1  a2 b2  ...  a8 b8]  (units: eV)
p = [ 17.95 ,  5.136 , -33.18 , -25.76 , -8.536 , ...
      -8.667 , 16.41 ,  8.545 ,   6.482 ,  3.706 , ...
      -3.49  , -1.038 , -1.577 , -0.4536, 0.1913 , ...
       0.01524 , 0.0677 ];

% convenience handles ------------------------------------------------
Vq = @(q) fourier_V(p,q,w,vshift,orderN);        % value in eV
dVdq = @(q) d_fourier_V(p,q,w,vshift,orderN);    % derivative in eV

lambda_q_eV = @(q) 0.025*q.*dVdq(q);          % SOC prefactor (eV)
E2RY        = 1/13.605693009;                      % eV → Ry conversion

%% ---------------- BUILD LINEAR SYSTEM  A·λ = b  -----------------------
G_max   = 8;                 %  (2π/a) units – same as your earlier loop
A   = [];
bRy = [];

for h = 0:4
    for k = 0:4
        for l = 0:4
            Gvec   = K(h,k,l);
            q      = norm(Gvec);
            if q < 1e-6 || q > G_max,  continue;  end   % skip G=0 and beyond cutoff

            % λ(q) from analytic derivative  (convert to Ry)
            lam_Ry = lambda_q_eV(q) * E2RY;

            % structure‑factor sums (dimensionless)
            phase  = Gvec * Tlist.';           % 1×3 dot products
            cos_sum = sum(cos(phase));
            sin_sum = sum(sin(phase));

            % accumulate row into least‑squares matrix
            A   = [A ; cos_sum , sin_sum];
            bRy = [bRy ; lam_Ry];
        end
    end
end


%% ---------------- BUILD FULL λ_S(G) TABLE -----------------------------
%  For high‑accuracy runs we feed the *exact* λ_S(q) into the Hamiltonian.
%  Format requested :  G : λ_S  ,  G : λ_S  ,  ...               [units: eV]

lambda_pairs = compose('%.7g:%.7g', results(2:end,1) , ...
                                          lambda_q_eV(results(2:end,1)));
fprintf('%s', strjoin(lambda_pairs,','));
fprintf('\n');
















function V = fourier_V(p,q,w,v,orderN)
    % p = [a0  a1 b1  a2 b2 ...], q vectorised
    V = p(1);
    for n = 1:orderN
        an = p(2*n);      bn = p(2*n+1);
        V  = V + an*cos(n*w*(q-v)) + bn*sin(n*w*(q-v));
    end
end

function dV = d_fourier_V(p,q,w,v,orderN)
    dV = 0;
    for n = 1:orderN
        an = p(2*n);      bn = p(2*n+1);
        dV = dV + (-an)*n*w .* sin(n*w*(q-v)) + bn*n*w .* cos(n*w*(q-v));
    end
end