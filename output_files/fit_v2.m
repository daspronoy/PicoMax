%=========================================================================
%  fit_pseudopotential_Te.m
%-------------------------------------------------------------------------
%  Purpose  :
%    •  Read a discrete V(q) dataset (vectors x, y) for trigonal Te.
%    •  Fit those points to a compact Fourier‑like analytic form that is
%       compatible with the symmetric + antisymmetric decomposition used in
%       empirical pseudopotential and SOC routines.
%    •  Regenerate the reciprocal‑lattice form‑factor table (|G|² , V(G)) in
%       the exact layout your solver expects and echo it to the console so
%       that you can paste it straight into the C++/MATLAB front end that
%       builds H(k).
%-------------------------------------------------------------------------
%  How to use :
%    1)  Put your column vectors  q  ("x")  and  Vq  ("y")  in the workspace
%        *in 2π/a units* and *in the same energy unit you feed to the solver*.
%        Example :  load('Te_Vq_dataset.mat');
%    2)  Run the script.  A non‑linear least‑squares fit (lsqcurvefit) will
%        optimise the Fourier coefficients.  The script prints the updated
%        form‑factor list and draws V(q)_fit on top of the raw data.
%-------------------------------------------------------------------------
%  NOTE      :
%    •  The analytic form keeps the 8th‑order cosine/sine family you already
%       use, which is standard for chiral‑axis semiconductors and meshes
%       cleanly with the λ_S, λ_A SOC terms that live elsewhere in your code.
%    •  If you have a separate dataset for the SOC form factor λ_A(q), copy
%       the *FIT BLOCK* a second time with its own parameter vector.
%=========================================================================

%% ---------------- USER INPUT ------------------------------------------
% 1)  Supply your (q , V(q)) dataset here.  Replace the two example lines.
%     Units : q in  2π/a  ;  V  in eV  (or whatever your solver expects)
%-------------------------------------------------------------------------
% EXAMPLE load:
% load('Te_Vq_dataset.mat');   % -> variables x , y
%-------------------------------------------------------------------------
% QUICK manual entry (delete in production):
% x = [0.000; 0.800; 1.600; 2.400; 3.200; 4.000];   %  2π/a  units
% y = [17.96;  5.02; -4.31; -8.55;  3.68;  0.07];    %  eV
%-------------------------------------------------------------------------

assert(exist('x','var')==1 && exist('y','var')==1, ...
    'Vectors  x  and  y  must exist in the workspace before running.');

x  = x(:);               % ensure column vectors
y  = y(:);

%% ---------------- LATTICE + RECIPROCAL VECTORS ------------------------
u = 0.263170;                % internal parameter  (a units)
f = 1.330207;                % c/a ratio          (a units)

vec1 = [ 1 , -1/sqrt(3) , 0     ];   % (2π/a) primitive reciprocal basis
vec2 = [ 1 ,  1/sqrt(3) , 0     ];
vec3 = [ 0 ,  0         , 1/f   ];
K     = @(h,k,l) h*vec1 + k*vec2 + l*vec3;   % |G(hkl)| in (2π/a)

%% ---------------- ANALYTIC PSEUDOPOTENTIAL FORM -----------------------
%   V(q) = a0 + Σ_{n=1..8} [ a_n cos(n w q) + b_n sin(n w q) ]
%   w    = 0.6636   (kept fixed – gives first zero near  q ≈ π/a )
%   v    = phase‑shift centre; keep v = 0  for Te three‑fold axis.
%-------------------------------------------------------------------------
w         = 0.6636;          % rad / (2π/a)
vshift    = 0;               % phase‑shift (kept zero)
orderN    = 8;               % highest harmonic

% Initial guess  (take your previous values; tweak a0 slightly to match y(1))
p0 = [ 17.95 ,  5.136 , -33.18 , -25.76 , -8.536 , ...
       -8.667 , 16.41 ,  8.545 ,  6.482 , 3.706 , -3.49 , ...
       -1.038 , -1.577 , -0.4536 , 0.1913 , 0.01524 , 0.0677 ];
%               a0      a1      b1      a2      b2
%               a3      b3      a4      b4     a5     b5
%               a6      b6      a7      b7     a8     b8

% Fit function handle ---------------------------------------------------
Vq_model = @(p,q) p(1) + ...
    p(2)*cos((q-vshift)*w)    + p(3)*sin((q-vshift)*w)    + ...
    p(4)*cos(2*(q-vshift)*w)  + p(5)*sin(2*(q-vshift)*w)  + ...
    p(6)*cos(3*(q-vshift)*w)  + p(7)*sin(3*(q-vshift)*w)  + ...
    p(8)*cos(4*(q-vshift)*w)  + p(9)*sin(4*(q-vshift)*w)  + ...
    p(10)*cos(5*(q-vshift)*w) + p(11)*sin(5*(q-vshift)*w) + ...
    p(12)*cos(6*(q-vshift)*w) + p(13)*sin(6*(q-vshift)*w) + ...
    p(14)*cos(7*(q-vshift)*w) + p(15)*sin(7*(q-vshift)*w) + ...
    p(16)*cos(8*(q-vshift)*w) + p(17)*sin(8*(q-vshift)*w);

%% ---------------- FIT BLOCK -------------------------------------------
opts = optimoptions('lsqcurvefit','Display','iter','MaxFunEvals',5e4);
[pfit,~,~,exitflag] = lsqcurvefit(Vq_model , p0 , x , y , [] , [] , opts);
assert(exitflag>0,'Fit did not converge – check initial guesses or data.');

fprintf('\n>>>  Optimised Fourier coefficients  (Te local potential)  <<<\n');
Cnames = {'a0','a1','b1','a2','b2','a3','b3','a4','b4','a5','b5','a6','b6','a7','b7','a8','b8'};
for ii=1:numel(Cnames)
    fprintf('%s = %12.5g\n',Cnames{ii},pfit(ii));
end

%% ---------------- REBUILD |G|² – V(G) TABLE ---------------------------
results = [];
for h = 0:4
    for k = 0:4
        for l = 0:4
            Gvec = K(h,k,l);
            q    = norm(Gvec);            % 2π/a  units
            if q <= 5                     % plane‑wave cutoff
                results = [results ; q^2 , Vq_model(pfit,q) ];
            end
        end
    end
end
results          = sortrows(results);
rounded_results  = round(results , 10);
results          = unique(rounded_results , 'rows');

%  Duplicate V(G) three times as per your solver's (|G|²,V,V,V) layout ----
output = cat(1 , results(2:end,1) , results(2:end,2) , ...
                   results(2:end,2) , results(2:end,2) );

fprintf('\noutput: [%s]\n\n', join(string(output), ','));

%% ---------------- VISUAL CHECK ----------------------------------------
figure(1);  clf;
scatter(x , y , 55 , 'k' , 'filled');  hold on;
qfine = linspace(min(x) , max(x) , 800);
plot(qfine , Vq_model(pfit,qfine) , 'LineWidth',1.4);
set(gca,'FontSize',14,'LineWidth',0.8);
xlabel('|q|   (2\pi/a)'); ylabel('V(q)  [unit]');
legend('data','fit','Location','best'); grid on;
xlim([0 6]);

title('Empirical Pseudopotential fit for Trigonal Te');

%% ---------------- NEXT STEPS  -----------------------------------------
% *  Pass  pfit  to your SOC module :
%      λ_S(q) = λ_S_eV * cos(q·τ)      (symmetric)
%      λ_A(q) = λ_A_eV * sin(q·τ)      (antisymmetric)
%   Keep the same harmonic w and cut‑off so that V(q) & λ(q) share a common
%   plane‑wave basis – this guarantees the correct splitting of valence
%   bands at Γ/L and gives you direct control over the Weyl‑point energy.
%-------------------------------------------------------------------------
