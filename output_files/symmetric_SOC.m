%primitive lattice (a units)
u=0.263170;
f=1.330207;


vecp1=[-1/2, -sqrt(3)/2, 0];
vecp2=[1, 0, 0];
vecp3=[0, 0, f];

%reciprocal lattice
vecp2timesvecp3= [0, f, 0];
vecp3timesvecp1= [0.8660*f, -0.5*f, 0];
vecp1timesvecp2=[0, 0 ,0.8660];

%reciprocal lattice vectors(2pi/a units)
vec1=[1, -1/sqrt(3), 0];
vec2=[1, 1/sqrt(3), 0];
vec3=[0, 0, 1/f];


%tau vectors:
vecT1 = [0, 0.5*sqrt(3)*u, -1.0/3.0*f];
vecT2 = [1.5, 0, 0]*u;
vecT3 = [0, -0.5*sqrt(3)*u, 1.0/3.0*f];

a_0 =       0.108;
a_1 =     -0.5043;
a_2 =     -0.1753;
b_0 =      0.5751;
b_1 =      -1.011;
b_2 =       1.412;



%G-vector:
K = @(h,k,l) h*vec1 + k*vec2 + l*vec3;

f      = @(B,X)  sin(B.*X)./(B.*X);                                   %  sin(bx)/(bx)
dfdx   = @(B,X)  B.*((B.*X).*cos(B.*X) - sin(B.*X))./((B.*X).^2);     %  d/dx of above

% ---------- NEW polyval and its derivative ------------------------------
polyval     = @(x)  a_0*f(b_0,x)                    + ...
                    a_1*(f(b_1,x)).^2               + ...
                    a_2*(f(b_2,x)).^3;

der_polyval = @(x)  a_0*dfdx(b_0,x)                 + ...
                    a_1*2*f(b_1,x).*dfdx(b_1,x)     + ...
                    a_2*3*(f(b_2,x)).^2.*dfdx(b_2,x);

results=[];

%---------------- loop over G-shells --------------------------------------
results = [];

for h = 0:4
    for k = 0:4
        for l = 0:4
            Gnorm = norm(K(h,k,l));
            if Gnorm <= 8
                results = [results ; Gnorm^2 , 0.025*(1/Gnorm)*der_polyval(Gnorm)];
            end
        end
    end
end


results=sortrows(results);
plot(results(2:end,1),results(2:end,2)); hold on;

rounded_results = round(results, 10);
results=unique(rounded_results,"rows");



a =   -0.0569*0.025;
b =     -0.2292;
c =      0.2347;
bess = @(x) a*(sin(b*x)./((b*x).^2)-cos(c*x)./(c*x));

plot(results(2:end,1),bess(results(2:end,1))); hold off;

rows = 2:size(results,1);              % use 1:size(...) if you want all
pairs = compose('%.10g:%.10g', ...
                results(rows,1), bess(results(rows,1)));   % cell array of strings

fprintf('output: [%s]\n', strjoin(pairs, ','));

