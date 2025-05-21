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

% a0 =       17.95;
a1 =       5.136;
b1 =      -33.18;
a2 =      -25.76;
b2 =      -8.536;
a3 =      -8.667;
b3 =       16.41;
a4 =       8.545;
b4 =       6.482;
a5 =       3.706;
b5 =       -3.49;
a6 =      -1.038;
b6 =      -1.577;
a7 =     -0.4536;
b7 =      0.1913;
a8 =     0.01524;
b8 =      0.0677;
% w =      0.6636;
w=0.6636;
v=0;
a0=17.955;



%G-vector:
K = @(h,k,l) h*vec1 + k*vec2 + l*vec3;

polyval = @(x) a0 + a1*cos((x-v)*w) + b1*sin((x-v)*w) + ...
               a2*cos(2*(x-v)*w) + b2*sin(2*(x-v)*w) + a3*cos(3*(x-v)*w) + b3*sin(3*(x-v)*w) + ...
               a4*cos(4*(x-v)*w) + b4*sin(4*(x-v)*w) + a5*cos(5*(x-v)*w) + b5*sin(5*(x-v)*w) + ...
               a6*cos(6*(x-v)*w) + b6*sin(6*(x-v)*w) + a7*cos(7*(x-v)*w) + b7*sin(7*(x-v)*w) + ...
               a8*cos(8*(x-v)*w) + b8*sin(8*(x-v)*w);

results=[];

for h = 0:4
    for k = 0:4
        for l = 0:4
            if norm(K(h,k,l)) <= 5
                results =[results;[norm(K(h,k,l))^2,polyval(norm(K(h,k,l)))]];
            end
        end
    end
end
results=sortrows(results);


rounded_results = round(results, 10);
results=unique(rounded_results,"rows");

output=cat(1, results(2:end,1), results(2:end,2), results(2:end,2), results(2:end,2));


fprintf('output: [%s]\n', join(string(output), ','));

a=a0 + a1*cos((x-v)*w) + b1*sin((x-v)*w) + ...
               a2*cos(2*(x-v)*w) + b2*sin(2*(x-v)*w) + a3*cos(3*(x-v)*w) + b3*sin(3*(x-v)*w) + ...
               a4*cos(4*(x-v)*w) + b4*sin(4*(x-v)*w) + a5*cos(5*(x-v)*w) + b5*sin(5*(x-v)*w) + ...
               a6*cos(6*(x-v)*w) + b6*sin(6*(x-v)*w) + a7*cos(7*(x-v)*w) + b7*sin(7*(x-v)*w) + ...
               a8*cos(8*(x-v)*w) + b8*sin(8*(x-v)*w);


scatter(x,y,'black'); hold on;
plot(x, a);
grid on;
xlim([0 6]);

