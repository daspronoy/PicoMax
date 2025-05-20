%primitive lattice (a units)
vecp1=[-1/2, -sqrt(3)/2, 0];
vecp2=[1, 0, 0];
vecp3=[0, 0, 1.33];

%reciprocal lattice
vecp2timesvecp3= [0, 1.33, 0];
vecp3timesvecp1= [0.8660*1.33, -0.5*1.33, 0];
vecp1timesvecp2=[0, 0 ,0.8660];

%reciprocal lattice vectors(2pi/a units)
vec2=[0, 2/sqrt(3), 0];
vec1=[1, -1/sqrt(3), 0] ;
vec3=[0, 0 ,1/1.33];


%tau vectors:
vecT1 = [-0.2632/2, -1.71*0.2632/2, -1.33/3];
vecT2 = [1, 0, 0]*0.2632;
vecT3 = [-0.2632/2, 1.71*0.2632/2, 1.33/3];

a1=    -2.4872;
a2=    -2.1289;
a3=    -1.2941;
a4=    1.2112;
a5=    1.3127;
a6=    -1.3769;
a7=    1.2635;
a8=    0.5433;
b1=    1.1642;
b2=    -1.5589;
b3=    1.6235;
b4=    1.5175;
b5=    1.5216;
b6=    3.4701;
b7=    -3.5268;
b8=    1.3766;
c1=    -4.8124;
c2=    0.1810;
c3=    -0.9587;
c4=    2.7676;
c5=    -10.6904;
c6=    -0.9569;
c7=    4.2242;
c8=    5.2054;
d=     1.3956;
e=     0.8983;




%G-vector:
K = @(h,k,l) h*vec1 + k*vec2 + l*vec3;

polyval = @(x) (a1*sin(b1*x + c1) + a2*sin(b2*x+ c2) + a3*sin(b3*x + c3) + a4*sin(b4*x + c4) + a5*sin(b5*x + c5) + a6*sin(b6*x + c6) + a7*sin(b7*x+ c7) + a8*sin(b8*x + c8)).*exp(-d.*x.^e);

results=[];

for h = 0:4
    for k = 0:4
        for l = 0:4
            if exp(dot(K(h,k,l), vecT1 + vecT2 + vecT3)) > 0
                results =[results;[norm(K(h,k,l))^2,polyval(norm(K(h,k,l)))]];
            end
        end
    end
end
results=sortrows(results);


rounded_results = round(results, 8);
results=unique(rounded_results,"rows");

output=cat(1, results(2:end,1), results(2:end,2), results(2:end,2), results(2:end,2));


fprintf('output: [%s]\n', join(string(output), ','));


a=(a1*sin(b1*x + c1) + a2*sin(b2*x+ c2) + a3*sin(b3*x + c3) + a4*sin(b4*x + c4) + a5*sin(b5*x + c5) + a6*sin(b6*x + c6) + a7*sin(b7*x+ c7) + a8*sin(b8*x + c8)).*exp(-d.*x.^e);
scatter(x,y,'black'); hold on;
plot(x, a);
xlim([0 6]);

