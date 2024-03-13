



function [] = dos(o)

NFREQ = o.num.NFREQ;
DFREQ = o.num.DFREQ;
% NBAND = o.num.NBAND;


lambda = o.dat.electron.D;
omega_min = min(lambda(1,:))-5;
omega = DFREQ*(0:1:NFREQ-1)+omega_min;
o.dat.electron.dos = zeros([1 NFREQ]);



for f=1:NFREQ
    o.dat.electron.dos(f) = o.dat.electron.dos(f) ...
        +sum(picomax.util.diracdelta_gauss(omega(f)-lambda,o.num.EPSILON),'all');
end
o.dat.electron.omega = omega;





