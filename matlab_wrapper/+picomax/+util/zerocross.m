function [fcross] = zerocross(y)



fpos = y .* (y >= 0);
% fneg = y .* (y <= 0);
fcross = find(diff(fpos>0));