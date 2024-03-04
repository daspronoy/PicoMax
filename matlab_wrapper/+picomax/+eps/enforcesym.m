function [epsmn] = enforcesym(epsmn,Q)



if abs(Q(1)-Q(2))<1e-7 ... [111] direction (G-L)
        && abs(Q(2)-Q(3))<1e-7 ...
        && abs(Q(3)-Q(1))<1e-7
    epsmn = enforcesym_111(epsmn);
elseif abs(Q(2))<1e-7 ... [100] direction (G-X)
        && abs(Q(3))<1e-7
    epsmn = enforcesym_100(epsmn);
else
    epsmn = enforcesym_xxx(epsmn);
end


end

function [epsmn] = enforcesym_100(epsmn)
epsmn([3 4 5]) = epsmn(2); % s01
epsmn([7 8 9]) = epsmn(6); % s02
epsmn([13 23 32]) = epsmn(12); % s16
epsmn(22) = epsmn(14); % s18
epsmn([25 35 45]) = epsmn(15); % s14
epsmn([17 24 27 33 36 43 44]) = epsmn(16); % s12
epsmn([26 34 42]) = epsmn(18); % s15
epsmn([21 31 41]) = epsmn(11); % s11
epsmn([61 71 81]) = epsmn(51); % s22
epsmn([53 63 72]) = epsmn(52); % s24
epsmn(62) = epsmn(54); % s23
end


function [epsmn] = enforcesym_111(epsmn)
epsmn([4 6],1) = epsmn(3,1); % s05
epsmn([7 8],1) = epsmn(5,1); % s01
epsmn([4 6],2) = epsmn(3,2); % s58
epsmn([7 8],2) = epsmn(5,2); % s18
epsmn(4,4) = epsmn(3,3); % s55
epsmn(6,6) = epsmn(3,3);
epsmn(7,7) = epsmn(5,5); % s11
epsmn(8,8) = epsmn(5,5);
epsmn([24 33]) = epsmn(22); % s56
epsmn([25 32 35 52 53]) = epsmn(23); % s16
epsmn([34 42]) = epsmn(26); % s15
epsmn([36 54]) = epsmn(27); % s45
epsmn([44 62]) = epsmn(43); % s12
epsmn([63 72]) = epsmn(45); % s14
end


function [epsmn] = enforcesym_xxx(epsmn)
epsmn(3) = epsmn(2); % s07
epsmn([5 6 7]) = epsmn(4); % s01
epsmn(9) = epsmn(8); % s03
epsmn(21) = epsmn(11); % s77
epsmn([25 23 25]) = epsmn(13); % s17
epsmn([16 22 24]) = epsmn(14); % s18
epsmn(27) = epsmn(17); % s38
epsmn(26) = epsmn(18); % s37
epsmn(52) = epsmn(32); % s16
epsmn(43) = epsmn(33); % s12
epsmn(42) = epsmn(34); % s15
epsmn([45 53 63]) = epsmn(35); % s14
epsmn([44 54 62]) = epsmn(36); % s13
epsmn([41 51 61]) = epsmn(31); % s11
epsmn(81) = epsmn(71); % s33




end




