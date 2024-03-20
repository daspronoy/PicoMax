function [epsmn] = enforcesym(epsmn,Q)
tol = 1e-7;
if abs(Q(1)-Q(2))<tol ... [111] direction (G-L)
        && abs(Q(2)-Q(3))<tol ...
        && abs(Q(3)-Q(1))<tol
    epsmn = enforcesym_111(epsmn);
elseif abs(Q(2))<tol ... [100] direction (G-X)
        && abs(Q(3))<tol
    epsmn = enforcesym_100(epsmn);
elseif abs(Q(1)-Q(2))<tol ... [110] direction (G-K)
        && abs(Q(3))<tol
    epsmn = enforcesym_110(epsmn);
else
    epsmn = enforcesym_110(epsmn);
end
end




% function [epsmn] = enforcesym_100(epsmn)
% epsmn([2 3 4 5]) = mean(epsmn([2])); % s01
% epsmn([6 7 8 9]) = mean(epsmn([6])); % s05
% epsmn([11 21 31 41]) = mean(epsmn([11])); % s11
% epsmn([12 13 23 32]) = mean(epsmn([12])); % s12
% epsmn([14 22]) = mean(epsmn([14])); % s14
% epsmn([15 25 35 45]) = mean(epsmn([15])); % s15
% epsmn([16 17 24 27 33 36 43 44]) = mean(epsmn([16])); % s16
% epsmn([18 26 34 42]) = mean(epsmn([18])); % s18
% epsmn([51 61 71 81]) = mean(epsmn([51])); % s55
% epsmn([52 53 63 72]) = mean(epsmn([52])); % s56
% epsmn([54 62]) = mean(epsmn([54])); % s58
% end
% function [epsmn] = enforcesym_111(epsmn)
% epsmn([3 4 6]) = mean(epsmn([3])); % s02
% epsmn([5 7 8]) = mean(epsmn([5])); % s04
% epsmn([12 13 15]) = mean(epsmn([12])); % s12
% epsmn([14 16 17]) = mean(epsmn([14])); % s14
% epsmn([21 31 51]) = mean(epsmn([21])); % s22
% epsmn([22 24 33]) = mean(epsmn([22])); % s23
% epsmn([23 25 32 35 52 53]) = mean(epsmn([23])); % s24
% epsmn([26 34 42]) = mean(epsmn([26])); % s27
% epsmn([27 36 54]) = mean(epsmn([27])); % s28
% epsmn([41 61 71]) = mean(epsmn([41])); % s44
% epsmn([43 44 62]) = mean(epsmn([43])); % s46
% epsmn([45 63 72]) = mean(epsmn([45])); % s48
% end
% function [epsmn] = enforcesym_110(epsmn)
% epsmn([2 3]) = mean(epsmn([2])); % s01
% epsmn([4 5 6 7]) = mean(epsmn([4])); % s03
% epsmn([8 9]) = mean(epsmn([8])); % s07
% epsmn([11 21]) = mean(epsmn([11])); % s11
% epsmn([13 15 23 25]) = mean(epsmn([13])); % s13
% epsmn([14 16 22 24]) = mean(epsmn([14])); % s14
% epsmn([17 27]) = mean(epsmn([17])); % s17
% epsmn([18 26]) = mean(epsmn([18])); % s18
% epsmn([31 41 51 61]) = mean(epsmn([31])); % s33
% epsmn([32 52]) = mean(epsmn([32])); % s34
% epsmn([33 43]) = mean(epsmn([33])); % s35
% epsmn([34 42]) = mean(epsmn([34])); % s36
% epsmn([35 45 53 63]) = mean(epsmn([35])); % s37
% epsmn([36 44 54 62]) = mean(epsmn([36])); % s38
% epsmn([71 81]) = mean(epsmn([71])); % s77
% end
% function [epsmn] = enforcesym_Rz(epsmn)
% epsmn([2 3]) = mean(epsmn([2])); % s01
% epsmn([4 5]) = mean(epsmn([4])); % s03
% epsmn([6 7]) = mean(epsmn([6])); % s05
% epsmn([8 9]) = mean(epsmn([8])); % s07
% epsmn([11 21]) = mean(epsmn([11])); % s11
% epsmn([13 23]) = mean(epsmn([13])); % s13
% epsmn([14 22]) = mean(epsmn([14])); % s14
% epsmn([15 25]) = mean(epsmn([15])); % s15
% epsmn([16 24]) = mean(epsmn([16])); % s16
% epsmn([17 27]) = mean(epsmn([17])); % s17
% epsmn([18 26]) = mean(epsmn([18])); % s18
% epsmn([31 41]) = mean(epsmn([31])); % s33
% epsmn([33 43]) = mean(epsmn([33])); % s35
% epsmn([34 42]) = mean(epsmn([34])); % s36
% epsmn([35 45]) = mean(epsmn([35])); % s37
% epsmn([36 44]) = mean(epsmn([36])); % s38
% epsmn([51 61]) = mean(epsmn([51])); % s55
% epsmn([53 63]) = mean(epsmn([53])); % s57
% epsmn([54 62]) = mean(epsmn([54])); % s58
% epsmn([71 81]) = mean(epsmn([71])); % s77
% end


function [epsmn] = enforcesym_100(epsmn)
epsmn([2 3 4 5]) = mean(epsmn([2 3 4 5])); % s01
epsmn([6 7 8 9]) = mean(epsmn([6 7 8 9])); % s05
epsmn([11 21 31 41]) = mean(epsmn([11 21 31 41])); % s11
epsmn([12 13 23 32]) = mean(epsmn([12 13 23 32])); % s12
epsmn([14 22]) = mean(epsmn([14 22])); % s14
epsmn([15 25 35 45]) = mean(epsmn([15 25 35 45])); % s15
epsmn([16 17 24 27 33 36 43 44]) = mean(epsmn([16 17 24 27 33 36 43 44])); % s16
epsmn([18 26 34 42]) = mean(epsmn([18 26 34 42])); % s18
epsmn([51 61 71 81]) = mean(epsmn([51 61 71 81])); % s55
epsmn([52 53 63 72]) = mean(epsmn([52 53 63 72])); % s56
epsmn([54 62]) = mean(epsmn([54 62])); % s58
end
function [epsmn] = enforcesym_111(epsmn)
epsmn([3 4 6]) = mean(epsmn([3 4 6])); % s02
epsmn([5 7 8]) = mean(epsmn([5 7 8])); % s04
epsmn([12 13 15]) = mean(epsmn([12 13 15])); % s12
epsmn([14 16 17]) = mean(epsmn([14 16 17])); % s14
epsmn([21 31 51]) = mean(epsmn([21 31 51])); % s22
epsmn([22 24 33]) = mean(epsmn([22 24 33])); % s23
epsmn([23 25 32 35 52 53]) = mean(epsmn([23 25 32 35 52 53])); % s24
epsmn([26 34 42]) = mean(epsmn([26 34 42])); % s27
epsmn([27 36 54]) = mean(epsmn([27 36 54])); % s28
epsmn([41 61 71]) = mean(epsmn([41 61 71])); % s44
epsmn([43 44 62]) = mean(epsmn([43 44 62])); % s46
epsmn([45 63 72]) = mean(epsmn([45 63 72])); % s48
end
function [epsmn] = enforcesym_110(epsmn)
epsmn([2 3]) = mean(epsmn([2 3])); % s01
epsmn([4 5 6 7]) = mean(epsmn([4 5 6 7])); % s03
epsmn([8 9]) = mean(epsmn([8 9])); % s07
epsmn([11 21]) = mean(epsmn([11 21])); % s11
epsmn([13 15 23 25]) = mean(epsmn([13 15 23 25])); % s13
epsmn([14 16 22 24]) = mean(epsmn([14 16 22 24])); % s14
epsmn([17 27]) = mean(epsmn([17 27])); % s17
epsmn([18 26]) = mean(epsmn([18 26])); % s18
epsmn([31 41 51 61]) = mean(epsmn([31 41 51 61])); % s33
epsmn([32 52]) = mean(epsmn([32 52])); % s34
epsmn([33 43]) = mean(epsmn([33 43])); % s35
epsmn([34 42]) = mean(epsmn([34 42])); % s36
epsmn([35 45 53 63]) = mean(epsmn([35 45 53 63])); % s37
epsmn([36 44 54 62]) = mean(epsmn([36 44 54 62])); % s38
epsmn([71 81]) = mean(epsmn([71 81])); % s77
end
function [epsmn] = enforcesym_Rz(epsmn)
epsmn([2 3]) = mean(epsmn([2 3])); % s01
epsmn([4 5]) = mean(epsmn([4 5])); % s03
epsmn([6 7]) = mean(epsmn([6 7])); % s05
epsmn([8 9]) = mean(epsmn([8 9])); % s07
epsmn([11 21]) = mean(epsmn([11 21])); % s11
epsmn([13 23]) = mean(epsmn([13 23])); % s13
epsmn([14 22]) = mean(epsmn([14 22])); % s14
epsmn([15 25]) = mean(epsmn([15 25])); % s15
epsmn([16 24]) = mean(epsmn([16 24])); % s16
epsmn([17 27]) = mean(epsmn([17 27])); % s17
epsmn([18 26]) = mean(epsmn([18 26])); % s18
epsmn([31 41]) = mean(epsmn([31 41])); % s33
epsmn([33 43]) = mean(epsmn([33 43])); % s35
epsmn([34 42]) = mean(epsmn([34 42])); % s36
epsmn([35 45]) = mean(epsmn([35 45])); % s37
epsmn([36 44]) = mean(epsmn([36 44])); % s38
epsmn([51 61]) = mean(epsmn([51 61])); % s55
epsmn([53 63]) = mean(epsmn([53 63])); % s57
epsmn([54 62]) = mean(epsmn([54 62])); % s58
epsmn([71 81]) = mean(epsmn([71 81])); % s77
end
