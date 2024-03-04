







function [] = permittivity(o)
o.tmp.PICOMAX_OPT = sprintf( ...
    "-switch eps " + ...
    "-GCUT %d -neps %d " + ...
    "-nband %d " + ...
    "-a %f -epm %s " + ...
    "-nfreq %d -dfreq %f -kk %d -epsilon %f ", ...
    o.num.GCUT,o.num.NEPS, ...
    o.num.NBAND, ...
    o.sys.a,regexprep(num2str(o.sys.epm),'\s+',','), ...
    o.num.NFREQ,o.num.DFREQ,o.num.KK,o.num.EPSILON);

% BZ integration (KPOINT)
if o.num.KPOINT_SCHEME==0 % Cohen-Chadi grid
    o.tmp.PICOMAX_OPT = sprintf( ...
        "%s -gridscheme %d -gridorder %d ", ...
        o.tmp.PICOMAX_OPT, ...
        o.num.KPOINT_SCHEME,o.num.KPOINT_ORDER);
elseif o.num.KPOINT_SCHEME==1 % Monkhorst-Pack grid
    o.tmp.PICOMAX_OPT = sprintf( ...
        "%s -gridscheme %d -kgridfile %s ", ...
        o.tmp.PICOMAX_OPT, ...
        o.num.KPOINT_SCHEME,o.num.KPOINT_FILE);
end

% q vectors
if (~isempty(o.num.QNUM) && ~isempty(o.num.QPATH))
    o.tmp.PICOMAX_OPT = ...
        sprintf("%s -numq %s -pathq %s", ...
        o.tmp.PICOMAX_OPT, ...
        regexprep(num2str(o.num.QNUM),'\s+',','), ...
        o.num.QPATH);
else
    o.tmp.PICOMAX_OPT = ...
        sprintf("%s -q %s", ...
        o.tmp.PICOMAX_OPT, ...
        regexprep(num2str(o.num.QVEC(:).'),'\s+',','));
end

o.tmp.SYSTEM_CMD = sprintf( ...
    "powershell -command ""wsl bash -c '%s %s'""", ...
    o.num.PICOMAX_CMD,o.tmp.PICOMAX_OPT);

tic
[status,o.tmp.CMDOUT] = system(o.tmp.SYSTEM_CMD);
toc

fid = fopen(o.num.PATH_DAT,'r');
% fid = fopen('\\wsl.localhost\ubuntu-22.04\tmp\tmpepsilon.dat','r');
eps = fread(fid,'double');
fclose(fid);


% tmp = splitlines(o.CMDOUT);
% tmp(1:o.NumHeaderLines) = [];
% tmp(end) = [];
% eps = cell2mat(cellfun(@str2num,tmp,'uniform',0));
NFREQ = o.num.NFREQ;

o.dat.omega = o.num.DFREQ*(0:1:NFREQ-1);
if (~isempty(o.num.QNUM) && ~isempty(o.num.QPATH))
    nq = sum(o.num.QNUM)+1;
    o.dat.qi = 0:nq-1;
else
    nq = 1;
    o.dat.qi = 0;
end


ne = o.num.NEPS*(o.num.NEPS+1)/2;

% N = ne*nq;
% o.dat.realepsL = permute(reshape(eps(0*N+1:1*N,:),[ne nq o.nfreq]),[3 2 1]);
% o.dat.imagepsL = permute(reshape(eps(1*N+1:2*N,:),[ne nq o.nfreq]),[3 2 1]);
% o.dat.realepsT = permute(reshape(eps(2*N+1:3*N,:),[ne nq o.nfreq]),[3 2 1]);
% o.dat.imagepsT = permute(reshape(eps(3*N+1:4*N,:),[ne nq o.nfreq]),[3 2 1]);

N = ne*nq*NFREQ;
o.dat.realepsL = permute(reshape(eps(0*N+1:1*N),[NFREQ ne nq]),[1 3 2]);
o.dat.imagepsL = permute(reshape(eps(1*N+1:2*N),[NFREQ ne nq]),[1 3 2]);
o.dat.realepsT = permute(reshape(eps(2*N+1:3*N),[NFREQ ne nq]),[1 3 2]);
o.dat.imagepsT = permute(reshape(eps(3*N+1:4*N),[NFREQ ne nq]),[1 3 2]);
end

