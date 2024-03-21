







function [] = epsij(o)
picomax_cmd = sprintf( ...
    "-switch epsij " + ...
    "-encut %d -neps %d " + ...
    "-nband %d -originshift %d " + ...
    "-a %f -epm %s " + ...
    "-nfreq %d -dfreq %f -kk %d -epsilon %f -delta %d ", ...
    o.num.encut,o.num.neps, ...
    o.num.nband,o.num.originshift, ...
    o.sys.a,regexprep(num2str(o.sys.epm),'\s+',','), ...
    o.num.nfreq,o.num.dfreq,o.num.kk,o.num.epsilon,o.num.delta);


% 
picomax_cmd = sprintf( ...
    "%s -wdir %s " + ...
    "-outputfile %s ", ...
    picomax_cmd, ...
    [o.num.base_remote o.num.wdir], ...
    o.num.outputfile);



% BZ integration (KPOINT)
if o.num.kpoint==0
    picomax_cmd = sprintf( ...
        "%s -kpoint %d -kpointfile %s ", ...
        picomax_cmd, ...
        o.num.kpoint,o.num.kpointfile);
elseif o.num.kpoint==3 % Cohen-Chadi grid
    picomax_cmd = sprintf( ...
        "%s -kpoint %d -kpointorder %d ", ...
        picomax_cmd, ...
        o.num.kpoint,o.num.kpointorder);
elseif o.num.kpoint==2 % Monkhorst-Pack grid
    picomax_cmd = sprintf( ...
        "%s -kpoint %d -kpointorder %d ", ...
        picomax_cmd, ...
        o.num.kpoint,o.num.kpointorder);
end

% q vectors
if (~isempty(o.num.qnum) && ~isempty(o.num.qpath))
    picomax_cmd = ...
        sprintf("%s -qnum %s -qpath %s", ...
        picomax_cmd, ...
        regexprep(num2str(o.num.qnum),'\s+',','), ...
        o.num.qpath);
else
    picomax_cmd = ...
        sprintf("%s -qvec %s", ...
        picomax_cmd, ...
        regexprep(num2str(o.num.qvec(:).'),'\s+',','));
end


% construct system cmd
picomax_cmd = sprintf("%s %s", ...
    o.num.picomax_exe,picomax_cmd);
switch o.num.system
    case 'linux'
        system_cmd = picomax_cmd;
    case 'wsl'
        system_cmd = sprintf( ...
            "powershell -command ""wsl bash -c '%s'""", ...
            picomax_cmd);
end

% run the system cmd
tic
[status,cmdout] = system(system_cmd);
toc

% save the commands for debugging purposes
o.tmp.system_cmd = system_cmd;
o.tmp.picomax_cmd = picomax_cmd;
o.tmp.cmdout = cmdout;


% import the result file (.dat)
datfile = [o.num.base_local o.num.wdir o.num.outputfile '.dat'];
% fid = fopen(datfile,'r');
% eps = fread(fid,'double');
% fclose(fid);


o.import(datfile,'epsij');

% nfreq = o.num.nfreq;
% o.dat.omega = o.num.dfreq*(0:1:nfreq-1);
% if (~isempty(o.num.qnum) && ~isempty(o.num.qpath))
%     nq = sum(o.num.qnum)+1;
%     o.dat.qi = 0:nq-1;
% else
%     nq = 1;
%     o.dat.qi = 0;
% end
% ne = o.num.neps*(o.num.neps+1)/2;
% N = ne*nq*nfreq;
% 
% 
% o.dat.realepsL = permute(reshape(eps(0*N+1:1*N),[nfreq ne nq]),[1 3 2]);
% o.dat.imagepsL = permute(reshape(eps(1*N+1:2*N),[nfreq ne nq]),[1 3 2]);
% o.dat.realepsT = permute(reshape(eps(2*N+1:3*N),[nfreq ne nq]),[1 3 2]);
% o.dat.imagepsT = permute(reshape(eps(3*N+1:4*N),[nfreq ne nq]),[1 3 2]);
end

