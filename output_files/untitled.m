filename = ['eband_141pts_KLHGAM_0.003eV.dat'];
fid = fopen(filename);
data= fread(fid,'double');
fclose(fid);

npts=141;
eband=reshape(data(1:end,:),[40 npts]);
set(gcf,'renderer','opengl', 'Position',  [500, 100, 500, 500])

plot(1:npts,eband(1:end,:), 'Color','#505050','LineWidth',1.5);

xticks([1 21 61 81 121 141]);
% xticklabels({"$\Gamma$", "A", "H", "K", "$\Gamma$"});
xticklabels({"K", "L", "H", "$\Gamma$", "A", "M"});

ylabel("Energy (eV)","Interpreter","latex");

ylim([-2 2]);
xlim([1 npts]);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 20, 'LineWidth', 0.5);