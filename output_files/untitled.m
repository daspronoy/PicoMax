filename = ['temp.dat'];
fid = fopen(filename);
data= fread(fid,'double');
fclose(fid);

eband=reshape(data(1:end,:),[40 31]);
set(gcf,'renderer','opengl', 'Position',  [500, 100, 500, 500])

plot(1:31,eband(1:end,:), 'Color','#0072BD','LineWidth',0.5);

set(gca,'fontsize',20,'linewidth',0.5);