

# Matlab wrapper

Matlab wrapper scripts for PicoMax software.

Windows SubLinux (WSL) 


Case 1: PicoMax is installed in WSL and MATLAB is installed in Windows.
Submits commmands to picomax installed in WSL through powershell
`base_local`
`base_remote`
`

```MATLAB
c = picomax;
...
c.setopt('system','wsl')
c.setopt('base_local','//wsl.localhost/ubuntu-22.04/')
c.setopt('base_remote','/')
c.setopt('picomax_exe','/path/to/picomax_exe/picomax')
```


Case 2: Both PicoMax and MATLAB are installed in linux (or WSL)
Submits commands to picomax through MATLAB `system` command
```MATLAB
c = picomax;
...
c.setopt('system','linux')
...
```




## Appendix A. Examples

### Electron band structure of silicon carbide (SiC) in WSL system

```MATLAB
clear;clc;

% material parameters of SiC
a = 4.35;
g2 = [3 4 8 11];
v1 = [-0.225 0.005 0.118 0.127];
v2 = [-0.545 -0.005 0.118 0.121];
vg = [g2 v1 v2];
crystal = 'zincblende';
nvalence = 4;

% numerical parameters
encut = 300;
nband = 20;

% q-point path
qnum = [20,20,10,0,20];
qpath = 'L,G,X,U,K,G';

% solve
c = picomax;
c.setopt('encut',encut,'nband',nband,'nvalence',nvalence, ...
    'qpath',qpath,'qnum',qnum);
c.set('crystal',crystal,'a',a,'vg',vg)
c.setopt('force_cpp',1)
c.setopt('picomax_exe','/path/to/picomax_exe/picomax')
c.setopt('base_remote','/')
c.setopt('base_local','//wsl.localhost/ubuntu-22.04/')
c.run('eband') % solve for bandstructure

% plot
figure
set(gcf,'units','centimeters','Position',[4 4 6 5])
cmap = colororder("gem12");
plot(0:sum(qnum),c.dat.electron.band,'-','LineWidth',1.0,'Color',cmap(1,:))
xlim([0 sum(qnum)])
colororder("gem12")
qpathlabel = {'L','\Gamma','X','U/K','\Gamma'};
set(gca,'XTick',unique([0 cumsum(qnum)]), ...
    'XTickLabel',qpathlabel, ...
    'FontName','Aptos', ...
    'FontSize',8.0, ...
    'Layer','top', ...
    'TickDir','out', ...
    'XMinorTick','on','YMinorTick','on', ...
    'Color','none')
ylabel('Energy (eV)','FontName','Aptos')
xline(cumsum(qnum(1:end-1)),':')
grid on
set(gca,'YTick',-40:5:40)
ylim([-20 15])
```


### Deep microscopic optical band structure of SiC

```matlab
clear;clc;

% material parameters of SiC
a = 4.35;
g2 = [3 4 8 11];
v1 = [-0.225 0.005 0.118 0.127];
v2 = [-0.545 -0.005 0.118 0.121];
vg = [g2 v1 v2];
crystal = 'zincblende';
nvalence = 4;

% numerical parameters
encut = 300;
nband = 20;
nchi = 51;
nomega = 2;
domega = 0.04135665538536; % 10 THz

% q-point path
qnum = [20,20,10,0,20];
qpath = 'L,G,X,U,K,G';

% solve
c = picomax;
c.setopt('encut',encut,'nchi',nchi,'nband',nband,'nvalence',nvalence, ...
    'qpath',qpath,'qnum',qnum, ...
    'nomega',nomega,'domega',domega,'kk_transform',0,'epsilon',0.3,'delta',0);
c.set('crystal',crystal,'a',a,'vg',vg)
c.setopt('kptopt',1,'kptorder',11)
c.setopt('picomax_exe','/path/to/picomax_exe/picomax')
c.setopt('base_remote','/')
c.setopt('base_local','//wsl.localhost/ubuntu-22.04/')
c.run('Xij')


nchi = 
Xij = c.unwrap()


```

### Visualize deep microscopic optical eigenwaves in SiC




### Wurtzite crystal CdSe

```MATLAB



```
