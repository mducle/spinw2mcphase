% Calculations of Gd2Pt2O7
% Original author: J. R. Stewart

yl = 1;

%% McPhase setup
% Some parameters to be defined for your system
global mcphasedir; mcphasedir='c:/mcphas5_3';  % This is the folder to McPhase - it should be something like c:/mcphas5_3
global workdir; workdir='mtmp';  % This is the working directory for the McPhase calculation. If it is a relative path it will be created off the current folder
global numproc; numproc=2;       % Number of processors to use for the parallel calculation of the dispersion
global emin; emin=-200;          % Minimum energy of modes to calculate the intensity for (modes below this energy will have zero intensity to save computation time)
global emax; emax=200;           % Maximum energy of modes to calculate the intensity for (modes above this energy will have zero intensity to save computation time)
global temperature; temperature = 0.1;

%%
% Generate GPO lattice
% This is the Palmer-Chalker Structure. 
symStr = ['-z,y+3/4,x+3/4; z+3/4,-y,x+3/4; z+3/4,y+3/4,-x; '...
          'y+3/4,x+3/4,-z; x+3/4,-z,y+3/4; -z,x+3/4,y+3/4'];
% 2nd origin setting for F-d3m

gpo = spinw;
a = 10.225;
x = 0.3348;
gJ = 2;                  % g-factor
S  = 3.5;
r_limit = 20;

gpo.genlattice('lat_const',[a a a],'angled',[90 90 90],'spgr',symStr,...
    'label','F d -3 m Z');
gpo.addatom('label','Gd3+','r',[1/2 1/2 1/2],'S',7/2);
%gpo.addatom('label','Pt4+','r',[0 0 0]);
%gpo.addatom('label','O2-','r',[x 1/8 1/8]);
%gpo.addatom('label','O2-','r',[3/8 3/8 3/8]);

gpo.gencoupling('maxDistance',r_limit)   

%% input Palmer-Chalker structure
S1 = [1 -1 0];      % (1/2 1/2 1/2) GD_1
S2 = [-1 -1 0];     % (1/2 1/4 1/4) GD_4
S3 = [1 1 0];       % (3/4 0   1/4) GD_3
S4 = [-1 1 0];      % (3/4 3/4 1/2) GD_2
S5 = [1 1 0];       % (1/4 1/2 1/4) GD_3
S6 = [-1 -1 0];     % (0   3/4 1/4) GD_4
S7 = [1 -1 0];      % (0   0   1/2) GD_1
S8 = [-1 1 0];      % (1/4 1/4 1/2) GD_2
S9 = [1 -1 0];      % (0   1/2 0  ) GD_1
S10= [-1 -1 0];     % (0   1/4 3/4) GD_4
S11= [1 1 0];       % (1/4 0   3/4) GD_3
S12= [-1 1 0];      % (1/4 3/4 0  ) GD_2
S13= [1 1 0];       % (3/4 1/2 3/4) GD_3
S14= [-1 -1 0];     % (1/2 3/4 3/4) GD_4
S15= [1 -1 0];      % (1/2 0   0  ) GD_1
S16= [-1 1 0];      % (3/4 1/4 0  ) GD_2
spins = cat(3,S1',S2',S3',S4',S5',S6',S7',S8',S9',S10',S11',S12',S13',S14',S15',S16');
spins = permute(spins,[1 3 2]);

gpo.genmagstr('k',[0 0 0],'S',spins);
magatoms = gpo.table('mag');

%%  Interaction model - J1, J2 and J3b + D-aniso and Dipolar
boltzmann = 0.08617;                % meV K^-1

J1 = 4.81 * boltzmann / S / (S + 1) ;
gpo.addmatrix('label','J1','value',J1);
gpo.addcoupling('mat','J1','bond',1);

J2 = 0.139 * boltzmann / S / (S + 1) ;
gpo.addmatrix('label','J2','value',J2);
gpo.addcoupling('mat','J2','bond',2);

fprintf('\n Heisenberg interactions:')
fprintf('\n J1  = %0.5f meV (%0.5f K); J1*S(S+1)  = %0.5f K', J1, J1/boltzmann, J1*S*(S+1)/boltzmann)
fprintf('\n J2  = %0.5f meV (%0.5f K); J2*S(S+1)  = %0.5f K', J2, J2/boltzmann, J2*S*(S+1)/boltzmann)

%% add easy plane (111) anisotropies 

D = boltzmann * 3.52 / S / S;        % from refs [22,23] - check with Phil
gpo.addmatrix('value',D*[1 1 1 ; 1 1 1 ; 1 1 1]/sqrt(3),'label','D1','color','r');
gpo.addaniso('D1');

fprintf('\n Single Ion Anisotropy:')
fprintf('\n D =%0.5f meV (%0.5f K); D*S^2=%0.5f K \n\n', D, D/boltzmann, D*S*S/boltzmann)

%% Dipolar interactions - Sandor

gpo.coupling.rdip=r_limit;

%% Optimise structure and minimise energy

energy1 = gpo.energy;
gpo.optmagsteep('nRun',1e3);
energy2 = gpo.energy;

fprintf('\n Energy before optimisation: %0.4f \n',energy1)
fprintf(' Energy after optimisation: %0.4f \n',energy2)

if (abs(energy1 - energy2) > 0.001)
    disp('Warning - magnetic structure changed by optmagsteep....')
end

%% Plot 

swpref.setpref('nmesh',3,'npatch',40);
plot(gpo,'atomMode','mag','atomColor','yellow','magColor',...
         'orange','bondMode','cylinder')

%% Calculate spectra
spec = sw_egrid(sw_neutron(gpo.spinwave({[0 1.5 0] [0 1 0] [0 0 0] [0.5 0.5 0.5] [1 1 0] [0 1 0] [0.5 0.5 0] 100},'hermit',true)),'Evect',0:0.001:yl);

figure
subplot(211);
sw_plotspec(spec,'mode','color','dE',0.025,'dashed',true)
legend('off'); ylim([0 yl]); title('SpinW - 0T');

%% Now do the same using McPhase
hkl = sw_qscan({[0 1.5 0] [0 1 0] [0 0 0] [0.5 0.5 0.5] [1 1 0] [0 1 0] [0.5 0.5 0] 100});
specm = sw_egrid(sw_neutron(mcphase_sqw(hkl(1,:), hkl(2,:), hkl(3,:), gpo)),'Evect',0:0.001:yl);

subplot(212);
sw_plotspec(specm,'mode','color','dE',0.025,'dashed',true)
legend('off'); ylim([0 yl]); title('McPhase - 0T');
