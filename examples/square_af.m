global mcphasedir; mcphasedir='c:/mcphas5_3';  % This is the folder to McPhase - it should be something like c:/mcphas5_3
global workdir; workdir='mtmp';  % This is the working directory for the McPhase calculation. If it is a relative path it will be created off the current folder
global numproc; numproc=2;       % Number of processors to use for the parallel calculation of the dispersion
global emin; emin=-200;          % Minimum energy of modes to calculate the intensity for (modes below this energy will have zero intensity to save computation time)
global emax; emax=200;           % Maximum energy of modes to calculate the intensity for (modes above this energy will have zero intensity to save computation time)
global temperature; temperature = 0.02;

swobj = sw_model('squareAF', [1 0.3]);
swobj.genmagstr('mode', 'direct', 'S', [1 0 0], 'n', [0 0 1], 'k', [0.5 0.5 0]);
qpts = {[0 1.5 0] [0 1 0] [0 0 0] [0.5 0.5 0.5] [1 1 0] [0 1 0] [0.5 0.5 0] 100};
spec = sw_egrid(sw_neutron(swobj.spinwave(qpts,'hermit',true)),'Evect',0:0.01:5);
hkl = sw_qscan(qpts);
specm = sw_egrid(sw_neutron(mcphase_sqw(hkl(1,:), hkl(2,:), hkl(3,:), swobj)),'Evect',0:0.01:5);

figure; 
subplot(211); sw_plotspec(spec,'mode','color','dE',0.1,'dashed',true); title('SpinW');
subplot(212); sw_plotspec(specm,'mode','color','dE',0.1,'dashed',true); title('McPhase');
