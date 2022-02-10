global mcphasedir; mcphasedir='c:/mcphas5_3';  % This is the folder to McPhase - it should be something like c:/mcphas5_3
global workdir; workdir='mtmp';  % This is the working directory for the McPhase calculation. If it is a relative path it will be created off the current folder
global numproc; numproc=2;       % Number of processors to use for the parallel calculation of the dispersion
global emin; emin=-200;          % Minimum energy of modes to calculate the intensity for (modes below this energy will have zero intensity to save computation time)
global emax; emax=200;           % Maximum energy of modes to calculate the intensity for (modes above this energy will have zero intensity to save computation time)
global temperature; temperature = 35;

lfpo = spinw();
lfpo.genlattice('lat_const',[10.33 6.01 4.69],'angled',[90 90 90],'spgr','P n m a');
lfpo.addatom('label','MFe2','r',[0.28222 0.25 0.97472],'S',2);
lfpo.gencoupling('maxDistance', 7)
lfpo.addmatrix('label', 'Jbc', 'value', 0.46); lfpo.addcoupling('mat', 'Jbc', 'bond', 1);
lfpo.addmatrix('label', 'Jb', 'value', 0.09);  lfpo.addcoupling('mat', 'Jb', 'bond', 6);
lfpo.addmatrix('label', 'Jc', 'value', 0.01);  lfpo.addcoupling('mat', 'Jc', 'bond', 2);
lfpo.addmatrix('label', 'Jab', 'value', 0.09); lfpo.addcoupling('mat', 'Jab', 'bond', 7);
lfpo.addmatrix('label', 'Jac', 'value', 0.01); lfpo.addcoupling('mat', 'Jac', 'bond', 4);
lfpo.addmatrix('label', 'D', 'value', diag([0.86 0 2.23])); lfpo.addaniso('D');

lfpo.genmagstr('mode', 'direct', 'S', [0 0 0 0; 1 -1 -1 1; 0 0 0 0], 'k', [0 0 0]);
plot(lfpo, 'range', [2 2 2])

qpts = {{[-1 1 0] [1 1 0] 50} {[0 0 0] [0 4 0] 50} {[0 0 0] [0 0 3] 50}};
for ii = 1:numel(qpts)
    spec(ii) = sw_egrid(sw_neutron(lfpo.spinwave(qpts{ii},'hermit',true)),'Evect',0:0.01:15);
    hkl = sw_qscan(qpts{ii});
    specm(ii) = sw_egrid(sw_neutron(mcphase_sqw(hkl(1,:), hkl(2,:), hkl(3,:), lfpo)),'Evect',0:0.01:15);
end

figure; sp = [1 2 3; 4 5 6]; 
for ii = 1:numel(qpts)
    subplot(2,3,sp(1,ii));
    sw_plotspec(spec(ii),'mode','color','dE',0.1,'dashed',true); title('SpinW'); legend('off');
    subplot(2,3,sp(2,ii));
    sw_plotspec(specm(ii),'mode','color','dE',0.1,'dashed',true); title('McPhase'); legend('off');
end
