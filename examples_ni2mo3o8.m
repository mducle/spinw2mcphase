% Some parameters to be defined for your system
global mcphasedir; mcphasedir='c:/mcphas5_3';  % This is the folder to McPhase - it should be something like c:/mcphas5_3
global workdir; workdir='mtmp';  % This is the working directory for the McPhase calculation. If it is a relative path it will be created off the current folder
global numproc; numproc=8;       % Number of processors to use for the parallel calculation of the dispersion
global emin; emin=-200;          % Minimum energy of modes to calculate the intensity for (modes below this energy will have zero intensity to save computation time)
global emax; emax=200;           % Maximum energy of modes to calculate the intensity for (modes above this energy will have zero intensity to save computation time)
 
%% Defines the model.
% Uses the bilinear-biquadratic model of Li et al., PRM 6 014405 (2022) for details
% Magnetic structure from Morey et al., PRM 3 014410

% Magnitude of the DM interaction in meV;
DM = 3;

% There are two models of magnetic structures, define them as cell arrays here.
mdls = {[-0.994 0.07; 0.473 -0.586] [-0.955 -0.559; 0.514 0.021]};   % From PRM 3 014410 Supp Mat
bv = [2 1 0; 0 0 2];  % From PRM 3 014410 Table IV
% Note that we manually change the c-axis moment between layers in code below instead of using difference basis vectors as in Table IV

% Now set up SpinW model with published magnetic structure
model = 1;   % Use model 1 or 2

h7 = spinw;
h7.genlattice('lat_const', [5.74683 5.74683 9.86260], 'angled', [90 90 120], 'spgr', 'P 63 m c')
h7.addatom('r', [2/3 1/3 0.4480], 'S', 1, 'label', 'Ni2+', 'color', 'blue') % Tetrahedral site
h7.addatom('r', [1/3 2/3 0.5116], 'S', 1, 'label', 'Ni2+', 'color', 'red')  % Octahedral site

% Recalculate angles and magnetic moments - should reproduce part of Table V of PRM 3 014410
sites = {'Tet', 'Oct'};
disp('azimuth(deg) elevation(deg) moment(ub)')
for imd=1:numel(mdls);
    for iss=1:2;
        cmg{imd,iss} = h7.basisvector(true) * (mdls{imd}(iss,:) * bv)'; % .basisvector(true) give the normalised basis vector
        smg = mat2cell(cmg{imd,iss},ones(1,3),1); 
        [smg{:}] = cart2sph(smg{:});
        disp(sprintf('%11.4f %11.4f %11.4f     model %i %s',cell2mat(smg)'.*[[1 1]*180/pi 1], imd, sites{iss}));
    end;
end
disp(sprintf('\n'));

% Define the moment vector in the structural unit cell in lattice units (from Morey paper - component * basis_vectors)
Svec = zeros(3, 2);
for iss = 1:2;
    % There are two layers in the structure, and the layers have spins antiparallel in the c-direction
    S0 = mdls{model}(iss,:) * bv;
    Svec(:,(iss-1)*2+1) = S0;
    Svec(:,(iss-1)*2+2) = S0 .* [1 1 -1];  % This would have been achieved by the different basis vectors in Table IV
end
Svec;
kvec = [1/2 0 0]; nExt = 1./kvec; nExt(isinf(nExt)) = 1;
h7.genmagstr('mode', 'helical', 'k', kvec, 'S', Svec, 'n', [0 1 0], 'unit', 'lu')

h7.gencoupling('maxDistance', 6, 'forceNoSym', true)

% Gets the angle between nearest neighbour bonds
disp(sprintf('\n'));
tb = h7.table('bond', 1);
Sv = h7.magstr('nExt', nExt).S; % This is the magnetic moment in Cartesian units in the extended cell
h7.genmagstr('mode', 'direct', 'k', [0 0 0], 'nExt', nExt, 'S', Sv, 'n', [0 1 0])
for ib = 1:6
    sdb(ib) = dot(Sv(:,tb.idx1(ib)), Sv(:,tb.idx2(ib) + tb.dl(ib,1)*4));
end
disp('Angle between spins')
sdbd = (acos(sdb) * 180 / pi);
disp(sdbd);
[usdbd, iu, ia] = unique(round(sdbd*100)/100);

% Gets the indices for the J+ / J- couplings
if numel(iu) > 2; error('Something went wrong in the calculation of the biquadratic term - there are more than two inequivalent nn couplings'); end
idxbqp = find(ia==2)';
idxbqm = find(ia==1)';
Jp = 2.62 - 1.13*(5*sdb(iu(2))/4 - 0.5);
Jm = 2.62 - 1.13*(5*sdb(iu(1))/4 - 0.5);
disp(['Effective nn couplings are J+=' num2str(Jp) ' J-=' num2str(Jm)]);
disp(sprintf('\n'));

% Defines the exchanges - parameters notation from Li paper
h7.addmatrix('label', 'J1+', 'value', Jp, 'color', 'gray');
h7.addmatrix('label', 'J1-', 'value', Jm, 'color', 'green');
h7.addmatrix('label', 'JTF', 'value', 0.09, 'color', 'orange');
h7.addmatrix('label', 'JTA', 'value', 0.96, 'color', 'cyan');
h7.addmatrix('label', 'JOF', 'value', 0.48, 'color', 'yellow');
h7.addmatrix('label', 'JOA', 'value', 0.25, 'color', 'blue');
v1 = -[-1 sqrt(3) 0]; v1 = v1 ./ norm(v1);
v2 = [0.5 1 0]; v2 = v2 ./ norm(v2);
v3 = [1 0 0]; v3 = v3 ./ norm(v3);
h7.addmatrix('label', 'DM1', 'value',  DM*v1, 'color', 'green')
h7.addmatrix('label', 'DM2', 'value',  DM*v2, 'color', 'green')
h7.addmatrix('label', 'DM3', 'value',  DM*v3, 'color', 'green')
h7.addmatrix('label', 'DM4', 'value', -DM*v2, 'color', 'green')

h7.addcoupling('mat', 'J1+', 'bond', 1, 'subIdx', idxbqp)
h7.addcoupling('mat', 'J1-', 'bond', 1, 'subIdx', idxbqm)

plot(h7, 'range', [0 0  0.2; 2 2 0.8]')
M = swplot.transform(gcf); 
M(1:3, 1:3) = [cosd(30) sind(30) 0; -sind(30) cosd(30) 0; 0 0 1];  % a-b plane, rotated 30deg
%M(1:3, 1:3) = [cosd(30) sind(30) 0; 0 0 1; sind(30) -cosd(30) 0];  % a-c plane, rotated 30deg
swplot.transform(M, gcf);
hcm = camlight('right')
%return;  % Stops script here to get just picture of structure with nearest neighbour couplings plotted only.

h7.addcoupling('mat', 'JTF', 'bond', 4, 'subIdx', [5 3])
h7.addcoupling('mat', 'JTA', 'bond', 4, 'subIdx', [4 11 1 10])
h7.addcoupling('mat', 'JOF', 'bond', 4, 'subIdx', [2 9])
h7.addcoupling('mat', 'JOA', 'bond', 4, 'subIdx', [6 7 8 12])
%h7.addcoupling('mat', 'JOF', 'bond', 4, 'subIdx', [7 12])    % For checking Zn calcs
%h7.addcoupling('mat', 'JOA', 'bond', 4, 'subIdx', [6 2 8 9]) % For checking Zn calcs

h7.addcoupling('mat', 'DM1', 'bond', 1, 'subIdx', [4 5])
h7.addcoupling('mat', 'DM2', 'bond', 1, 'subIdx', [1])
h7.addcoupling('mat', 'DM3', 'bond', 1, 'subIdx', [3 6])
h7.addcoupling('mat', 'DM4', 'bond', 1, 'subIdx', [2])

%h7.addmatrix('label', 'JC1', 'value', 2.0); h7.addcoupling('mat', 'JC1', 'bond', 2);
%h7.addmatrix('label', 'JC2', 'value', 0.0); h7.addcoupling('mat', 'JC2', 'bond', 3);

% SIA using parameters defined in Li et al. rather than Morey et al.
h7.addmatrix('value', diag([0 0 30.41]), 'label', 'DT', 'color', 'green')
h7.addmatrix('value', diag([0 0 -0.53]), 'label', 'DO', 'color', 'blue')
h7.addaniso('DT', 1)
h7.addaniso('DO', 2)

%% Runs pure SpinW calculation
S0 = h7.magstr.S';
optm = h7.optmagsteep('nRun', 1000)
h7 = optm.obj;
[S0 S0(:, 1)*0 h7.magstr.S']
[az0, el0, r] = cart2sph(S0(:,1), S0(:,2), S0(:,3)); 
[az1, el1, r] = cart2sph(h7.magstr.S(1,:), h7.magstr.S(2,:), h7.magstr.S(3,:));
[az0 el0 az0*0 az1(:) el1(:)].*(180/pi)

plot(h7, 'range', [0 0 0.8; 4 2 1.3]');
M = swplot.transform(gcf); 
%M(1:3, 1:3) = [cosd(30) sind(30) 0; -sind(30) cosd(30) 0; 0 0 1];    % a-b plane, rotated 30deg
M(1:3, 1:3) = [-cosd(30) -sind(30) 0; 0 0 1; -sind(30) cosd(30) 0];  % a-c plane, rotated 30deg
swplot.transform(M, gcf);

% Define twins
h7.twin.rotc(:,:,2) = [cosd(120) sind(120) 0; -sind(120) cosd(120) 0; 0 0 1];
h7.twin.rotc(:,:,3) = [cosd(240) sind(240) 0; -sind(240) cosd(240) 0; 0 0 1];
h7.twin.vol = [1/3 1/3 1/3];

hkl = sw_qscan({[-1 -1 0] [1 1 0] [1 1 1] 400});
spec_sw = h7.spinwave(hkl, 'hermit', true, 'formfact', true);
spec_sw = sw_neutron(spec_sw);

figure
spec_sw = sw_egrid(spec_sw, 'Evect', linspace(0, 60, 800), 'component', 'Sperp', 'T', 1.6);
subplot(2, 1, 1)
sw_plotspec(spec_sw, 'mode', 'disp', 'axLim', [0 60], 'colormap', [0 0 0], 'colorbar', false)
subplot(2, 1, 2)
sw_plotspec(spec_sw, 'mode', 'color', 'axLim', [0 4], 'dE', 1)
set(gcf, 'PaperPosition', [1.2813 0.4792 5.9375 10.0417]);
set(gcf, 'Position', [675 4 570 964]);

%% Runs the McPhase calculations
MoC = 9; OxC = -4; NiC = 1; rr = 3.5;   % Theory paper parameters (Shuyi+Andriy)

global temperature; temperature = 0.02;

h7.cache.mcphase = struct('pointcharge', struct('maxdist', rr, ...
          'atoms', struct('label', {'Ni1', 'Ni2', 'Mo', 'O1', 'O2', 'O3', 'O4'}, ...
                          'valence', {2 2 4 -2.0 -2.0 -2.0 -2.0}, ...  % Formal valences, for form factor lookup and electronic configuration
                          'charges', {NiC NiC MoC OxC OxC OxC OxC}, ...  % Effective charge, for point charge calculation (optional)
                          'r', {[2/3 1/3 0.448] [1/3 2/3 0.5116] [0.144 -0.144 0.2489] [0 0 0.6839] [1/3 2/3 0.1461] [0.488 -0.488 0.3659] [0.1688 -1/6 0.6342]})))

% Using CEF parameters directly
%{
% Pars in Pengcheng's paper's SI
h7.unit_cell.label = {'Ni1', 'Ni2'};
h7.cache.mcphase = struct('Ni1', struct('module', 'ic1ion', 'iontype', 'Ni2+', 'L20', 735, 'L40', -882, 'L43', 117.138), ... % [0 18 18 288 288]
                          'Ni2', struct('module', 'ic1ion', 'iontype', 'Ni2+', 'L20', 105, 'L40', -630, 'L43', 1597.3));     % [0 0.79 0.79 743 743]
%}

hh = hkl(1,:); kk = hkl(2,:); ll = hkl(3,:);
% This will take about 5 min to run
spec_mc = mcphase_sqw(hh,kk,ll,h7)

figure
spec_mc = sw_egrid(spec_mc, 'Evect', linspace(0, 60, 800), 'component', 'Sperp', 'T', 1.6);
subplot(2, 1, 1)
sw_plotspec(spec_mc, 'mode', 'disp', 'axLim', [0 60], 'colormap', [0 0 0], 'colorbar', false)
subplot(2, 1, 2)
sw_plotspec(spec_mc, 'mode', 'color', 'axLim', [0 2], 'dE', 1)
set(gcf, 'PaperPosition', [1.2813 0.4792 5.9375 10.0417]);
set(gcf, 'Position', [675 4 570 964]);

% Prints out the magnetic structure found from the mean-field self-consistent loop
% compared with that determined by the minimization in SpinW
[az0, el0, r] = cart2sph(S0(:,1), S0(:,2), S0(:,3)); 
[az1, el1, r] = cart2sph(spec_mc.obj.magstr.S(1,:), spec_mc.obj.magstr.S(2,:), spec_mc.obj.magstr.S(3,:));
disp(' -------------- ')
disp('Azi_ref Polar_ref | Azi_mf Polar_mf')
[az0 el0 az0*NaN az1(:) el1(:)].*(180/pi)
disp('ref is the experimental value; mf is the mean-field calculation')
disp(sprintf('Canting angle = %0.3f deg', max(abs(el1))*180/pi));
disp(' -------------- ')

% Plots out the resulting structure
plot(spec_mc.obj, 'range', [0 0 0.8; 4 2 1.3]');
M = swplot.transform(gcf); 
M(1:3, 1:3) = [cosd(30) sind(30) 0; -sind(30) cosd(30) 0; 0 0 1];     % a-b plane, rotated 30deg
%M(1:3, 1:3) = [-cosd(30) -sind(30) 0; 0 0 1; -sind(30) cosd(30) 0];  % a-c plane, rotated 30deg
swplot.transform(M, gcf);

% Plots the spectra from SpinW and McPhase together
figure
spec_mc = sw_egrid(spec_mc, 'Evect', linspace(0, 60, 800), 'component', 'Sperp', 'T', 1.6);
subplot(2, 1, 1)
sw_plotspec(spec_sw, 'mode', 'color', 'axLim', [0 2], 'dE', 1); title('SpinW'); legend('off'); ylim([0 40])
subplot(2, 1, 2)
sw_plotspec(spec_mc, 'mode', 'color', 'axLim', [0 2], 'dE', 1); title('McPhase'); legend('off'); ylim([0 40])
set(gcf, 'PaperPosition', [1.2813 0.4792 5.9375 10.0417]);
set(gcf, 'Position', [675 4 570 964]);
