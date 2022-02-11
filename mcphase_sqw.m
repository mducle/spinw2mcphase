function spec = mcphase_sqw(h,k,l,i,o,nmf)
% Function to calculate the dispersion relation using McPhase from a SpinW model
%
% NB. Before using, you should set the McPhase path in a global variable, e.g.:
%   >> global mcphasedir; mcphasedir = 'D:\mcphase\';
%
% If this is not set, the function assumes all McPhase programs are in the current
%   directory - this functions needs the programs "mcphasit", "mcdispit", "spins".
%
% Syntax: [w,s] = mcphase_sqw(h,k,l,o)
%   w is a cell of the energies (hbar-omega) in meV at particular (hkl) points
%   s is a cell of the corresponding cross-sections (in barns/Sr/f.u.)
%   h is a vector of the input Qh (in r.l.u.)
%   k is a vector of the input Qk (in r.l.u.)
%   l is a vector of the input Ql (in r.l.u.) - h,k,l must be the same size!
%   o is a SpinW object
%  
% The exchange parameters and crystal structure for the McPhase model will be 
% obtained from the SpinW object.
% If you don't want to use the SIA parameters in SpinW, but instead use
% crystal field parameters directly or a point charge model, you should set
% a "mcphase" field in the SpinW .cache field as described below.
%
% Example:
%   >> proj.u=[1,0,0]; proj.v=[0,0,1]; proj.type='rrr';proj.uoffset=[0,0,0,0];
%   >> sqw_file='/data/home/ctf81281/duc/ei25-ortho.sqw';
%   >> w3=cut_sqw(sqw_file,proj,[-6,3],[-6,0.15,6],[],[],'-nopix');
%   >> w3calc=sqw_eval(w3,@mcphase_sqw,1,{swobj})
%
% The function can also use the Matlab provided perl interpreter to launch parallel
% processes using Perls fork() emulation on Windows (and native behaviour on Linux). 
%
% In order to do this, you need to set up a global variable:
%
%   >> global numproc; numproc = 4;    % To run 4 processes
%
% Windows is not very good at scheduling computationally intensive processes, so
% note that if you set numproc to the actual number of processors you have, you 
% might find that all other applications become unresponsive.
%
% This function sets the priority of the McPhase calculations to "BelowNormal" but this
% might not work so well.
%
% There are also three other global parameters which can be set:
%
%   magfield - a 3-vector which sets the applied magnetic field in Tesla (default: [0 0 0])
%   temperature - a scalar which sets the measurement temperature in Kelvin (defulat: 2K)
%   workdir - a string with the path to the working directory (default: current folder)
%
% In order to use single-ion properties you should create a "mcphase" field in the
% SpinW object's .cache property.
% If you want to input the crystal field parameters directly you should 

% This field may be a cell array or struct.
% If it is a cell array it should contain n_at elements where n_at is the number of
% magnetic atoms in the unit cell. If it is a struct, the fields of the struct
% should be the atom labels used in the SpinW model.
% The elements/fields should be a struct with 
%
% duc.le@stfc.ac.uk - Oct 13 12:03:48 BST 2021

% Example script:
%
% global mcphasedir; mcphasedir='c:/mcphas5_3';  % This is the folder to McPhase - it should be something like c:/mcphas5_3
% global workdir; workdir='mt';    % This is the working directory for the McPhase calculation. If it is a relative path it will be created off the current folder
% global numproc; numproc=2;       % Number of processors to use for the parallel calculation of the dispersion
% global emin; emin=-200;          % Minimum energy of modes to calculate the intensity for (modes below this energy will have zero intensity to save computation time)
% global emax; emax=200;           % Maximum energy of modes to calculate the intensity for (modes above this energy will have zero intensity to save computation time)
% global temperature; temperature = 0.02;
% 
% swobj = sw_model('squareAF', [1 0.3]);
% swobj.genmagstr('mode', 'direct', 'S', [1 0 0], 'n', [0 0 1], 'k', [0.5 0.5 0]);
% qpts = {[0 1.5 0] [0 1 0] [0 0 0] [0.5 0.5 0.5] [1 1 0] [0 1 0] [0.5 0.5 0] 100};
% spec = sw_egrid(sw_neutron(swobj.spinwave(qpts,'hermit',true)),'Evect',0:0.01:5);
% hkl = sw_qscan(qpts);
% specm = sw_egrid(sw_neutron(mcphase_sqw(hkl(1,:), hkl(2,:), hkl(3,:), swobj)),'Evect',0:0.01:5);
%
% figure; 
% subplot(211); sw_plotspec(spec,'mode','color','dE',0.1,'dashed',true); title('SpinW');
% subplot(212); sw_plotspec(specm,'mode','color','dE',0.1,'dashed',true); title('McPhase');

% Output struct
spec = struct('datestart', string(datetime), 'formfact', 1);

% Some definitions
global mcphasedir magfield temperature workdir emin emax; 
if ~ispc && isempty(mcphasedir); mcphasedir = '.'; end
if ~ispc && isempty(workdir); workdir = '.'; end
if isempty(magfield); magfield=[0 0 0]; end; 
if isempty(temperature); temperature = 2; end
if isempty(emin); emin = 0; end
if isempty(emax); emax = 100; end
nl = sprintf('\n'); if ispc; sp = '\'; else sp = '/'; end

if isa(i, 'spinw')
    if exist('o', 'var')
        if ~isa(o, 'spinw')
            if ~exist('nmf', 'var')
                nmf = o;
            end
            o = i;
        end
    else
        o = i;
    end
end
if ~exist('nmf', 'var'); nmf = false; else nmf = true; end

% Change to working direction
if exist(workdir, 'dir') ~= 7; mkdir(workdir); end
cwd = pwd;
cleanup = onCleanup(@()cd(cwd));

cd(workdir);
if exist('results','dir')~=7; mkdir('results'); end

% Generates the required input files
[infiles, tstinfo] = generate_mcphase_inputs(o);

% Runs the calculation
if nmf && exist('tmpmcphasobj.mat', 'file')
    load('tmpmcphasobj.mat');
else
    o = run_mcphase(copy(o), tstinfo);
    save('tmpmcphasobj.mat', 'o');
end
if numel(o.twin.vol) > 1
    tQ = o.twinq([h(:) k(:) l(:)]');
    for it = 1:numel(tQ)
        [w{it}, s{it}, sab{it}] = run_mcdisp(tQ{it}(1,:), tQ{it}(2,:), tQ{it}(3,:), infiles);
    end
else
    [w, s, sab] = run_mcdisp(h, k, l, infiles);
end

cd(cwd);

spec.omega = w;
spec.Sab = sab;
spec.hkl = [h(:) k(:) l(:)]';
spec.hklA = 2*pi*(spec.hkl'/o.basisvector)';
spec.incomm = 0;
spec.helical = 0;
spec.norm = 0;
spec.nformula = 0;
spec.param = struct('notwin', 0, 'sortMode', 1, 'tol', 1.0000e-04, 'omega_tol', 1.0000e-05, 'hermit', 1, 'n', [0 0 1]);
spec.param.uv = {}; spec.param.pol = 0;
spec.title = 'McPhase RPA spectrum';
spec.gtensor = 0;
spec.dateend = string(datetime);
spec.obj = o;
spec.intP = {[]  []  []};
spec.Pab = {[]  []  []};
spec.Mab = {[]  []  []};
spec.Sperp = s;

end

% % Writes the input parameters to a mcphase exchange constants input file [mcphas.j]
% if exist('makejfromtemplate.pl','file')~=2  
%   template=['#!/usr/bin/perl' nl ...
%             '# Reads in template file and converts into a single line' nl ...
%             'open(IFH,"mcphas.j-template"); @template = <IFH>; close(IFH); $templine=join("",@template);' nl ...
%             '# Creates the parameter list' nl ...
%             '@parlist = split('','',$ARGV[0]);' nl ...
%             '' nl ...
%             '# Goes through input list and replace [IJK]# with input.' nl ...
%             '$count=1; $pref="I"; $con=1;' nl ...
%             'foreach(@parlist)' nl ...
%             '{' nl ...
%             '  $templine =~ s/$pref$con/$parlist[$count-1]/g;' nl ...
%             '  if($pref=~/I/) {$pref="J";} elsif($pref=~/J/) {$pref="K";} elsif($pref=~/K/) {$pref="I";$con++;}' nl ...
%             '  $count=$count+1;' nl ...
%             '}' nl ...
%             '' nl ...
%             '# Writes output' nl ...
%             'open(OFH,">mcphas.j"); $count=1; print OFH "$templine"; close(OFH);' nl];
%   fid=fopen('makejfromtemplate.pl','w'); fprintf(fid,template); fclose(fid);
% end
% [rs,rv]=perl('makejfromtemplate.pl',sprintf('%s%f',sprintf('%f,',p(1:(end-1))),p(end)));

function swobj = run_mcphase(swobj, tstinfo)
    global mcphasedir temperature magfield emin emax;
    nl = sprintf('\n'); if ispc; sp = '\'; psp = ';'; else; sp = '/'; psp = ':'; end
    % Runs the mean field program in order to get a magnetic moment and spin direction at the 
    %    field, temperature point to be calculated (default: 2K, 0T)
    % To set different values, define global variables: temperature and/or magfield. E.g.:
    %   >> global temperature magfield; temperature=10; magfield=[1 1 0]/sqrt(2);
    if exist('results','dir')~=7; mkdir('results'); end
    if exist('run_env.pl','file')~=2  
      runmulti=['#!/usr/bin/perl' nl ...
                '' nl ...
                'if ($#ARGV<1) { $cmdline = "mcdispit -t"; } else { $cmdline = $ARGV[1]; }' nl ...
                'if ($#ARGV<0) { $mcphasedir = "C:/mcphas5_3" } else { $mcphasedir = $ARGV[0]; }' nl ...
                '$ENV{"MCPHASE_DIR"} = $mcphasedir;' nl ...
                '$ENV{"PATH"} = $mcphasedir."/bin' psp '".$ENV{"PATH"};' nl ...
                'system("%s $cmdline");' nl];
      if ispc; multicmdprefix = 'start /B /WAIT /BELOWNORMAL'; else multicmdprefix=''; end
      fid=fopen('run_env.pl','w'); fprintf(fid,runmulti,multicmdprefix); fclose(fid);
    end
    global no_mf
    if isempty(no_mf) || ~no_mf
        [rs,rv]=perl('run_env.pl',mcphasedir,'mcphasit -v'); if(length(rs)>8e3); rs=rs((end-8e3):end); end
        if rv~=0; error(sprintf('Failed to run meanfield program mcphasit. Output is:\n%s',rs)); end
    end
    spinscmd = sprintf('spins -f results%smcphas.mf %.4g %.4g %.4g %.4g > mcdisp.mf',sp,temperature,magfield);
    [rs,rv]=perl('run_env.pl',mcphasedir,spinscmd);
    if rv~=0; error(sprintf('Failed to run spins program to extract magnetic moment. Output is:\n%s',rs)); end
    global scale_mf
    if ~isempty(scale_mf) && any(scale_mf); rescale_mf(scale_mf); end
    if exist('mcdisp.par','file')~=2; fid=fopen('mcdisp.par','w'); 
    fprintf(fid,[...
      '#<!--mcdisp.mcdisp.par>' nl ...
      '#*********************************************************************' nl ...
      '# mcdisp - program to calculate the dispersion of magnetic excitations' nl ...
      '# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751' nl ...
      '#*********************************************************************' nl ...
      '#!ki=99999' nl ...
      '#!emin=%f' nl ...
      '#!emax=%f' nl ...
      '#!calculate_magmoment_oscillation=0  creates mcdisp.qem' nl ...
      '#!calculate_phonon_oscillation=0  creates mcdisp.qep' nl ...
      '#!calculate_spinmoment_oscillation=0  creates mcdisp.qes' nl ...
      '#!calculate_orbmoment_oscillation=0  creates mcdisp.qeo' nl ...
      '#!calculate_chargedensity_oscillation=0  creates mcdisp.qee' nl ...
      '#!calculate_spindensity_oscillation=0  creates mcdisp.qsd' nl ...
      '#!calculate_orbmomdensity_oscillation=0  creates mcdisp.qod' nl ...
      '#!outS=1' nl ...
      '#!hklfile=mcdisp.hkl' nl], emin, emax);
    fclose(fid); end
    if exist('mcdisp.hkl','file')~=2; fid=fopen('mcdisp.hkl','w'); fprintf(fid,'0.5 0.5 0.5\n'); fclose(fid); end
    S2 = 2*max(swobj.matom.S); e2 = max(abs([emin emax])); e1 = -e2;
    [rs,rv]=perl('run_env.pl',mcphasedir,sprintf('mcdispit -c -minE %f -maxE %f -max %f', e1, e2, S2*10)); if(length(rs)>8e3); rs=rs((end-8e3):end); end
    %[rs,rv]=perl('run_env.pl',mcphasedir,'mcdispit -c'); if(length(rs)>8e3); rs=rs((end-8e3):end); end
    if rv~=0; error(sprintf('Failed to run mcdispit program to calculate single ion levels. Output is:\n%s',rs)); end
    [rs,rv]=perl([mcphasedir sp 'bin' sp 'range.pl'],'7','0.01','1e99','results/mcdisp.trs');
    [rs,rv]=perl([mcphasedir sp 'bin' sp 'range.pl'],'6',num2str(emin),num2str(emax),'results/mcdisp.trs');
    if rv~=0; error(sprintf('Failed to run range program to remove zero weight modes. Output is:\n%s',rs)); end
    fid = fopen('results/mcphas.sps'); while true; ll=strip(fgetl(fid)); if ll(1)~='#'; break; end; end
    tst = textscan(fid, '%f'); fclose(fid); tst = tst{1};
    if numel(tst) == numel(tstinfo.tst)
        [na, nb, nc] = size(tstinfo.tst); S0 = []; tstm = zeros(na,nb,nc);
        nblk = tstinfo.n_at * tstinfo.n_comp; nExt = tstinfo.nExt;
        for ia = 1:na; for ib = 1:nb; for ic = 1:nc; 
            tstm(ia, ib, ic) = tst((ic-1)*na*nb+(ia-1)*nb+ib);
        end; end; end;
        for ia = 1:nExt(1); for ib = 1:nExt(2); for ic = 1:nExt(3); 
            S0 = [S0 reshape(tstm(((ib-1)*nblk+1):(ib*nblk), ia, ic), [tstinfo.n_comp tstinfo.n_at])];
        end; end; end;
        if tstinfo.n_comp == 6
            S0 = S0([1 3 5],:) + 2*S0([2 4 6],:);
        end
        swobj.genmagstr('mode', 'direct', 'nExt', nExt, 'S', S0);
    else
        try; swobj.genmagstr('mode', 'direct', 'nExt', swobj.magstr.N_ext, 'S', repmat([0 0 1]',1,size(swobj.magstr.S,1))); end;
    end
end
    
function [w, s, sab] = run_mcdisp(h, k, l, infiles)
    global mcphasedir;
    % Rounds to 4 d.p.
    dp=4; rhkl=str2num(sprintf(sprintf('%%.%df %%.%df %%.%df\\n',dp,dp,dp),[h(:) k(:) l(:)]'));
    uhkl = unique(rhkl,'rows');
    nl = sprintf('\n'); if ispc; sp = '\'; psp = ';'; else; sp = '/'; psp = ':'; end
    qeiformatstr = '%f %f %f %f %f %f %f %f %s %s %s %s  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
    % If user wants to use more than one processor
    global numproc;
    
    if ~isempty(numproc) && numproc>1
      % Splits the input between <numproc> input files in different directories
      hklcnt = 1; hklblock = max([1 floor(size(uhkl,1)/numproc)]); hklnextcnt = hklcnt+hklblock;
      for np = 1:numproc; 
        if np == numproc && (hklnextcnt < size(uhkl,1)); hklnextcnt=size(uhkl,1); end
        idxmem{np} = hklcnt:hklnextcnt;
        hklcnt=hklnextcnt+1; hklnextcnt = hklcnt+hklblock-1; 
        if (hklnextcnt>size(uhkl,1)); hklnextcnt=size(uhkl,1); end
      end
      idxmem = idxmem(find(~cellfun(@isempty, idxmem))); npp = numel(idxmem);
      w = {}; s = {}; sab = {};
      for np = 1:npp
        dirname = ['mcdisp_temp-p' num2str(np)];
        if exist(dirname,'dir')~=7; mkdir(dirname); end
        if exist([dirname sp 'results'],'dir')~=7; mkdir([dirname sp 'results']); end
        for nf = 1:length(infiles); copyfile(infiles{nf},[dirname sp infiles{nf}]); end
        copyfile('results/mcdisp.trs',[dirname sp 'results/mcdisp.trs']);
        copyfile('mcdisp.par',[dirname sp 'mcdisp.par']);
        fid = fopen([dirname sp 'mcdisp.hkl'],'w'); 
        if np > 1
            % McPhase uses the first two Q-vectors to determine the UVW coordinate system
            % in which the Sab spin susceptibility tensor is calculated.
            % Thus we ensure here that these are the same for all processes so the Sab
            % are all in the same coordinate system.
            fprintf(fid,sprintf('%%.%df %%.%df %%.%df\\n',dp,dp,dp),uhkl(idxmem{1}(1:2),:)');
        end
        fprintf(fid,sprintf('%%.%df %%.%df %%.%df\\n',dp,dp,dp),uhkl(idxmem{np},:)');
        fclose(fid);
      end
      % Run parallel processes.
      if exist('run_multi.pl','file')~=2  
        % Creates a perl file in order to use the Perl Fork() function to spawn child processes
        runmulti=[ '#!/usr/bin/perl' nl ...
                   '' nl ...
                   'if ($#ARGV<3) { $cmdline = "mcdispit -t"; } else { $cmdline = $ARGV[3]; }' nl ...
                   'if ($#ARGV<2) { $mcphasedir = "C:/mcphas5_3" } else { $mcphasedir = $ARGV[2]; }' nl ...
                   'if ($#ARGV<1) { $numproc = 2; } else { $numproc = $ARGV[1]; }' nl ...
                   'if ($#ARGV<0) { $dirprefix = "mcd-p"; } else { $dirprefix = $ARGV[0]; }' nl ...
                   '$ENV{"MCPHASE_DIR"} = $mcphasedir;' nl ...
                   '$ENV{"PATH"} = $mcphasedir."/bin' psp '".$ENV{"PATH"};' nl ...
                   '' nl ...
                   '# Parallel fork code shameless stolen from: http://hell.jedicoder.net/?p=82' nl ...
                   '' nl ...
                   'for(1..$numproc) {' nl ...
                   '  my $pid = fork();' nl ...
                   '  if ($pid) {' nl ...
                   '    push(@children, $pid); # parent' nl ...
                   '  } ' nl ...
                   '  elsif ($pid == 0) {' nl ...
                   '    system("cd $dirprefix$_ && %s $cmdline");' nl ...
                   '    exit(0);' nl ...
                   '  } ' nl ...
                   '  else {' nl ...
                   '    die ''couldn\\''t fork: $!\\n'';' nl ...
                   '  }' nl ...
                   '}' nl ...
                   '' nl ...
                   'foreach (@children) { waitpid($_, 0); }'];
        if ispc; multicmdprefix = 'start /B /WAIT /BELOWNORMAL'; else multicmdprefix=''; end
        fid=fopen('run_multi.pl','w'); fprintf(fid,runmulti,multicmdprefix); fclose(fid);
      end
      [rs,rv]=perl('run_multi.pl','mcdisp_temp-p',num2str(npp),mcphasedir,'mcdispit -t');
      if(length(rs)>8e3); rs=rs((end-8e3):end); end
      if rv~=0; error(sprintf('Failed to run mcdispit program to calculate dispersion. Output is:\n%s',rs)); end
      % Combines all the output files
      for np = 1:numproc
        filename = ['mcdisp_temp-p' num2str(np) sp 'results' sp 'mcdisp.qei'];
        [rs,rv]=perl([mcphasedir sp 'bin' sp 'range.pl'],'10','-0.0001','1e99',filename);
        if rv~=0; error(sprintf('Failed to run range program to remove uncalculated modes. Output is:\n%s',rs)); end
        [ha,hb,hc,t,qh,qk,ql,qmod,en,cs,~,~,...
         sr{1},si{1},sr{2},si{2},sr{3},si{3},sr{4},si{4},sr{5},si{5},sr{6},si{6},sr{7},si{7},sr{8},si{8},sr{9},si{9}] = ...
            textread(filename,qeiformatstr,'commentstyle','shell');
        if isempty(ha); error(sprintf('Failed to calculate dispersion. Output was %s',rs)); end
        en=cell2mat(cellfun(@str2num,en,'UniformOutput',0)); cs=cell2mat(cellfun(@str2num,cs,'UniformOutput',0));
        if(length(cs)~=length(en)); error(sprintf('Failed to calculate dispersion. Output was %s',rs)); end
        qhkl = [qh(:) qk(:) ql(:)]; uqhkl = unique(qhkl,'rows');
        nchunk = size(qhkl,1)/size(uqhkl,1);
        if mod(nchunk, 1) ~= 0
            error(sprintf(['Number of modes returned by McPhase is inconsistent over different Q points.\n' ...
                   'This is likely due to some modes having imaginary eigenvalues or falling outside the energy range.\n' ...
                   'Please adjust the exchange / SIA parameters and/or change the emin/emax parameters and try again.']));
        end
        ismem = zeros([size(rhkl,1) nchunk]); 
        for iq=1:length(rhkl); iin=find(ismember(qhkl,rhkl(iq,:),'rows')); if ~isempty(iin); ismem(iq,:)=iin; end; end
        out_sz = [size(rhkl,1) size(qhkl,1)/size(uqhkl,1)];
        if isempty(w); w = zeros(out_sz); end
        if isempty(s); s = zeros(out_sz); end
        if isempty(sab); sab = zeros([out_sz 3 3]); blksz = prod(out_sz); end
        id_mem = find(ismem); ismem = ismem(id_mem);
        w(id_mem) = en(ismem); s(id_mem) = cs(ismem);
        for jj = 1:3
            for ii = 1:3
                idx = (jj-1)*3+ii;
                rp = cell2mat(cellfun(@str2num, sr{idx}, 'UniformOutput', 0));
                ip = cell2mat(cellfun(@str2num, si{idx}, 'UniformOutput', 0));
                sab((idx-1)*blksz + id_mem) = rp(ismem) + 1j*ip(ismem);
            end
        end
      end
    else
      % Writes the list of (hkl) to a single mcdisp input file [mcdisp.par]
      fid = fopen('mcdisp.hkl','w');
      fprintf(fid,sprintf('%%.%df %%.%df %%.%df\\n',dp,dp,dp),uhkl');
      fclose(fid);
      [rs,rv]=perl('run_env.pl',mcphasedir,'mcdispit -t'); if(length(rs)>8e3); rs=rs((end-8e3):end); end
      if rv~=0; error(sprintf('Failed to run mcdispit program to calculate dispersion. Output is:\n%s',rs)); end
      [ha,hb,hc,t,qh,qk,ql,qmod,en,cs,~,~,...
       sr{1},si{1},sr{2},si{2},sr{3},si{3},sr{4},si{4},sr{5},si{5},sr{6},si{6},sr{7},si{7},sr{8},si{8},sr{9},si{9}] = ...
          textread(['results' sp 'mcdisp.qei'],qeiformatstr,'commentstyle','shell');
      if isempty(ha); error(sprintf('Failed to calculate dispersion. Output was %s',rs)); end
      en=cell2mat(cellfun(@str2num,en,'UniformOutput',0)); cs=cell2mat(cellfun(@str2num,cs,'UniformOutput',0));
      if(length(cs)~=length(en)); error(sprintf('Failed to calculate dispersion. Output was %s',rs)); end
      qhkl = [qh(:) qk(:) ql(:)]; for iq=1:length(rhkl); ismem(iq,:)=find(ismember(qhkl,rhkl(iq,:),'rows')); end
      w = en(ismem);
      s = cs(ismem);
      sab = zeros(size(w,1),size(w,2),3,3);
      for jj = 1:3
          for ii = 1:3
              rp = cellfun(@str2num, sr{(jj-1)*3+ii}, 'UniformOutput', 0);
              ip = cellfun(@str2num, si{(jj-1)*3+ii}, 'UniformOutput', 0);
              id0 = cellfun(@isempty, rp); rp(id0) = {[0]}; rp = cell2mat(rp);
              id0 = cellfun(@isempty, ip); ip(id0) = {[0]}; ip = cell2mat(ip);
              sab(:,:,ii,jj) = rp(ismem) + 1j*ip(ismem);
          end
      end
    end
    w = permute(w, [2 1]);
    s = permute(s, [2 1]);
    sab = permute(sab, [3 4 2 1]);
end


function [infiles, tstinfo] = generate_mcphase_inputs(swobj)
    if isfield(swobj.cache, 'mcphase')
        if isfield(swobj.cache.mcphase, 'pointcharge')
            [sipf, n_comp] = generate_sipf_pc(swobj);
        else
            [sipf, n_comp] = generate_sipf_cf(swobj);
        end
    else
        % Check for symmetry modified SIA tensors
        [SS, SI, RR] = intmatrix(swobj,'plotmode',true,'extend',false,'sortDM',false,'zeroC',false,'nExt',[1 1 1]);
        anisos = reshape(SI.aniso, 9, size(SI.aniso, 3))'; [ua,ia,ja] = unique(anisos, 'rows');
        n_at_uniq = unique(swobj.matom.idx);
        if numel(ia) == numel(n_at_uniq)
            [sipf, n_comp] = generate_sipf_sia(swobj);
        else
            [sipf, n_comp] = generate_sipf_sia_sym(swobj);
        end
    end
    for ii = 1:numel(sipf)
        infiles{ii} = [sipf{ii} '.sipf'];
    end
    spinw2mcphasj(swobj, sipf, n_comp);
    tstinfo = generate_mcphas_ini_tst(swobj, sipf, n_comp);
    infiles = {infiles{:} 'mcphas.j' 'mcphas.ini' 'mcphas.tst' 'mcdisp.par' 'mcdisp.hkl' 'mcdisp.mf'};
end

function [sipf, n_comp] = generate_sipf_pc(swobj)
    pcstruct = swobj.cache.mcphase.pointcharge;
    if isfield(pcstruct, 'cif')
        cifstr = pcstruct.cif;
    else
        cifstr = sprintf([...
        '_cell_length_a                     %f\n' ...
        '_cell_length_b                     %f\n' ...
        '_cell_length_c                     %f\n' ...
        '_cell_angle_alpha                  %f\n' ...
        '_cell_angle_beta                   %f\n' ...
        '_cell_angle_gamma                  %f\n' ...
        '_symmetry_space_group_name_H-M     %s\n' ...
        'loop_\n' ...
        '_atom_site_label\n' ...
        '_atom_type_oxidation_number\n' ...
        '_atom_site_fract_x\n' ...
        '_atom_site_fract_y\n' ...
        '_atom_site_fract_z\n'], swobj.lattice.lat_const, swobj.lattice.angle * 180 / pi, swobj.lattice.label);
        if isfield(pcstruct, 'atoms')
            atoms = pcstruct.atoms;
            for ii = 1:numel(atoms)
                cifstr = sprintf('%s%s %f %f %f %f\n', cifstr, atoms(ii).label, atoms(ii).valence, atoms(ii).r);
            end
        else
            for ii = 1:numel(swobj.unit_cell.label)
                label = swobj.unit_cell.label{ii};
                oxy = swobj.unit_cell.ox(ii);
                rr = swobj.unit_cell.r(:, ii);
                cifstr = sprintf('%s%s %f %f %f %f\n', cifstr, label, oxy, rr);
            end
        end
    end
    fid = fopen('mcphas.cif', 'w');
    fprintf(fid, cifstr);
    fclose(fid);
    nl = sprintf('\n'); if ispc; sp = '\'; psp = ';'; else; sp = '/'; psp = ':'; end
    global mcphasedir;
    setenv('MCPHASE_DIR', mcphasedir);
    if isempty(strfind(getenv('PATH'), mcphasedir))
        setenv('PATH',[getenv('PATH') psp mcphasedir sp 'bin']);
    end
    if isempty(strfind(getenv('PERL5LIB'), mcphasedir))
        setenv('PERL5LIB',[getenv('PERL5LIB') psp mcphasedir sp 'bin'])
    end
    [rs,rv]=perl([mcphasedir sp 'bin' sp 'cif2mcphas.pl'],'-pc',num2str(pcstruct.maxdist),'-sp','mcphas.cif');
    if rv~=0; error(sprintf('Failed to run cif2mcphas program to calculate point charges. Output is:\n%s',rs)); end
    fid = fopen('mcphas.j', 'r');
    sipfstr = fscanf(fid, '%s');
    fclose(fid);
    id0 = strfind(sipfstr, 'sipffilename=') + 13;
    n_comp = 3;
    for ii = 1:numel(id0)
        id1 = strfind(sipfstr(id0(ii):end), '.sipf');
        sipf{ii} = sipfstr(id0(ii):(id0(ii)+id1(1)-2));
        fid = fopen(sprintf('%s.sipf', sipf{ii}));
        modul = fscanf(fid, '#!MODULE=%s');
        fclose(fid);
        if strcmp(modul, 'ic1ion')
            n_comp = 6;
        end
    end
    if isfield(pcstruct, 'atoms') && isfield(pcstruct.atoms(1), 'charges')
        % Convert valence values to labels
        for ii = 1:numel(sipf)
            for ia = 1:numel(atoms)
                %disp(['|' sprintf('% 8.4g ',atoms(ia).valence) '|' sprintf('% 8.4g ',atoms(ia).charges), sprintf('results%s%s.pc',sp,sipf{ii})])
                [rs,rv] = perl([mcphasedir sp 'bin' sp 'substitute.pl'], ...
                        sprintf('% 8.4g ',atoms(ia).valence), atoms(ia).label, sprintf('results%s%s.pc',sp,sipf{ii}));
                if rv~=0; warning(sprintf('Failed to set point charges. Output is:\n%s',rs)); end
            end
        end
        % Then convert labels to desired effective charges
        % This avoids the error where the desired charge of one ion is the same as the valence of another.
        for ii = 1:numel(sipf)
            for ia = 1:numel(atoms)
                [rs,rv] = perl([mcphasedir sp 'bin' sp 'substitute.pl'], ...
                        atoms(ia).label, sprintf('% 8.4g ',atoms(ia).charges), sprintf('results%s%s.pc',sp,sipf{ii}));
                if rv~=0; warning(sprintf('Failed to set point charges. Output is:\n%s',rs)); end
            end
        end
        [rs,rv]=perl([mcphasedir sp 'bin' sp 'cif2mcphas.pl'],'-pc',num2str(pcstruct.maxdist),'-rp','mcphas.cif');
    end
    [rs,rv]=perl([mcphasedir sp 'bin' sp 'substitute.pl'],'B00','"#B00"','*.sipf');
    [rs,rv]=perl([mcphasedir sp 'bin' sp 'substitute.pl'],'L00','"#L00"','*.sipf');
    if isfield(pcstruct, 'atoms') && isfield(pcstruct.atoms(1), 'zeta')
        for ia = 1:numel(atoms)
            [rs,rv]=perl([mcphasedir sp 'bin' sp 'setvariable.pl'],'zeta',num2str(atoms(ia).zeta),[atoms(ia).label '_*.sipf']);
        end
    end
end

function [sipf, n_comp] = generate_sipf_cf(swobj, parstruct)
    nl = sprintf('\n'); if ispc; sp = '\'; else sp = '/'; end
    fstr = [...
    '#!MODULE=%s' nl ...
    'IONTYPE=%s' nl ...
    'DWF=0' nl ...
    'SCATTERINGLENGTHREAL=%f' nl ...
    'SCATTERINGLENGTHIMAG=%f' nl ...
    'units=meV' nl];
    blen = swobj.unit_cell.b(1, :) / 10;
    ffstr = generate_ff(swobj);
    if exist('parstruct', 'var')
        u_ions = parstruct.ions;
        blen = blen(parstruct.idx);
        ffstr = ffstr(parstruct.idx);
    else
        [u_ions, id_uions, id_ions] = unique(swobj.unit_cell.label);
        blen = blen(id_uions);
        ffstr = ffstr(id_uions);
        sipf = swobj.unit_cell.label(id_ions(swobj.matom.idx));
    end
    for io = 1:numel(u_ions)
        if exist('parstruct', 'var') || isfield(swobj.cache.mcphase, u_ions{io})
            if exist('parstruct', 'var')
                cfstruct = parstruct.cfpars(io);
            else
                cfstruct = swobj.cache.mcphase.(u_ions{io});
            end
            non_magnetic = false;
            if isfield(cfstruct, 'non_magnetic') && cfstruct.non_magnetic
                module = 'kramer';
                n_comp = 3;
                fstr = strrep(fstr, 'IONTYPE', '#IONTYPE');
                if ~isfield(cfstruct, 'iontype')
                    cfstruct.iontype = 'S=0';
                end
                ffstr{io} = sprintf('GJ=2\nA=0\nB=0\nC=0\n%s\n', ffstr{io});
                non_magnetic = true;
                ffstr{io} = strrep(regexprep(ffstr{io},'[0-9]','0'), '00.', ' 0.');
            elseif ~isfield(cfstruct, 'module')
                module = 'so1ion';
                n_comp = 3;
            else
                module = cfstruct.module;
                if strcmp(cfstruct.module, 'ic1ion')
                    n_comp = 6;
                else
                    n_comp = 3;
                end
            end
            sipfstr = sprintf(fstr, module, cfstruct.iontype, real(blen(io)), imag(blen(io)));
            pars = fieldnames(cfstruct);
            for ii = 1:numel(pars)
                if ~strcmp(pars{ii}, 'iontype') && ~strcmp(pars{ii}, 'module') && ~non_magnetic
                    sipfstr = [sipfstr sprintf(['%s = %f' nl], pars{ii}, cfstruct.(pars{ii}))];
                end
            end
            fid = fopen([u_ions{io} '.sipf'], 'w');
            fprintf(fid, sipfstr);
            fprintf(fid, ffstr{io});
            fclose(fid);
        else
            error(sprintf('Ion label %s not found in mcphase field', u_ions{io}));
        end
    end
end

function [sipf, n_comp] = generate_sipf_sia(swobj)
    [u_mat, id_mat, r_idx] = unique(swobj.single_ion.aniso);
    idx = swobj.matom.idx(id_mat);
    ions = swobj.unit_cell.label(idx);
    for ii = 1:numel(ions)
        ions{ii} = sprintf('%s_%i', ions{ii}, ii);
        if u_mat(ii) > 0
            dmat = swobj.matrix.mat(:,:,u_mat(ii));
        else
            dmat = zeros(3);
        end
        Sval = swobj.matom.S(id_mat(ii));
        cfstruct(ii).iontype = sprintf('S=%2.1f', Sval);
        cfstruct(ii).gJ = 2;
        if Sval < 0.1
            cfstruct(ii).non_magnetic = true;
        else
            cfstruct(ii).Dx2 = dmat(1, 1);
            cfstruct(ii).Dy2 = dmat(2, 2);
            cfstruct(ii).Dz2 = dmat(3, 3);
            cfstruct(ii).B21 = dmat(1, 3);
            cfstruct(ii).B21S = dmat(2, 3);
            cfstruct(ii).B22S = dmat(1, 2)/2;
            cfstruct(ii).non_magnetic = false;
        end
    end
    parstruct.idx = idx;
    parstruct.ions = ions;
    parstruct.cfpars = cfstruct;
    generate_sipf_cf(swobj, parstruct);
    sipf = ions(r_idx);
    n_comp = 3;
end

function [sipf, n_comp] = generate_sipf_sia_sym(swobj)
    if any(swobj.single_ion.aniso > 0)
        [SS, SI, RR] = intmatrix(swobj,'plotmode',true,'extend',false,'sortDM',false,'zeroC',false,'nExt',[1 1 1]);
        swobj2 = copy(swobj);
        Svals = double(swobj.atom.mag); Svals(find(Svals==1)) = swobj.matom.S;
        lab = {}; ox = []; occ = []; b = []; ff = []; A = []; Z = []; biso = []; cols = [];
        kk = 1;
        cell2 = swobj2.unit_cell;
        for ii = 1:numel(swobj2.unit_cell.S)
            n_at = numel(find(swobj.atom.idx == ii));
            for jj = 1:n_at
                lab = [lab sprintf('%s_%d', swobj.unit_cell.label{ii}, kk)];
                kk = kk + 1;
            end
            ox = [ox repmat(swobj.unit_cell.ox(ii), 1, n_at)];
            occ = [occ repmat(swobj.unit_cell.occ(ii), 1, n_at)];
            b = [b repmat(swobj.unit_cell.b(:,ii), 1, n_at)];
            ff = cat(3, ff, repmat(swobj.unit_cell.ff(:,:,ii), 1, 1, n_at));
            A = [A repmat(swobj.unit_cell.A(ii), 1, n_at)];
            Z = [Z repmat(swobj.unit_cell.Z(ii), 1, n_at)];
            biso = [biso repmat(swobj.unit_cell.biso(ii), 1, n_at)];
            cols = [cols repmat(swobj.unit_cell.color(:,ii), 1, n_at)];
        end
        swobj2.unit_cell = struct('r', swobj.atom.r, 'S', Svals, 'label', {lab}, 'color', cols, ...
            'ox', ox, 'occ', occ, 'b', b, 'ff', ff, 'A', A, 'Z', Z, 'biso', biso);
        swobj2.lattice.sym = swobj2.lattice.sym(:,:,1);
        % Recompute SIA
        [ua, ia, ja] = unique(round(transpose(reshape(SI.aniso, 9, size(SI.g,3)))*1e9)/1e9, 'rows');
        mat_id = swobj2.single_ion.aniso(find(swobj2.single_ion.aniso>0));
        cc = swobj2.matrix.color(:, mat_id(1));
        ua = transpose(ua);
        i0 = size(swobj2.matrix.mat, 3);
        n_aniso = size(ua, 2);
        swobj2.matrix.mat = cat(3, swobj2.matrix.mat, reshape(ua, 3, 3, n_aniso));
        swobj2.matrix.label = cat(2, swobj2.matrix.label, cellfun(@(ii)sprintf('K%02d',ii), num2cell(1:n_aniso), 'UniformOutput', false));
        swobj2.matrix.color = cat(2, swobj2.matrix.color, repmat(cc, 1, n_aniso));
        swobj2.single_ion.aniso = ja(:)' + i0;
    end
    [sipf, n_comp] = generate_sipf_sia(swobj2);
end

function ffstr =  generate_ff(swobj)
    nl = sprintf('\n'); if ispc; sp = '\'; else sp = '/'; end
    for ii = 1:size(swobj.unit_cell.ff, 3)
        ff = swobj.unit_cell.ff(1,[1:6 end],ii);
        ffstr{ii} = sprintf([...
        'FFj0A=%7.4f FFj0a=%7.4f FFj0B=%7.4f FFj0b=%7.4f FFj0C=%7.4f FFj0c=%7.4f FFj0D=%7.4f' nl ...
        'FFj2A= 0.0000 FFj2a= 0.0000 FFj2B= 0.0000 FFj2b= 0.0000 FFj2C= 0.0000 FFj2c= 0.0000 FFj2D= 0.0000' nl ...
        'FFj4A= 0.0000 FFj4a= 0.0000 FFj4B= 0.0000 FFj4b= 0.0000 FFj4C= 0.0000 FFj4c= 0.0000 FFj4D= 0.0000' nl], ff);
    end
end

function spinw2mcphasj(swobj, sipf, n_comp)
    if n_comp ~= 3 && n_comp ~= 6
        n_comp = 3;
    end
    n_at = numel(swobj.matom.S);
    jfile = 'mcphas.j';
    
    header = sprintf([...
    '#<!--mcphase.mcphas.j-->\n' ...
    '#***************************************************************\n' ...
    '# Lattice and Exchange Parameter file for\n' ...
    '# mcphas version 5.2\n' ...
    '# - program to calculate static magnetic properties\n' ...
    '# reference: M. Rotter JMMM 272-276 (2004) 481\n' ...
    '# mcdisp version 5.2\n' ...
    '# - program to calculate the dispersion of magnetic excitations\n' ...
    '# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n' ...
    '#***************************************************************\n' ...
    '#\n' ...
    '# SpinW model\n' ...
    '#\n' ...
    '# Lattice Constants (A)    Spacegroup=%s\n' ...
    '#\n' ...
    '#! a= %f b= %f c= %f  alpha= %f beta= %f gamma= %f\n' ...
    '#\n' ...
    '#! r1a=   1 r2a= 0 r3a=  0\n' ...
    '#! r1b=   0 r2b= 1 r3b=  0   primitive lattice vectors [a][b][c]\n' ...
    '#! r1c=   0 r2c= 0 r3c=  1\n' ...
    '#\n' ...
    '#! nofatoms= %i  nofcomponents=%i  number of atoms in primitive unit cell/number of components of each spin\n' ...
    '#*************************************************************************'], ...
    swobj.lattice.label, swobj.lattice.lat_const, swobj.lattice.angle * 180 / pi, n_at, n_comp);
    
    atom_header = '#! da= %f [a] db= %f [b] dc= %f [c] nofneighbours=%i diagonalexchange=2 sipffilename= %s.sipf';
    comp_header{3} = '#! symmetricexchange=0 indexexchange= JaJa JbJb JcJc JaJb JbJa JaJc JcJa JbJc JcJb\n';
    comp_header{6} = '#! symmetricexchange=0 indexexchange= JaJa JcJc JeJe JaJc JcJa JaJe JeJa JcJe JeJc\n';
    
    fid = fopen(jfile, 'w');
    fprintf(fid, '%s\n', header);

    [SS, SI, RR] = swobj.intmatrix('extend', false);
    if ~isempty(SS.bq)
        warning('Biquadratic interactions defined - McPhase does not support biquadratic interactions. These will be ignored.');
    end

    dl = double(swobj.coupling.dl);
    at1 = swobj.coupling.atom1;
    at2 = swobj.coupling.atom2;
    id0 = [dl; at1; at2];
    dr = swobj.matom.r(:,at2) + dl - swobj.matom.r(:,at1);
    len = sqrt(sum((swobj.basisvector * dr).^2, 1));
    atstr = cell(1, n_at);
    for iat = 1:n_at
        % Get list of neighbouring bonds and interactions matrices
        idx = find(at1 == iat | at2 == iat);
        pm = ones(1, numel(idx)); pm(at1(idx) ~= iat) = -1;
        dr1 = repmat(pm, 3, 1) .* dr(:, idx);
        % For same atom in adjacent cells, use also -dl
        ids = find(at1 == iat & at2 == iat);
        [ln1, isr] = sort([len(idx) len(ids)]);
        pm = [pm ids*0+1]; pm = pm(isr);
        dr1 = [dr1 -dr(:, ids)]; dr1 = round(dr1(:, isr)*1e6)/1e6;
        idd = [idx ids]; idd = idd(isr);
        n_neigh = numel(ln1);
        jmat = zeros(9, n_neigh);
        if isfield(swobj.cache, 'mcphase') && isfield(swobj.cache.mcphase, 'self_interaction')
            fprintf(fid, [atom_header '\n'], swobj.matom.r(:,iat), n_neigh + 1, sipf{iat});
            fprintf(fid, comp_header{n_comp});
            fprintf(fid, '% f\t', [0 0 0 swobj.cache.mcphase.self_interaction 0 0 0 0 0 0]);
            fprintf(fid, '\n');
        else
            fprintf(fid, [atom_header '\n'], swobj.matom.r(:,iat), n_neigh, sipf{iat});
            fprintf(fid, comp_header{n_comp});
        end
        %idm = swobj.coupling.mat_idx(:, idd);
        for inn = 1:n_neigh
            %jmat(:, inn) = reshape(sum(swobj.matrix.mat(:, :, idm((idm(:, inn) > 0), inn)), 3), 9, 1);
            %jm = -jmat([1 5 9 2 4 3 7 6 8], inn);
            %% Reverse sign for j->i opposite sense pair
            %if pm(inn) == -1
            %    jm([4:9]) = -jm([4:9]);
            %end
            idall = find(ismember(SS.all(1:5,:)', id0(:, idd(inn))', 'rows'));
            % Remove the biquadratic interactions
            idall(find(SS.all(15, idall)==1)) = [];
            % Sum all the defined matrices
            if ~isempty(idall)
                jm = -SS.all([1 5 9 2 4 3 7 6 8]+5, idall); % McPhase uses opposite sign convention to SpinW
                if size(jm, 2) > 1
                    jm = sum(jm, 2);
                end
                % Reverse sign for j->i opposite sense pair
                if pm(inn) == -1
                    jm([4:9]) = -jm([4:9]);
                end
            else
                jm = zeros(9, numel(idd(inn)));
            end
            % Adds the dipolar interaction if it exists
            if ~isempty(SS.dip)
                jm = jm - SS.dip([1 5 9 2 4 3 7 6 8]+5, idd(inn));
            end
            fprintf(fid, '% f\t', [dr1(:, inn); jm]);
            fprintf(fid, '\n');
        end
        %atstr{iat}{1} = sprintf(atom_header, swobj.matom.r(:,iat), n_neigh, sipf{iat});
        %atstr{iat}{2} = num2str([dr1; -jmat([1 5 9 2 4 3 7 6 8], :)]');
        %disp(atstr{iat}{1}); disp(atstr{iat}{2});
    end
    fclose(fid);
end

function tstinfo = generate_mcphas_ini_tst(swobj, sipf, n_comp)
    n_at = numel(sipf);
    tstheader = sprintf('#! nofatoms=%i\n#! nofcomponents=%i\n', n_at, n_comp);
    magStr = swobj.magstr;
    if any(abs(magStr.k-round(magStr.k)) > 1e-4)
        if any(magStr.N_ext ~= [1 1 1])
            error('Cannot handle both incommensurate and extended unit cell cases');
        end
        k_mc = magStr.k;
        nExt = floor(1 ./ k_mc);
        nExt(isinf(nExt)) = 1;
        nExt(isnan(nExt)) = 1;
        S0 = swobj.magstr('nExt', nExt).S;
    else
        k_mc = 1 ./ magStr.N_ext;
        nExt = magStr.N_ext;
        S0 = magStr.S;
    end
    if n_comp == 6
        S0 = [S0; S0*2]; 
        S0 = S0([1 4 2 5 3 6],:);
    end
    tst = zeros([n_comp*n_at 1 1] .* nExt([2 1 3]));
    for ic = 1:nExt(3)
        for ib = 1:nExt(2)
            b0 = (ib-1) * n_comp * n_at + 1;
            b1 = ib * n_comp * n_at;
            for ia = 1:nExt(1)
                tstcol = [];
                for isp = 1:n_at
                    tstcol = [tstcol; S0(:,(ic-1)*nExt(1)*nExt(2)*n_at+(ib-1)*nExt(1)*n_at+(ia-1)*n_at+isp)];
                end
                tst(b0:b1,ia,ic) = tstcol;
            end
        end
    end
% {
    fid = fopen('mcphas.tst', 'w');
    fprintf(fid, tstheader);
    for ic = 1:size(tst,3)
        for ia = 1:size(tst,1)
            for ib = 1:size(tst,2)
                fprintf(fid, ' %6.4f', tst(ia, ib, ic));
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
% }
    global magfield temperature; 
    tt = temperature;
    hh = norm(magfield);
    if hh > 0
        hd = magfield ./ hh;
    else
        hd = magfield;
    end
    fid = fopen('mcphas.ini', 'w');
    fprintf(fid, [...
        'xT=1\nxHa=0\nxHb=0\nxHc=0\nxmin=%f\nxmax=%f\nxstep=1\n' ...
        'yT=0\nyHa=%f\nyHb=%f\nyHc=%f\nymin=%f\nymax=%f\nystep=1\n' ...
        'hmin=%f\nhmax=%f\ndeltah=%f\n' ...
        'kmin=%f\nkmax=%f\ndeltak=%f\n' ...
        'lmin=%f\nlmax=%f\ndeltal=%f\n' ...
        'maxqperiod=%d\n' ...
        'maxnofspins=%d\n' ...
        'nofrndtries=0\n' ...
        'maxnofmfloops=1000000\n' ...
        'maxspinchange=10000\n' ...
        'maxstamf=1e-6\n' ...
        'bigstep=5.8\n'], ...
        temperature, temperature, hd, hh, hh, ...
        k_mc(1)-0.05, k_mc(1)+0.01, 0.05, ...
        k_mc(2)-0.05, k_mc(2)+0.01, 0.05, ...
        k_mc(3)-0.05, k_mc(3)+0.01, 0.05, ...
        max(nExt)+1, n_at*prod(nExt)+1);
    fclose(fid);
    tstinfo = struct('n_comp', n_comp, 'n_at', n_at, 'nExt', nExt, 'tst', tst);
end

function rescale_mf(scale_fac)
    if numel(scale_fac) < 2; scale_fac = [1 1]*scale_fac; end
    scale_fac
    fid = fopen('mcdisp.mf'); 
    txt=[]; lnn=1; 
    while ~feof(fid);
        ll = fgets(fid); 
        if ll(1)~='#'; 
            vv=sscanf(ll, '%f'); 
            fstr = repmat('%0.5g ', [1 numel(vv)]);
            if lnn>12; idx=1; else; idx=2; end
            ll=sprintf([fstr '\n'], vv*scale_fac(idx)); 
            lnn=lnn+1; 
        end;
        txt=[txt ll];
    end;
    fclose(fid);
    copyfile('mcdisp.mf', 'mcdisp.mf-orig');
    fid = fopen('mcdisp.mf', 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
    %global mcphasedir;
    %spinscmd = sprintf('factcol 1 %f mcdisp.mf', scale_fac(1)); [rs,rv]=perl('run_env.pl',mcphasedir,spinscmd); if rv~=0; error('bad'); end
    %spinscmd = sprintf('factcol 2 %f mcdisp.mf', scale_fac(1)); [rs,rv]=perl('run_env.pl',mcphasedir,spinscmd); if rv~=0; error('bad'); end
end
