function rg_makespace(inputfolder, outputfolder, loadopts, ncomp_params, clus_params)
% function rg_makespacegen(inputfolder, outputfolder, ncomp_params, clus_params)
%   This function will get called from the command line, with the input and
%   output folder as parameters. Fourier data and parameters are automatically 
%   loaded from the inputfolder (fourfull.mat, space_configs.mat), and the
%   results are saved out to the outputfolder.
% Inputs:
%   inputfolder: input folder with data and config files (fourfull.mat, space_configs.mat).
%   outputfulder: output folder for data saveout
%   loadopts: [dosplit, kthrootnorm, loadfromdisk]
%       flags for loading/preprocessing fourier data
%       - dosplit: 1 for splitting the data in this process, 0 if split
%       - kthrootnorm: 0 for not norming, otherwise norm by power at this root
%       - loadfromdisk: 1 to use filenames as input, 0 to load into memory
%   ncomp_params: [ncompstart, ncompend, ncompstep]; defaults to [10, 10, 1]
%       start, stop, and step size of # of components. If start equals to
%       stop size, then do preset search
%   clus_params: [timreq, nstacks]; defaults to [10*60*60, 1]
%       parameters for torque: computation time required for each init and 
%       number of inits to stack per core

comp_local = 0;

% set default params
if isempty(ncomp_params)
    ncomp_params = [10, 10, 1];
end

if isempty(clus_params)
    clus_params = [10*60*60, 1];
end

% load configs
space_config = load([inputfolder 'space_configs.mat']);

% options for loading/preprocessing datafile
dosplit = loadopts(1);
kthrootnorm = loadopts(2);
loadfromdisk = loadopts(3);

% figure out the control flow for norming/splitting/loading data
if dosplit
    % need to split
    % first load the data
    FC_full_orig = load([inputfolder 'fourfull.mat']);    
    if kthrootnorm
        % need to norm by kth root of eigen values (for spiking data)
        disp(sprintf('Norming Fourier coefficients by %ith root.', kthrootnorm));
        [fourier, scaling] = rmr_normfourierbyroot(FC_full_orig.fourier, kthrootnorm);
    else
        % no need to norm
        fourier = FC_full_orig.fourier;
    end
    % split data
    fourier1 = fourier(:,:,1:2:end,:);
    fourier2 = fourier(:,:,2:2:end,:);
    if loadfromdisk
        % save split files on disk first and use filenames as input
        splitfolder = [inputfolder 'split_data/'];
        if ~exist(splitfolder, 'dir')
           mkdir(splitfolder);
        end
        disp(['Saving split data in: ', splitfolder]);
        % save full data
        save([splitfolder 'fourfull.mat'], 'fourier', '-v7.3');
        % save halves
        fourier = fourier1;
        save([splitfolder 'fourprt1.mat'], 'fourier','-v7.3');
        fourier = fourier2;
        save([splitfolder 'fourprt2.mat'], 'fourier','-v7.3');
        
        % the following variables will be used for the fourier fields
        disp(['Loading data from: ' splitfolder]);
        fourfull = [splitfolder 'fourfull.mat'];
        fourprt1 = [splitfolder 'fourprt1.mat'];
        fourprt2 = [splitfolder 'fourprt2.mat'];       
    else
        % directly use the matrices
        disp('Using data already loaded in memory.');
        fourfull = fourier;
        fourprt1 = fourier1;
        fourprt2 = fourier2;
    end
else
    % dont need to split, so also dont need to norm. load directly from
    % file in the current directory, assuming the same filename convention
    disp(['Loading data from: ' inputfolder]);
    fourfull = [inputfolder 'fourfull.mat'];
    fourprt1 = [inputfolder 'fourprt1.mat'];
    fourprt2 = [inputfolder 'fourprt2.mat'];
end

% define output file
if ~exist(outputfolder, 'dir')
    mkdir(outputfolder);
end
outputpath = [outputfolder 'nwaydecomp.mat'];
    
% build fourierdata
fourierdata = [];
fourierdata.fourier        = fourfull;
fourierdata.fourierpart  = {fourprt1, fourprt2};
fourierdata.freq           = space_config.freqoi; %loaded from space_config
fourierdata.label          = space_config.labels; %loaded from space_config
fourierdata.dimord         = 'chan_freq_epoch_tap';
fourierdata.cfg             =    []; %empty config

% nway settings
cfg = [];
cfg.model              = 'spacetime';
cfg.datparam           = 'fourier';
cfg.Dmode              = 'identity';
cfg.numiter            = 3000;
cfg.convcrit           = 1e-8;
cfg.randstart          = 50;
cfg.outputfile         = outputpath;

% odd-even split settings
cfg.ncompestrandstart   = cfg.randstart;
cfg.ncompestsrdatparam  = 'fourierpart';
cfg.ncompest            = 'splitrel';

% grab input parameters for here, either set n-comp or search
cfg.ncompeststart       = ncomp_params(1);
cfg.ncompestend         = ncomp_params(2);
cfg.ncompeststep        = ncomp_params(3);

if cfg.ncompeststart == cfg.ncompestend
    % manually set component number , set critval to nans to fail
    cfg.ncompestsrcritval   = [nan nan 0 nan 0];
else
    % do search with these critvals
    cfg.ncompestsrcritval   = [.7 .7 0 .7 0];
end
cfg.ncompestsrcritjudge = 'meanoversplitscons';

if ~comp_local
    % compute on torque, set up core distribution and grab parameter for
    % core stack and time requirement
    cfg.distcomp.system          = 'torque';
    cfg.distcomp.timreq          = clus_params(1); % in seconds
    cfg.distcomp.matlabcmd       = '/opt/matlab/2015a/bin/matlab';
    cfg.distcomp.torquequeue     = 'hotel';
    cfg.distcomp.torquestack     = clus_params(2); % number of stacks per core
    cfg.distcomp.qsuboptions     = ' -k oe ';
end

disp(sprintf('Starting job with %i inits, with %i-init stacks each, allocating %i seconds per core.', cfg.randstart, clus_params(2), clus_params(1)));
disp(sprintf('Starting at %i components, to %i, at steps of %i',cfg.ncompeststart,cfg.ncompestend,cfg.ncompeststep));

% blast off!
nwaycomp = nd_nwaydecomposition(cfg,fourierdata);