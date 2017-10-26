function rg_spacedishlfp

comp_local = 0; %0 for cluster, 1 for local

if comp_local
    % local storage
    datafolder = '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/space_inputs/';
else
    % server storage
    datafolder = '/projects/ps-voyteklab/ricgao/data/spacedish/';
end

% load configs
space_config = load([datafolder 'space_configs.mat']);

% define filepath for fourier 
fourfull = [datafolder 'fourfull.mat'];
fourprt1 = [datafolder 'fourprt1.mat'];
fourprt2 = [datafolder 'fourprt2.mat'];

% define output file
outputpath = [datafolder 'nwaydecomp_set10comps.mat'];

% nwaydecomp settings
nwayalg      = 'spacetime'; % 'spacefsp'
nwaynmethod  = 'splitrel'; %'ncomp10sr'; %'ncomp15sr'; % 'splitrel'
nwaynrand    = 30;
nwayconvcrit = 1e-8;
normmethod   = 'none'; % coh none 16throotpower
nwaysplit    = 'oddeventrials'; %  oddeventrials

% build fourierdata
fourierdata = [];
fourierdata.fourier        = fourfull;
fourierdata.fourierpart  = {fourprt1, fourprt2};
fourierdata.freq           = space_config.freqoi; %loaded from space_config
fourierdata.label          = space_config.labels; %loaded from space_config
fourierdata.dimord         = 'chan_freq_epoch_tap';
%fourierdata.trialinfo      = data.trialinfo;
%fourierdata.cfg            = data.cfg;
fourierdata.cfg             =    []; %empty config

% nway settings
cfg = [];
cfg.model              = nwayalg;
cfg.datparam           = 'fourier';
cfg.Dmode              = 'identity';
cfg.numiter            = 3000;
cfg.convcrit           = nwayconvcrit;
cfg.randstart          = nwaynrand;
cfg.outputfile         = outputpath;

% odd-even split settings
cfg.ncompestrandstart   = nwaynrand;
cfg.ncompestsrdatparam  = 'fourierpart';
cfg.ncompest            = 'splitrel';
cfg.ncompeststart       = 10;
cfg.ncompeststep        = 1;
cfg.ncompestend         = 10;
cfg.ncompestsrcritval   = [nan nan 0 nan 0];
cfg.ncompestsrcritjudge = 'meanoversplitscons';

if ~comp_local
    % compute on torque, set up core distribution
    cfg.distcomp.system          = 'torque';
    cfg.distcomp.timreq          = 60*60*40; %
    cfg.distcomp.matlabcmd       = '/opt/matlab/2015a/bin/matlab';
    cfg.distcomp.torquequeue     = 'hotel';
    cfg.distcomp.torquestack     = 1;    
    cfg.distcomp.qsuboptions     = ' -k oe ';
end

% blast off!
nwaycomp = nd_nwaydecomposition(cfg,fourierdata);


