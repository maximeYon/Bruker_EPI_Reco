%% reconstruction of a Bruker EPI dataset
%% Coded by Maxime Yon in 2026
%% maxime.yon@univ-rennes.fr
%
% With help of Valéry Ozenne for the regridding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;

home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
% addpath(genpath('functions'));
addpath('functions');
addpath(genpath('pvtools'));

%% Data path
% expno_dir = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\EPI_test\9\';
% expno_dir = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\EPI_test\18\';
% expno_dir = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\ON84_invivo3d\14\';
expno_dir = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\Prostate15\12\';
% expno_dir = 'E:\Documents\shared_folder\20230423_101023_ON_79_mtbi_pilot_3d_1_1\12\';

do_pocs = 0;
% navigator correction (echo time shifting: segment s starts PVM_EpiEchoShiftA(s) later)
%   'fit'     = weighted linear fit of phase and log-amplitude of each full navigator (edges excluded)
%               vs absolute time (sample time + EpiEchoShiftA), all segments evaluated at the same
%               reference time -> works whatever the number of segments / overlap
%   'aligned' = navigators compared sample by sample on the time window shared by all segments
nav_mode = 'fit';

%% retrieve Paravision version & import FID
Param.PV_version = ReadPV360Param(expno_dir, 'ACQ_sw_version');
if isfile([expno_dir 'fid'])
    [~,fid] = readPV360Raw([expno_dir 'fid']);           % PV <= 7
else
    [~,fid] = readPV360Raw([expno_dir 'rawdata.job0']);  % PV360
end
fid = squeeze(fid);
Param.Npoints = numel(fid);

%% get PV parameters
Param = getPVparams(expno_dir,Param);

%% get reference image, used for comparison only
imageObj = ImageDataObject([expno_dir filesep, 'pdata', filesep, '1']);
imgBruker = squeeze(imageObj.data);

%% split navigator + navigator correction (RecoStabCorrFilter, drift=false)
% amplitude & zero-order phase of each segment aligned to segment 1,
% separately for each receiver / slice / repetition (or kz step in 3D)
if strcmp(Param.EpiPrefixNavYes,'yes') && strcmp(Param.EpiSingleNav,'no')
    scanSize = Param.ACQ_size(1)/2;   % complex points per scan: nav + junk + echo train
    assert(mod(numel(fid), scanSize*Param.EncNReceivers*Param.NSlice*Param.NSegments)==0, 'Unexpected scan size');
    % scan order: 2D [rx, slice, seg, rep...]  /  3D (NSlice=1) [rx, seg, kz, rep...]
    fid = reshape(fid, scanSize, Param.EncNReceivers, Param.NSlice, Param.NSegments, []);
    Nav = fid(1:Param.EpiPrefixNavSize,:,:,:,:);
    if strcmp(Param.DoubleSampling,'yes')
        echotrainSize = Param.EncMatrix(1,1).*Param.EncMatrix(1,2)./Param.NSegments*2; % *2 for double sampling
    else
        echotrainSize = Param.EncMatrix(1,1).*Param.EncMatrix(1,2)./Param.NSegments;
    end
    fid = fid(end-echotrainSize+1:end,:,:,:,:);           % drop nav + junk
    % timing of the navigator samples
    dwNav = Param.EpiPrefixNavDur/Param.EpiPrefixNavSize;                  % ms per nav sample
    A = zeros(1,Param.NSegments);                                          % start delay of each segment (ms)
    if strcmpi(Param.EpiEchoTimeShifting,'yes') && numel(Param.EpiEchoShiftA)==Param.NSegments
        A = Param.EpiEchoShiftA(:).';
    end
    navLo = 16; navHi = Param.EpiPrefixNavSize-6;                          % drop filter transients at both ends
    switch nav_mode
        case 'fit'
            % phase(t) and log|nav(t)| fitted as a + b*(t - Tref) for each receiver/slice/segment/rep,
            % t = absolute time after excitation. Weights |nav|^2 (variance of phase and log-amplitude
            % ~ 1/|nav|^2). Tref = centre of the time span covered by all navigators; the intercept a
            % is the value at Tref, so segments are compared at the same time.
            n    = (navLo:navHi).';
            Tref = ((navLo-1)*dwNav + min(A) + (navHi-1)*dwNav + max(A))/2;
            t    = (n-1)*dwNav + reshape(A,1,1,1,[]) - Tref;                % [Nn 1 1 NSeg]
            Y    = Nav(n,:,:,:,:);
            w    = abs(Y).^2;
            S0 = sum(w,1); S1 = sum(w.*t,1); S2 = sum(w.*t.^2,1); den = S0.*S2 - S1.^2;
            atTref = @(y) (S2.*sum(w.*y,1) - S1.*sum(w.*t.*y,1))./den;     % weighted LS intercept
            phT = atTref(unwrap(angle(Y),[],1));
            laT = atTref(log(abs(Y)+realmin));
            navPhase = angle(exp(1i*(phT(:,:,:,1,:) - phT)));            % segment 1 = reference
            navAmp   = exp(laT(:,:,:,1,:) - laT);
            clearvars Y w t S0 S1 S2 den phT laT
        case 'aligned'
            % nav sample n of segment s is at the same time as sample n+sh(s) of segment 1:
            % compare segments only on the time window they all share
            sh = round((A - A(1))/dwNav);
            u  = (navLo+max(sh)):(navHi+min(sh));                          % common window, segment-1 time base
            assert(numel(u)>=20, 'Navigator window common to all segments is too short (%d samples), use nav_mode = ''fit''', numel(u));
            NavAl = zeros([numel(u) size(Nav,2) size(Nav,3) size(Nav,4) size(Nav,5)], 'like', Nav);
            for s = 1:Param.NSegments
                NavAl(:,:,:,s,:) = Nav(u-sh(s),:,:,s,:);
            end
            Nav = NavAl; clearvars NavAl
            NavRef   = Nav(:,:,:,1,:);                    % segment 1 = reference
            navAmp   = sqrt(sum(abs(NavRef).^2,1)./sum(abs(Nav).^2,1));
            navPhase = angle(sum(NavRef.*conj(Nav),1));
        otherwise
            error('Unknown nav_mode ''%s''', nav_mode);
    end
    fid = fid.*(navAmp.*exp(1i*navPhase));                % segment 1 -> factor 1
    fid = fid(:);
    clearvars Nav NavRef
end

%% fid reshape
if isempty(Param.DwNDiffExp)==1 % allows also the reco of non classical diff experiment
    Param.DwNDiffExp =1;
end
Param.total_rep = Param.NRepetitions * Param.NEchoImages * Param.DwNDiffExp;

if strcmp(Param.SpatDimEnum,'<2d>')
    if strcmp(Param.DoubleSampling,'yes')
        New_readout_size = size(fid,1)/(2*(Param.EncMatrix(1,2)/Param.NSegments)*Param.EncNReceivers*Param.NSlice*Param.NSegments*Param.total_rep);
        Param.EncMatrix(1,1) = New_readout_size; % for ramp sampling compensation cases
        fid = reshape(fid,[Param.EncMatrix(1,1), 2, Param.EncMatrix(1,2)/Param.NSegments, Param.EncNReceivers, Param.NSlice, Param.NSegments,  Param.total_rep]);
        fid = permute(fid,[1 3 6 5 7 2 4]);

    else
        New_readout_size = size(fid,1)/((Param.EncMatrix(1,2)/Param.NSegments)*Param.EncNReceivers*Param.NSlice*Param.NSegments*Param.total_rep);
        Param.EncMatrix(1,1) = New_readout_size;% for ramp sampling compensation cases
        fid = reshape(fid,[Param.EncMatrix(1,1), Param.EncMatrix(1,2)/Param.NSegments, Param.EncNReceivers, Param.NSlice, Param.NSegments, Param.total_rep]);
        fid = permute(fid,[1 2 5 4 6 3]); %dim1 dim2 NSegments dim3 rep 1 receivers, might need special case for Nslice =1
    end

    if strcmp(Param.DoubleSampling,'no')
        %% use dim5 to store odd and even echoes
        fid = reshape(fid,size(fid,1),size(fid,2),size(fid,3),size(fid,4),1,size(fid,5),size(fid,6),size(fid,7));
        fid_odd = fid(:,1:2:end,:,:,:,:,:);
        fid_even = fid(:,2:2:end,:,:,:,:,:);
        add_line_to_even =0;
        if size(fid_odd,2)>size(fid_even,2)
            add_line_to_even =1;
            fid_even = cat(2,fid_even,zeros(size(fid_even,1),1,size(fid_even,3),size(fid_even,4),size(fid_even,5),size(fid_even,6),size(fid_even,7)));
        end
        fid = cat(5,fid_odd,fid_even);
    end

    %% regridding
    if strcmp(Param.DoubleSampling,'yes')
        polDim = 6; % [read, lines/seg, seg, slice, rep, DS, rx]
    else
        polDim = 5; % odd/even echoes stored in dim 5
    end
    [fid] = my_EPI_reconX_trapezoid_compute_apply(fid,Param,polDim);

    if strcmp(Param.DoubleSampling,'no')
        % (ghost correction RecoEpiGhostFilter is done inside the regridding)
        %% re-create correct dim 2 by combining odd and even echoes
        fid = permute(fid,[1 5 2 3 4 6 7]);
        fid = reshape(fid,size(fid,1),size(fid,2)*size(fid,3),size(fid,4),size(fid,5),size(fid,6),size(fid,7)); % combine odd even

        if add_line_to_even ==1
            fid = fid(:,1:end-1,:,:,:,:);
        end
    end

    %% combine segments (all other dimensions kept: slice, rep, [DS], rx)
    sz = size(fid); sz(end+1:7) = 1;
    fid = reshape(fid,[sz(1) sz(2)*sz(3) sz(4:end)]);      % [read, phase, slice, rep, (DS), rx]
    fid(:,Param.EncSteps1_pos,:) = fid(:,:,:);
    fid(:,:,Param.SliceOrder,:) = fid(:,:,:,:);

    %% Zero filling Partial FT reco
    szZF = size(fid); szZF(2) = Param.Matrix(1,2);
    fidZF = zeros(szZF);
    Nline = Param.Matrix(1,2)-size(fid,2);
    fidZF(:,Nline+1:end,:) = fid(:,:,:);

    %% POCS
    if do_pocs ==1
        for dim3 = 1:size(fidZF,3)
            for dim4 = 1:size(fidZF,4)
                for dim5 = 1:size(fidZF,5)
                    [~, fidZF(:,:,dim3,dim4,dim5)] = my_pocs(squeeze(fidZF(:,:,dim3,dim4,dim5)), 500, 0);
                end
            end
        end
    end

else
    if strcmp(Param.DoubleSampling,'yes')
        New_readout_size = size(fid,1)/(2*(Param.EncMatrix(1,2)/Param.NSegments)*Param.EncNReceivers*Param.NSlice*Param.NSegments*Param.EncMatrix(1,3)*Param.total_rep);
        Param.EncMatrix(1,1) = New_readout_size; % for ramp sampling compensation cases
        fid = reshape(fid,[Param.EncMatrix(1,1), 2,  Param.EncMatrix(1,2)/Param.NSegments, Param.EncNReceivers, Param.NSegments,  Param.EncMatrix(1,3), Param.total_rep]);
        fid = permute(fid,[1 3 5 6 7 2 4]); % dim1 dim2 NSegments dim3 rep double_sampling receivers
    else
        New_readout_size = size(fid,1)/((Param.EncMatrix(1,2)/Param.NSegments)*Param.EncNReceivers*Param.NSlice*Param.NSegments*Param.EncMatrix(1,3)*Param.total_rep);
        Param.EncMatrix(1,1) = New_readout_size; % for ramp sampling compensation cases
        fid = reshape(fid,[Param.EncMatrix(1,1), 1, Param.EncMatrix(1,2)/Param.NSegments, Param.EncNReceivers, Param.NSegments,  Param.EncMatrix(1,3), Param.total_rep ]);
        fid = permute(fid,[1 3 5 6 7 2 4]); %dim1 dim2 NSegments dim3 rep 1 receivers
    end

    if strcmp(Param.DoubleSampling,'no')
        %% use dim5 to store odd and even echoes
        if size(fid,5)~=1
            % fid = reshape(fid,size(fid,1),size(fid,2),size(fid,3),size(fid,4),1,size(fid,5));
            fid = permute(fid,[1 2 3 4 6 5 7]);
        end
        fid_odd = fid(:,1:2:end,:,:,:,:,:);
        fid_even = fid(:,2:2:end,:,:,:,:,:);
        add_line_to_even =0;
        if size(fid_odd,2)>size(fid_even,2)
            add_line_to_even =1;
            fid_even = cat(2,fid_even,zeros(size(fid_even,1),1,size(fid_even,3),size(fid_even,4),size(fid_even,5),size(fid_even,6),size(fid_even,7)));
        end
        fid = cat(5,fid_odd,fid_even);
    end

    %% regridding
    if strcmp(Param.DoubleSampling,'yes')
        polDim = 6; % [read, lines/seg, seg, dim3, rep, DS, rx]
    else
        polDim = 5; % odd/even echoes stored in dim 5
    end
    [fid] = my_EPI_reconX_trapezoid_compute_apply(fid,Param,polDim);

    % (ghost correction RecoEpiGhostFilter is done inside the regridding)

    if strcmp(Param.DoubleSampling,'no')
        %% re-create correct dim 2 by combining odd and even echoes
        fid = permute(fid,[1 5 2 3 4 6 7]);
        fid = reshape(fid,size(fid,1),size(fid,2)*size(fid,3),size(fid,4),size(fid,5),size(fid,6),size(fid,7)); % combine odd even

        if add_line_to_even ==1
            fid = fid(:,1:end-1,:,:,:,:);
        end
    end

    fid = reshape(fid,size(fid,1),size(fid,2)*size(fid,3),size(fid,4),size(fid,5),size(fid,6),size(fid,7)); % combine segments -> [read, phase, kz, rep, (DS), rx]
    fid(:,Param.EncSteps1_pos,:,:,:,:) = fid;

    %% Zero filling Partial FT reco (phase) + kz placement from PVM_EncSteps2 (acquisition order, partial FT in kz)
    kzPos = Param.EncSteps2(:).' + floor(Param.Matrix(1,3)/2) + 1;    % kz = 0 at floor(N/2) (0-based), as for ky
    assert(numel(kzPos)==size(fid,3) && all(kzPos>=1 & kzPos<=Param.Matrix(1,3)), ...
        'PVM_EncSteps2 (%d steps) does not match the kz dimension (%d) / matrix (%d)', numel(kzPos), size(fid,3), Param.Matrix(1,3));
    szZF = size(fid); szZF(end+1:6) = 1;
    szZF(2) = Param.Matrix(1,2); szZF(3) = Param.Matrix(1,3);
    fidZF = zeros(szZF);
    Nline = Param.Matrix(1,2)-size(fid,2);
    fidZF(:,Nline+1:end,kzPos,:) = fid(:,:,:,:);
end
clearvars fid;

%% receivers: Paravision RecoScaleFilter per channel (scale + constant phase offset)
% layout of fidZF: [read, phase, slice|kz, rep, (DS), rx]
if strcmp(Param.DoubleSampling,'yes'), rxDim = 6; else, rxDim = 5; end
chanFactor = Param.RecoScaleChan(:).*exp(1i*pi/180*Param.RecoPhaseChan(:));
assert(numel(chanFactor)==Param.EncNReceivers, 'RecoScaleChan/RecoPhaseChan do not match PVM_EncNReceivers');
assert(size(fidZF,rxDim)==Param.EncNReceivers, 'Receiver dimension %d has size %d, expected %d', rxDim, size(fidZF,rxDim), Param.EncNReceivers);
fidZF = fidZF.*reshape(chanFactor, [ones(1,rxDim-1) numel(chanFactor)]);

%% Fourier transform, reproducing Paravision orientation and position shifts
% read (dim 1): exact inverse of the k-space transform at the end of the regridding.
%   The position (RECO_offset, RECO_rotate) is already in the regridding grid.
imgCmp = ifftshift(ifft(ifftshift(fidZF,1),[],1),1);                 % Paravision scale (see fftscale in the regridding)
% phase (dim 2) [and slab (dim 3) in 3D]: RecoFTShiftFilter with RECO_rotate + RecoFTFilter (exponent +1)
%   RECO_rotate=0.5 is the centring, other values add the off-centre shift. The ramp uses a signed
%   k index (ky=0 at floor(N/2)) + ifftshift: same magnitude as Paravision, no (-1)^n phase checkerboard.
%   ifft includes the 1/N normalisation of Paravision's inverse FT (RecoFTFilter exponent +1).
ftDims = 2;
if strcmp(Param.SpatDimEnum,'<3d>'), ftDims = [2 3]; end
for d = ftDims
    N   = size(imgCmp,d);
    rot = Param.RecoRotate(d,:);                % one value per frame, acquisition order
    if strcmp(Param.SpatDimEnum,'<2d>')
        rot(Param.SliceOrder) = rot;            % same slice reordering as the data
        rot = reshape(rot,1,1,[]);              % slices are in dim 3
    end
    kc = reshape((0:N-1) - floor(N/2), [ones(1,d-1) N]);   % signed k index, ky=0 at floor(N/2)
    imgCmp = ifft(ifftshift(imgCmp.*exp(-1i*2*pi*rot.*kc), d), [], d);     % same values as Paravision, smooth phase
end

%% display
for nrep_disp = 1:2:min(20,size(imgCmp,4))
    % nrep_disp =1;
    slice_disp =round(size(fidZF,3)/2);
    if strcmp(Param.DoubleSampling,'yes')
        if Param.EncNReceivers >1
            figure(2)
            subplot(2,4,1)
            imagesc(abs(imgCmp(:,:,slice_disp,nrep_disp,1,1)))
            subplot(2,4,4+1)
            imagesc(angle(imgCmp(:,:,slice_disp,nrep_disp,1,1)))
            subplot(2,4,2)
            imagesc(abs(imgCmp(:,:,slice_disp,nrep_disp,2,1)))
            subplot(2,4,4+2)
            imagesc(angle(imgCmp(:,:,slice_disp,nrep_disp,2,1)))
            subplot(2,4,3)
            imagesc(abs(imgCmp(:,:,slice_disp,nrep_disp,1,2)))
            subplot(2,4,4+3)
            imagesc(angle(imgCmp(:,:,slice_disp,nrep_disp,1,2)))
            subplot(2,4,4)
            imagesc(abs(imgCmp(:,:,slice_disp,nrep_disp,2,2)))
            subplot(2,4,4+4)
            imagesc(angle(imgCmp(:,:,slice_disp,nrep_disp,2,2)))
            linkaxes
        else
            figure(2)
            subplot(2,4,1)
            imagesc(abs(imgCmp(:,:,slice_disp,nrep_disp,1,1)))
            subplot(2,4,4+1)
            imagesc(angle(imgCmp(:,:,slice_disp,nrep_disp,1,1)))
            subplot(2,4,2)
            imagesc(abs(imgCmp(:,:,slice_disp,nrep_disp,2,1)))
            subplot(2,4,4+2)
            imagesc(angle(imgCmp(:,:,slice_disp,nrep_disp,2,1)))
            linkaxes
        end
    else
        if Param.EncNReceivers >1
            figure(2)
            subplot(2,2,1)
            imagesc(abs(imgCmp(:,:,slice_disp,nrep_disp,1)))
            subplot(2,2,2)
            imagesc(angle(imgCmp(:,:,slice_disp,nrep_disp,1)))
            subplot(2,2,3)
            imagesc(abs(imgCmp(:,:,slice_disp,nrep_disp,2)))
            subplot(2,2,4)
            imagesc(angle(imgCmp(:,:,slice_disp,nrep_disp,2)))
            linkaxes
        else
            figure(2)
            subplot(2,1,1)
            imagesc(abs(imgCmp(:,:,slice_disp,nrep_disp,1)))
            subplot(2,1,2)
            imagesc(angle(imgCmp(:,:,slice_disp,nrep_disp,1)))
            linkaxes
        end
    end

    drawnow
        pause(1)
end

img_abs = sqrt(sum(abs(imgCmp).^2, rxDim));     % RecoSumOfSquaresFilter over receivers (= abs for 1 receiver)
if strcmp(Param.DoubleSampling,'yes')
    img_abs = mean(img_abs, 5);                 % Paravision: (|pos| + |neg|)/2
end

%% Display
slice_disp =round(size(fidZF,3)/2);
figure(3)
subplot(1,4,1)
imagesc(squeeze(abs(fidZF(:,:,slice_disp,1,1,1))))
title('Matlab k-space')
subplot(1,4,2)
imagesc(img_abs(:,:,slice_disp,1,1,1));
title('Matlab image')
imgPV = imgBruker(:,:,slice_disp,1);
if Param.RecoTransposition(1)==1, imgPV = imgPV.'; end   % 2dseq stored [phase, read]
subplot(1,4,3)
imagesc(imgPV)
% imagesc(flip(flip(imgBruker(:,:,slice_disp,1),1),2))
title('Paravision image')
subplot(1,4,4)
imagesc(imgPV-img_abs(:,:,slice_disp,1))
title('difference')
colormap gray


function [Param] = getPVparams(expno_dir,Param)
%% Get parameter from procno
%Encoding
Param.SpatDimEnum=ReadPV360Param(expno_dir, 'PVM_SpatDimEnum') ;
Param.EncMatrix = ReadPV360Param(expno_dir, 'PVM_EncMatrix') ;
Param.Matrix = ReadPV360Param(expno_dir, 'PVM_Matrix') ;
Param.NSlice = ReadPV360Param(expno_dir, 'PVM_SPackArrNSlices') ;
Param.NRepetitions = ReadPV360Param(expno_dir, 'PVM_NRepetitions') ;
Param.NAverages = ReadPV360Param(expno_dir, 'PVM_NAverages') ;
Param.NEchoImages = ReadPV360Param(expno_dir, 'PVM_NEchoImages') ;
Param.NSegments = ReadPV360Param(expno_dir, 'NSegments') ;
Param.ObjOrderList = ReadPV360Param(expno_dir, 'PVM_ObjOrderList') ;
Param.DwNDiffExp = ReadPV360Param(expno_dir, 'PVM_DwNDiffExp') ;
Param.ACQ_size = ReadPV360Param(expno_dir, 'ACQ_size') ; % [2*complex points per scan, ...]
%Receivers
Param.EncNReceivers = ReadPV360Param(expno_dir, 'PVM_EncNReceivers') ;
Param.EncChanScaling = ReadPV360Param(expno_dir, 'PVM_EncChanScaling') ;
Param.EncChanPhase = ReadPV360Param(expno_dir, 'PVM_EncChanPhase') ;
%Offset and acceleration
Param.EncPftOverscans = ReadPV360Param(expno_dir, 'PVM_EncPftOverscans') ;
Param.Param.MissingLines = Param.EncMatrix-2*Param.EncPftOverscans;
Param.Read_offset = ReadPV360Param([expno_dir filesep], 'PVM_SPackArrReadOffset') ;
Param.Phase_1offset = ReadPV360Param([expno_dir filesep], 'PVM_SPackArrPhase1Offset') ;
Param.Phase_2offset = ReadPV360Param([expno_dir filesep], 'PVM_SPackArrPhase2Offset') ;
Param.EffSliceOffset = ReadPV360Param([expno_dir filesep], 'PVM_EffSliceOffset') ;
Param.FoV = ReadPV360Param([expno_dir filesep], 'PVM_Fov') ;
Param.EncSteps1 = ReadPV360Param([expno_dir filesep], 'PVM_EncSteps1') ;
Param.EncSteps1_pos = Param.EncSteps1+abs(min(Param.EncSteps1))+1;
Param.EncSteps2 = ReadPV360Param([expno_dir filesep], 'PVM_EncSteps2') ; % kz encoding steps (3D), acquisition order
Param.PVM_SpatResol = ReadPV360Param([expno_dir filesep], 'PVM_SpatResol') ;
Param.PVM_SliceThick = ReadPV360Param([expno_dir filesep], 'PVM_SliceThick') ;
%Options
Param.DoubleSampling = ReadPV360Param(expno_dir, 'PVM_EpiCombine') ;
Param.AntiAlias = ReadPV360Param(expno_dir, 'PVM_AntiAlias') ;
Param.SliceOrder = ReadPV360Param(expno_dir, 'PVM_ObjOrderList')+1 ;

%navigator
Param.EpiPrefixNavYes = ReadPV360Param(expno_dir, 'PVM_EpiPrefixNavYes') ;
Param.EpiSingleNav = ReadPV360Param(expno_dir, 'PVM_EpiSingleNav') ;
Param.EpiPrefixNavSize = ReadPV360Param(expno_dir, 'PVM_EpiPrefixNavSize') ;
Param.EpiPrefixNavDur = ReadPV360Param(expno_dir, 'PVM_EpiPrefixNavDur') ;     % ms
Param.EpiEchoTimeShifting = ReadPV360Param(expno_dir, 'PVM_EpiEchoTimeShifting') ;

% combination
Param.RecoScaleChan = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'RecoScaleChan') ;
Param.RecoPhaseChan = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'RecoPhaseChan') ;
Param.EpiPhaseCorrection = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'PVM_EpiPhaseCorrection') ;
% Paravision FT shift / cutoff (rows = FT dims, columns = frames in acquisition order)
Param.RecoRotate = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'RECO_rotate') ;
Param.RecoOffset = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'RECO_offset') ;
Param.RecoSize   = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'RECO_size') ;
Param.RecoFtSize = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'RECO_ft_size') ;
Param.RecoTransposition = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'RECO_transposition') ;

% trajectories
Param.EpiTrajAdjkx = ReadPV360Param(expno_dir, 'PVM_EpiTrajAdjkx') ;
Param.EpiReadCenter = ReadPV360Param(expno_dir, 'PVM_EpiReadCenter') ;
Param.EpiReadCenterPV = Param.EpiReadCenter; % raw value (index of k=0 in the regridded readout), empty in PV6
Param.EpiTrajAdjb0 = ReadPV360Param(expno_dir, 'PVM_EpiTrajAdjb0') ; % measured trajectory phase (regridding)
Param.EpiReadAsym = ReadPV360Param(expno_dir, 'PVM_EpiReadAsym') ;
if isempty(Param.EpiReadCenter)==1
    Param.EpiReadCenter = round(size(Param.EpiTrajAdjkx,2)/2)-size(Param.EpiTrajAdjkx,2)*Param.EpiReadAsym;
end

Param.EpiPhaseCorrection = ReadPV360Param(expno_dir, 'PVM_EpiPhaseCorrection') ;
Param.EpiEchoShiftA = ReadPV360Param(expno_dir, 'PVM_EpiEchoShiftA') ; % echo time shifting: delay before nav + echo train, per segment (ms)
Param.EpiEchoShiftB = ReadPV360Param(expno_dir, 'PVM_EpiEchoShiftB') ;


end

function [fidR] = my_EPI_reconX_trapezoid_compute_apply(fid,Param,polDim)
% Regrids the EPI readout (dim 1) with the positive / negative trajectory operators.
% polDim : dimension holding the readout polarity (index 1 = positive, 2 = negative).
%          All other dimensions (segments, slices, reps, receivers...) are arbitrary.

% Paravision RecoRegridFilter (rampForm=5), identified exactly on a Paravision fidRegrid file:
%   out(j) = sum_t s(t) * sinc(k_t - j) * cos^2(pi*(k_t - j)/10) * [|k_t - j| < 5] * dk/dt(t)
%   k_t = PVM_EpiTrajAdjkx (ktraject), j = 0..Nk-1, Nk = floor(max(k_t))+1, dk/dt = gradient(k_t)
%   positive echoes * exp(-i*b0), negative echoes * exp(+i*(b0 - b0(end))) (b0traject = PVM_EpiTrajAdjb0),
%   negative echoes regridded with the same forward trajectory, then reversed (RecoReverseFilter).
traj = Param.EpiTrajAdjkx(:);                        % raw measured trajectory (as Paravision)
numSamples = numel(traj);
assert(size(fid,1)==numSamples, 'Readout has %d points but PVM_EpiTrajAdjkx has %d', size(fid,1), numSamples);
Nk   = floor(max(traj)) + 1;                         % regridded size
cen  = Param.EpiReadCenterPV;                        % index of k = 0 (PVM_EpiReadCenter)
if isempty(cen), cen = floor((Nk-1)/2); end          % e.g. PV6 (no PVM_EpiReadCenter)
keven = (0:Nk-1) - cen;                              % evenly spaced k-space locations
Ne    = numel(keven);

D  = traj.' - (0:Nk-1).';                            % [Nk x numSamples]: k_t - j
Kc = sinc(D) .* (abs(D)<5) .* cos(pi*D/10).^2 .* gradient(traj).';

b0 = Param.EpiTrajAdjb0(:);
if isempty(b0), b0 = zeros(numSamples,1); end
pQp = Kc .* exp(-1i*b0).';                           % positive echoes
pQn = Kc(end:-1:1,:) .* exp(1i*(b0 - b0(end))).';    % negative echoes: same regrid, then REV

% Image-domain positions = Paravision grid, in units of the encoded FOV:
% RecoZfill to RECO_ft_size, RecoFTShift with RECO_rotate, RecoFT (exponent +1),
% RecoCutoff with RECO_offset / RECO_size  ->  x = ((0:RECO_size-1)+RECO_offset)/RECO_ft_size - RECO_rotate
% This includes the centring (RECO_rotate=0.5), the read off-centre and the anti-alias crop.
% One grid per frame (2D: slices in raw acquisition order, dim 4 of fid).
nFrames = size(Param.RecoRotate,2);
Nout    = Param.RecoSize(1);
fftscale = 1.0 / Param.RecoFtSize(1);               % Paravision: zero filling to RECO_ft_size + inverse FT normalised by 1/N

sz = size(fid); sz(end+1:max([polDim,2,4,7])) = 1;   % receivers are in dim 7
assert(sz(polDim)==2, 'Polarity dimension %d must have size 2 (pos/neg)', polDim);
assert(nFrames==1 || sz(4)==nFrames, 'RECO_rotate has %d frames but fid dim 4 has %d', nFrames, sz(4));

% Ghost correction (Paravision RecoEpiGhostFilter, not used with double sampling):
% PVM_EpiPhaseCorrection(receiver, [2*frame-1, 2*frame]) = [ph0 ph1] (frames in acquisition order).
% Negative echoes get the relative phase ph0 - ph1*n, n = x*RECO_ft_size (read pixel, 0 at the centre).
ghost = strcmp(Param.DoubleSampling,'no') && ~isempty(Param.EpiPhaseCorrection);
if ghost
    assert(size(Param.EpiPhaseCorrection,2) >= 2*nFrames, 'PVM_EpiPhaseCorrection has %d values per receiver, expected %d', size(Param.EpiPhaseCorrection,2), 2*nFrames);
end

fidR = zeros([Nout sz(2:end)]);
idx  = repmat({':'}, 1, numel(sz));
for f = 1:nFrames
    x = ((0:Nout-1).' + Param.RecoOffset(1,f))/Param.RecoFtSize(1) - Param.RecoRotate(1,f);
    F = fftscale*exp(1i*2*pi*x*keven);      % DFT operator, inverse-FT sign convention (as Paravision)
    if nFrames > 1, idx{4} = f; end
    for r = 1:sz(7)
        if ghost
            row = min(r, size(Param.EpiPhaseCorrection,1));
            ph  = Param.EpiPhaseCorrection(row, 2*f-1:2*f);          % [ph0 ph1]
            gN  = exp(1i*(ph(1) - ph(2)*x*Param.RecoFtSize(1)));     % ghost correction, negative echoes
        else
            gN  = 1;
        end
        M = {F*pQp, (gN.*F)*pQn};           % recon operators, positive / negative readout
        idx{7} = r;
        for pol = 1:2
            idx{polDim} = pol;
            A   = fid(idx{:});
            szA = size(A); szA(end+1:numel(sz)) = 1;
            fidR(idx{:}) = reshape(M{pol} * reshape(A, sz(1), []), [Nout szA(2:end)]);
        end
    end
end

fidR = fftshift(fft(fftshift(fidR,1),[],1),1);   % back to k-space (inverted exactly before display)
end