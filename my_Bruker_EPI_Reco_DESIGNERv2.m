%% reconstruction of a Bruker EPI dataset
%% VARIANT with DESIGNER-v2 complex denoising (MP-PCA modules mpcomplex.py / mpdenoise_adaptive.py called from
%% Python, 2-pass scheme of designer -denoise -phase) on the complex images, before coil combination / magnitude,
%% and DESIGNER-v2 magnitude denoising of the final image + Rician bias correction (designer -denoise -rician);
%% the image without denoising still reproduces Paravision and is compared with it.
%% DESIGNER-v2: NYU, non-commercial research use only (license_DESIGNER-v2.txt in Designer_denoising_v2).
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
expno_dir = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\ON84_invivo3d\14\';
% expno_dir = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\Prostate15\12\';
% expno_dir = 'E:\Documents\shared_folder\20230423_101023_ON_79_mtbi_pilot_3d_1_1\12\';

do_pocs = 0;
% navigator correction (echo time shifting: segment s starts PVM_EpiEchoShiftA(s) later)
%   'fit'     = weighted linear fit of phase and log-amplitude of each full navigator (edges excluded)
%               vs absolute time (sample time + EpiEchoShiftA), all segments evaluated at the same
%               reference time -> works whatever the number of segments / overlap
%   'aligned' = navigators compared sample by sample on the time window shared by all segments
nav_mode = 'fit';

%% DESIGNER-v2 complex denoising (before coil combination and magnitude)
% DESIGNER-v2 MP-PCA modules (copied in Designer_denoising_v2) run with Anaconda, 2-pass complex scheme:
%   pass 1: complex MP-PCA (15x15x1 box, 'threshold', 'cordero-grande') -> denoised phase phi_dn
%   pass 2: MP-PCA of real(S*exp(-1i*phi_dn)) with the settings below -> signal*exp(1i*phi_dn)
do_designer = 1;
DESIGNER2PY.module_dir = 'C:\Users\User\Mon_Drive\Python\Designer_denoising_v2';
DESIGNER2PY.conda_root = 'C:\Users\User\anaconda3';       % Anaconda installation (fallback: C:\ProgramData\anaconda3)
DESIGNER2PY.conda_env  = 'base_python310';   % needs numpy, scipy, joblib, nibabel (+ scikit-image, numba, numpy-groupies for mpcomplex)
DESIGNER2PY.variant    = 'adaptive';         % 'mpcomplex' (designer -denoise) or 'adaptive' (designer -denoise -adaptive_patch)
DESIGNER2PY.shrinkage  = 'frob';              % 'threshold' (designer default) or 'frob'
DESIGNER2PY.algorithm  = 'jespersen';    % 'cordero-grande' (default), 'veraart' or 'jespersen'
DESIGNER2PY.n_cores    = -1;                  % joblib cores (-1 = all)
% MP-PCA kernel: set automatically from PVM_SpatDimEnum (3D -> [7 7 7], 2D -> [15 15 1], each slice on its own)
% topup: repetitions with opposite blips are denoised separately (TopUpYesNo / TopUpFirstOnly, as LaunchProcessing_3D)
% magnitude denoising of the final image (same DESIGNER2PY settings) + Rician bias correction sqrt(|S^2 - sigma^2|),
% sigma = MP-PCA noise map x sqrt(2) with double sampling (mean of the two DS images)
do_magnitude = 1;
keep_designer_files = 0;             % 1: keep the temporary NIfTI files

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

%% final image without denoising (Paravision-equivalent): sum of squares over receivers, mean of the DS images
img_abs = combine_rx_ds(imgCmp, rxDim, Param.DoubleSampling);
nDS = 1 + strcmp(Param.DoubleSampling,'yes');

%% denoising groups (topup) and MP-PCA kernel
% topup detection (as LaunchProcessing_3D): blip polarity of each repetition -> separate denoising groups
%   TopUpYesNo = yes & TopUpFirstOnly = yes : only repetition 2 is blip down (kept noisy, as the b0s of the pipeline)
%   TopUpYesNo = yes otherwise              : blip up / down alternate (odd / even repetitions), denoised separately
%   TopUpYesNo = no or absent               : all repetitions denoised together
blip = ones(1, size(imgCmp,4));                              % +1 blip up, -1 blip down, per repetition
if strcmp(Param.TopUpYesNo,'yes')
    if strcmp(Param.TopUpFirstOnly,'yes')
        blip(2) = -1;  topupMode = 'topup, reversed blip on repetition 2 only';
    else
        blip(2:2:end) = -1;  topupMode = 'topup, alternating blip up/down (odd/even repetitions)';
    end
else
    topupMode = 'no topup';
end
if do_designer || do_magnitude
    fprintf('Denoising groups: %s -> %d blip-up and %d blip-down repetitions\n', topupMode, sum(blip>0), sum(blip<0));
    switch Param.SpatDimEnum
        case '<3d>', DESIGNER2PY.kernel = [7 7 7];
        case '<2d>', DESIGNER2PY.kernel = [15 15 1];
        otherwise,   error('PVM_SpatDimEnum = %s not handled for the denoising kernel', Param.SpatDimEnum);
    end
    fprintf('MP-PCA kernel (%s): [%d %d %d]\n', Param.SpatDimEnum, DESIGNER2PY.kernel);
end

%% 1) DESIGNER-v2 complex denoising: each receiver and each double-sampling image is a 4D series [read, phase, slice|kz, rep]
if do_designer
    imgCmp_C = denoise_groups(@(x) designer2py_denoise(x, DESIGNER2PY, keep_designer_files), imgCmp, blip, 'DESIGNER v2 complex');
    img_C    = combine_rx_ds(imgCmp_C, rxDim, Param.DoubleSampling);
end

%% 2) DESIGNER-v2 magnitude denoising of the final image (designer -denoise) + Rician bias correction
% correction sqrt(|S^2 - sigma^2|) (designer -rician) with sigma = MP-PCA noise map x sqrt(2) with double sampling
% (the final image is the mean of the two DS images)
if do_magnitude
    [img_M, sig_M] = mag_denoise_groups(img_abs, blip, DESIGNER2PY, keep_designer_files);
    img_MR = rician_correct(img_M, sqrt(nDS)*sig_M);
end

%% display of the complex channel images (complex-denoised if do_designer)
if do_designer, imgDisp = imgCmp_C; else, imgDisp = imgCmp; end

%% display
for nrep_disp = 1:2:min(20,size(imgDisp,4))
    % nrep_disp =1;
    slice_disp =round(size(fidZF,3)/2);
    if strcmp(Param.DoubleSampling,'yes')
        if Param.EncNReceivers >1
            figure(2)
            subplot(2,4,1)
            imagesc(abs(imgDisp(:,:,slice_disp,nrep_disp,1,1)))
            subplot(2,4,4+1)
            imagesc(angle(imgDisp(:,:,slice_disp,nrep_disp,1,1)))
            subplot(2,4,2)
            imagesc(abs(imgDisp(:,:,slice_disp,nrep_disp,2,1)))
            subplot(2,4,4+2)
            imagesc(angle(imgDisp(:,:,slice_disp,nrep_disp,2,1)))
            subplot(2,4,3)
            imagesc(abs(imgDisp(:,:,slice_disp,nrep_disp,1,2)))
            subplot(2,4,4+3)
            imagesc(angle(imgDisp(:,:,slice_disp,nrep_disp,1,2)))
            subplot(2,4,4)
            imagesc(abs(imgDisp(:,:,slice_disp,nrep_disp,2,2)))
            subplot(2,4,4+4)
            imagesc(angle(imgDisp(:,:,slice_disp,nrep_disp,2,2)))
            linkaxes
        else
            figure(2)
            subplot(2,4,1)
            imagesc(abs(imgDisp(:,:,slice_disp,nrep_disp,1,1)))
            subplot(2,4,4+1)
            imagesc(angle(imgDisp(:,:,slice_disp,nrep_disp,1,1)))
            subplot(2,4,2)
            imagesc(abs(imgDisp(:,:,slice_disp,nrep_disp,2,1)))
            subplot(2,4,4+2)
            imagesc(angle(imgDisp(:,:,slice_disp,nrep_disp,2,1)))
            linkaxes
        end
    else
        if Param.EncNReceivers >1
            figure(2)
            subplot(2,2,1)
            imagesc(abs(imgDisp(:,:,slice_disp,nrep_disp,1)))
            subplot(2,2,2)
            imagesc(angle(imgDisp(:,:,slice_disp,nrep_disp,1)))
            subplot(2,2,3)
            imagesc(abs(imgDisp(:,:,slice_disp,nrep_disp,2)))
            subplot(2,2,4)
            imagesc(angle(imgDisp(:,:,slice_disp,nrep_disp,2)))
            linkaxes
        else
            figure(2)
            subplot(2,1,1)
            imagesc(abs(imgDisp(:,:,slice_disp,nrep_disp,1)))
            subplot(2,1,2)
            imagesc(angle(imgDisp(:,:,slice_disp,nrep_disp,1)))
            linkaxes
        end
    end

    drawnow
        pause(1)
end

%% Figures 3 to 12: repetitions 1:2:19 (blip up with alternating topup), middle slice
% row 1: Paravision, MATLAB without denoising, complex denoising, magnitude denoising + Rician correction (same scale)
% row 2: histogram of each image (same bins)
% row 3: Paravision - MATLAB (no denoising) and complex - magnitude denoising, with their histograms
slice_disp = round(size(imgCmp,3)/2);
bAll = rep_bvalues(expno_dir, size(imgCmp,4), blip);        % b-value of each repetition ([] if not available)
names = {'Paravision', 'MATLAB, no denoising'};
if do_designer,  names{end+1} = 'complex denoising'; end
if do_magnitude, names{end+1} = 'magnitude denoising + Rician corr.'; end
reps_disp = 1:2:min(19, size(imgCmp,4));
for k = 1:numel(reps_disp)
    r = reps_disp(k);
    if r <= size(imgBruker,4)
        imgPV = double(imgBruker(:,:,slice_disp,r));
        if Param.RecoTransposition(1)==1, imgPV = imgPV.'; end   % 2dseq stored [phase, read]
    else
        imgPV = nan(size(img_abs,1), size(img_abs,2));            % repetition not in the Paravision image
    end
    I = {imgPV, img_abs(:,:,slice_disp,r)};
    if do_designer,  I{end+1} = img_C(:,:,slice_disp,r);  end
    if do_magnitude, I{end+1} = img_MR(:,:,slice_disp,r); end
    Dn = {'Paravision - MATLAB (no denoising)'};  Dv = {I{1} - I{2}};
    if do_designer && do_magnitude
        Dn{end+1} = 'complex - magnitude denoising';  Dv{end+1} = I{3} - I{4};
    end
    nI = numel(I); nC = max(nI, 2*numel(Dv));
    figure(2 + k); clf; set(gcf, 'Position', [60 60 1500 850]); colormap gray
    tl = tiledlayout(3, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
    if isempty(bAll), bTxt = ''; else, bTxt = sprintf(', b = %.0f s/mm^2', bAll(r)); end
    title(tl, sprintf('repetition %d%s, slice %d', r, bTxt, slice_disp))
    cl = [0 prctile(I{2}(:), 99.5)];                               % common scale: MATLAB image without denoising
    edges = linspace(0, max(cellfun(@(x) prctile(x(:), 99.9), I)), 101);
    hAx = gobjects(1, nI);
    for m = 1:nI
        nexttile(tl, m); imagesc(I{m}, cl); axis image off; colorbar; title(names{m})
        hAx(m) = nexttile(tl, nC + m);
        histogram(I{m}(:), edges, 'EdgeColor', 'none'); xlabel('intensity'); ylabel('voxels'); title(['histogram: ' names{m}], 'FontWeight', 'normal')
    end
    linkaxes(hAx, 'xy');                                           % same axes for the 4 histograms
    for d = 1:numel(Dv)
        a = prctile(abs(Dv{d}(:)), 99.5); if ~(a > 0), a = 1; end   % symmetric scale for each difference
        nexttile(tl, 2*nC + 2*d - 1); imagesc(Dv{d}, [-a a]); axis image off; colorbar; title(Dn{d})
        nexttile(tl, 2*nC + 2*d); histogram(Dv{d}(:), linspace(-a, a, 101), 'EdgeColor', 'none')
        xlabel('difference'); ylabel('voxels'); title(['histogram: ' Dn{d}], 'FontWeight', 'normal')
    end
end


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
%topup (custom methods): reversed phase-encoding blips
Param.TopUpYesNo     = ReadPV360Param(expno_dir, 'TopUpYesNo') ;      % 'yes' / 'no' / [] (absent)
Param.TopUpFirstOnly = ReadPV360Param(expno_dir, 'TopUpFirstOnly') ;  % 'yes' / 'no' / [] (absent)


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


function out = denoise_groups(fun, img, grp, name)
% apply a denoising function separately to each group of repetitions (dim 4) sharing the same blip polarity;
% groups with fewer than 6 repetitions are left noisy (e.g. the single reversed-blip b0 of 'TopUpFirstOnly')
out = img;
for g = unique(grp(:)).'
    idx = find(grp == g);
    if numel(idx) < 6
        fprintf('%s: blip group %+d has %d repetition(s) -> not denoised\n', name, g, numel(idx));
        continue
    end
    fprintf('%s: denoising blip group %+d (%d repetitions)\n', name, g, numel(idx));
    out(:,:,:,idx,:,:) = fun(img(:,:,:,idx,:,:));
end
end


function imgOut = designer2py_denoise(img, opt, keepFiles)
% DESIGNER-v2 MP-PCA modules (mpcomplex.py or mpdenoise_adaptive.py) run with Anaconda, on each complex 4D series
% [read, phase, slice|kz, rep]. Same 2-pass scheme as DESIGNER-v2 denoise(..., phase=...), but calling the MP class
% directly on complex arrays (no antspyx, no phase file, no abs()):
%   pass 1: complex MP-PCA, 15x15x1 box, 'threshold', 'cordero-grande' -> denoised phase phi_dn
%   pass 2: MP-PCA of real(S*exp(-1i*phi_dn)) with the user settings -> signed real signal, put back as signal*exp(1i*phi_dn)
sz = size(img); sz(end+1:6) = 1; nv = sz(4);
assert(nv >= 6, 'DESIGNER v2 needs several repetitions (here %d)', nv);
work = [tempname '_designer2py']; mkdir(work);                % unique temporary folder
imgR = reshape(img, sz(1), sz(2), sz(3), nv, []); nser = size(imgR,5);
for v = 1:nser
    niftiwrite(single(real(imgR(:,:,:,:,v))), fullfile(work, sprintf('series%d_real.nii', v)));
    niftiwrite(single(imag(imgR(:,:,:,:,v))), fullfile(work, sprintf('series%d_imag.nii', v)));
end
py = {
'import os, sys'
'import numpy as np'
'import nibabel as ni'
'def main():'
'    sys.path.insert(0, sys.argv[1])'
'    work = sys.argv[2]; nser = int(sys.argv[3]); variant = sys.argv[4]'
'    kernel = [int(k) for k in sys.argv[5].split(",")]'
'    shrink = sys.argv[6]; algo = sys.argv[7]; ncores = int(sys.argv[8])'
'    if variant == "adaptive":'
'        import mpdenoise_adaptive as mpd'
'    else:'
'        import mpcomplex as mpc'
'    for v in range(1, nser + 1):'
'        nre = ni.load(os.path.join(work, "series%d_real.nii" % v))'
'        nim = ni.load(os.path.join(work, "series%d_imag.nii" % v))'
'        img = np.asarray(nre.dataobj, dtype=np.float64) + 1j * np.asarray(nim.dataobj, dtype=np.float64)'
'        nv = img.shape[-1]'
'        kp = [15, 15, 1]'
'        if variant == "adaptive":'
'            mp1 = mpd.MP(img, kp, patchtype="box", patchsize=None, shrinkage="threshold", algorithm="cordero-grande", crop=nv // 2, n_cores=ncores)'
'        else:'
'            mp1 = mpc.MP(img, kp, [kp[0] // 2 - 1, kp[1] // 2 - 1, 1], "threshold", "cordero-grande", min(nv // 2, int(np.prod(kp)) // 2), ncores)'
'        phi_dn = np.angle(mp1.process()[0])'
'        x = np.real(img * np.exp(-1j * phi_dn))'
'        if variant == "adaptive":'
'            mp2 = mpd.MP(x, kernel, patchtype="adaptive", patchsize=None, shrinkage=shrink, algorithm=algo, crop=0, n_cores=ncores)'
'        else:'
'            mp2 = mpc.MP(x, kernel, None, shrink, algo, 0, ncores)'
'        sig = np.real(mp2.process()[0])'
'        out = sig * np.exp(1j * phi_dn)'
'        ni.save(ni.Nifti1Image(out.real.astype(np.float32), nre.affine), os.path.join(work, "series%d_real_den.nii" % v))'
'        ni.save(ni.Nifti1Image(out.imag.astype(np.float32), nre.affine), os.path.join(work, "series%d_imag_den.nii" % v))'
'        print("DESIGNER v2 (%s): series %d / %d done" % (variant, v, nser), flush=True)'
'if __name__ == "__main__":'
'    main()'
};
fp = fopen(fullfile(work, 'run_designer2py.py'), 'w'); fprintf(fp, '%s\n', py{:}); fclose(fp);
conda_on_path(opt.conda_root);
cmd = sprintf('conda activate %s && cd /d "%s" && python "%s" "%s" "%s" %d %s %s %s %s %d', opt.conda_env, work, ...
    fullfile(work, 'run_designer2py.py'), opt.module_dir, work, nser, opt.variant, sprintf('%d,%d,%d', opt.kernel), ...
    opt.shrinkage, opt.algorithm, opt.n_cores);
fprintf('DESIGNER v2 (%s): %d complex series, %d repetitions\n', opt.variant, nser, nv);
status = system(cmd, '-echo');
if status ~= 0, error('DESIGNER v2 (python) failed.\nCommand: %s', cmd); end
for v = 1:nser
    Rr = double(niftiread(fullfile(work, sprintf('series%d_real_den.nii', v))));
    Ri = double(niftiread(fullfile(work, sprintf('series%d_imag_den.nii', v))));
    assert(isequal(size(Rr), [sz(1:3) nv]) && isequal(size(Ri), [sz(1:3) nv]), 'DESIGNER v2 output size mismatch (series %d)', v);
    imgR(:,:,:,:,v) = complex(Rr, Ri);
end
imgOut = reshape(imgR, size(img));
if ~keepFiles, rmdir(work, 's'); end
end


function conda_on_path(pyRoot)
% put Anaconda on the PATH of MATLAB (as in my_function_python_denoising_3D.m)
if ~isfolder(fullfile(pyRoot, 'Scripts')), pyRoot = 'C:\ProgramData\anaconda3'; end
p = strsplit(getenv('PATH'), ';');
add = {pyRoot; fullfile(pyRoot,'Library','mingw-w64','bin'); fullfile(pyRoot,'Library','usr','bin'); ...
       fullfile(pyRoot,'Library','bin'); fullfile(pyRoot,'Scripts'); fullfile(pyRoot,'bin')};
setenv('PATH', strjoin(unique([add; p(:)], 'stable'), ';'));
end


function img = combine_rx_ds(img, rxDim, DoubleSampling)
% Paravision combination: sum of squares over receivers, then mean of the two double-sampling images
img = sqrt(sum(abs(img).^2, rxDim));
if strcmp(DoubleSampling,'yes'), img = mean(img, 5); end
end


function b = rep_bvalues(expno_dir, nRep, blip)
% b-value of each repetition from BValueList (custom diffusion methods) or PVM_DwEffBval; [] if not available
b = [];
for name = {'BValueList', 'PVM_DwEffBval'}
    v = ReadPV360Param(expno_dir, name{1});
    if isnumeric(v) && ~isempty(v), b = double(v(:)).'; break; end
end
if isempty(b) || numel(b) == nRep, return; end
if 2*numel(b) == nRep && all(blip(2:2:end) == -1)
    b = repelem(b, 2);                    % alternating topup: each experiment acquired with both blips
elseif numel(b) + 1 == nRep && blip(2) == -1
    b = [b(1) b];                         % TopUpFirstOnly: repetition 2 is the reversed-blip b0
else
    warning('%d b-values for %d repetitions: not used', numel(b), nRep); b = [];
end
end


function out = rician_correct(img, sigma)
% DESIGNER-v2 Rician bias correction (run_rice_bias_correct): sqrt(|S^2 - sigma^2|), non-finite values set to 0
out = sqrt(abs(img.^2 - sigma.^2));
out(~isfinite(out)) = 0;
end


function [out, sig] = mag_denoise_groups(img, grp, opt, keepFiles)
% DESIGNER-v2 magnitude denoising of each blip group of the final image [x, y, z, rep];
% sig: MP-PCA noise map of each group, repeated over its repetitions. Groups with fewer than 6 repetitions stay
% noisy and get the noise map of the largest group (used only by the Rician correction).
out = img; sig = zeros(size(img)); sigMain = []; nMain = 0;
for g = unique(grp(:)).'
    idx = find(grp == g);
    if numel(idx) < 6
        fprintf('DESIGNER v2 magnitude: blip group %+d has %d repetition(s) -> not denoised\n', g, numel(idx));
        continue
    end
    fprintf('DESIGNER v2 magnitude: denoising blip group %+d (%d repetitions)\n', g, numel(idx));
    [out(:,:,:,idx), s] = designer2py_mag_denoise(img(:,:,:,idx), opt, keepFiles);
    sig(:,:,:,idx) = repmat(s, [1 1 1 numel(idx)]);
    if numel(idx) > nMain, sigMain = s; nMain = numel(idx); end
end
for g = unique(grp(:)).'
    idx = find(grp == g);
    if numel(idx) < 6 && ~isempty(sigMain), sig(:,:,:,idx) = repmat(sigMain, [1 1 1 numel(idx)]); end
end
end


function [imgOut, sigma] = designer2py_mag_denoise(img, opt, keepFiles)
% DESIGNER-v2 magnitude MP-PCA (mpcomplex / mpdenoise_adaptive denoise(...) without phase, as designer -denoise
% [-adaptive_patch]) of a magnitude series [read, phase, slice|kz, rep] with Anaconda;
% sigma: MP-PCA noise map [read, phase, slice|kz] (non-finite values set to 0, as designer)
sz = size(img); sz(end+1:4) = 1; nv = sz(4);
assert(isreal(img) && ndims(img) <= 4, 'magnitude series [x y z rep] expected');
assert(nv >= 6, 'DESIGNER v2 needs several repetitions (here %d)', nv);
work = [tempname '_designer2mag']; mkdir(work);              % unique temporary folder
niftiwrite(single(img), fullfile(work, 'mag.nii'));
py = {
'import os, sys'
'import numpy as np'
'import nibabel as ni'
'def main():'
'    sys.path.insert(0, sys.argv[1])'
'    work = sys.argv[2]; variant = sys.argv[3]'
'    kernel = [int(k) for k in sys.argv[4].split(",")]'
'    shrink = sys.argv[5]; algo = sys.argv[6]; ncores = int(sys.argv[7])'
'    nii = ni.load(os.path.join(work, "mag.nii"))'
'    img = np.asarray(nii.dataobj, dtype=np.float64)'
'    zero = img == 0'
'    img[zero] = np.finfo(img.dtype).eps'
'    if variant == "adaptive":'
'        import mpdenoise_adaptive as mpd'
'        mp = mpd.MP(img.copy(), kernel, patchtype="adaptive", patchsize=None, shrinkage=shrink, algorithm=algo, crop=0, n_cores=ncores)'
'        sig, sigma, npars = mp.process()'
'        sig = np.abs(sig)'
'    else:'
'        import mpcomplex as mpc'
'        mp = mpc.MP(img, kernel, None, shrink, algo, 0, ncores)'
'        sig, sigma, npars = mp.process()'
'        sig = np.real(sig)'
'    sig[zero] = 0'
'    sigma = np.real(sigma)'
'    sigma[~np.isfinite(sigma)] = 0'
'    ni.save(ni.Nifti1Image(sig.astype(np.float32), nii.affine), os.path.join(work, "mag_den.nii"))'
'    ni.save(ni.Nifti1Image(sigma.astype(np.float32), nii.affine), os.path.join(work, "sigma.nii"))'
'    print("DESIGNER v2 magnitude (%s) done" % variant, flush=True)'
'if __name__ == "__main__":'
'    main()'
};
fp = fopen(fullfile(work, 'run_designer2mag.py'), 'w'); fprintf(fp, '%s\n', py{:}); fclose(fp);
conda_on_path(opt.conda_root);
cmd = sprintf('conda activate %s && cd /d "%s" && python "%s" "%s" "%s" %s %s %s %s %d', opt.conda_env, work, ...
    fullfile(work, 'run_designer2mag.py'), opt.module_dir, work, opt.variant, sprintf('%d,%d,%d', opt.kernel), ...
    opt.shrinkage, opt.algorithm, opt.n_cores);
fprintf('DESIGNER v2 magnitude (%s): %d repetitions\n', opt.variant, nv);
status = system(cmd, '-echo');
if status ~= 0, error('DESIGNER v2 magnitude (python) failed.\nCommand: %s', cmd); end
imgOut = reshape(double(niftiread(fullfile(work, 'mag_den.nii'))), size(img));
sigma  = reshape(double(niftiread(fullfile(work, 'sigma.nii'))), sz(1:3));
if ~keepFiles, rmdir(work, 's'); end
end
