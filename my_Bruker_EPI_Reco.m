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
addpath(genpath([home_path filesep 'MATLABfunctions']));
addpath(genpath('functions'));

%% Data path
% expno_dir = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\EPI_test\9\';
expno_dir = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\EPI_test\18\';

do_pocs = 0;

%% retrieve Paravision version & import FID
Param.PV_version = ReadPV360Param(expno_dir, 'ACQ_sw_version');
if strcmp(Param.PV_version,'<pv 6.0.1>')
    [~,fid] = readPV360Raw([expno_dir 'fid']);
else
    [~,fid] = readPV360Raw([expno_dir 'rawdata.job0']);
end
fid = squeeze(fid);
Param.Npoints = numel(fid);

% figure()
% plot(abs(fid))

%% get PV parameters
Param = getPVparams(expno_dir,Param);

%% get reference image, used for comparison only
imageObj = ImageDataObject([expno_dir filesep, 'pdata', filesep, '1']);
imgBruker = squeeze(imageObj.data);

%% split navigator data
if strcmp(Param.EpiPrefixNavYes,'yes')
    if strcmp(Param.EpiSingleNav,'no')
        Param.Nnav_rep = Param.NRepetitions * Param.NEchoImages*Param.NSegments*Param.NSlice;
        fid = reshape(fid, numel(fid)/Param.Nnav_rep,Param.Nnav_rep);
        Nav = fid(1:Param.EpiPrefixNavSize,:);
        Nav = mean(abs(Nav(round(Param.EpiPrefixNavSize*0.1):round(Param.EpiPrefixNavSize*0.9),:))); % from 10 to 90%  arbitrary
        if strcmp(Param.SpatDimEnum,'<2d>')
            Param.total_rep = Param.NRepetitions * Param.NEchoImages;
            Nav =  reshape(Nav,[1,1, Param.EncNReceivers, Param.NSlice, Param.NSegments,  Param.total_rep]);
            Nav = permute(Nav,[1 2 5 4 6 3]);
        end
        if strcmp(Param.DoubleSampling,'yes')
            echotrainSize = Param.EncMatrix(1,1).* Param.EncMatrix(1,2)./ Param.NSegments * 2; %last *2 is for Double sampling
        else
            echotrainSize = Param.EncMatrix(1,1).* Param.EncMatrix(1,2)./ Param.NSegments;
        end
        fid = fid(end-echotrainSize+1:end,:);
        fid = fid(:);
    end
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

        % Navigator correction of 2 segments acquisition
        if strcmp(Param.EpiPrefixNavYes,'yes')
            if strcmp(Param.EpiSingleNav,'no')
                nav_norm = Nav(:,:,1,:,:,:)./Nav(:,:,2,:,:,:);
                fid(:,:,2,:,:,:) =fid(:,:,2,:,:,:).*repmat(nav_norm,size(fid,1),size(fid,2),1,1,1,2);
            end
        end
    else
        New_readout_size = size(fid,1)/((Param.EncMatrix(1,2)/Param.NSegments)*Param.EncNReceivers*Param.NSlice*Param.NSegments*Param.total_rep);
        Param.EncMatrix(1,1) = New_readout_size;% for ramp sampling compensation cases
        fid = reshape(fid,[Param.EncMatrix(1,1), Param.EncMatrix(1,2)/Param.NSegments, Param.EncNReceivers, Param.NSlice, Param.NSegments, Param.total_rep]);
        fid = permute(fid,[1 2 5 4 6 3]); %dim1 dim2 NSegments dim3 rep 1 receivers, might need special case for Nslice =1
    end
    % %% correction of EpiEchoShift, not working
    %     % RecoStabCorrFilter NAV0{innerSegments=1;navSize=0;echoSize=96;nEchoes=12;junkSize=0;nSlices=<NI>;nSegments=<PVM_EpiNShots>;sliding=true;useNav=true;dc=false;select=0;drift=false;dyncor=0;junkSize0=0;navSize2=0;multiecho=1};
    %     ind_odd = 1:2:size(fid,2);
    %     ind_even = 2:2:size(fid,2);
    %     x= linspace(-0.5, 0.5,size(fid,1));
    %     % How to use Param.EpiEchoShiftA Param.EpiEchoShiftB ?
    %         for ind_seg = 1:size(fid,3)
    %         % --- scaling  ---
    %         delta_k1 = Param.EpiEchoShiftB(ind_seg) * size(fid,1) / (4*pi);
    %         delta_k2 = Param.EpiEchoShiftA(ind_seg) * size(fid,1) / (4*pi);
    %
    %         % Phase ramps
    %         corrpos = exp(-1i * 2*pi * delta_k1 * x)';
    %         corrneg = exp(-1i * 2*pi * delta_k2 * x)';
    %
    %         fid(:,ind_odd,ind_seg,:,:,:) = fid(:,ind_odd,ind_seg,:,:,:).*repmat(corrpos,1,numel(ind_odd),1,size(fid,4),size(fid,5),size(fid,6));
    %         fid(:,ind_even,ind_seg,:,:,:) = fid(:,ind_even,ind_seg,:,:,:).*repmat(corrneg,1,numel(ind_even),1,size(fid,4),size(fid,5),size(fid,6));
    %         end
    %     %% End of correction of EpiEchoShift

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
    [fid] = my_EPI_reconX_trapezoid_compute_apply(fid,Param);

    if strcmp(Param.DoubleSampling,'no')
        %% ghost correction : might not match exactly the Bruker one
        % RecoEpiGhostFilter EG0{innerSegments=1;nechoes=64;phc=<PVM_EpiPhaseCorrection[0]>;nslices=<NI>;adjust=0;refecho=33;time=true;center=<PVM_EpiReadCenter>;channel=0;maxorder=1;doubleshotadj=0;dyncor=0;usedyn=1};
        disp('Ghost correction')
        [~,center_coordinate] = min(abs(Param.EpiTrajAdjkx-Param.EpiReadCenter));
        center_coordinate = center_coordinate/size(fid,1)-0.5;
        x= linspace(-0.5, 0.5,size(fid,1))+center_coordinate;

        for ind_slice = 1:size(fid,4)
            epiCorr_all =Param.EpiPhaseCorrection(1,:);
            epiCorr_all = reshape(epiCorr_all,2,size(fid,4));
            epiCorr = epiCorr_all(:,ind_slice);

            % --- scaling  ---
            delta_k1 = epiCorr(1) * size(fid,1) / (4*pi) / Param.NSegments; % / NSegments does not make sense
            delta_k2 = epiCorr(2) * size(fid,1) / (4*pi) / Param.NSegments;

            % Phase ramps
            corrpos = exp(-1i * 2*pi * delta_k1 * x)';
            corrneg = exp(-1i * 2*pi * delta_k2 * x)';

            fid(:,:,:,ind_slice,1,:,:) = FFTXSpace2KSpace(FFTKSpace2XSpace(fid(:,:,:,ind_slice,1,:,:),1).*repmat(corrpos,1,size(fid,2),size(fid,3),1,1,size(fid,6),size(fid,7)),1);
            fid(:,:,:,ind_slice,2,:,:) = FFTXSpace2KSpace(FFTKSpace2XSpace(fid(:,:,:,ind_slice,2,:,:),1).*repmat(corrneg,1,size(fid,2),size(fid,3),1,1,size(fid,6),size(fid,7)),1);
        end

        %% re-create correct dim 2 by combining odd and even echoes
        fid = permute(fid,[1 5 2 3 4 6 7]);
        fid = reshape(fid,size(fid,1),size(fid,2)*size(fid,3),size(fid,4),size(fid,5),size(fid,6),size(fid,7)); % combine odd even

        if add_line_to_even ==1
            fid = fid(:,1:end-1,:,:,:,:);
        end
    end

    %% combine segments
    fid = reshape(fid,size(fid,1),size(fid,2)*size(fid,3),size(fid,4),size(fid,5),size(fid,6));
    fid(:,Param.EncSteps1_pos,:,:,:) = fid;
    fid(:,:,Param.SliceOrder,:,:) = fid;

    %% Zero filling Partial FT reco
    fidZF = zeros(size(fid,1),Param.Matrix(1,2),size(fid,3),size(fid,4),size(fid,5));
    Nline = Param.Matrix(1,2)-size(fid,2);
    fidZF(:,Nline+1:end,:,:,:) = fid;

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
    % %% correction of EpiEchoShift , not working
    %     % RecoStabCorrFilter NAV0{innerSegments=1;navSize=0;echoSize=96;nEchoes=12;junkSize=0;nSlices=<NI>;nSegments=<PVM_EpiNShots>;sliding=true;useNav=true;dc=false;select=0;drift=false;dyncor=0;junkSize0=0;navSize2=0;multiecho=1};
    %     ind_odd = 1:2:size(fid,2);
    %     ind_even = 2:2:size(fid,2);
    %     x= linspace(-0.5, 0.5,size(fid,1));
    %     % How to use Param.EpiEchoShiftA Param.EpiEchoShiftB ?
    %         for ind_seg = 1:size(fid,3)
    %         % --- scaling  ---
    %         delta_k1 = Param.EpiEchoShiftA(ind_seg) * size(fid,1) / (4*pi);
    %         delta_k2 = Param.EpiEchoShiftB(ind_seg) * size(fid,1) / (4*pi);
    %
    %         % Phase ramps
    %         corrpos = exp(-1i * 2*pi * delta_k1 * x)';
    %         corrneg = exp(-1i * 2*pi * delta_k2 * x)';
    %
    %         fid(:,ind_odd,ind_seg,:,:,:) = fid(:,ind_odd,ind_seg,:,:,1).*repmat(corrpos,1,numel(ind_odd),1,size(fid,4),size(fid,5),size(fid,6));
    %         fid(:,ind_even,ind_seg,:,:,:) = fid(:,ind_even,ind_seg,:,:,2).*repmat(corrneg,1,numel(ind_even),1,size(fid,4),size(fid,5),size(fid,6));
    %         end
    %     %% End of correction of EpiEchoShift

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
    [fid] = my_EPI_reconX_trapezoid_compute_apply(fid,Param);

    if strcmp(Param.DoubleSampling,'no')
        %% ghost correction : might not match exactly the Bruker one
        disp('Ghost correction')
        [~,center_coordinate] = min(abs(Param.EpiTrajAdjkx-Param.EpiReadCenter));
        center_coordinate = center_coordinate/size(fid,1)-0.5;
        x= linspace(-0.5, 0.5,size(fid,1))+center_coordinate;
        coord_PhaseCorr = [1 2];
        epiCorr = Param.EpiPhaseCorrection(1,coord_PhaseCorr);

        % --- scaling  ---
        delta_k1 = epiCorr(1) * size(fid,1) / (4*pi) / Param.NSegments; % / NSegments does not make sense
        delta_k2 = epiCorr(2) * size(fid,1) / (4*pi) / Param.NSegments; % / NSegments does not make sense

        % Phase ramps
        corrpos = exp(-1i * 2*pi * delta_k1 * x)';
        corrneg = exp(-1i * 2*pi * delta_k2 * x)';

        fid(:,:,:,:,1,:,:) = FFTXSpace2KSpace(FFTKSpace2XSpace(fid(:,:,:,:,1,:,:),1).*repmat(corrpos,1,size(fid,2),size(fid,3),size(fid,4),1,size(fid,6),size(fid,7)),1);
        fid(:,:,:,:,2,:,:) = FFTXSpace2KSpace(FFTKSpace2XSpace(fid(:,:,:,:,2,:,:),1).*repmat(corrneg,1,size(fid,2),size(fid,3),size(fid,4),1,size(fid,6),size(fid,7)),1);
    end

    %% re-create correct dim 2 by combining odd and even echoes
    fid = permute(fid,[1 5 2 3 4 6 7]);
    fid = reshape(fid,size(fid,1),size(fid,2)*size(fid,3),size(fid,4),size(fid,5),size(fid,6),size(fid,7)); % combine odd even

    if add_line_to_even ==1
        fid = fid(:,1:end-1,:,:,:,:);
    end


    fid = reshape(fid,size(fid,1),size(fid,2)*size(fid,3),size(fid,4),size(fid,5),size(fid,6),size(fid,7)); % combine segments
    fid(:,Param.EncSteps1_pos,:,:,:,:) = fid;

    %% Zero filling Partial FT reco
    fidZF = zeros(size(fid,1),Param.Matrix(1,2),size(fid,3),size(fid,4),size(fid,5),size(fid,6));
    Nline = Param.Matrix(1,2)-size(fid,2);
    fidZF(:,Nline+1:end,:,:,:,:) = fid;
end
clearvars fid;

% reshape if Nslice == 1
if Param.NSlice ==1 && numel(size(fidZF))==5 && strcmp(Param.SpatDimEnum,'<2d>')==1
    fidZF = reshape(fidZF,size(fidZF,1),size(fidZF,2),1,size(fidZF,3),size(fidZF,4),size(fidZF,5));
end

% apply scaling on coils
for ind_coil = 1:Param.EncNReceivers
    if numel(size(fidZF))==5
        fidZF(:,:,:,:,ind_coil) = fidZF(:,:,:,:,ind_coil).*Param.RecoScaleChan(ind_coil);
    end
    if numel(size(fidZF))==6
        fidZF(:,:,:,:,:,ind_coil) = fidZF(:,:,:,:,:,ind_coil).*Param.RecoScaleChan(ind_coil);
    end
end

% apply phase offset between receivers
if Param.EncNReceivers > 1
    for ind_coil = 1:Param.EncNReceivers
        phase_mat{ind_coil} = -(1.*exp(1i*pi*(Param.RecoPhaseChan(ind_coil)/180)*(0:1/(size(fidZF,1)-1):1)))';

        if numel(size(fidZF))==5
            phase_mat{ind_coil} = repmat(phase_mat{ind_coil},1,size(fidZF,2),size(fidZF,3),size(fidZF,4),1);
            fidZF(:,:,:,:,ind_coil) = fidZF(:,:,:,:,ind_coil).*phase_mat{ind_coil};
        end
        if numel(size(fidZF))==6
            phase_mat{ind_coil} = repmat(phase_mat{ind_coil},1,size(fidZF,2),size(fidZF,3),size(fidZF,4),size(fidZF,5),1);
            fidZF(:,:,:,:,:,ind_coil) = fidZF(:,:,:,:,:,ind_coil).*phase_mat{ind_coil};
        end
    end
end

% apply phase shift for position correction
dims=[size(fidZF,1), size(fidZF,2), size(fidZF,3), size(fidZF,4), size(fidZF,5), size(fidZF,5)];
Offset_in_pixel = zeros(3,1);
Offset_in_pixel(2,1) =Param.Read_offset/(Param.FoV(1,1)).*((Param.Matrix(1,1)));
Offset_in_pixel(2,1) =-Param.Phase_1offset/(Param.FoV(1,2)).*((Param.Matrix(1,2)));
if strcmp(Param.SpatDimEnum,'<3d>')
    Offset_in_pixel(3,1) =-Param.Phase_1offset/(Param.FoV(1,3)).*((Param.Matrix(1,3)));
end


phase_matrix=ones(size(fidZF));
if strcmp(Param.DoubleSampling,'yes')
    if strcmp(Param.SpatDimEnum,'<2d>')
        for Offset_dim = 2:2 % shift in dim 1 is included in gridding
            f=0:1/(dims(Offset_dim)):1-(1/(dims(Offset_dim)));
            phase_vector=exp(1i*2*pi*Offset_in_pixel(Offset_dim)*f);
            switch Offset_dim
                case 1
                    phase_matrix=phase_matrix.*repmat(phase_vector, [1, dims(2), dims(3), dims(4), dims(5), dims(6)]);
                case 2
                    phase_matrix=phase_matrix.*repmat(phase_vector, [dims(1), 1, dims(3), dims(4), dims(5), dims(6)]);
            end
        end
    end

    if strcmp(Param.SpatDimEnum,'<3d>')
        for Offset_dim = 2:3 % shift in dim 1 is included in gridding
            f=0:1/(dims(Offset_dim)):1-(1/(dims(Offset_dim)));
            phase_vector=exp(1i*2*pi*Offset_in_pixel(Offset_dim)*f);
            switch Offset_dim
                case 1
                    phase_matrix=phase_matrix.*repmat(phase_vector, [1, dims(2), dims(3), dims(4), dims(5), dims(6)]);
                case 2
                    phase_matrix=phase_matrix.*repmat(phase_vector, [dims(1), 1, dims(3), dims(4), dims(5), dims(6)]);
                case 3
                    phase_matrix=phase_matrix.*repmat(permute(phase_vector,[1 3 2]), [dims(1), dims(2),1 , dims(4), dims(5), dims(6)]);
            end
        end
    end
    fidZF = fidZF.*phase_matrix;
else
    if strcmp(Param.SpatDimEnum,'<2d>')
        for Offset_dim = 2:2 % shift in dim 1 is included in gridding
            f=0:1/(dims(Offset_dim)):1-(1/(dims(Offset_dim)));
            phase_vector=exp(1i*2*pi*Offset_in_pixel(Offset_dim)*f);
            switch Offset_dim
                case 1
                    phase_matrix=phase_matrix.*repmat(phase_vector, [1, dims(2), dims(3), dims(4), dims(5)]);
                case 2
                    phase_matrix=phase_matrix.*repmat(phase_vector, [dims(1), 1, dims(3), dims(4), dims(5)]);
            end
        end
    end

    if strcmp(Param.SpatDimEnum,'<3d>')
        for Offset_dim = 2:3 % shift in dim 1 is included in gridding
            f=0:1/(dims(Offset_dim)):1-(1/(dims(Offset_dim)));
            phase_vector=exp(1i*2*pi*Offset_in_pixel(Offset_dim)*f);
            switch Offset_dim
                case 1
                    phase_matrix=phase_matrix.*repmat(phase_vector, [1, dims(2), dims(3), dims(4), dims(5)]);
                case 2
                    phase_matrix=phase_matrix.*repmat(phase_vector, [dims(1), 1, dims(3), dims(4), dims(5)]);
                case 3
                    phase_matrix=phase_matrix.*repmat(permute(phase_vector,[1 3 2]), [dims(1), dims(2),1 , dims(4), dims(5)]);
            end
        end
    end
    fidZF = fidZF.*phase_matrix;
end

% Fourier transform
if strcmp(Param.SpatDimEnum,'<2d>')
    imgCmp= FFTKSpace2XSpace(FFTKSpace2XSpace(fidZF,1),2);
else
    imgCmp= FFTKSpace2XSpace(FFTKSpace2XSpace(FFTKSpace2XSpace(fidZF,1),2),3);
end

%% display
nrep_disp =1;
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

if strcmp(Param.DoubleSampling,'yes')
    img_abs = sum(sqrt(abs(imgCmp.^2)),5); % sum over double sampling
else
    img_abs = abs(imgCmp); % abs
end

if Param.EncNReceivers > 1
    img_abs = sum(sqrt(abs(img_abs.^2)),6); % sum over coils
end

%% Display
slice_disp =round(size(fidZF,3)/2);
figure(3)
subplot(1,3,1)
imagesc(squeeze(abs(fidZF(:,:,slice_disp,1,1,1))))
title('Matlab k-space')
subplot(1,3,2)
imagesc(img_abs(:,:,slice_disp,1,1,1));
title('Matlab image')
subplot(1,3,3)
% imagesc(flip(flip(imgBruker(:,:,slice_disp,1),1),2)')
imagesc(flip(flip(imgBruker(:,:,slice_disp,1),1),2))
title('Paravision image')
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

% combination
Param.RecoScaleChan = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'RecoScaleChan') ;
Param.RecoPhaseChan = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'RecoPhaseChan') ;
Param.EpiPhaseCorrection = ReadPVParamRECO([expno_dir 'pdata' filesep '1' filesep], 'PVM_EpiPhaseCorrection') ;

% trajectories
Param.EpiTrajAdjkx = ReadPV360Param(expno_dir, 'PVM_EpiTrajAdjkx') ;
Param.EpiReadCenter = ReadPV360Param(expno_dir, 'PVM_EpiReadCenter') ;
Param.EpiReadAsym = ReadPV360Param(expno_dir, 'PVM_EpiReadAsym') ;
if isempty(Param.EpiReadCenter)==1
    Param.EpiReadCenter = round(size(Param.EpiTrajAdjkx,2)/2)-size(Param.EpiTrajAdjkx,2)*Param.EpiReadAsym;
end

Param.EpiPhaseCorrection = ReadPV360Param(expno_dir, 'PVM_EpiPhaseCorrection') ;
Param.EpiEchoShiftA = ReadPV360Param(expno_dir, 'PVM_EpiEchoShiftA') ;
Param.EpiEchoShiftB = ReadPV360Param(expno_dir, 'PVM_EpiEchoShiftB') ;


%Smotting of the trajectory
xTraj = (1:size(Param.EpiTrajAdjkx,2));
xqTraj = -4:size(Param.EpiTrajAdjkx,2)+5;
EpiTrajAdjkxPChip = pchip(xTraj,Param.EpiTrajAdjkx,xqTraj);
[~,Values] = spaps(xqTraj,EpiTrajAdjkxPChip,0.1);% 0.3
Param.EpiTrajAdjkxSodd = Values(1,6:end-5);
Param.EpiTrajAdjkxSeven = Param.EpiTrajAdjkxSodd(end)-(Param.EpiTrajAdjkxSodd);

% figure();
% plot(diff(Param.EpiTrajAdjkx),'k');hold on;
% plot(diff(Param.EpiTrajAdjkxSodd),'r');hold off;
% figure();
% plot([diff(Param.EpiTrajAdjkxSodd) diff(Param.EpiTrajAdjkxSeven)]);
end

function [fidR] = my_EPI_reconX_trapezoid_compute_apply(fid,Param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

maxK = round(max(Param.EpiTrajAdjkxSodd));

reconx.encodeNx_  = maxK; %Param.Matrix(1); % hdr.encoding.encodedSpace.matrixSize.x;
reconx.encodeFOV_ = Param.FoV(1); %hdr.encoding.encodedSpace.fieldOfView_mm.x;
reconx.reconNx_   = Param.Matrix(1); % hdr.encoding.reconSpace.matrixSize.x;
reconx.reconFOV_  = Param.FoV(1); %hdr.encoding.reconSpace.fieldOfView_mm.x;

trajectoryPos_ = (Param.EpiTrajAdjkxSodd'-max(abs(Param.EpiTrajAdjkxSodd/2)));
trajectoryNeg_ = (Param.EpiTrajAdjkxSeven'-(trajectoryPos_(end)));
reconx.numSamples_ = Param.EncMatrix(1);

% Compute the reconstruction operator
Km = floor(reconx.encodeNx_ / 2.0);
Ne = 2*Km + 1;

% resize the reconstruction operator
Mpos_=zeros(reconx.reconNx_,reconx.numSamples_);
Mneg_=zeros(reconx.reconNx_,reconx.numSamples_);

% evenly spaced k-space locations
keven = linspace(-Km, Km, Ne);
%keven.print("keven =");

% image domain locations [-0.5,...,0.5)
x = linspace(-0.5,(reconx.reconNx_-1.)/(2.*reconx.reconNx_),reconx.reconNx_);
x = x./Param.AntiAlias(1);
%x.print("x =");

% DFT operator
% Going from k space to image space, we use the IFFT sign convention
F=zeros(reconx.reconNx_, Ne);
fftscale = 1.0 / sqrt(Ne);

for p=1:1:reconx.reconNx_
    for q=1:1:Ne
        F(p,q) = fftscale * exp(complex(0.0,1.0*2*pi*keven(q)*x(p)));
    end
end
%F.print("F =");

% forward operators
Qp=zeros(reconx.numSamples_, Ne);
Qn=zeros(reconx.numSamples_, Ne);

for p=1:1:reconx.numSamples_
    %GDEBUG_STREAM(trajectoryPos_(p) << "    " << trajectoryNeg_(p) << std::endl);
    for q=1:1:Ne
        Qp(p,q) = sinc(trajectoryPos_(p)-keven(q));
        Qn(p,q) = sinc(trajectoryNeg_(p)-keven(q));
    end
end

%Qp.print("Qp =");
%Qn.print("Qn =");
% figure(1)
% subplot(2,2,1)
% imagesc(Qp)
% subplot(2,2,2)
% imagesc(Qn)
% subplot(2,2,3)
% plot(trajectoryNeg_)
% subplot(2,2,4)
% plot(trajectoryPos_)

% recon operators
Mp=zeros(reconx.reconNx_,reconx.numSamples_);
Mn=zeros(reconx.reconNx_,reconx.numSamples_);
Mp = F * pinv(Qp);
Mn = F * pinv(Qn);

% Compute the off-center correction:     /////

% Compute the off-center distance in the RO direction:
% roOffCenterDistance = calcOffCenterDistance( hdr_in );
roOffCenterDistance =Param.Read_offset; % is unit OK ?
%GDEBUG_STREAM("roOffCenterDistance_: " << roOffCenterDistance << ";       encodeFOV_: " << encodeFOV_);

my_keven = linspace(0, reconx.numSamples_ -1, reconx.numSamples_);
%         % find the offset:
%         % PV: maybe find not just exactly 0, but a very small number?
trajectoryPosArma = trajectoryPos_;

% n = find( trajectoryPosArma==0, 1, 'first');
[~,n] = min(abs(trajectoryPosArma));

my_keven = my_keven - (n-1);  %%attention 128

% my_keven
%         % Scale it:
%         % We have to find the maximum k-trajectory (absolute) increment:
Delta_k = abs( trajectoryPosArma(2:reconx.numSamples_,1) - trajectoryPosArma(1:reconx.numSamples_-1,1));
% max(Delta_k(:))

my_keven = my_keven *max(Delta_k(:));
%
%         % off-center corrections:


clear offCenterCorrN offCenterCorrP

myExponent = zeros(reconx.numSamples_,1);


for l=1:size(trajectoryPosArma,1)
    myExponent(l,1)=imag( 2*pi*roOffCenterDistance/reconx.encodeFOV_*(trajectoryPosArma(l,1)-my_keven(1,l)) );
end

offCenterCorrN = exp( myExponent );

for l=1:size(trajectoryPosArma,1)
    myExponent(l,1)=imag( 2*pi*roOffCenterDistance/reconx.encodeFOV_*(trajectoryNeg_(l,1)+my_keven(1,l)) );
end

offCenterCorrP = exp( myExponent );

%         Finally, combine the off-center correction with the recon operator:
Mp = Mp * diag(offCenterCorrP);
Mn = Mn * diag(offCenterCorrN);

%          and save it into the NDArray members:
for p=1:1:reconx.reconNx_
    for q=1:1:reconx.numSamples_
        Mpos_(p,q) = Mp(p,q);
        Mneg_(p,q) = Mn(p,q);
    end
end

fidR = zeros(Param.Matrix(1),size(fid,2),size(fid,3),size(fid,4),size(fid,5),size(fid,6),size(fid,7));
for dim2 =1:size(fid,2)
    for dim3 =1:size(fid,3)
        for dim4 =1:size(fid,4)
            for dim6 =1:size(fid,6)
                for dim7 =1:size(fid,7)
                    fidR(:,dim2,dim3,dim4,1,dim6,dim7) = Mpos_*squeeze(fid(:,dim2,dim3,dim4,1,dim6,dim7));
                    fidR(:,dim2,dim3,dim4,2,dim6,dim7) = Mneg_*squeeze(fid(:,dim2,dim3,dim4,2,dim6,dim7));
                end
            end
        end
    end
end

fidR = fftshift(fft(fftshift(fidR,1),[],1),1);
end