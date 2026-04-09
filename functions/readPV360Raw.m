function [Acqp,data] = readPV360Raw(path_to_jobFile)
%% Open PV360 Job

% Retrieve ACQP parameters
directory_path = split(path_to_jobFile,filesep);
directory_path = join(directory_path(1:end-1,1),filesep,1); directory_path = directory_path{1};
Acqp = readBrukerParamFile([directory_path filesep 'acqp']);
if Acqp.ACQ_jobs_size ==0
   Acqp.ACQ_jobs_size = 1; 
end
endian = Acqp.BYTORDA(1,1);
numSelectedReceivers = bruker_getSelectedReceivers(Acqp);
single_bool = 1; % double;
isComplexRaw = 1; % complex data 

switch(Acqp.ACQ_word_size)
    case ('_32_BIT')
        format='int32';
        bits=32;
    case ('_16_BIT')
        format='int16';
        bits=16;
    case ('_32_BIT_FLOAT')
        format='float32';
        bits=32;
    otherwise
        format='int32';
        disp('Data-Format not correct specified! Set to int32')
        bits=32;
end


try
    fileID = fopen(path_to_jobFile,'r');
catch
    fileID = -1;
end
if fileID == -1
    error('Cannot open parameter file. Problem opening file %s.',path_to_jobFile);
end

%read File to Variable X with fread() and make it single-precission with single():
%  Attention: Matlab ACQ_size(1) = ACQ_size(0)
if(single_bool)
    X=( fread(fileID, [numSelectedReceivers*Acqp.ACQ_jobs_size, inf], [format, '=>single'], 0, endian) );
else
    X=fread(fileID, [numSelectedReceivers*Acqp.ACQ_jobs_size, inf], format, 0, endian);
end

dim1=numel(X) / (Acqp.ACQ_jobs_size*numSelectedReceivers);
X=reshape(X,[numSelectedReceivers, Acqp.ACQ_jobs_size,dim1]);

% Save to output variable:
if isComplexRaw
    % convert to complex:
    X=complex(X(:,:,1:2:end,:), X(:,:,2:2:end,:));
    % else: don't convert, only save
end

% file close
fclose(fileID);

data = X;
end