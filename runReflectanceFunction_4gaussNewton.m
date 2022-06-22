%% ----- Run the Reflectance Function calculation over a range of re and tau -----


% the input names file must chagne tau across the column space, and must
% change r across row space

% By Andrew J. Buggee

%%

function [R,Rl] = runReflectanceFunction_4gaussNewton(names,INP_folderName,saveCalculations_fileName)

% what computer are we using?

userName = whatComputer;

if strcmp(userName,'anbu8374')
    
    libRadTran_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4'];
    
elseif strcmp(userName,'andrewbuggee')
    
    libRadTran_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4'];
    
else
    error('I dont recognize this computer user name')
end


addpath(libRadTran_path);

% ----- extract inputs -----

inp_folder = [libRadTran_path,'/',INP_folderName]; % where the newly created .inp files will be saved
inputFileNames = names.inp;
outputFileNames = names.out;



R = zeros(size(inputFileNames)); % each value here is integrated over the band provided
Rl = cell(size(inputFileNames)); % each value here is the spectral reflectance over the entire band




% step through the band dimension
parfor bb = 1:length(inputFileNames)
    
    % --- For now, calculate inputSettings every time ---
    
    % start by running uvspec
    [inputSettings] = runUVSPEC(inp_folder,inputFileNames{bb},outputFileNames{bb});
    
    
    
    % read in the data structure and calculate reflectance function
    
    
    
    
    [ds,~,~] = readUVSPEC(inp_folder,outputFileNames{bb},inputSettings(2,:)); % headers don't change per iteration
    [R(bb),~] = reflectanceFunction(inputSettings(2,:),ds);
    
end




% save the relfectance calculation
% save(saveCalculations_fileName,"R",'-append'); % save inputSettings to the same folder as the input and output file






end

