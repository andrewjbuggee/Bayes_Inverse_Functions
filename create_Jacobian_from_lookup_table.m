%% -- Create Jacobian Matrix from pre computed look up tables --


% By Andrew J. Buggee
%%


function [K] = create_Jacobian_from_lookup_table(inputs,data_inputs,folderName)

% Load the look-up table that is needed for this data set
if strcmp(inputs.data_type,'aviris')==true
    
    error('There are no look-up tables yet!')
    
    
elseif strcmp(inputs.data_type,'modis')==true
    
    num_pixels = data_inputs.inputs.pixels.num_2calculate;
    re = data_inputs.inputs.re; % - microns - effective radius - value is in the file name
    tau_c = data_inputs.inputs.tau_c; % - cloud optical depth
    interpGridScaleFactor = data_inputs.inputs.interpGridScaleFactor; % scale factor the will be used to increase the grid size for interpolation.
    
    R = data_inputs.R;
    
    num_bands = size(R,4); % number of bands the look-up table spans
    % lets interpolate the model data to increase our grid
    
    % in our model data, the column space, whixh spans the x direction, varies
    % with tau. The row space, which spans the y direction, vaires with re
    
    [X,Y] = meshgrid(tau_c,re);
    
    % now lets expand the two variables in question to the preferred size
    
    newTau_c = linspace(min(tau_c),max(tau_c),interpGridScaleFactor*length(tau_c));
    new_re = linspace(min(re),max(re),interpGridScaleFactor*length(re));
    
    % set up the new grid to interpolate on
    [Xq,Yq] = meshgrid(newTau_c,new_re);
    
    
    % the jacobian is calculate for each wavelength bin, and for each model
    % parameter
    K = cell(1,num_pixels);
    
    for pp = 1:num_pixels
        
        
        % lets shape our data in an easy to use format
        % first extract the data from the bands of interest
        
        % lets compare with non-interpolated data and use the 1km resolution
        % modis reflected data, right now, is in 500 meter resolution. So we
        % have to use the 500 meter pixel indexes
        
        
        R_pix = R(pp,:,:,:);
        R_pix = reshape(R_pix,size(R,2),size(R,3),size(R,4));
        
        % preform 2D interpolation
        newR = zeros(length(new_re),length(newTau_c),size(R_pix,3));
        
        for bb = 1:size(R_pix,3)
            newR(:,:,bb) = interp2(X,Y,R_pix(:,:,bb),Xq,Yq); % new interpolated look up table
        end
        
        % ----- Find Derivatives with Respect to Each Retrieval Variable -----
        % always taking the first derivatinve here, hence the 1 in diff
        derivative_newR_1 = diff(newR,1,1)./diff(Yq,1,1); % derivative with respect to row space (y -space) which is r
        derivative_newR_2 = diff(newR,1,2)./diff(Xq,1,2); % derivative with respect to column space (x-space) with is tau
        % we need these two derivatives to be the same size matrix. So we
        % drop the boundary values for both
        
        derivative_newR_1(:,end,:) = [];
        derivative_newR_2(end,:,:) = [];
        
        for ii = 1:size(derivative_newR_1,1)
            for jj = 1:size(derivative_newR_1,2)
                
                K{pp}(:,1,ii*jj) = reshape(derivative_newR_1(ii,jj,:),num_bands,1);
                K{pp}(:,2,ii*jj) = reshape(derivative_newR_2(ii,jj,:),num_bands,1);
                
            end
        end
        
        
        
    end
    
    
end











end