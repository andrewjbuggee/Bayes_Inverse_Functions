% create measurement prior


% By Andrew J. Buggee
%%


function [bayes_inputs] = create_measurement_covariance(bayes_inputs,data_struct,data_inputs)


covariance_type = bayes_inputs.measurement.covariance_type;


if strcmp(bayes_inputs.data_type,'aviris')==true
    
    bayes_inputs.spectral_bins = length(data_struct.wavelengths); % number of spectral channels, for this application
    
    
            
        % ---**--- Important Quantity ---**---
        bayes_inputs.measurement.uncertainty = 0.03; % percentage of measurement uncertainty according To Coddington et al. (2012) paragraph 24
        
          % measurement (spectral channel) is independent of one another
        bayes_inputs.measurement.covariance = diag(bayes_inputs.measurement.variance);
        
        % lets create our measurement pdf at each wavelength. We do this by
        % sampling a gaussian pdf with a  mean and variance as defined above
        
        
        
    end
    
    
elseif strcmp(bayes_inputs.data_type,'modis') == true
    
    % This defines the number of spectral channels used to create the
    % measurement covariance
    bayes_inputs.spectral_bins = size(data_struct,3);

    
        
        % ---**--- Important Quantity ---**---
        bayes_inputs.measurement.uncertainty = 0.02; % percentage of measurement uncertainty for reflectance according To "VALIDATION OF MODIS-DERIVED TOP-OF-ATMOSPHERE SPECTRAL RADIANCES BY MEANS OF VICARIOUS CALIBRATION"
        
         
        bayes_inputs.measurement.mean = radiance;
        bayes_inputs.measurement.variance = ones(length(wavelengths),1)*bayes_inputs.measurement.uncertainty; % Radiometric variance is some percentage of the mean
        
        % --------------------------------------------------------
        % Create the covariance matrix by taking the cross
        % correlation between spectral channels with the modis data
        % ---------------------------------------------------------
        
        if strcmp(covariance_type,'computed') == true
            data = cat(3,data_struct.EV.m250.reflectance,data_struct.EV.m500.reflectance);
            data = data(:,:,bayes_inputs.spectral_bins);
            for bb = 1:length(bayes_inputs.spectral_bins)
                for ii = 1:length(data_inputs.pixels2use.res500m.row)
                    data2run(ii,bb) = data(data_inputs.pixels2use.res500m.row(ii),data_inputs.pixels2use.res500m.col(ii),bb);
                end
            end
            
            bayes_inputs.measurement.covariance = cov(data2run);
            
        elseif strcmp(covariance_type,'independent') == true
            % create the covaraince matrix of the model parameters
            % if the covariance matrix is diagonal, then we are assuming each
            % measurement (spectral channel) is independent of one another
            bayes_inputs.measurement.covariance = diag(bayes_inputs.measurement.variance);
            
        else
            
            error('I dont understand how you want me to compute the covariance')
            
        end
        
        
        
        
    end
    
    



end

