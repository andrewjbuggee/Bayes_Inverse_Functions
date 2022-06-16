function [bayes_inputs] = create_model_prior_and_covariance(bayes_inputs,truthTable)

% -------------------------------------------------------------
% -------------------------------------------------------------
% THIS CODE ASSUMES THE RETIREVAL VARIABLES ARE ALL INDEPENDENT
% -------------------------------------------------------------
% -------------------------------------------------------------


% defien the model variance and mean using the Truth Table found by my
% TBLUT algorithm

if bayes_inputs.numPixels2Calculate<=size(truthTable,1)
    
    % Pick the first n pixels in the table
    n = bayes_inputs.numPixels2Calculate;
    
    % the order of the values below: (r_top, r_bottom, tau_c)
    % The model mean is the a priori; our first guess
    bayes_inputs.model.mean = [truthTable.estR17(1:n), truthTable.estR17(1:n), truthTable.estT17(1:n)]; % expected values for the effective radius (microns) and the optical depth
    
    
    % lets create the variance and mean for each model parameter
    % Using the same values defined by King and Vaughn (2012)
    bayes_inputs.model.variance = [linspace(3,3,n)',linspace(10,10,n)',0.1*truthTable.estT17(1:n)]; % variance for the effective radius (microns squared) and optical thickness respectively
    
    
    % Define the initial guess based off of the variance
    bayes_inputs.model.initialGuess = bayes_inputs.model.mean - sqrt(bayes_inputs.model.variance);
    
    % For now lets claim the desired variables are independent
    for ii = 1:n
        bayes_inputs.model.covariance(:,:,ii) = diag(bayes_inputs.model.variance(ii,:));
    end
    
else
    
    % We will only use what is in the truth table!
    bayes_inputs.numPixels2Calculate = size(truthTable,1);
    n = size(truthTable,1);
    % the order of the values below: (r_top, r_bottom, tau_c)
    % The model mean is the a priori; our first guess
    bayes_inputs.model.mean = [truthTable.estR17, truthTable.estR17, truthTable.estT17]; % expected values for the effective radius (microns) and the optical depth
    
    
    % lets create the variance and mean for each model parameter
    % Using the same values defined by King and Vaughn (2012)
    bayes_inputs.model.variance = [linspace(3,3,n)',linspace(10,10,n)',0.1*truthTable.estT17]; % variance for the effective radius (microns squared) and optical thickness respectively
    
    
    % For now lets claim the desired variables are independent
    for ii = 1:n
        
        bayes_inputs.model.covariance(:,:,ii) = diag(bayes_inputs.model.variance(ii,:));
    end
    
    
end


end

