function [bayes_inputs] = create_model_prior_and_covariance(bayes_inputs,truthTable, use_MODIS_estimates, modis, vocalsRex)

% -------------------------------------------------------------
% -------------------------------------------------------------
% THIS CODE ASSUMES THE RETIREVAL VARIABLES ARE ALL INDEPENDENT
% -------------------------------------------------------------
% -------------------------------------------------------------


% define the model variance and mean using the Truth Table found by my
% TBLUT algorithm

if bayes_inputs.numPixels2Calculate<=size(truthTable,1)

    % Pick the first n pixels in the table
    n = bayes_inputs.numPixels2Calculate;

    % the order of the values below: (r_top, r_bottom, tau_c)
    % The model mean is the a priori, which is separate from our initial
    % guess.


    if use_MODIS_estimates==false

        %--------------------------------------------
        % use my own TBLUT estimates to define priors
        %--------------------------------------------




        %----------------------------------------------------
        % ----------- Set the a priori value ----------------
        %----------------------------------------------------

        bayes_inputs.model.apriori = [1.5*truthTable.estR17(1:n), 0.5*truthTable.estR17(1:n), truthTable.estT17(1:n)]; % expected values for the effective radius (microns) and the optical depth
        %bayes_inputs.model.apriori = [truthTable.estR17(1:n), truthTable.estR17(1:n), truthTable.estT17(1:n)];

        % lets create the variance and mean for each model parameter
        % Using the same values defined by King and Vaughn (2012)
        % King and Vaughn define the standard deviation of each variable
        % The first two values are the standard deviation of the effective
        % radius at the top of the cloud and the bottom of the cloud, measured
        % in microns. The third value is the percentage of the optical depth
        % that defines the standard deviation.
        stdev_variables = [sqrt(3), sqrt(10), sqrt(0.3)];

        

        
        bayes_inputs.model.variance = [linspace(stdev_variables(1)^2,stdev_variables(1)^2,n)',...
            linspace(stdev_variables(2)^2,stdev_variables(2)^2,n)',...
            stdev_variables(3)^2 *truthTable.estT17(1:n)]; % variance for the effective radius (microns squared) and optical thickness respectively

    

        %----------------------------------------------------
        % ----------- Set the Initial guess  ----------------
        %----------------------------------------------------

        % Define the initial guess as being similar to the a priori except that
        % we define the initial guess as having the same value for the
        % effective radius at cloud top and bottom
        %bayes_inputs.model.initialGuess = [1.5*truthTable.estR17(1:n), 0.5*truthTable.estR17(1:n), truthTable.estT17(1:n)];
        bayes_inputs.model.initialGuess = [truthTable.estR17(1:n), truthTable.estR17(1:n), truthTable.estT17(1:n)];



        %----------------------------------------------------------
        % ----------- Define the Covariance Matrix ----------------
        %----------------------------------------------------------

        % For now lets claim the desired variables are independent
        for ii = 1:n
            bayes_inputs.model.covariance(:,:,ii) = diag(bayes_inputs.model.variance(ii,:));
        end


    



    
    else

        %---------------------------------------------------
        % use MODIS retrievals for initial guess and priori
        %---------------------------------------------------


        %----------------------------------------------------
        % ----------- Set the a priori value ----------------
        %----------------------------------------------------

        bayes_inputs.model.apriori = [1.5*truthTable.modisR17(1:n), 0.5*truthTable.modisR17(1:n), truthTable.modisT17(1:n)]; % expected values for the effective radius (microns) and the optical depth
        %bayes_inputs.model.apriori = [truthTable.modisR17(1:n), truthTable.modisR17(1:n), truthTable.modisT17(1:n)];


        % lets create the variance and mean for each model parameter
        % Using the same values defined by King and Vaughn (2012)
        % King and Vaughn define the standard deviation of each variable...
        % The first two values are the standard deviation of the effective
        % radius at the top of the cloud and the bottom of the cloud, measured
        % in microns. The third value is the percentage of the optical depth
        % that defines the standard deviation.
        %stdev_variables = [sqrt(3), sqrt(10), sqrt(0.1 *truthTable.modisT17(1:n))];
        stdev_variables = [sqrt(0.75), sqrt(2.25), sqrt(0.05 *truthTable.modisT17(1:n))];

        bayes_inputs.model.variance = [linspace(stdev_variables(1)^2,stdev_variables(1)^2,n)',...
            linspace(stdev_variables(2)^2,stdev_variables(2)^2,n)',...
            stdev_variables(3).^2 ]; % variance for the effective radius (microns squared) and optical thickness respectively





        %----------------------------------------------------
        % ----------- Set the Initial guess  ----------------
        %----------------------------------------------------

        % Define the initial guess as being similar to the a priori except that
        % we define the initial guess as having the same value for the
        % effective radius at cloud top and bottom
        
        %bayes_inputs.model.initialGuess = [1.5*truthTable.modisR17(1:n), 0.5*truthTable.modisR17(1:n), truthTable.modisT17(1:n)];
        bayes_inputs.model.initialGuess = [truthTable.modisR17(1:n), truthTable.modisR17(1:n), truthTable.modisT17(1:n)];


        %----------------------------------------------------------
        % ----------- Define the Covariance Matrix ----------------
        %----------------------------------------------------------
        
        % For now lets claim the desired variables are independent
        for ii = 1:n
            bayes_inputs.model.covariance(:,:,ii) = diag(bayes_inputs.model.variance(ii,:));
        end


        %------------------------------------------------------
        % ----------- Define cloud models inputs --------------
        %------------------------------------------------------

        % Define a custom cloud depth
        %bayes_inputs.model.cloudDepth = 0.5;            % km

        % Define cloud depth using Vocals Rex
        bayes_inputs.model.cloudDepth = (vocalsRex.altitude(end) - vocalsRex.altitude(1))/1e3;      % km

        % Define cloud top height using Vocals Rex
        bayes_inputs.model.cloudTop_height = vocalsRex.altitude(end)/1e3;           % km

        % Define cloud top height using MODIS data
        %bayes_inputs.model.cloudTop_height = modis.cloud.topHeight(pixel_row, pixel_col)/1e3;        % km


        % Define number of layers to use in libRadTran when defining
        % vertically inhomogenous clouds
        bayes_inputs.model.cloud_layers = 10;



    end





else

    % We will only use what is in the truth table!
    bayes_inputs.numPixels2Calculate = size(truthTable,1);
    n = size(truthTable,1);
    % the order of the values below: (r_top, r_bottom, tau_c)

    %----------------------------------------------------
    % ----------- Set the a priori value ----------------
    %----------------------------------------------------

    % The model mean is the a priori; our first guess
    bayes_inputs.model.apriori = [truthTable.estR17, truthTable.estR17, truthTable.estT17]; % expected values for the effective radius (microns) and the optical depth


    % lets create the variance and mean for each model parameter
    % Using the same values defined by King and Vaughn (2012)
    r_top_var = 3;
    r_bot_var = 10;
    percent_tau = 0.05;

    bayes_inputs.model.variance = [linspace(r_top_var,r_top_var,n)',linspace(r_bot_var,r_bot_var,n)',percent_tau*truthTable.estT17]; % variance for the effective radius (microns squared) and optical thickness respectively


    % For now lets claim the desired variables are independent
    for ii = 1:n

        bayes_inputs.model.covariance(:,:,ii) = diag(bayes_inputs.model.variance(ii,:));
    end


end


end

