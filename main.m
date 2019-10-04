%% The script corresponds to the Algorithm 4 with r = (1 3)
clc
clear all
close all
%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon.
% Also, since the ACOSQ is an iterative algorithm, the distirtion
% value for every iteration is provided in this file for a given epsilon and delta.
FileID = fopen ('Results.txt' , 'a') ;

%% Source distributions
sigma = 1 ;  % the standard devision of the source
u = 0 ; % the source's mean value
alpha = 300000  ; % the size of the training set to find the initial codebook using the splitting algorithm

%% Channel's cross-over probability epsilon
epsilon = unique ([10 ^-6 10^-5 : 2 * 10^-5 : 10^-4 , 10^-4 : 10^-4 : 10^-3 , 10^-3 , 0.005 0.01 0.05 0.1]);
% Since the design converges to a locally optimal solution, to avoid
% bad local optimums, we use a increase-decrease method.

SIZE = length(epsilon) ;
i_noise =  [1 : SIZE , SIZE : -1 : 1 , 1 : SIZE , SIZE : -1 : 1 , 1 : SIZE , SIZE : -1 : 1 ] ;

% The variable resolution determines the accuracy of the Riemann summation
resolution = 2 ^ 11 ;


%% Initialize the parameters
SDR_1 = zeros (length(i_noise) , 1) ;
final_SDR_1 = zeros(SIZE , 1) ;

SDR_4 = zeros (2 , 1) ;
Final_SDR_4 = zeros (length(i_noise) , 1) ;
final_SDR_4 = zeros(SIZE , 1) ;

D_4 = zeros(2 , 1) ;


Probability_y_1 = zeros (2 , 1) ;

[Training_set , T , delta_u] = initialization (sigma , u , alpha , resolution) ;

% Compute the source pdf. We herein consider a zero-mean
% unit-variance Gaussian source distribution.
u = T(: , 1) ;
f =  1 ./ (sqrt (2 .* pi)) .* exp (-u .^ 2 ./ 2) ;
f = f ./( sum(f) * delta_u ) ;

%% noise correlation
for delta = [0 5 10]
    D_4 = zeros(2 , 1) ;
    for k = 1: length(i_noise)
        i = i_noise(k) ;
        
        Pr_1 = [1 - epsilon(i) , epsilon(i) ;
            epsilon(i) , 1 - epsilon(i)] ;
        
        %% COSQ for rate 1
        if (k == 1)
            % Using the splitting algorithm to find the initial codebook.
            [~  , codebook] = kmeans (Training_set , 2 , 'MaxIter',1000 , 'OnlinePhase','on') ;
        else
            % We slightly increase the channel's cross-over probability,
            % setting the codebook from the system with small epsilon as
            % the initial state of the system with new epsilon.
            codebook = load ('codebook_rate_1' , 'codebook_1') ;
            codebook = struct2array(codebook) ;
        end
        % The first step of the ACOSQ described in Section 4.1 and
        % Algorithm 4. In this step a 1-bit COSQ is designed for the source pdf f as computed in line 47.
        [SDR_1(k) , D_1 , T , codebook_1] = COSQ_1(f , Pr_1 , T(: , 1) , codebook , delta_u ) ;
        
        % save the codebook to initialize the system with the next value of
        % epsilon. In other words, the codebook obtained at every step is used to initialize
        % the quantizer in the SAME step for the new epsilon.
        save ('codebook_rate_1' , 'codebook_1' ) ;
        
        %% COSQ for rate 4
        codebook_4 = [] ;
        for y_1 = 1 : 2
            % Based on (4.8) compute the source conditional pdf given the received
            % sequence y_1 where y_1 is the channel output corresponding to
            % the transmitted sequence in the first step.
            [f_u_given_y_1] = generate_pdf_rate_4(Pr_1 , f , T , y_1 , delta_u ) ;
            for inner_noise = 1 : k
                i = i_noise(inner_noise) ;
                Pr_z = [(1 - epsilon(i) + delta) / (1 + delta)  , epsilon(i) / (1 + delta) ;
                    (1 - epsilon(i)) / (1 + delta)  , (epsilon(i) + delta) / (1 + delta)] ;
                
                if (inner_noise == 1)
                    % Find the initial codebook for the conditional source
                    % pdf using splitting algorithm.
                    codebook = init_codebook(f_u_given_y_1 , delta_u , T , alpha) ;
                else
                    % We slightly increase the channel's cross-over probability,
                    % setting the codebook from the system with small epsilon as
                    % the initial state of the system with new epsilon.
                    Data = ['codebook_rate_4_y_1_' num2str(y_1)] ;
                    load (Data) ;
                end
                
                % The second step of the ACOSQ described in Section 4.1 and
                % Algorithm 4. In this step a 3-bit COSQ is designed for the conditional source pdf f_u_given_y_1
                % as computed in line 84.
                [SDR_4(y_1) ,  D_4(y_1) , hold_T , codebook] = COSQ_4(Pr_z , f_u_given_y_1 , T(: , 1 : 2) ,  y_1 , codebook , delta_u) ;
            end
            Data = ['codebook_rate_4_y_1_' num2str(y_1)] ;
            save (Data , 'codebook') ;
            
            T(: , 2 + y_1) = hold_T(: , 3) ;
            % Compute the probability of P(Y_1 = y_1) for y_1 = 0 , 1.
            Probability_y_1(y_1) = Pr_y_1(y_1, Pr_1 , f , T , delta_u) ;
            codebook_4 = [codebook_4 ; codebook] ;
        end
        Probability_y_1 = Probability_y_1 ./ (sum(Probability_y_1)) ;
        % As mentioned in the Thesis, the ultimate distortion at every step is the
        % weighted sum of the conditional distortion functions given the
        % received sequence y_1.
        Final_D_4 = sum(D_4 .* Probability_y_1) ;
        % Compute the SDR value in the following way as the source is zero
        % mean and unit variance.
        Final_SDR_4(k) = 10 * log10(1 / Final_D_4) ;
        fprintf (FileID , '\nFinal SDR_4 = %7.4f\n' , Final_SDR_4(i)) ;
        % Store the ultimate partition indexes and codebook to compute
        % experimental results.
        Data = ['T\T_k_' num2str(k) '_noise_Delta_' num2str(delta)] ;
        save(Data, 'T' , 'codebook_4' , 'codebook_1') ;
    end
    % Pick the best SDR value after the end of the so-called
    % increase-decrease method.
    fprintf (FileID , '\nSDR for rate 1\n') ;
    clear D_1
    D_1 = zeros (SIZE , 1) ;
    for i = 1 : SIZE
        index = find (i_noise == i) ;
        hold_var = SDR_1(index) ;
        hold_var = hold_var(:) ;
        final_SDR_1(i) = max(hold_var) ;
        fprintf (FileID , '\ni = %d' , i) ;
        D_1(i) = 10 ^(- final_SDR_1(i) / 10) ;
        fprintf (FileID , '\nD = %f' , D_1(i)) ;
        fprintf (FileID , '\nSDR_1 = %7.4f' , final_SDR_1(i)) ;
    end
    
    fprintf (FileID , '\nSDR for rate 4\n') ;
    clear D_4
    D_4 = zeros(SIZE , 1) ;
    for i = 1 : SIZE
        index = find (i_noise == i) ;
        hold_var = Final_SDR_4(index) ;
        hold_var = hold_var(:) ;
        final_SDR_4(i) = max(hold_var) ;
        fprintf (FileID , '\ni = %d' , i) ;
        D_4(i) = 10 ^(- final_SDR_4(i) / 10) ;
        fprintf (FileID , '\nD = %f' , D_4(i)) ;
        fprintf (FileID , '\nSDR_4 = %7.4f\n' , final_SDR_4(i)) ;
    end
    %  Save the best SDR values for every channel parameters.
    clear D_4 ;
    Data = ['ACOSQ_1_3_delta_' num2str(delta)] ;
    save(Data , 'final_SDR_4' , 'Final_SDR_4' , 'final_SDR_1' , 'epsilon' , 'SDR_1') ;
end
