function [SDR , Distortion , T , codebook] = COSQ_4 (Pr_z , f , T ,  y_1 , codebook , delta)

FileID = fopen ('Results.txt' , 'a') ;
D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal partitions
    parfor u_index = 1 : length(T)
        summation = 0 ;
        d_4 = zeros (8 , 1) ;
        u = T(u_index , 1) ;
        x_1 = T(u_index , 2) ;
        for x_2 = 1 : 2
            for x_3 = 1 : 2
                for x_4 = 1 : 2
                    x_prime = (x_2 - 1) * 4 + (x_3 - 1) * 2 + x_4 ;
                    for y_2 = 1 : 2
                        for y_3 = 1 : 2
                            for y_4 = 1 : 2
                                y_prime = (y_2 - 1) * 4 + (y_3 - 1) * 2 + y_4 ;
                                
                                summation = summation + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) ...
                                    * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                                    * Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) * (u - codebook (y_prime)) ^ 2 ;
                            end
                        end
                    end
                    d_4 (x_prime) = summation ;
                    summation = 0 ;
                end
            end
        end
        [~ , partition_index] = min (d_4) ;
        T_u (u_index) = partition_index ;
    end
    T (: , 3) = T_u ;
    
    %% Optimal Centroids
    for y_2 = 1 : 2
        for y_3 = 1 : 2
            for y_4 = 1 : 2
                y_prime = (y_2 - 1) * 4 + (y_3 - 1) * 2 + y_4 ;
                numerator = 0 ;
                denominator = 0 ;
                for x_2 = 1 : 2
                    for x_3 = 1 : 2
                        for x_4 = 1 : 2
                            x_prime = (x_2 - 1) * 4 + (x_3 - 1) * 2 + x_4 ;
                            u_index = find (T (: , 3) == x_prime) ;
                            for u_i = 1 : length(u_index)
                                x_1 = T(u_index(u_i) , 2) ;
                                u = T(u_index(u_i) , 1) ; 
                                numerator = numerator + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) ...
                                    * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                                    * Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) ...
                                    * u * f(u_index(u_i)) ;
                                
                                denominator = denominator + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) ...
                                    * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                                    * Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) ...
                                    * f(u_index(u_i)) ;
                            end
                        end
                    end
                end
                codebook(y_prime) = numerator / denominator ;
            end
        end
    end
    %% Distortion
    D(2) = distortion_4 (f , y_1 , codebook , delta , Pr_z , T ) ;
    fprintf (FileID , 'Overall D_4 = %f\n' ,D(2)) ;
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR_4 = %f\n' , SDR) ;
fclose (FileID) ;
end