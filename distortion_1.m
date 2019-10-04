function Distortion = distortion_1(f , codebook , delta , Pr_1 , T)
summation = 0 ; 
parfor x_1 = 1 : 2 
    for y_1 = 1 : 2 
        u_index = find (T(: , 2) == x_1) ; 
        u = T(u_index , 1) ; 
                
        summation = summation + Pr_1(x_1 , y_1) * delta * sum(f(u_index).* (u - codebook(y_1)) .^ 2) ; 
        
    end 
end 
Distortion = summation ; 
end 