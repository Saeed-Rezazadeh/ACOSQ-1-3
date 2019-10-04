function Probability = Pr_y_1(y_1, Pr , f , T , delta)
summation = 0 ; 
for x_1 = 1 : 2 
    u_index_x_1 = find(T(: , 2) == x_1) ; 
    summation = summation + Pr(y_1 , x_1) * delta * sum(f(u_index_x_1)) ; 
end 
Probability = summation ;
end 