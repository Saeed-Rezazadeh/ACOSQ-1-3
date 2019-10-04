function [f_u_given_y_1] = generate_pdf_rate_4(Pr_1 , f , T , y_1 , delta )
%% numerator
numerator = zeros (length(T) , 1) ;
for u_index = 1 : length(T)
    x_1 = T(u_index , 2) ;
    numerator(u_index) = Pr_1(x_1 , y_1) * f(u_index) ;
end

%% denominator
denominator = 0 ;
for x_1 = 1 : 2
    u_index_x_1 = find (T(: , 2) == x_1) ;
    denominator = denominator + Pr_1(x_1 , y_1) * delta * sum(f(u_index_x_1)) ;
end
f_u_given_y_1 = numerator ./ denominator ;
f_u_given_y_1 = f_u_given_y_1 ./ (sum(f_u_given_y_1) .* delta) ;
end