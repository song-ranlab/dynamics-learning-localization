function B = gauss_newton_sin(x,yn,B)

damp = .01;

% x = 0:.1:6;
% 
% y=0.5*sin(1.2*x+0.3)+0.6;
% 
% yn = y' + (rand(length(y),1) - .5)/10;

%B =  [0.4; 1; 0.4; 0.5];

rows = length(yn);

for interation = 1:10000

    %build J and r
    for row = 1:rows
         xk = x(row);

        J(row,1) = sin(B(2) * xk + B(3));

        J(row,2) = B(1) * xk * cos(B(2) * xk + B(3)); 

        J(row,3) = B(1) * cos(B(2) * xk + B(3)); 

        J(row,4) = 1.0;  

        r(row) =  -yn(row) +  B(1) * sin(B(2) * xk + B(3)) + B(4);
    end   

    Jt = transpose(J);    
    delta = (inv(Jt*J))*Jt*r';   
    B = B - damp*delta;

end