function sum = trapez(y,x)
% Trapezoidal integration 3/8
% Only for any sort of grid spacing

sum = 0;
Nx = max(size(x));

for i=1:Nx-1;
    
    sum = sum + (x(i+1)-x(i))*(y(i+1)+y(i))/2;

end

