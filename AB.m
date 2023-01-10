function z = AB(x,theta)

y = x;
y(abs(y)>=theta) = 0;


Rxx = 1/size(x,2)*(x*x');
Ryx = 1/size(x,2)*(y*x');

B = Ryx*(Rxx)^-1;
z = B*x;


