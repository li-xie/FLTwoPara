function y=mu_factor50(n)
u=rand(n,1);
smax=0.14;
y=-log(1-u)/50;
y(y>smax)=0;








