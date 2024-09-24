pf= 1600+(30*48);
Ef= (40+81)*10^8;
SigmaF= (200+(100*7))*10^6;
pc= 20+(5*37);
SigmaC= (0.5+(6.5*((37/100)^1.5)))*10^6;
Tc= (0.5+(4.5*((37/100)^1.5)))*10^6;


syms x y lambda
f = y*x +((pc/pf)*x) ;
g = 4*(x^2)*y + 4*(y^2)*(x^2);
L=f-lambda*g
Lx=diff(L,x)==0;
Ly=diff(L,y)==0;
system=[Lx,Ly,g==0]
[xs,ys,lambdas]= solve(system,[x y lambda], 'Real', true)
f(x,y)= y*x +((pc/pf)*x) ;
f(xs,ys)
[min(f(xs,ys)), max(f(xs,ys))]