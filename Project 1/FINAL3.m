%%Question 4
c=0:0.000001:0.000225;
t=0.000070248-(0.325.*c);
plot(c,t)
title('Optimal Design of Example Beam c vs t')
xlabel c(m)
ylabel t(m)