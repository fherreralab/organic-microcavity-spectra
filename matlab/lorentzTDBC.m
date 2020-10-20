function f = lorentzTDBC(w)
eps_b_organic3 = 2;
fo_organic3 = 0.2;
f1_organic3 = 0.03;

wo_organic3 = 3.236098450319052e+15;
w1_organic3 = 3.5706e15;

y1_organic3 = 3.038590094196293e+14;
yo_organic3 = 7.596475235490733e+13;
f = eps_b_organic3 + fo_organic3*wo_organic3^2/(wo_organic3^2 - w^2 - 1i*yo_organic3*w) + f1_organic3*w1_organic3^2/(w1_organic3^2 - w^2 - 1i*y1_organic3*w);
end