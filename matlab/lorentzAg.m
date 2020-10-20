function f = lorentzAg(w)

eps_inf_Ag = 1;
f_Ag1 = 0.845;
f_Ag2 = 0.065;
f_Ag3 = 0.124;
f_Ag4 = 0.011;
f_Ag5 = 0.84;
f_Ag6 = 5.646;

wo_Ag1 = 0;
wo_Ag2 = 1.2397e+15;
wo_Ag3 = 6.8078e+15;
wo_Ag4 = 1.2435e+16;
wo_Ag5 = 1.38e+16;
wo_Ag6 = 3.0826e+16;
wp_Ag = 1.3689e+16;

y_Ag1 = 7.2925e+13;
y_Ag2 = 5.9039e+15;
y_Ag3 = 6.8671e+14;
y_Ag4 = 9.8752e+13;
y_Ag5 = 1.3916e+15;
y_Ag6 = 3.6751e+15;

f = eps_inf_Ag + f_Ag1*wp_Ag^2/(wo_Ag1^2-w^2-1i*y_Ag1*w) + f_Ag2*wp_Ag^2/(wo_Ag2^2-w^2-1i*y_Ag2*w) + f_Ag3*wp_Ag^2/(wo_Ag3^2-w^2-1i*y_Ag3*w) + f_Ag4*wp_Ag^2/(wo_Ag4^2-w^2-1i*y_Ag4*w) + f_Ag5*wp_Ag^2/(wo_Ag5^2-w^2-1i*y_Ag5*w) + f_Ag6*wp_Ag^2/(wo_Ag6^2-w^2-1i*y_Ag6*w);
end
%eps_inf_Ag,f_Ag1,f_Ag2,f_Ag3,f_Ag4,f_Ag5,f_Ag6,wo_Ag1,wo_Ag2,wo_Ag3,wo_Ag4,wo_Ag5,wo_Ag6,wp_Ag,y_Ag1,y_Ag2,y_Ag3,y_Ag4,y_Ag5,y_Ag6