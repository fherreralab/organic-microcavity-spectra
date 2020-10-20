function f = sellmeier(lambda)
B1 = 0.6961663;
B2 = 0.4079426;
B3 = 0.8974794;
C1 = 0.0684043;
C2 = 0.1162414;
C3 = 9.896161;
f = (1+(B1*lambda^2/(lambda^2-C1^2))+(B2*lambda^2/(lambda^2-C2^2))+(B3*lambda^2/(lambda^2-C3^2)))^0.5;
end