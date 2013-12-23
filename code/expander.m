function [phi epsilon] = expander(rp,rpm,p)

x_0_0=2.098;			% Value of the pressure ratio when epsilon=0, in the reference conditions
dydx_0=0.8886;		    % Slope at epsilon=0, in the reference conditions
shape_0=1.617;          % Shape parameter of the epsilon vs. rp curve
x_m_0=4.45;             % rp corresponding the maximum epsilon. in the reference conditions
y_m_0=0.7037;           % Max epsilon in the reference conditions
rpm_m=2437;             % rpm corresponding to epsilon_max
a0=0.7;
a1=0.6373;
a2=1.628;
a3=0.9768;
a4=2.162;
a5=0.1906;
a6=0.3256;

rpm_star = (rpm - 3000)/3000;
rpm_star_m = (rpm_m - 3000)/3000;
p_star = (p - 10)/10;
r_star_p = (rp-4)/4;

x_0=x_0_0 + a0 * rpm_star;
dydx= dydx_0 + a1*p_star - a2*rpm_star;
shape=shape_0;
x_m=x_m_0 -  a3 * p_star + a4 * rpm_star;
y_max=y_m_0 +a5*p_star -a6*(rpm_star - rpm_star_m).^2;

x=rp;
A = x_0;
C = shape;
D = y_max;
B = dydx./(C*D);
E = (B .* (x_m-x_0) - tan(pi/(2*C)))./(B .* (x_m-x_0) - atan(B .* (x_m-x_0)));
epsilon = D .* sin(C * atan(B.*(x-A) - E .* (B .* (x-A) - atan(B.*(x-A)))));

phi=0.751399733-0.336838575*log(rpm/3000)+0.0239211190*r_star_p-0.0366144375*p_star;

end



