IFwing = 1.00;
sHT = 1.6145 ; % Horizontal wing area
sVT = 0.430556417 ; % Vertical wing area
sHTwet = 1.05 * 2 * sHT;
sVTwet = 1.05 * 2 * sVT;

IFtail = 1.05;
FFHT = (1+(0.6/xsem)*atrsem+100*(atrsem^4))*(1.34*(M^0.18)*(cos(0)^0.28));
FFVT = (1+(0.6/xsem)*atrsem+100*(atrsem^4))*(1.34*(M^0.18)*(cos(0)^0.28));
T = 509 ; % Temperature
airspeed = 50; % cruise speed
soundspd = sqrt(1.4*1716*T);
M = airspeed/soundspd;
xcamber = 29.43/100; % Location of the maximum thickness (MH114)
atrcamber = t/c; %airfoil thickness ratio (MH114)
xsem = 30.84/100;% Location of the maximum thickness (NACA 0012)
atrsem = 12/100; %airfoil thickness ratio (NACA 0012)
sWing = 6.13542894 ; % Wing area
wDf = 0.001981 * IFwing * FFwing * sWingwet ; % Wing drag factor
FFwing = (1+(0.6/xcamber)*atrcamber+100*(atrcamber^4))*(1.34*(M^0.18)*(cos(0)^0.28));
x_wing = ((0.012*sWing)-2.7045e-04)/(IFwing * FFwing * sWingwet)
x_VT = ((0.0044*sWing)-2.7045e-04)/(IFtail * FFVT * sVTwet)
x_HT = ((0.0096*sWing)-2.7045e-04)/(IFtail * FFHT * sHTwet)