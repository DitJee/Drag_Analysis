%% Drag Estimation

%{
 ALL UNITS ARE IN IMPERIAL 
 WRITTEN BY THANAPOL SUPITAYAKUL
 AND DIT DEJPHACHON
 THIS CODE IS DEDICATE TO DBF(2020)
%}

%% USER INPUT
airspeed = 50; % cruise speed
viscosity = 3.745*10^-7; % viscosity of air
density = 0.002378; % air density
dpressure = 0.5 * density * airspeed^2;
T = 509 ; % Temperature
sWing = 6.13542894 ; % Wing area
sHT = 1.6145 ; % Horizontal wing area
sVT = 0.430556417 ; % Vertical wing area
c = 1.24671916  ; % Chord length[ft]
t = 13.06*c/100; % Maximum thickness of the airfoil[ft]
xcamber = 29.43/100; % Location of the maximum thickness (MH114)
atrcamber = t/c; %airfoil thickness ratio (MH114)
xsem = 30.84/100;% Location of the maximum thickness (NACA 0012)
atrsem = 12/100; %airfoil thickness ratio (NACA 0012)
D = 0.45/1.5     ; %Fuselage Diameter
NL = 0.83333 ; %Nose section Length
CL = 1.66666667 ; %Center section Length
TL = 1.66666667/2 ; %Tail section Length
dWheel = 0.377627953; %wheel diameter[ft]
wWheel = 0.123645013; %wheel width[ft]
CDs = 0.615;

D_pass = (3/12)*2; %passenger compartment diameter[ft]
CL_pass = 13/12; %passenger compartment center section lenght[ft]
NL_pass = 2/12; %passenger Nose section Length[ft]

CMGC = 1.2467191601;%wing mean geometric chord[ft]
CmW = 0.05;%wing pitching moment coefficient
AR_wing = 3.95;%wing aspect ratio
e = 1.78*(1 - 0.045*(AR_wing)^(0.68))-0.64;%estimation of Oswald's span efficiency
k = 1/(pi*AR_wing*e);%lift-induced drag constant
hAC = 0.25*(CMGC);%distance between the wing LE and AC
h = 0.25*(CMGC); %distance between wing and CG
lHT = 3.467847769; %distance between wing and HT[ft]
W = 7.7161791765;%weight at condition[lbs]
T = 8.3555197368;%engine thrust[lbf]
zT = 0.16404199475;%distance between the CG and thrustline.
Pitch_prop = 0.91666666667;%propeller pitch
D_prop = 1.3333333333;%prop diameter
volume = ((pi*D_prop^2)/4)*Pitch_prop;
KV = 560;
m_dot = volume*density*KV/60;

%% Atmospheric data
u = symunit;
soundspd = sqrt(1.4*1716*T);
M = airspeed/soundspd;

%% Wetted area
sWingwet = 1.07 * 2 * sWing;
sHTwet = 1.05 * 2 * sHT;
sVTwet = 1.05 * 2 * sVT;

%% Fuselage
sFuse = pi * D /4 * (( 1/(3*(NL^2)) * ((4*NL^2+((D^2)/4))^1.5) - ((D^3)/8))- D + 4*CL + 2*sqrt((TL^2)+((D^2)/4))) 
sFusewet = sFuse*1.25;

%% Passenger
sPass = pi * D_pass /4 * (( 1/(3*(NL_pass^2)) * ((4*NL_pass^2+((D_pass^2)/4))^1.5) - ((D_pass^3)/8)- D_pass + 4*CL_pass + 2*sqrt((TL^2)+(D_pass^2)/4))) ;
sPasswet = sPass*1.25;

%% Interference Factor
IFfuse = 1.0;
IFtail = 1.05;
IFLandinggear = 1.10;
IFwing = 1.00;

%% Landing gear
cdwheel = ((dWheel*wWheel)/sWing)*CDs*IFLandinggear;

%% Reynold Numbers
Reroot = (density * airspeed * c)/viscosity;
ReCutoff = 38.21*(4.875/1.7*10^-6)^1.053;

%% Form Factor
FFwing = (1+(0.6/xcamber)*atrcamber+100*(atrcamber^4))*(1.34*(M^0.18)*(cos(0)^0.28));
FFHT = (1+(0.6/xsem)*atrsem+100*(atrsem^4))*(1.34*(M^0.18)*(cos(0)^0.28));
FFVT = (1+(0.6/xsem)*atrsem+100*(atrsem^4))*(1.34*(M^0.18)*(cos(0)^0.28));
f = (NL+CL+TL)/D; %finese ratio
FFfuse = 1+(60/f^3)+(f/(400));

%% CRUD trim drag

Mw = dpressure*sWing*CMGC*CmW;
A = k/(dpressure*sWing)^2;
B = 1-(h-hAC)/(lHT);
C = W-((Mw-(T*zT))/(lHT-h+hAC));
DeltaCDtrim = A*(B^2)*(C^2) - A*(W^2);


%% Weighted drag factor
 wDf =  0.0056 * IFwing * FFwing * sWingwet ; % Wing drag factor
 fDf = 0.002112 * IFfuse * FFfuse * sFusewet ; % Fuselage drag factor
 hDf = 0.0118* IFtail * FFHT * sHTwet ; % Horizontal drag factor
 vDf = 0.0202 * IFtail * FFVT * sVTwet ; % Vertical drag factor
 pDf = 0.002112 * IFfuse * FFfuse * sPasswet; % Passenger drag factor
 
 disp('pDf')
 %%disp((pDf+dragDelta)/sWing)
 
%% Combined drag
drag = (wDf+fDf+hDf+vDf);
dragDelta = (wDf+fDf+hDf+vDf)*DeltaCDtrim;
disp('DeltaCD')
disp(DeltaCDtrim)
dragT = drag+dragDelta;
cd = (dragT/sWing)+cdwheel;

%% Display
disp('Wing')
disp((wDf+dragDelta)/sWing)
disp('Fuselage')
disp((fDf+dragDelta)/sWing)
disp('Horizontal')
disp((hDf+dragDelta)/sWing)
disp('Vertical')
disp((vDf+dragDelta)/sWing)
disp('Overall Cd')
disp(cd)
disp('Cd induced')
disp('0.014')
Cdt = cd + 0.014;
disp('Total drag coeff')
disp(Cdt)
disp('drag?')
disp(drag)
disp('sWingwet')
disp(sWingwet)
disp(xcamber)
disp(cdwheel)



%% Banner sizing
sThrust_kg = 3.790;
sThrust_N = sThrust_kg*9.81;
disp(sThrust_N)
dThrust_N = sThrust_N*0.5;
disp('dThrust_N')
disp(dThrust_N)
aThrust_N = (dThrust_N)-((Cdt*dpressure*sWing)/2);
disp('aThrust_N')
disp(aThrust_N)
aThrust_lblF = aThrust_N*0.224808943;
AR = 5;

cdBanner = 0.405*(AR^(-0.494))
disp('cdbanner')
disp(cdBanner)
SBanner = (aThrust_lblF/((cdBanner)*dpressure)); %{ft^2}
bBanner = (5*SBanner)^(1/2);
disp('SBanner')
disp(SBanner)
disp('bBanner')
disp(bBanner)

%% Total CD0

M1_CDT0 = ((wDf+dragDelta)/sWing) + ((fDf+dragDelta)/sWing) + ((hDf+dragDelta)/sWing) + ((vDf+dragDelta)/sWing)
M2_CDT0 = ((wDf+dragDelta)/sWing) + ((fDf+dragDelta)/sWing) + ((hDf+dragDelta)/sWing) + ((vDf+dragDelta)/sWing) + ((pDf+dragDelta)/sWing)
M3_CDT0 = ((wDf+dragDelta)/sWing) + ((fDf+dragDelta)/sWing) + ((hDf+dragDelta)/sWing) + ((vDf+dragDelta)/sWing) + (cdBanner)