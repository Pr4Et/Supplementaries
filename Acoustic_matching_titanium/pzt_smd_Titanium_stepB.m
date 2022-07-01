%Script to find best thickness of PZT ceramic plate for optimum transmission of ultrasound, given the thickness and properties of the protective layer
%Layers assumed in the current calculation: PZT- Titanium - water
%See the paper: Titanium made ultrasonic reactor with 500KHz resonance
%Written by Shahar Seifer
%(C) copyright 2018, Shahar Seifer


%Code run once to solve symbolically the system of equations
%coefficients of A B C D E and right side of equation
%sym mat matsq d
%mat(1,:)=sym('[-i*k_pzt i*k_pzt,  0,  0,  0,  V*d33/t_pzt]');
%mat(2,:)=sym('[exp(-i*k_pzt*t_pzt),  exp(i*k_pzt*t_pzt),  -exp(-i*k_smd*t_pzt),  -exp(i*k_smd*t_pzt),  0,  0]');
%mat(3,:)=sym('[-i*k_pzt*exp(-i*k_pzt*t_pzt)/s33_pzt,  i*k_pzt*exp(i*k_pzt*t_pzt)/s33_pzt,  i*k_smd*exp(-i*k_smd*t_pzt)/s33_smd,  -i*k_smd*exp(i*k_smd*t_pzt)/s33_smd,  0,  d33*V/t_pzt]'); 
%mat(4,:)=sym('[0,  0,  -i*k_smd*exp(-i*k_smd*(t_pzt+t_smd))/s33_smd,  i*k_smd*exp(i*k_smd*(t_pzt+t_smd))/s33_smd,  +i*k_water*exp(-i*k_water*(t_pzt+t_smd))/s33_water,  0]');
%mat(5,:)=sym('[0,  0,  exp(-i*k_smd*(t_pzt+t_smd)),  exp(i*k_smd*(t_pzt+t_smd)),  -exp(-i*k_water*(t_pzt+t_smd)),  0]');
%matsq=mat(1:5,1:5);
%d=mat(1:5,6);
%ABCDE=matsq^(-1)*d;
%abcde=simplify(ABCDE);

%found abcde =
% 
% (V*d33*(k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i + k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i - k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i)*i - k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i)*i + k_pzt*k_water*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i)*i + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i - k_pzt*k_water*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i)*i + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i)*i - k_pzt*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*i)*i - k_pzt*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i))/(k_pzt*t_pzt*(k_pzt*k_water*s33_smd^2 - k_smd^2*s33_pzt*s33_water - k_pzt*k_smd*s33_smd*s33_water + k_smd*k_water*s33_smd*s33_pzt + k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) - k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i) - k_pzt*k_water*s33_smd^2*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i) + k_smd^2*s33_pzt*s33_water*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i) - k_pzt*k_smd*s33_smd*s33_water*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_smd*t_smd*2*i)))
% -(V*d33*(k_pzt*k_water*s33_smd^2*i - k_smd^2*s33_pzt*s33_water*i - k_pzt*k_smd*s33_smd*s33_water*i + k_smd*k_water*s33_smd*s33_pzt*i - k_pzt*k_water*s33_smd^2*exp(k_smd*t_smd*2*i)*i + k_smd^2*s33_pzt*s33_water*exp(k_smd*t_smd*2*i)*i - k_pzt*k_water*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i)*i + k_pzt*k_water*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i - k_pzt*k_smd*s33_smd*s33_water*exp(k_smd*t_smd*2*i)*i + k_smd*k_water*s33_smd*s33_pzt*exp(k_smd*t_smd*2*i)*i + k_pzt*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*i)*i + k_pzt*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i))/(k_pzt*t_pzt*(k_pzt*k_water*s33_smd^2 - k_smd^2*s33_pzt*s33_water - k_pzt*k_smd*s33_smd*s33_water + k_smd*k_water*s33_smd*s33_pzt + k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) - k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i) - k_pzt*k_water*s33_smd^2*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i) + k_smd^2*s33_pzt*s33_water*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i) - k_pzt*k_smd*s33_smd*s33_water*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_smd*t_smd*2*i)))
% -(V*d33*s33_smd*(k_smd*s33_water + k_water*s33_smd)*(- 2*exp(k_pzt*t_pzt*i + k_smd*t_pzt*i + k_smd*t_smd*2*i) + s33_pzt*exp(k_smd*(t_pzt*i + t_smd*2*i)) + s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_pzt*i + k_smd*t_smd*2*i))*i)/(k_pzt*k_water*s33_smd^2*t_pzt - k_smd^2*s33_pzt*s33_water*t_pzt - k_pzt*k_water*s33_smd^2*t_pzt*exp(k_pzt*t_pzt*2*i) - k_pzt*k_water*s33_smd^2*t_pzt*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i) + k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_smd*t_smd*2*i) - k_pzt*k_smd*s33_smd*s33_water*t_pzt + k_smd*k_water*s33_smd*s33_pzt*t_pzt + k_pzt*k_water*s33_smd^2*t_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i) - k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_pzt*t_pzt*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i))
% (V*d33*k_water*s33_smd^2*s33_pzt*i - V*d33*k_water*s33_smd^2*exp(k_pzt*t_pzt*i)*2*i + V*d33*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*i)*2*i + V*d33*k_water*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*2*i)*i - V*d33*k_smd*s33_smd*s33_pzt*s33_water*i - V*d33*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i)*i)/(k_pzt*k_water*s33_smd^2*t_pzt*exp(k_smd*t_pzt*i) - k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_smd*t_pzt*i) - k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_smd*t_pzt*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_smd*t_pzt*i) - k_pzt*k_water*s33_smd^2*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_pzt*i) - k_pzt*k_water*s33_smd^2*t_pzt*exp(k_smd*t_pzt*i)*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_pzt*i) + k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_smd*t_pzt*i)*exp(k_smd*t_smd*2*i) + k_pzt*k_water*s33_smd^2*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_pzt*i)*exp(k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_pzt*i)*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_pzt*i) - k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_smd*t_pzt*i)*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_pzt*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_smd*t_pzt*i)*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_pzt*i)*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_pzt*i)*exp(k_smd*t_smd*2*i))
%  -(V*d33*k_smd*s33_smd*s33_pzt*s33_water*exp(k_smd*t_smd*i)*exp(k_water*t_pzt*i)*exp(k_water*t_smd*i)*2*i - V*d33*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*i)*exp(k_smd*t_smd*i)*exp(k_water*t_pzt*i)*exp(k_water*t_smd*i)*4*i + V*d33*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*i)*exp(k_water*t_pzt*i)*exp(k_water*t_smd*i)*2*i)/(k_pzt*k_water*s33_smd^2*t_pzt - k_smd^2*s33_pzt*s33_water*t_pzt - k_pzt*k_water*s33_smd^2*t_pzt*exp(k_pzt*t_pzt*2*i) - k_pzt*k_water*s33_smd^2*t_pzt*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i) + k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_smd*t_smd*2*i) - k_pzt*k_smd*s33_smd*s33_water*t_pzt + k_smd*k_water*s33_smd*s33_pzt*t_pzt + k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i) - k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_pzt*t_pzt*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_smd*t_smd*2*i) + k_pzt*k_water*s33_smd^2*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*2*i));

%all variables in mks units
t_pzt_vect=0.001*(3.5:0.025:4.5);%search range of PZT thickness in meters
grand_n=0;
grand_acoustic=zeros(size(t_pzt_vect));
grand_pzt=zeros(size(t_pzt_vect));
grand_peak=zeros(size(t_pzt_vect));
grand_bandwidth=zeros(size(t_pzt_vect));
for t_pzt=0.0040;%t_pzt_vect
    grand_n=grand_n+1;
    V=4; %amplitude of electrical voltage
    d33=265*1e-12; %PZT strngth
    eps33=1200*8.85e-12; %dielectric strength in PZT
    e33=sqrt(0.46*16.6*1e10*eps33);%e33=sqrt(kt*c33D*eps33S)
    s33_pzt=14.2*1e-12;
    t_smd=0.0065; %6.5mm Titanium thickness fixed
    s33_water=1/(2.15e9);  %=1/[bulk modulus], =1/K , in units m^2/Newton.  c=sqrt(K/ro).
    s33_air=1/(1.4e5);
    s33_smd=1/(100e9);  %=1/[bulk modulus], =1/K , in units m^2/Newton. (also the same for young modulus) 
    c_pzt=4220; %Nt*2
    c_water=1480;
    c_air=340;
    ro_smd=4500; %Titanium
    c_smd=6000;%Titanium (Engineering Toolbox)

    As=1e-4; %transducer area in m^2
    ro_water=1000;  %kg/m^3
    ro_air=1.2;%kg/m^3
    ro_pzt=7800;

    n=0;
    freq_vect=(0.4:0.002:0.6);
    Yvect=zeros(size(freq_vect));
    Yvect_air=zeros(size(freq_vect));
    acoustic_power=zeros(size(freq_vect));
    pzt_forward_reverse_power=zeros(size(freq_vect));
    for freq=freq_vect
    %freq- frequency
    n=n+1;    
    k_pzt=(1-0.1i)*2*pi*freq*1000000/c_pzt;
    k_smd=(1+.1i)*2*pi*freq*1000000/c_smd;
    k_water=2*pi*freq*1000000/c_water;
    k_air=2*pi*freq*1000000/c_air;
    omega=2*pi*freq*1000000;
    E= -(V*d33*k_smd*s33_smd*s33_pzt*s33_water*exp(k_smd*t_smd*i)*exp(k_water*t_pzt*i)*exp(k_water*t_smd*i)*2*i - V*d33*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*i)*exp(k_smd*t_smd*i)*exp(k_water*t_pzt*i)*exp(k_water*t_smd*i)*4*i + V*d33*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*i)*exp(k_water*t_pzt*i)*exp(k_water*t_smd*i)*2*i)/(k_pzt*k_water*s33_smd^2*t_pzt - k_smd^2*s33_pzt*s33_water*t_pzt - k_pzt*k_water*s33_smd^2*t_pzt*exp(k_pzt*t_pzt*2*i) - k_pzt*k_water*s33_smd^2*t_pzt*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i) + k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_smd*t_smd*2*i) - k_pzt*k_smd*s33_smd*s33_water*t_pzt + k_smd*k_water*s33_smd*s33_pzt*t_pzt + k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i) - k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_pzt*t_pzt*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_smd*t_smd*2*i) + k_pzt*k_water*s33_smd^2*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*t_pzt*exp(k_pzt*t_pzt*2*i)*exp(k_smd*t_smd*2*i));
    A=(V*d33*(k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i + k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i - k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i)*i - k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i)*i + k_pzt*k_water*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i)*i + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i - k_pzt*k_water*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i)*i + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i)*i - k_pzt*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*i)*i - k_pzt*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i))/(k_pzt*t_pzt*(k_pzt*k_water*s33_smd^2 - k_smd^2*s33_pzt*s33_water - k_pzt*k_smd*s33_smd*s33_water + k_smd*k_water*s33_smd*s33_pzt + k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) - k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i) - k_pzt*k_water*s33_smd^2*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i) + k_smd^2*s33_pzt*s33_water*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i) - k_pzt*k_smd*s33_smd*s33_water*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_smd*t_smd*2*i)));
    B=-(V*d33*(k_pzt*k_water*s33_smd^2*i - k_smd^2*s33_pzt*s33_water*i - k_pzt*k_smd*s33_smd*s33_water*i + k_smd*k_water*s33_smd*s33_pzt*i - k_pzt*k_water*s33_smd^2*exp(k_smd*t_smd*2*i)*i + k_smd^2*s33_pzt*s33_water*exp(k_smd*t_smd*2*i)*i - k_pzt*k_water*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i)*i + k_pzt*k_water*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i - k_pzt*k_smd*s33_smd*s33_water*exp(k_smd*t_smd*2*i)*i + k_smd*k_water*s33_smd*s33_pzt*exp(k_smd*t_smd*2*i)*i + k_pzt*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*i)*i + k_pzt*k_smd*s33_smd*s33_pzt*s33_water*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i))/(k_pzt*t_pzt*(k_pzt*k_water*s33_smd^2 - k_smd^2*s33_pzt*s33_water - k_pzt*k_smd*s33_smd*s33_water + k_smd*k_water*s33_smd*s33_pzt + k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) - k_pzt*k_water*s33_smd^2*exp(k_pzt*t_pzt*2*i) - k_pzt*k_water*s33_smd^2*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_water*exp(k_pzt*t_pzt*2*i) + k_smd^2*s33_pzt*s33_water*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_water*exp(k_pzt*t_pzt*2*i) - k_pzt*k_smd*s33_smd*s33_water*exp(k_smd*t_smd*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i) + k_smd*k_water*s33_smd*s33_pzt*exp(k_smd*t_smd*2*i)));
    Y=1i*omega*As*(eps33/t_pzt+(e33/V)*(-1i*k_pzt*A*exp(-1i*k_pzt*t_pzt)+1i*k_pzt*B*exp(1i*k_pzt*t_pzt)));
    Yvect(n)=Y;
    acoustic_power(n)=0.5*As*ro_water*c_water*abs(E*conj(E))*(omega)^2;
    pzt_forward_reverse_power(n)=0.5*As*ro_pzt*c_pzt*(abs(A*conj(A))-abs(B*conj(B)))*(omega)^2; %forward-reverse of waves inside pzt
    %special case in air
    A_air=(V*d33*(k_pzt*k_air*s33_smd^2*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i + k_smd^2*s33_pzt*s33_air*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i - k_pzt*k_air*s33_smd^2*exp(k_pzt*t_pzt*2*i)*i - k_smd^2*s33_pzt*s33_air*exp(k_pzt*t_pzt*2*i)*i + k_pzt*k_air*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i)*i + k_pzt*k_smd*s33_smd*s33_air*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i + k_smd*k_air*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i)*i - k_pzt*k_air*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i + k_pzt*k_smd*s33_smd*s33_air*exp(k_pzt*t_pzt*2*i)*i + k_smd*k_air*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i)*i - k_pzt*k_smd*s33_smd*s33_pzt*s33_air*exp(k_pzt*t_pzt*i)*i - k_pzt*k_smd*s33_smd*s33_pzt*s33_air*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i))/(k_pzt*t_pzt*(k_pzt*k_air*s33_smd^2 - k_smd^2*s33_pzt*s33_air - k_pzt*k_smd*s33_smd*s33_air + k_smd*k_air*s33_smd*s33_pzt + k_pzt*k_air*s33_smd^2*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_air*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) - k_pzt*k_air*s33_smd^2*exp(k_pzt*t_pzt*2*i) - k_pzt*k_air*s33_smd^2*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_air*exp(k_pzt*t_pzt*2*i) + k_smd^2*s33_pzt*s33_air*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_air*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd*k_air*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_air*exp(k_pzt*t_pzt*2*i) - k_pzt*k_smd*s33_smd*s33_air*exp(k_smd*t_smd*2*i) + k_smd*k_air*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i) + k_smd*k_air*s33_smd*s33_pzt*exp(k_smd*t_smd*2*i)));
    B_air=-(V*d33*(k_pzt*k_air*s33_smd^2*i - k_smd^2*s33_pzt*s33_air*i - k_pzt*k_smd*s33_smd*s33_air*i + k_smd*k_air*s33_smd*s33_pzt*i - k_pzt*k_air*s33_smd^2*exp(k_smd*t_smd*2*i)*i + k_smd^2*s33_pzt*s33_air*exp(k_smd*t_smd*2*i)*i - k_pzt*k_air*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i)*i + k_pzt*k_air*s33_smd^2*s33_pzt*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i - k_pzt*k_smd*s33_smd*s33_air*exp(k_smd*t_smd*2*i)*i + k_smd*k_air*s33_smd*s33_pzt*exp(k_smd*t_smd*2*i)*i + k_pzt*k_smd*s33_smd*s33_pzt*s33_air*exp(k_pzt*t_pzt*i)*i + k_pzt*k_smd*s33_smd*s33_pzt*s33_air*exp(k_pzt*t_pzt*i + k_smd*t_smd*2*i)*i))/(k_pzt*t_pzt*(k_pzt*k_air*s33_smd^2 - k_smd^2*s33_pzt*s33_air - k_pzt*k_smd*s33_smd*s33_air + k_smd*k_air*s33_smd*s33_pzt + k_pzt*k_air*s33_smd^2*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd^2*s33_pzt*s33_air*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) - k_pzt*k_air*s33_smd^2*exp(k_pzt*t_pzt*2*i) - k_pzt*k_air*s33_smd^2*exp(k_smd*t_smd*2*i) - k_smd^2*s33_pzt*s33_air*exp(k_pzt*t_pzt*2*i) + k_smd^2*s33_pzt*s33_air*exp(k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_air*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_smd*k_air*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i + k_smd*t_smd*2*i) + k_pzt*k_smd*s33_smd*s33_air*exp(k_pzt*t_pzt*2*i) - k_pzt*k_smd*s33_smd*s33_air*exp(k_smd*t_smd*2*i) + k_smd*k_air*s33_smd*s33_pzt*exp(k_pzt*t_pzt*2*i) + k_smd*k_air*s33_smd*s33_pzt*exp(k_smd*t_smd*2*i)));
    Y_air=1i*omega*As*(eps33/t_pzt+(e33/V)*(-1i*k_pzt*A_air*exp(-1i*k_pzt*t_pzt)+1i*k_pzt*B_air*exp(1i*k_pzt*t_pzt)));
    Yvect_air(n)=Y_air;
    end %for freq

    figure(1)
    plot(freq_vect,acoustic_power,'-')
    xlabel('Frequency [MHz]');
    ylabel('Acoustic intensity [W/cm^2]');
    title(sprintf('t_{pzt}=%.2f  mm',t_pzt*1000));

    figure(11)
    plot(freq_vect,abs(real(Yvect)),'g-',freq_vect,abs(real(Yvect_air)),'r-')
    xlabel('Frequency [MHz]');
    ylabel('|Real(Y)|  [Simense]');
    legend('water','air')
    title(sprintf('t_{pzt}=%.2f  mm',t_pzt*1000));

    figure(13)
    [ax h1 h2]=plotyy(freq_vect,acoustic_power,freq_vect,abs(real(Yvect)));
    xlabel('Frequency [MHz]');
    axes(ax(1)); ylabel('Acoustic power [W]');
    axes(ax(2)); ylabel('Electrical G  [Simense]');
    title(sprintf('Simulation %g',t_pzt*1000));


    pause (1)
    close(13)
    grand_acoustic(grand_n)=max(acoustic_power);
    grand_pzt(grand_n)=max(pzt_forward_reverse_power);
    grand_peak(grand_n)=mean(freq_vect(acoustic_power==max(acoustic_power)));
    tentative_vector1=freq_vect(freq_vect>grand_peak(grand_n) &  (acoustic_power<=0.5*max(acoustic_power)));
    tentative_vector2=freq_vect(freq_vect<grand_peak(grand_n) &  (acoustic_power<=0.5*max(acoustic_power)));
    if ~isempty(tentative_vector1)
        df1=min(tentative_vector1)-grand_peak(grand_n);
    else
        df1=10;
    end
    if ~isempty(tentative_vector2)
        df2=grand_peak(grand_n)-max(tentative_vector2);
    else
        df2=10;
    end
    grand_bandwidth(grand_n)=min([df1 df2]);
end % for t_smd


figure(5)
plot([t_pzt_vect*1000],[grand_acoustic],'-')
xlabel('PZT thickness [mm]');
ylabel('Peak acoustic intensity [W /cm^2]');

figure(55)
plot([t_pzt_vect]*k_smd/(2*pi),[grand_acoustic],'-')
xlabel('PZT thickness [wavelength number]');
ylabel('Peak acoustic intensity [W /cm^2]');

figure(6)
plot(t_pzt_vect*1000,grand_peak,'-')
xlabel('PZT thickness [mm]');
ylabel('Peak frequency [MHz]');

%figure(9)
%plot(t_pzt_vect*1000,grand_bandwidth,'-')
%xlabel('PZT thickness [mm]');
%ylabel('Bandwidth [MHz]');
