function [w,e,de] = stkerr()

%Load data
[fname,pname] = uigetfile('.mat','Select Stark Spectrum');
if fname == 0, return; end

S = load([pname,fname]);

w = S.i(1,:)';
I = S.i(2,:)';
I0 = S.i0(2,:)';
dI = S.istd(2,:)';
dI0 = S.i0std(2,:)';

%Spec parameters
prompt = {'Field (Vrms):','Pathlength (um):','Concentration (mM):','Incidence (deg):'};
dlg_title = 'Stark parameters';
num_lines = 1;
def = {'1300','55','0.845','45'};
answer = inputdlg(prompt,dlg_title,num_lines,def);


Vrms = str2double(answer{1});
b0 = str2double(answer{2})*1e-4;
c = str2double(answer{3})*1e-3;
inc = str2double(answer{4});


% Vrms = 1300;
% b0 = 55e-4;
% c = .835*1e-3;
% inc = 45;

%pathlength correction
n1 = 1.2; %LN2
n2 = 1.4; %MeTHF

theta1 = inc*pi/180;
theta2 = asin(n1/n2*sin(theta1));
b = b0/cos(theta2);


%Calculate extinction
e = (2*sqrt(2).*I)./((2.303*b*c).*I0);


%Propagation of error
dedI = (2*sqrt(2))./(2.303*b*c.*I0);
dedI0 = -(2*sqrt(2).*I)./(2.303*b*c.*I0.^2);
de = sqrt((dedI.^2.*dI.^2)+(dedI0.^2.*dI0.^2));

%Normalize to 1MV/cm
V = (Vrms*2.828)/2;
E2 = (V/b0).^2;


e = (e./E2)*1e6^2;
de = (de./E2)*1e6^2;






errorbar(w,e,de)




end