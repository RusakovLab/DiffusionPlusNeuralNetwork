function dy = NMDA(t,y,f,ft, TimeControl)
   f = interp1(ft, f, t); % Interpolate the data set (ft, f) at times t
 
% if (t < 0.1*TimeControl) || (t > 0.102*TimeControl)
%     Glu=0;
% end

 
	Rb	= 25*10^-3;    %(/uM /ms)	: binding        (5e-3   x5)
	Ru	= 47.5*10^-3;  %(/ms)	: unbinding      (9.5e-3  x5)
	Ro	= 80*10^-3;    %(/ms)	: opening        (16e-3   x5) -> open state
	Rc	= 65*10^-3;    %(/ms)	: closing        (13e-3   x5)
	Rd	= 125*10^-3;   %(/ms)	: desensitization(25e-3   x5)
	Rr	= 295*10^-3;   %(/ms)	: resensitization(59e-3   x5)
 
    rb1 = Rb * (1e3) * f;
dy = zeros(5,1);    % a column vector
NormParam = 1;
 
dy(1) =  NormParam*(Ru*y(2) - rb1*y(1));
dy(2) =  NormParam*(rb1*y(1) - (Ru+rb1)*y(2) + Ru*y(3));
dy(3) =  NormParam*(rb1*y(2) - (Ru+Rd+Ro)*y(3) + Rc*y(4) + Rr*y(5));
dy(4) =  NormParam*(Ro*y(3) - Rc*y(4));   % index 4 = OPEN (strong)
dy(5) =  NormParam*(Rd*y(3) - Rr*y(5));   % index 5 = DESENSITIZED (weak)
end