function [Z, AT, AW, W] = gnlse(t, V, A, w0, gamma, betas, loss, fr, Tr, dlength,RT,nplot)

% === propagatation of an optical field using the generalised NLSE
nt = length(t); dt = t(2)-t(1); 		% grid parameters		
alpha = log(10.^(loss/10));  			% attenuation coefficient
HOD = 0;
for ii = 1:length(betas)        		% Higher Order Dispersion (HOD)
  HOD = HOD + betas(ii)/factorial(ii).*V.^(ii);
end

L = 1i*HOD-(alpha)/2;           		% linear operator
if abs(w0) > eps               			% if w0>0 then include shock
    gamma = gamma/w0;    
    W = V + w0;                			% for shock W is true freq
else
    W = 1;                     			% set W to 1 when no shock
end
RW = nt*ifft(fftshift(RT.'));   		% frequency domain Raman
L = fftshift(L); W = fftshift(W); 		% shift to fft space


    function R = rhs(z, AW)
    AT = fft(AW.*exp(L*z)); % time domain field with variable change
    IT = abs(AT).^2; % time domain intensity
    if abs(fr) < eps % No Raman case
        M = ifft(AT.*IT); % Kerr + TPA response
    else
        RS = dt*fr*fft(ifft(IT).*RW); % Raman convolution
        M = ifft(AT.*((1-fr).*IT + RS + Tr)); % Full response (Kerr + TPA + Raman)
    end
    R = 1i*gamma*W.*M.*exp(-L*z); % full RHS of GNLSE
end

% === define function to print ODE integrator status
function status = report(z, y, flag) 
  status = 0;
  if isempty(flag)
    fprintf('%05.1f %% complete\n', z/dlength*100);
  end
end



% === setup and run the ODE integrator
Z = linspace(0, dlength, nplot);  		% select output z points
%Z = (flength/(nplot-1))*(0:nplot-1);



% === set error control options
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-12,'NormControl', 'on','OutputFcn', @report);
[Z, AW] = ode45(@rhs, Z, ifft(A), options); 	% run integrator


% === process output of integrator
AT = zeros(size(AW(1,:)));
for ii = 1:length(AW(:,1))
  AW(ii,:) = AW(ii,:).*exp(L.'*Z(ii)); 		% change variables
  AT(ii,:) = fft(AW(ii,:));           		% time domain output
  AW(ii,:) = fftshift(AW(ii,:))./dt;  		% scale
end
W = V + w0; 					% the absolute frequency grid
end

