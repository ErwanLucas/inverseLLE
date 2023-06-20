function [Psi]  = LLE_ode_para(psi_0,FT, Ft, zeta, D)

%[~,i_peak] = max(abs(psi_0));

F = FT*Ft;
dd = D-1i*zeta;

t = linspace(0,.1,3);

figure
options = odeset('NormControl','off', 'OutputFcn', @(t,y,flag) report(t,y,flag));

[T,Psi] = ode45(@(t,y) rhs(t,y,dd,F), t, psi_0.', options);


end

%% Other sub functions
function [dy] = rhs(t,psi,D,S)
if ~isrow(psi)
    psi = psi.';
end
Disp = fft( D.*ifft(psi) );
Nlin = 1i * conj(psi) .* psi.^2;
dy = Disp + Nlin + S;
dy = dy.';
end

function status = report(t,psi, flag)
status = 0;
if isempty(flag)
    plot(log(abs(psi(:,end))));
    drawnow()
    clc
    fprintf('Time %f\n', t(end));
end

end