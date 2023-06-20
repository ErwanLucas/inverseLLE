function [deltaKerr,gainKerr] = computeKerrResponse(Aa, alpha)

if iscolumn(Aa)
    Aa = Aa.';
end

E3f = ifft(abs(Aa.').^2.*Aa.').';
E1f = ifft(Aa.').';
Kerr_resp = E3f./E1f;
if ~isrow(Aa)
    Kerr_resp = ifftshift(Kerr_resp, 2)*alpha;
else
    Kerr_resp = ifftshift(Kerr_resp)*alpha;
end
deltaKerr = real(Kerr_resp);
gainKerr = -imag(Kerr_resp);

end

