function spectr = spectrumF(X)
spectr = ifftshift(abs(ifft(X)).^2,1);
end