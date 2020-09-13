function ns = N(~,yh,~)
%N Returns ODE nonlinear function
y  = ifft(yh);
ns = 2*1i*fft(y.^2 .* conj(y));
end