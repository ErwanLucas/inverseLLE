function [Hom, R] = init_hom(D,X,root_selct, varargin)

if isempty(varargin)
    verbose = true;
else
    verbose = varargin{1};
end

if abs(X)==0
    Hom = 0;
    R = 0;
    return
end

R = roots([1, -2*D, (D^2+1), -X]);

real_soln = abs(imag(R)./real(R))<1e-2;
R = R(real_soln);
R = real(R);

if length(R)>2
    R = [min(R), max(R)];
    Hom = sqrt(X)/(1 + 1i*(D-R(root_selct))); % solution in the bistable region
    if verbose
        opt = {'lower', 'upper'};
        fprintf('Choosing %s branch\n', opt{root_selct})
    end
else
    Hom = sqrt(X)/(1 + 1i*(D-R(1)));
    if verbose
        disp('Single branch background')
    end
end
