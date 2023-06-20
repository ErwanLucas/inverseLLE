function y = shapingFunctions(x, selctStr, varargin)

switch selctStr
    case 'raisedCos'
        T = 1/varargin{1};
        beta = varargin{2};
        x1 = abs(x) <= (1-beta)/(2*T);
        x2 = abs(x) >= (1+beta)/(2*T);
        x3 = ~x2 & ~x1;
        if any(x3)
            x3 = double(x3) .* ( 1 + cos(pi*T/beta * (abs(x) - (1-beta)/(2*T)) ) )/2;
        end
        y = 1 * double(x1) + 0 * double(x2) + x3;
        
    case ''
        
    otherwise
        error('Unknown pulse type: %s', selctStr)
end


end

