function varargout = main(varargin)

varargout{1} = nargout;
varargout{2} = nargin;
varargout{3} = 'You did it';
if isempty(varargout{2})
    disp('2nd output omitted')
end