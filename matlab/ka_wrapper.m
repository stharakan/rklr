function [KA] = ka_wrapper(ka_type,varargin)
%KA_WRAPPER wraps all existing KernApprox types. 
% By passing in the ka_type and the proper args (as listed)
% in that KernApprox's constructor, you can use this 
% to construct any kernel approximation

switch ka_type
    case 'OneShot'
	    KA = OneShot(varargin{:});
    case 'DiagNyst'
	    KA = DiagNyst(varargin{:});
    case 'EnsNyst'
	    KA = EnsNyst(varargin{:});
    
    otherwise
        disp(' KernApprox type not found!!! - EXITING KA_WRAPPER EARLY');
        return;
end
