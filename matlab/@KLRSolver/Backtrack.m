function [stepLength,varargout] = Backtrack(klr, curr_vec,curr_err,step_dir,residual)
%function v = backtack(vc,q,dv,residual, st, tol)
% vc: current, q: current residual inf norm, dv:direction, 
% residual: eval function  
%	 OPTIONAL st:steps, tol: stop tolerance
  
  if nargin == 1, selftest(klr); return; end;  

  %/*---*//*---*/
  %  DEFAULTS
  st = klr.bktrk_its; 
  tol = 0.99;
  ratio=2;
  stepLength = 1.;
  
  %/*---*//*---*/
  % MAIN ALGORITHM
  v  = curr_vec + stepLength*step_dir;
  %rp = max(abs(   feval( residual, v)  ))
  rp = feval(residual,v);

  rm=rp;
  for m=2:st
	 if rm < tol*curr_err, break;	 end;
	 if rm > rp,    
		v = curr_vec + ratio*stepLength*step_dir; m=m-1;
		break;   
	 end;

	 stepLength = stepLength/ratio;
	 v  = curr_vec + stepLength*step_dir;
	 rp=rm; 
     %rm = max(abs(   feval( residual, v)  ))	 
     rm = feval(residual,v);
  end
  
  if nargout > 1, varargout{1}=m; varargout{2}=v; end;
end
  
    
%/*
%**********************************************************************
%*/
function selftest(klr)
  vc =  1;
  dv = -1.8;
  q  =  1;
  v = klr.Backtrack(vc,q,dv,@(x) x.^2);

  v
end

