function [y,klr] = FindDirEN( klr,rhs,Pr )
%FINDDIRDPCG Solves the problem H y = rhs for a given rhs, where H is the
%hessian defined by KLR_SmHessMult, along with a diagonal preconditioner
%defined by klr.Make

% Project out of nullspace
prhs = klr.ProjectNullSq(rhs);

% Define tolerance
tol = 1E-4; 

% make anon funcs
H = @(x) klr.KLR_SmHessMult(Pr,x);
d = klr.MakeDiagPrec(Pr);
M = @(x) klr.KLR_ApplyDiagPrec(d,x);

% Call appropriate pcg
[ysq,~,pcg_res,its] = pcg(H,reshape(prhs,klr.mm*klr.cc,1) ...
            ,tol,klr.inner_its,M);
        
% Reshape
y = reshape(-ysq,klr.mm,klr.cc);

if klr.pr_flag
    disp(['PCG tol: ',num2str(tol)]);
    disp(['PCG res: ',num2str(pcg_res)]);
    disp(['PCG its: ',num2str(its)]);
end
klr.in_steps = klr.in_steps + its;
 
end

