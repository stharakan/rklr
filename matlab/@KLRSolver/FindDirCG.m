function [y,klr] = FindDirCG( klr,rhs,Pr )
%FINDDIRDPCG Solves the problem H y = rhs for a given rhs, where H is the
%hessian defined by KLR_SmHessMult, along with a diagonal preconditioner
%defined by klr.Make

% Project out of nullspace
prhs = klr.ProjectNullSq(rhs);
tol_min = 1E-3;

% Define tolerance
tol = min(max(klr.relgrad,tol_min),0.5)/10;

% anon func
H = @(x) klr.KLR_SmHessMult(Pr,x);

% Call appropriate pcg
[ysq,~,cg_res,its] = pcg(H,reshape(prhs,klr.mm*klr.cc,1) ...
            ,tol,klr.inner_its);
        
% Reshape
y = reshape(-ysq,klr.mm,klr.cc);

if klr.pr_flag
    disp(['CG tol: ',num2str(tol)]);
    disp(['CG res: ',num2str(cg_res)]);
    disp(['CG its: ',num2str(its)]);
end
klr.in_steps = klr.in_steps + its;
end

