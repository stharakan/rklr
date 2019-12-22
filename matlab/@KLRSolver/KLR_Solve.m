function [ klr ] = KLR_Solve( klr )
%KLR_SOLVE Solves for a given klr object, outputting theta and all relevant
%statistics within the object itself

% Initial probabilities needed
tic;
[f,g,Pr] = klr.KLR_Obj();
obj_time = toc;
%g = klr.g0;
%f = klr.f0;
%size(klr.Kt)

if klr.ws
    klr = klr.KLR_WarmStart();
end


% Loop over outer iterations
for ii = 1:klr.outer_its
    klr.iter = klr.iter + 1;
    if klr.pr_flag
        disp('---------------------------');
        fprintf('Outer iteration %d \n',klr.iter);
    end
    
    % Find the direction
    tic;
    [dir,klr] = klr.FindDir(g,Pr);
    inv_time = toc;
    if klr.pr_flag
        fprintf('Inv time: %f\n',inv_time);
    end
    
    % Backtrack
    tic;
    [mu,bst] = klr.Backtrack(klr.theta,f,dir,@(t) klr.KLR_Obj(t));
    klr.bt_steps = klr.bt_steps + bst;
    bck_time =  toc;
    if klr.pr_flag
        fprintf('Bck time: %f | Bck steps: %d | Bck mu: %f\n',bck_time,bst,mu);
    end
    
    % Update theta
    tic;
    klr.theta = klr.theta + mu.*dir;
    [f,g,Pr] = klr.KLR_Obj();
    gn = norm(g,'fro');
    klr.relgrad = gn/(klr.gnorm0);
    obj_time = obj_time + toc;
    
    % Calculate stats
    [klr,cnv_flag,cnv_time,errs] = klr.CheckConv();
		tot_time = obj_time + cnv_time + bck_time + inv_time;
    klr.it_times(klr.iter + klr.ws + 1) = tot_time;
    if klr.pr_flag
        fprintf('Obj time: %f\n',obj_time);
        fprintf('Cnv time: %f\n',cnv_time);
        fprintf('It time: %f\n',tot_time);
        disp('        ------          ');
        fprintf('Tot time: %f\n',sum(klr.it_times));
        fprintf('Trn err: %f\n',errs.trn);
        fprintf('Tst err: %f\n',errs.tst);
        fprintf('Abs grd: %f\n',gn);
        fprintf('Rel grd: %f\n',klr.relgrad);
    end
    
    
    % Update times;
    klr.times.cnv = klr.times.cnv + cnv_time;
    klr.times.obj = klr.times.obj + obj_time;
    klr.times.bck = klr.times.bck + bck_time;
    klr.times.inv = klr.times.inv + inv_time;
    klr.times.tot = klr.times.tot + tot_time;
    obj_time = 0;
    
    % Break if necessary
    if cnv_flag
        break;
    end
    
end
if klr.pr_flag
    disp('---------------------------');
    disp('   Outer iters complete!   ');
    fprintf('Tot steps: %d\n',klr.iter);
    fprintf('Inv time: %f | In steps %d\n',klr.times.inv,klr.in_steps);
    fprintf('Bck time: %f | Bck steps: %d\n',klr.times.bck,klr.bt_steps);
    fprintf('Obj time: %f\n',klr.times.obj);
    fprintf('Cnv time: %f\n',klr.times.cnv);
    fprintf('Tot time: %f\n',klr.times.tot);
    disp('---------------------------');
    disp('    Final conv metrics     ');
    fprintf('Trn err: %f\n',errs.trn);
    fprintf('Tst err: %f\n',errs.tst);
    fprintf('Abs grd: %f\n',gn);
    fprintf('Rel grd: %f\n',klr.relgrad);
    disp('---------------------------');
end

