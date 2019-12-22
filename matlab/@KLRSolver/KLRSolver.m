classdef KLRSolver
    %KLRSOLVER A class for the many different solving methods, assumes that
    % all classes are being used.
    
    properties
        % needed for any solve
        KA        % Kernel approximation, object of class KernApprox
        data      % Struct with Xtrain, Ytrain, Xtest, Ytest
        labels    % labels for the classes
        lambda    % regularization
        nn        % number of training points
        nt        % number of testing points
        mm        % size of small rank system
        dd        % dimension of data
        cc        % number of classes
        
        % specify types of solves, all have default vals
        ws = 0;   % number of warm start iterations
        inv_meth = 'dpcg'; % method of finding direction 
        tol_meth = 'tst'; % method for determining Newton convergence
        pr_flag  = true; % iteration print flag
        
        % iteration maxes, tol
        outer_its = 10;  % outer iteration max
        inner_its = 50;  % inner iteration max
        bktrk_its = 10;  % backtrack iter max
        tst_tol = 1e-2;  % Newton tolerance -- tst err
        trn_tol = 5e-3;  % Newton tolerance -- trn err
        grd_tol = 1e-3;  % Newton tolerance -- relgrad
        
        % Input, output
        theta0    % initial guess
        f0        % initial f val
        g0        % initial gradient
        gnorm0    % initial gradient norm
        theta     % current solve
        relgrad   % current rel_grad
        
        % Stats  
        times     % struct with tot,bck,inv,obj,cnv,con
        grd_errs  % rel grad vals across iters
        tst_errs  % tst errs across iters
        trn_errs  % trn errs across iters
        it_times  % iteration times
        iter = 0  % current iter
        in_steps = 0 % how many cg/inner steps we've taken
        bt_steps = 0 % how many backtrack steps
        
        % Test kernel -- if we want to store it
        Kt = 0;
    end
    
    methods
        % Constructor
        function obj = KLRSolver(ka,dat,lam,th0,options)
            tic;
            % the essentials
            obj.KA = ka;
            obj.data = dat;
            obj.lambda = lam;
            obj.nn = ka.nn;
            obj.nt = size(dat.Xtest,1);
            obj.mm = ka.mm;
            obj.cc = ka.cc;
            obj.dd = ka.dd;
            obj.labels = unique(dat.Ytrain);
            
            % check on data
            if norm(obj.KA.Xtrain - obj.data.Xtrain) ~= 0
                disp(['Warning!! data and KA do not match, using ', ...
                    'KA for training and data for testing']);
                obj.data.Xtrain = obj.KA.Xtrain;
                obj.data.Ytrain = obj.KA.Ytrain;
            end
            
            % Set options
            obj = obj.SetOptions(options);
            
            % initialize theta
            if isempty(th0)
                obj.theta0 = zeros(obj.mm,obj.cc);
            else
                obj.theta0 = th0;
            end
            obj.theta = obj.theta0;
            
            % Initial eval -- errs TODO?
            [obj.f0,obj.g0] = obj.KLR_Obj();
            obj.gnorm0 = norm(obj.g0,'fro');
            obj.relgrad = 1;
            
            % initialize times tot,bck,inv,obj,cnv,con
            obj.times.tot = 0;
            obj.times.bck = 0;
            obj.times.inv = 0;
            obj.times.obj = 0;
            obj.times.cnv = 0;
            obj.times.ws  = 0;
            obj.times.con = toc;
            obj.it_times(1) = 0;
            
            % initialize errs
						[terr,obj] = obj.TestError();
            obj.tst_errs(1) = terr;
            obj.grd_errs(1) = 1;
            obj.trn_errs(1) = obj.TrainError();
            
            if obj.pr_flag
                disp('-------------------------');
                disp(['Initialized KLRSolver object with ',obj.inv_meth]);
                disp(['Sigma: ',num2str(obj.KA.sigma)]);
                disp(['Lambda: ',num2str(obj.lambda)]);
                disp(['Tot rank: ',num2str(obj.KA.mm)]);
                disp(['Tol meth: ', obj.tol_meth]);
                disp(['Inv meth: ', obj.inv_meth]);
                disp('-------------------------');
                disp(['Initial fobj: ',num2str(obj.f0)]);
                disp(['Initial g0: ', num2str(obj.gnorm0)]);
                disp(['Initial tst err: ', num2str(obj.tst_errs(1))]);
                disp(['Initial trn err: ', num2str(obj.trn_errs(1))]);
                disp('-------------------------');
            end
            
        end
        
        % Returns options, if you want to use on another guy
        function options = GetOptions(obj)
            options.ws = obj.ws;
            options.tol_meth = obj.tol_meth;
            options.pr_flag = obj.pr_flag;
            options.outer_its = obj.outer_its;
            options.inner_its = obj.inner_its;
            options.bktrk_its = obj.bktrk_its;
            options.tst_tol = obj.tst_tol;
            options.trn_tol = obj.trn_tol;
            options.grd_tol = obj.grd_tol;        
        end
        
        % Set given options
        function obj = SetOptions(obj,options)
            % Set FindDir
            if isfield(options,'inv_meth')
                obj.inv_meth = options.inv_meth;
            end
            switch obj.inv_meth
                case 'cg'
                    obj.inner_its = 100; % adjust for cg
            end
            
            % Set up all the other options
            if isfield(options,'ws')
                obj.ws = options.ws;
            end
            if isfield(options,'tol_meth')
                obj.tol_meth = options.tol_meth;
            end
            if isfield(options,'pr_flag')
                obj.pr_flag = options.pr_flag;
            end
            if isfield(options,'outer_its')
                obj.outer_its = options.outer_its;
            end
            if isfield(options,'inner_its')
                obj.inner_its = options.inner_its;
            end
            if isfield(options,'bktrk_its')
                obj.bktrk_its = options.bktrk_its;
            end
            if isfield(options,'tst_tol')
                obj.tst_tol = options.tst_tol;
            end
            if isfield(options,'trn_tol')
                obj.trn_tol = options.trn_tol;
            end
            if isfield(options,'grd_tol')
                obj.grd_tol = options.grd_tol;
            end
        end
        
        % Methods needed, if they're in files leave them here 
        
        % Apply diag prec
        function out = KLR_ApplyDiagPrec(obj,d,vec)
            pvec = obj.ProjectNullSq(vec);
            pout = obj.KA.ApplyDiagPrec(d, pvec);
            out = obj.ProjectNullSq(pout);
        end
        
        % Create diag prec
        function d = MakeDiagPrec(obj,Pr)
            d = obj.KA.CompDiagPrec(Pr,obj.lambda);
        end
        
        
        % Make own files for these
        % [f,g,Pr] = KLR_Obj(obj) -> evaluates objective, grad, probs (X)
        % [Pr,G] = KLR_Pr(obj) -> evaluates probabilities (X)
        % [out] = KLR_SmHessMult(obj,Pr,vec) -> Small Hess multiply (X)
        % [mu,bst] = Backtrack(sdkfjsldk ) -> Backtracking (X) 
        % Solve(obj) -> returns solution in obj.theta  ( )
        % All the diff FindDirs DPCG (X), CG (X)
        
        
        % Project out of nullspace (vec is square)
        function [proj] = ProjectNullSq(obj,vec)
            % Q^T * theta
            th_sum = sum(vec,2)./obj.cc;
            
            % Q * theta not needed
            
            % I - Q Q^T theta
            proj = bsxfun(@minus, vec, th_sum);
            
        end
        
        % Project out of nullspace (vec is 1 long entry)
        function [proj] = ProjectNull(obj,vec)
            s = size(vec,2);
            
            if s ~= 1
                % Proj matrix
                Q = repmat(eye(obj.mm)./sqrt(obj.cc),obj.cc,1); % mem cost: c * m^2
                
                proj = vec - (Q * (Q' * vec));
            else
                % Do reshapes instead
                vsq = reshape(vec,obj.mm,obj.cc);
                psq = obj.ProjectNullSq(vsq);
                proj = reshape(psq,obj.mm*obj.cc,1);
            end
            
        end
        
        % Check tol_meth for convergence
        function [obj,flag,varargout] = CheckConv(obj,varargin)
            % process theta if we don't want to use obj's theta
            ll = length(varargin);
            th_flag = false;
            if ll == 1
                th = varargin{1};
                th_flag = true;
            end
            
            
            if th_flag
                % Compute train error
                tic;
                trn_err = obj.TrainError(th);
                trn_time = toc;
                
                % Compute test error
                tic;
                tst_err = obj.TestError(th);
                tst_time = toc;
                
                % Relgrad
                [~,g,~] = obj.KLR_Obj(th);
                grd_err = norm(g,'fro')/obj.gnorm0;
            else
                % Compute train error
                tic;
                trn_err = obj.TrainError();
                obj.trn_errs(obj.iter + obj.ws + 1) = trn_err;
                trn_time = toc;
                
                % Compute test error
                tic;
                [tst_err,obj] = obj.TestError();
                obj.tst_errs(obj.iter + obj.ws + 1) = tst_err;
								tst_time = toc;
                
                % Relgrad
                grd_err = obj.relgrad;
                obj.grd_errs(obj.iter + obj.ws + 1) = grd_err;
            end
            
            errs.grd = grd_err;
            errs.tst = tst_err;
            errs.trn = trn_err;
            
            if ~th_flag
                switch obj.tol_meth
                    case 'grd'
                        flag = grd_err < obj.grd_tol;
                        out_time = 0;
                    case 'tst'
                        t0 = obj.tst_errs(max(obj.iter + obj.ws,1));
                        flag = abs(tst_err - t0)/t0 < obj.tst_tol || ...
                            (tst_err - t0)/t0 > obj.tst_tol || ...
														tst_err == 0;
												flag = flag && ~(~obj.ws && obj.iter == 1); 

                        out_time = tst_time;
                    case 'trn'
                        t0 = obj.trn_errs(max(obj.iter + obj.ws,1));
                        flag = abs(trn_err - t0)/t0 < obj.trn_tol || ...
														trn_err == 0;
                        out_time = trn_time;
                end
            else
                flag = 0;
                out_time = tst_time;
            end
            
            varargout{1} = out_time;
            varargout{2} = errs;
        end
        
        % Evaluate train err at current theta (or alternate)
        function err = TrainError(obj,varargin)
            ll = length(varargin);
            th = obj.theta;
            if ll == 1
                th = varargin{1};
            end
            
            % K * th
            Kw = obj.KA.SM_MULT(th);
            
            % Max across classes
            [~,ycl] = max(Kw,[],2);
            
            % Evaluate
            err = sum(obj.labels(ycl) ~= obj.data.Ytrain)/obj.nn;
            
        end
            
        % Evaluate test err at current theta (or alternate) 
        function [err,obj] = TestError(obj,varargin)
            % Process extra input
            ll = length(varargin);
            th = obj.theta;
            if ll == 1
                th = varargin{1};
            end
            
            if ~obj.Kt
                % Need to calculate Kt
                obj.Kt = obj.KA.SKernel(obj.data.Xtest);
								%disp('kernel calc ...');
            end
            
            % Need to calculate Kt
            Kw = obj.KA.SM_TMULT(obj.Kt,th);
            
            % Max across classes
            [~,ycl] = max(Kw,[],2);
            
            % Evaluate
            err = sum(obj.labels(ycl) ~= obj.data.Ytest)/obj.nt;
            
        end
        
        % Evaluate error on other data
        function err = OthError(obj,Xoth,Yoth,varargin)
            % Process extra input
            ll = length(varargin);
            th = obj.theta;
            if ll == 1
                th = varargin{1};
            end
            
            % Need to calculate Kt
            Koth = obj.KA.SKernel(Xoth);
            
            % Need to calculate Kt
            Kw = obj.KA.SM_TMULT(Koth,th);
            
            % Max across classes
            [~,ycl] = max(Kw,[],2);
            
            % Evaluate
            noth = size(Xoth,1);
            err = sum(obj.labels(ycl) ~= Yoth)/noth;
            
        end
        
        % Finding direction function
        function [y,obj] = FindDir(obj,b,Pr)
            switch obj.inv_meth
                case 'dpcg'
                    [y,obj] = obj.FindDirDPCG(b,Pr);
                case 'lpcg'
                    [y,obj] = obj.FindDirLPCG(b,Pr);
                case 'cg'
                    [y,obj] = obj.FindDirCG(b,Pr);
                case 'en'
                    [y,obj] = obj.FindDirEN(b,Pr);
                otherwise
                    error(sprintf('Inversion method "%s" not supported',obj.inv_meth));
            end
        end
        
        % Truncate kernel approximation
        function [trun] = Truncate(obj,trank,Th,varargin)
            ll = length(varargin);
            if ll == 1
                options = varargin{1};
            else
                options = obj.GetOptions();
            end
            
            tTh = obj.KA.TruncateVec(trank,Th);
            tKA = obj.KA.Truncate(trank);
            trun = KLRSolver(tKA,obj.data,obj.lambda,tTh,options);       
        end

        % assemble table from outputs
        function [T] = AssembleTable(obj)
            iters = [-obj.ws:obj.iter]';
            cum_times = cumsum(obj.it_times(:));
            T = table(iters,obj.tst_errs(:),cum_times + obj.KA.decomp_time,'VariableNames', ...
                {'Iters','Errs','Time'});
        end
    end
    
end

