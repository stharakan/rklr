classdef OneShot < KernApprox
    %ONESHOT Subclass of KernApprox, builds the orthogonal 1 shot
    %approximation to a kernel matrix K.
    
    properties
        U 
        l
        Us
        smp_idx
        mm
        smp
        rnk
        err
        decomp_time
    end
    
    methods
        % Constructor
        function obj = OneShot(Xtr,Ytr,samp,rank,sig,varargin)
            % call superconstructor
            obj@KernApprox(Xtr,Ytr,sig,varargin{:});
            
            % deal with samp being given as smp_idx, rank not there
            ll = length(samp);
            if ll == 1
                idx = randperm(obj.nn);
                idx = idx(1:samp);
            else
                idx = samp;
                samp = length(samp);
                
            end
            if isempty(rank)
                rank = samp;
            end
            obj.mm = rank;
            obj.smp = samp;
            obj.smp_idx = idx;
            obj.rnk = rank;
            
            % Call actual one shot decomp
            tic;
            [UU,ll,Uss] = one_shot(obj.Xtrain, ...
                obj.smp_idx,obj.kernelfun);
            obj.decomp_time = toc;
            obj.U = UU(:,1:rank);
            obj.l = ll(1:rank);
            obj.Us = Uss(:,1:rank);
            
        end
        
        
        % Overload all KernApprox methods
        
        % Big multiply U * L * U' * vec
        function out = BG_MULT(obj,vec)
            s = size(vec,2);
            if s == 1
                out = obj.U * (obj.l .* (obj.U' * vec) );
            else
                out = obj.U * (diag(obj.l) * (obj.U' * vec) );
            end
        end
        
        % Inner multiply U * L^1/2 * vec
        function out = SM_MULT(obj,vec)
            s = size(vec,2);
            if s == 1
                out = obj.U * (sqrt(obj.l) .* vec);
            else
                out = obj.U * (diag(sqrt(obj.l)) * vec);
            end
        end
        
        % Convert from big to small L^1/2 * U' * vec
        function out = BG_SM(obj,vec)
            s = size(vec,2);
            if s == 1
                out = sqrt(obj.l) .* ( obj.U' * vec);
            else
                out = diag(sqrt(obj.l)) * (obj.U' * vec);
            end
        end
        
        % Approx inverse K^-1 * vec
        function out = BG_INV(obj,vec)
            s = size(vec,2);
            if s == 1
                out = obj.U * ((obj.U' * vec) ./ obj.l );
            else
                out = obj.U * (diag(1./obj.l) * (obj.U' * vec) );
            end
        end
        
        % Approx small inverse l \ vec 
        function out = SM_INV(obj,vec)
            out = bsxfun(@rdivide,vec,obj.l);
        end
        
        % Map small to big U * L^-1/2 TODO
        function out = SM_BG(obj,vec)
            s = size(vec,2);
            if s == 1
                out = obj.U * (vec ./ sqrt(obj.l));
            else
                out = obj.U * (diag(1./sqrt(obj.l)) * vec);
            end
        end
        
        % Small test multiply Ft * vec
        function out = SM_TMULT(obj,Kt,vec)
            s = size(vec,2);
            if s == 1
                out = Kt * (obj.Us * (sqrt(obj.l) .* vec) );
            else
                out = Kt * (obj.Us * (diag(sqrt(obj.l)) * vec) );
            end
        end
        
        % Test multiply Ft * F * vec
        function out = BG_TMULT(obj,Kt,vec)
            s = size(vec,2);
            if s == 1
                out = Kt * (obj.Us * (obj.l .* (obj.U' * vec) ) );
            else
                out = Kt * (obj.Us * (diag(obj.l) * (obj.U' * vec) ) );
            end
        end
        
        % Form kernel with samp_indices
        function Kt = SKernel(obj,Xt)
            Kt = obj.kernelfun(Xt,obj.Xtrain(obj.smp_idx,:));
        end
        
        % Apply diag prec
        function [out] = ApplyDiagPrec(obj, d, vec)
            out = bsxfun(@rdivide,vec,d);
        end
        
        % Return a truncated rank-m version of full decomposition
        function [trunc] = Truncate(obj, trank,th)
            % First set it equal, then modify
            trunc = obj;
            
            trunc.U = obj.U(:,1:trank);
            trunc.l = obj.l(1:trank);
            trunc.Us = obj.Us(1:trank,1:trank);
            trunc.mm = trank;
            trunc.smp_idx = obj.smp_idx(1:trank);
            trunc.smp = trank;
            trunc.rnk = trank;
        end
        
        % Return a truncated rank-m version of full decomposition
        function [tvec] = TruncateVec(obj, trank,vec)
            % Assume it is of small size
            tvec = vec(1:trank,:);
        end
        
        % Return a resized rank-m version of a truncated vector
        function [rvec] = ResizeVec(obj,vec)
            % Assume it is of small size
            [sm,cc] = size(vec);
            dif = obj.mm - sm;
            rvec = [vec;zeros(dif,cc)];
        end
        
        % Compute diag prec
        function [d] = CompDiagPrec(obj, Pr,lambda)
            
            [n,cl] = size(Pr);
            d = zeros(cl*obj.rnk,1);
            
            %reg = spdiags(double(lambda).*l,0,k,k);
            reg = lambda.*single(ones(obj.rnk,1));
            %Ul = obj.U* diag(sqrt(obj.l));
            Ul = bsxfun(@times,obj.U, sqrt(obj.l)');
            Ul = Ul .* Ul;

            
            for i = 1:cl
                idx = ((i-1)*obj.rnk +1):(i*obj.rnk);
                pi = Pr(:,i);
                
                % do diagonal block first
                dc = pi.*(1-pi);
                
                % scale Ul to get W U L
                %A = repmat(dc,1,obj.rnk) .* Ul;
                A = bsxfun(@times,dc,Ul);
                
                % dot multiply by U and sum to get diag(L U^T W UL)
                Bd = sum(A);
                
                d(idx) = reg + Bd'./n;
                
            end
        end
        
    end
    
end

