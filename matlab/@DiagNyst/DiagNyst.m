classdef DiagNyst < KernApprox
    %DIAGNYST Samples b batches from Xtrain, and forms rank rnk
    %approximations to each. The system is then solved as a diagonal system
    %of kernel matrices.
    
    properties
        U
        l
        Us
        smp_idx
        bb
        bb_size;
        mm
        smp
        rnk
        err
        decomp_time
        reorder_idx
    end
    
    methods
        % Constructor
        function obj = DiagNyst(Xtr,Ytr,samp,rank,sig,b,varargin)
            % First subsample Xtr and Ytr as needed
            n = size(Xtr,1);
            bsize = floor(n/b);
            nn = bsize * b;
            
            % call superconstructor
            if nn ~= n 
                bidx = get_idx(Ytr, nn);
                Xtr = Xtr(bidx,:);
                Ytr = Ytr(bidx);
            end
            % call superconstructor
            obj@KernApprox(Xtr,Ytr,sig,varargin{:});
           
            % deal with rank not there,init
            if isempty(rank)
                rank = samp;
            end
            obj.mm = rank * b;
            obj.smp = samp;
            obj.rnk = rank;
            obj.bb = b;
            obj.bb_size = bsize;
            obj.U = zeros(bsize,obj.mm);
            obj.l = zeros(obj.mm,1);
            obj.Us = zeros(rank,obj.mm);
            
            % Call actual one shot decomp
            tic;
            oth_idx = [];
            rem_idx = 1:nn;
            rnk_idx = 1:rank;
            bbs_idx = 1:bsize;
            for i = 1:b
                disp(['Computing decomposition ',num2str(i),' ...']);
                % Get smaller idx
                if b == 1
                    cur_idx = rem_idx;
                else
                    rem_idx = setdiff(rem_idx,oth_idx);
                    cur_idx = get_idx(Ytr(rem_idx),bsize);
                    oth_idx = [oth_idx,cur_idx];
                end
                
                % Subsample
                Xi = Xtr(cur_idx,:);
                Yi = Ytr(cur_idx);
                b_idx = bbs_idx + (i-1)*bsize; 
                obj.Xtrain(b_idx,:) = Xi;
                obj.Ytrain(b_idx) = Yi;
                
                % Decomposition
                [uu,ll,uss,sidx] = one_shot(Xi,rank,obj.kernelfun);
                
                % Load into obj
                m_idx = rnk_idx + (i-1)*rank;
                obj.U(:,m_idx) = uu;
                obj.l(m_idx) = ll;
                obj.Us(:,m_idx) = uss;
                obj.smp_idx(m_idx) = (i-1)*bsize + sidx;  
            end
            obj.reorder_idx = oth_idx;
            obj.decomp_time = toc;
        end
        
        % Big multiply U * L * U' * vec
        function out = BG_MULT(obj,vec)
            % Set up
            rnk_idx = 1:obj.rnk;
            bbs_idx = 1:obj.bb_size;
            out = zeros(obj.nn,size(vec,2));
            
            for i = 1:obj.bb
                m_idx = rnk_idx + (i-1)*obj.rnk;
                n_idx = bbs_idx + (i-1)*obj.bb_size;
                
                % Ut mult
                utprod = (obj.U(:,m_idx)' * vec(n_idx,:));
                
                % L mult
                lprod = bsxfun(@times,utprod,obj.l(m_idx));
                
                % U mult
                out(n_idx,:) = obj.U(:,m_idx) * lprod;
            end
            
            
        end
        
        % Inner multiply U * L^1/2 * vec
        function out = SM_MULT(obj,vec)
            % init
            out = zeros(obj.nn,size(vec,2));
            rnk_idx = 1:obj.rnk;
            bbs_idx = 1:obj.bb_size;
            
            % L^1/2 mult
            lprod = bsxfun(@times,vec,sqrt(obj.l));
            
            % U mult
            for i = 1:obj.bb
                m_idx = rnk_idx + (i-1)*obj.rnk;
                out(bbs_idx + (i-1)*obj.bb_size,:) = ...
                    obj.U(:,m_idx) * lprod(m_idx,:);
            end
        end
        
        % Convert from big to small L^1/2 * U' * vec
        function out = BG_SM(obj,vec)
            rnk_idx = 1:obj.rnk;
            bbs_idx = 1:obj.bb_size;
            dummy = zeros(obj.rnk,size(vec,2));
            
            % U mult
            for i = 1:obj.bb
                m_idx = rnk_idx + (i-1)*obj.rnk;
                n_idx = bbs_idx + (i-1)*obj.bb_size;
                
                dummy(m_idx,:) = obj.U(:,m_idx)' * vec(n_idx,:);
            end
            
            out = bsxfun(@times,dummy,sqrt(obj.l));
        end
        
        % Approx inverse K^-1 * vec
        function out = BG_INV(obj,vec)
            rnk_idx = 1:obj.rnk;
            bbs_idx = 1:obj.bb_size;
            out = zeros(obj.nn,size(vec,2));
            
            for i = 1:obj.bb
                m_idx = rnk_idx + (i-1)*obj.rnk;
                n_idx = bbs_idx + (i-1)*obj.bb_size;
                
                % Ut mult
                utprod = (obj.U(:,m_idx)' * vec(n_idx,:));
                
                % L mult
                lprod = bsxfun(@rdivide,utprod,obj.l(m_idx));
                
                % U mult
                out(n_idx,:) = obj.U(:,m_idx) * lprod;
            end
        end
        
        % Approx small inverse l \ vec
        function out = SM_INV(obj,vec)
            out = bsxfun(@rdivide,vec,obj.l);
        end
        
        % Approx inner inverse U * L^-1/2
        function out = SM_BG(obj,vec)
            % init
            out = zeros(obj.nn,size(vec,2));
            rnk_idx = 1:obj.rnk;
            bbs_idx = 1:obj.bb_size;
            
            % L^1/2 mult
            lprod = bsxfun(@rdivide,vec,sqrt(obj.l));
            
            % U mult
            for i = 1:obj.bb
                m_idx = rnk_idx + (i-1)*obj.rnk;
                out(bbs_idx + (i-1)*obj.bb_size,:) = ...
                    obj.U(:,m_idx) * lprod(m_idx,:);
            end
        end
        
        % Small test multiply Ft * vec
        function out = SM_TMULT(obj,Kt,vec)
            % Set up
            rnk_idx = 1:obj.rnk;
            %out = zeros(size(Kt,1)*obj.bb,size(vec,2));
            out = zeros(size(Kt,1),size(vec,2));
            bbt_idx = 1:size(Kt,1);
            
            % L^1/2 mult
            lprod = bsxfun(@times,vec,sqrt(obj.l));
            
            for i = 1:obj.bb
                m_idx = rnk_idx + (i-1)*obj.rnk;
                
                % Us mult
                usprod = obj.Us(:,m_idx) * lprod(m_idx,:);
                
                % Kt mult
                out = out + Kt(:,m_idx) * usprod;
                
                %t_idx = bbt_idx + (i-1)*size(Kt,1);
                %out(t_idx,:) = Kt(:,m_idx) * usprod;
            end
            %out = out./obj.bb;
        end
        
        % Test multiply Ft * F * vec
        function out = BG_TMULT(obj,Kt,vec)
            % Set up
            rnk_idx = 1:obj.rnk;
            bbs_idx = 1:obj.bb_size;
            %out = zeros(size(Kt,1)*obj.bb,size(vec,2));
            out = zeros(size(Kt,1),size(vec,2));
            bbt_idx = 1:size(Kt,1);
            
            for i = 1:obj.bb
                m_idx = rnk_idx + (i-1)*obj.rnk;
                n_idx = bbs_idx + (i-1)*obj.bb_size;
                
                % Ut mult
                utprod = (obj.U(:,m_idx)' * vec(n_idx,:));
                
                % L mult
                lprod = bsxfun(@times,utprod,obj.l(m_idx));
                
                % Us mult
                usprod = obj.Us(:,m_idx) * lprod;
                
                % Kt mult
                out = out + Kt(:,m_idx) * usprod;
                
                %t_idx = bbt_idx + (i-1)*size(Kt,1);
                %out(t_idx,:) = Kt(:,m_idx) * usprod;
            end
            %out = out./obj.bb;
        end
        
        % Form kernel with samp_indices
        function Kt = SKernel(obj,Xt)
            % Get using kernelfun
            Kt = obj.kernelfun(Xt, obj.Xtrain(obj.smp_idx,:));
            
            % Split up?
            if 0
                xte = size(Xt,1);
                bbt = xte/obj.bb;
                te_idx = 1:bbt;
                Kt = zeros(bbt,obj.mm);
                rnk_idx = 1:obj.rnk;
                
                for i = 1:obj.bb
                    m_idx = rnk_idx + (i-1)*obj.rnk;
                    t_idx = te_idx + (i-1)*bbt;
                    
                    Kt(:,m_idx) = obj.kernelfun(Xt(t_idx,:), ...
                        obj.Xtrain(obj.smp_idx(m_idx),:));
                end
            end
            
        end
        
        % Apply the diagonal preconditioner
        function [out] = ApplyDiagPrec(obj, d, vec)
            out = bsxfun(@rdivide,vec,d);
        end
        
        % Compute the diagonal preconditioner
        function [d] = CompDiagPrec(obj, Pr,lambda)
            [n,cl] = size(Pr);
            d = zeros(obj.mm*cl,1);
            
            %reg = spdiags(double(lambda).*l,0,k,k);
            reg = lambda.*single(ones(obj.mm,1));
            
            if 1 % old
            % form Ul
            Ul = bsxfun(@times,obj.U,sqrt(obj.l'))';
            Ul = reshape(Ul,obj.rnk,obj.nn);
            Uln = bsxfun(@times,Ul,Ul);
            
            % m idx
            rnk_idx = 1:obj.mm;
            
            % Probs
            Pr = bsxfun(@times,Pr,1-Pr);
            
            % Loop over cl
            for i = 1:cl
                % extract pi
                pi = Pr(:,i);
                %pi = reshape(Pr(:,i),obj.bsize,obj.bb);
                
                % Uln * diag(pi)
                A = reshape(bsxfun(@times,Uln,pi'),obj.mm,obj.bb_size);
                
                % sum up to get d
                d(rnk_idx + (i-1)*obj.mm) = sum(A,2) + reg;
            end
            
            else
                
                
            end
            
            
        end
        
        % Return a truncated rank-m version of full decomposition
        function [trunc] = Truncate(obj, trank)
            % First set it equal, then modify
            trunc = obj;
            trc = trank/obj.bb;
            trc_idx = repmat((1:trc)',1,obj.bb) + repmat(obj.rnk * (0:obj.bb-1),trc,1);
            trc_idx = reshape(trc_idx,trank,1);
            
            trunc.U = obj.U(:,trc_idx);
            trunc.l = obj.l(trc_idx);
            trunc.Us = obj.Us(1:trc,trc_idx);
            trunc.mm = trank;
            trunc.smp_idx = obj.smp_idx(trc_idx);
            trunc.smp = trc;
            trunc.rnk = trc;
        end
        
        % Return a truncated rank-m version of a vector
        function [tvec] = TruncateVec(obj, trank,vec)
            % Assume it is of small size
            trc = trank/obj.bb;
            trc_idx = repmat((1:trc)',1,obj.bb) + repmat(obj.rnk * (0:obj.bb-1),trc,1);
            trc_idx = reshape(trc_idx,trank,1);
            tvec = vec(trc_idx,:);
        end
        
        % Return a resized rank-m version of a truncated vector
        function [rvec] = ResizeVec(obj, vec)
            % Assume it is of small size
            [trank,cc] = size(vec);
            trc = trank/obj.bb;
            trc_idx = repmat((1:trc)',1,obj.bb) + repmat(obj.rnk *(0:obj.bb-1),trc,1);
            trc_idx = reshape(trc_idx,trank,1);
            
            rvec = zeros(obj.mm,cc);
            rvec(trc_idx,:) = vec;
        end
        
        
    end
    
end

