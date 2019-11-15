classdef EnsNyst < KernApprox
    %ONESHOT Subclass of KernApprox, builds the orthogonal 1 shot
    %approximation to a kernel matrix K.
    
    properties
        U 
        l
        Us
        smp_idx
        bb
        mm
        smp
        rnk
        decomp_time
    end
    
    methods
        % Constructor
        function obj = EnsNyst(Xtr,Ytr,samp,rank,sig,b)
            % call superconstructor
            obj@KernApprox(Xtr,Ytr,sig);
            
            % deal with rank not there
            ll = length(samp);
            if ll == 1
                if obj.yflag
                    idx = get_idx(Ytr,samp);
                else
                    idx = randperm(N);
                    idx = idx(1:samp);
                end
            else
                idx = samp;
                samp = length(samp);
                
            end
            if isempty(rank)
                rank = samp;
            end
            obj.mm = rank * b;
            obj.smp = samp;
            obj.smp_idx = idx;
            obj.rnk = rank;
            obj.bb = b;
            
            % Call actual one shot decomp
            tic;
            [uu,ll,uss,sidx] = obj.Decomp();
            obj.U = uu;
            obj.l = ll;
            obj.Us = uss;
            obj.smp_idx = sidx;
            obj.decomp_time = toc;
        end
        
        
        % Overload all KernApprox methods 
        
        % Big multiply U * L * U' * vec
        function out = BG_MULT(obj,vec)
            % U' multiply
            Utv = obj.U' * vec;
            
            % l multiply
            lUtv = bsxfun(@times,Utv, obj.l);
            
            % U multiply
            out = obj.U * lUtv;
            
            % Scale
            out = out./obj.bb;
        end
        
        % Inner multiply U * L^1/2 * vec
        function out = SM_MULT(obj,vec)
            % L^1/2 mult
            lprod = bsxfun(@times,vec,sqrt(obj.l));
            
            % U mult
            out = obj.U * lprod;
            
            % scale
            out = out./sqrt(obj.bb);
        end
        
        % Convert from big to small L^1/2 * U' * vec
        function out = BG_SM(obj,vec)
            % U' mult
            Utv = obj.U' * vec;
            
            % l^1/2 mult
            out = bsxfun(@times,Utv,sqrt(obj.l));
            
            % Scale
            out = out./sqrt(obj.bb);
            
        end
        
        % Approx inverse K^-1 * vec
        function out = BG_INV(obj,vec)
            % U' multiply
            Utv = obj.U' * vec;
            
            % l multiply
            lUtv = bsxfun(@rdivide,Utv, obj.l);
            
            % U multiply
            out = obj.U * lUtv;
            
            % Scale
            out = out./obj.bb;
        end
        
        % Approx small inverse l \ vec
        function out = SM_INV(obj,vec)
            % L^-1
            out = bsxfun(@rdivide,vec,obj.l);
        end
        
        % Map small to big U * L^-1/2
        function out = SM_BG(obj,vec)
            % Apply l^-1/2
            lv = bsxfun(@rdivide,vec,sqrt(obj.l));
            
            % U multiply
            out = obj.U * lv;
            
            % Scale
            out = out./sqrt(obj.bb);
        end
        
        % Small test multiply Kt * Us * L^1/2 * vec 
        function out = SM_TMULT(obj,Kt,vec)
            % init
            s = size(vec,2);
            m_idx = 1:obj.rnk;
            
            % L^1/2 mult
            lv = bsxfun(@times, vec, sqrt(obj.l));
            
            % Us mult, need loop
            Ulv = zeros(obj.rnk * obj.bb, s);
            for i = 1:obj.bb
                mc_idx = m_idx + (i-1) * obj.rnk;
                
                Ulv(mc_idx,:) = obj.Us(:,mc_idx) * lv(mc_idx,:);
            end
            
            % Kt mult
            out = Kt * Ulv;
            
            out = out./sqrt(obj.bb);
        end
        
        % Test multiply Kt * Us * L * U^T * vec TODO
        function out = BG_TMULT(obj,Kt,vec)
            % init
            s = size(vec,2);
            m_idx = 1:obj.rnk;
            
            % U' mult
            Uv = obj.U' * vec;
            
            % L^1/2 mult
            lv = bsxfun(@times, Uv, obj.l);
            
            % Us mult, need loop
            Ulv = zeros(obj.rnk * obj.bb, s);
            for i = 1:obj.bb
                mc_idx = m_idx + (i-1) * obj.rnk;
                
                Ulv(mc_idx,:) = obj.Us(:,mc_idx) * lv(mc_idx,:);
            end
            
            % Kt mult
            out = Kt * Ulv;
            
            
            out = out./obj.bb;
        end
        
        % Form kernel with samp_indices
        function Kt = SKernel(obj,Xt)
            % Get using kernelfun
            Kt = obj.kernelfun(Xt, obj.Xtrain(obj.smp_idx,:));
        end
        
        % Apply diag prec
        function [out] = ApplyDiagPrec(obj, d, vec)
            out = bsxfun(@rdivide,vec,d);
        end
        
        % Compute diag prec
        function [d] = CompDiagPrec(obj, Pr,lambda)
            [n,cl] = size(Pr);
            d = zeros(obj.mm*cl,1);
            
            %reg = spdiags(double(lambda).*l,0,k,k);
            reg = lambda.*single(ones(obj.mm,1));
            
            % form Ul
            Ul = bsxfun(@times,obj.U,sqrt(obj.l.')./obj.bb); %O(10m^2b^2)
            %Ul = Ul./sqrt(obj.bb);
            Uln = bsxfun(@times,Ul,Ul); %O(10m^2b^2)
            clear Ul
            
            % m idx
            m_idx = 1:obj.mm;
            
            Pr = bsxfun(@times,Pr,1-Pr); %O(20mC)
            %Q = zeros(size(Uln));
            
            for j = 1:cl %O(20 m^2 b^2 C)
                % idx
                mc_idx = m_idx + (j-1) * obj.mm;
                
                % Extract pi * 1 - pi
                dc = Pr(:,j);
                %Q = repmat(dc,1,size(Uln,2));
                
                % scale Ul to get W U L
                A = bsxfun(@times,Uln,dc); %O(10m^2b^2)
                %A = Q.*Uln;
                
                % dot multiply by U and sum to get diag(L U^T W UL)
                Bd = sum(A); %O(10m^2b^2)
                
                % load into d
                d(mc_idx) = reg + Bd'./n; 
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

