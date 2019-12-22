classdef KernApprox
    %KernApprox Parent class for all kernel deocmpositions
    %   Never called explicitly, simply call the appropriate kernel
    %   decomposition, and it will call this class. All the methods listed
    %   are overwritten in the particular subclass
    
    properties
        Xtrain
        Ytrain
        cc
        nn
        dd
        kernelfun
        sigma
        yflag
        kerr
    end
    
    methods
        % Constructors
        function obj = KernApprox(Xtr,Ytr,sig,varargin)
            % do nothing with zero arg case
            if nargin > 0
                % initialize
                obj.Xtrain = Xtr;
                obj.Ytrain = Ytr;
                [N,d] = size(Xtr);
                obj.nn = N;
                obj.dd = d;
                
                % deal with no Ytrain case
                yfl = ~isempty(Ytr);
                obj.yflag = yfl;
               
                if obj.yflag
                    obj.cc = length(unique(Ytr));
                end
                
                % deal with sigma
                obj.sigma = sig;
                
                
                %obj.kernelfun = @(x,y) gaussiankernel(x,y,sig);
                if isempty(varargin)
                    obj.kernelfun = @(x,y) gaussiankernel(x,y,sig);
                else
                    obj.kernelfun = varargin{1};
                end
            end
            
        end
        
        % Big multiply K * vec
        function out = BG_MULT(obj,vec)
            error('Function BG_MULT not overwritten');
            out = zeros(obj.nn,1);
        end
        
        % Inner multiply F * vec
        function out = SM_MULT(obj,vec)
            error('Function SM_MULT not overwritten');
            out = zeros(obj.nn,1);
        end
        
        % Convert from big to small F' * vec
        function out = BG_SM(obj,vec)
            error('Function BG_SM not overwritten');
            out = zeros(obj.samp,1);
        end
        
        % Approx inverse K^-1 * vec
        function out = BG_INV(obj,vec)
            error('Function BG_INV not overwritten');
            out = zeros(obj.nn,1);
        end
        
        % Approx small inverse l \ vec 
        function out = SM_INV(obj,vec)
            error('Function SM_INV not overwritten');
            out = zeros(obj.mm,1);
        end
         
        % Approx inner inverse F \ vec
        function out = SM_BG(obj,vec)
            error('Function SM_BG not overwritten');
            out = zeros(obj.nn,1);
        end
        
        % Small test multiply Ft * vec
        function out = SM_TMULT(obj,Kt,vec)
            error('Function SM_TMULT not overwritten');
            nt = size(Kt,1);
            out = zeros(nt,1);
        end
        
        % Test multiply Ft * F * vec
        function out = BG_TMULT(obj,Kt,vec)
            error('Function BG_TMULT not overwritten');
            nt = size(Kt,1);
            out = zeros(nt,1);
        end
        
        % Form kernel with samp_indices
        function Kt = SKernel(obj,Xt)
            error('Function SKernel not overwritten');
            nt = size(Xt,1);
            Kt = zeros(nt,1);
        end
        
        % Apply the diagonal preconditioner
        function [out] = ApplyDiagPrec(obj, d, vec)
            out = 0;
            error('Function ApplyDiagPrec not overwritten');
        end
        
        % Compute the diagonal preconditioner
        function [d] = CompDiagPrec(obj, Pr,lambda)
            d = 0;
            error('Function CompDiagPrec not overwritten');
        end
        
        % Return a truncated rank-m version of full decomposition
        function [trunc] = Truncate(obj, trank)
            trunc = 0;
            error('Function Truncate not overwritten');
        end
        
        % Return a truncated rank-m version of a vector
        function [tvec] = TruncateVec(obj, trank,vec)
            % Assume it is of small size
            tvec = 0;
            error('Function TruncateVec not overwritten');
        end
        
        % Return a resized rank-m version of a truncated vector
        function [rvec] = ResizeVec(obj, m,vec)
            % Assume it is of small size
            rvec = 0;
            error('Function ResizeVec not overwritten');
        end
        
        % Compute the error of the approximation
        function [rel_err,W] = matvec_errors(obj,runs)
            norm_sample_size = 1000;
            if nargin == 1
                runs = 50;
            end
            
            [x1,x2] = size(runs);
            if x1 == 1 && x2 == 1 
                % runs is a number
                W = normrnd(0,1,[obj.nn,runs]);
                W2 = W.*W;
                normw = sqrt(abs(sum(W2,1)));
                normw = repmat(1./normw,obj.nn,1);
                W = W.*normw;
            else
                % evaluate at given W
                W = runs;
                runs = size(W,2);
            end
            
            % get test idx
            norm_sample_size = min(norm_sample_size, obj.nn);
            tst_idx = floor(1:(obj.nn/norm_sample_size):obj.nn);
            
            % real kernel
            truK = obj.kernelfun(obj.Xtrain(tst_idx,:),obj.Xtrain);
            
            % multiplies
            truKw = truK * W;
            estKw = obj.BG_MULT(W);
            estKw = estKw(tst_idx,:);
            
            % difference norms for each vector
            diff = truKw - estKw;
            diff = diff.^2;
            abserrs = sqrt(abs(sum(diff,1)));
            truKw = truKw.^2;
            relerrs = sqrt(abs(sum(truKw,1)));
            
            % get final errors
            abs_err = sum(abserrs)/runs;
            rel_err = sum(abserrs./relerrs)/runs;
        end

        function [Pr, G, varargout] = klr_probs(obj,theta)
            theta_sm = obj.BG_SM(theta);

            if (nargout > 2)
                [Pr, G , varargout] = obj.klr_reduced_probs(theta_sm);
            else
                [Pr,G] = obj.klr_reduced_probs(theta_sm);
            end
        end

        function [Pr, G,varargout] = klr_reduced_probs(obj, theta)
            [rank, cc] = size(theta);

            norm_flag = false;
            if(nargout > 2)
                norm_flag = true;
                norms = zeros(cc - 1,1,'single');
            end

            % rewrite multiply
            F = obj.SM_MULT(theta);

            maxF = max(F,[],2);
            G = exp(bsxfun(@minus,F,maxF));

            if norm_flag
                norms = sqrt(sum(F.*F));
            end

            Pr = bsxfun(@rdivide,G,sum(G,2));


            pr_sum = sum(Pr,2);
            Pr = Pr./repmat(pr_sum,1,cc);

            if(norm_flag)
                varargout{1} = norms;
            end
        end

    end
    
end

