%%STOCHPLS(dataname,method,k,numpasses,ITERS,RESTART,eta_0,cap) runs stochastic PLS
%  The inputs are as follows
%
%  dataname   - one of these strings: 'VidTIMIT', 'VidTIMIT2', 'SIM'
%  method     - one of these strings: 'sgd', 'brand', 'batch', 'truth',
%                                       'msg', 'cmsg'
%  k          - positive integer - denotes the desired RANK
%  numpasses  - a positive integer - denotes number of passes over the data
%  ITERS      - an array of positve integers represening which random
%  splits - use only (1-1000)
%  RESTART    - a boolean flag - set to 1 if you want to restart from a
%  previous run - set to 0 if you want to discard previous runs
%  The output containing the population objective, empirical objective,
%  runtime, and the singular value decomposition (U,S,V) is written to a
%  MATFILE in ../PAGE/PROFILE/PLS/METHOD/DATANAME
%  The filename follows the format: 'method_progress[rank=1,pass=1,iter=1].mat'
%  eta_0      - initial learning rate
%  If possible combine random splits in the same run as each separate
%  run loads the data and increase I/O on the cluster
%  CITATION:
%  @inproceedings{arora2016stochastic,
%   title={Stochastic optimization for multiview representation learning using partial least squares},
%   author={Arora, Raman and Mianjy, Poorya and Marinov, Teodor},
%   booktitle={Proceedings of The 33rd International Conference on Machine Learning},
%   pages={1786--1794},
%   year={2016}
%  }

function score_tune = stochPLS(dataname,method,k,numpasses,ITERS,RESTART,eta_0,cap)

addpath(genpath('../'));

%% Default is cap = intmax
if(nargin<8)
    cap=intmax;
end

%% Default is eta_0 = 1
if(nargin<7)
    eta_0=1;
end

%% Default is over-write the previous runs
if(nargin<6)
    RESTART=0;
end

%% Default is just one run
if(nargin<5)
    ITERS=1;
end

%% Set PAGE directories
pagepath=sprintf('../PAGE/PROFILE/PLS/%s/%s/',method,dataname);
pageprefix=@(method,rank,pass,numiter)[pagepath,...
    sprintf('%s_progress[rank=%d,pass=%d,iter=%d].mat',method,rank,pass,numiter)];
if(~exist(pagepath,'dir'))     % Check if the PAGE directory is structured properly
    flag=createpath(pagepath); % If not create the desired directory structure
    if(~flag)                  % If the directory structure could not be created
        error('Could not create path for result files');% Display error message and quit
    end
end

%% Load data
if(strcmp(dataname,'SYN'))
    load('../DATA/SYN.mat','data');
    load('../DATA/permSYN.mat','perm');
elseif(strcmp(dataname,'XRMB'))
    load('../DATA/XRMB.mat','data');
    load('../DATA/permXRMB.mat','perm');
elseif(strcmp(dataname,'JW11'))
    load('../DATA/JW11.mat','data');
    load('../DATA/permJW11.mat','perm');
end

M1=size(data.view1.training,1); %#ok<NODEF>     dimension of the first view
M2=size(data.view2.training,1); % dimension of the second view
d=min(M1,M2); % PLS dim could not be larger than min of two views
N=size(data.view1.training,2); % number of training points

score_tune = 0;
for ITER=ITERS % each iteration uses a random shuffling over the data, as stored in perm
    %% Display the run
    fprintf('Starting run: (%s,%s,%s,%d,%d)\n',...
        dataname,'PLS',method,k,ITER);
    
    %% Get a random permutation
    iperm=perm(ITER,:); %#ok<NODEF> iperm is one specific permutation corresponding to the current ITER
    
    %% Random permutation of the data
    X=[data.view1.training data.view1.testing];
    data.view1.training=X(:,iperm(1,1:N)); % First half of data goes to traning
    data.view1.testing=X(:,iperm(1,N+1:2*N)); % Second half goes to testing
    clear('X');
    
    Y=[data.view2.training data.view2.testing];
    data.view2.training=Y(:,iperm(1,1:N));
    data.view2.testing=Y(:,iperm(1,N+1:2*N));
    clear('Y');
    
    %% Compute covariance matrix for the training and the held-out set
    
    CXY=(1/(size(data.view1.testing,2)-1))*... % CXY is TRUE covariance (TEST)
        (data.view1.testing-repmat(mean(data.view1.testing,2),1,size(data.view1.testing,2)))*...
        (data.view2.testing-repmat(mean(data.view2.testing,2),1,size(data.view1.testing,2)))';
    
    CXYhat=(1/(size(data.view1.training,2)-1))*... % CXY is EMPIRICAL covariance (TRAIN)
        (data.view1.training-repmat(mean(data.view1.training,2),1,size(data.view1.training,2)))*...
        (data.view2.training-repmat(mean(data.view2.training,2),1,size(data.view1.training,2)))';
    
    %% method=truth: Compute the maximum covariance that can be captured by any k dimensions
    if(strcmp(method,'truth'))
        [U,S,V]=svd(CXY); %#ok<ASGLU>
        fname=[pagepath,sprintf('truth[rank=%d,iter=%d].mat',k,ITER)];
        dS=diag(S); 
        dS=sort(dS,'descend'); %#ok<NASGU>
        save(fname,'dS');
        
        [U,S,V]=svd(CXYhat);
        fname=[pagepath,sprintf('truth_emp[rank=%d,iter=%d].mat',k,ITER)];
        dS=diag(S);
        dS=sort(dS,'descend'); %#ok<NASGU>
        save(fname,'dS');
        continue;
    end
    
    %% Sequence of iterations on which to compute objective
    [seq,L]=equilogseq(N,numpasses);
    % L is number of times objective will be evaluated
    % seq is the iterates on which objective will be evaluated
    offset=1;
    
    %% Initialize the basis for sgd and incremental PLS
    S=0;
    if(any(strcmp(method,'sgd')))
        U=orth(randn(M1,k));
        V=orth(randn(M2,k));
    elseif(strcmp(method,'inc'))
        U=zeros(M1,1);
        V=zeros(M2,1);
        S=0;
    elseif(any(strcmp(method,{'meg','msg'})))
        U=orth(random('normal',0,1,M1,k));
        V=orth(random('normal',0,1,M2,k));
        S=diag(ones(k,1)/k);
    end
    
    %% Initialize objective value and runtime
    objV=zeros(L(numpasses+1),1);
    objVe=zeros(L(numpasses+1),1);
    runtime=zeros(L(numpasses+1),1);
    
    %% Check if we can start from a previous run
    firstpass=1;
    initsamp=1;
    if(RESTART)
        for ipass=numpasses:-1:0
            if(exist(pageprefix(method,k,ipass,ITER),'file'))
                fprintf('Loading state from a previous run at pass number: %d\n',ipass);
                load(pageprefix(method,k,ipass,ITER),'runtime',...
                    'U','S','V','objV','objVe','seq');
                firstpass=ipass;
                if(isempty(find(objV>0,1,'last')))
                    offset=1;
                else
                    offset=find(objV>0,1,'last');
                end
                if(offset>=L(ipass+1))
                    firstpass=firstpass+1;
                    offset=1;
                end
                break;
            end
        end
    end
    
    %% Loop over the passes
    for ipass=firstpass:numpasses % ipass is the number of current pass over data
        
        fprintf('Starting pass number: %d\n',ipass);
        
        %% Output filename
        fname=pageprefix(method,k,ipass,ITER);
        
        %% Loop over data
        for iter=L(ipass)+offset:L(ipass+1)
            fprintf('Sequence number %d...\n',seq(iter));
            switch(method)
                case 'batch'
                    %% BATCH PLS
                    isamp=seq(iter);
                    if isamp==1 % otherwise computing Ctrain will be a devison by zero
                        continue;
                    end
                    if ipass>1
                        isamp=N;
                    end
                    tcounter=tic;
                    Ctrain=(1/(isamp-1))*((data.view1.training(:,1:isamp)-...
                        repmat(mean(data.view1.training(:,1:isamp),2),1,isamp))*...
                        (data.view2.training(:,1:isamp)-...
                        repmat(mean(data.view2.training(:,1:isamp),2),1,isamp))');
                    [U,S,V]=svds(Ctrain,k);
                    k2=min(k,size(U,2));
                    S=diag(S);
                    [~,idx]=sort(S,'descend');
                    idx=idx(1:k2);
                    S=S(idx);
                    U=U(:,idx);
                    V=V(:,idx);
                    runtime(iter)=toc(tcounter);
                    
                case 'sgd'
                    %% Stochastic gradient descent
                    for isamp=initsamp:seq(iter)
                        modisamp=1+mod(isamp-1,N);
                        etax=eta_0/(sqrt((ipass-1)*N+isamp));
                        etay=eta_0/(sqrt((ipass-1)*N+isamp));
                        tcounter=tic;
                        if size(U,2)>size(V,2)
                            zeropad=zeros(size(V,1),size(U,2)-size(V,2));
                            V=[V zeropad]; %#ok
                        elseif size(V,2)>size(U,2)
                            zeropad=zeros(size(U,1),size(V,2)-size(U,2));
                            U=[U zeropad]; %#ok
                        end
                        [U,V]=updateSGD_Cxy(U,V,data.view1.training(:,modisamp),...
                            data.view2.training(:,modisamp),etax,etay);
                        runtime(iter)=runtime(iter)+toc(tcounter);
                        if(mod(modisamp,1000)==0) %projection needed infrequently
                            U=Gram_Schmidt(U);
                            V=Gram_Schmidt(V);
                        end
                    end
                    
                case 'inc'
                    %% Incremental PLS method
                    for isamp=initsamp:seq(iter)
                        modisamp=1+mod(isamp-1,N);
                        tcounter=tic;
                        [U,S,V]=updateSVD_Cxy(U,S,V,data.view1.training(:,modisamp),...
                            data.view2.training(:,modisamp),k);
                        k2=min(k,size(U,2));
                        [~,idx]=sort(S,'descend');
                        idx=idx(1:k2);
                        S=S(idx);
                        U=U(:,idx);
                        V=V(:,idx);
                        runtime(iter)=runtime(iter)+toc(tcounter);
                    end
                    
                case 'msg'
                    %% Capped MSG-PLS
                    for isamp=initsamp:seq(iter)
                        modisamp=1+mod(isamp-1,N);
                        eta_t=eta_0/(sqrt((ipass-1)*N+isamp));
                        tcounter=tic;
                        [U,S,V]=updateSVD_Cxy(U,S,V,...
                            sqrt(eta_t)*data.view1.training(:,modisamp),...
                            sqrt(eta_t)*data.view2.training(:,modisamp),d);
                        S=msgproject(S,k);
                        [~,idx]=sort(S,'descend');
                        idx=idx(1:min(cap,length(idx)));
                        S=S(idx);
                        U=U(:,idx);
                        V=V(:,idx);
                        runtime(iter)=runtime(iter)+toc(tcounter);
                    end
                    
                    
                case 'meg'
                    %% MEG PLS
                    for isamp=initsamp:seq(iter)
                        modisamp=1+mod(isamp-1,N);
                        eta_t=eta_0/(sqrt((ipass-1)*N+isamp));
                        tcounter=tic;
                        [U,S,V]=meg_pls(k,U,S,V,data.view1.training(:,modisamp),...
                            data.view2.training(:,modisamp),eta_t);
                        [S,idx]=sort(S);
                        U=U(:,idx);
                        V=V(:,idx);
                        runtime(iter)=runtime(iter)+toc(tcounter);
                    end
                    
            end
            
            if(strcmp(method,'sgd'))
                U=Gram_Schmidt(U);
                V=Gram_Schmidt(V);
            end
            initsamp=seq(iter)+1;
            k2=min(size(U,2),size(V,2));
            keff=min(k,k2);
            if(any(strcmp(method,'msg'))) %JUST LIKE WARMUTH, NOT LIKE MSG PCA
                [U_msg,V_msg]=top_subspace(keff,U,S,V);
                objV(iter)=trace(abs(U_msg(:,1:keff)'*CXY*V_msg(:,1:keff)));
                objVe(iter)=trace(abs(U_msg(:,1:keff)'*CXYhat*V_msg(:,1:keff)));
            elseif(strcmp(method,'meg'))
                [U_meg,V_meg]=meg_pls_solution(keff,U,S,V);
                objV(iter)=trace(abs(U_meg(:,1:keff)'*CXY*V_meg(:,1:keff)));
                objVe(iter)=trace(abs(U_meg(:,1:keff)'*CXYhat*V_meg(:,1:keff)));
            elseif(strcmp(method,'sgd'))
                [U_sgd,V_sgd]=top_subspace(keff,U,ones(1,keff),V);
                objV(iter)=trace(abs(U_sgd(:,1:keff)'*CXY*V_sgd(:,1:keff)));
                objVe(iter)=trace(abs(U_sgd(:,1:keff)'*CXYhat*V_sgd(:,1:keff)));
            else
                [U_else,V_else]=top_subspace(keff,U,S,V);
                objV(iter)=trace(abs(U_else(:,1:keff)'*CXY*V_else(:,1:keff)));
                objVe(iter)=trace(abs(U_else(:,1:keff)'*CXYhat*V_else(:,1:keff)));
            end
            save(fname,'runtime','U','S','V','objV','objVe','seq');            
        end
        save(fname,'runtime','U','S','V','objV','objVe','seq');
    end
    
    %% tuning, currently using test
    if nargout > 0
        CXY_tun=(1/(size(data.view1.testing,2)-1))*... % CXY_tun is covariance on dev set (tune)
            (data.view1.tuning-repmat(mean(data.view1.tuning,2),1,size(data.view1.tuning,2)))*...
            (data.view2.tuning-repmat(mean(data.view2.tuning,2),1,size(data.view1.tuning,2)))';
        
        keff=min([k,size(U,2),size(V,2)]);
        if(strcmp(method,'sgd'))
            U=Gram_Schmidt(U);
            V=Gram_Schmidt(V);
            [U,V]=top_subspace(keff,U,ones(1,keff),V);
        elseif(strcmp(method,'meg'))
            [U,V]=meg_pls_solution(keff,U,S,V);
        else
            [U,V]=top_subspace(keff,U,S,V);
        end
        score_tune=score_tune+trace(abs(U(:,1:keff)'*CXY_tun*V(:,1:keff)));
    end
    
    
end

score_tune=score_tune/length(ITERS);

end
