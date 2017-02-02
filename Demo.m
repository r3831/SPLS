% DEMO runs several PLS algorithms, tunes the parameters on a held out set,
% and outputs the accuracy over a test set.
% methods = {truth, batch, incremental, msg, meg, sgd}
% datasets = {SYN, XRMB}
% eta is a structure which keeps the track of best learning rate over the
% held out set.
%     N     should be set to the size of the dataset being used. For example, for
%           our 'SYN' dataset, we have N=1000 data points.
%     k     Desired rank
%     rst   1: start from a previous session, 0: reset
%     run   number of runs over the dataset
%     itr   number of experiments, to be averaged over

%% initialize experimental variables
k=4; itr=2; dataset='SYN'; N=1000; rst=0; runs=5;
methods={'batch','inc','sgd','meg','msg'};

%% Running {truth, batch, incremental}. No parameters to be tunned.
stochPLS(dataset,'batch',k,runs,1:itr,rst);
stochPLS(dataset,'truth',k,runs,1:itr,rst);
stochPLS(dataset,'inc',k,runs,1:itr,rst);

%% Running {msg, meg, sgd}. The initial learning rate eta to be tunned.
eta.msg.value=.001;
eta.msg.best=-Inf;
eta.meg.value=.001;
eta.meg.best=-Inf;
eta.sgd.value=.001;
eta.sgd.best=-Inf;

for eta_0=[.001 .01 .1 1 10]
    %% msg
    current_msg=stochPLS(dataset,'msg',k,1,1:itr,rst,eta_0);
    if current_msg > eta.msg.best
        eta.msg.best=current_msg;
        eta.msg.value=eta_0;
    end
        
    %% meg
    current_meg=stochPLS(dataset,'meg',k,1,1:itr,rst,eta_0);
    if current_meg > eta.meg.best
        eta.meg.best=current_meg;
        eta.meg.value=eta_0;
    end
    
    %% sgd
    current_sgd=stochPLS(dataset,'sgd',k,1,1:itr,rst,eta_0);
    if current_sgd > eta.sgd.best
        eta.sgd.best=current_sgd;
        eta.sgd.value=eta_0;
    end
end

%% Testing the methods using the best learning rate eta
stochPLS(dataset,'msg',k,runs,1:itr,rst,eta.msg.value);
stochPLS(dataset,'meg',k,runs,1:itr,rst,eta.meg.value);
stochPLS(dataset,'sgd',k,runs,1:itr,rst,eta.sgd.value);

%% Plotting the PLS objective on the test set
plotobjV(dataset,methods,N,k,runs,1:itr,'avg',0)
plotobjV(dataset,methods,N,k,runs,1:itr,'avg',1)

