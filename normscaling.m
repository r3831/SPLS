function view=normscaling(view)
%% Get dimensions
[~,N1]=size(view.training);
[~,N3]=size(view.testing);

%% Mean center everything
avg=mean(view.training,2);           % Mean of View1 training data
view.training=view.training-...      % Center View1 training data
  repmat(avg,1,N1);
view.testing=view.testing-...        % Center View testing data
  repmat(avg,1,N3);

%% Scale zero-mean data by standard deviation (N(0,1) scaling)

stdtrg1=std(view.training,[],2);     % STD DEV of  View1 train data
scmult1=1;
view.training=view.training./...     % Normalizing View1 train data
  (scmult1*repmat(stdtrg1,1,N1));
view.testing=view.testing./...       % Normalizing View1 test data
  (scmult1*repmat(stdtrg1,1,N3));

%% Only if data needed to be tuned
if isfield(view,'tuning')
    [~,N2]=size(view.tuning);
    view.tuning=view.tuning-...          % Center View1 tuning data
        repmat(avg,1,N2);
    view.tuning=view.tuning./...         % Normalizing View1 tune data
        (scmult1*repmat(stdtrg1,1,N2));
end

end

