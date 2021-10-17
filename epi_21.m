clear all
data=edfread('/Users/rohanmoney/git_rohan/NonLinearTISO/chb01_03.edf');
sizure_start_time=seconds(2996)
sizure_end_time=seconds(3036)
sizure_timerange = timerange(sizure_start_time,sizure_end_time,'closed');
sizure_data= data(sizure_timerange,:);
m_X_sizure_data=cell2mat(sizure_data(:,:).Variables);
pre_sizure_start_time=seconds(2000);
pre_sizure_end_time=seconds(2995);
pre_sizure_timerange = timerange(pre_sizure_start_time,pre_sizure_end_time,'closed');
pre_sizure_data= data(pre_sizure_timerange,:);
m_X_pre_sizure_data=cell2mat(pre_sizure_data(:,:).Variables);
post_sizure_start_time=seconds(3037);
post_sizure_end_time=seconds(4000);
post_sizure_timerange = timerange(post_sizure_start_time,post_sizure_end_time,'closed');
post_sizure_data= data(post_sizure_timerange,:);
m_X_post_sizure_data=cell2mat(post_sizure_data(:,:).Variables);
save('X_chb01_03.mat','m_X_sizure_data','m_X_pre_sizure_data','m_X_post_sizure_data')
%    v_avgNorms(1:(noOfNodes +1):end)=0;  % avoiding selfloops to plot
       