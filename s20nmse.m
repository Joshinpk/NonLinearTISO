clear all
% close all
filename = './data/s20.mat';
Data_struct = load(filename,'-mat');
% M_x=Data_struct.m_X;
temp_m_x=Data_struct.tt_2;
M_x=temp_m_x.Variables';
noOfNodes        = 12;
filtOrder        =12;
noOfObservations = 300
nTimeInstants=noOfObservations;
D=5
step=12;
order=filtOrder;
gamma_RF=1e-4;
gamma_NLTISO=1e-3;
lamda=70/50000;
m_X=M_x(1:24,1:noOfObservations);
m_X=M_x(1:24,1:noOfObservations);
% m_X(6,1:noOfObservations)=M_x(15,1:noOfObservations);
% m_X(12,1:noOfObservations)=M_x(16,1:noOfObservations);
% m_X(13,1:noOfObservations)=M_x(20,1:noOfObservations);
% m_X(14,1:noOfObservations)=M_x(19,1:noOfObservations);
%    m_X=normalize(m_X,2)

%% TIRSO
addpath '/Users/rohanmoney/git_rohan/STInference/gsim/ProcessingBlocks/STFunctionEstimators/TirsoObjects/'


[n_row_mx,nTimeInstants]=size(m_X);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = 24;
tirsoObj.order     = 12; % we can try a higher order later
tirsoObj.regPar    = 1e-2;
tirsoObj.b_shrinkSelfLoops  = 1; % Bolstad
tirsoObj.forgettingFactor   = 0.98;
tirsoObj.h_stepsize         = @(ts)1/eigs(ts.m_Phi,1);
% initialize
tState_in = tirsoObj.initialize(0, m_X( :,1:tirsoObj.order)');
e_tirso=zeros(24,noOfObservations);
for t = tirsoObj.order+1:nTimeInstants
    tic
    mtemp= m_X(:, t);
    tState_in = tirsoObj.update(tState_in, mtemp);
    m_predic(:,:)=tState_in.predictManyFromBuffer(12)';
    m_prediction(1:24,t)= m_predic(:,12);
    if t==50
        flag=1;
    end
      
end
%% RF_NLtiso
RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = 24;
RFObj.filtOrder = 12; % we can try a higher order later
RFObj.lambda    = 1900/900;
RFObj.NoOfRF    =5;
RFObj.vsigma(1:6)    =1*ones(6,1);
RFObj.vsigma(7:16)    =1*ones(10,1);
RFObj.vsigma(17:19)    =1*ones(3,1);
RFObj.vsigma(20:24)    =1*ones(5,1);
RFObj.forgettingFactor=.98 ;
RFObj.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj.eta       =4000;
%initialize
RFState_in = RFObj.initialize(1, m_X( :,1:RFObj.filtOrder)');
 for t = RFObj.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in = RFObj.update(RFState_in, mtemp);
    RF_m_predic(:,:)=RFState_in.predictManyFromBuffer(12)';
    RF_m_prediction(1:24,t)= RF_m_predic(:,12);
%     if t==700
%       RFObj.eta       =500  
%     e
 end
%  sensor_to_plot=15;
% %  plot(m_X(sensor_to_plot,12:noOfObservations),'LineWidth',2)
% %  hold on
% % % 
% %   plot(RF_m_prediction(sensor_to_plot,:),'LineWidth',2)
% %  hold on
% %   plot(m_prediction(sensor_to_plot,:),'LineWidth',2)
% % %   hold on
% %   plot((RF_m_prediction(sensor_to_plot,1:end-12)-m_X(sensor_to_plot,13:end)).^2)
%  error_temp=(RF_m_prediction(sensor_to_plot,1:end-12)-m_X(sensor_to_plot,13:end)).^2;
%  for ii=1:length(error_temp)
%      nmse(ii)=sum(error_temp(1:ii))/sum((m_X(sensor_to_plot,13:ii+12)).^2);
%  end
%  hold on
%  plot(nmse)
%   error_temp=(m_prediction(sensor_to_plot,1:end-12)-m_X(sensor_to_plot,13:end)).^2;
%  for ii=1:length(error_temp)
%      nmse(ii)=sum(error_temp(1:ii))/sum((m_X(sensor_to_plot,13:ii+12)).^2);
%  end
%  hold on
%  plot(nmse)


 error_temp=sum((RF_m_prediction(:,1:end-12)-m_X(:,13:end)).^2,1);
 for ii=1:length(error_temp)
     nmse(ii)=sum(error_temp(1:ii))/sum(sum((m_X(:,13:ii+12)).^2,1));
 end
 hold on
 plot(nmse)
 error_temp=sum((m_prediction(:,1:end-12)-m_X(:,13:end)).^2,1);
 for ii=1:length(error_temp)
     nmse(ii)=sum(error_temp(1:ii))/sum(sum((m_X(:,13:ii+2)).^2,1));
 end
 hold on
 plot(nmse)