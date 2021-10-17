clc
clear all
close all
filename = './data/sytem20.mat';
Data_struct = load(filename,'-mat');
M_x=importdata('./epilepsy/sz5_ict.dat');
noOfNodes        = 76;
filtOrder        = 2;
noOfObservations = 4000

D=10
step=12;
order=filtOrder;
gamma_RF=1e-4;
gamma_NLTISO=1e-3;
lamda=70/50000;
m_X=M_x(:,1:noOfObservations);
m_X=M_x(:,1:noOfObservations);
% m_X(6,1:noOfObservations)=M_x(15,1:noOfObservations);
% m_X(12,1:noOfObservations)=M_x(16,1:noOfObservations);
% m_X(13,1:noOfObservations)=M_x(20,1:noOfObservations);
% m_X(14,1:noOfObservations)=M_x(19,1:noOfObservations);
 m_X=normalize(m_X,2)

%% TIRSO
addpath '/Users/rohanmoney/git_rohan/STInference/gsim/ProcessingBlocks/STFunctionEstimators/TirsoObjects/'


[n_row_mx,nTimeInstants]=size(m_X);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = 76;
tirsoObj.order     = 2; % we can try a higher order later
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
    m_predic(:,:)=tState_in.predictManyFromBuffer(1)';
    m_prediction(1:76,t)= m_predic(:,1);
    if t==50
        flag=1;
    end
      
end
%% RF_NLtiso
RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = 76;
RFObj.filtOrder = 2; % we can try a higher order later
RFObj.lambda    = 70/50000000;
RFObj.NoOfRF    =10;
RFObj.vsigma    =100;
RFObj.forgettingFactor=.98;
RFObj.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj.eta       =500;
%initialize
RFState_in = RFObj.initialize(1, m_X( :,1:RFObj.filtOrder)');
 for t = tirsoObj.order+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in = RFObj.update(RFState_in, mtemp);
    RF_m_predic(:,:)=RFState_in.predictManyFromBuffer(1)';
    RF_m_prediction(1:76,t)= RF_m_predic(:,1);
      
 end
sensor_to_plot=11;
plot(m_X(sensor_to_plot,24:noOfObservations),'LineWidth',2)
hold on

 plot(RF_m_prediction(sensor_to_plot,:),'LineWidth',2)
 hold on
  plot(m_prediction(sensor_to_plot,:),'LineWidth',2)
       