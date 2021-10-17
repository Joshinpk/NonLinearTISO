clc
clear all
% close all
filename = './data/13hrs.mat';
Data_struct = load(filename,'-mat');
% M_x=Data_struct.m_X;
temp_m_x=Data_struct.tt_2;
M_x=temp_m_x.Variables';
noOfNodes        = 12;
filtOrder        = 12;
noOfObservations = 200
nTimeInstants=noOfObservations;
D=10
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
 m_X=normalize(m_X,2)
m_X(16,:)=0;

%% RF_NLtiso
RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = 24;
RFObj.filtOrder = 12; % we can try a higher order later
RFObj.lambda    = 700/50000;
RFObj.NoOfRF    =5;
RFObj.vsigma    =1*ones(24,1);

RFObj.forgettingFactor=.98 ;
RFObj.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj.eta       =900;
%initialize
RFState_in = RFObj.initialize(1, m_X( :,1:RFObj.filtOrder)');
 for t = RFObj.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in = RFObj.update(RFState_in, mtemp);
    RF_m_predic(:,:)=RFState_in.predictManyFromBuffer(12)';
    RF_m_prediction(1:24,t)= RF_m_predic(:,12);
%     if t==700
%       RFObj.eta       =500  
%     end
      
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


img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/system20.png');
 %image('CData',img,'XData',[-2 9],'YData',[-1 7]);
 image(img);
 hold on

Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5
v_avgNorms= mean(Psuedo_Adj(:,:,12), 3, 'omitNaN');
 v_avgNorms=v_avgNorms-diag(diag(v_avgNorms));
%         v_avgNorms(1:(noOfNodes +1):end)=0;  % avoiding selfloops to plot
        [~, order] = sort(v_avgNorms(:), 'descend');
        n_toPreserve = 2*noOfNodes; % displaying the edges in the order of the number of nodes 
        threshold = v_avgNorms(order(n_toPreserve));
        AdjacencyMatrix=(v_avgNorms>threshold).*v_avgNorms;
nodes_labels={'S1 (P)','S2 (P)','S3 (P)','S4 (P)','S5 (P)','S6 (P)','S7 (L)','S8 (L)','S9 (L)','S10(L)','S11(T)',...
            'S12 (L)','S13 (L)','S14 (L)','S15 (L)','S16 (L)','S17 (P)','S18 (P)','S19 (P)','S20 (T)','S21 (T)','S22 (T)','S23 (T)','S24 (T)'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
        x=[ 1000,1300,1750,1700,700,900,1000,900,1200,1900,800,800,700,800,900,1800,1000,800,1100,1800,1750,1100,700,750]
        y=[ 590,1000,1100,1050,170,250,720,720,1100,1100,600,290,360,330,330,1100,700,300,600,900,950,600,190,210]
       
%         imagesc(B)
                plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',16,'LineWidth',2)
                colorbar
set(gca,'XTick',[], 'YTick', [])
