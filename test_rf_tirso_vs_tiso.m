clc
clear all
%close all
% filename = './data/sytem20.mat';
% Data_struct = load(filename,'-mat');
% M_x=Data_struct.m_X;
noOfNodes        = 5;
filtOrder        = 2;
T=250;
noOfObservations = 3000
order=2;
edge_probability=.001;
m_X=zeros(noOfNodes,noOfObservations);
 m_X_estimated=zeros(noOfNodes,noOfObservations);
%% data generation
    m_X(1,:)=.5*randn(1,noOfObservations);
    m_X(2,:)= .5*randn(1,noOfObservations);
    m_X(3,:)= .5*randn(1,noOfObservations);
    error_optimal=m_X;
    a(1:noOfObservations)=1;
    b(1:noOfObservations)=.5;
    c(1:noOfObservations)=.5;
    step_co=.2;
    w_T=zeros(noOfObservations,1);
    for ii=5:noOfObservations
        
         m_X(4,ii)= a(ii)*exp(-1*(m_X(2,ii-1)-m_X(2,ii-2))^2)+b(ii)*exp(-1*(m_X(2,ii-1)-m_X(2,ii-3))^2)+c(ii)*exp(-1*(m_X(2,ii-1)-m_X(2,ii-4))^2);
        m_X(5,ii)= a(ii)*exp(-1*(m_X(3,ii-1)-m_X(3,ii-2))^2)+b(ii)*exp(-1*(m_X(3,ii-1)-m_X(3,ii-3))^2)+c(ii)*exp(-1*(m_X(3,ii-1)-m_X(3,ii-4))^2);
       error_optimal(4,ii)=0;
       error_optimal(5,ii)=0;
       
       a(ii+1)=a(ii)+step_co*sin(ii);
       
       b(ii+1)=b(ii)+step_co*sin(ii);
       c(ii+1)=c(ii)+step_co*sin(ii);
       w_T(ii)=w_T(ii-1)+((a(ii+1)-a(ii))^2+(b(ii+1)-b(ii))^2+(c(ii+1)-c(ii))^2);
    end
% m_X=normalize(m_X,2);
[n_row_mx,nTimeInstants]=size(m_X);
RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = 5;
RFObj.filtOrder = 2; % we can try a higher order later
RFObj.lambda    = 50/75000;
RFObj.NoOfRF    =50;
RFObj.vsigma    =1;
  RFObj.forgettingFactor=.98;
 RFObj.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj.eta       =1000
RFState_in = RFObj.initialize( 100,m_X( :,1:RFObj.filtOrder)');
cum_orginal=0;
cum_estimated=0;
 for t = RFObj.filtOrder+1:nTimeInstants-1
    mtemp= m_X(:, t);
    RFState_in = RFObj.update(RFState_in, mtemp);
    if t==2000
        fffff=4;
    end
     RF_m_predic_error=m_X(:,t+1)-RFState_in.predictManyFromBuffer(1)';
     m_X_estimated(:,t+1)=RFState_in.predictManyFromBuffer(1)';
    
  Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
  % Psuedo_Adj_real=(sum(adjacency_real(:,:,:,1:T).^2,4)).^0.5;
   orginal_loss =sum(error_optimal(4:5,t+1).^2,1)+RFObj.lambda *.75;
    estimated_loss=sum(RF_m_predic_error(4:5).^2)+RFObj.lambda *sum(Psuedo_Adj(:));
   cum_orginal=cum_orginal+orginal_loss;
   cum_estimated=cum_estimated+estimated_loss;
   regret(t)=(cum_estimated-cum_orginal);
% temp_eig(t)=RFState_in.eig_value_phi;
 end
  Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        figure
        imagesc(B/max(B(:)))
        tit=title('RF-NLTISO-Estimate at t= 2990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
figure
plot(regret)
hold on
plot(w_T)