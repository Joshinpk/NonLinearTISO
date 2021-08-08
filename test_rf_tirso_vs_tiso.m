clc
clear all
close all
% filename = './data/sytem20.mat';
% Data_struct = load(filename,'-mat');
% M_x=Data_struct.m_X;
noOfNodes        = 5;
filtOrder        = 2;
T=250;
noOfObservations = 3000
order=2;
edge_probability=.1;
m_X=zeros(noOfNodes,noOfObservations);
%% data generation
    m_X(1,:)=.5*randn(1,noOfObservations);
    m_X(2,:)= .5*randn(1,noOfObservations);
    m_X(3,:)= .5*randn(1,noOfObservations);
    
    for ii=3:noOfObservations
        m_X(5,ii)= exp(-1*m_X(1,ii-1)*m_X(1,ii-1))+exp(-1*m_X(3,ii-2)*m_X(3,ii-2));
        m_X(4,ii)= exp(-1*m_X(2,ii-1)*m_X(2,ii-1));
    end
% m_X=normalize(m_X,2);
[n_row_mx,nTimeInstants]=size(m_X);
RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = 5;
RFObj.filtOrder = 2; % we can try a higher order later
RFObj.lambda    = 50/50;
RFObj.NoOfRF    =50;
RFObj.vsigma    =1;
  RFObj.forgettingFactor=.98;
 RFObj.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj.eta       =1000
RFState_in = RFObj.initialize( 10,m_X( :,1:RFObj.filtOrder)');
 for t = RFObj.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in = RFObj.update(RFState_in, mtemp);
    if t==2000
        fffff=4;
    end
   
 end
  Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        
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