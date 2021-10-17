clc
clear all
%% Data generation
noOfNodes        = 5; %no of nodes
filtOrder        = 2; % process order
noOfObservations = 3000 % no of time samples
edge_probability=.1;    % probaility of non-zeros edge in graph


A=[];
for i=1:filtOrder
    temp1=.1*rand(noOfNodes,noOfNodes).*(rand(noOfNodes,noOfNodes)<edge_probability);
    A=[A;temp1];
end
%plot intial true dependency
figure
A1=A/max(max(A));
subplot(4,3,1)
imagesc(A/max(A(:)))
tit=title('True dependencies intial')
set(tit,'fontSize',24)
ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);

% Intialize kernel paramters
y_initial=1*randn(noOfNodes ,1);
y=[y_initial];
Ksigma=10;
gsigma=Ksigma;
gscale=(1/sqrt(2*pi)/Ksigma);
GausDenom=-1/2/Ksigma^2;
T_Beta=50*randn(filtOrder *noOfNodes,noOfNodes,noOfObservations);
M_x=zeros(noOfNodes ,noOfObservations );
M_x(:,1)=y_initial+.2*randn;
for ii=2:filtOrder
    temp=A(1:noOfNodes ,:)*y_initial;
    M_x(:,ii)=temp+.2*rand;
    if ii==2990
        M_x(:,ii)=0;
    end
    y_initial= M_x(:,ii);
    
end

for ii=filtOrder +1:noOfObservations
    
    
    for k=1:noOfNodes
        temp=0;
        for m=1:filtOrder 
            for l=1:noOfNodes
                temp=temp+ A(noOfNodes *(m-1)+k,l)*((1*gscale*exp(-1000*(M_x(l,1:ii-m-1)-M_x(l,ii-m)).^2/(10*gsigma*gsigma))*reshape(T_Beta(noOfNodes *(m-1)+k,l,1:ii-m-1),[ii-m-1,1])));
                
            end
        end
        
        M_x(k,ii)=temp+.2*randn;
        
    end
    if ii==1000
        no_non_zero=nnz(A);
        no_to_change=round(no_non_zero*.3);%30 percentgae of edges change
        non_zero_entries=find(A);
        entries_to_change= randperm(numel(non_zero_entries),no_to_change)
        temp=non_zero_entries(entries_to_change)
        A(temp)=0;
        subplot(4,3,2)
        imagesc(A/max(A(:)))
        tit=title('True dependencies intial')
        set(tit,'fontSize',24)
        ax=gca;
        ax.XAxis.FontSize = 24;
        ax.YAxis.FontSize = 24;
        c = colorbar;
        w = c.LineWidth;
        c.LineWidth = 1.5;
        set(c,'fontSize',24);
        %
        %             end
    end
    if ii==2000
        no_non_zero=nnz(A);
        no_to_change=round(no_non_zero*.3);%30 percentgae of edges change
        non_zero_entries=find(A);
        entries_to_change= randperm(numel(non_zero_entries),no_to_change)
        temp=non_zero_entries(entries_to_change)
        A(temp)=0;
        subplot(4,3,3)
        imagesc(A/max(A(:)))
        tit=title('True dependencies intial')
        set(tit,'fontSize',24)
        ax=gca;
        ax.XAxis.FontSize = 24;
        ax.YAxis.FontSize = 24;
        c = colorbar;
        w = c.LineWidth;
        c.LineWidth = 1.5;
        set(c,'fontSize',24);
        %             end
    end
    
end
    m_X=M_x; % Data matrix
    
    
%% RFNL-TIRSO ESTIMATION AND PLOTTING   
addpath '/Users/rohanmoney/git_rohan/NonLinearTISO'

    [n_row_mx,nTimeInstants]=size(m_X);
RFObj_tirso = RF_nltirso; % set tirso object up
RFObj_tirso.noOfNodes = 5;
RFObj_tirso.filtOrder = 2; % we can try a higher order later
RFObj_tirso.lambda    = 100/50;
RFObj_tirso.NoOfRF    =30;
RFObj_tirso.vsigma    =1*ones(5,1);
  RFObj_tirso.forgettingFactor=.98;
 RFObj_tirso.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj_tirso.eta       =500
RFState_in_tirso = RFObj_tirso.initialize( 10,m_X( :,1:RFObj_tirso.filtOrder)');
 for t = RFObj_tirso.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in_tirso = RFObj_tirso.update(RFState_in_tirso, mtemp);
    RF_m_predic_tirso(:,:)=RFState_in_tirso.predictManyFromBuffer(1)';
    RF_m_prediction_tirso(1:5,t)= RF_m_predic_tirso(:,1);
    if t== 990
        Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5
        B=[ Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1)];
       subplot(4,3,4)
        imagesc(B/max(B(:)))
        tit=title('RFNL-TIRSO-Estimate at t= 990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
    if t== 1990
        Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5
        B=[ Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1)];
        subplot(4,3,5)
        imagesc(B/max(B(:)))
        tit=title('RFNL-TIRSO-Estimate at t= 1990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
    if t== 2990
        Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5
        B=[ Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1)];
        subplot(4,3,6)
        imagesc(B/max(B(:)))
        tit=title('RFNL-TIRSO-Estimate at t= 2990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
 end
%  Psuedo_Adj_real=(sum(adjacency_real(:,:,:,1:T).^2,4)).^0.5;
  Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5
        B=[Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1) ];
       
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
Y_n=zeros(1000,1);
for ii=filtOrder:1000
    Y_n(ii)=m_X(1,ii)+Y_n(ii-1)*RFObj_tirso.forgettingFactor;
end
nn=10
cvx_begin
    variable x(nn)
    variable x_r(10)
    minimize(  (Y_n(ii) - x'*ones(nn,1))*(Y_n(ii) - x'*ones(nn,1))+norm(x_r,1) ) b
 
