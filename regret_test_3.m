% Dictionary based NL-TISO is implemented in this code. This implementation
% builds dictionary based on the simalarity measured using kernal
% functions. The length of the dictionary is not fixed in this
% implementation.

% Comment#2 to test git branch
clc
clear all
close all
% filename = './data/sytem20.mat';
% Data_struct = load(filename,'-mat');
% M_x=Data_struct.m_X;
noOfNodes        = 3;
filtOrder        = 2;
noOfObservations = 10000
order=2;
A=[];
edge_probability=.2;
kernal_forget(1)=1;
T=1000;
for ii=2:T
    kernal_forget(ii)= kernal_forget(ii-1)*kernal_forget(1);
end
for ii=1:order
    temp1=.1*rand(noOfNodes,noOfNodes).*(rand(noOfNodes,noOfNodes)<edge_probability);
    A(:,:,ii)=temp1;
end
 
Ksigma=10;
gsigma=Ksigma;
Gscale=(1/sqrt(2*pi)/Ksigma);
gscale=Gscale;
GausDenom=-1/2/Ksigma^2;

T_Beta=50*randn(noOfNodes,noOfNodes,order,T);
adjacency_real=zeros(noOfNodes,noOfNodes,order,T);
M_x=zeros(noOfNodes ,noOfObservations );
M_x_optimal=zeros(noOfNodes ,noOfObservations );
error_optimal=zeros(noOfNodes ,noOfObservations );
M_x(:,1:T)=randn(noOfNodes,T);
W_k=zeros(noOfObservations,1);
adj_norm=zeros(noOfObservations,1);
for kk=1:noOfNodes
    for ii=1:noOfNodes
        for jj=1:order
            adjacency_real(kk,ii,jj,:)=A(kk,ii,jj)*(T_Beta(kk,ii,jj,:).*reshape(kernal_forget,[1,1,1,T])) ;
        end
    end
end
for ii=order+1+T:noOfObservations
    
    temp_norm=0;
    for k=1:noOfNodes
        temp=0;
        for m=1:order
           
            for l=1:noOfNodes
                temp=temp+ ((exp(-1000*(flip(M_x(l,ii-m-1-T+1:ii-m-1))-M_x(l,ii-m)).^2/(1*gsigma*gsigma))))*vec(adjacency_real(k,l,m,:));
                adjacency_prev=adjacency_real(k,l,m,:);
                if norm(vec(adjacency_real(k,l,m,:)))>0
                    
%                    adjacency_real(k,l,m,:)=adjacency_real(k,l,m,:)+.01*sin(ii);
                end
                temp_norm=temp_norm+(norm(vec(adjacency_real(k,l,m,:))-vec(adjacency_prev)))^2;
            end
        end
        adj_norm(ii)= (temp_norm)^.5;
        W_k(ii)=W_k(ii-1)+adj_norm(ii);
        M_x_optimal(k,ii)=temp;
        M_x(k,ii)=temp+.5*randn;
        error_optimal(k,ii)=M_x(k,ii)-M_x_optimal(k,ii);
        
    end
    
end

m_X=M_x;
    [n_row_mx,nTimeInstants]=size(m_X);
RFObj_tirso = RF_nltirso; % set tirso object up
RFObj_tirso.noOfNodes = 3;
RFObj_tirso.filtOrder = 2; % we can try a higher order later
RFObj_tirso.lambda    = 50/50;
RFObj_tirso.NoOfRF    =30;
RFObj_tirso.vsigma    =1*ones(3,1);
  RFObj_tirso.forgettingFactor=.98;
 RFObj_tirso.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj_tirso.eta       =1000
RFState_in_tirso = RFObj_tirso.initialize( 10,m_X( :,1:RFObj_tirso.filtOrder)');
cum_orginal=0;
cum_estimated=0;
 for t = RFObj_tirso.filtOrder+1:nTimeInstants-1
    mtemp= m_X(:, t);
    RFState_in_tirso = RFObj_tirso.update(RFState_in_tirso, mtemp);
    if t==5000
        fffff=4;
    end
   %temp_eig(t)= RFState_in_tirso.eig_value_phi;
   RF_m_predic_error=m_X(:,t+1)-RFState_in_tirso.predictManyFromBuffer(1)';
    
   Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5;
   Psuedo_Adj_real=(sum(adjacency_real(:,:,:,1:T).^2,4)).^0.5;
   orginal_loss =sum(error_optimal(:,t+1).^2,1)+RFObj_tirso.lambda *sum(Psuedo_Adj_real(:));
    estimated_loss=sum(RF_m_predic_error.^2)+RFObj_tirso.lambda *sum(Psuedo_Adj(:));
   cum_orginal=cum_orginal+orginal_loss;
   cum_estimated=cum_estimated+estimated_loss;
   regret(t)=(cum_estimated-cum_orginal);
  
 end
  Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5
     Psuedo_Adj_real=(sum(adjacency_real(:,:,:,1:T).^2,4)).^0.5
        B=[Psuedo_Adj_real(:,:,1);Psuedo_Adj_real(:,:,2) ];
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
plot(m_X')
figure
 plot(regret)
 hold on
 plot(W_k)
 figure
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