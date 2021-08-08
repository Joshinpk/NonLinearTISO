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
noOfNodes        = 5;
filtOrder        = 2;
T=1;
noOfObservations = 3000
D=30;
order=2;
gamma_RF=1e-3;
gamma_NLTISO=1e-1;
lamda=1/10;
edge_probability=.1;
A=[];
for kk=1:1
for i=1:order
    temp1=10*rand(noOfNodes,noOfNodes).*(rand(noOfNodes,noOfNodes)<edge_probability);
    A=[A;temp1];
end
% A=zeros(10,5);
% [len_A,bred_A]=size(A)
% ad_in=randi([1,len_A*bred_A],1,1);
% A(ad_in)=1*rand;
%  ad_in=randi([1,len_A*bred_A],1,1);
 A(5)=10*rand;
%  A(14)=8;

subplot(4,3,1)
imagesc(A/max(A(:)))
title('True dependencies intial')
colorbar;
alpha=.1*rand(noOfNodes,noOfObservations);
y_initial=1*randn(noOfNodes ,1);
y=[y_initial];
Ksigma=10;
gsigma=Ksigma;
Gscale=(1/sqrt(2*pi)/Ksigma);
gscale=Gscale;
GausDenom=-1/2/Ksigma^2;
T=250
T_Beta=3*rand(order*noOfNodes,noOfNodes,T);

M_x=zeros(noOfNodes ,noOfNodes );
M_x(:,1)=y_initial+.1*randn;
for ii=2:T+order
    temp=A(1:noOfNodes ,:)*y_initial;
    M_x(:,ii)=temp+.1*rand;
    y_initial= M_x(:,ii);
end

for ii=T+order+1:noOfObservations
    A5(ii)=A(5);
%    list_nz=find(A);
%         len=length(list_nz);
%         rem_in=randi([1,len],1,1);
%         rem=list_nz(rem_in);
        A(5)=A(5)+1*sin(ii*.003);
    
    for k=1:noOfNodes
        temp=0;
        for m=1:order
            for l=1:noOfNodes
                temp=temp+ A(noOfNodes *(m-1)+k,l)*((1*gscale*exp(-10000000*(M_x(l,ii-m-T:ii-m-1)-M_x(l,ii-m)).^2/(10*gsigma*gsigma))*reshape(T_Beta(noOfNodes *(m-1)+k,l,:),[T,1]))+(0*gscale*exp(-1000000*(M_x(l,ii-m-T:ii-m-1)-M_x(l,ii-m)).^2/(2*gsigma*gsigma))*reshape(T_Beta(noOfNodes *(m-1)+k,l,:),[T,1])));
                
            end
        end
        M_x(k,ii)=temp+.05*randn;
        
    end
    
end
m_X=M_x;
% m_X=normalize(m_X,2);
Ksigma=10;
vsigma=[10 10 10 10 10];
for i=1:noOfNodes
    v(i,:,:)=vsigma(i)*vsigma(i)*randn(1,filtOrder,D);
end
Gscale=(1/sqrt(2*pi)/Ksigma);
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal=zeros(noOfNodes,filtOrder,noOfObservations);

[pred_RF,alpha_RF,mse_RF(kk,:),B15]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,10);
[pred_RF,alpha_RF,mse_RF(kk,:),B25]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,30);
[pred_RF,alpha_RF,mse_RF(kk,:),B35]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,40);
[pred_NLTISO,alpha_NLTISO,mse_NL(kk,:),BNL5]=NLTISOF(order,lamda,gamma_NLTISO,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D);
%%%%Moving average
for t=12:noOfObservations
    movavg_pred(:,t)=sum(m_X(:,t-11:t-1),2)/10;
end
%% Tiso
addpath '/Users/rohanmoney/git_rohan/STInference/gsim/ProcessingBlocks/STFunctionEstimators/TirsoObjects/'


[n_row_mx,nTimeInstants]=size(M_x);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = 5;
tirsoObj.order     = 2; % we can try a higher order later
tirsoObj.regPar    = 1e-1;
tirsoObj.b_shrinkSelfLoops  = 0; % Bolstad
tirsoObj.forgettingFactor   = 0.8;
tirsoObj.h_stepsize         = @(ts)1/eigs(ts.m_Phi,1);
% initialize
tState_in = tirsoObj.initialize(0, M_x( :,1:tirsoObj.order)');
for t = tirsoObj.order+1:nTimeInstants
    if t==990
        subplot(4,3,4)
        imagesc(tState_in.m_A/max(tState_in.m_A(:)))
        colorbar;
        title('TIRSO estimate at t=990')
        
    end
    if t==1990
        subplot(4,3,5)
        imagesc(tState_in.m_A/max(tState_in.m_A(:)))
        colorbar;
        title('TIRSO estimate at t=1990')
        
    end
    if t==2990
         subplot(4,3,6)
        imagesc(tState_in.m_A/max(tState_in.m_A(:)))
        colorbar;
        title('TIRSO estimate at t=2990')
        
    end
    mtemp= M_x(:, t);
    tState_in = tirsoObj.update(tState_in, mtemp);
    m_predic(:,:)=tState_in.predictManyFromBuffer(10)';
     m_prediction(1:5,t)= m_predic(:,1);
    error_tiso=(m_X(:,t)-m_prediction(:,t)).^2;
    e_tirso(:,t)=error_tiso;
    mse_tirso(:,t)=sum(error_tiso)/sum(m_X(:,t).^2);
    Bt=tState_in.m_A/max(tState_in.m_A(:));
    Bt5(t)=Bt(5);
end
% mse_t(kk,:)=mse_tirso;
end
figure

plot(A5/max(A5),'LineWidth',2)
hold on
plot(Bt5/max(Bt5),'LineWidth',2)
hold on
 plot(BNL5/max(BNL5),'LineWidth',2)
 hold on
% plot(B15/max(B15),'LineWidth',2)
% hold on
% plot(B25/max(B25),'LineWidth',2)
% hold on
plot(B35/max(B35),'LineWidth',2)
title('')
lgd2=legend('Orginal','TIRSO','NL-TISO','RF-NLTISO')
set(lgd2,'fontsize',18)
ylabel('Magnitude','FontSize',18)
xlabel('t','FontSize',18)
% hold on
% plot(B25)
% hold on
% plot(B35)
figure
plot(sum(mse_RF(:,300:end))/10,'LineWidth',2)
 hold on
 plot(sum(mse_NL(:,300:end))/10,'LineWidth',2)
 hold on
%  plot( sum(mse_t(:,300:end))/10,'LineWidth',2)
% legend('RF-NLTISO','NLTISO','TIRSO')
% Full_Ker_time=toc;
% Psuedo_Adj=(sum(alpha_RF.^2,4)).^0.5
% plot(pred_RF(1,:),'LineWidth',2)
% hold on
% plot(m_X(1,1:noOfObservations),'LineWidth',2)
% hold on
% plot(movavg_pred(1,:))
% hold on
% % plot(pred_NLTISO(15,:),'LineWidth',2)
% legend('prediction (RF)','orginal','Moving Average')
% alpha_full=alpha(:);
% Num_samples_Dict=Kernal_Size_Dict(3);
% Num_samples_Ful=Kernal_Size_Full(3);
% title_str=strcat('NoSamplesDict=',num2str(Num_samples_Dict),', ',...
%     'NoSamplesFull=',num2str(Num_samples_Ful))
% title(title_str)
%
% Kernal_Size_Dict
% Kernal_Size_Full
% Full_Ker_time
% Dict_time
function [predt,alpha,mse,B5]= NLTISOF(order,lamda,gamma,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D)
Gscale=(1/sqrt(2*pi)/Ksigma);
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal=zeros(noOfNodes,filtOrder,noOfObservations);

eta=1/gamma;
%% new code
for t=filtOrder+1:noOfObservations
    
        Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        B5(t)=B(5);
        
        
    for n1=1:noOfNodes
        m_X_P_pre=m_X(n1,t-filtOrder:t-1);
        m_X_P=flip(m_X_P_pre);
        m_X_P_big=repmat(m_X_P',1,t);
        Kernal(n1,:,1:t)=Gscale*exp(-1000*(m_X_P_big-m_X(n1,1:t)).^2);
    end
    K_vec=vec(Kernal(:,:,1:t));
    for n1=1:noOfNodes
        predt(n1,t)=(vec(alpha(n1,:,:,1:t)))'*K_vec;
        grad_vec=(vec(alpha(n1,:,:,1:t)))'*K_vec*K_vec'-K_vec'*m_X(n1,t);
        error_temp(n1)=(predt(n1,t)-m_X(n1,t))^2;
        grad=reshape(grad_vec,noOfNodes,filtOrder,t);
        g=eta*squeeze(alpha(n1,:,:,1:t))-grad;
        alpha(n1,:,:,1:t)= SoftThresold(g,t,noOfNodes,eta,lamda,n1,order);
    end
     error(:,t)=error_temp;
     mse(t)=sum(error_temp)/(sum(m_X(:,t).^2));
    (sum(alpha.^2,4)).^0.5;
end
end
function [predt,alpha,mse,B5]= NLTISORFF(order,lamda,gamma,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,D)
Kernal=zeros(noOfNodes,filtOrder,2*D);

for i=1:noOfNodes
    v(i,:,:)=vsigma(i)*vsigma(i)*randn(1,filtOrder,D);
end

eta=1/gamma;

alpha=zeros(noOfNodes,noOfNodes,filtOrder,2*D);
tic
%%RF-TISO
for t=filtOrder+1:noOfObservations
    
        Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        
       B5(t)=B(5);
      
    for n1=1:noOfNodes
        
        for tau=1:filtOrder
            m_X_P_pre=m_X(n1,t-tau);
            v_vec=v(n1,tau,:);
            
            
            
            Kernal(n1,tau,1:2*D)=[reshape((sin(m_X_P_pre*v_vec)),[1,D]) reshape((cos(m_X_P_pre*v_vec)),[1,D])];
        end
    end
    
    Kernal_Size_Full=size(Kernal);
    K_vec=vec(Kernal(:,:,1:2*D));
    for n1=1:noOfNodes
        predt(n1,t)=(vec(alpha(n1,:,:,1:2*D)))'*K_vec;
        grad_vec=(vec(alpha(n1,:,:,1:2*D)))'*K_vec*K_vec'-K_vec'*m_X(n1,t);
        error_temp(n1)=(predt(n1,t)-m_X(n1,t))^2;
        grad=reshape(grad_vec,noOfNodes,filtOrder,2*D);
        g=eta*squeeze(alpha(n1,:,:,1:2*D))-grad;
        alpha(n1,:,:,1:2*D)= SoftThresold(g,2*D,noOfNodes,eta,lamda,n1,order);
    end
    error(:,t)=error_temp;
     mse(t)=sum(error_temp)/sum(m_X(:,t).^2);
end
end
function [alpha_temp] =SoftThresold(g,t,noOfNodes,eta,lamda,n1,order)
alpha_temp=zeros(size(g));
for n2=1:noOfNodes
    for m=1:order
        alpha_temp(n2,m,:)=g(n2,m,:)/eta;
        alpha_temp(n2,m,:)=g(n2,m,:)/eta*max(0,(1-lamda/norm(reshape(g(n2,m,:),1,t))));
    end
end
end