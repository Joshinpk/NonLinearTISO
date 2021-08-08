% Dictionary based NL-TISO is implemented in this code. This implementation
% builds dictionary based on the simalarity measured using kernal
% functions. The length of the dictionary is not fixed in this
% implementation.

% Comment#2 to test git branch
clc 
clear all
close all
filename = './data/sytem20.mat';
Data_struct = load(filename,'-mat');
M_x=Data_struct.m_X;
noOfNodes        = 24;
filtOrder        = 2;
noOfObservations = 1000
D=20
order=filtOrder;
gamma_RF=1e-4;
gamma_NLTISO=1e-2;
lamda=1/100000;
m_X=M_x(:,1:noOfObservations);
m_X=normalize(m_X,2)
Ksigma=10;
vsigma=1;
Gscale=(1/sqrt(2*pi)/Ksigma);
[pred_RF,alpha_RF,mse_RF,time_RF]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,D);
[pred_NLTISO,alpha_NLTISO,mse_NLTISO,time_NLTISO]=NLTISOF(order,lamda,gamma_NLTISO,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D);
%%%%Moving average
for t=12:noOfObservations
    tic
    movavg_pred(:,t)=sum(m_X(:,t-11:t-1),2)/10;
     timov(t)=toc;
     error_mov=(m_X(:,t)- movavg_pred(:,t)).^2;
    e_mov(:,t)=(error_mov)
    mse_mov(:,t)=sum(e_mov,2)/t;
   
end
%%% TIRSO
 addpath '/Users/rohanmoney/git_rohan/STInference/gsim/ProcessingBlocks/STFunctionEstimators/TirsoObjects/'
 

[n_row_mx,nTimeInstants]=size(m_X);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = 24;
tirsoObj.order     = 2; % we can try a higher order later
tirsoObj.regPar    = 1e-2;
tirsoObj.b_shrinkSelfLoops  = 0; % Bolstad
tirsoObj.forgettingFactor   = 0.8;
tirsoObj.h_stepsize         = @(ts)1/eigs(ts.m_Phi,1);
% initialize
tState_in = tirsoObj.initialize(0, m_X( :,1:tirsoObj.order)');
for t = tirsoObj.order+1:nTimeInstants
    tic
    mtemp= m_X(:, t);
    tState_in = tirsoObj.update(tState_in, mtemp);
    m_predic(:,:)=tState_in.predictManyFromBuffer(10)';
    m_prediction(1:24,t)= m_predic(:,1);
    error_tiso=(m_X(:,t)-m_prediction(:,t)).^2;
    e_tirso(:,t)=error_tiso;
    mse_tirso(:,t)=sum(e_tirso,2)/t;
    titirso(t)=toc;
end
Full_Ker_time=toc;
Psuedo_Adj=(sum(alpha_RF.^2,4)).^0.5;
plot(m_X(17,1:noOfObservations),'LineWidth',2)
hold on
plot(pred_RF(17,:),'LineWidth',2)
hold on

plot(movavg_pred(17,:))
hold on
plot(pred_NLTISO(17,:),'LineWidth',2) 
hold on
plot(m_prediction(17,:),'LineWidth',2)
legend('Orginal','RF-TISO','Moving Average','NL-TISO','TIRSO')
figure
% plot(mse_mov(18,:),'LineWidth',2)
% hold on
plot(mse_NLTISO(18,:),'LineWidth',2)
hold on
plot(mse_RF(18,:),'LineWidth',2)
hold on
% plot(mse_tirso(18,:),'LineWidth',2)
% legend('NL-TISO','RF-TISO')
% legend('Moving average','NL-TISO','RF-TISO','TIRSO')
figure
% plot(timov,'LineWidth',2)
% hold on
plot(time_NLTISO,'LineWidth',2)
legend('NL-TISO','RF-TISO')
hold on
plot(time_RF,'LineWidth',2)
% hold on
% plot(titirso,'LineWidth',2)
legend('NL-TISO','RF-TISO')
% legend('Moving average','NL-TISO','RF-TISO','TIRSO')
xlabel('t')
ylabel('computation time(s)')
title('Time taken for each prediction')

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
function [predt,alpha,mse,ti]= NLTISOF(order,lamda,gamma,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D)
Gscale=(1/sqrt(2*pi)/Ksigma);
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal=zeros(noOfNodes,filtOrder,noOfObservations);

eta=1/gamma;
%% new code
for t=filtOrder+1:noOfObservations
    tic
   for n1=1:noOfNodes
       m_X_P_pre=m_X(n1,t-filtOrder:t-1);
       m_X_P=flip(m_X_P_pre);
       m_X_P_big=repmat(m_X_P',1,t);
              Kernal(n1,:,1:t)=Gscale*exp(-.001*(m_X_P_big-m_X(n1,1:t)).^2);
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
   ti(t)=toc;
   error(:,t)=error_temp;
     mse(:,t)=sum(error,2)/t;
   
   (sum(alpha.^2,4)).^0.5;
end
end
function [predt,alpha,mse,ti]= NLTISORFF(order,lamda,gamma,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,D)
Kernal=zeros(noOfNodes,filtOrder,2*D);
v=vsigma*vsigma*randn(noOfNodes,filtOrder,D);

eta=1/gamma;

alpha=zeros(noOfNodes,noOfNodes,filtOrder,2*D);
tic
%%RF-TISO
for t=filtOrder+1:noOfObservations
    tic
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
    ti(t)=toc;
    error(:,t)=error_temp;
     mse(:,t)=sum(error,2)/t;
 
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