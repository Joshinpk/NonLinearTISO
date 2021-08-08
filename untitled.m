clc 
clear all
close all
filename = './data/sytem20.mat';
Data_struct = load(filename,'-mat');
M_x=Data_struct.m_X;
noOfNodes        = 24;
filtOrder        = 2;
noOfObservations = 800
order=2;
gamma=1e-2;
lamda=1/10000000;
m_X=M_x(:,1:noOfObservations);
m_X=normalize(m_X,2)
Ksigma=10;
Gscale=(1/sqrt(2*pi)/Ksigma);
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal=zeros(noOfNodes,filtOrder,noOfObservations);                                     

eta=1/gamma;
D=100;
sigmapi=10;
v=sigmapi*sigmapi*randn(noOfNodes,noOfNodes,filtOrder,D)
% Without dictionary
%%%%%%%%%%%%%%%
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
tic
for t=filtOrder+1:noOfObservations
    for n1=1:noOfNodes  
        m_X_P_pre=m_X(n1,t-filtOrder:t-1);
        m_X_P=flip(m_X_P_pre);
        m_X_P_big=repmat(m_X_P',1,t);
        Kernal(n1,:,1:t)=Gscale*exp(-.001*(m_X_P_big-m_X(n1,1:t)).^2);
    end 
    Kernal_Size_Full=size(Kernal);
    K_vec=vec(Kernal(:,:,1:t));
    for n1=1:noOfNodes
        predt(n1,t)=(vec(alpha(n1,:,:,1:t)))'*K_vec; 
        grad_vec=(vec(alpha(n1,:,:,1:t)))'*K_vec*K_vec'-K_vec'*m_X(n1,t);
        error_temp(n1)=(predt(n1,t)-m_X(n1,t))^2;
        grad=reshape(grad_vec,noOfNodes,filtOrder,t);
        g=eta*squeeze(alpha(n1,:,:,1:t))-grad;
        alpha(n1,:,:,1:t)= SoftThresold(g,t,noOfNodes,eta,lamda,n1,order);
    end
%         error(t)=sum(error_temp(1))/5;
end
Full_Ker_time=toc;
Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
plot(predt(7,:),'LineWidth',2)
legend('original','prediction (dict)','prediction (full)')
alpha_full=alpha(:);
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

function [alpha_temp] =SoftThresold(g,t,noOfNodes,eta,lamda,n1,order)
alpha_temp=zeros(size(g));
    for n2=1:noOfNodes 
       for m=1:order
           alpha_temp(n2,m,:)=g(n2,m,:)/eta;
           alpha_temp(n2,m,:)=g(n2,m,:)/eta*max(0,(1-lamda/norm(reshape(g(n2,m,:),1,t)))); 
       end 
    end
end