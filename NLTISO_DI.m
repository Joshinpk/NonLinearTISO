
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
noOfObservations = 250 
order=2;
gamma=1e-1;
lamda=1/10000000;
m_X=M_x(:,1:noOfObservations);
Ksigma=10;
Gscale=(1/sqrt(2*pi)/Ksigma);
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal=zeros(noOfNodes,filtOrder,noOfObservations);                                     

eta=1/gamma;
%% Dictionery
m_X_Dict(:,1:filtOrder)=m_X(:,1:filtOrder);
    [~  ,T]=size(m_X_Dict);
tic
for t=filtOrder+1:noOfObservations   
    flag=0;
    for n11=1:noOfNodes
        temp= Gscale*exp(-100*(m_X(n11,t)-m_X(n11,1:t-1)).^2);%find similarity with previous samples
        if max(temp)<.035 %check weather similarity is grater than treshold
            if flag==0
               T=T+1;% In crease legth of m_X_Dict
               flag=1;
               m_X_Dict(:,T)=m_X(:,t);
               break
            end
        m_X_Dict(:,T)=  m_X(:,t);     
        end
     end
   for n1=1:noOfNodes
       m_X_P_pre=m_X(n1,t-filtOrder:t-1);
       m_X_P=flip(m_X_P_pre);
       m_X_P_big=repmat(m_X_P',1,T);
       Kernal_Dict(n1,:,1:T)=Gscale*exp(-.001*(m_X_P_big-m_X_Dict(n1,1:T)).^2);
   end 
   Kernal_Size_Dict=size(Kernal_Dict);
   K_vec_Dict=vec(Kernal_Dict(:,:,1:T));  
   for n1=1:noOfNodes
       try
           predt_dict(n1,t)=(vec(alpha(n1,:,:,1:T)))'*K_vec_Dict;
           grad_vec=(vec(alpha(n1,:,:,1:T)))'*K_vec_Dict*K_vec_Dict'-K_vec_Dict'*m_X(n1,T);
           error_temp(n1)=(predt_dict(n1,t)-m_X(n1,t))^2;
           grad=reshape(grad_vec,noOfNodes,filtOrder,T);
           g=eta*squeeze(alpha(n1,:,:,1:T))-grad;
           alpha(n1,:,:,1:T)= SoftThresold(g,T,noOfNodes,eta,lamda,n1,order);
       catch
           debug=1;
       end
   end
%    error(t)=sum(error_temp(1))/5;
end
Dict_time=toc;
Psuedo_Adj_Dict=(sum(alpha.^2,4)).^0.5;
plot(m_X(18,1:noOfObservations),'LineWidth',2)
hold on
plot(predt_dict(18,:),'LineWidth',2)
hold on

%%%%%%%%%%%%%%%
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
plot(predt(18,:),'LineWidth',2)
legend('original','prediction (dict)','prediction (full)')
alpha_full=alpha(:);
Num_samples_Dict=Kernal_Size_Dict(3);
Num_samples_Ful=Kernal_Size_Full(3);
title_str=strcat('NoSamplesDict=',num2str(Num_samples_Dict),', ',...
    'NoSamplesFull=',num2str(Num_samples_Ful))
title(title_str)

Kernal_Size_Dict
Kernal_Size_Full
Full_Ker_time
Dict_time

function [alpha_temp] =SoftThresold(g,t,noOfNodes,eta,lamda,n1,order)
alpha_temp=zeros(size(g));
    for n2=1:noOfNodes 
       for m=1:order
           alpha_temp(n2,m,:)=g(n2,m,:)/eta;
           alpha_temp(n2,m,:)=g(n2,m,:)/eta*max(0,(1-lamda/norm(reshape(g(n2,m,:),1,t)))); 
       end 
    end
end