% Dictionary based NL-TISO is implemented in this code. This implementation
% builds dictionary based on the simalarity measured using kernal
% functions. The length of the dictionary is not fixed in this
% implementation.

% Comment#2 to test git branch
clc 
clear all
close all
T = readtable('time_series_covid19_confirmed_US.csv');
 header_colum_num=11;
 Alabama_start_index=2;
 Alabama_end_index=300;
 
 
 start_index=Alabama_start_index;
 end_index=Alabama_end_index;
 T_selected=T(start_index:end_index,header_colum_num+1:end);
 TMx_cell=table2array(T_selected);
 M_xtotal=(str2double(TMx_cell))';
 M_x= M_xtotal(:,6:21)';

noOfNodes        = 15;
filtOrder        = 24;
noOfObservations = 299
D=50
order=filtOrder;
gamma_RF=1e-4;
gamma_NLTISO=1e-2;
lamda=1/5;
m_X=M_x(:,1:noOfObservations);
 m_X=normalize(m_X,2)
Ksigma=10;
vsigma=100;
Gscale=(1/sqrt(2*pi)/Ksigma);
[pred_RF,alpha_RF,timerf]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,D);
% [pred_NLTISO,alpha_NLTISO,timenl]=NLTISOF(order,lamda,gamma_NLTISO,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D)
%%%%Moving average
for t=12:noOfObservations
    movavg_pred(:,t)=sum(m_X(:,t-11:t-1),2)/10;
end
    
Full_Ker_time=toc;
Psuedo_Adj=(sum(alpha_RF.^2,4)).^0.5
B=[];

for ii=1:24
    B=[B;Psuedo_Adj(1,:,ii)];
end
plot(pred_RF(1,:),'LineWidth',2)
hold on
plot(m_X(1,1:noOfObservations),'LineWidth',2)
hold on
plot(movavg_pred(1,:))
hold on
% plot(pred_NLTISO(15,:),'LineWidth',2) 
legend('prediction (RF)','orginal','Moving Average')
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
function [predt,alpha,ti]= NLTISOF(order,lamda,gamma,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D)
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
   ti(t)=toc
     error(t)=sum(error_temp(1))/5;
   (sum(alpha.^2,4)).^0.5;
end
end
function [predt,alpha,ti]= NLTISORFF(order,lamda,gamma,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,D)
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
    ti(t)=toc
%         error(t)=sum(error_temp(1))/5;
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