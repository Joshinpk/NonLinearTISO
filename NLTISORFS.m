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
noOfNodes        = 5;
filtOrder        = 2;
T=1;
noOfObservations = 9000
D=60;
order=2;
gamma=1e+1;
lamda=1/10;
edge_probability=.05;
A=zeros(10,5);
    A(1,5)=.9;
      A(9,3)=.7;
    A(1,2)=.8;
     A(8,4)=.9;
     A=[];
for i=1:order
    temp1=1*rand(noOfNodes,noOfNodes).*(rand(noOfNodes,noOfNodes)<edge_probability);
    A=[A;temp1];
end
  
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
    
    for k=1:noOfNodes 
        temp=0;
        for m=1:order
        for l=1:noOfNodes 
              temp=temp+ A(noOfNodes *(m-1)+k,l)*((gscale*exp(-10000000*(M_x(l,ii-m-T:ii-m-1)-M_x(l,ii-m)).^2/(10*gsigma*gsigma))*reshape(T_Beta(noOfNodes *(m-1)+k,l,:),[T,1]))+(0*gscale*exp(-1000000*(M_x(l,ii-m-T:ii-m-1)-M_x(l,ii-m)).^2/(2*gsigma*gsigma))*reshape(T_Beta(noOfNodes *(m-1)+k,l,:),[T,1])));
              
        end
        end
         M_x(k,ii)=temp+.05*randn;

    end

end
m_X=M_x;
m_X=normalize(m_X,2);
Ksigma=10;
vsigma=[10 10 10 10 10];
Gscale=(1/sqrt(2*pi)/Ksigma);
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal=zeros(noOfNodes,filtOrder,noOfObservations);
for i=1:noOfNodes
v(i,:,:)=vsigma(i)*vsigma(i)*randn(1,filtOrder,D);
end

eta=1000;

alpha=zeros(noOfNodes,noOfNodes,filtOrder,2*D);
tic
for t=filtOrder+1:noOfObservations
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
%         error(t)=sum(error_temp(1))/5;
end
Full_Ker_time=toc;
Psuedo_Adj=(sum(alpha.^2,4)).^0.5
B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
A
plot(predt(1,:),'LineWidth',2)
hold on
plot(m_X(1,1:noOfObservations),'LineWidth',2)
legend('prediction (RF)','orginal')
figure
subplot(2,2,1)
imagesc(A/max(A(:)))
colorbar;
title('True causal dependencies')
subplot(2,2,2)
imagesc(B/max(B(:)))
colorbar;
title('NL-TISO estimate')

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

function [alpha_temp] =SoftThresold(g,t,noOfNodes,eta,lamda,n1,order)
alpha_temp=zeros(size(g));
    for n2=1:noOfNodes 
        for m=1:order
           alpha_temp(n2,m,:)=g(n2,m,:)/eta;
           alpha_temp(n2,m,:)=g(n2,m,:)/eta*max(0,(1-lamda/norm(reshape(g(n2,m,:),1,t)))); 
       end 
    end
end