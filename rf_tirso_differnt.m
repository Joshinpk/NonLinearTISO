



clear all


%% Initiazlization rf-nltirso
noOfNodes = 5;
filtOrder = 2; % we can try a higher order later
lambda    = 10/10;
D    =30;
NoOfRF=D;
vsigma    =10;
forgettingFactor=.98;
h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
eta       =1000;
edge_probability=.1;
noOfObservations=2000;
m_Phi=eye(noOfNodes*filtOrder*2*NoOfRF);
m_R=zeros(noOfNodes*filtOrder*2*NoOfRF,noOfNodes);
m_X=zeros(noOfNodes,noOfObservations);
%% data generation
    m_X(1,:)=randn(1,noOfObservations);
    m_X(2,:)= .2*randn(1,noOfObservations);
    m_X(3,:)= .5*randn(1,noOfObservations);
    m_X(4,:)= .05*randn(1,noOfObservations);
    for ii=3:noOfObservations
        m_X(5,ii)= sin(m_X(1,ii-1))+sin(m_X(3,ii-2));
    end

for i=1:noOfNodes
    v(i,:,:)=vsigma*vsigma*randn(1,filtOrder,D);
end
alpha=zeros(noOfNodes,noOfNodes,filtOrder,2*D);
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
m_Phi=eta*m_Phi+(1-eta)*K_vec*K_vec';
for n1=1:noOfNodes
    m_R(:,n1)=eta*m_R(:,n1)+(1-eta)*K_vec*m_X(n1,t);
    grad_vec=m_Phi*(vec(alpha(n1,:,:,1:2*D)))-m_R(:,n1);
    grad=reshape(grad_vec,noOfNodes,filtOrder,2*D);
    g=eta*squeeze(alpha(n1,:,:,1:2*D))-grad;
    alpha(n1,:,:,1:2*D)= SoftThresold(g,2*D,noOfNodes,eta,lambda,n1,filtOrder);
end
end
figure
Psuedo_Adj=(sum(alpha.^2,4)).^0.5
B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
imagesc(B/max(B(:)))
function [alpha_temp] =SoftThresold(g,t,noOfNodes,eta,lamda,n1,order)
alpha_temp=zeros(size(g));
for n2=1:noOfNodes
    for m=1:order
        alpha_temp(n2,m,:)=g(n2,m,:)/eta;
        alpha_temp(n2,m,:)=g(n2,m,:)/eta*max(0,(1-lamda/norm(reshape(g(n2,m,:),1,t))));
    end
end
end