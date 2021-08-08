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
filtOrder        = 12;
noOfObservations = 1000;
D=20
step=20;
order=filtOrder;
gamma_RF=1e-4;
gamma_NLTISO=1e-3;
lamda=70/500000;
m_X=M_x(1:24,1:noOfObservations);
m_X=normalize(m_X,2)
Ksigma=50;
vsigma=100;
Gscale=(1/sqrt(2*pi)/Ksigma);
 [pred_RF1,alpha_RF1,mse1,ti1,npred1]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,10,step);
 [pred_RF2,alpha_RF2,mse2,mset2,ti2,npred2]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,10,step);
 [pred_RF3,alpha_RF3,mse3,ti3,npred3]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,100,step);
  [pred_NLTISO,alpha_NLTISO,mse_nltiso,ti_nltiso]=NLTISOF(order,lamda,gamma_NLTISO,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D);
%%%%Moving average
for t=12:noOfObservations
    tic
    movavg_pred(:,t)=sum(m_X(:,t-11:t-1),2)/10;
    timov(t)=toc;
    error_mov=(m_X(:,t)- movavg_pred(:,t)).^2;
    
    mse_mov(t)=sum( error_mov)/sum(m_X(:,t).^2);
    
end
%%% TIRSO
addpath '/Users/rohanmoney/git_rohan/STInference/gsim/ProcessingBlocks/STFunctionEstimators/TirsoObjects/'


[n_row_mx,nTimeInstants]=size(m_X);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = 24;
tirsoObj.order     = 12; % we can try a higher order later
tirsoObj.regPar    = 1e-1;
tirsoObj.b_shrinkSelfLoops  = 1; % Bolstad
tirsoObj.forgettingFactor   = 0.9;
tirsoObj.h_stepsize         = @(ts)1/eigs(ts.m_Phi,1);
% initialize
tState_in = tirsoObj.initialize(0, m_X( :,1:tirsoObj.order)');
e_tirso=zeros(24,noOfObservations);
for t = tirsoObj.order+1:nTimeInstants
    tic
    mtemp= m_X(:, t);
    tState_in = tirsoObj.update(tState_in, mtemp);
    m_predic(:,:)=tState_in.predictManyFromBuffer(20)';
    m_prediction(1:24,t)= m_predic(:,20);
      if t<800
    error_tiso=(m_X(:,t+20)-m_prediction(:,t)).^2;
  
        e_tirso(:,t)=error_tiso;
        mset(:,t)=sum(e_tirso(:,t),2)/t;
        mse_tiso(t)=sum(error_tiso)/sum(m_X(:,t).^2);
    end
    titirso(t)=toc;
end
Psuedo_Adj=(sum(alpha_RF2.^2,4)).^0.5
v_avgNorms= mean(Psuedo_Adj(:,:,12), 3, 'omitNaN');

        v_avgNorms(1:(noOfNodes +1):end)=0;  % avoiding selfloops to plot
        [~, order] = sort(v_avgNorms(:), 'descend');
        n_toPreserve = noOfNodes; % displaying the edges in the order of the number of nodes 
        threshold = v_avgNorms(order(n_toPreserve));
        AdjacencyMatrix=(v_avgNorms>threshold).*v_avgNorms;
nodes_labels={'S1 (P)','S2 (P)','S3 (P)','S4 (P)','S5 (P)','S6 (P)','S7 (L)','S8 (L)','S9 (L)','S10(L)','S11(T)',...
            'S12 (L)','S13 (L)','S14 (L)','S15 (L)','S16 (L)','S17 (P)','S18 (P)','S19 (P)','S20 (T)','S21 (T)','S22 (T)','S23 (T)','S24 (T)'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
        
        plot(G,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
        figure
        Psuedo_Adj=(sum(alpha_RF1.^2,4)).^0.5
v_avgNorms= mean(Psuedo_Adj(:,:,12), 3, 'omitNaN');

        v_avgNorms(1:(noOfNodes +1):end)=0;  % avoiding selfloops to plot
        [~, order] = sort(v_avgNorms(:), 'descend');
        n_toPreserve = noOfNodes; % displaying the edges in the order of the number of nodes 
        threshold = v_avgNorms(order(n_toPreserve));
        AdjacencyMatrix=(v_avgNorms>threshold).*v_avgNorms;
nodes_labels={'S1 (P)','S2 (P)','S3 (P)','S4 (P)','S5 (P)','S6 (P)','S7 (L)','S8 (L)','S9 (L)','S10(L)','S11(T)',...
            'S12 (L)','S13 (L)','S14 (L)','S15 (L)','S16 (L)','S17 (P)','S18 (P)','S19 (P)','S20 (T)','S21 (T)','S22 (T)','S23 (T)','S24 (T)'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
        
        plot(G,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
        figure
        Psuedo_Adj=(sum(alpha_RF3.^2,4)).^0.5
v_avgNorms= mean(Psuedo_Adj(:,:,12), 3, 'omitNaN');

        v_avgNorms(1:(noOfNodes +1):end)=0;  % avoiding selfloops to plot
        [~, order] = sort(v_avgNorms(:), 'descend');
        n_toPreserve = noOfNodes; % displaying the edges in the order of the number of nodes 
        threshold = v_avgNorms(order(n_toPreserve));
        AdjacencyMatrix=(v_avgNorms>threshold).*v_avgNorms;
nodes_labels={'S1 (P)','S2 (P)','S3 (P)','S4 (P)','S5 (P)','S6 (P)','S7 (L)','S8 (L)','S9 (L)','S10(L)','S11(T)',...
            'S12 (L)','S13 (L)','S14 (L)','S15 (L)','S16 (L)','S17 (P)','S18 (P)','S19 (P)','S20 (T)','S21 (T)','S22 (T)','S23 (T)','S24 (T)'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
        
        plot(G,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
        figure
        colorbar
        Psuedo_Adj=(sum(alpha_NLTISO.^2,4)).^0.5
v_avgNorms= mean(Psuedo_Adj(:,:,12), 3, 'omitNaN');

        v_avgNorms(1:(noOfNodes +1):end)=0;  % avoiding selfloops to plot
        [~, order] = sort(v_avgNorms(:), 'descend');
        n_toPreserve = noOfNodes; % displaying the edges in the order of the number of nodes 
        threshold = v_avgNorms(order(n_toPreserve));
        AdjacencyMatrix=(v_avgNorms>threshold).*v_avgNorms;
nodes_labels={'S1 (P)','S2 (P)','S3 (P)','S4 (P)','S5 (P)','S6 (P)','S7 (L)','S8 (L)','S9 (L)','S10(L)','S11(T)',...
            'S12 (L)','S13 (L)','S14 (L)','S15 (L)','S16 (L)','S17 (P)','S18 (P)','S19 (P)','S20 (T)','S21 (T)','S22 (T)','S23 (T)','S24 (T)'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
        
        plot(G,'layout','circle','EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
%          axes('FontSize',18)
% % Full_Ker_time=toc;
% % % Psuedo_Adj=(sum(alpha_RF.^2,4)).^0.5;
% plot(m_X(8,20:noOfObservations),'LineWidth',2)
% hold on
% % plot(npred1(10,:),'LineWidth',2)
% hold on
% plot(npred2(8,:),'LineWidth',2)
% hold on
% % plot(npred3(10,:),'LineWidth',2)
% % step1=(npred1(:,1:991)-m_X(8,10:noOfObservations)).^2;
% % step2=npred2(:,1:991)-m_X(8,10:noOfObservations).^2;
% % step3=npred3(:,1:991)-m_X(8,10:noOfObservations).^2;
% % % for ii=1:991
% %     ms1(ii)=sum(step1(:,ii),2);
% %     ms2(ii)=sum(step2(:,ii),2);
% %     ms3(ii)=sum(step3(:,ii),2);
% % end
% % hold on
% % % plot(movavg_pred(7,:))
% % % hold on
%  plot(m_prediction(8,:),'LineWidth',2)
% % legend('Orginal','RF D= 10','RF D= 50','RF D= 100','TIRSO')
%  figure
% % 
% %  plot(mse_nltiso(10,:),'LineWidth',2)
% % hold on
% % plot(mse1,'LineWidth',2)
% % hold on
%  plot(mse2,'LineWidth',2)
% % hold on
% % plot(mse3,'LineWidth',2)
%  hold on
%  plot(mse_tiso,'LineWidth',2)
% 
% 
%  legend('RF D= 200','RF D= 50','RF D= 500','TIRSO')
% figure
% % plot(timov,'LineWidth',2)
% % hold on
% % plot(ti_nltiso,'LineWidth',2)
%
% hold on
% plot(ti1,'LineWidth',2)
% hold on
% plot(ti2,'LineWidth',2)
% hold on
% plot(ti3,'LineWidth',2)
% legend('RF D= 10','RF D= 50','RF D= 100')
% % legend('Moving average','NL-TISO','RF-TISO','TIRSO')
% xlabel('t')
% ylabel('computation time(s)')
% title('Time taken for each prediction')

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
function [predt,alpha,mse,mset,ti]= NLTISOF(order,lamda,gamma,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D)
Gscale=(1/sqrt(2*pi)/Ksigma);
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal=zeros(noOfNodes,filtOrder,noOfObservations);

eta=1/gamma;
%% new code
error=zeros(24,noOfObservations);
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
        alpha(n1,:,:,t)= alpha(n1,:,:,t-1);
        predt(n1,t)=(vec(alpha(n1,:,:,1:t)))'*K_vec;
        grad_vec=(vec(alpha(n1,:,:,1:t)))'*K_vec*K_vec'-K_vec'*m_X(n1,t);
        error_temp(n1)=(predt(n1,t)-m_X(n1,t))^2;
        grad=reshape(grad_vec,noOfNodes,filtOrder,t);
        g=eta*squeeze(alpha(n1,:,:,1:t))-grad;
        
        alpha(n1,:,:,1:t)= SoftThresold(g,t,noOfNodes,eta,lamda,n1,order);
    end
    ti(t)=toc;
    if t>0
        error(:,t)=error_temp;
        mset(:,t)=sum(error(:,t),2);
        mse(t)=sum(error_temp)/sum(m_X(:,t).^2);
    end
    
    (sum(alpha.^2,4)).^0.5;
end
end
function [predt,alpha,mse,mset,ti,npred]= NLTISORFF(order,lamda,gamma,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,D,step)
Kernal=zeros(noOfNodes,filtOrder,2*D);
v=vsigma*vsigma*randn(noOfNodes,filtOrder,D);

eta=1/gamma;

alpha=zeros(noOfNodes,noOfNodes,filtOrder,2*D);
error=zeros(24,noOfObservations);

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
    npred(:,t)=pred_nahead(order,lamda,gamma,vsigma,m_X(:,1:t),noOfNodes,filtOrder,noOfObservations,D,alpha(:,:,:,1:2*D),t,step,v);
    ti(t)=toc;
    if t<800
        error_temp=(npred(:,t)-m_X(:,t+20)).^2;
        error(:,t)=error_temp;
        mset(:,t)=sum(error(:,t),2)/t;
        mse(:,t)=sum(error_temp)/sum(m_X(:,t).^2);
    end
    
end
end
function [out]=pred_nahead(order,lamda,gamma,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,D,alpha,t,step,v)
for ii=1:step
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
        predt(n1)=(vec(alpha(n1,:,:,1:2*D)))'*K_vec;
        m_X(n1,t)=predt(n1);
    end
    t=t+1;
end
out=m_X(:,t-1);
end
function [alpha_temp] =SoftThresold(g,t,noOfNodes,eta,lamda,n1,order)
alpha_temp=zeros(size(g));
for n2=1:noOfNodes
    for m=1:order
%          if n2==n1
%                           alpha_temp(n2,m,:)=g(n2,m,:)/eta;
% % %             
%         else
            %            
            alpha_temp(n2,m,:)=g(n2,m,:)/eta*max(0,(1-lamda/norm(reshape(g(n2,m,:),1,t))));
%           end
    end
end
end