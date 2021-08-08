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
lamda=10/10;
edge_probability=.1;

for kk=1:1
    A=[];
    for i=1:order
        temp1=10*rand(noOfNodes,noOfNodes).*(rand(noOfNodes,noOfNodes)<edge_probability);
        A=[A;temp1];
    end
        A=zeros(10,5);
        [len_A,bred_A]=size(A)
        ad_in=randi([1,len_A*bred_A],1,1);
        A(ad_in)=10*rand;
        ad_in=randi([1,len_A*bred_A],1,1);
        A(ad_in)=10*rand;
    A1=A/max(max(A));
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
    T_Beta=5*rand(order*noOfNodes,noOfNodes,T);
    
    M_x=zeros(noOfNodes ,noOfNodes );
    M_x(:,1)=y_initial+.1*randn;
    for ii=2:T+order
        temp=A(1:noOfNodes ,:)*y_initial;
        M_x(:,ii)=temp+.1*rand;
        y_initial= M_x(:,ii);
    end
    
    for ii=T+order+1:noOfObservations
        if ii==1000
            list_nz=find(A);
            len=length(list_nz);
            rem_in=randi([1,len],1,1);
            rem=list_nz(rem_in);
            A(rem)=0;
            [len_A,bred_A]=size(A)
            ad_in=randi([1,len_A*bred_A],1,1)
            A(ad_in)=10*rand;
            %                 A=zeros(10,5);
            %         [len_A,bred_A]=size(A)
            %                 ad_in=randi([1,len_A*bred_A],1,1)
            %                 A(ad_in)=8;
            A2=A/max(max(A));
            subplot(4,3,2)
            imagesc(A/max(A(:)))
            title('True dependencies change at t=1000')
            colorbar;
            
        end
        if ii==2000
            list_nz=find(A);
            len=length(list_nz);
            rem_in=randi([1,len],1,1);
            rem=list_nz(rem_in);
            A(rem)=0;
            [len_A,bred_A]=size(A)
            ad_in=randi([1,len_A*bred_A],1,1)
            A(ad_in)=10*rand;
            %                 A=zeros(10,5);
            %         [len_A,bred_A]=size(A)
            %                 ad_in=randi([1,len_A*bred_A],1,1)
            %                 A(ad_in)=6;
            A3=A/max(max(A));
            subplot(4,3,3)
            imagesc(A/max(A(:)))
            title('True  dependencies change at t=2000')
            colorbar;
            
        end
        
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
    %  m_X=normalize(m_X,2);
    Ksigma=10
    vsigma=[10 10 10 10 10];
    for i=1:noOfNodes
        v(i,:,:)=vsigma(i)*vsigma(i)*randn(1,filtOrder,D);
    end
    Gscale=(1/sqrt(2*pi)/Ksigma);
    alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
    Kernal=zeros(noOfNodes,filtOrder,noOfObservations);
    [pred_RF1,alpha_RF1,mse_RF1(kk,:),pmdRF1(kk,:),pfaRF1(kk,:)]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,40,A1,A2,A3);
     [pred_RF2,alpha_RF2,mse_RF2(kk,:),pmdRF2(kk,:),pfaRF2(kk,:)]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,30,A1,A2,A3);
    [pred_RF3,alpha_RF3,mse_RF3(kk,:),pmdRF3(kk,:),pfaRF3(kk,:)]=NLTISORFF(order,lamda,gamma_RF,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,100,A1,A2,A3);
        [pred_NLTISO,alpha_NLTISO,mse_NL(kk,:),pmd1(kk,:),pfa1(kk,:)]=NLTISOF(order,lamda,gamma_NLTISO,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D,A1,A2,A3);
%     
    %%%%Moving average
    for t=12:noOfObservations
        movavg_pred(:,t)=sum(m_X(:,t-11:t-1),2)/10;
    end
    
end
% figure
% plot(sum(mse_RF(:,300:end))/10,'LineWidth',2)
% hold on
% plot(sum(mse_NL(:,300:end))/10,'LineWidth',2)
% hold on
% plot( sum(mse_t(:,300:end))/10,'LineWidth',2)
% lgd1=legend('RF-NLTISO','NLTISO','TIRSO')
% set(lgd1,'fontsize',18)
% ylable('MSE','FontSize',18)
% xlabel('t','FontSize',18)
figure
plot((sum(pmd1))/10,'LineWidth',2);
hold on
plot((sum(pmdRF1))/10,'LineWidth',2);
hold on
plot((sum(pmdRF2))/10,'LineWidth',2);
hold on
plot((sum(pmdRF3))/10,'LineWidth',2);
title('Probability of miss detection ')
lgd2=legend('NL-TISO','RF-TISO-D=10','RF-TISO-D=30','RF-TISO-D=100')
set(lgd2,'fontsize',18)
ylabel('P_{MD}','FontSize',18)
xlabel('t','FontSize',18)

 figure
plot(sum(pfa1)/10,'LineWidth',2);
hold on
 plot(sum(pfaRF1)/10,'LineWidth',2);
hold on
plot(sum(pfaRF2)/10,'LineWidth',2);
hold on
plot(sum(pfaRF3)/10,'LineWidth',2);
title('Probability of false alarm detection ')
lgd3=legend('NL-TISO','RF-TISO-D=10','RF-TISO-D=30','RF-TISO-D=100')
set(lgd3,'fontsize',18)
ylabel('P_{FA}','FontSize',18)
xlabel('t','FontSize',18)


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
function [predt,alpha,mse,pmd,pfa]= NLTISOF(order,lamda,gamma,Ksigma,m_X,noOfNodes,filtOrder,noOfObservations,D,A1,A2,A3)
Gscale=(1/sqrt(2*pi)/Ksigma);
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal=zeros(noOfNodes,filtOrder,noOfObservations);

eta=1/gamma;
%% new code
for t=filtOrder+1:noOfObservations
    if t<1000
        Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        C=B/max(B(:));
        
        if t==990
            
            subplot(4,3,7)
            imagesc(B/max(B(:)))
            title('NLTISO-Estimate at t=990')
            colorbar;
            
        end
        countmd=0;
        MF=wthresh(C-A1,'s',.5);
         MF(isnan(MF))=0;
        
           indmd=find(A1);
           MF(indmd)=0;
           pfa(t)=(sum(sum(MF)))/50;
        
        cv=wthresh(1-abs(A1(indmd)-C(indmd)),'s',.4);
        cv(isnan(cv))=0;
        cunt=nnz(cv);
        tcount=nnz(indmd);
        
       pmd(t)=(tcount-cunt)/ tcount;
    elseif t<2000
        Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        C=B/max(B(:));
        if t==1990
            
            
            subplot(4,3,8)
            imagesc(B/max(B(:)))
            title('NLTISO-Estimate at t=1990');
            colorbar;
            
        end
        countmd=0;
        MF=wthresh(C-A2,'s',.5);
        MF(isnan(MF))=0;
       
          indmd=find(A2);
          MF(indmd)=0;
           pfa(t)=(sum(sum(MF)))/50;
        
        cv=wthresh(1-abs(A2(indmd)-C(indmd)),'s',.4);
         cv(isnan(cv))=0;
        cunt=nnz(cv);
        tcount=nnz(indmd);
        
       pmd(t)=(tcount-cunt)/ tcount;
    else
        Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        C=B/max(B(:));
        if t==2990
            
            
            subplot(4,3,9)
            imagesc(B/max(B(:)))
            title('NLTISO-Estimate at t=2000');
            colorbar;
            
        end
        countmd=0;
        MF=wthresh(C-A3,'s',.5);
        MF(isnan(MF))=0;
       
         indmd=find(A3);
         MF(indmd)=0;
         pfa(t)=(sum(sum(MF)))/50;
        cv=wthresh(1-abs(A3(indmd)-C(indmd)),'s',.4);
         cv(isnan(cv))=0;
        cunt=nnz(cv);
        tcount=nnz(indmd);
        
       pmd(t)=(tcount-cunt)/ tcount;
    end
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
function [predt,alpha,mse,pmd,pfa]= NLTISORFF(order,lamda,gamma,vsigma,m_X,noOfNodes,filtOrder,noOfObservations,D,A1,A2,A3)
Kernal=zeros(noOfNodes,filtOrder,2*D);

for i=1:noOfNodes
    v(i,:,:)=vsigma(i)*vsigma(i)*randn(1,filtOrder,D);
end

eta=1/gamma;

alpha=zeros(noOfNodes,noOfNodes,filtOrder,2*D);
tic
%%RF-TISO
for t=filtOrder+1:noOfObservations
    if t<1000
        Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        C=B/max(B(:));
        
        if t==990
            
            subplot(4,3,7)
            imagesc(B/max(B(:)))
            title('RFTISO-Estimate at t=990')
            colorbar;
            
        end
        countmd=0;
        indmd=find(A1);
        MF=wthresh((C-A1),'s',.5);
        MF(isnan(MF))=0;
        
        indmd=find(A1);
        MF(indmd)=0;
        pfa(t)=(sum(sum(MF)))/50;
        cv=wthresh(1-abs(A1(indmd)-C(indmd)),'s',.4);
         cv(isnan(cv))=0;
        cunt=nnz(cv);
        tcount=nnz(indmd);
        
       pmd(t)=(tcount-cunt)/ tcount;
        
    elseif t<2000
        Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        C=B/max(B(:));
        if t==1990
            
            
            subplot(4,3,8)
            imagesc(B/max(B(:)))
            title('RFTISO-Estimate at t=1990');
            colorbar;
            
        end
        countmd=0;
        indmd=find(A2);
        MF=wthresh(C-A1,'s',.5);
        MF(isnan(MF))=0;
       
        
        indmd=find(A1);
        MF(indmd)=0;
         pfa(t)=(sum(sum(MF)))/50;
        
        cv=wthresh(1-abs(A1(indmd)-C(indmd)),'s',.4);
         cv(isnan(cv))=0;
        cunt=nnz(cv);
        tcount=nnz(indmd);
        
        pmd(t)=(tcount-cunt)/ tcount;
    else
        Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,1);Psuedo_Adj(:,:,2) ];
        C=B/max(B(:));
        if t==2990
            
            
            subplot(4,3,9)
            imagesc(B/max(B(:)))
            title('RFTISO-Estimate at t=2000');
            colorbar;
            
        end
        countmd=0;
        indmd=find(A3);
        MF=wthresh(C-A1,'s',.49);
        MF(isnan(MF))=0;
        
        
        indmd=find(A1);
        MF(indmd)=0;
        pfa(t)=(sum(sum(MF)))/50;
        
        cv=wthresh(1-abs(A1(indmd)-C(indmd)),'s',.4);
         cv(isnan(cv))=0;
        cunt=nnz(cv);
        tcount=nnz(indmd);
        
        pmd(t)=(tcount-cunt)/ tcount;
    end
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