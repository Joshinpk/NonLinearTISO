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
filtOrder        = 4;
T=1;
noOfObservations = 2000
D=30;
order=4;
gamma_RF=1e-3;
gamma_NLTISO=1e-1;
lamda1=10/10;
lamda2=100/10;
edge_probability=.1;


    A=[];
    for i=1:order
        temp1=.1*rand(noOfNodes,noOfNodes).*(rand(noOfNodes,noOfNodes)<edge_probability);
        A=[A;temp1];
    end
    if all(A(:)==0)==1
            A=zeros(10,5);
        [len_A,bred_A]=size(A)
        ad_in=randi([1,len_A*bred_A],1,1);
         A(ad_in)=rand;
        ad_in=randi([1,len_A*bred_A],1,1);
        A(ad_in)=rand;
         ad_in=randi([1,len_A*bred_A],1,1);
        A(ad_in)=rand;
    end

%         A=zeros(10,5);
%         [len_A,bred_A]=size(A)
%         ad_in=randi([1,len_A*bred_A],1,1);
%          A(ad_in)=10+5*randn;
%         ad_in=randi([1,len_A*bred_A],1,1);
%         A(ad_in)=10+5*randn;
%          ad_in=randi([1,len_A*bred_A],1,1);
%         A(ad_in)=10+5*randn;
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
    T=20
    T_Beta=50*randn(order*noOfNodes,noOfNodes,noOfObservations);
    
    M_x=zeros(noOfNodes ,noOfObservations );
    M_x_optimal=zeros(noOfNodes ,noOfObservations );
    error_optimal=zeros(noOfNodes ,noOfObservations );
    M_x(:,1)=y_initial+.1*randn;
    for ii=2:order+T
        temp=A(1:noOfNodes ,:)*y_initial;
        M_x(:,ii)=temp+.1*rand;
        if ii==2990
         M_x(:,ii)=0;
        end
        y_initial= M_x(:,ii);
    end
    
    for ii=order+1+T:noOfObservations
        
        
        for k=1:noOfNodes
            temp=0;
            for m=1:order
                for l=1:noOfNodes
                    temp=temp+ A(noOfNodes *(m-1)+k,l)*((1*gscale*exp(-1000*(M_x(l,1:ii-m-1)-M_x(l,ii-m)).^2/(10*gsigma*gsigma))*reshape(T_Beta(noOfNodes *(m-1)+k,l,1:ii-m-1),[ii-m-1,1])));
                    
                end
            end
            M_x_optimal(k,ii)=temp;
            M_x(k,ii)=temp+.5*randn;
            error_optimal(k,ii)=M_x(k,ii)-M_x_optimal(k,ii);
            
        end
        
    end
    m_X=M_x;
    [n_row_mx,nTimeInstants]=size(m_X);
RFObj_tirso = RF_nltirso; % set tirso object up
RFObj_tirso.noOfNodes = 5;
RFObj_tirso.filtOrder = 4; % we can try a higher order later
RFObj_tirso.lambda    = 50/500000;
RFObj_tirso.NoOfRF    =30;
RFObj_tirso.vsigma    =1;
  RFObj_tirso.forgettingFactor=.98;
 RFObj_tirso.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj_tirso.eta       =1000
RFState_in_tirso = RFObj_tirso.initialize( 10,m_X( :,1:RFObj_tirso.filtOrder)');
 for t = RFObj_tirso.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in_tirso = RFObj_tirso.update(RFState_in_tirso, mtemp);
    if t==2000
        fffff=4;
    end
   %temp_eig(t)= RFState_in_tirso.eig_value_phi;
 end
%  Psuedo_Adj_real=(sum(adjacency_real(:,:,:,1:T).^2,4)).^0.5;
  Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5
        B=[Psuedo_Adj(:,:,4);Psuedo_Adj(:,:,3) ];
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
RFObj_tirso_2=RFObj_tirso;
 RFObj_tirso_2.forgettingFactor=.001;
 RFState_in_tirso_2 = RFObj_tirso_2.initialize( 10,m_X( :,1:RFObj_tirso.filtOrder)');

  for t = RFObj_tirso.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in_tirso_2 = RFObj_tirso_2.update(RFState_in_tirso_2, mtemp);
    if t==2000
        fffff=4;
    end
   %temp_eig(t)= RFState_in_tirso.eig_value_phi;
  end
  Psuedo_Adj=(sum(RFState_in_tirso_2.coeff.^2,4)).^0.5
        B=[Psuedo_Adj(:,:,4);Psuedo_Adj(:,:,3) ];
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