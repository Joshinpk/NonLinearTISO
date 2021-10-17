clc
clear all

noOfNodes        = 10;
filtOrder        = 4;
T=1;
noOfObservations = 3000
D=30;
order=4;
gamma_RF=1e-3;
gamma_NLTISO=1e-1;
lamda1=10/10;
lamda2=100/10;
edge_probability=.25;


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
figure
    A1=A/max(max(A));
    subplot(4,3,1)
    imagesc(A/max(A(:)))
    tit=title('True dependencies intial')
    set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
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
    M_x(:,1)=y_initial+.15*randn;
    for ii=2:order+T
        temp=A(1:noOfNodes ,:)*y_initial;
        M_x(:,ii)=temp+.15*rand;
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
            M_x(k,ii)=temp+.15*randn;
            error_optimal(k,ii)=M_x(k,ii)-M_x_optimal(k,ii);
            
        end
        if ii==1000
            no_non_zero=nnz(A);
            no_to_change=round(no_non_zero*.3);%30 percentgae of edges change
            non_zero_entries=find(A);
            entries_to_change= randperm(numel(non_zero_entries),no_to_change)
            temp=non_zero_entries(entries_to_change)
            A(temp)=0;
%             for jj=1:no_to_change
%               r = randi([1 noOfNodes*noOfNodes*order],1,1)  
%              
%                   A(r)=.1*rand;
%               A1=A/max(max(A));
    subplot(4,3,2)
    imagesc(A/max(A(:)))
    tit=title('True dependencies intial')
    set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
% 
%             end
        end
         if ii==2000
            no_non_zero=nnz(A);
            no_to_change=round(no_non_zero*.3);%30 percentgae of edges change
            non_zero_entries=find(A);
            entries_to_change= randperm(numel(non_zero_entries),no_to_change)
            temp=non_zero_entries(entries_to_change)
            A(temp)=0;
%             for jj=1:no_to_change
%               r = randi([1 noOfNodes*noOfNodes*order],1,1)  
%              
%                   A(r)=.1*rand;
%              
%  A1=A/max(max(A));
    subplot(4,3,3)
    imagesc(A/max(A(:)))
    tit=title('True dependencies intial')
   set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
%             end
        end
        
    end
    m_X=M_x;
    [n_row_mx,nTimeInstants]=size(m_X);
RFObj_tirso = RF_nltirso; % set tirso object up
RFObj_tirso.noOfNodes = 10;
RFObj_tirso.filtOrder = 4; % we can try a higher order later
RFObj_tirso.lambda    = 50/5000;
RFObj_tirso.NoOfRF    =30;
RFObj_tirso.vsigma    =1;
  RFObj_tirso.forgettingFactor=.98;
 RFObj_tirso.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj_tirso.eta       =500
RFState_in_tirso = RFObj_tirso.initialize( 10,m_X( :,1:RFObj_tirso.filtOrder)');
 for t = RFObj_tirso.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in_tirso = RFObj_tirso.update(RFState_in_tirso, mtemp);
    RF_m_predic_tirso(:,:)=RFState_in_tirso.predictManyFromBuffer(1)';
    RF_m_prediction_tirso(1:10,t)= RF_m_predic_tirso(:,1);
    if t== 990
        Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5
        B=[Psuedo_Adj(:,:,4);Psuedo_Adj(:,:,3); Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1)];
       subplot(4,3,4)
        imagesc(B/max(B(:)))
        tit=title('RFNL-TIRSO-Estimate at t= 990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
    if t== 1990
        Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,4);Psuedo_Adj(:,:,3); Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1)];
        subplot(4,3,5)
        imagesc(B/max(B(:)))
        tit=title('RFNL-TIRSO-Estimate at t= 1990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
    if t== 2990
        Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,4);Psuedo_Adj(:,:,3); Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1)];
        subplot(4,3,6)
        imagesc(B/max(B(:)))
        tit=title('RFNL-TIRSO-Estimate at t= 2990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
 end
%  Psuedo_Adj_real=(sum(adjacency_real(:,:,:,1:T).^2,4)).^0.5;
  Psuedo_Adj=(sum(RFState_in_tirso.coeff.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,4);Psuedo_Adj(:,:,3);Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1) ];
       
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

RFObj_tiso = RF_nltiso_test; % set tirso object up
RFObj_tiso.noOfNodes = 10;
RFObj_tiso.filtOrder = 4; % we can try a higher order later
RFObj_tiso.lambda    = 50/500000;
RFObj_tiso.NoOfRF    =30;
RFObj_tiso.vsigma    =1;
RFObj_tiso.eta       =1000

 RFState_in_tiso = RFObj_tiso.initialize( m_X( :,1:RFObj_tirso.filtOrder)');

  for t = RFObj_tiso.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in_tiso = RFObj_tiso.update(RFState_in_tiso, mtemp);
     RF_m_predic_tiso(:,:)=RFState_in_tiso.predictManyFromBuffer(1)';
    RF_m_prediction_tiso(1:10,t)= RF_m_predic_tiso(:,1);
   %temp_eig(t)= RFState_in_tirso.eig_value_phi;
    if t== 990
        Psuedo_Adj=(sum(RFState_in_tiso.coeff.^2,4)).^0.5
        B=[Psuedo_Adj(:,:,4);Psuedo_Adj(:,:,3); Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1)];
       subplot(4,3,7)
        imagesc(B/max(B(:)))
        tit=title('RFNL-TISO-Estimate at t= 990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
    if t== 1990
        Psuedo_Adj=(sum(RFState_in_tiso.coeff.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,4);Psuedo_Adj(:,:,3); Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1)];
        subplot(4,3,8)
        imagesc(B/max(B(:)))
        tit=title('RFNL-TISO-Estimate at t= 1990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
    if t== 2990
        Psuedo_Adj=(sum(RFState_in_tiso.coeff.^2,4)).^0.5;
        B=[Psuedo_Adj(:,:,4);Psuedo_Adj(:,:,3); Psuedo_Adj(:,:,2);Psuedo_Adj(:,:,1)];
        subplot(4,3,9)
        imagesc(B/max(B(:)))
        tit=title('RFNL-TISO-Estimate at t= 2990');
%         colorbar;
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
  end
  %% TIRSO
addpath '/Users/rohanmoney/git_rohan/STInference/gsim/ProcessingBlocks/STFunctionEstimators/TirsoObjects/'


[n_row_mx,nTimeInstants]=size(m_X);
tirsoObj = Tirso; % set tirso object up
tirsoObj.noOfNodes = 10;
tirsoObj.order     = 4; % we can try a higher order later
tirsoObj.regPar    = 1e-2;
tirsoObj.b_shrinkSelfLoops  = 1; % Bolstad
tirsoObj.forgettingFactor   = 0.98;
tirsoObj.h_stepsize         = @(ts)1/eigs(ts.m_Phi,1);
% initialize
tState_in = tirsoObj.initialize(0, m_X( :,1:tirsoObj.order)');
e_tirso=zeros(24,noOfObservations);
for t = tirsoObj.order+1:nTimeInstants
    tic
    mtemp= m_X(:, t);
    tState_in = tirsoObj.update(tState_in, mtemp);
    m_predic(:,:)=tState_in.predictManyFromBuffer(1)';
    m_prediction(1:10,t)= m_predic(:,1);
    if t==990
         subplot(4,3,10)
        imagesc(tState_in.m_A/max(tState_in.m_A(:)))
%         colorbar;
       tit= title('TIRSO estimate at t=1990')
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
        if t==1990
         subplot(4,3,11)
        imagesc(tState_in.m_A/max(tState_in.m_A(:)))
%         colorbar;
       tit= title('TIRSO estimate at t=1990')
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
        end
        if t==2990
         subplot(4,3,12)
        imagesc(tState_in.m_A/max(tState_in.m_A(:)))
%         colorbar;
       tit= title('TIRSO estimate at t=1990')
       set(tit,'fontSize',24)
       ax=gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
 c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
set(c,'fontSize',24);
    end
      
end
figure
sensor_to_plot=1;
%  plot(m_X(sensor_to_plot,12:noOfObservations),'LineWidth',2)
%  hold on
% % 
%   plot(RF_m_prediction(sensor_to_plot,:),'LineWidth',2)
%  hold on
%   plot(m_prediction(sensor_to_plot,:),'LineWidth',2)
% %   hold on
%   plot((RF_m_prediction(sensor_to_plot,1:end-12)-m_X(sensor_to_plot,13:end)).^2)
 error_temp=sum((RF_m_prediction_tirso(:,1:end-1)-m_X(:,2:end)).^2,1);
 for ii=1:length(error_temp)
      nmse(ii)=sum(error_temp(1:ii))/sum(sum((m_X(:,2:ii+1)).^2,1)) ;
 end
 hold on
 plot(nmse)
 hold on
  error_temp=sum((RF_m_prediction_tiso(:,1:end-1)-m_X(:,2:end)).^2,1);
 for ii=1:length(error_temp)
     nmse(ii)=sum(error_temp(1:ii))/sum(sum((m_X(:,2:ii+1)).^2,1)) ;
 end
 hold on
 plot(nmse)
  error_temp=sum((m_prediction(:,1:end-1)-m_X(:,2:end)).^2,1);
 for ii=1:length(error_temp)
     nmse(ii)=sum(error_temp(1:ii))/sum(sum((m_X(:,2:ii+1)).^2,1)) ;
 end
 hold on
 plot(nmse)