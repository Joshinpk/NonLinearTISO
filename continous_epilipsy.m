clear all
data1=cell2mat(struct2cell(load('X_chb01_03.mat','m_X_pre_sizure_data')));
data2=cell2mat(struct2cell(load('X_chb01_03.mat','m_X_sizure_data')));
data3=cell2mat(struct2cell(load('X_chb01_03.mat','m_X_post_sizure_data')));
M_x1=resample(data1,2,5);
M_x2=resample(data2,2,5);
M_x3=resample(data3,2,5);
m_X1=M_x1(end-4000:end,:);
m_X3=M_x3(1:4000,:);
m_X2=M_x2;
M_x=[m_X1;m_X2;m_X3];
m_X=M_x(:,1:23)';
[n_row_mx,nTimeInstants]=size(m_X);
m_X=normalize(m_X,2);
% nTimeInstants=1500
%% RF_NLtiso
RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = 23;
RFObj.filtOrder = 2; % we can try a higher order later
RFObj.lambda    = 70/500000;
RFObj.NoOfRF    =20;
RFObj.vsigma    =3*ones(23,1);
RFObj.forgettingFactor=.98;
RFObj.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj.eta       =500;
flag=NaN(nTimeInstants,1);
%initialize
RFState_in = RFObj.initialize(1, m_X( :,1:RFObj.filtOrder)');
for t = RFObj.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in = RFObj.update(RFState_in, mtemp);
    RF_m_predic(:,:)=RFState_in.predictManyFromBuffer(1)';
    RF_m_prediction(1:23,t)= RF_m_predic(:,1);
     Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/Brain.png');
        theta = linspace(0,4*pi,200);
%        image('CData',img,'XData',[-2 9],'YData',[-1 7])
        hold on
        [~, order] = sort(B(:), 'descend');
        n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes
%         threshold = B(order(RFObj.noOfNodes));
        %         threshold = B(order(100));
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
        binaryad=AdjacencyMatrix;
        binaryad(AdjacencyMatrix~=0)=1;
         nnodes(t)=nnz(AdjacencyMatrix)/(23*23);
         degree_count(:,t)=sum(binaryad,2);
         if t>10
        flag(t)=nnodes(t)/mean(nnodes(t-1));
         end
    if t==3000
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain.png');
        theta = linspace(0,4*pi,200);
        subplot(2,3,1)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain.png')
theta = linspace(0,4*pi,200);
image('CData',img,'XData',[-2 9],'YData',[7 -1])
hold on
 [~, order] = sort(B(:), 'descend');
%         n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes 
%         threshold = B(order(RFObj.noOfNodes));
          
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
        
nodes_labels={'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11',...
            'C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
         % to add lables to the nodes
        x=[ 1.7,0.5,1.1,1.5,2.2,2.0,2.0,2.5,4.8,5.0,5.0,4.7,5.5,6.4,6.2,5.5,3.6,3.6,0.2,0.0,4.0,7.,6.8]
        y=[ 4.5,3.4,1.7,0.7,4.5,3.2,2.2,0.9,4.3,3.3,2.2,0.9,4.5,3.5,1.7,0.8,3.3,1.8,2.0,2.9,-.6,3.0,1.7]
        
        plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',16,'LineWidth',6)
        title('P1:Before sizure')
        set(gca,'XTick',[], 'YTick', [])
    end
    if t==4990
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain.png');
        theta = linspace(0,4*pi,200);
        subplot(2,3,2)
         img =imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain.png');
theta = linspace(0,4*pi,200);
image('CData',img,'XData',[-2 9],'YData',[7 -1])
hold on
 [~, order] = sort(B(:), 'descend');
%         n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes 
%         threshold = B(order(RFObj.noOfNodes));
          
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
        
nodes_labels={'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11',...
            'C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
         x=[ 1.7,0.5,1.1,1.5,2.2,2.0,2.0,2.5,4.8,5.0,5.0,4.7,5.5,6.4,6.2,5.5,3.6,3.6,0.2,0.0,4.0,7.1,6.8]
        y=[ 4.5,3.4,1.7,0.7,4.5,3.2,2.2,0.9,4.3,3.3,2.2,0.9,4.5,3.5,1.7,0.8,3.3,1.8,2.0,2.9,-.6,3.0,1.7]
        
        plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',16,'LineWidth',6)
         title('P1:During sizure')
         set(gca,'XTick',[], 'YTick', [])
    end
    if t==11000
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        subplot(2,3,3)
       img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain.png');
theta = linspace(0,4*pi,200);
image('CData',img,'XData',[-2 9],'YData',[7 -1])
hold on
 [~, order] = sort(B(:), 'descend');
%         n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes 
%         threshold = B(order(RFObj.noOfNodes));
          
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
        
nodes_labels={'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11',...
            'C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
          x=[ 1.7,0.5,1.1,1.5,2.2,2.0,2.0,2.5,4.8,5.0,5.0,4.7,5.5,6.4,6.2,5.5,3.6,3.6,0.2,0.0,4.0,7.1,6.8]
        y=[ 4.5,3.4,1.7,0.7,4.5,3.2,2.2,0.9,4.3,3.3,2.2,0.9,4.5,3.5,1.7,0.8,3.3,1.8,2.0,2.9,-.6,3.0,1.7]
        
        plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',16,'LineWidth',6)
         title('P1:After sizure')
         set(gca,'XTick',[], 'YTick', [])
    end
    if t==1000000
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        subplot(2,2,4)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain.png');
theta = linspace(0,4*pi,200);
image('CData',img,'XData',[-2 9],'YData',[7 -1])
hold on
 [~, order] = sort(B(:), 'descend');
%         n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes 
%         threshold = B(order(RFObj.noOfNodes));
%         AdjacencyMatrix=(B>threshold).*B;
          
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
        
nodes_labels={'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11',...
            'C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
         x=[ 1.7,0.5,1.1,1.5,2.2,2.0,2.0,2.5,4.8,5.0,5.0,4.7,5.5,6.4,6.2,5.5,3.6,3.6,0.2,0.0,4.0,7.1,6.8]
        y=[ 4.5,3.4,1.7,0.7,4.5,3.2,2.2,0.9,4.3,3.3,2.2,0.9,4.5,3.5,1.7,0.8,3.3,1.8,2.0,2.9,-.6,3.0,1.7]
        
        plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',16,'LineWidth',6)
    end
    if t==110000
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain_plain.png');
        theta = linspace(0,4*pi,200);
%         image('CData',img,'XData',[-2 9],'YData',[-1 7])
        hold on
        [~, order] = sort(B(:), 'descend');
%         n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes
%         threshold = B(order(RFObj.noOfNodes));
%                 
%         AdjacencyMatrix=(B>threshold).*B;
  
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
        
%         threshold = B(order(RFObj.noOfNodes));
        nodes_labels={'FP1:F7','F7:T7','T7:P7','P7:O1','FP1:F3','F3:C3','C3:P3','P3:O1','FP2:F4','F4:C4','C4:P4',...
            'P4:O2','FP2:F8','F8:T8','T8:P8','P8:O2','FZ:CZ','CZ:PZ','P7:T7','T7:FT9','FT9:FT10','FT10:T8','2T8:P8'};
        G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes
        x=[ 2,3,4,1,2,3,4,5,0,1,2,3,4,5,6,1,2,3,4,5,2,3,4]
        y=[ 1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5]
        subplot(2,2,4)
        imagesc(B)
        %         plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
    end
end
