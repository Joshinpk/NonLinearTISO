clear all
close all
data1=cell2mat(struct2cell(load('X_chb07_12.mat','m_X_pre_sizure_data')));
data2=cell2mat(struct2cell(load('X_chb07_12.mat','m_X_sizure_data')));
data3=cell2mat(struct2cell(load('X_chb07_12.mat','m_X_post_sizure_data')));
M_x1=resample(data1,2,5);
M_x2=resample(data2,2,5);
M_x3=resample(data3,2,5);
m_X1=M_x1(end-4000:end,:);
m_X3=M_x3(1:40,:);
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
RFObj.lambda    =70/80
RFObj.NoOfRF    =10;
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
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain_plain.png');
        theta = linspace(0,4*pi,200);
%        image('CData',img,'XData',[-2 9],'YData',[-1 7])
        hold on
        [~, order] = sort(B(:), 'descend');
        n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes
%         threshold = B(order(RFObj.noOfNodes));
        %         threshold = B(order(100));
        threshold = max(B(:))/3;
        AdjacencyMatrix=(B>threshold).*B;
        binaryad=AdjacencyMatrix;
        binaryad(AdjacencyMatrix~=0)=1;
         nnodes(t)=nnz(AdjacencyMatrix)/(23*23);
         degree_count(:,t)=sum(binaryad,2);
         if t>10
        flag(t)=nnodes(t)/mean(nnodes(t-1));
         end
    if t==2500
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain_plain.png');
        theta = linspace(0,4*pi,200);
        subplot(2,2,1)
        % image('CData',img,'XData',[-2 9],'YData',[-1 7])
        hold on
        [~, order] = sort(B(:), 'descend');
        n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes
%         threshold = B(order(RFObj.noOfNodes));
        %         threshold = B(order(100));
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
       
        nodes_labels={'FP1:F7','F7:T7','T7:P7','P7:O1','FP1:F3','F3:C3','C3:P3','P3:O1','FP2:F4','F4:C4','C4:P4',...
            'P4:O2','FP2:F8','F8:T8','T8:P8','P8:O2','FZ:CZ','CZ:PZ','P7:T7','T7:FT9','FT9:FT10','FT10:T8','2T8:P8'};
        G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes
        x=[ 2,3,4,1,2,3,4,5,0,1,2,3,4,5,6,1,2,3,4,5,2,3,4]
        y=[ 1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5]
        
        imagesc(B)
        
%                 plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
    end
    if t==4900
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain_plain.png');
        theta = linspace(0,4*pi,200);
        subplot(2,2,2)
       % image('CData',img,'XData',[-2 9],'YData',[-1 7])
        hold on
        [~, order] = sort(B(:), 'descend');
        n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes
%         threshold = B(order(RFObj.noOfNodes));
        %         threshold = B(order(100));
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
        
        
        nodes_labels={'FP1:F7','F7:T7','T7:P7','P7:O1','FP1:F3','F3:C3','C3:P3','P3:O1','FP2:F4','F4:C4','C4:P4',...
            'P4:O2','FP2:F8','F8:T8','T8:P8','P8:O2','FZ:CZ','CZ:PZ','P7:T7','T7:FT9','FT9:FT10','FT10:T8','2T8:P8'};
        G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes
        x=[ 2,3,4,1,2,3,4,5,0,1,2,3,4,5,6,1,2,3,4,5,2,3,4]
        y=[ 1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5]
        
        imagesc(B)
            %   plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
    end
    if t==6600
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain_plain.png');
        theta = linspace(0,4*pi,200);
        subplot(2,2,3)
       % image('CData',img,'XData',[-2 9],'YData',[-1 7])
        hold on
        [~, order] = sort(B(:), 'descend');
        n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes
%         threshold = B(order(RFObj.noOfNodes));
        %         threshold = B(order(100));
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
       
        nodes_labels={'FP1:F7','F7:T7','T7:P7','P7:O1','FP1:F3','F3:C3','C3:P3','P3:O1','FP2:F4','F4:C4','C4:P4',...
            'P4:O2','FP2:F8','F8:T8','T8:P8','P8:O2','FZ:CZ','CZ:PZ','P7:T7','T7:FT9','FT9:FT10','FT10:T8','2T8:P8'};
        G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes
        x=[ 2,3,4,1,2,3,4,5,0,1,2,3,4,5,6,1,2,3,4,5,2,3,4]
        y=[ 1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5]
        
         imagesc(B)
                 %plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
    end
    if t==10000
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain_plain.png');
        theta = linspace(0,4*pi,200);
         subplot(2,2,4)
        %image('CData',img,'XData',[-2 9],'YData',[-1 7])
        hold on
        [~, order] = sort(B(:), 'descend');
        n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes
%         threshold = B(order(RFObj.noOfNodes));
        %         threshold = B(order(100));
        threshold = max(B(:))/2;
        AdjacencyMatrix=(B>threshold).*B;
        
        nodes_labels={'FP1:F7','F7:T7','T7:P7','P7:O1','FP1:F3','F3:C3','C3:P3','P3:O1','FP2:F4','F4:C4','C4:P4',...
            'P4:O2','FP2:F8','F8:T8','T8:P8','P8:O2','FZ:CZ','CZ:PZ','P7:T7','T7:FT9','FT9:FT10','FT10:T8','2T8:P8'};
        G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes
        x=[ 2,3,4,1,2,3,4,5,0,1,2,3,4,5,6,1,2,3,4,5,2,3,4]
        y=[ 1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5]
       
        imagesc(B)
 %               plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
    end
    if t==110000
        Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
        B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
        %imagesc(B)
        img = imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain.png');
        theta = linspace(0,4*pi,200);
%         image('CData',img,'XData',[-2 9],'YData',[-1 7])
        hold on
        [~, order] = sort(B(:), 'descend');
        n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes
        threshold = B(order(RFObj.noOfNodes));
                
        AdjacencyMatrix=(B>threshold).*B;
%         threshold = B(order(RFObj.noOfNodes));
        nodes_labels={'FP1:F7','F7:T7','T7:P7','P7:O1','FP1:F3','F3:C3','C3:P3','P3:O1','FP2:F4','F4:C4','C4:P4',...
            'P4:O2','FP2:F8','F8:T8','T8:P8','P8:O2','FZ:CZ','CZ:PZ','P7:T7','T7:FT9','FT9:FT10','FT10:T8','2T8:P8'};
        G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes
        x=[ 4.5,0.3,1.2,1.5,2.0,2.0,2.0,2.5,5.0,5.0,5.0,4.7,5.5,6.2,6.2,5.5,3.5,3.5,0.5,0.3,4.0,7.0,6.7]
        y=[ 1.7,3.4,1.5,0.5,4.5,3.2,2.0,0.7,4.5,3.5,2.5,0.7,4.5,3.5,1.7,0.5,3.5,1.5,2.0,2.9,-.5,3.2,1.7]
        subplot(2,2,4)
        imagesc(B)
        %         plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',10,'LineWidth',2)
    end
end
figure
t=degree_count([3,15,19,20,22,23],:)/6;
f=degree_count([1,5,6,9,10,13,17],:)/7;
o=degree_count([4,12,16],:)/3;
plot(sum(t))
hold on
plot(sum(f))
plot(sum(o))