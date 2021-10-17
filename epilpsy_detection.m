clear all
data=cell2mat(struct2cell(load('X_chb01_03.mat','m_X_post_sizure_data')));
M_x=resample(data,2,5);
   m_X=M_x(end-1500:end,:)';
%      m_X=M_x';
[n_row_mx,nTimeInstants]=size(m_X);
m_X=normalize(m_X,2)
nTimeInstants=1500
%% RF_NLtiso
RFObj = RF_nltirso; % set tirso object up
RFObj.noOfNodes = 23;
RFObj.filtOrder = 2; % we can try a higher order later
RFObj.lambda    = 70/5000;
RFObj.NoOfRF    =5;
RFObj.vsigma    =3;
RFObj.forgettingFactor=.98;
RFObj.h_stepsize= @(RF_ts)1/eigs(RF_ts.m_Phi,1);
RFObj.eta       =500;
%initialize
RFState_in = RFObj.initialize(1, m_X( :,1:RFObj.filtOrder)');
 for t = RFObj.filtOrder+1:nTimeInstants
    mtemp= m_X(:, t);
    RFState_in = RFObj.update(RFState_in, mtemp);
    RF_m_predic(:,:)=RFState_in.predictManyFromBuffer(1)';
    RF_m_prediction(1:23,t)= RF_m_predic(:,1);
      
 end
figure
Psuedo_Adj=(sum(RFState_in.coeff.^2,4)).^0.5;
B=Psuedo_Adj(:,:,2)-diag(diag(Psuedo_Adj(:,:,2)));
%imagesc(B) 
img = imrotate(imread('/Users/rohanmoney/git_rohan/NonLinearTISO/brain.png'),180);
theta = linspace(0,4*pi,200);
image('CData',img,'XData',[-2 9],'YData',[-1 7])
hold on
 [~, order] = sort(B(:), 'descend');
        n_toPreserve = RFObj.noOfNodes; % displaying the edges in the order of the number of nodes 
        threshold = B(order(RFObj.noOfNodes));
        AdjacencyMatrix=(B>threshold).*B;
nodes_labels={'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11',...
            'C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23'};
		G=digraph( AdjacencyMatrix', nodes_labels); % transposed
        Weights=G.Edges.Weight;
        Names= G.Nodes;
        % to add lables to the nodes 
        x=[ 1.5,0.5,0.3,0.5,2.1,2.1,2.1,.75,1,2,3,4,5,6,1,2,3,4,5,2,3,4,4]
        y=[ 4.5,3.5,1.5,1.5,4.5,3.5,2.0,0.8,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5]
        
        plot(G,'XData',x,'YData',y,'EdgeCData', Weights,'ArrowSize',16,'LineWidth',6)