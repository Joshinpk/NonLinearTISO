<<<<<<< HEAD
% Dictionary based NL-TISO is implemented in this code. This implementation
% builds dictionary based on the simalarity measured using kernal
% functions. The length of the dictionary is not fixed in this
% implementation.


clc 
clear all
close all
filename = './data/sytem20.mat';
Data_struct = load(filename,'-mat');
M_x=Data_struct.m_X;
noOfNodes        = 24;
=======
clc; clear all
addpath(genpath('./Execution/'));
initializeGsim.initializePath;
noOfNodes        = 4;%edited by joshin, it was 12
>>>>>>> e36bcfe75a67f56ed31d101d9d2fd84488dea84d
filtOrder        = 2;
noOfObservations = 100; %edited by joshin, it was 1000 before| % 300 when testing plots; 10000 normally
transitionTime   = noOfObservations/3;
transitionPar    = TimeVaryingVARSTFunctionGenerator.stParameter(...
0.99, transitionTime);%0.5E-6;
% determines the speed of the transition(the higher, thefaster)
noOfExperiments  = 1; %Monte Carlo iterations %edited by joshin, it was 10
erdosRenyiEdgeProbability = 0.2;%edited by jpk, was 0.1 before
% 		erdosRenyiEdgeProbability2 = 0.1;
syntheticSigma            = 0.5; %0.02
regPar       = 1e-5; % regPar==1E-4 -> proxRecRDA=0.05, proxRDA=0.005
f_factor     = 0.99;
prox_rdaRec  = 0.01;
prox_rda     = 0.1;
initialSigma = 0.01; %0.01

rng(3);
m_adjacency = (rand(noOfNodes)<erdosRenyiEdgeProbability)...
    .*not(eye(noOfNodes)); %directed graph, noSelfLoops
%     		figure(1); clf; bg = biograph(m_adjacency); bg.view;
myGraph = Graph('m_adjacency',m_adjacency);

myCoefficients1 = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, filtOrder);
myCoefficients2 = VARSTFunctionGenerator.randomCoefficientsFromGraph(myGraph, filtOrder);
gen = TimeVaryingVARSTFunctionGenerator;
gen.N = noOfNodes;
gen.Q = filtOrder;
gen.T = noOfObservations;
gen.sigma = syntheticSigma;
gen.Gsigma=1;

[oldCoefficients, t_A, transitionCurve] = gen.smoothTransitionModel(myCoefficients1,myCoefficients2, ...
    transitionPar, transitionTime, noOfObservations);
reverseCoefficients = gen.smoothTransitionModel(myCoefficients2, myCoefficients1, ...
    transitionPar, transitionTime, noOfObservations);

dum_pre=zeros(size(reverseCoefficients(:,:,1,1)));
dum_pre=dum_pre+0*diag(ones(noOfNodes,1));
dum_pre(1,4,:,:)=1;

% dum_pre(5,8,:,:)=0.5;
% dum_pre(7,1,:,:)=0.5;
% dummyCoe=zeros(size(reverseCoefficients));dummyCoe(2,3,:,:)=0.9;dummyCoe(3,2,:,:)=0.9;
dummyCoe=repmat(dum_pre,1,1,filtOrder,noOfObservations);
% dummyCoe(:,:,2,:)=0*dummyCoe(:,:,1,:);

myCoefficients = dummyCoe;%reverseCoefficients;%dummyCoe;%
gen.t_coefficients = myCoefficients;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD gamma genereation logic
gamma=1e-3;
lamda=1e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sq_error = zeros(noOfExperiments, noOfObservations, 4);
predict_error= zeros(noOfExperiments, noOfObservations, 4);
power= zeros (noOfExperiments, noOfObservations, 4);
generator=gen;
m_X1 = 1*generator.realization(); % time series generator
m_X=m_X1;%(m_X1'-mean(m_X1'))';
% r12=corr(m_X(1,:)',m_X(2,:)')
% r23=corr(m_X(2,:)',m_X(3,:)')

%m_X=rand(3,6);
%% joshin
% initialization
Gsigma=gen.Gsigma;
Gscale=1/sqrt(2*pi)/Gsigma;
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
%   alpha=randi([1,10],size(alpha))/1000;

% Kernal=zeros(noOfNodes,filtOrder,noOfObservations);
GausDenom=-1/2/Gsigma^2;
eta=1/gamma;
for t=filtOrder+1:noOfObservations
   for n1=1:noOfNodes
       m_X_P_pre=m_X(n1,t-filtOrder:t-1);
       m_X_P=flip(m_X_P_pre);
       m_X_P_big=repmat(m_X_P',1,t);
       Kernal(n1,:,1:t)=Gscale*exp((m_X_P_big-m_X(n1,1:t)).^2*GausDenom);
%        Kernal(n1,:,1:t)=1*((m_X_P_big-m_X(n1,1:t)).^2*1);

   end 
   K_vec=vec(Kernal(:,:,1:t));  
   for n1=1:noOfNodes
       grad_vec=(vec(alpha(n1,:,:,1:t)))'*K_vec*K_vec'-K_vec'*m_X(n1,t);
       grad=reshape(grad_vec,noOfNodes,filtOrder,t);
       g=eta*squeeze(alpha(n1,:,:,1:t))-grad;
       alpha(n1,:,:,1:t)= SoftThresold(g,t,noOfNodes,eta,lamda,n1);
   end
   (sum(alpha.^2,4)).^0.5;
end
 
Psuedo_Adj=(sum(alpha(:,:,:,60:100).^2,4)).^0.5

<<<<<<< HEAD
%%%%%%%%%%%%%%%
% Without dictionary
%%%%%%%%%%%%%%%
alpha=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
tic
for t=filtOrder+1:noOfObservations
    for n1=1:noOfNodes  
        m_X_P_pre=m_X(n1,t-filtOrder:t-1);
        m_X_P=flip(m_X_P_pre);
        m_X_P_big=repmat(m_X_P',1,t);
        Kernal(n1,:,1:t)=Gscale*exp(-.001*(m_X_P_big-m_X(n1,1:t)).^2);
    end 
    Kernal_Size_Full=size(Kernal);
    K_vec=vec(Kernal(:,:,1:t));
    for n1=1:noOfNodes
        predt(n1,t)=(vec(alpha(n1,:,:,1:t)))'*K_vec; 
        grad_vec=(vec(alpha(n1,:,:,1:t)))'*K_vec*K_vec'-K_vec'*m_X(n1,t);
        error_temp(n1)=(predt(n1,t)-m_X(n1,t))^2;
        grad=reshape(grad_vec,noOfNodes,filtOrder,t);
        g=eta*squeeze(alpha(n1,:,:,1:t))-grad;
        alpha(n1,:,:,1:t)= SoftThresold(g,t,noOfNodes,eta,lamda,n1,order);
    end
%         error(t)=sum(error_temp(1))/5;
end
Full_Ker_time=toc;
Psuedo_Adj=(sum(alpha.^2,4)).^0.5;
plot(predt(18,:),'LineWidth',2)
legend('original','prediction (dict)','prediction (full)')
alpha_full=alpha(:);
Num_samples_Dict=Kernal_Size_Dict(3);
Num_samples_Ful=Kernal_Size_Full(3);
title_str=strcat('NoSamplesDict=',num2str(Num_samples_Dict),', ',...
    'NoSamplesFull=',num2str(Num_samples_Ful))
title(title_str)

Kernal_Size_Dict
Kernal_Size_Full
Full_Ker_time
Dict_time

function [alpha_temp] =SoftThresold(g,t,noOfNodes,eta,lamda,n1,order)
=======
function [alpha_temp] =SoftThresold(g,t,noOfNodes,eta,lamda,n1)
if t==90
   stopflag=1; 
end
>>>>>>> e36bcfe75a67f56ed31d101d9d2fd84488dea84d
alpha_temp=zeros(size(g));
    for n2=1:noOfNodes
       for m=1:t
           if n2==n1
               alpha_temp(n2,:,m)=g(n2,:,m)/eta;
           else
                alpha_temp(n2,:,m)=g(n2,:,m)/eta*max(0,(1-lamda/norm(g(n2,:,m))));
           end
       end
    end
end