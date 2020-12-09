clc
clear all
close all
noOfNodes        = 5;%edited by joshin, it was 12
filtOrder        = 2;
noOfObservations = 2000; 
n_observation=2000;
T=400;
order=2;
number_of_nodes=5;
edge_probability=.05;
gamma=1e-2;
lamda=1/10000;
A=[];
T_Beta=.3*rand(order*number_of_nodes,number_of_nodes,T);
alpha=.1*rand(number_of_nodes,n_observation);
for i=1:order
temp1=5*rand(number_of_nodes,number_of_nodes).*(rand(number_of_nodes,number_of_nodes)<edge_probability);
A=[A;temp1];
end
y_initial=1*randn(number_of_nodes,1);
y=[y_initial];
Ksigma=1;
gsigma=Ksigma;
Gscale=(1/sqrt(2*pi)/Ksigma);
gscale=Gscale;
GausDenom=-1/2/Ksigma^2;
T_Beta=.5*rand(order*number_of_nodes,number_of_nodes,T);
% 
   A=zeros(10,5);
   A(1,5)=.6;
   A(3,2)=.9;
   
%    
   y_initial=1*randn(number_of_nodes,1);
y=[y_initial];
Ksigma=10;
gsigma=Ksigma;
Gscale=(1/sqrt(2*pi)/Ksigma);
gscale=Gscale;
GausDenom=-1/2/Ksigma^2;
for t=2:n_observation
    if t<=order
        y=[y zeros(number_of_nodes,1)];
        for j=1:t-1
            temp2=A(1+(j-1)*number_of_nodes:j*number_of_nodes,:)* y(:,t-j);
            y(:,t)=temp2-y(:,t);
        end
    end
    if t>order
         y=[y zeros(number_of_nodes,1)];
    for j=1:order
          temp2=A(1+(j-1)*number_of_nodes:j*number_of_nodes,:)*sum(alpha(:,1:t-j-1).*exp(GausDenom*(y(:,t-j)-y(:,1:t-j-1)).^2)*Gscale,2);
            y(:,t)=temp2+y(:,t);
    end
    end
end
M_x=zeros(number_of_nodes,n_observation);
M_x(:,1)=y_initial+.1*randn;
for ii=2:T+order
   temp=A(1:number_of_nodes,:)*y_initial;
   M_x(:,ii)=temp+.1*rand;
   y_initial= M_x(:,ii);
end
for ii=T+order+1:n_observation
    
    for k=1:number_of_nodes
        temp=0;
        for m=1:order
        for l=1:number_of_nodes
              temp=temp+ A(number_of_nodes*(m-1)+k,l)*(1*sin(M_x(l,ii-m)^6*3.14)+0*.3*sin(M_x(l,ii-m)*2*3.14)+0*.3*sin(M_x(l,ii-m)*3*3.14));
        end
        end
         M_x(k,ii)=temp+.01*randn;
          M_x(5,ii)=1*sin(ii*.03);
    end
    
    A(1,5)=1+1*sin(ii*.02)+.1*sin(ii*.010);
end
m_X=M_x;
Ksigma=10;
Gscale=(1/sqrt(2*pi)/Ksigma);
alpha1=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal1=zeros(noOfNodes,filtOrder,noOfObservations);
alpha2=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal2=zeros(noOfNodes,filtOrder,noOfObservations);
alpha3=zeros(noOfNodes,noOfNodes,filtOrder,noOfObservations);
Kernal3=zeros(noOfNodes,filtOrder,noOfObservations);
alpha=zeros(noOfNodes,noOfNodes,filtOrder,3*noOfObservations);


eta=1/gamma;
%% new code
m_Xnew(:,1:filtOrder)=m_X(:,1:filtOrder)
for t=filtOrder+1:noOfObservations
   
    [~  ,T]=size(m_Xnew);
     flag=0;
      for n11=1:noOfNodes
          
       temp= Gscale*exp(-1000*(m_X(n11,t)-m_Xnew(n11,1:T-1)).^2);
      
       if max(temp)<.034
           if flag==0
               T=t+1;
               flag=1;
               m_Xnew(:,T)=100000000000000000000000;
           end
         m_Xnew(n11,T)=  m_X(n11,t);
           
       end
      end
        alpha1=alpha(:,:,:,1:T);
    alpha2=alpha(:,:,:,T+1:2*T);
    alpha3=alpha(:,:,:,2*T+1:3*T);
   for n1=1:noOfNodes
     
       m_X_P_pre=m_X(n1,t-filtOrder:t-1);
       m_X_P=flip(m_X_P_pre);
       m_X_P_big=repmat(m_X_P',1,T);
              Kernal1(n1,:,1:T)=Gscale*exp(-100000*(m_X_P_big-m_Xnew(n1,1:T)).^2);
              Kernal2(n1,:,1:T)=Gscale*exp(-100*(m_X_P_big-m_Xnew(n1,1:T)).^2);
              Kernal3(n1,:,1:T)=Gscale*exp(-1000000000*(m_X_P_big-m_Xnew(n1,1:T)).^2);
   end 
   K_vec1=vec(Kernal1(:,:,1:T));  
   K_vec2=vec(Kernal2(:,:,1:T));  
   K_vec3=vec(Kernal3(:,:,1:T));  
   K_vec=[ K_vec1;K_vec2;K_vec3];
   for n1=1:noOfNodes
       alphavec=[vec(alpha1(n1,:,:,1:T));vec(alpha1(n1,:,:,1:T));vec(alpha1(n1,:,:,1:T))];
      predt(n1,t)=(alphavec)'*K_vec; 
       grad_vec=(alphavec)'*K_vec*K_vec'-K_vec'*m_X(n1,T);
        error_temp(n1)=(predt(n1,t)-m_X(n1,t))^2;
       grad=reshape(grad_vec,noOfNodes,filtOrder,3*T);
       g=eta*cat(3,squeeze(alpha1(n1,:,:,1:T)),squeeze(alpha2(n1,:,:,1:T)),squeeze(alpha3(n1,:,:,1:T)))-grad;
       alpha(n1,:,:,1:3*T)= SoftThresold(g,3*T,noOfNodes,eta,lamda,n1,order);
   end
     error(t)=sum(error_temp(1))/5;
   (sum(alpha.^2,4)).^0.5;
end
Psuedo_Adj=(sum(alpha.^2,4)).^0.5
 plot(m_X(1,1:2000))
hold on
plot(predt(1,:))
function [alpha_temp] =SoftThresold(g,t,noOfNodes,eta,lamda,n1,order)
if t==90
   stopflag=1; 
end
alpha_temp=zeros(size(g));
    for n2=1:noOfNodes
       
       for m=1:order
               alpha_temp(n2,m,:)=g(n2,m,:)/eta;
%                 
        
                alpha_temp(n2,m,:)=g(n2,m,:)/eta*max(0,(1-lamda/norm(reshape(g(n2,m,:),1,t))));
        
       end
      
    end
end

