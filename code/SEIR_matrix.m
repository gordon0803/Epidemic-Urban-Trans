  sigma=0.1; %1/sigma is the latent period of the disease
  gamma=0.1; %1/gamma is the recovery rate of the disease
  betainn=0.1; %within zone contact rate
  betatra=[0.3,0.2,0.1]; %contactrate in travel mode m, m=1,2,3, which represent large, medium and small capacity
  d=[0.6,0.8,1]; %level of controls for mode m.
%   c12=[0.4,0.3,0.3];  %ratio of people move from i to j who choose travel mode m, i,j=1,2,3
%   c13=[0.4,0.4,0.2];
%   c23=[0.4,0.3,0.3];
%   c21=[0.4,0.3,0.3];
%   c31=[0.4,0.4,0.2];
%   c32=[0.4,0.3,0.3];
  c1=[0,0.4,0.4;0.4,0,0.4;0.4,0.4,0]; %mode1 c1(i,j)
  c2=[0,0.3,0.4;0.3,0,0.3;0.4,0.3,0]; %mode2
  c3=[0,0.3,0.2;0.3,0,0.3;0.2,0.3,0]; %mode3
  C=cell(1,3);
  C{1}=c1;
  C{2}=c2;
  C{3}=c3;
  Np=[1,1,1]; % population of patch i
  N=[0.5,0.3,0.2;0.3,0.4,0.3;0.2,0.3,0.5];  %the amount of people currently at patch j who are the residents of patch i
  g=[0.5,0.6,0.5]; %total departure rate of patch i
  m=[0,0.6,0.4;0.5,0,0.5;0.4,0.6,0]; %the rate of movement from patch i to patch j
  r=[0,0.6,0.3;0.6,0,0.3;0.3,0.6,0]; %return rate of people who have traveled from patch i to patch j
  M=zeros(36,117); %SEIR matrix
  R=[2,3;1,3;1,2];
  
  Path=[0,0,0,0,0,0,0,0,0;
        0,1,1,0,0,0,0,0,0; 
        0,1,1,0,0,1,0,0,0;
        0,0,0,1,0,0,1,0,0;
        0,0,0,0,0,0,0,0,0;
        0,0,1,0,0,1,0,0,0;
        0,0,0,1,0,0,1,1,0;
        0,0,0,0,0,0,1,1,0;
        0,0,0,0,0,0,0,0,0];

%S
n_zone=3;
index=1;
I=cell(1,9);
I{1}=[1,1];I{2}=[1,2];I{3}=[1,3];I{4}=[2,1];I{5}=[2,2];I{6}=[2,3];I{7}=[3,1];I{8}=[3,2];I{9}=[3,3];
P=cell(1,2);
P{1}=zeros(9,3);
P{2}=zeros(9,3);
L=zeros(9,3);
Q=zeros(9,3);
K=zeros(9,3);
for i=1:n_zone
    for j=1:n_zone
        M(n_zone*(i-1)+j,36+27*(i-1)+9*(j-1)+3*(j-1)+j)=-betainn/(Np(j)); %Sij*Ijj
        for k=1:2
        M(n_zone*(i-1)+j,36+27*(i-1)+9*(j-1)+3*(R(i,k)-1)+j)=-betainn*[C{1}(R(i,k),j),C{2}(R(i,k),j),C{3}(R(i,k),j)]*d'; %Sij*Ikj
        end
        if i==j %Sii_dot
            M(n_zone*(i-1)+j,n_zone*(i-1)+j)=-g(i); %Sii
            for l=1:n_zone-1
            M(n_zone*(i-1)+j,n_zone*(i-1)+R(i,l))=r(i,R(i,l)); %Sij
            end
        else %Sij_dot
            M(n_zone*(i-1)+j,n_zone*(i-1)+j)=-r(i,j); %Sij
            M(n_zone*(i-1)+j,n_zone*(i-1)+i)=g(i)*m(i,j); %Sii
            %infected via transportation
            for x=1:9 
                if Path(n_zone*(i-1)+j,x)==1
                    L(n_zone*(i-1)+j,index)=x;
                    P{1}(n_zone*(i-1)+j,index)=I{x}(1);
                    P{2}(n_zone*(i-1)+j,index)=I{x}(2);
                    Q(n_zone*(i-1)+j,index)=N(I{x}(1),I{x}(2)); %Nkj
                    index=index+1;
                else
                end
            end
            for l=1:3 %different modes of Nkj
                for y=1:n_zone
                    K(n_zone*(i-1)+j,l)=K(n_zone*(i-1)+j,l)+C{l}(P{1}(n_zone*(i-1)+j,y),P{2}(n_zone*(i-1)+j,y))*Q(n_zone*(i-1)+j,y);
                end
            end
            for z1=1:n_zone %different Ikj
                for z2=1:3 %summation over mode
                    if Q(n_zone*(i-1),z1)==0
                    else
                    M(n_zone*(i-1)+j,36+27*(i-1)+9*(j-1)+L(n_zone*(i-1)+j,z1))=...
                    M(n_zone*(i-1)+j,36+27*(i-1)+9*(j-1)+L(n_zone*(i-1)+j,z1))+...
                    C{z2}(i,j)*betatra(z2)*C{z2}(P{1}(n_zone*(i-1)+j,z1),P{2}(n_zone*(i-1)+j,z1))*d(z2)/K(n_zone*(i-1)+j,z2);
                    end
                end
            end
        end
    end
end



   
%E
for i=1:3
    for j=1:3
        M(9+3*(i-1)+j,36+27*(i-1)+9*(j-1)+3*(j-1)+j)=-betainn/(Np(j)); %Sij*Ijj
        for k=1:2
        M(9+3*(i-1)+j,36+27*(i-1)+9*(j-1)+3*(R(i,k)-1)+j)=-betainn*[C{1}(R(i,k),j),C{2}(R(i,k),j),C{3}(R(i,k),j)]*d'; %Sij*Ikj
        end
        if i==j %Eii_dot
           M(9+3*(i-1)+j,9+3*(i-1)+j)=-g(i)-sigma; %Eii
           for l=1:2
           M(9+3*(i-1)+j,9+3*(i-1)+R(i,l))=r(i,R(i,l));%Eij
           end
        else %Eij_dot
           M(9+3*(i-1)+j,9+3*(i-1)+j)=-r(i,j)-sigma; %Eij
           M(9+3*(i-1)+j,3*(i-1)+j)=g(i)*m(i,j); %Eii
           for l=1:3 %infected via transportaion
               M(9+n_zone*(i-1)+j,36+27*(i-1)+9*(j-1)+L(n_zone*(i-1)+j,l))= M(n_zone*(i-1)+j,36+27*(i-1)+9*(j-1)+L(n_zone*(i-1)+j,l));
           end
        end
    end
end

%I
for i=1:3
    for j=1:3
        M(18+3*(i-1)+j,9+3*(i-1)+j)=sigma; %Eii/Eij
        if i==j %Iii_dot
           M(18+3*(i-1)+j,18+3*(i-1)+j)=-g(i)-gamma; %Iii
           for k=1:2
           M(18+3*(i-1)+j,18+3*(i-1)+R(i,k))=r(i,R(i,k))*[C{1}(i,R(i,k)),C{2}(i,R(i,k)),C{3}(i,R(i,k))]*d'; %Iij
           end
        else %Iij_dot
           M(18+3*(i-1)+j,18+3*(i-1)+j)=-r(i,j)-[C{1}(i,j),C{2}(i,j),C{3}(i,j)]*(gamma*d'+(1-d')); %Iij
           M(18+3*(i-1)+j,18+3*(i-1)+i)=g(i)*m(i,j)*[C{1}(i,j),C{2}(i,j),C{3}(i,j)]*d';  %Iii
        end
    end
end

%R
for i=1:3
    for j=1:3
        if i==j %Rii_dot
            M(27+3*(i-1)+j,27+3*(i-1)+j)=-g(i); %Rii
            M(27+3*(i-1)+j,18+3*(i-1)+j)=gamma; %Iii
            for k=1:2
            M(27+3*(i-1)+j,27+3*(i-1)+R(i,k))=r(i,R(i,k)); %Rij
            M(27+3*(i-1)+j,18+3*(i-1)+R(i,k))=r(i,R(i,k))*[C{1}(i,R(i,k)),C{2}(i,R(i,k)),C{3}(i,R(i,k))]*(1-d'); %Iij
            end
        else
            M(27+3*(i-1)+j,27+3*(i-1)+j)=-r(i,j); %Rij
            M(27+3*(i-1)+j,27+3*(i-1)+i)=g(i)*m(i,j); %Rii
            M(27+3*(i-1)+j,18+3*(i-1)+j)=[C{1}(i,j),C{2}(i,j),C{3}(i,j)]*(gamma*d'+(1-d')); %Iij
        end
    end
end


