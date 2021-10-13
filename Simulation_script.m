    %function lineage_deterministic(q)
rng('shuffle');

Tg0=120;
mu0=log(2)/Tg0;

totaltime=5*Tg0;

%% number of proteins
ntot=2e3;
np=ntot-2;


taum0=10; %in the unit of the number per micrometer qubic
CVtaum=0;
sigmay=sqrt(log(1+CVtaum^2));
meany=log(taum0)-sigmay^2/2;
taum=exp(meany+sigmay*randn(1,ntot));
taum(1:2)=taum0;

decayp=zeros(1,ntot);
% decayp(10)=1/10;
% decayp(11)=1/10;
% decayp(12)=1/10;

ndecay=0;
%decayp(ntot-ndecay+1:end)=1/10;

vr=12*60; % number of amnino acids per min
vn=12*60; 

% The length is measured in the number of codons
len=500*ones(1,np+2);
len(1)=1e3;
len(2)=1e4;
Lr=len(2);
Ln=len(1);


% critical number of RNAP and ribosomes
nc=1e4;

% typical cell mass at birth
Mb=1e9;

phi_r=mu0*len(2)/vr;
phi_n=len(1)/len(2)*phi_r*0.1;

%% MM constant for binding of RNAP and ribosome
% the density is in the unit of amino acids per micrometer qubic
Kn0=10*1e-6*6e23*10^-15; %in the unit of the number per micrometer qubic

Kn_mean=5e3;

CVKn=0;
sigmay=sqrt(log(1+CVKn^2));
meany=log(Kn0)-sigmay^2/2;
Kn=exp(meany+sigmay*randn(1,ntot));

% Kn(1)=Kn_mean;
% Kn(2)=Kn_mean;

% Kn=Kn0*ones(1,ntot);
% Kn(10)=20*Kn0;
% Kn(11)=0.2*Kn0;
% Kn(1:2)=Kn0;


Kr=10*1e-6*6e23*10^-15;

rho=1e10;
%ratio between total protein mass and cell volume (num of aa per micro qubic)

alpha=20; %ratio between nuclear volume and cytoplasmic volume, assumed to be constant


% An approximate expression of gene copy number assuming cn*Fn << Kn
% Once g_tot is chosen, it is fixed
gr=5; %all the other genes have gene copy number 1
g_tot=ones(1,ntot);
g_tot(2)=gr;

y=phi_n/phi_r*(g_tot(2)*taum(2)*len(2)/Kn(2))/(g_tot(1)*taum(1)*len(1)/Kn(1));

A=g_tot(2)*taum(2).*len(2)./Kn(2);
B=g_tot(1)*taum(1).*len(1)./Kn(1);
C=sum(g_tot(3:end).*taum(3:end).*len(3:end)./Kn(3:end));
    

x=(1/(phi_r)-1-y*B/A)*A/C;

kcatn_r=(nc-sum(g_tot))*vn/(g_tot(2)*len(2)+y*g_tot(1)*len(1)+...
    x*sum(g_tot(3:end).*len(3:end)));

kcatn_n=y*kcatn_r;

kcatn=zeros(1,ntot);

kcatn(1)=kcatn_n;
kcatn(2)=kcatn_r;
kcatn(3:end)=kcatn_r*x;


kcatr=10;


Lambda_n=kcatn.*len/vn;
Lambda_r=kcatr*len/vr;


%% we first compute the fraction of free RNAP and ribosome
n=Mb*phi_n/len(1)
r=Mb*phi_r/len(2)
n0=n;
r0=r;

V=Mb/rho;
Vn=V/alpha; % nuclear volume

cn=n/Vn;

Gn=g_tot.*(1+Lambda_n);
temp=0;
test_0=1e8;

while(1)
    temp1=n*(1-temp);
    temp3=sum(Gn.*(cn*temp./(cn*temp+Kn)));
    test=abs(temp1-temp3);
    if test<test_0
        test_0=test;
        temp=temp+1e-5;
    else
        break;
    end
end
Fn=temp;

m_tot=kcatn.*g_tot.*(cn*Fn./(cn*Fn+Kn)).*taum;

Gr=sum(m_tot.*(1+Lambda_r));
cr=r/(V-Vn);
temp=0;
test_0=1e8;
while(1)
    temp2=r*(1-temp);
    temp4=Gr*(cr*temp./(cr*temp+Kr));
    test=abs(temp2-temp4);
    if test<test_0
        test_0=test;
        temp=temp+1e-5;
    else
        break;
    end
end
Fr=temp;

phi_tot=m_tot.*len/sum(m_tot.*len);

P_tot=Mb*phi_tot./len;

P_tot=P_tot*n/P_tot(1);% we make sure the number of RNAP is continuous in the beginning

M=sum(P_tot.*len);
Mb=M;

%Kn(1:2)=sum(phi_tot.*Kn);


tra=zeros(totaltime,np+3);
tram=zeros(totaltime,np+3);
tram1=zeros(totaltime,np+3);
Mt=zeros(totaltime,3);
Ft=zeros(totaltime,5);


next=1;
deltaM=0;
%%
record=1;
tnext=0;
t=0;

count=1;
while(1)
    
    if t>=tnext
        tra(record,:)=[t P_tot];

        tram(record,:)=[t m_tot];
        
        temp=kcatn.*g_tot.*cn*Fn./(cn*Fn+Kn);
        tram1(record,:)=[t temp];
        
        Mt(record,:)=[t M deltaM];
        Ft(record,:)=[t Fn Fr cn*Fn./(cn*Fn+Kn0) cr*Fr/(cr*Fr+Kr)];
        record=record+1;
        tnext=t+1;
    end
    
    % we first compute the fraction of free RNAP and ribosome
    n=P_tot(1);
    r=P_tot(2);
    Gr=sum(m_tot.*(1+Lambda_r));
    
    
    V=M/rho;
    Vn=V/alpha;
    cn=n/Vn;
    cr=r/(V-Vn);
    
    test_0=1e8;
    
    temp=0;
    dtemp=1e-5;
    while(1)
        temp1=n*(1-temp);
        temp3=sum(Gn.*(cn*temp./(cn*temp+Kn)));
        test=abs(temp1-temp3);
        if test<test_0
            test_0=test;
            temp=temp+dtemp;
        else
            break;
        end
    end
    Fn=temp;
    
    temp=0;
    test_0=1e8;
    dtemp=1e-5;
    while(1)
        temp2=r*(1-temp);
        temp4=Gr*(cr*temp./(cr*temp+Kr));
        test=abs(temp2-temp4);
        if test<test_0
            test_0=test;
            temp=temp+dtemp;
        else
            break;
        end
    end
    
    Fr=temp;
    
    
    rate_transcription=kcatn.*g_tot.*(cn*Fn./(cn*Fn+Kn));
    
    rate_mRNAdecay=m_tot./taum;
    
    rate_translation=kcatr*m_tot.*(cr*Fr/(cr*Fr+Kr));
    
    rate_proteindecay=P_tot.*decayp;
    
    
    delta_t=5/max(rate_translation);  
    
    t=t+delta_t;
    
    m_tot=m_tot+rate_transcription*delta_t-rate_mRNAdecay*delta_t;
    
    P_tot=P_tot+rate_translation*delta_t-rate_proteindecay*delta_t;
    
    Mold=M;
    M=sum(P_tot.*len);
    deltaM=deltaM+M-Mold;
    
    
    if  deltaM>=8*Mb
        
        break;
     
    end
end

tra=tra(1:record-1,:);
tram=tram(1:record-1,:);
tram1=tram1(1:record-1,:);
Mt=Mt(1:record-1,:);
Ft=Ft(1:record-1,:);
mtot=sum(tram(:,1:end),2);

phi=P_tot.*len/M;
