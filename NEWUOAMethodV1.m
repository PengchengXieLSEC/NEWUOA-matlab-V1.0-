function [Fopt,xopt,NF] = NEWUOAMethodV1(Fobj,M,N,xbeg,rhobeg,rhoend,Max)
%NEWUOAMETHOD NEWUOA算法的函数
%   TRSAPP By us
%   F表示目标函数,m是插值点个数,n是问题维数,rho_beg是初始的rho,rho_end是终止的rho,Max是迭代的最大次数
global Xn Fn n m F_times rho_beg rho_end x0 opt xb c g Gamma gamma H F rho delta Krho D3 QF3 CRVMIN d NORMD Qnew NFMAX RATIO MOVE w Hw beta Fnew DIST XXopt NXX flag
F=Fobj;
m=M;
n=N;
xb=xbeg;
rho_beg=rhobeg;
rho_end=rhoend;
NFMAX=Max;
flag=1;
while flag<16
    switch flag
        case 1
            NEWUOAStep1();
        case 2
            NEWUOAStep2();
        case 3
            NEWUOAStep3();
        case 4
            NEWUOAStep4();
        case 5
            NEWUOAStep5();
        case 6
            NEWUOAStep6();
        case 7
            NEWUOAStep7();
        case 8
            NEWUOAStep8();
        case 9
            NEWUOAStep9();
        case 10
            NEWUOAStep10();
        case 11
            NEWUOAStep11();
        case 12
            NEWUOAStep12();
        case 13
            NEWUOAStep13();
        case 14
            NEWUOAStep14();
        case 15
            NEWUOAStep15();
    end
end
Fopt=Fn(opt);
xopt=Xn(:,opt);
NF=F_times;
end


function [] = NEWUOAStep1()
%NEWUOAStep1 对应NEWUOA算法的第一步
%   在第一步中，会根据初始点生成初始的插值节点，以及计算它们处的函数值，同时计算此时的最小值的Fopt以及最小值点xopt
%   Xn是插值点矩阵，W是插值问题的矩阵，H是W的逆（不用求逆算法得）,Fn是插值处目标函数的函数值
%   最初始的插值函数Q(x)=c+g'*(x-x0)+0.5*(x-x0)'*G*(x-x0);
%   W=[A,X';X,0]; H=[s.*Z*Z';Theta';Theta,Upsilon];
global Xn Fn n m F_times rho_beg x0 opt xb c g Gamma gamma H F rho delta Krho D3 QF3 XX0 flag
global Steps %测试用
Steps=1;
gamma=zeros(m,1);
Xn=zeros(n,m);
Fn=zeros(m,1);
x0=xb;
Xn(:,1)=x0;
Fn(1)=F(Xn(:,1));
I=eye(n);
Gamma=zeros(n,n);
c=Fn(1);
g=zeros(n,1);
%%
%设置参数
rho=rho_beg;
delta=rho;
Krho=0;
D3=zeros(3,1);
QF3=zeros(3,1);
%%
%生成插值点
if m>=2*n+1
    Delta=zeros(n,1);
    for i=1:n
        Xn(:,i+1)=x0+rho*I(:,i);
        Fn(i+1)=F(Xn(:,i+1));
        Xn(:,i+n+1)=x0-rho*I(:,i);
        Fn(i+n+1)=F(Xn(:,i+n+1));
        if Fn(i+n+1)>=Fn(i+1)
            Delta(i)=1;
        else
            Delta(i)=-1;
        end
        g(i)=(Fn(i+1)-Fn(i+n+1))/(2*rho);
        Gamma(i,i)=(Fn(i+1)+Fn(i+n+1)-2*Fn(1))/(rho*rho);
    end
    
    for i=2*n+2:1:m
        j=floor((i-n-2)/n);
        p=i-n-1-j*n;
        if p+j>n
            q=p+j-n;
        else
            q=p+j;
        end
        Xn(:,i)=x0+rho*Delta(p)*I(:,p)+rho*Delta(q)*I(:,q);
        Fn(i)=F(Xn(:,i));
        Gamma(p,q)=( Fn(1)-Fn(1+p+(1-Delta(p))*n/2)-Fn(1+q+(1-Delta(q))*n/2)+Fn(i) )/(Delta(p)*Delta(q)*rho*rho);
        Gamma(q,p)=Gamma(p,q);
    end
else
    for i=1:n
        Xn(:,i+1)=x0+rho*I(:,i);
        Fn(i+1)=F(Xn(:,i+1));
    end
    for i=n+1:m-1
        Xn(:,i+1)=x0-rho*I(:,i-n);
        Fn(i+1)=F(Xn(:,i+1));
        g(i-n)=(Fn(i-n+1)-Fn(i+1))/(2*rho);
        Gamma(i-n,i-n)=(Fn(i-n+1)+Fn(i+1)-2*Fn(1))/(rho*rho);
    end
    for i=m:2*n
        g(i-n)=(Fn(i-n+1)-Fn(1))/(rho);
    end
end
F_times=m;
%%
%计算W
% XX0=Xn-x0;
% A=((XX0'*XX0).^2)/2;
% X=[ones(1,m);XX0];
% W=zeros(m+n+1,m+n+1);
% W(1:m,1:m)=A;
% W(m+1:m+n+1,1:m)=X;
% W(1:m,m+1:m+n+1)=X';
%%
%计算W的逆H
Theta=zeros(n+1,m);
Theta(1,1)=1;
Upsilon=zeros(n+1,n+1);
if m>=2*n+1
    for i=2:n+1
        Theta(i,i)=1/(2*rho);
        Theta(i,i+n)=-1/(2*rho);
    end
else
    for i=2:m-n
        Theta(i,i)=1/(2*rho);
        Theta(i,i+n)=-1/(2*rho);
    end
    for i=m-n+1:n+1
        Theta(i,1)=-1/rho;
        Theta(i,i)=1/rho;
        Upsilon(i,i)=-0.5*rho*rho;
    end
end
Z=zeros(m,m-n-1);
s=ones(1,m-n-1);
if m<=2*n+1
    for k=1:m-n-1
        Z(1,k)=-sqrt(2)/(rho*rho);
        Z(k+1,k)=-Z(1,k)/2;
        Z(k+n+1,k)=Z(k+1,k);
    end
else
    for k=1:n
        Z(1,k)=-sqrt(2)/(rho*rho);
        Z(k+1,k)=-Z(1,k)/2;
        Z(k+n+1,k)=Z(k+1,k);
    end
    for k=n+1:m-n-1
        i=k+n+1;
        j=floor((i-n-2)/n);
        p=i-n-1-j*n;
        if p+j>n
            q=p+j-n;
        else
            q=p+j;
        end
        if Delta(p)==1
            pt=p+1;
        else
            pt=p+1+n;
        end
        if Delta(q)==1
            qt=q+1;
        else
            qt=q+1+n;
        end
        Z(1,k)=1/(rho*rho);
        Z(pt,k)=-Z(1,k);
        Z(qt,k)=Z(pt,k);
        Z(k+n+1,k)=Z(1,k);
    end
end
H=zeros(m+n+1,m+n+1);
H(1:m,1:m)=Z*Z';
H(1:m,m+1:m+n+1)=Theta';
H(m+1:m+n+1,1:m)=Theta;
H(m+1:m+n+1,m+1:m+n+1)=Upsilon;
%%
%计算二次模型系数c,g,lambda
% R=zeros(m+n+1,1);
% R(1:m)=Fn;
% Lcg=H*R;
% lambda=Lcg(1:m);
% c=Lcg(m+1);
% g=Lcg(m+2:m+n+1);
%%
[~,opt]=min(Fn);
XX0=Xn-x0;
% %% 用于测试数组
% global Test
% Test(F_times-m+1,:)=[Xn(:,opt)',Fn(opt)];
%%
flag=2;
%NEWUOAStep2();
end

function [] = NEWUOAStep2()
%NEWUOAStep2 对应NEWUOA算法的第二步
%   插值多项式为Q(x)=c+g'*(x-x0)+0.5*(x-x0)'*G*(x-x0)
%   G=Gamma+sum(lambda(j)(Xn(:,j)-x0)*(Xn(:,j)-x0)');
%   目标函数为Q(xopt+d)
%   约束条件为 ||d||<=delta
%   方法为截断共轭梯度法
%   XX=Xn-x0;
global n g x0 Xn Fn opt delta CRVMIN d NORMD Qnew flag
global Steps %测试用
Steps=[Steps;2];
xopt=Xn(:,opt);
Fopt=Fn(opt);
d=zeros(n,1);
NORMD=0;
ITERC=0;
ITERMAX=n;
ITERSW=ITERMAX;
DELSQ=delta*delta;
%% 截断共轭梯度法
QRED=0;
SS=0;
s=xopt-x0;
HS=qHD(s);
d=zeros(n,1);
DD=0;
HD=zeros(n,1);
G=g+HS;%Q'(xopt)=g+H*(xopt-x0);
s=-G;
SS=s'*s;
CRVMIN=0;%默认设置为0
DS=0;
DD=0;
GG=SS;
GGBEG=GG;
K1=0;%(5.13)(2)
K2=0;%(5.13)(1)
K3=0;
K4=0;
K5=0;
if (SS < 10e-15)
    Qnew=Fopt;
    DD=0;
else
    while(ITERC<ITERMAX)
        ITERC=ITERC+1;
        TEMP=DELSQ-DD;
        BSTEP=TEMP/(DS+sqrt(DS*DS+SS*TEMP));
        HS=qHD(s);
        SHS=s'*HS;
        ALPHA=BSTEP;
        if(SHS> 0)
            TEMP=SHS/SS;%公式(5.14)
            if (ITERC == 1)
                CRVMIN=TEMP;
            end
            CRVMIN=min([CRVMIN,TEMP]);
            ALPHA=min([ALPHA,GG/SHS]);%公式(5.10)
        end
        QADD=ALPHA*(GG-0.5*ALPHA*SHS);% 公式(5.8)
        QRED=QRED+QADD;
        GGSAV=GG;
        GG=0;
        d=d+ALPHA*s;
        HD=HD+ALPHA*HS;
        GG=(G+HD)'*(G+HD);
        if (ALPHA < BSTEP)
            if (QADD <=0.01D0*QRED)%公式(5.13)下面那个
                K1=1;
                break;
            end
            if (GG <= 1.0D-4*GGBEG)%公式(5.13)上面那个
                K2=1;
                break;
            end
            TEMP=GG/GGSAV; %公式(5.5)中的BETA_j
            DD=0;
            DS=0;
            SS=0;
            s=TEMP*s-G-HD;
            SS=s'*s;
            DS=d'*s;
            DD=d'*d;
            %             NORMD=sqrt(d);
            if (DS <= 0)
                K3=1;
                break;
            end
            if(DD>=DELSQ)
                K4=1;
                break;
            end
        else
            %             NORMD=delta;
            DD=DELSQ;
            K5=1;
            break;
        end
    end
    if K5==1 ||(K1==0&&K2==0&&K3==0&&K4==1&&(ITERC<ITERMAX))
        CRVMIN=0;% DD=DELSQ,且前面几个不等式都不成立
        ITERSW=ITERC;
    end
    %     xnew=xopt+d;
    Qnew=Fopt+d'*G+0.5*d'*HD;%Qnew=Fopt+Q'(xopt)*d+0.5d*H*d
    %% 额外部分
    if(  (K1+K2+K3<1)&&(ITERC<ITERMAX) )
        while(GG>1.0e-4*GGBEG)%进入额外部分的条件1
            DG=d'*G;
            DHD=d'*HD;
            DGK=DG+DHD;
            ANGTEST=DGK/sqrt(GG*DELSQ);
            if (ANGTEST > -0.99D0) %进入额外部分的条件2
                %49个辅助量
                CTH=zeros(49,1);
                STH=zeros(49,1);
                ANGLE=zeros(49,1);
                dPi=2*pi/50;
                for i=1:49
                    ANGLE(i)=dPi*i;
                    CTH(i)=cos(ANGLE(i));
                    STH(i)=sin(ANGLE(i));
                end
                %
                ITERC=ITERC+1;
                TEMP=sqrt(DELSQ*GG-DGK*DGK);
                TEMPA=DELSQ/TEMP;
                TEMPB=DGK/TEMP;
                s=TEMPA*(G+HD)-TEMPB*d;
                HS=qHD(s);
                SG=0;
                SHS=0;
                SHD=0;
                SG=s'*G;
                SHS=s'*HS;
                SHD=HS'*d;
                QBEG=Fopt+DG+0.5*DHD;%Q(xopt+d);Q'(xopt+d)-Q'(xopt)=g+H(xopt+d-x0)-g-H(xopt-x0)=HD;
                QMIN=QBEG;
                QSAV=QBEG;
                ISAVE=0;
                for i=1:49
                    QNEW=Fopt+CTH(i)*DG+STH(i)*SG+0.5*(CTH(i)^2*DHD+STH(i)^2*SHS)+CTH(i)*STH(i)*SHD;
                    if (QNEW < QMIN)
                        QMIN=QNEW;
                        ISAVE=i;
                        TEMPA=QSAV;
                    else
                        if (i == ISAVE+1)
                            TEMPB=QNEW;
                        end
                    end
                    QSAV=QNEW;
                end
                if (ISAVE == 0)
                    TEMPA=QNEW;
                end
                if (ISAVE == 49)
                    TEMPB=QBEG;
                end
                angle=0;
                if (TEMPA ~= TEMPB)
                    TEMPA=TEMPA-QMIN;
                    TEMPB=TEMPB-QMIN;
                    angle=0.5*(TEMPA-TEMPB)/(TEMPA+TEMPB);
                end
                angle=dPi*(ISAVE+angle);
                cosA=cos(angle);
                sinA=sin(angle);
                Qnew=Fopt+cosA*DG+sinA*SG+0.5*(cosA^2*DHD+sinA^2*SHS)+cosA*sinA*SHD;
                d=cosA*d+sinA*s;
                HD= cosA*HD+sinA*HS;
                GG=(G+HD)'*(G+HD);
                RATIO=(QBEG-Qnew)/(Fopt-Qnew);
                if( (RATIO<=0.01) || (ITERC >= ITERMAX) )
                    break;
                end
            else
                break;
            end
        end
        DD=DELSQ;
    end
end
%%
NORMD=sqrt(DD);
%%
%NEWUOAStep3();
flag=3;
end
function [hs]=qHD(s)
%计算矩阵G乘向量s的快速算法，这里s是n*1的
global Gamma gamma XX0 m
eta=gamma.*(XX0'*s);
hs=Gamma*s;
for i=1:m
    hs=hs+eta(i)*XX0(:,i);%%公式(5.3)
end
end


function [] = NEWUOAStep3()
%NEWUOAStep3 对应NEWUOA算法的第三步

global NORMD rho flag
global Steps %测试用
Steps=[Steps;3];
if NORMD>0.5*rho
    %NEWUOAStep4();%To NEWUOAStep4
    flag=4;
else
    %NEWUOAStep14();%To NEWUOAStep14
    flag=14;
end

end

function [] = NEWUOAStep4()
%NEWUOAStep4 对应NEWUOA算法的第四步
% 计算比例系数RATIO，以此作为迭代delta的依据。
% 把MOVE设置为需要被替换的点
global m n RATIO Fn Xn F F_times delta opt d NORMD rho MOVE Qnew H x0 w Hw beta Fnew Krho XX0 NFMAX flag
global DandR Steps Test %测试用
Steps=[Steps;4];
Fopt=Fn(opt);
xopt=Xn(:,opt);
xnew=xopt+d;%新生成的点
F_times=F_times+1;
if F_times>NFMAX
    F_times=F_times-1;
    flag=16;
else
    Fnew=F(xnew);
    Krho=Krho+1;
    %%
    %计算RATIO
    % QNEW=Q(xnew);
    dQ=Fopt-Qnew;
    dF=Fopt-Fnew;
    % DDQ=Q(Xn(:,opt))-Q(xnew);
    % R1=dF/dQ;
    % R2=dF/DDQ;
    if dQ<=0 %%二次模型没有下降
        if dF>0 %但是实际值是下降的
            RATIO=0;%去修改模型
        else%实际值也没有下降
            RATIO=-1;
        end
    else
        RATIO=dF/dQ;
    end
    %% 测试
    [E,en] = InterpolationError();
    DandR=[DandR;delta,RATIO,rho,E];
    %Test(F_times-m,:)=[Xn(:,opt)',Fn(opt)];
    Test=[Test;Xn(:,opt)',Fn(opt),F_times];
    %%
    %更新delta
    if RATIO<=0.1
        deltaint=0.5*NORMD;
    else
        if RATIO<=0.7
            deltaint=max([NORMD,0.5*delta]);
        else
            deltaint=max([2*NORMD,0.5*delta]);
        end
    end
    if deltaint>1.5*rho
        delta=deltaint;
    else
        delta=rho;
    end
    %%
    %挑选出被替代的点
    T=1:m;%待挑选的点
    if dF>0
        Case=1;
        xstar=xopt+d;%之后的最优值点
    else
        xstar=xopt;%之后的最优值点
        T(opt)=[];
        Case=-1;
    end
    w=zeros(m+n+1,1);
    dxx0=xnew-x0;
    for i=1:m
        w(i)=0.5*(((XX0(:,i))'*(dxx0))^2);%XX0(:,i)=Xn(:,i)-x0
    end
    w(m+1)=1;
    w(m+2:m+n+1)=dxx0;
    Hw=H*w;
    beta=0.5*((dxx0'*dxx0)^2)-w'*Hw;
    M1=max([0.1*delta,rho]);
    if Case>0 %新的点让函数值下降
        alpha=H(1,1);
        tau=Hw(1);
        Sigma=alpha*beta+tau*tau;
        Weight=max([1,(norm(Xn(:,1)-xstar)/M1)^6]);
        WDOLD=Weight*abs(Sigma);
        TSAVE=1;
        for i=2:m
            alpha=H(i,i);
            tau=Hw(i);
            Sigma=alpha*beta+tau*tau;
            Weight=max([1,(norm(Xn(:,i)-xstar)/M1)^6]);
            WDNEW=Weight*abs(Sigma);
            if WDNEW>WDOLD
                WDOLD=WDNEW;
                TSAVE=i;
            end
        end
        MOVE=TSAVE;
    else %新的点没有让函数值下降
        WDOLD=0;
        TSAVE=0;
        for i=1:m
            if i~=opt
                alpha=H(i,i);
                tau=Hw(i);
                Sigma=alpha*beta+tau*tau;
                Weight=max([1,(norm(Xn(:,i)-xstar)/M1)^6]);
                WDNEW=Weight*abs(Sigma);
                if WDNEW>WDOLD
                    WDOLD=WDNEW;
                    TSAVE=i;
                end
            end
        end
        if WDOLD<=1
            MOVE=0;
        else
            MOVE=TSAVE;
        end
    end
    
    %%
    %NEWUOAStep5();%To Setp5
    flag=5;
end
end

function [] = NEWUOAStep5()
%NEWUOAStep5 对应NEWUOA算法的第五步
% 根据前面的计算更新插值节点
global m n d MOVE Krho D3 QF3 NORMD Xn Fn opt x0 Qnew Fnew c g Gamma gamma H w Hw beta RATIO XX0 flag
global Steps %测试用
Steps=[Steps;5];
%%
%计算步骤14的参变量
K=mod(Krho,3)+1;
D3(K)=NORMD;
QF3(K)=abs(Qnew-Fnew);
%%
t=MOVE;
if t>0
    %%
    %判断是否该改变x0
    XX=XX0;%XX=Xn-x0;
    % s=Xn(:,opt)-x0;
    s=XX(:,opt);
    NS2=s'*s;
    xnew=Xn(:,opt)+d;
    if (NORMD*NORMD<=(0.001*NS2 ))
        %此时替换x0为xopt
        %由此，需要修改的还有 H，c,g,Gamma,gamma,w,Hw ,beta
        xav=0.5*(x0+Xn(:,opt));
        %s=Xn(:,opt)-x0;
        %     [HE1,W1,Hreal1] = HError();
        %     [E1,en1] = InterpolationError();
        %%
        %修改H
        %     %%
        %     Y=XX-0.5*s;
        %     Z=zeros(n,m);
        %     for i=1:m
        %        Z(:,i)=(s'*Y(:,i))* Y(:,i);
        %     end
        %     Z=Z+0.25*(NS2)*s;
        %     G1=eye(m+n+1);
        %     G2=eye(m+n+1);
        %     G1(m+1,m+2:m+n+1)=0.5*s;
        %     G2(m+2:m+n+1,1:m)=Z;
        %     H=G1*G2*G1*H*G1'*G2'*G1';
        %%
        %修改其他参数
        x0=Xn(:,opt);
        v=zeros(n,1);
        for i=1:m
            v=v+gamma(i)*(Xn(:,i)-xav);
        end
        vs=v*s';
        Gamma=Gamma+vs+vs';
        XX0=Xn-x0;
        XX=XX0;
        eta=gamma.*(XX'*s);
        DDQs=Gamma*s;
        for i=1:m
            DDQs=DDQs+eta(i)*XX(:,i);
        end
        g=g+DDQs;
        c=Fn(opt);
        %%
        X=zeros(n+1,m);
        X(1,:)=ones(1,m);
        X(2:n+1,:)=XX;
        A=zeros(m,m);
        for i=1:m
            for j=1:m
                A(i,j)=0.5*(XX(:,i)'*XX(:,j))^2;
            end
        end
        W=zeros(m+n+1,m+n+1);
        W(1:m,1:m)=A;
        W(1:m,m+1:m+n+1)=X';
        W(m+1:m+n+1,1:m)=X;
        H=inv(W);
        %     [HE,W,Hreal] = HError();
        %     HE
        %     [E,en] = InterpolationError();
        %     E
        %%
        %新计算w,Hw以及beta
        w=zeros(m+n+1,1);
        dxx0=xnew-x0;
        for i=1:m
            w(i)=0.5*(((Xn(:,i)-x0)'*(dxx0))^2);
        end
        w(m+1)=1;
        w(m+2:m+n+1)=dxx0;
        Hw=H*w;
        beta=0.5*((dxx0'*dxx0)^2)-w'*Hw;
        
    else
        %%
        %不修改x0，可以直接调用H,w,Hw,beta,
    end
    %%
    %开始更新模型
    alpha=H(t,t);
    tau=Hw(t);
    % sigma=max([0,alpha])*max([0,beta])+tau^2;
    sigma=alpha*beta+tau^2;
    et=zeros(n+m+1,1);
    et(t)=1;
    eHw=et-Hw;
    H=H+(    alpha*(eHw*eHw')-beta*H(:,t)*H(t,:) +tau*(H(:,t)*eHw'+eHw*H(t,:) ) )/sigma;
    if (RATIO>=0.5) %%如果模型预测较好，只修改插值点的值
        C=Fnew-Qnew;
        Lcg=C*H(:,t);
    else %%如果模型不好，则修改所有插值点处的值
        R=zeros(m+n+1,1);
        
        for i=1:m
            if i~=t
                R(i)=Fn(i)-Q(Xn(:,i));
            else
                % R(i)=Fnew-Qnew;
                R(i)=Fnew-Q(xnew);
            end
        end
        Lcg=H*R;
    end
    lambda=Lcg(1:m);
    dc=Lcg(m+1);
    dg=Lcg(m+2:m+n+1);
    %%
    c=c+dc;
    g=g+dg;
    Gamma=Gamma+gamma(t)*XX(:,t)*XX(:,t)';
    for i=1:m
        if i~=t
            gamma(i)=gamma(i)+lambda(i);
        else
            gamma(i)=lambda(i);
        end
    end
    
    if Fnew<Fn(opt)
        opt=t;%更新最优的位置，并更新插值点
        Fn(t)=Fnew;
        Xn(:,t)=xnew;
    else
        %只更新点
        Fn(t)=Fnew;
        Xn(:,t)=xnew;
    end
    XX0=Xn-x0;
    
else
    %模型不变
end

%% 用于测试
% [HE] = HError()
%Test(F_times-m+1,:)=[Xn(:,opt)',Fn(opt)];
% [E,en] = InterpolationError();
%  E
% [HE,W,Hreal] = HError();
% HE
%%
%NEWUOAStep6();%To Setp6
flag=6;
end

function [] = NEWUOAStep6()
%NEWUOAStep6 对应NEWUOA算法的第六步

global RATIO flag
global Steps %测试用
Steps=[Steps;6];
if RATIO>=0.1
    %NEWUOAStep2();%To NEWUOAStep2
    flag=2;
else
    %NEWUOAStep7();%To NEWUOAStep7
    flag=7;
end

end

function [] = NEWUOAStep7()
%NEWUOAStep7 对应NEWUOA算法的第七步

global DIST MOVE opt Xn n m XXopt NXX flag
global Steps %测试用
Steps=[Steps;7];
DIST=0;
XXopt=zeros(n,m);
NXX=zeros(1,m);
for i=1:m
    XXopt(:,i)=Xn(:,i)-Xn(:,opt);
    NXX(i)=sqrt(XXopt(:,i)'*XXopt(:,i));
    if NXX(i)>DIST
        DIST=NXX(i);
        MOVE=i;
    end
end
%NEWUOAStep8();%To Setp8
flag=8;

end

function [] = NEWUOAStep8()
%NEWUOAStep8 对应NEWUOA算法的第八步

global DIST delta flag
global Steps %测试用
Steps=[Steps;8];
if DIST>=2*delta
    %NEWUOAStep9();%To NEWUOAStep9
    flag=9;
else
    %NEWUOAStep10();%To NEWUOAStep10
    flag=10;
end

end

function [] = NEWUOAStep9()
%NEWUOAStep9 对应NEWUOA算法的第九步
%BIGLAGandBIGDENOfNEWUOA 移动插值点的方法
%   目标函数为|l_t(xopt+d)|
%   约束条件为 ||d||<=deltah
%%
global F d NORMD beta w Hw rho delta DIST H Xn MOVE opt x0 n m RATIO Fnew Qnew F_times Krho XXopt NXX XX0 NFMAX flag
global Steps %测试用
Steps=[Steps;9];
%BIGLAG
deltah=max([min([0.1*DIST,0.5*delta]),rho]);
DELSQ=deltah*deltah;
N=n;
t=MOVE;
xt=Xn(:,t);
xopt=Xn(:,opt);
% X=Xn-x0;=XX0
%%
%l_t(x)=c+g*(x-x0)+0.5*(x-x0)'*G*(x-x0)
Lambda=diag(H(1:m,t));
c=H(m+1,t);
g=H(m+2:m+n+1,t);
G=XX0*Lambda*XX0';
%%
Dl=g+G*XX0(:,opt);%Dlt(xopt)=g+G*(xopt-x0);
eta=@(u)(H(1:m,t).*(XX0'*u));
DDlu=@(u)( XX0*eta(u) );
xx=XXopt(:,t);%xx=xt-xopt;
d=deltah*(xx/norm(xx));
%d2=-d1;
GD=DDlu(d);
DD=d'*d;
GG=Dl'*Dl;
SP=d'*Dl;
DHD=d'*GD;
if  (SP*DHD <= 0)
    d=-d;
    GD=-GD;
    lxd=abs(-SP+0.5*DHD);
else
    lxd=abs(SP+0.5*DHD);
end
TEMP=0;
if(SP*SP >= 0.99*DD*GG)||(GG*DELSQ <= 0.01*lxd*lxd) %%是否平行
    TEMP=1;
end
S=Dl+TEMP*GD;
K1=0;%计算次数
CTH=zeros(49,1);
STH=zeros(49,1);
ANGLE=zeros(49,1);
dPi=2*pi/50;
for i=1:49
    ANGLE(i)=dPi*i;
    CTH(i)=cos(ANGLE(i));
    STH(i)=sin(ANGLE(i));
end
while K1<=N
    K1=K1+1;
    DD=d'*d;
    SP=d'*S;
    SS=S'*S;
    TEMP=DD*SS-SP*SP;
    if (TEMP <= 1.0D-8*DD*SS)
        break;
    end
    DENOM=sqrt(TEMP);
    S=(DD*S-SP*d)/DENOM;%更新S
    W=DDlu(S);
    CF1=0.5*S'*W;
    CF2=d'*Dl;
    CF3=S'*Dl;
    CF4=0.5*d'*GD-CF1;
    CF5=S'*GD;
    TAUBEG=CF1+CF2+CF4;
    TAUMAX=TAUBEG;
    TAUOLD=TAUBEG;
    ISAVE=0;
    for i=1:49
        TAU=CF1+(CF2+CF4*CTH(i))*CTH(i)+(CF3+CF5*CTH(i))*STH(i);
        if (abs(TAU) >= abs(TAUMAX))
            TAUMAX=TAU;
            ISAVE=i;
            TEMPA=TAUOLD;
        else
            if ( i ==(ISAVE+1) )
                TEMPB=TAU;
            end
        end
        TAUOLD=TAU;
    end
    if (ISAVE == 0)
        TEMPA=TAU;
    end
    if (ISAVE == 49)
        TEMPB=TAUBEG;
    end
    STEP=0;
    if (TEMPA ~= TEMPB)
        TEMPA=TEMPA-TAUMAX;
        TEMPB=TEMPB-TAUMAX;
        STEP=0.5*(TEMPA-TEMPB)/(TEMPA+TEMPB);
    end
    angle=dPi*(ISAVE+STEP);
    cosA=cos(angle);
    sinA=sin(angle);
    TAU=CF1+(CF2+CF4*cosA)*cosA+(CF3+CF5*cosA)*sinA;
    d=cosA*d+sinA*S;
    GD=cosA*GD+sinA*W;
    S=Dl+GD;
    if (abs(TAU) <= 1.1D0*abs(TAUBEG))
        break;
    end
end
%%
NORMD=deltah;
xnew=xopt+d;
alpha=H(t,t);
w=zeros(m+n+1,1);
xx0=xnew-x0;
for i=1:m
    w(i)=0.5*((XX0(:,i)'*xx0)^2);
end
w(m+1)=1;
w(m+2:m+n+1)=xx0;
Hw=H*w;
beta=0.5*(norm(xx0)^4)-w'*Hw;
tau=Hw(t);
tau2=tau*tau;
%%
%BIGDEN
if abs(alpha*beta+tau2)<=0.8*tau2
    S=XXopt(:,t);
    DD=d'*d;
    DS=d'*S;
    SS=S'*S;
    v=zeros(m+n+1,1);
    for i=1:m
        v(i)=0.5*((XX0(:,i)'*XX0(:,opt))^2);
    end
    v(m+1)=1;
    v(m+2:m+n+1)=XX0(:,opt);
    NCN=0.5*((XX0(:,opt)'*XX0(:,opt))^2);
    DEN=@(d)(DENFunction(d,v,t,NCN));
    Weight=DS*DS/(DD*SS);
    Kt=t;
    if (Weight>= 0.99)
        for i=1:m
            if i~=opt
                DSTEMP=XXopt(:,i)'*d;
                SSTEMP=NXX(i)^2;
                Weightnew=DSTEMP*DSTEMP/(DD*SSTEMP);
                if Weightnew<=Weight
                    Kt=i;
                    Weight=Weightnew;
                    DS=DSTEMP;
                    SS=SSTEMP;
                end
            end
        end
        S=XXopt(:,Kt);
    end
    SSDEN=DD*SS-DS*DS;
    K2=0;
    TEMP=1/sqrt(SSDEN);
    S=TEMP*(DD*S-DS*d);
    while(K2<=N)
        DENTEST=zeros(1,50);
        DENTEST(1)=abs(DEN(d));
        for i=1:49
            dtest=CTH(i)*d+STH(i)*S;
            DENTEST(i+1)=abs(DEN(dtest));
        end
        [DENMAX,ISAVE]=max(DENTEST);
        if ISAVE==1
            TEMPA=DENTEST(50);
            TEMPB=DENTEST(2);
        else
            if ISAVE==50
                TEMPA=DENTEST(49);
                TEMPB=DENTEST(1);
            else
                TEMPA=DENTEST(ISAVE-1);
                TEMPB=DENTEST(ISAVE+1);
            end
        end
        if (TEMPA ~= TEMPB)
            TEMPA=TEMPA-DENMAX;
            TEMPB=TEMPB-DENMAX;
            STEP=0.5*(TEMPA-TEMPB)/(TEMPA+TEMPB);
        end
        angle=dPi*(ISAVE+STEP);
        d=cos(angle)*d+sin(angle)*S;
        [DENNEW,w,Hw,beta,S]=DENFALL(d,v,t,NCN);
        if K2==0
            DENOLD=abs(DENNEW);
        else
            if (abs(DENNEW) <= 1.1*DENOLD)
                DENOLD=abs(DENNEW);
                break;
            else
                DENOLD=abs(DENNEW);
            end
        end
        DD=d'*d;
        DS=d'*S;
        SS=S'*S;
        SSDEN=DD*SS-DS*DS;
        if (SSDEN >= 1.0D-8*DD*SS)
            K2=K2+1;
            TEMP=1/sqrt(SSDEN);
            S=TEMP*(DD*S-DS*d);
        else
            break;
        end
    end
end

%%
xopt=Xn(:,opt);
xnew=xopt+d;%新生成的点
F_times=F_times+1;
if F_times>NFMAX
    F_times=F_times-1;
    flag=16;
else
    Fnew=F(xnew);
    Qnew=Q(xnew);
    Krho=Krho+1;
    RATIO=1;
    %NEWUOAStep5();%To Setp5
    flag=5;
end
end
function [value] = Q(x)
%Q 插值函数
%测试用函数
global c g Gamma gamma x0 Xn m
G=Gamma;
for i=1:m
    G=G+gamma(i)*((Xn(:,i)-x0)*(Xn(:,i)-x0)');
end
value=c+(x-x0)'*g+0.5*(x-x0)'*G*(x-x0);

end
function DEN=DENFunction(d,v,t,C)
% C=0.5*||xopt-x0||^4
global Xn m n XX0 opt x0 H
Xnew=Xn(:,opt)+d;
xx0=Xnew-x0;
w=zeros(m+n+1,1);
for i=1:m
    w(i)=0.5*((XX0(:,i)'*xx0)^2);
end
w(m+1)=1;
w(m+2:m+n+1)=xx0;
alpha=H(t,t);
wv=w-v;
Hwv=H*wv;
DEN=alpha*(0.5*((xx0'*xx0)^2)-(XX0(:,opt)'*xx0)^2 + C  )-alpha*(wv'*Hwv)+Hwv(t);

end
function [DEN,w,Hw,beta,DDEN]=DENFALL(d,v,t,C)
% 返回[DEN,w,DDEN]
% C=0.5*||xopt-x0||^4
global Xn m n XX0 opt x0 H
Xnew=Xn(:,opt)+d;
xx0=Xnew-x0;
w=zeros(m+n+1,1);
for i=1:m
    w(i)=0.5*((XX0(:,i)'*xx0)^2);
end
w(m+1)=1;
w(m+2:m+n+1)=xx0;
alpha=H(t,t);
wv=w-v;
Hwv=H*wv;
NX0=(xx0'*xx0);
DEN=alpha*(0.5*(NX0^2)-(XX0(:,opt)'*xx0)^2 + C  )-alpha*(wv'*Hwv)+Hwv(t);

%%
%计算DDEN
Hw=H*w;
beta=0.5*(NX0^2)-w'*Hw;
tau=Hw(t);
eta1=Hwv(1:m);
eta2=Hwv(m+2:m+n+1);
DDEN2=zeros(n,1);
DDEN3=zeros(n,1);
for i=1:m
    DDEN2=DDEN2+((tau*H(t,i)-alpha*eta1(i))*(xx0'*XX0(:,i)))*XX0(:,i);
end
DDEN2=2*DDEN2;
for i=1:n
    DDEN3(i)=2*(tau*H(t,i+m+1)-alpha*eta2(i));
end
DDEN=2*alpha*(NX0*d+d'*xx0*XX0(:,opt))+(DDEN2)+(DDEN3);
end

function [] = NEWUOAStep10()
%NEWUOAStep10 对应NEWUOA算法的第十步

global delta RATIO NORMD rho flag
global Steps %测试用
Steps=[Steps;10];
M1=max(NORMD,delta);
if ( (M1-rho)/rho<10^-6 )&&( RATIO<=0  )
    %NEWUOAStep11();%To NEWUOAStep11
    flag=11;
else
    % NEWUOAStep2();%To NEWUOAStep2
    flag=2;
end

end

function [] = NEWUOAStep11()
%NEWUOAStep11 对应NEWUOA算法的第十一步

global rho rho_end flag
global Steps %测试用
Steps=[Steps;11];
if log10(rho)<=log10(rho_end)
    %NEWUOAStep13();%To NEWUOAStep13
    flag=13;
else
    %NEWUOAStep12();%To NEWUOAStep12
    flag=12;
end

end

function [] = NEWUOAStep12()
%NEWUOAStep12 对应NEWUOA算法的第十二步

global rho delta Krho flag
global Steps %测试用
Steps=[Steps;12];
delta=0.5*rho;
rho=rho/10;
Krho=0;
%NEWUOAStep2();%To Setp2
flag=2;
end

function [] = NEWUOAStep13()
%NEWUOAStep13 对应NEWUOA算法的第十三步

global d D3 rho Xn opt Fnew F F_times flag
global Steps %测试用
Steps=[Steps;13];
if D3(1)<0.5*rho
    xopt=Xn(:,opt);
    Fnew=F(xopt+d);
    F_times=F_times+1;
end
flag=16;
end

function [] = NEWUOAStep14()
%NEWUOAStep14 对应NEWUOA算法的第十四步

global CRVMIN D3 QF3 rho Krho flag
global Steps %测试用
Steps=[Steps;14];
if Krho>=3
    
    TEMP=0.125*CRVMIN*rho*rho;
    if max(QF3)<=TEMP&&max(D3)<=rho
        %NEWUOAStep11();%To Setp11
        flag=11;
    else
        %NEWUOAStep15();%To Setp15
        flag=15;
    end
else
    %NEWUOAStep15();%To Setp15
    flag=15;
end

end

function [] = NEWUOAStep15()
%NEWUOAStep15 对应NEWUOA算法的第十五步

global RATIO delta rho flag
global Steps %测试用
Steps=[Steps;15];
delta=0.1*delta;
RATIO=-1;
if delta<=1.5*rho
    delta=rho;
end
%NEWUOAStep7();%To Setp7
flag=7;
end

