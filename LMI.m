setlmis([]);
Q=lmivar(1,[4,1]);
P=lmivar(1,[4,1]);
G=lmivar(2,[1,4]);
alpha1=lmivar(1,[1,0]);
alpha2=lmivar(1,[1,0]);



mu1=3;mu2=1;I4=eye(4);eta1=0.6513;eta2=1.961;var=1.8;
k1=2.5;k2=2;m1=3.5;m2=1.2;c1=0.5;c2=3;d=1.3;
a1=(k1+k2)/m1;a2=(k2)/m1;a3=(c1+c2)/m1;a4=(c2)/m1;a5=(k2)/m2;a6=(c2)/m2;
b1=1/m1;b2=1/m2;
A=[0,0,1,0;0,0,0,1;-a1,a2,-a3,a4;a5,-a5,a6,-a6];B=[0,0,0,0;0,0,0,0;0,0,b1,0;0,0,0,b2];C=[0;0;b1;0];
D=[1,0,0,0;0,1,0,0;0,0,d,0;0,0,0,1];
lmiterm([-1 1 1 Q],1,1);
lmiterm([1 1 1 0],1/mu1*I4);
lmiterm([2 1 1 Q],1,1);
lmiterm([-2 1 1 0],1/mu2*I4);
lmiterm([3 1 1 Q],1,-var);
lmiterm([3 1 2 Q],1,D');
lmiterm([3 2 2 Q],-1,1);
lmiterm([4 1 1 Q],A,1,'s');
lmiterm([4 1 1 G],C,1,'s');
lmiterm([4 1 1 Q],eta1,1);
lmiterm([4 1 2 0],B);
lmiterm([4 2 2 alpha1],-1,1);
lmiterm([5 1 1 Q],A,1,'s');
lmiterm([5 1 1 Q],-eta2,1);
lmiterm([5 1 2 0],B);
lmiterm([5 2 2 alpha2],-1,1);
lmiterm([-6 1 1 Q],1,1);
lmiterm([7 1 1 alpha1],1,1);
lmiterm([-7 1 1 0],1);
lmiterm([7 1 1 alpha2],1,1);
lmiterm([-7 1 1 0],3);






  
  
  
  lmisys = getlmis; %这都是范式了，获取上边描述的LMI系统
[tmin,xfeas] = feasp(lmisys); %验证 LMI 的可行性
if tmin<0
Q=dec2mat(lmisys,xfeas,Q)
P=vpa(inv(Q),4)
G=dec2mat(lmisys,xfeas,G)
H=vpa(G*P,4)

alpha1=dec2mat(lmisys,xfeas,alpha1)
alpha2=dec2mat(lmisys,xfeas,alpha2)
elseif tmin>=0
    Q=NAN;
end
