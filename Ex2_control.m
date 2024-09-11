clear; clc; clf;
t=[0];
dt=1e-4;
tau=0.2;
T=15;
h=0.0001;
t_total=[0:h:T];
x_total=zeros(4,length(t_total));
y_total=zeros(1,length(t_total));
y_total(1,1)=2.1678;
u=zeros(1,length(t_total));
u(1,1)=0.9333;
a=0.1;
mu1=3;mu2=1;I4=eye(4);eta1=0.6513;eta2=1.961;var=1.8;
k1=2.5;k2=2;m1=3.5;m2=1.2;c1=0.5;c2=3;d=1.3;
a1=(k1+k2)/m1;a2=(k2)/m1;a3=(c1+c2)/m1;a4=(c2)/m1;a5=(k2)/m2;a6=(c2)/m2;
b1=1/m1;b2=1/m2;
A=[0,0,1,0;0,0,0,1;-a1,a2,-a3,a4;a5,-a5,a6,-a6];B=[0,0,0,0;0,0,0,0;0,0,b1,0;0,0,0,b2];C=[0;0;b1;0];
D=[1,0,0,0;0,1,0,0;0,0,d,0;0,0,0,1];
P=[   1.992,  -0.1415,   0.7554, -0.2287;-0.1415,     1.67, -0.01754,  0.5689;  0.7554, -0.01754,    1.684, -0.1489; -0.2287,   0.5689,  -0.1489,    1.59];
H=[ -2.005, -4.391, -5.702, -9.21];


x_total(1,1)=1.8;
x_total(2,1)=0.8;
x_total(3,1)=-0.1;
x_total(4,1)=-0.9;
d=3.5;
m=5;
for k=1:1:7
     for i=2:length(t_total)
         t_now=(i-1)*h;
         M=round((t(k)+d)/h);
       if t_now>t(k)+d && x_total(:,i-1)'*P*x_total(:,i-1)>=exp(0.1)*x_total(:,M-1)'*P*x_total(:,M-1)
           t(k+1)=t_now;
       
       elseif t_now<t(k)+d
           t(k+1)=t(k)+m;
       end
       if t(k+1)>=t(k)+m
           t(k+1)=t(k)+m;
       end
     if t(k)<=t_now&&t_now<t(k)+d
         taut=t_now-tau;
           if taut<=0
               M1=0;
           elseif taut>0
           M1=round(taut/h);
           end
         x_total(:,i)=x_total(:,i-1)+(A*x_total(:,i-1)+B*[0;0;sin(a*x_total(1,i-1));sin(a*x_total(2,M1+1))]+C*H*x_total(:,i-1))*h;
         y_total(1,i)=sqrt(x_total(1,i)^2+x_total(2,i)^2+x_total(3,i)^2+x_total(4,i)^2);
         u(1,i)=H*x_total(:,i);
             end
  if abs(t_now-t(k))<=1e-6
      x_total(:,i)=D*x_total(:,i-1);
  
  y_total(1,i)=sqrt(x_total(1,i)^2+x_total(2,i)^2+x_total(3,i)^2+x_total(4,i)^2);
  u(1,i)=H*x_total(:,i);
   end
       
       if t(k)+d<=t_now&&t_now<t(k+1)
           taut=t_now-tau;
           if taut<=0
               M1=0;
           elseif taut>0
           M1=round(taut/h);
           end
           x_total(:,i)=x_total(:,i-1)+(A*x_total(:,i-1)+B*[0;0;sin(a*x_total(1,i-1));sin(a*x_total(2,M1+1))])*h;
         y_total(1,i)=sqrt(x_total(1,i)^2+x_total(2,i)^2+x_total(3,i)^2+x_total(4,i)^2);
         u(1,i)=H*x_total(:,i);
   
        end
  if abs(t_now-t(k))<=1e-6
  x_total(:,i)=D*x_total(:,i-1);
  
  y_total(1,i)=sqrt(x_total(1,i)^2+x_total(2,i)^2+x_total(3,i)^2+x_total(4,i)^2);
  u(1,i)=H*x_total(:,i);
  end  
     y1(1,i)=0.45*t_now/t_now;
     
end
end

for i=1:6
    t1(i)=t(i+1);
y(1,i)=0.5*(t1(i)/t1(i));
end
       
    
    

     
   m=-5:0.01:20; l=-5:0.01:20;
   h=1.85*m;   j=4.8.*l./l;
   n=h./m;     
a=-5:0.01:20;
b=6.*a./a;

 plot(t_total,x_total(1,:),t_total,x_total(2,:),t_total,x_total(3,:),t_total,x_total(4,:),t_total,y_total(1,:),'g',t1,y(1,i),'x k',b,a,'--',m,n,'--',j,l,'--');%,t_total,y1(1,i),'--'
 xlim([0,10]);ylim([-2,3]);
 legend('\iota_1(t)','\iota_2(t)','\iota_3(t)','\iota_4(t)','||\iota(t)||','t_k');
xlabel('t'),ylabel('\iota_1(t),\iota_2(t),\iota_3(t),\iota_4(t),||\iota(t)||')
%title('(a)')
annotation('arrow',[0.4 0.5],[0.6,0.6])
 str='$$\mathcal{T}-\pi=4.8$$';
 text(3.1,1.2,str,'interpreter','latex')
 annotation('arrow',[0.7 0.6],[0.3,0.3])
 str='$$\mathcal{T}=6$$';
 text(6.3,-0.6,str,'interpreter','latex')
 annotation('arrow',[0.26 0.26],[0.62,0.72])
  text(1.9,1.3,'\rho=1.85')
 figure
  plot(t_total,u(1,:),'m',t1,y(1,i),'x');
  xlim([0,10]);ylim([-1,5]);
  legend('u(t)','t_k');
xlabel('t'),ylabel('u(t)')
%title('(c)')