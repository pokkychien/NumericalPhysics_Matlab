clc;clear;close;
lambda = 0.00063; % (mm)
x1=linspace(-0.05,0.05,100);
x2=linspace(-50,50,100);
y1=linspace(-0.05,0.05,100);
y2=linspace(-50,50,100);
L=2000;
L_x2=length(x2);
L_y2=length(y2);
I=zeros([L_x2 L_y2]);
f=@(r) exp(1i*(pi/lambda*(r^2)/L));
for i=1:length(x1)
    for j=1:length(x2)
        for k=1:length(y1)
            for l=1:length(y2)
                x=x2(j)-x1(i);
                y=y2(l)-y1(k);
                r=sqrt(x^2+y^2);
                if (x1(i)^2+y1(k)^2<0.05^2)
                I(j,l)=I(j,l)+ f(r);
                end
            end
        end
    end
end
I=abs(I);
% plot(x2,I)
[A B]=meshgrid(x2,y2);
mesh(A,B,I)