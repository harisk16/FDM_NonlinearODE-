% ASR NODES 2020
% The non-linear mid-point method
% SOLVES y'(x)=q(x,y) 
% g(y(a),g(b)) = 0 

clear all
format long
clc

N = 200; % no of discritization points
mm = 4;
h = 1/N; % mesh size
F = 0.000001;
Kn = 0.5;  
y = (0:h:1)'; % nodes
M = zeros(mm*(N+1),mm*(N+1));
b = zeros(mm*(N+1),1);

Done = 0;
Max_iter = 50;

y0 = zeros(mm*(N+1),1);
for i = 1:N+1
   y0(mm*(i-1)+1) = 0;
   y0(mm*(i-1)+2) = 1;
   y0(mm*(i-1)+3) = 1;
   y0(mm*(i-1)+4) = 1;
end

iter = 0;

while (~Done)
    for i = 1:N
        yip1 = 0.5*(y(i)+y(i+1));
        Pip1 = 0.5*(y0(mm*(i-1)+1:mm*i)+y0(mm*i+1:mm*(i+1))); % dP/dy = A(y)Y + f(y)
        Vip1 = Pip1(1); % X velocity
        Sigmaip1 = Pip1(2); % Shear stress
        Thetaip1 = Pip1(3); % Temperature
        qip1 = Pip1(4); % Flux
        Aip1 = [0,-1/Kn,0,0;
                0,0,-F/((Thetaip1)^2),0;
                0,0,0,-4/(15*Kn);
                0,2*Sigmaip1/(Kn),0,0];
        Qip1 = [-Sigmaip1/Kn ; F/Thetaip1 ; -4*qip1/(15*Kn) ; (Sigmaip1)^2/Kn];
        fip1 = Qip1-Aip1*Pip1;
        M(mm*(i-1)+1:mm*i,mm*(i-1)+1:mm*i) = -(0.5*Aip1+eye(mm)/h);
        M(mm*(i-1)+1:mm*i,mm*i+1:mm*(i+1)) = (-0.5*Aip1+eye(mm)/h);
        b(mm*(i-1)+1:mm*i) = fip1;
    end 
    
    Ba = [sqrt(2/pi), 1, 0, 0;
        0, 0, 0, 0;
        0, 0,2*sqrt(2/pi), 1;
        0, 0, 0, 0];
    Bb = [0, 0, 0, 0;
        -sqrt(2/pi),1,0,0;
        0, 0, 0, 0;
        0, 0, -2*sqrt(2/pi), 1];
    Ya =[y0(1);y0(2);y0(3);y0(4)];
    Yb =[y0(N*mm+1);y0(N*mm+2);y0(N*mm+3);y0(N*mm+4)];
    g_at_a_b = [y0(2); y0(N*mm+2); y0(4); y0(N*mm+4)];
    M(mm*((N+1)-1)+1:mm*(N+1),mm*((N+1)-1)+1:mm*(N+1)) = Bb;
    M(mm*((N+1)-1)+1:mm*(N+1),1:mm) = Ba;
    alpha = Ba*Ya + Bb*Yb - g_at_a_b ;
    b(mm*((N+1)-1)+1:mm*(N+1)) = alpha;
    
    Y = M\b;
    
    error = max(abs(Y-y0));% max |y_i-y_i^*|
    iter = iter+1;
    disp(['Error in step ',num2str(iter),' is : ',num2str(error)])
    y0 = Y;

    Done = (iter>=Max_iter)||(error<0.000001);
    
    
end

figure(1);
y1 = zeros(N+1,1);
y2 = zeros(N+1,1);
for i = 1:N+1
    y1(i) = Y(mm*(i-1)+1);
    y2(i) = Y(mm*(i-1)+2);
    y3(i) = Y(mm*(i-1)+3);
    y4(i) = Y(mm*(i-1)+4);
end

figure(1)
plot(y1,y,'ok','LineWidth',1),grid on;
xlabel('$V(y)$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$y$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')

figure(2)
plot(y2,y,'ok','LineWidth',1),grid on;
xlabel('$Shear$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$y$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')

figure(3)
plot(y3,y,'ok','LineWidth',1),grid on;
xlabel('$Temperature$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$y$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')

figure(4)
plot(y4,y,'ok','LineWidth',1),grid on;
xlabel('$Heat Flux$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$y$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
