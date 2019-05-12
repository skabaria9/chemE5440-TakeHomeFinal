%Sneha Kabaria, srk97
%Take Home Final - ChemE 5440 
%Problem 2 

%-------------------------------------------------------------------------
%Part B - Plot nullclines
%-------------------------------------------------------------------------

%Define Variables
alpha = 10;  
syms u v n

%Define equations
f(u,v,n) = alpha/(1+v^n)-u;
g(u,v,n) = alpha/(1+u^n)-v;

%Define space to plot 
u_plot = linspace(-10,10,200);
u_plot_small = linspace(-20,20,200);
u_plot_pos = linspace(0.0001,10,200);
v_plot = linspace(-10,10,200);

%Plot nullclines when n = 1; 
n=1; 
f_nullclines_n1(u) = solve(f(u,v,n)==0,[v])
g_nullclines_n1(u) = solve(g(u,v,n)==0,[v])
figure(1); 
%Plot f nullclines
subplot(2,3,1); plot(u_plot,f_nullclines_n1(u_plot),'-k')
xlabel("u"); ylabel("v"); title("f nullcline, n=1");
%Plot g nullclines
subplot(2,3,2); plot(u_plot,g_nullclines_n1(u_plot),'-r')
xlabel("u"); ylabel("v"); title("g nullcline, n=1");
%Plot together 
subplot(2,3,3);plot(u_plot,f_nullclines_n1(u_plot),'-k',u_plot,g_nullclines_n1(u_plot),'-r');
xlabel("u"); ylabel("v"); title("nullclines, n=1");
axis([-5 5 -35 35])
xL = xlim;
yL = ylim;
line([0 0], yL);  %y-axis
line(xL, [0 0]);  %x-axis
legend("f nullcline", "g nullcline","x-axis","y-axis")
%Plot nullclines when n = 2; 
n=2; 
f_nullclines_n2(u) = solve(f(u,v,n)==0,[v])
g_nullclines_n2(u) = solve(g(u,v,n)==0,[v])
f_nullcline_n2_1(u) = (-(u - 10)/u)^(1/2); 
f_nullcline_n2_2(u) = -(-(u - 10)/u)^(1/2); 
figure(1); 
%Plot f nullclines
subplot(2,3,4); plot(u_plot_pos,f_nullcline_n2_1(u_plot_pos),'-k',u_plot_pos,f_nullcline_n2_2(u_plot_pos),'--k')
xlabel("u"); ylabel("v"); title("f nullcline, n=2");
legend("f nullcline 1", "f nullcline 2")
%Plot g nullclines
subplot(2,3,5); plot(u_plot_small,g_nullclines_n2(u_plot_small),'-r')
xlabel("u"); ylabel("v"); title("g nullcline, n=2");
%Plot together 
subplot(2,3,6);
plot(u_plot_pos,f_nullcline_n2_1(u_plot_pos),'-k',u_plot_pos,f_nullcline_n2_2(u_plot_pos),'--k',u_plot_small,g_nullclines_n2(u_plot_small),'-r')
xlabel("u"); ylabel("v"); title("nullclines, n=2");
axis([-5 15 -20 20])
xL = xlim;
yL = ylim;
line([0 0], yL);  %y-axis
line(xL, [0 0]);  %x-axis
legend("f nullcline 1","f nullcline 2", "g nullcline","x-axis","y-axis")

%Find Intersections in each case
%n =1
u_inter_n1 = double(solve(f_nullclines_n1(u) == g_nullclines_n1(u),u))
u_inter_n2_1 = double(solve(f_nullcline_n2_1(u) == g_nullclines_n2(u),u))
u_inter_n2_2 = double(solve(f_nullcline_n2_2(u) == g_nullclines_n2(u),u))


%-------------------------------------------------------------------------
%Part C - Plot StreamPlot
%-------------------------------------------------------------------------

%Define Functions with different n
%n = 1
f_n1(u,v) = f(u,v,1);
g_n1(u,v) = g(u,v,1);
%n = 2
g_n2(u,v) = f(u,v,2);
f_n2(u,v) = f(u,v,2);

%Calculate the Steady-States
%n=1
SS_n1 = solve(f_n1(u,v)==0, g_n1(u,v)==0,[u,v]);
u_SS_n1 = double(SS_n1.u)
v_SS_n1 = double(SS_n1.u)
%n=2
SS_n2 = solve(f_n2(u,v)==0, g_n2(u,v)==0,[u,v]);
u_SS_n2 = SS_n2.u
v_SS_n2 = SS_n2.v
%Also Calculate v's from previous u's found from nullclines
%only for n=2
u_inter_n2_1
v_inter_n2_1 = double(f_nullcline_n2_1(u_inter_n2_1))



