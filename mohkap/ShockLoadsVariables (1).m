%% Shock loading calcs an ddynamic response of the rocket body and parachute system after deployment


%% 
clear
clc
close all

%% Define variables and set up the simulation
 

Shock_Cord__Widths = 0.01:0.005:0.06;
%%
Lengths = 2:0.1:15;
% %for i = 1:length(Lengths)
%     for j = 1:length(Shock_Cord__Widths)
%         L = 12;
%         L_sc = Shock_Cord__Widths(j);%0.01524; %0.015875;
%         Es = 2.7e9;
%         rho = 0.9161;
%         g = 9.81;
%         Dp = 6.09;

%         A_p = pi/4 * (Dp^2);
%         A_b = 0.01824;
%         Cdp = 0.97;
%         Cdb = 0.54;
%         mp = 1.355;
%         mb = 38;
%         V_0 = 35;
%         D_sc = 0.0127;
%         W_sc = 0.002;
%         A_sc = L_sc*W_sc;
%         %A_sc = pi*D_sc^2/4;
%         k = Es*A_sc/L;
%         out = sim("PostDeploymentResponse.slx", 20);
%         maxShock(j) = max(out.ShockLoad.Data);
%     %end
%     end
% for i = 1:length(Shock_Cord__Widths)
%     for j = 1:length(Lengths)
        L = 8; %length of shockcord
        L_sc = 0.015875; %flat length of shockcord
        Es = 2.7e9; %youngs modulus
        height = 300; %height in metres
        [~,~,~,rho] = atmoscoesa(height);
        g = 9.81; 
        Dp = ; %diameter of parachute
        A_p = pi/4 * (Dp^2); %area of parachute
        %Db = 0.018; 
        %A_b = pi/4 * (Db^2); %area of body
        A_b = 0.018;
        Cdp = 1.75; %coefficent of drag of parachute
        Cdb = 0.5; %coefficient of drag of body
        density_nylon = 0.06;
        mp = density_nylon * A_p; %mass of parachute
        mb = 0.9; %mass of body
        V_0 = 10; %deployment speed
        D_sc = 0.04; %diameter of shockcord (not use now)
        W_sc = 0.002; %width of shockcord
        A_sc = L_sc*W_sc; %cross sectional area of shockcord
        %A_sc = pi*D_sc^2/4;
        k = Es*A_sc/L; %stiffness ratio
        simulation_time = 10; %adjust until converges
        out = sim("PostDeploymentResponse.slx", simulation_time);
        %%
        loads = out.ShockLoad.Data;
        
        t = out.tout;
        maxShock = max(out.ShockLoad.Data);
        figure()

        plot(t, out.ShockLoad.Data, LineWidth=2)
        xlabel("Time (s)", "Interpreter","latex", "FontSize", 20)
        ylabel("Load (N)", "Interpreter","latex", "FontSize", 20)
        set(gca,'FontSize', 20, "TickLabelInterpreter", "Latex")
        grid on
        grid minor
        %end
%     end
%     filename = strcat('Shocks_', num2str(L_sc), 'm_vary_lengths.mat');
%     save(filename, 'maxShock', 'Lengths', 'L_sc')
% end
% %% Not needed for sims for individual event - DELETE IF NOT NEEDED
% figure()
% for i = [2, 4, 6, 8, 10, 12]
%     filename = strcat('Shocks_', int2str(i), 'm_vary_L_sc.mat');
%     load(filename)
%     plot(Shock_Cord__Widths, maxShock, LineWidth=2)
%     xlabel("Shock cord width (m)", "Interpreter","latex", "FontSize", 20)
%     ylabel("Shock load (N)", "Interpreter","latex", "FontSize", 20)
%     set(gca,'FontSize', 20, "TickLabelInterpreter", "Latex")
% 
%     hold on
% end
% legend("L = 2", "L = 4", "L = 6", "L = 8", "L = 10", "L = 12", "Interpreter","latex", "FontSize", 20, "location", "northwest")
% Shock_Cord__Widths = 0.01:0.005:0.06;
% figure()
% for i = 1:length(Shock_Cord__Widths)-4
%     filename = strcat('Shocks_', num2str(Shock_Cord__Widths(i)), 'm_vary_lengths.mat');
%     load(filename)
%     plot(Lengths, maxShock, LineWidth=2)
%     xlabel("Shock cord length (m)", "Interpreter","latex", "FontSize", 20)
%     ylabel("Shock load (N)", "Interpreter","latex", "FontSize", 20)
%     set(gca,'FontSize', 20, "TickLabelInterpreter", "Latex")
%     xlim([2 15]);
%     hold on
% end
% for i = 8:length(Shock_Cord__Widths)
%     filename = strcat('Shocks_', num2str(Shock_Cord__Widths(i)), 'm_vary_lengths.mat');
%     load(filename)
%     plot(Lengths, maxShock, "--", LineWidth=2)
%     xlabel("Shock cord length (m)", "Interpreter","latex", "FontSize", 20)
%     ylabel("Shock load (N)", "Interpreter","latex", "FontSize", 20)
%     set(gca,'FontSize', 20, "TickLabelInterpreter", "Latex")
%     xlim([2 15]);
%     hold on
% end
% legendStrings = "Cord Width = " + string(Shock_Cord__Widths);
% legend(legendStrings, "Interpreter","latex", "FontSize", 20, "location", "northeast")
% 
% %% Plot surface
% 
% Shock_Cord__Widths = 0.01:0.001:0.06;
% Lengths = 2:0.1:15;
% 
% for i = 1:length(Lengths)
%     for j = 1:length(Shock_Cord__Widths)
%         L_sc = Shock_Cord__Widths(j);
%         L = Lengths(i);
%         Es = 2.7e9;
%         rho = 0.9161;
%         g = 9.81;
%         Dp = 6.09;
%         A_p = pi/4 * (Dp^2);
%         A_b = 0.01824;
%         Cdp = 0.97;
%         Cdb = 0.54;
%         mp = 1.355;
%         mb = 38;
%         V_0 = 35;
%         D_sc = 0.0127;
%         W_sc = 0.002;
%         A_sc = L_sc*W_sc;
%         %A_sc = pi*D_sc^2/4;
%         k = Es*A_sc/L;
%         out = sim("PostDeploymentResponse.slx", 20);
%         maxShock(i,j) = max(out.ShockLoad.Data);
% 
%     end
% end
% %% Plot the Data
% [Ls, W] = meshgrid(Lengths, Shock_Cord__Widths);
% figure()
% surf(Ls, W, maxShock')
% 
% xslice = [];
% yslice = [];
% zslice = 6800;
% figure()
% 
% %% Plot Drogue and main plots
% 
% load('main_lower_deploy.mat')
% figure()
% plot(t, loads, "LineWidth", 2)
% hold on 
% load('drogue_lower_deploy.mat')
% plot(t, loads, "LineWidth", 2)
% load('drogue_nose_deploy.mat')
% plot(t, loads, "LineWidth", 2)
% xlim([0 8])
% grid on
% grid minor
% xlabel("Time (s)", "Interpreter","latex", "FontSize", 20)
% ylabel("Force (N)", "Interpreter","latex", "FontSize", 20)
% set(gca,'FontSize', 20, "TickLabelInterpreter", "Latex")
% legend("Load on lower body - main deploy", "Load on lower body - drogue deploy", "Load on nose - drogue deploy", "Interpreter","latex", "FontSize", 20, "location", "northeast")