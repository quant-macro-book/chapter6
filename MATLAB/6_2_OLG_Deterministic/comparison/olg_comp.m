
clear all; clc; close all;


ind=1; % GENERATE COMPARISON FIGURES
if ind==1
    
    a1=load('asset_1.txt');
    a2=load('asset_2.txt');
    a3=load('asset_3.txt');
    a4=load('asset_4.txt');

    c1=load('cons_1.txt');
    c2=load('cons_2.txt');
    c3=load('cons_3.txt');
    c4=load('cons_4.txt');
    
    norm=c1(1);

    ageA1=20:81;
    ageA2=20:86;
    
    
    ageC1=20:80;
    ageC2=20:85;


    minJ=20;
    maxJ=80;

    figure(3)
    hold on
    plot(ageA1,a1/norm,'-k', 'LineWidth',3)
    plot(ageA1,a2/norm,'-.k','LineWidth',3)
    plot(ageA1,a3/norm,':k', 'LineWidth',3)
    plot(ageA2,a4/norm,'--k', 'LineWidth',3)
    hold off
    legend('ベースライン','シナリオ 1','シナリオ 2','シナリオ 3','Location','NW')
    %legend('Baseline','Scenario 1','Scenario 2','Scenario 3','Location','NW')
    %title('Asset : a_{j}')
    xlim([minJ 86])
    xlabel('年齢')
    %ylim([0 12])
    grid on 
    box on 
    set(gca,'Fontsize',14,'FontName','Times New Roman');
    saveas(gcf,'fig_a_all.eps','epsc2')
    
    figure(4)
    hold on
    plot(ageC1,c1/norm,'- k', 'LineWidth',3)
    plot(ageC1,c2/norm,'-. k','LineWidth',3)
    plot(ageC1,c3/norm,': k', 'LineWidth',3)
    plot(ageC2,c4/norm,'-- k', 'LineWidth',3)
    hold off
    legend('ベースライン','シナリオ 1','シナリオ 2','シナリオ 3','Location','NW')
    %legend('Baseline','Scenario 1','Scenario 2','Scenario 3','Location','NW')
    %title('Asset : a_{j}')
    xlim([minJ 86])
    xlabel('年齢')
    %ylim([0 12])
    grid on 
    box on 
    set(gca,'Fontsize',14,'FontName','Times New Roman');
    saveas(gcf,'fig_c_all.eps','epsc2')
    
end