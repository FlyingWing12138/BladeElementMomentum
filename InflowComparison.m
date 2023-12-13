clear
clc
close all

data.NASA = load('Data\lambdaMeasuredNASA.mat');
data.Coleman = load('Data\ColemanDynamics.mat');
data.Drees = load('Data\DreesDynamics.mat');
data.Payne = load('Data\PayneDynamics.mat');
data.WhiteBlake = load('Data\WhiteBlakeDynamics.mat');
data.PittPeters = load('Data\PittPetersDynamics.mat');
data.Howlett = load('Data\HowlettDynamics.mat');
data.Uniform = load('Data\UniformDynamics.mat');

modelCell = {'NASA','Coleman','Drees','Payne','WhiteBlake','PittPeters','Howlett','Uniform'};
lineStyles = {'.','-','-','-','-','-','-','-'};
plotColor = {'#FF0000','#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#000000'};

%% Inflow Model Plotting

h = figure;
hold on
box on
grid on
set(gca, 'YDir', 'reverse')

angle1 = pi/2;
angle2 = pi+angle1;

% plot(data.NASA.r, data.NASA.lambdaMeasured(:,4),'.r','MarkerSize',8)

for i = 1:size(modelCell,2)

    model = modelCell{i};
    pos = find(abs(data.(model).psi-angle1)<0.001);
    plot(data.(model).r, data.(model).lambda(:,pos),lineStyles{i},'Color',plotColor{i},'LineWidth',1.5)

end

% plot(-data.NASA.r, data.NASA.lambdaMeasured(:,6),'.r','MarkerSize',8)

for i = 1:size(modelCell,2)

    model = modelCell{i};
    pos = find(abs(data.(model).psi-angle2)<0.001);
    if size(data.(model).lambda(:,pos),2)==0
        continue
    else
    plot(-data.(model).r, data.(model).lambda(:,pos),lineStyles{i},'Color',plotColor{i},'LineWidth',1.5)
    end

end

title('Inflow Model Comparison with Experiment Data', 'FontSize',12,'FontWeight','bold')
xlabel('Nondimensional Rotor Position')
ylabel('Inflow Ratio')
axis([-1.26 1.26 -0.025 0.085])
xticks(-1.2:0.4:1.2)
yticks(-0.02:0.02:0.08)
% legend([{'NASA Measurement'},modelCell],'Location','southwest')
legend(modelCell,'Location','best')

if angle1 == 0
    figureName = 'InflowComparisonLongitudinal';
elseif angle1 == pi/2
    figureName = 'InflowComparisonLateral';
else
    figureName = ['InflowComparison',num2str(angle1)];
end

saveas(h, ['Figures\',figureName],'fig')
saveas(h, ['Figures\',figureName],'png')
