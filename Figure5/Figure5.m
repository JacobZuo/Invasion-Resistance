clear
clc

% Parameters of the dynamics
NumberSpecies = 15; % number of samples
steps = 5000;
dt = 0.01;

% Set seed for reproducibility
rng(123);
GrowthRate = normrnd(1,0.2,1,15); % growth rate
GrowthRate(1) = 1;


MeanAAll = [-0.1,-0.3,-0.5];


for MeanAIndex = 1:3
    MeanA = MeanAAll(MeanAIndex); % mean of pairwise matrix A
    A_STD = 0.1; % standard deviation of pairwise matrix A

    rng(123);
    A = normrnd(MeanA, A_STD, NumberSpecies, NumberSpecies); % Initial interaction matrix

    for i = 1:NumberSpecies
        A(i, i) = -1;
    end


    Collection = [];
    Abundance = [];

    FinalAbundance=[];
    for SynComSize = 1:14
        Combination = nchoosek(2:NumberSpecies, SynComSize);
        if size(Combination, 1) > 50
            Combination = Combination(randsample(1:size(Combination, 1), 50),:);
        end
        for j = 1:size(Combination, 1)
            Abundance = zeros(steps+1, NumberSpecies);
            Abundance(1, :) = 0;
            % rng(SynComSize);
            Abundance(1, [1 Combination(j,:)]) = 0.01;

            for Time = 1:steps
                Abundance(Abundance < 0) = 0; % Ensure abundance cannot be negative
                pair_term = sum(A .* Abundance(Time, :),2)'; % Pairwise interaction term
                Abundance(Time + 1, :) = Abundance(Time, :) + dt * Abundance(Time, :) .* (GrowthRate + pair_term);
            end

            FinalAbundance = [FinalAbundance; Abundance(end, 1)]; % Final abundance
            Collection_sub = zeros(1, NumberSpecies);
            Collection_sub(Combination(j,:)) = 1;
            Collection = [Collection; Collection_sub]; % Species combination
        end
    end


    MatrixATable = A;
    AbundanceGroup{MeanAIndex} = cell(0);


    for i = 1:15
        AbundanceGroup{MeanAIndex}{i}=FinalAbundance(sum(Collection,2)==i);
    end
end

BoxColor2=[0.6:(-0.45/20):0.15;0.8:(-0.65/20):0.15;1:(-0.45/20):0.55]';


figure('Position',[100,100,1000,750],'Color',[1,1,1])
hold on

xOffsets = [-0.25, 0, 0.25];
GroupColors = [0.6, 0.75, 1.0; 0.45, 0.6, 0.9; 0.3, 0.45, 0.75]; 

for MeanAIndex = 1:3
    for i = 1:15
        xPos = ones(size(AbundanceGroup{MeanAIndex}{i})) * (i + xOffsets(MeanAIndex));
        boxchart(xPos, log10(AbundanceGroup{MeanAIndex}{i}), 'BoxFaceColor', GroupColors(MeanAIndex,:), 'MarkerColor', GroupColors(MeanAIndex,:), 'MarkerSize', 4, 'LineWidth', 1.5, 'BoxEdgeColor', [0, 0, 0], 'WhiskerLineColor', [0, 0, 0], 'BoxWidth', 0.18, 'BoxFaceAlpha', 1);
    end
end

axis([0 15 -0.9 0.1])
set(gca,'FontSize',24,'LineWidth',1.5)
set(gca,'layer','top')
xticks(1:4:14)
xlabel('SynCom Strain Number')
ylabel('log_{10} Biomass(a.u.)')
box on
print('Figure5A.pdf','-dpdf')


save('Figure5A.mat')
load('Figure5A.mat')


MeanAAll = 0.1:-0.025:-0.5;
A_STDAll = 0.02:0.02:0.4;


for  MeanAIndex=1:size(MeanAAll,2) % mean of pairwise matrix A
    for A_STDIndex=1:size(A_STDAll,2) % deviation of pairwise matrix A
        for TestIndex=1:3
        MeanA = MeanAAll(MeanAIndex); % mean of pairwise matrix A
        A_STD = A_STDAll(A_STDIndex); % standard deviation of pairwise matrix A
        FileName=['Simulation\Results\gLV-',num2str(MeanA),'-',num2str(A_STD),'-',num2str(TestIndex),'.mat'];
        load(FileName,'R2La2','Lasso2ParaNum')
        R2La2All(MeanAIndex,A_STDIndex,TestIndex)=R2La2;
        Lasso2ParaNumAll(MeanAIndex,A_STDIndex,TestIndex)=Lasso2ParaNum;
        end
    end
end


figure('Position',[100,100,800,700],'Color',[1,1,1])
imagesc(mean(R2La2All,3)',[0,1])
colorbar off
axis equal
axis([0.5 25.5 0.5 20.5])
yticks(5:5:20)
xticks(1:4:25)
yticklabels(A_STDAll(5:5:20))
xticklabels(MeanAAll(1:4:25))
colorbar
set(gca,'YDir','normal')
xlabel('')
ylabel('')
box on
set(gca,'FontSize',20,'LineWidth',1.5)
set(gcf,'PaperType','A3')
set(gca,'layer','top')
set(gcf,'Renderer','painters')
colorbar
print('Figure5B.pdf','-dpdf')



LassoOrder1=Collection(15:end,2:end);
LassoOrder1(LassoOrder1==0)=-1;

LassoOrder2=[];

for i = 1:size(LassoOrder1,1)
    Order2Matrix=nchoosek(LassoOrder1(i,:),2);
    LassoOrder2(i,:)=[LassoOrder1(i,:),(Order2Matrix(:,1).*Order2Matrix(:,2))'];
end

close all

figure('Position',[100,100,1050,1000],'Color',[1,1,1])
hold on

for MeanAIndex = 1:3
FinalAbundance=cell2mat(AbundanceGroup{MeanAIndex}');


LassoOrderY=log10(FinalAbundance(15:end));

[B1,FitInfo1]=lasso(LassoOrder1, LassoOrderY,'cv',10);
[B2,FitInfo2]=lasso(LassoOrder2, LassoOrderY,'cv',10);

R2La1=fitlm(sum(LassoOrder1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),LassoOrderY).Rsquared.Ordinary;
R2La2=fitlm(sum(LassoOrder2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),LassoOrderY).Rsquared.Ordinary;

plot(sum(LassoOrder1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),LassoOrderY,'o','MarkerSize',10,'LineWidth',1.5,'Color',GroupColors(MeanAIndex,:))

end
plot([-3.5 0.5],[-3.5 0.5],'--','Color',[0.5 0.5 0.5])

axis([-1 0.1 -1 0.1])

axis equal
box on
axis([-1 0.1 -1 0.1])
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Predicted by Regression (log 10)')
ylabel('gLV model (log 10)')
print('FigureS13A.pdf','-dpdf','-r300')



figure('Position',[100,100,1050,1000],'Color',[1,1,1])
hold on
for MeanAIndex = 1:3
FinalAbundance=cell2mat(AbundanceGroup{MeanAIndex}');


LassoOrderY=log10(FinalAbundance(15:end));

[B1,FitInfo1]=lasso(LassoOrder1, LassoOrderY,'cv',10);
[CorrLa1,PLa1]=corr(sum(LassoOrder1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),LassoOrderY);

[B2,FitInfo2]=lasso(LassoOrder2, LassoOrderY,'cv',10);
[CorrLa2,PLa2]=corr(sum(LassoOrder2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),LassoOrderY);


R2La1=fitlm(sum(LassoOrder1.*B1(:,FitInfo1.Index1SE)',2)+FitInfo1.Intercept(FitInfo1.Index1SE),LassoOrderY).Rsquared.Ordinary;
R2La2=fitlm(sum(LassoOrder2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),LassoOrderY).Rsquared.Ordinary;

plot(sum(LassoOrder2.*B2(:,FitInfo2.Index1SE)',2)+FitInfo2.Intercept(FitInfo2.Index1SE),LassoOrderY,'o','MarkerSize',10,'LineWidth',1.5,'Color',GroupColors(MeanAIndex,:))
end
plot([-3.5 0.5],[-3.5 0.5],'--','Color',[0.5 0.5 0.5])
box on
axis([-1 0.1 -1 0.1])

axis equal
set(gca,'FontSize',28,'LineWidth',1.5)
set(gca,'layer','top')
xlabel('Predicted by Regression (log 10)')
ylabel('gLV model (log 10)')

text(-3.3,0.2,['R^2 = ',num2str(R2La2,3)],'FontSize',24)
print('FigureS13B.pdf','-dpdf','-r300')



IndexOrder=nchoosek(1:14,2);
Order2CorrIndex=zeros(14);
Order2CorrValue=zeros(14);

NonZeroOrder2Index=IndexOrder(B2(15:end,FitInfo2.Index1SE)'~=0,:);
NonZeroOrder2Index2=find(B2(15:end,FitInfo2.Index1SE)'~=0);


for j=1:size(NonZeroOrder2Index,1)
    Order2CorrIndex(NonZeroOrder2Index(j,1),NonZeroOrder2Index(j,2))=Order2CorrIndex(NonZeroOrder2Index(j,1),NonZeroOrder2Index(j,2))+1;
    Order2CorrValue(NonZeroOrder2Index(j,1),NonZeroOrder2Index(j,2))=B2(14+NonZeroOrder2Index2(j),FitInfo2.Index1SE);

end

cmap=[];
cmap(1:128,1)=0.85:0.15/127:1;
cmap(1:128,2)=0.23:0.77/127:1;
cmap(1:128,3)=0.10:0.90/127:1;
cmap(129:256,1)=1:-1.00/127:0.00;
cmap(129:256,2)=1:-0.65/127:0.35;
cmap(129:256,3)=1:-0.15/127:0.85;




close all

figure('Position',[100,100,1050,1000],'Color',[1,1,1])
hold on
imagesc(A,[-1 1])
colorbar off
axis equal
axis([0.5 15.5 0.5 15.5])
colormap(cmap)
xticks(1:15)
yticks(1:15)
box on
set(gca,'YDir','reverse')
set(gca,'yaxislocation','right')
set(gca,'xaxislocation','top')

colorbar
for i=1:16
    plot([i-0.5,i-0.5],[0.5,15.5],'-k','LineWidth',1)
    plot([0.5,15.5],[i-0.5,i-0.5],'-k','LineWidth',1)
end

set(gca,'FontSize',16,'LineWidth',1.5)
box off
set(gcf,'PaperType','A3')
set(gca,'layer','top')
set(gcf,'Renderer','painters')
colorbar
print('FigureS15A.pdf','-dpdf')
close all


figure('Position',[100,100,1050,1000],'Color',[1,1,1])
hold on
imagesc(Order2CorrValue,[-0.01 0.01])
colorbar off
axis equal
axis([0.5 14.5 0.5 14.5])
colormap(cmap)
xticks(1:14)
yticks(1:14)
box on
set(gca,'YDir','reverse')
set(gca,'yaxislocation','right')
set(gca,'xaxislocation','top')

for i=1:13
    plot([i+0.5,i+0.5],[0.5,i+1.5],'-k','LineWidth',1)
    plot([i-0.5,14.5],[i+0.5,i+0.5],'-k','LineWidth',1)
end
plot([0.5,0.5],[0.5,1.5],'-k','LineWidth',1)
plot([13.5,14.5],[14.5,14.5],'-k','LineWidth',1)

set(gca,'FontSize',16,'LineWidth',1.5)
box off
set(gcf,'PaperType','A3')
set(gca,'layer','top')
set(gcf,'Renderer','painters')
colorbar
print('Figure15B.pdf','-dpdf')