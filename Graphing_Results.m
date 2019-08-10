%% Load all the Fringe Length data
f = load('FringesMethaneNoDilution.mat');
FringesMethaneNoDilution = f.FringesMethaneNoDilution;

f = load('FringesMethane50Dilution.mat');
FringesMethane50Dilution = f.FringesMethaneNoDilution;

f = load('FringesMethane67Dilution.mat');
FringesMethane67Dilution = f.FringesMethaneNoDilution;

f = load('FringesEthyleneNoDilution.mat');
FringesEthyleneNoDilution = f.FringesMethaneNoDilution;

f = load('FringesEthylene67Dilution.mat');
FringesEthylene67Dilution = f.FringesMethaneNoDilution;

f = load('FringesEthylene90Dilution.mat');
FringesEthylene90Dilution = f.FringesMethaneNoDilution;

f = load('FringesEthaneNoDilution.mat');
FringesEthaneNoDilution = f.FringesMethaneNoDilution;

f = load('FringesEthane50Dilution.mat');
FringesEthane50Dilution = f.FringesMethaneNoDilution;

f = load('FringesEthane85Dilution.mat');
FringesEthane85Dilution = f.FringesMethaneNoDilution;

%% Load all the Fringe Tortuosity data
f = load('TortuosityMethaneNoDilution.mat');
TortuosityMethaneNoDilution = f.TortuosityMethaneNoDilution;

f = load('TortuosityMethane50Dilution.mat');
TortuosityMethane50Dilution = f.TortuosityMethaneNoDilution;

f = load('TortuosityMethane67Dilution.mat');
TortuosityMethane67Dilution = f.TortuosityMethaneNoDilution;

f = load('TortuosityEthyleneNoDilution.mat');
TortuosityEthyleneNoDilution = f.TortuosityMethaneNoDilution;

f = load('TortuosityEthylene67Dilution.mat');
TortuosityEthylene67Dilution = f.TortuosityMethaneNoDilution;

f = load('TortuosityEthylene90Dilution.mat');
TortuosityEthylene90Dilution = f.TortuosityMethaneNoDilution;

f = load('TortuosityEthaneNoDilution.mat');
TortuosityEthaneNoDilution = f.TortuosityMethaneNoDilution;

f = load('TortuosityEthane50Dilution.mat');
TortuosityEthane50Dilution = f.TortuosityMethaneNoDilution;

f = load('TortuosityEthane85Dilution.mat');
TortuosityEthane85Dilution = f.TortuosityMethaneNoDilution;

%% Load all the Separation Distance data
f = load('SeparationDistanceMethaneNoDilution.mat');
SeparationDistanceMethaneNoDilution = f.SeparationDistanceMethaneNoDilution;

f = load('SeparationDistanceMethane50Dilution.mat');
SeparationDistanceMethane50Dilution = f.SeparationDistanceMethaneNoDilution;

f = load('SeparationDistanceMethane67Dilution.mat');
SeparationDistanceMethane67Dilution = f.SeparationDistanceMethaneNoDilution;

f = load('SeparationDistanceEthyleneNoDilution.mat');
SeparationDistanceEthyleneNoDilution = f.SeparationDistanceMethaneNoDilution;

f = load('SeparationDistanceEthylene67Dilution.mat');
SeparationDistanceEthylene67Dilution = f.SeparationDistanceMethaneNoDilution;

f = load('SeparationDistanceEthylene90Dilution.mat');
SeparationDistanceEthylene90Dilution = f.SeparationDistanceMethaneNoDilution;

f = load('SeparationDistanceEthaneNoDilution.mat');
SeparationDistanceEthaneNoDilution = f.SeparationDistanceMethaneNoDilution;

f = load('SeparationDistanceEthane50Dilution.mat');
SeparationDistanceEthane50Dilution = f.SeparationDistanceMethaneNoDilution;

f = load('SeparationDistanceEthane85Dilution.mat');
SeparationDistanceEthane85Dilution = f.SeparationDistanceMethaneNoDilution;

%% Load all the Percent of Stacked Layers data
f = load('PercentStackedFringesMethaneNoDilution.mat');
PercentStackedFringesMethaneNoDilution = f.PSL;

f = load('PercentStackedFringesMethane50Dilution.mat');
PercentStackedFringesMethane50Dilution = f.PSL;

f = load('PercentStackedFringesMethane67Dilution.mat');
PercentStackedFringesMethane67Dilution = f.PSL;

f = load('PercentStackedFringesEthyleneNoDilution.mat');
PercentStackedFringesEthyleneNoDilution = f.PSL;

f = load('PercentStackedFringesEthylene67Dilution.mat');
PercentStackedFringesEthylene67Dilution = f.PSL;

f = load('PercentStackedFringesEthylene90Dilution.mat');
PercentStackedFringesEthylene90Dilution = f.PSL;

f = load('PercentStackedFringesEthaneNoDilution.mat');
PercentStackedFringesEthaneNoDilution = f.PSL;

f = load('PercentStackedFringesEthane50Dilution.mat');
PercentStackedFringesEthane50Dilution = f.PSL;

f = load('PercentStackedFringesEthane85Dilution.mat');
PercentStackedFringesEthane85Dilution = f.PSL;

%% Fringe Label -- character array with all of the labels

A = 'Methane No Dilution ';
B = 'Methane 50 Dilution ';
C = 'Methane 67 Dilution ';
D = 'Ethylene No Dilution';
E = 'Ethylene 67 Dilution';
F = 'Ethylene 90 Dilution';
G = ' Ethane No Dilution ';
H = ' Ethane 50 Dilution ';
I = ' Ethane 85 Dilution ';


%% 
FL=zeros(length(FringesMethaneNoDilution),length('Methane No Dilution '));
for i=1:length(FringesMethaneNoDilution)
    FL(i,:)= char(A) ;
end

for j=i:i+length(FringesMethane50Dilution)
    FL(j,:)= char(B) ;
end

for i=j:j+length(FringesMethane67Dilution)
    FL(i,:)= char(C) ;
end

for j=i:i+length(FringesEthyleneNoDilution)
    FL(j,:)= char(D) ;
end

for i=j:j+length(FringesEthylene67Dilution)
    FL(i,:)= char(E) ;
end

for j=i:i+length(FringesEthylene90Dilution)
    FL(j,:)= char(F) ;
end

for i=j:j+length(FringesEthaneNoDilution)
    FL(i,:)= char(G) ;
end

for j=i:i+length(FringesEthane50Dilution)
    FL(j,:)= char(H) ;
end

for i=j:j+length(FringesEthane85Dilution)
    FL(i,:)= char(I) ;
end
FringeLabel = char(FL);


%%
FLSD=zeros(length(FringesMethaneNoDilution),length('Methane No Dilution '));

for i=1:length(SeparationDistanceMethaneNoDilution)
    FLSD(i,:)= char(A) ;
end

for j=i:i+length(SeparationDistanceMethane50Dilution)
    FLSD(j,:)= char(B) ;
end

for i=j:j+length(SeparationDistanceMethane67Dilution)
    FLSD(i,:)= char(C) ;
end

for j=i:i+length(SeparationDistanceEthyleneNoDilution)
    FLSD(j,:)= char(D) ;
end

for i=j:j+length(SeparationDistanceEthylene67Dilution)
    FLSD(i,:)= char(E) ;
end

for j=i:i+length(SeparationDistanceEthylene90Dilution)
    FLSD(j,:)= char(F) ;
end

for i=j:j+length(SeparationDistanceEthaneNoDilution)
    FLSD(i,:)= char(G) ;
end

for j=i:i+length(SeparationDistanceEthane50Dilution)
    FLSD(j,:)= char(H) ;
end

for i=j:j+length(SeparationDistanceEthane85Dilution)
    FLSD(i,:)= char(I) ;
end

FLPSL = [A; B; C; D; E; F; G; H; I];
FringeLabelPSL=char(FLPSL);
FringeLabelSD = char(FLSD);


Fringes = [FringesMethaneNoDilution FringesMethane50Dilution FringesMethane67Dilution FringesEthyleneNoDilution FringesEthylene67Dilution FringesEthylene90Dilution FringesEthaneNoDilution FringesEthane50Dilution FringesEthane85Dilution];
Tortuosity = [TortuosityMethaneNoDilution TortuosityMethane50Dilution TortuosityMethane67Dilution TortuosityEthyleneNoDilution TortuosityEthylene67Dilution TortuosityEthylene90Dilution TortuosityEthaneNoDilution TortuosityEthane50Dilution TortuosityEthane85Dilution];
SeparationDistance = [SeparationDistanceMethaneNoDilution; SeparationDistanceMethane50Dilution; SeparationDistanceMethane67Dilution; SeparationDistanceEthyleneNoDilution; SeparationDistanceEthylene67Dilution; SeparationDistanceEthylene90Dilution; SeparationDistanceEthaneNoDilution; SeparationDistanceEthane50Dilution; SeparationDistanceEthane85Dilution];
PercentStackedFringes=[PercentStackedFringesMethaneNoDilution PercentStackedFringesMethane50Dilution PercentStackedFringesMethane67Dilution PercentStackedFringesEthyleneNoDilution PercentStackedFringesEthylene67Dilution PercentStackedFringesEthylene90Dilution PercentStackedFringesEthaneNoDilution PercentStackedFringesEthane50Dilution PercentStackedFringesEthane85Dilution];
%% What I want to plot
%Fringe Length
boxplot(Fringes,FringeLabel)
ylim([.45 1.5])

%Tortuosity
boxplot(Tortuosity,FringeLabel)
ylim([1 1.5])

boxplot(SeparationDistance,FringeLabelSD)
 
boxplot(PercentStackedFringes,FringeLabelPSL)


%% 4)Plot Methane fringe length and tortuosity Recirculating and Non recirculating case
x1=0:.01:2.5;   %x and such for the Fringe Length Plot
pdFringesMethaneNoDilution = fitdist(transpose(FringesMethaneNoDilution),'Kernel','BandWidth',.1);
yFringesMethaneNoDilution = pdf(pdFringesMethaneNoDilution,x1);
pdFringesMethane50Dilution = fitdist(transpose(FringesMethane50Dilution),'Kernel','BandWidth',.1);
yFringesMethane50Dilution = pdf(pdFringesMethane50Dilution,x1);
pdFringesMethane67Dilution = fitdist(transpose(FringesMethane67Dilution),'Kernel','BandWidth',.1);
yFringesMethane67Dilution = pdf(pdFringesMethane67Dilution,x1);
x2=1:.01:2;   %x and such for the Fringe Length Plot
pdTortuosityMethaneNoDilution = fitdist(transpose(TortuosityMethaneNoDilution),'Kernel','BandWidth',.1);
yTortuosityMethaneNoDilution = pdf(pdTortuosityMethaneNoDilution,x2);
pdTortuosityMethane50Dilution = fitdist(transpose(TortuosityMethane50Dilution),'Kernel','BandWidth',.1);
yTortuosityMethane50Dilution = pdf(pdTortuosityMethane50Dilution,x2);
pdTortuosityMethane67Dilution = fitdist(transpose(TortuosityMethane67Dilution),'Kernel','BandWidth',.1);
yTortuosityMethane67Dilution = pdf(pdTortuosityMethane67Dilution,x2);
x3=.3:.01:0.6;   %x and such for the Fringe Length Plot
pdSeparationDistanceMethaneNoDilution = fitdist(SeparationDistanceMethaneNoDilution,'Kernel','BandWidth',.1);
ySeparationDistanceMethaneNoDilution = pdf(pdSeparationDistanceMethaneNoDilution,x3);
pdSeparationDistanceMethane50Dilution = fitdist(SeparationDistanceMethane50Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceMethane50Dilution = pdf(pdSeparationDistanceMethane50Dilution,x3);
pdSeparationDistanceMethane67Dilution = fitdist(SeparationDistanceMethane67Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceMethane67Dilution = pdf(pdSeparationDistanceMethane67Dilution,x3);

a3=subplot(1,3,1); plot1Recirculation=plot(x1,yFringesMethaneNoDilution,'k-','LineWidth',[2]); 
hold on
a1=subplot(1,3,1); plot1Shell=plot(x1,yFringesMethane50Dilution,'k--','LineWidth',[2]);
a2=subplot(1,3,1); plot1Shell=plot(x1,yFringesMethane67Dilution,'b--','LineWidth',[2]);
xlabel('Fringe Length (nm)'); ylabel('Probability Density'); 
hold off

a4=subplot(1,3,2); plot2=plot(x2,yTortuosityMethaneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x2,yTortuosityMethane50Dilution,'k--','LineWidth',[2]);
a4=subplot(1,3,2); plot2=plot(x2,yTortuosityMethane67Dilution,'b--','LineWidth',[2]);

xlabel('Tortuosity'); ylabel('Probability Density');
legend('Methane No Dilution','Methane 50% Dilution','Methane 67% Dilution'); 
legend('boxoff');
hold off

a4=subplot(1,3,3); plot5 = plot(x3,ySeparationDistanceMethaneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x3,ySeparationDistanceMethane50Dilution,'k--','LineWidth',[2]);
a4=subplot(1,3,3); plot2=plot(x3,ySeparationDistanceMethane67Dilution,'b--','LineWidth',[2]);

xlabel('Separation Distance (nm)'); ylabel('Probability Density');
hold off

%% 4)Plot Ethane fringe length and tortuosity Recirculating and Non recirculating case
x1=0:.01:2.5;   %x and such for the Fringe Length Plot
pdFringesEthaneNoDilution = fitdist(transpose(FringesEthaneNoDilution),'Kernel','BandWidth',.1);
yFringesEthaneNoDilution = pdf(pdFringesEthaneNoDilution,x1);
pdFringesEthane50Dilution = fitdist(transpose(FringesEthane50Dilution),'Kernel','BandWidth',.1);
yFringesEthane50Dilution = pdf(pdFringesEthane50Dilution,x1);
pdFringesEthane85Dilution = fitdist(transpose(FringesEthane85Dilution),'Kernel','BandWidth',.1);
yFringesEthane85Dilution = pdf(pdFringesEthane85Dilution,x1);
x2=1:.01:2;   %x and such for the Fringe Length Plot
pdTortuosityEthaneNoDilution = fitdist(transpose(TortuosityEthaneNoDilution),'Kernel','BandWidth',.1);
yTortuosityEthaneNoDilution = pdf(pdTortuosityEthaneNoDilution,x2);
pdTortuosityEthane50Dilution = fitdist(transpose(TortuosityEthane50Dilution),'Kernel','BandWidth',.1);
yTortuosityEthane50Dilution = pdf(pdTortuosityEthane50Dilution,x2);
pdTortuosityEthane85Dilution = fitdist(transpose(TortuosityEthane85Dilution),'Kernel','BandWidth',.1);
yTortuosityEthane85Dilution = pdf(pdTortuosityEthane85Dilution,x2);
x3=.3:.01:0.6;   %x and such for the Fringe Length Plot
pdSeparationDistanceEthaneNoDilution = fitdist(SeparationDistanceEthaneNoDilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthaneNoDilution = pdf(pdSeparationDistanceEthaneNoDilution,x3);
pdSeparationDistanceEthane50Dilution = fitdist(SeparationDistanceEthane50Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthane50Dilution = pdf(pdSeparationDistanceEthane50Dilution,x3);
pdSeparationDistanceEthane85Dilution = fitdist(SeparationDistanceEthane85Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthane85Dilution = pdf(pdSeparationDistanceEthane85Dilution,x3);

a3=subplot(1,3,1); plot1Recirculation=plot(x1,yFringesEthaneNoDilution,'k-','LineWidth',[2]); 
hold on
a1=subplot(1,3,1); plot1Shell=plot(x1,yFringesEthane50Dilution,'k--','LineWidth',[2]);
a2=subplot(1,3,1); plot1Shell=plot(x1,yFringesEthane85Dilution,'b--','LineWidth',[2]);
xlabel('Fringe Length (nm)'); ylabel('Probability Density'); 
hold off

a4=subplot(1,3,2); plot2=plot(x2,yTortuosityMethaneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x2,yTortuosityEthane50Dilution,'k--','LineWidth',[2]);
a4=subplot(1,3,2); plot2=plot(x2,yTortuosityEthane85Dilution,'b--','LineWidth',[2]);

xlabel('Tortuosity'); ylabel('Probability Density');
legend('Ethane No Dilution','Ethane 50% Dilution','Ethane 67% Dilution'); 
legend('boxoff');
hold off

a4=subplot(1,3,3); plot5=plot(x3,ySeparationDistanceEthaneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x3,ySeparationDistanceEthane50Dilution,'k--','LineWidth',[2]);
a4=subplot(1,3,3); plot2=plot(x3,ySeparationDistanceEthane85Dilution,'b--','LineWidth',[2]);

xlabel('Separation Distance (nm)'); ylabel('Probability Density');
hold off
%% 4)Plot Ethylene fringe length and tortuosity Recirculating and Non recirculating case
x1=0:.01:2.5;   %x and such for the Fringe Length Plot
pdFringesEthyleneNoDilution = fitdist(transpose(FringesEthyleneNoDilution),'Kernel','BandWidth',.1);
yFringesEthyleneNoDilution = pdf(pdFringesEthyleneNoDilution,x1);
pdFringesEthylene67Dilution = fitdist(transpose(FringesEthylene67Dilution),'Kernel','BandWidth',.1);
yFringesEthylene67Dilution = pdf(pdFringesEthylene67Dilution,x1);
pdFringesEthylene90Dilution = fitdist(transpose(FringesEthylene90Dilution),'Kernel','BandWidth',.1);
yFringesEthylene90Dilution = pdf(pdFringesEthylene90Dilution,x1);
x2=1:.01:2;   %x and such for the Fringe Length Plot
pdTortuosityEthyleneNoDilution = fitdist(transpose(TortuosityEthyleneNoDilution),'Kernel','BandWidth',.1);
yTortuosityEthyleneNoDilution = pdf(pdTortuosityEthyleneNoDilution,x2);
pdTortuosityEthylene67Dilution = fitdist(transpose(TortuosityEthylene67Dilution),'Kernel','BandWidth',.1);
yTortuosityEthylene67Dilution = pdf(pdTortuosityEthylene67Dilution,x2);
pdTortuosityEthylene90Dilution = fitdist(transpose(TortuosityEthylene90Dilution),'Kernel','BandWidth',.1);
yTortuosityEthylene90Dilution = pdf(pdTortuosityEthylene90Dilution,x2);

x3=.3:.01:0.6;   %x and such for the Fringe Length Plot
pdSeparationDistanceEthyleneNoDilution = fitdist(SeparationDistanceEthyleneNoDilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthyleneNoDilution = pdf(pdSeparationDistanceEthyleneNoDilution,x3);
pdSeparationDistanceEthylene67Dilution = fitdist(SeparationDistanceEthylene67Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthylene67Dilution = pdf(pdSeparationDistanceEthylene67Dilution,x3);
pdSeparationDistanceEthylene90Dilution = fitdist(SeparationDistanceEthylene90Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthylene90Dilution = pdf(pdSeparationDistanceEthylene90Dilution,x3);

a3=subplot(1,3,1); plot1Recirculation=plot(x1,yFringesEthyleneNoDilution,'k-','LineWidth',[2]); 
hold on
a1=subplot(1,3,1); plot1Shell=plot(x1,yFringesEthylene67Dilution,'k--','LineWidth',[2]);
a2=subplot(1,3,1); plot1Shell=plot(x1,yFringesEthylene90Dilution,'b--','LineWidth',[2]);
xlabel('Fringe Length (nm)'); ylabel('Probability Density'); 
hold off

a4=subplot(1,3,2); plot2=plot(x2,yTortuosityEthyleneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x2,yTortuosityEthylene67Dilution,'k--','LineWidth',[2]);
a4=subplot(1,3,2); plot2=plot(x2,yTortuosityEthylene90Dilution,'b--','LineWidth',[2]);

xlabel('Tortuosity'); ylabel('Probability Density');
legend('Ethylene No Dilution','Ethylene 67% Dilution','Ethylene 90% Dilution'); 
legend('boxoff');
hold off

a4=subplot(1,3,3); plot5=plot(x3,ySeparationDistanceEthyleneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x3,ySeparationDistanceEthylene67Dilution,'k--','LineWidth',[2]);
a4=subplot(1,3,3); plot2=plot(x3,ySeparationDistanceEthylene90Dilution,'b--','LineWidth',[2]);

xlabel('Separation Distance (nm)'); ylabel('Probability Density');
hold off


%% 4)Plot Fringe Length, Tortuosity, and Separation Distance

x1 = 0: .01: 2.5;   %x and such for the Fringe Length Plot
pdFringesMethaneNoDilution = fitdist(transpose(FringesMethaneNoDilution),'Kernel','BandWidth',.1);
yFringesMethaneNoDilution = pdf(pdFringesMethaneNoDilution,x1);
pdFringesMethane50Dilution = fitdist(transpose(FringesMethane50Dilution),'Kernel','BandWidth',.1);
yFringesMethane50Dilution = pdf(pdFringesMethane50Dilution,x1);
pdFringesMethane67Dilution = fitdist(transpose(FringesMethane67Dilution),'Kernel','BandWidth',.1);
yFringesMethane67Dilution = pdf(pdFringesMethane67Dilution,x1);
x2 = 1: .01: 2;   %x and such for the Fringe Length Plot
pdTortuosityMethaneNoDilution = fitdist(transpose(TortuosityMethaneNoDilution),'Kernel','BandWidth',.1);
yTortuosityMethaneNoDilution = pdf(pdTortuosityMethaneNoDilution,x2);
pdTortuosityMethane50Dilution = fitdist(transpose(TortuosityMethane50Dilution),'Kernel','BandWidth',.1);
yTortuosityMethane50Dilution = pdf(pdTortuosityMethane50Dilution,x2);
pdTortuosityMethane67Dilution = fitdist(transpose(TortuosityMethane67Dilution),'Kernel','BandWidth',.1);
yTortuosityMethane67Dilution = pdf(pdTortuosityMethane67Dilution,x2);
x3=.3:.01:0.6;   %x and such for the Fringe Length Plot
pdSeparationDistanceMethaneNoDilution = fitdist(SeparationDistanceMethaneNoDilution,'Kernel','BandWidth',.1);
ySeparationDistanceMethaneNoDilution = pdf(pdSeparationDistanceMethaneNoDilution,x3);
pdSeparationDistanceMethane50Dilution = fitdist(SeparationDistanceMethane50Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceMethane50Dilution = pdf(pdSeparationDistanceMethane50Dilution,x3);
pdSeparationDistanceMethane67Dilution = fitdist(SeparationDistanceMethane67Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceMethane67Dilution = pdf(pdSeparationDistanceMethane67Dilution,x3);

a3=subplot(3,3,1); plot1Recirculation=plot(x1,yFringesMethaneNoDilution,'k-','LineWidth',[2]); 
hold on
a1=subplot(3,3,1); plot1Shell=plot(x1,yFringesMethane50Dilution,'k--','LineWidth',[2]);
a2=subplot(3,3,1); plot1Shell=plot(x1,yFringesMethane67Dilution,'b--','LineWidth',[2]);
xlabel('Fringe Length (nm)'); ylabel('Probability Density'); 
hold off

a4=subplot(3,3,2); plot2=plot(x2,yTortuosityMethaneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x2,yTortuosityMethane50Dilution,'k--','LineWidth',[2]);
a4=subplot(3,3,2); plot2=plot(x2,yTortuosityMethane67Dilution,'b--','LineWidth',[2]);

xlabel('Tortuosity'); ylabel('Probability Density');
legend('Methane No Dilution','Methane 50% Dilution','Methane 67% Dilution'); 
legend('boxoff');
hold off

a4=subplot(3,3,3); plot5=plot(x3,ySeparationDistanceMethaneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x3,ySeparationDistanceMethane50Dilution,'k--','LineWidth',[2]);
a4=subplot(3,3,3); plot2=plot(x3,ySeparationDistanceMethane67Dilution,'b--','LineWidth',[2]);

xlabel('Separation Distance (nm)'); ylabel('Probability Density');
hold off

x1=0:.01:2.5;   %x and such for the Fringe Length Plot
pdFringesEthaneNoDilution = fitdist(transpose(FringesEthaneNoDilution),'Kernel','BandWidth',.1);
yFringesEthaneNoDilution = pdf(pdFringesEthaneNoDilution,x1);
pdFringesEthane50Dilution = fitdist(transpose(FringesEthane50Dilution),'Kernel','BandWidth',.1);
yFringesEthane50Dilution = pdf(pdFringesEthane50Dilution,x1);
pdFringesEthane85Dilution = fitdist(transpose(FringesEthane85Dilution),'Kernel','BandWidth',.1);
yFringesEthane85Dilution = pdf(pdFringesEthane85Dilution,x1);
x2=1:.01:2;   %x and such for the Fringe Length Plot
pdTortuosityEthaneNoDilution = fitdist(transpose(TortuosityEthaneNoDilution),'Kernel','BandWidth',.1);
yTortuosityEthaneNoDilution = pdf(pdTortuosityEthaneNoDilution,x2);
pdTortuosityEthane50Dilution = fitdist(transpose(TortuosityEthane50Dilution),'Kernel','BandWidth',.1);
yTortuosityEthane50Dilution = pdf(pdTortuosityEthane50Dilution,x2);
pdTortuosityEthane85Dilution = fitdist(transpose(TortuosityEthane85Dilution),'Kernel','BandWidth',.1);
yTortuosityEthane85Dilution = pdf(pdTortuosityEthane85Dilution,x2);
x3=.3:.01:0.6;   %x and such for the Fringe Length Plot
pdSeparationDistanceEthaneNoDilution = fitdist(SeparationDistanceEthaneNoDilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthaneNoDilution = pdf(pdSeparationDistanceEthaneNoDilution,x3);
pdSeparationDistanceEthane50Dilution = fitdist(SeparationDistanceEthane50Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthane50Dilution = pdf(pdSeparationDistanceEthane50Dilution,x3);
pdSeparationDistanceEthane85Dilution = fitdist(SeparationDistanceEthane85Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthane85Dilution = pdf(pdSeparationDistanceEthane85Dilution,x3);

a3=subplot(3,3,4); plot1Recirculation=plot(x1,yFringesEthaneNoDilution,'k-','LineWidth',[2]); 
hold on
a1=subplot(3,3,4); plot1Shell=plot(x1,yFringesEthane50Dilution,'k--','LineWidth',[2]);
a2=subplot(3,3,4); plot1Shell=plot(x1,yFringesEthane85Dilution,'b--','LineWidth',[2]);
xlabel('Fringe Length (nm)'); ylabel('Probability Density'); 
hold off

a4=subplot(3,3,5); plot2=plot(x2,yTortuosityMethaneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x2,yTortuosityEthane50Dilution,'k--','LineWidth',[2]);
a4=subplot(3,3,5); plot2=plot(x2,yTortuosityEthane85Dilution,'b--','LineWidth',[2]);

xlabel('Tortuosity'); ylabel('Probability Density');
legend('Ethane No Dilution','Ethane 50% Dilution','Ethane 67% Dilution'); 
legend('boxoff');
hold off

a4=subplot(3,3,6); plot5=plot(x3,ySeparationDistanceEthaneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x3,ySeparationDistanceEthane50Dilution,'k--','LineWidth',[2]);
a4=subplot(3,3,6); plot2=plot(x3,ySeparationDistanceEthane85Dilution,'b--','LineWidth',[2]);

xlabel('Separation Distance (nm)'); ylabel('Probability Density');
hold off

x1=0: .01: 2.5;   %x and such for the Fringe Length Plot
pdFringesEthyleneNoDilution = fitdist(transpose(FringesEthyleneNoDilution),'Kernel','BandWidth',.1);
yFringesEthyleneNoDilution = pdf(pdFringesEthyleneNoDilution,x1);
pdFringesEthylene67Dilution = fitdist(transpose(FringesEthylene67Dilution),'Kernel','BandWidth',.1);
yFringesEthylene67Dilution = pdf(pdFringesEthylene67Dilution,x1);
pdFringesEthylene90Dilution = fitdist(transpose(FringesEthylene90Dilution),'Kernel','BandWidth',.1);
yFringesEthylene90Dilution = pdf(pdFringesEthylene90Dilution,x1);

x2=1:.01:2;   %x Tortuosity Plot
pdTortuosityEthyleneNoDilution = fitdist(transpose(TortuosityEthyleneNoDilution),'Kernel','BandWidth',.1);
yTortuosityEthyleneNoDilution = pdf(pdTortuosityEthyleneNoDilution,x2);
pdTortuosityEthylene67Dilution = fitdist(transpose(TortuosityEthylene67Dilution),'Kernel','BandWidth',.1);
yTortuosityEthylene67Dilution = pdf(pdTortuosityEthylene67Dilution,x2);
pdTortuosityEthylene90Dilution = fitdist(transpose(TortuosityEthylene90Dilution),'Kernel','BandWidth',.1);
yTortuosityEthylene90Dilution = pdf(pdTortuosityEthylene90Dilution,x2);

x3=.3:.01:0.6;   %x and such for the Fringe Length Plot
pdSeparationDistanceEthyleneNoDilution = fitdist(SeparationDistanceEthyleneNoDilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthyleneNoDilution = pdf(pdSeparationDistanceEthyleneNoDilution,x3);
pdSeparationDistanceEthylene67Dilution = fitdist(SeparationDistanceEthylene67Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthylene67Dilution = pdf(pdSeparationDistanceEthylene67Dilution,x3);
pdSeparationDistanceEthylene90Dilution = fitdist(SeparationDistanceEthylene90Dilution,'Kernel','BandWidth',.1);
ySeparationDistanceEthylene90Dilution = pdf(pdSeparationDistanceEthylene90Dilution,x3);

a3=subplot(3,3,7); plot1Recirculation=plot(x1,yFringesEthyleneNoDilution,'k-','LineWidth',[2]); 
hold on
a1=subplot(3,3,7); plot1Shell=plot(x1,yFringesEthylene67Dilution,'k--','LineWidth',[2]);
a2=subplot(3,3,7); plot1Shell=plot(x1,yFringesEthylene90Dilution,'b--','LineWidth',[2]);
xlabel('Fringe Length (nm)'); ylabel('Probability Density'); 
hold off

a4=subplot(3,3,8); plot2=plot(x2, yTortuosityEthyleneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x2,yTortuosityEthylene67Dilution,'k--','LineWidth',[2]);
a4=subplot(3,3,8); plot2=plot(x2, yTortuosityEthylene90Dilution,'b--','LineWidth',[2]);

xlabel('Tortuosity'); ylabel('Probability Density');
legend('Ethylene No Dilution','Ethylene 67% Dilution','Ethylene 90% Dilution'); 
legend('boxoff');
hold off

a4=subplot(3,3,9); plot5=plot(x3,ySeparationDistanceEthyleneNoDilution,'k','LineWidth',[2]);  %!!!!
hold on
plot2coreRecirc=plot(x3,ySeparationDistanceEthylene67Dilution,'k--','LineWidth',[2]);
a4=subplot(3,3,9); plot2=plot(x3,ySeparationDistanceEthylene90Dilution,'b--','LineWidth',[2]);

xlabel('Separation Distance (nm)'); ylabel('Probability Density');
hold off
