% %%% Load the image that needs to be processed
% clear all
pic = 'TEM1054.tif';
A = imread(pic);
resolution = .0208;    %nm/pixel   .0208   .0231   0.0261
%% 1)Define a region of interest (ROI) to analyze the particle. Draw freehand, 
%double click or exit out of image after drawing to continue script
imshow(A)
roi = createMask(imfreehand);
position = wait(imfreehand);

for i = 1: size(A,1)
    for j = 1: length(A)
        if roi(i,j)==0
            A(i,j)=0;
        end      
    end
end

%% 2) Adjust the pictures contrast with a histogram equalization. 
% B = histeq(A);             %Increase Contrast
B = adapthisteq(A);        %Increase Contrast

%% 3)Remove Noise with a low pass Gaussian filter. Botero uses a 7 pixel
%%filter size and 3 pixel deviation
Gaus = fspecial('gaussian',[11 11],5);
C = imfilter(B,Gaus);
% P = imgaussfilt(B ,  5,'FilterSize',[11 11]);

%% 4) Bottom Hat transformation. Botero uses a 2x2 disk element
se = strel('disk',4);                   
D = imbothat(C,se);

%% 5) Sets a threshold and binarizes the image using Otsu's method
level = graythresh(D); 
E = imbinarize(D,level);    

%% 6) skeletonization of the image
F = bwmorph(E,'skel',Inf);
F = bwmorph(F,'fill');
F = bwmorph (F,'skel');

%% 7) Further morphological operations of cleaning isolated pixels and breaking fringes in the h form
G = bwmorph (F,'clean');
G = bwmorph ( F,'hbreak',8);


%% 8) breaks apart fringes with >3 connections
I=G;
% I( (conv2(I,[1 1 1;1 0 1;1 1 1],'same') .* I) > 2 ) = 0;

pts = bwlookup(G,makelut(@(x) sum(x(:))>=4 & x(5)==1,3)); %branch points?
I = G&~pts;
I = bwmorph(I,'clean');
I = I+pts/2;


%% do it again
i = 1;
for i=1:3
I(I>0)=1;
pts2 = bwlookup(I,makelut(@(x) sum(x(:))>=4 & x(5)==1,3)); %branch points?
I = I&~pts2;
I = bwmorph(I,'clean');
I = I+pts2/2;
end

I(I > 0) = 1;
% imshowpair(G,A,'montage')

   %% Measure the length of fringes in path walked start to end
% measurements = regionprops(G,'area');    % counts the total number of pixels
[l, k] = find(bwmorph(I,'endpoints'));
fringemeasurement=zeros(1,length(l));
H = I;
for a = 1: length(l)
    i = l(a);
    j = k(a);
    fringemeasurement(a) = 1;
    while H(i, j) == 1
        if bwarea(H(i + 1, j)) == 1            %% if the location TTR is true, then add 1
            fringemeasurement(a) = fringemeasurement(a) + 1;
            H(i,j)=0;
            i=i+1;                                  % Follow that path
        elseif bwarea(H(i,j+1))==1          %%% if the location Above is true, then add 1
            fringemeasurement(a)= fringemeasurement(a) +1;
            H(i,j)=0;
            j=j+1;
        elseif bwarea(H(i-1,j))==1          %% if the location TTL is true, then add 1
            fringemeasurement(a)= fringemeasurement(a) +1;
            H(i,j)=0;
            i=i-1;
        elseif bwarea(H(i,j-1))==1          %%% if the location Below is true, then add 1
            fringemeasurement(a)= fringemeasurement(a) +1;
            H(i,j)=0;
            j=j-1;
        elseif bwarea(H(i+1,j+1))==1          %%% if the location TTR and Above is true, then add 1
            fringemeasurement(a)= fringemeasurement(a) + sqrt(2);
            H(i,j)=0;
            i=i+1;
            j=j+1;
        elseif bwarea(H(i+1,j-1))==1          %%% if the location TTR and Below is true, then add 1
            fringemeasurement(a)= fringemeasurement(a) +sqrt(2);
            H(i,j)=0;
            i=i+1;
            j=j-1;
         elseif bwarea(H(i-1,j-1))==1          %%% if the location TTL and Below is true, then add 1
            fringemeasurement(a)= fringemeasurement(a) +sqrt(2);
             H(i,j)=0;
              i=i-1;
              j=j-1;
         elseif bwarea(H(i-1,j+1))==1          %%% if the location TTL and above is true, then add 1
            fringemeasurement(a)= fringemeasurement(a) +sqrt(2);
            H(i,j)=0;
            j=j+1;
            i=i-1;
        else
            H(i,j)=0;
        end
    end
end

%% Measure the distance between fringe endpoints
[l,k] = find(bwmorph(I,'endpoints'));
euclength=zeros(1,length(l));
H=I;
for a = 1:length(l)
    i = l(a);
    j = k(a);
    while H(i,j)== 1
        if bwarea(H(i+1,j)) == 1            %% if the location TTR is true, then add 1
            H(i,j)=0;
            i=i+1;                                  % Follow that path
        elseif bwarea(H(i,j+1))==1          %%% if the location Above is true, then add 1
            H(i,j)=0;
            j=j+1;
        elseif bwarea(H(i-1,j))==1          %% if the location TTL is true, then add 1
            H(i,j)=0;
            i=i-1;
        elseif bwarea(H(i,j-1))==1          %%% if the location Below is true, then add 1
            H(i,j)=0;
            j=j-1;
        elseif bwarea(H(i+1,j+1))==1          %%% if the location TTR and Above is true, then add 1
            H(i,j)=0;
            i=i+1;
            j=j+1;
        elseif bwarea(H(i+1,j-1))==1          %%% if the location TTR and Below is true, then add 1
            H(i,j)=0;
            i=i+1;
            j=j-1;
         elseif bwarea(H(i-1,j-1))==1          %%% if the location TTL and Below is true, then add 1
             H(i,j)=0;
              i=i-1;
              j=j-1;
         elseif bwarea(H(i-1,j+1))==1          %%% if the location TTL and above is true, then add 1
            H(i,j)=0;
            j=j+1;
            i=i-1;
        else
            H(i,j)=0;
        end
        euclength(a)=sqrt((i-l(a))^2+(j-k(a))^2);%i and j are then the end points. So sqrt((i-l(a))^2+(j-k(a))^2)
    end
end

%% Remove the uneeded ones, and measure the tortuosity
indices = fringemeasurement < .483/resolution;    % < 17.9 pixels (.483 nm) deemed unrealistic becasue size of two rings
fringemeasurement(indices) = [];
euclength(indices) = [];
tortuosity= fringemeasurement./euclength;
fringelength=fringemeasurement*resolution;

%% Calculate the angle of each point
[l,k] = find(bwmorph(I,'endpoints'));
l(indices)=[];          %removing the ones that are shorter than the required distance
k(indices)=[];          %removing the ones that are shorter than the required distance
H=double(I);
angle=zeros(2048:2048);
for a = 1:length(l)
    i = l(a);
    j = k(a);
    while H(i,j)== 1
        if bwarea(H(i+1,j)) == 1            %% if the location above is true, then angle is 90 degrees
            angle(i,j)=270;
            H(i,j)=0;
            i=i+1;                                  % Follow that path
        elseif bwarea(H(i,j+1))==1          %%% if the location Above is true, then add 1
            angle(i,j)=360;
            H(i,j)=0;
            j=j+1;
        elseif bwarea(H(i-1,j))==1          %% if the location TTL is true, then add 1
            angle(i,j)=90;
            H(i,j)=0;
            i=i-1;
        elseif bwarea(H(i,j-1))==1          %%% if the location Below is true, then add 1
            angle(i,j)=180;
            H(i,j)=0;
            j=j-1;
        elseif bwarea(H(i+1,j+1))==1          %%% if the location TTR and Above is true, then add 1
            angle(i,j)=315;
            H(i,j)=0;
            i=i+1;
            j=j+1;
        elseif bwarea(H(i+1,j-1))==1          %%% if the location TTR and Below is true, then add 1
            angle(i,j)=225;
            H(i,j)=0;
            i=i+1;
            j=j-1;
         elseif bwarea(H(i-1,j-1))==1          %%% if the location TTL and Below is true, then add 1
             angle(i,j)=135;
             H(i,j)=0;
              i=i-1;
              j=j-1;
         elseif bwarea(H(i-1,j+1))==1          %%% if the location TTL and above is true, then add 1
            angle(i,j)=45;
            H(i,j)=0;
            j=j+1;
            i=i-1;
        else
            H(i,j)=0;
            angle(i,j)=sum(sum(angle(i-1:i+1,j-1:j+1)));
        end
    end
end

%% Label each fringe with a unique number
[l,k] = find(bwmorph(I,'endpoints'));
l(indices)=[];          %removing the ones that are shorter than the required distance
k(indices)=[];          %removing the ones that are shorter than the required distance
Unique = double(I);
for a = 1:length(l)
    i = l(a);
    j = k(a);
    while Unique(i,j)== 1
        Unique(i,j)= a+1;
        if sum(sum(roi((i-1:i+1),(j-1:j+1))))~=9
            Unique(i,j)=0;
        elseif Unique(i+1,j) == 1            %% if the location below is true, then add 1
            i=i+1;
        elseif Unique(i,j+1)== 1          %%% if the location Above is true, then add 1
            j=j+1;
        elseif Unique(i-1,j)==1          %% if the location TTL is true, then add 1
            i=i-1;
        elseif Unique(i,j-1)==1          %%% if the location Below is true, then add 1
            j=j-1;
        elseif Unique(i+1,j+1)==1          %%% if the location TTR and Above is true, then add 1
            i=i+1;
            j=j+1;
        elseif Unique(i+1,j-1)==1          %%% if the location TTR and Below is true, then add 1
            i=i+1;
            j=j-1;
         elseif Unique(i-1,j-1)==1          %%% if the location TTL and Below is true, then add 1
              i=i-1;
              j=j-1;
         elseif Unique(i-1,j+1)==1          %%% if the location TTL and above is true, then add 1
            j=j+1;
            i=i-1;
        end
    end
end
%% Remove all points for short fringes
for i=1:2048
    for j=1:2048
        if Unique(i,j)==1
           Unique(i,j)=0;
        end
    end
end

%% 5 nearest neighbors average angle
[m,n]=find(angle);
angleave=zeros(2048,2048);
x=zeros(2048,2048);
y=zeros(2048,2048);

for i=1:length(m)
    x(m(i),n(i))=cosd(angle(m(i),n(i)));
    y(m(i),n(i))=sind(angle(m(i),n(i)));
end
for i=1:length(m)
    a = sum(sum(x(m(i)-2:m(i)+2,n(i)-2:n(i)+2)));
    b = sum(sum(y(m(i)-2:m(i)+2,n(i)-2:n(i)+2)));
    angleave(m(i),n(i))= atan2d(b,a);
end

for i=1:2048
    for j=1:2048
        if Unique(i,j)==0
            angleave(i,j)=0;
        end
    end
end
[m,n]=find(bwmorph(Unique,'endpoints'));

for i=1:length(m)
    a = sum(sum(x(m(i)-1:m(i)+1,n(i)-1:n(i)+1)));
    b = sum(sum(y(m(i)-1:m(i)+1,n(i)-1:n(i)+1)));
    angleave(m(i),n(i))= atan2d(b,a);
end
%% then need to say if there are any other fringes with the same vertical position and log if they are
% [l,n] = find(Unique);   %points of G that are on the fringes greater than a certain distance
verticali1=zeros(length(l),length(l));
verticali2=zeros(length(l),length(l));
verticalj=zeros(length(l),length(l));
stackedfringesi=[0];
stackedfringesj=[0];%zeros(length(l),length(l));
StackedDist =  0.8;   %(in nanometers)
UniqueFringes=length(unique(Unique));
angldev=10;
separation=0;
for iter = 1 : UniqueFringes
     [i,j] = find(Unique==iter);       %      i = l(iter);      %points on Unique matrix
       for q=1:length(i)     %      j = n(iter);       %points on Unique matrix
                if angleave(i(q),j(q))>0
                theta=(angleave(i(q),j(q)))-90;
                else
                theta=(angleave(i(q),j(q)))+90;
                end
            mi = i(q) + round(StackedDist/resolution*sind(theta));             
            mj = j(q) - round(StackedDist/resolution*cosd(theta));
            for z=1:round(2*StackedDist/resolution)  %while i < 2048 if the fringe numbers are equal do nothing, else log fringe
            if resolution*((i(q)-(mi-round(z*sind(theta))))^2+(j(q)-(mj+ round(z*cosd(theta))))^2)^.5  > StackedDist+.1
                meta
            elseif  Unique(i(q),j(q)) == Unique (mi-round(z*sind(theta)),mj+ round(z*cosd(theta)))      %if on same fringe, NEXT
            elseif Unique(mi-round(z*sind(theta)),mj+ round(z*cosd(theta))) > 1  & angleave(i(q),j(q))-angldev <= angleave(mi-round(z*sind(theta)),mj+ round(z*cosd(theta))) & angleave(i(q),j(q))+angldev>=angleave(mi-round(z*sind(theta)),mj+ round(z*cosd(theta))) 
             stackedfringesi=[stackedfringesi; Unique(i(q),j(q))];
             stackedfringesj=[stackedfringesj;Unique(mi-round(z*sind(theta)),mj+ round(z*cosd(theta)))]; %log the two fringe numbers  (iter,j)
             separation=[separation;((i(q)-(mi-round(z*sind(theta))))^2+(j(q)-(mj+ round(z*cosd(theta))))^2)^.5];
            else
            end
           end
     end
end
separation(separation<.00001)=[];
stackedfringesi(stackedfringesi<1)=[];
stackedfringesj(stackedfringesj<1)=[];
separation=separation*resolution;
stackedfringe=[stackedfringesi stackedfringesj];
stackedfringe = unique(stackedfringe);
% length(stackedfringe)/max(Unique(:))
% [mean(fringelength) mean(tortuosity) length(stackedfringesfin)/max(Unique(:)) max(Unique(:))]

%% Need to change script in order to account for the distance between pixels. Record the distance, then if they are 'stacked', add to the separation array. If not, delete the array
%Is there a way to index which specific pixels are stacked? The current
%process goes 1) find alll the stacked pixels, and label them by saying
%that fringe number 2) see how many repitions there are for that given
%pair, say fringe1 and fringe2. So this strategy is no bueno. I can log the
%distance as script finds the ones close to each other, then delete the
%inputs that are removed?
stackedfringes=[0 0];
stckdpxl=0.20/resolution;     %The number is nm
separationtot=0;
separationmean=0;
for i=1: length (stackedfringe)
    repititions=0;
    points = find(stackedfringesi == stackedfringe(i));
    sep=0;
    for j = 1:length (points)
        if stackedfringesj(points(1)) == stackedfringesj(points(j))
            if stackedfringe(i) == stackedfringesi(points(j)) & stackedfringesj(points(1)) == stackedfringesj(points(j))
            repititions = repititions + 1;   
            sep=[sep;separation(points(j))];
            end
            if repititions > stckdpxl
            stackedfringes = [stackedfringes; stackedfringesi(points(j)) stackedfringesj(points(j))];
            separationtot=[separationtot;sep];
            separationmean=[mean(sep);separationmean];
            repititions=0;
            end
        elseif stackedfringesj(points(1)) ~= stackedfringesj(points(j))
            if stackedfringe(i) == stackedfringesi(points(j)) & stackedfringesj(points(1)) == stackedfringesj(points(j))
            repititions = repititions + 1;
            sep=[sep;separation(points(j))];
            end
            if repititions > stckdpxl
            stackedfringes = [stackedfringes; stackedfringesi(points(j)) stackedfringesj(points(j))];
            separationtot=[separationtot;sep];
            separationmean=[mean(sep);separationmean];
            end
        end
    end
end
separationtot(separationtot < .0001)=[];
separationmean(separationmean < .00001)=[];
separationtot(separationtot > 0.6) = []; %0.6)=[];
separationmean(separationmean > 0.6) = [];     %2*min(separationtot)
stackedfringesfin= [stackedfringes(:, 1) ; stackedfringes(:, 2)];
stackedfringesfin(stackedfringesfin < 1) = [];
stackedfringesfin = unique(stackedfringesfin);
% length(stackedfringesfin)/length(unique(Unique))
%% Plot Stuff Plot fringe length and tortuosity for each individual trial
% plot1=histogram(fringelength,50,'binlimits',[0.4,6],'Normalization','probability');
% plot2=histogram(tortuosity,50,'binlimits',[1,5],'Normalization','probability');
% 
% subplot(1,2,1); plot1=histogram(fringelength,50,'binlimits',[0.4,4],'Normalization','probability');
% subplot(1,2,2); plot2=histogram(tortuosity,50,'binlimits',[1,2],'Normalization','probability');
% subplot(2,2,3); plot1=histogram(dishor,20,'binlimits',[0.3,0.9],'Normalization','probability');
% subplot(2,2,4); plot2=histogram(disvert,20,'binlimits',[0.3,0.9],'Normalization','probability');

% HighTor= sum(tortuosity > 1.5)/sum(tortuosity>0);  %HighTor is the percent of fringes with tortuosity>1.5
c = unique(stackedfringe);
%% Plot it
% yyaxis left
% plot(distance.*resolution,results(:,3),distance.*resolution,results(:,2))
% xlabel('Distance From Center (nm)')
% ylabel('Length (nm)                        Ratio')
% yyaxis right
% plot(distance.*resolution,results(:,4))
% ylabel('Percent Stacked')
% legend('Fringe Length','Tortuosity','Percent of Stacked Layers','Location','southeast')
% legend('boxoff')

%% After Loading the excel file in that contains Standard Deviation of: Tortuosity, Fringe Length, Percent stacked as 'results.mat'
% errtor = (results(:,1));
% errFL=(results(:,2));
% errstacked=results(:,3);
%%
ColorMap= zeros(2048: 2048);
stackedfringesfin(stackedfringesfin < 1) = [];
for i = 1: length(Unique)
    for j = 1: length(Unique)
        if Unique(i,j) ~= 0
            ColorMap(i, j) = 10;
        end
    end
end

for i=1:length(stackedfringesfin)
    [l,m] = find(Unique == stackedfringesfin(i));
    for k=1:length(l)
           ColorMap(l(k),m(k))=500;
    end
end
image(ColorMap)
% imagesc(ColorMap)
%% Plot it
% yyaxis left
% errorbar(distance.*resolution,results(:,2),errFL)
% hold on
% errorbar(distance.*resolution,results(:,1),errtor)
% xlabel('Distance From Center (nm)')
% ylabel('Length (nm)                        Ratio')
% yyaxis right
% yyaxis right
% errorbar(distance.*resolution,results(:,3),errstacked)
% ylabel('Percent Stacked')
% legend('Fringe Length','Tortuosity','Percent of Stacked Layers','Location','southeast')
% legend('boxoff')
% hold off
Results = ["FL" "T" "PSL" "Separation" "Fringes"; mean(fringelength) mean(tortuosity) length(stackedfringesfin)/length(unique(Unique)) mean(separationmean) length(unique(stackedfringe))]

%% Thicken the Unique into Unique2
Unique2 = Unique;


StackedYes= zeros(2048:2048);
for i=1:length(stackedfringesfin)
    [m,n]=find(Unique2==stackedfringesfin(i));
    for j=1:numel(m)
    StackedYes(m(j),n(j)) = 1;
    end
end

[m,n] = find (StackedYes);
for i=1:numel(m)
    StackedYes(m(i),n(i)) = 255;
    StackedYes(m(i)+1,n(i)) = 255;
    StackedYes(m(i),n(i)+1) = 255;
    StackedYes(m(i)-1,n(i)) = 255;
    StackedYes(m(i),n(i)-1) = 255;
end
% 
%% StackedNo = Unique&~StackedYes;

[m,n] = find (Unique2);
for i=1:numel(m)
    Unique2(m(i),n(i)) = 255;
    Unique2(m(i)+1,n(i)) = 255;
    Unique2(m(i),n(i)+1) = 255;
    Unique2(m(i)-1,n(i)) = 255;
    Unique2(m(i),n(i)-1) = 255;
end

A=imread(pic);
imshow(Unique2, 'InitialMag', 'fit')

% % Make a truecolor all-green image.
green = cat(3, 255*ones(size(StackedYes)), 255*ones(size(StackedYes)), zeros(size(StackedYes)));
hold on
h = imshow(green);

hold off
% 
set(h, 'AlphaData', StackedYes)

%%
% fnum=101;
% [m,n] = find(Unique==fnum);
% thisfringe=Unique(min(m)-2:max(m)+2,min(n)-2:max(n)+2);
% thisfringe=Unique(525:555,565:590);
% imwrite(thisfringe,'thisfringe.tif')
% % % imshow(thisfringe)
% [fringelength(fnum-1) tortuosity(fnum-1)]
%%
% [m,n] = find (Unique2);
% for i=1:numel(m)
%     Unique2(m(i),n(i)) = 255;
%     Unique2(m(i)+1,n(i)) = 255;
%     Unique2(m(i),n(i)+1) = 255;
%     Unique2(m(i)-1,n(i)) = 255;
%     Unique2(m(i),n(i)-1) = 255;
% end


%% boxplot historam of the average lattice spacing
% boxplot(separationmean)
% h = histogram(SeparationDistanceMethaneNoDilution, 19)
% set(h , 'FaceColor', 'k')