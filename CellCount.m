%CellCount.m
clearvars;
 
av1=1; %width of the narrow gaussian
av2=20; %width of the wide gaussian
thresh=1.8; %difference between the gaussians for detection 
cellsize = 150; %the minimum size of the cell 
cutoff=150; 
name='example.png';
 
%box 
x = 70;
y = 25;
xsize = 290;
ysize = 295; 
 
xbkg = 1;
ybkg = 1;  
textcolor = 'r';
 
A = 256*uint16(imread(name));
A_box = A;
A_box(y:y+ysize,x:x+1) = 65536;
A_box(y:y+1,x:x+xsize) = 65536;
A_box(y:y+ysize,x+xsize:x+xsize+1) = 65536;
A_box(y+ysize:y+ysize+1,x:x+xsize+1) = 65536;
A_box(ybkg:ybkg+20,xbkg:xbkg+1) = 65536;
A_box(ybkg:ybkg+1,xbkg:xbkg+20) = 65536;
A_box(ybkg:ybkg+20,xbkg+20:xbkg+21) = 65536;
A_box(ybkg+20:ybkg+21,xbkg:xbkg+21) = 65536;
figure(1); imshow(A_box)
 
H1=fspecial('gaussian',40,av1);
H2=fspecial('gaussian',40,av2);
a1=imfilter(A,H1,'replicate');
a2=imfilter(A,H2,'replicate');
a=a1-a2;
 
a_bw=medfilt2(a>(thresh*256));
a_bw2=bwmorph(a_bw,'open',inf);
a_bw3=imfill(a_bw2,'holes');
 
Amasked1=uint16(a_bw3).*A;
Am_box1 = Amasked1;
Am_box1(y:y+ysize,x:x+1) = 65536;
Am_box1(y:y+1,x:x+xsize) = 65536;
Am_box1(y:y+ysize,x+xsize:x+xsize+1) = 65536;
Am_box1(y+ysize:y+ysize+1,x:x+xsize+1) = 65536;
Am_box1(ybkg:ybkg+20,xbkg:xbkg+1) = 65536;
Am_box1(ybkg:ybkg+1,xbkg:xbkg+20) = 65536;
Am_box1(ybkg:ybkg+20,xbkg+20:xbkg+21) = 65536;
Am_box1(ybkg+20:ybkg+21,xbkg:xbkg+21) = 65536;
 
figure(2), imshow(Am_box1)
 
bkg=mean2(A(ybkg:20+ybkg,xbkg:20+xbkg))/256;
meanpix1=sum(sum(Amasked1(y:y+ysize,x:x+xsize)))/sum(sum(a_bw3(y:y+ysize,x:x+xsize)))/256;
meanpix1_bkgcorr=meanpix1-bkg
 
number_bins=20;
L=bwlabel(a_bw3);
objects=max(max(L));
a_bw4=false(size(a_bw3));
for m=1:objects
    av(m)=sum(sum(uint16(L==m).*A))./sum(sum(uint16(L==m)))/256;
    if av(m)<cutoff
        a_bw4=a_bw4+(L==m);
    end
end
 
bins=(floor(0.9*min(av))):(ceil(1.1*max(av)));
figure(3)
hist(av,bins), axis tight
avg=mean(av);
stdev=std(av);
av_2stdev=avg+2*stdev;
 
Amasked2=uint16(a_bw4).*A;
Am_box2 = Amasked2;
Am_box2(y:y+ysize,x:x+1) = 65536;
Am_box2(y:y+1,x:x+xsize) = 65536;
Am_box2(y:y+ysize,x+xsize:x+xsize+1) = 65536;
Am_box2(y+ysize:y+ysize+1,x:x+xsize+1) = 65536;
Am_box2(ybkg:ybkg+20,xbkg:xbkg+1) = 65536;
Am_box2(ybkg:ybkg+1,xbkg:xbkg+20) = 65536;
Am_box2(ybkg:ybkg+20,xbkg+20:xbkg+21) = 65536;
Am_box2(ybkg+20:ybkg+21,xbkg:xbkg+21) = 65536;
figure(4), imshow(Am_box2)
meanpix2=sum(sum(Amasked2(y:y+ysize,x:x+xsize)))/sum(sum(a_bw4(y:y+ysize,x:x+xsize)))/256;
meanpix2_bkgcorr=meanpix2-bkg
 
a_interest = a_bw4(y:y+ysize,x:x+xsize);
 
[L cells] = bwlabel(a_interest); 
cells
STATS = regionprops(L,'basic');
smallcells = 0; 
scidenity = 0;
for m=1:cells
    centroid(m,1) = STATS(m).Centroid(1,1);
    centroid(m,2) = STATS(m).Centroid(1,2);
    area(m) = STATS(m).Area;
    if(area(m) < cellsize)
        smallcells = smallcells + 1;
        scidenity(smallcells) = m;
    end
end
cells_threshold = cells - smallcells
figure(5)
imshow(uint16(a_interest).*A(y:y+ysize,x:x+xsize))
for m=1:cells
    text(centroid(m,1),centroid(m,2),num2str(m),'Color',textcolor)
end
text(0,-10,['Ex. Cells: ' num2str(scidenity)],'Color','k')
