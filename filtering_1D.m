clear
img = imread("evert.pgm");
dim = size(img);
lowFilter = zeros(dim(1),dim(2));
borne = 10;
lowFilter(borne+1:dim(1)-borne,:) = 1;
lowFilter(:,borne+1:dim(2)-borne) = 1;
filter = lowFilter;

for i=1:size(img,1)
  fft1Dimg(i,:)=fft(img(i,:));
endfor
for j=1:size(img,2)
   fft2Dimg(:,j)=fft(fft1Dimg(:,j));
endfor

fft2DFiltered = filter.*fft2Dimg;

for i=1:size(img,1)
   filterImg1D(i,:)=ifft(fft2DFiltered(i,:));
endfor
for j=1:size(img,2)
   filterImg(:,j)=ifft(filterImg1D(:,j));
endfor
imagesc(real(filterImg));
axis equal off
colormap(gray);
