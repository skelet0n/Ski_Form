A = imread(filename);
%29 x 8 = 232  = #boxes

%Only compute with subset of image for now
% N =900;
% A = A(1:N,30:1000,:);

%Is transform invariant under rottations?
% A = imrotate(A,0,'crop');

Abw = rgb2gray(A);


%get size of image
[mtmp,ntmp] = size(Abw);

%find where first edges are and crop accordingly
thresh = .8;
Abw_e = edge(Abw,'canny',thresh);
    
[i,j] = find(Abw_e);

imin = min(i);jmin = min(j);

imax = max(i);jmax = max(j);

bw = Abw_e;
% Abw = imcrop(Abw,[imin jmin imax jmax]);
% bw = imcrop(Abw_e,[imin jmin imax jmax]);

%get new size of image
[mtmp,ntmp] = size(Abw);


%params for 'ski.jpg'
% alpha1 = 1/6;alpha2 =  3/4;
% alpha3 = 1/6;alpha4 = 4/5;

%params for 'all_xs.jp' and scans/.*
alpha1 = 1/4; alpha2 = 3.95/5;
alpha3 = 2/8; alpha4 = 3.5/5;

% alpha1 = 0.01;alpha2 = 1;
% alpha3 = 0.01;alpha4 = 1;

%params for 'stef_scans/.*'
% alpha1 = 1/4 + 1/16; alpha2 = 3.95/5;
% alpha3 = 2/8; alpha4 = 3.5/5;
% 
%params for cropping th etop left co rner first
% alpha1 = 4/12;alpha2 = 12/12;
% alpha3 = 1/3; alpha4 = 3/4;
% 
m1 = ceil(alpha1*mtmp);m2 = ceil(alpha2*mtmp);
n1 = ceil(alpha3*ntmp);n2 = ceil(alpha4*ntmp);
% 
%do some ad hoc cropping
Abw = Abw(m1:m2,n1:n2);

Abw_e = Abw_e(m1:m2,n1:n2);

%Do edge detection
% thresh = 1/2^3;
% bw = edge(Abw,'canny');

bw = bw(m1:m2,n1:n2);






% bw = bw(m1:m2,n1:n2);

%imshow(bw);


%acceptable angles for lines:
t = 1/2^3;tw = 10;
theta = [-90:t:-90+tw -tw:t:tw 90-tw:t:90-t];

%apply hough transform
[H,T,R] = hough(bw,'theta',theta);

%grab most used lines
numpeaks = 100;
P = houghpeaks(H,numpeaks);

%grab lines found by houghpeaks
minlength = 50;fillgap = 200;
lines = houghlines(bw,T,R,P,'minlength',minlength,'fillgap',fillgap);


% lines are found, seperate into horizontal and vertical lines


hind = false(length(lines),1);
vind = ~hind;
lines_rho = zeros(length(lines),1);

for k = 1:length(lines)
    %for sorting purposes
    lines_rho(k) = abs(lines(k).rho);
end

%sort lines by translation
[lines_rho_sorted,sortind] =  sort(lines_rho);

%reorder lines
lines = lines(sortind);

for k=1:length(lines)
    %line info is stored in lines,
    line_angle = lines(k).theta;
    %angles around 0 are horizontal lines, get indices of horizontal and
    %vertical lines
    hthresh =  10;
    bool = abs(line_angle)<hthresh;
    hind(k) = ~bool;
    vind(k) = bool;
end


hlines = lines(hind);
vlines = lines(vind);


%need to detect if lines are really close together... (will be of the same
%class vertical/horizontal)

% take a point from line 1, and project it onto line 2. If this projection
% is relatively small then the lines are 'close'. Tolerance should be <7
% pixels ish

%tolerance for checking if lines are too close together.
linetol = 20;line_ind = [];
for k = length(hlines)-2:length(hlines)-1
   %Check if horizontal lines are 'close'
   
   %compute normal of line k
   x0 = hlines(k).point1;
   xend = hlines(k).point2;
   
   len = norm(xend-x0);
   uhat = (xend- x0)/len;
   nhat = fliplr(uhat);nhat(1) = -nhat(1);
   
   %take point on k +1 line and treat x0 as the origin 
   xp = hlines(k+1).point1 - x0;
   
   
   %decompose xp into new basis
   %xp = alpha uhat + beta nhat
   
   beta = dot(xp,nhat);
   
   if abs(beta)<linetol
       %UH OH, lines are too close
       line_ind = [line_ind k];
       disp('found lines that are too close');
   end
   
   
   
   
end

%delete lines that are too close to other lines
hlines(line_ind) = [];



%tolerance for checking if lines are too close together.
linetol = 20;line_ind = [];
for k = 1:length(vlines)-1
   %Check if horizontal lines are 'close'
   
   %compute normal of line k
   x0 = vlines(k).point1;
   xend = vlines(k).point2;
   
   len = norm(xend-x0);
   uhat = (xend- x0)/len;
   nhat = fliplr(uhat);nhat(1) = -nhat(1);
   
   %take point on k +1 line and treat x0 as the origin 
   xp = vlines(k+1).point1 - x0;
   
   
   %decompose xp into new basis
   %xp = alpha uhat + beta nhat
   
   beta = dot(xp,nhat);
   if abs(beta)<linetol
       %UH OH, lines are too close
       line_ind = [line_ind k];
       disp('found lines that are too close');
   end
   
   
   
   
end

vlines(line_ind) = [];





%format of form requires you to throw out the first three horiztonal lines,
%and the first three vertical lines (top/side of form)
% hlines = hlines(1:end-2);
% vlines = vlines(2:end);


%for each horizontal line, find point of intersection with all vertical
%lines
m = length(hlines);
n = length(vlines);

inters_x = zeros(m,n);
inters_y = zeros(m,n);

if doplot
    figure(1);imshow(Abw_e); hold on;
end

for i = 1:m
    
    for j = 1:n
        %find point of intersection of:
        %horizontal lines i and vertical line j
        hline = hlines(i);vline = vlines(j);
        
        
        h1 = hline.point1;h2 = hline.point2;
        v1 = vline.point1;v2 = vline.point2;
        
        % want a*h1 + (1-a)h2 = b v1 + (1-b) v2, %for some a,b
        % h2 + a(h1 - h2) = v2 + b(v1-v2);
        
        %build system
        S = [(h1-h2)' -(v1-v2)'];
        rhs = v2'-h2';
        
        soln = S\rhs;
        
        a = soln(1);
        xstar = a*h1 + (1-a)*h2;
        
        %store point of intersection
        inters_x(i,j) = xstar(1);
        inters_y(i,j) = xstar(2);
    end
    
end



%show image overlayed with lines

if doplot
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    
    % plot beginnings and ends of lines
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    
end
% plot(inters_x(:),inters_y(:),'bp','linewidth',2)

hold off;
end

%% Get subimages
%(inters_x,inters_y) give the points of intersections of the horizontal and
%vertical lines

%subimage given given by:
% [ inters_x(i,j): inters_x(i+1,j+1)] \times [ inters_y(i,j) :
% inters__y(i+1,j+1)]


%NOTE: ORIGIN IS BOTTOM LEFT
subimages = cell((m-1)*(n-1),1);

%SUBTRACT OFF DETECTED LINES FROM IMAGE
%this way subimages won't have black lines going through them, or at least
%will be less noticable.


% define width of line,
w = 8;

A_nolines = Abw;%zeros(size(Abw));% Abw;

for k = 1:length(lines)
    %get endpoints of line (are row vectors)
    x0 = lines(k).point1;xend = lines(k).point2;
    
    len = norm(xend - x0);
    
    %find normal vector to line
    vhat = (xend-x0)/len;
    nhat = fliplr(vhat);nhat(1) = -nhat(1);
    
    %discretize one dimensional space (finitely many convex combos of
    %x0,xend)
    alpha = 0:1/(len):1; alphanum = numel(alpha);
    
    %make matrix out of alphas (for matrix multiplication)
    alpha_big = repmat(alpha,2,1);
    
    x0_big = repmat(x0',1,alphanum);
    xend_big = repmat(xend',1,alphanum);
    
    for j = 1:(2*w+1)
        
        y = (j-w-1)/2*nhat;
        y_big =repmat(y',1,alphanum);
        
        %compute convex combo of endpoints
        xi = x0_big.*(1 - alpha_big) + alpha_big.*xend_big + y_big;
        
        %no width for now
        xi_round = round(xi);
        
        %store subscripts
        jsub = xi_round(1,:);
        isub = xi_round(2,:);
        
        %make sure subscripts are in image
        jsub(jsub>size(Abw,2)) = size(Abw,2);
        jsub(jsub<1) = 1;
        isub(isub>size(Abw,1)) = size(Abw,1);
        isub(isub<1) = 1;
        
        
        %convert to indices (and not subscripts)
        ind = sub2ind(size(Abw),isub,jsub);
        
        %set pixels to white
        A_nolines(ind) = 255;
        
    end
    
    
    
end


k = 1;

emptystuff = [];
for i = 1:m-1
    for j  =1:n-1
        indx = round(inters_x(i,j)):round(inters_x(i+1,j+1));
        indy = round(inters_y(i,j)):round(inters_y(i+1,j+1));
        
        %         tmp = zeros(size(Abw));
        %         tmp(indy,indx) = Abw(indy,indx);
%         subimages{k} = Abw(indy,indx);
        
        %use subimages of A_nolines (grid is removed)
        subimages{k} = A_nolines(indy,indx);
        
        if isempty(indy) || isempty(indx)
           emptystuff = [emptystuff k]; 
           indy
           indx
        end
        
        k =  k+1;
    end
end

%remove empty subimages
subimages(emptystuff) = [];


%% load images into large matrix
%M IS USED PREVIOUSLY, OVERWRITING
% M = 16;

%put each subimage into a row vector (a single observation)
X = zeros(numel(subimages),storeM^2);

for i = 1:numel(subimages)
    tmp_img = subimages{i};
    
    
    
    %resize image so they're all the same
    imres = imresize(tmp_img,[storeM,storeM]);
    
    %store observation
    X(i,:) = double(imres(:)');
end