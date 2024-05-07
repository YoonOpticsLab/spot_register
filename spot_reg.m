% These can be used instead of file dialog below
% stack_input = 'c:/code/spot_register/stack.avi_reg.avi';
% stack_output = 'c:/code/spot_register/registered100.avi';

% All pixels below this set to zero
NOISE_THRESHOLD=15;

% Finding spots to guess lenslet spacing
MAX_FRAC=0.5;
MIN_PEAK_DISTANCE=20;

% Minimum # of pixels in box to qualify as spot (for binarization)
BOX_THRESHOLD=100;

%%
[file,location] = uigetfile("*.avi");
stack_input = [location file ];
stack_output = [location 'reg_' file ];

v_in = VideoReader(stack_input);

clear frames;
nframe=1;

while hasFrame(v_in)
    frame = readFrame(v_in);
    frames(nframe,:,:) = frame(:,:,2);
    nframe=nframe+1
end

frames_orig = frames;
%close(v_in);

nframes=size(frames,1);
%% THRESHOLD

frames(frames<NOISE_THRESHOLD)=0;
%% Find peak in mean img: this should be a good spot center
mean_img = squeeze( mean(frames,1));
[val,loc]=max(mean_img(:) );

[maxy,maxx] = ind2sub(size(mean_img), loc);

%% Find peaks in mean x and mean y directions to determine spot spacing

meansigX=squeeze( mean( frames, [1 2] ));
meansigY=squeeze( mean( frames, [1 3] ));

minheight = max(meansigX)*MAX_FRAC;
[pks,locsX] = findpeaks(meansigX,'MinPeakHeight',minheight,'MinPeakDistance',MIN_PEAK_DISTANCE);
[pks,locsY] = findpeaks(meansigY,'MinPeakHeight',minheight,'MinPeakDistance',MIN_PEAK_DISTANCE);
spacing=mean( [diff(locsX)', diff(locsY)] );

%% Position boxes so that the good spot (max pixel) is centered in a box
corners_y=[1:floor( size(frames,2)/spacing) ] * spacing-spacing;
corners_x=[1:floor( size(frames,3)/spacing) ] * spacing-spacing;

corners_x = int32( floor(corners_x + mod(maxx,spacing) + spacing/2) );
corners_y = int32( floor(corners_y + mod(maxy,spacing) + spacing/2) );

lenslets_x = size(corners_x,2);
lenslets_y = size(corners_y,2);

figure; imagesc( mean_img); hold on;
plot(maxx,maxy, 'ro'); hold on;

global X
global Y
global box1

[X,Y]=meshgrid( 1:lenslets_x, 1:lenslets_y );

for y1=1:size(corners_y,2)
    plot( [1,size(frames,3)], [corners_y(y1), corners_y(y1)], 'r-'); hold on;
end

for x1=1:size(corners_x,2)
    plot(  [corners_x(x1), corners_x(x1)], [1,size(frames,2)],'r-'); hold on;
end

opts_all=zeros( nframes, 3 );

%% Make sum within each lenslet grid
stride=ceil(spacing);
box1=ones( [stride,stride] );

for nframe=1:size(frames,1)
    for x1=1:lenslets_x
        for y1=1:lenslets_y
            bottom_y = min(corners_y(y1)+stride, size(frames,2)); % clamp to maxy
            right_x = min(corners_x(x1)+stride, size(frames,3)); % clamp to maxx
         
            box1=frames(nframe,corners_y(y1):bottom_y,corners_x(x1):right_x);
            box_sums(y1,x1)=sum(box1(:));
        end
    end

    boxes_all(nframe,:,:)=box_sums;
end

binarized_boxes = boxes_all > BOX_THRESHOLD*1.0;

%% Fit
opts = zeros( size(frames,1), 3);
for nframe = [1:size(frames,1)]
    box1=squeeze(binarized_boxes(nframe,:,:));
    
    % As initial guess, use center of mass
    guess_x = sum( box1 .* X, [1 2] ) / sum(box1(:));
    guess_y = sum( box1 .* Y, [1 2] ) / sum(box1(:));
    guess_rad=sqrt(sum(box1(:))/3.1416);
    guess=[guess_x, guess_y, guess_rad];

    % Minimize "mycirc" function (LSQERR circle fitting vs. "box1")
    fitted=fminsearch(@circle,guess);
    err= guess - fitted;

    opts(nframe,:) = fitted;
end

%% Shift and write video
v_out = VideoWriter(stack_output,'Grayscale AVI');
v_out.FrameRate = v_in.FrameRate;
open(v_out);

for iframe=1:size(frames,1)
    yshift=round( (lenslets_y/2 - opts(iframe,2))) * stride;
    xshift=round( (lenslets_x/2 - opts(iframe,1))) * stride ;

    frame1=squeeze(frames(iframe,:,:));

    if (abs(xshift)<size(frames,3)) && (abs(yshift)<size(frames,2))
        frame1 = circshift(frame1,[yshift xshift]);
    end

    frame1=uint8(frame1);
    writeVideo(v_out, frame1);
    iframe
end
close(v_out);
disp(['OK, wrote ' num2str(nframes) ' frames to ' stack_output ' from ' stack_input ]);
%%
function y=mycirc(pees) %cx,cy,rad)
    global X
    global Y
    global box1
    cx=pees(1);
    cy=pees(2);
    rad=pees(3);
    %[X,Y]=meshgrid( [1:size(x,2)], [1:size(x,1)] );
    %r=sqrt((X-cx-0.0*sign(X))^2+(Y-cy-0.0*sign(Y))^2); % 0.5 fudge
    r=sqrt((X-cx).^2+(Y-cy).^2);
    y=(r<rad).*1.0;
    y=sum( (y(:)-box1(:)).^2 );
end

