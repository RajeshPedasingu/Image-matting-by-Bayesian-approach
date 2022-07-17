clc
clearvars
clear all
close all
imgo=imread('gandalf-input.png');    %%%original Image
trimap=imread('gandalf-trimap.png');  %%% trimap image same size as original image 
figure;
imshow(imgo);
figure
imshow(trimap);
[M,N]=size(trimap);

%%%Foreground, Background, Unknown

imgo=im2double(imgo);
trimap=im2double(trimap);
fg=zeros(M,N,3);bg=zeros(M,N,3);ug=zeros(M,N,3);
for i=1:M
    for j=1:N
        if trimap(i,j) == 1
            fg(i,j,:)=imgo(i,j,:);
        elseif trimap(i,j) == 0
            bg(i,j,:)=imgo(i,j,:);
        else
            ug(i,j,:)=imgo(i,j,:);
        end
    end
end

% figure('Name', 'FG, BG, UNKNOWN'); 
imshow([uint8(fg*255) uint8(bg*255) uint8(ug*255)]);


%%%%%Compute statistics from known foreground & background 

%Finding mean of known foreground & background
f_mean = [0 0 0];
b_mean = [0 0 0];
for i=1:3
    f = fg(:,:,i);
    b = bg(:,:,i);
    f_mean(i) = mean(f(find(f)));
    b_mean(i) = mean(b(find(b)));
end

%Finding variance of known foreground
f_div(:,:,1) = fg(:,:,1) - f_mean(1);
f_div(:,:,2) = fg(:,:,2) - f_mean(2);
f_div(:,:,3) = fg(:,:,3) - f_mean(3);
sumF = [0 0 0; 0 0 0; 0 0 0];
count = 0;
for c=1:size(f_div, 2)
    for r=1:size(f_div, 1)
        pixF = f_div(r,c, :);
        pixF = reshape(pixF,3,1);
        if(any(fg(r,c,:)))
            sumF = sumF + (pixF * pixF');
            count = count + 1;
        end
    end
end
f_var = sumF / count;

%%%finding variance of known background
b_div(:,:,1) = bg(:,:,1) - b_mean(1);
b_div(:,:,2) = bg(:,:,2) - b_mean(2);
b_div(:,:,3) = bg(:,:,3) - b_mean(3);
sumB = [0 0 0; 0 0 0; 0 0 0];
count = 0;
for c=1:size(b_div, 2)
    for r=1:size(b_div, 1)
        pixB = b_div(r,c, :);
        pixB = reshape(pixB,3,1);
        if(any(bg(r,c,:)))
            sumB = sumB + (pixB * pixB');
            count = count + 1;
        end
    end
end
b_var = sumB / count;

%%%Solve for unknown pixel values and alpha
alpha_un = trimap;
for i=1:M
    for j=1:N
        if(~any(ug(i,j,:)))
            continue
        end
        un_c = reshape(ug(i,j,:), 3,1);
        %initiating alpha as avg of neighbors
        alpha = 0;
        count = 0;
        try
        alpha = alpha + alpha_un(r-1, c);
        count = count + 1;
        end
        try
        alpha = alpha + alpha_un(r+1, c);
        count = count + 1;
        end
        try
        alpha = alpha + alpha_un(r, c-1);
        count = count + 1;
        end
        try
        alpha = alpha + alpha_un(r, c+1);
        count = count + 1;
        end
        alpha = alpha / count;

        %%%%iteratively solve for color of pixel
        for k=1:10
            alpha_prev = alpha;
            [F B] = FandB(f_var, b_var, f_mean, b_mean, un_c, alpha);
            alpha = dot((un_c - B), (F - B)) / norm(F-B).^2;
            if(abs(alpha -alpha_prev) <= 0.0001)
                break;
            end
        end
        alpha_un(i,j) = alpha;
        fg(i,j,:) = F;
        bg(i,j,:) = B;

    end
end

figure('Name', 'Alpha Matte'); imshow(alpha_un);

img3 = imread('gandalf-background.png');
img3=im2double(img3);
img_final=img3;
img_final(:,:,1) = fg(:,:,1).*alpha_un(:,:) + img3(:,:,1).*(1 - alpha_un(:,:)); 
img_final(:,:,2) = fg(:,:,2).*alpha_un(:,:) + img3(:,:,2).*(1 - alpha_un(:,:)); 
img_final(:,:,3) = fg(:,:,3).*alpha_un(:,:) + img3(:,:,3).*(1 - alpha_un(:,:)); 

figure('Name', 'Composed Image'); imshow(uint8(img_final*255));

function [out1,out2]=FandB(f_var, b_var, f_mean, b_mean, un_c, alpha); 
SigmaC_Sq = 64;
f_inv = inv(f_var);
b_inv = inv(b_var);
A_fl = f_inv + (eye(3)*alpha*alpha)/(SigmaC_Sq);
A_fr = (eye(3)*alpha*(1 - alpha))/(SigmaC_Sq);
A_sl = (eye(3)*alpha*(1-alpha))/(SigmaC_Sq);
A_sr = b_inv + (eye(3)*(1-alpha)*(1-alpha))/ (SigmaC_Sq);
A = [A_fl A_fr; A_sl A_sr];
b_l = (f_inv*f_mean') + (un_c*alpha)/(SigmaC_Sq);
b_r = (b_inv*b_mean') + (un_c*(1-alpha))/(SigmaC_Sq);
b = [b_l; b_r];

x = A\b;
F = x(1:3);
B = x(4:6);
out1=F;
out2=B;
end



