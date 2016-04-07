% test 3
clear all; close all; clc;

%jazz
[jazz1, Fs]=audioread('Billie Holiday - Lover Man.mp3');
[jazz2, Fs2]=audioread('Georgia on my Mind- Ray Charles.mp3');
[jazz3, Fs3]=audioread('Count Basie - April In Paris.mp3');
jazz1=decimate(jazz1(1:Fs*110,1),100);
jazz2=decimate(jazz2(1:Fs*110,1),100);
jazz3=decimate(jazz3(1:Fs*110,1),100);

%classical
[classical1, Fs4]=audioread('pathetic Beethoven.mp3');
[classical2, Fs5]=audioread('Chopin - Nocturne.mp3');
[classical3, Fs6]=audioread('Beethoven symphony 9 Movement2.mp3');
classical1=decimate(classical1(1:Fs*110,1),100);
classical2=decimate(classical2(1:Fs*110,1),100);
classical3=decimate(classical3(1:Fs*110,1),100);

%rock
[rock1, Fs7]=audioread('AC-DC - Thunderstruck.mp3');
[rock2, Fs8]=audioread('AC-DC - Whole Lotta Rosie.mp3');
[rock3, Fs9]=audioread('AC-DC - For Those About To Rock.mp3');
rock1=decimate(rock1(1:Fs*110,1),100);
rock2=decimate(rock2(1:Fs*110,1),100);
rock3=decimate(rock3(1:Fs*110,1),100);
%%
% spectrogram
n=5*Fs/100;
t2=linspace(0,5,n+1); t=t2(1:n);
k=(2*pi/5)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);
tslide=0:0.5:5;

jazz1_total=[]; %build empty spectrogram matrix
jazz2_total=[];
jazz3_total=[];
classical1_total=[]; 
classical2_total=[];
classical3_total=[];
rock1_total=[]; 
rock2_total=[];
rock3_total=[];
for j=3:22 %5-second/sample, 20 samples for each song, start from 11s
    %for jazz
    jazz1_spec=[];%spectrogram for each sample
    jazz2_spec=[];
    jazz3_spec=[];
    jazz1_sample=jazz1(5*(j-1)*Fs/100+1:+5*j*Fs/100)';
    jazz2_sample=jazz2(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    jazz3_sample=jazz3(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    for j2=1:length(tslide)
        g=exp(-10*(t-tslide(j2)).^2); %gabor filter
        jazz1_vg=fft(g.*jazz1_sample);
        jazz1_spec=[jazz1_spec; abs(fftshift(jazz1_vg))];
        jazz2_vg=fft(g.*jazz2_sample);
        jazz2_spec=[jazz2_spec; abs(fftshift(jazz2_vg))];
        jazz3_vg=fft(g.*jazz3_sample);
        jazz3_spec=[jazz3_spec; abs(fftshift(jazz3_vg))];
    end
    %put all spactrograms into one matrix for each song
    col_size=size(jazz1_spec,1)*size(jazz1_spec,2);
    jazz1_spec=reshape(jazz1_spec,[col_size,1]);
    jazz1_total=[jazz1_total jazz1_spec];
    col_size=size(jazz2_spec,1)*size(jazz2_spec,2);
    jazz2_spec=reshape(jazz2_spec,[col_size,1]);
    jazz2_total=[jazz2_total jazz2_spec];
    col_size=size(jazz3_spec,1)*size(jazz3_spec,2);
    jazz3_spec=reshape(jazz3_spec,[col_size,1]);
    jazz3_total=[jazz3_total jazz3_spec];
    %for classical
    classical1_spec=[];
    classical2_spec=[];
    classical3_spec=[];
    classical1_sample=classical1(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    classical2_sample=classical2(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    classical3_sample=classical3(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    for j2=1:length(tslide)
        g=exp(-10*(t-tslide(j2)).^2); 
        classical1_vg=fft(g.*classical1_sample);
        classical1_spec=[classical1_spec; abs(fftshift(classical1_vg))];
        classical2_vg=fft(g.*classical2_sample);
        classical2_spec=[classical2_spec; abs(fftshift(classical2_vg))];
        classical3_vg=fft(g.*classical3_sample);
        classical3_spec=[classical3_spec; abs(fftshift(classical3_vg))];
    end
    col_size=size(classical1_spec,1)*size(classical1_spec,2);
    classical1_spec=reshape(classical1_spec,[col_size,1]);
    classical1_total=[classical1_total classical1_spec];
    col_size=size(classical2_spec,1)*size(classical2_spec,2);
    classical2_spec=reshape(classical2_spec,[col_size,1]);
    classical2_total=[classical2_total classical2_spec];
    col_size=size(classical3_spec,1)*size(classical3_spec,2);
    classical3_spec=reshape(classical3_spec,[col_size,1]);
    classical3_total=[classical3_total classical3_spec];
    %for rock
    rock1_spec=[];
    rock2_spec=[];
    rock3_spec=[];
    rock1_sample=rock1(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    rock2_sample=rock2(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    rock3_sample=rock3(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    for j2=1:length(tslide)
        g=exp(-10*(t-tslide(j2)).^2);
        rock1_vg=fft(g.*rock1_sample);
        rock1_spec=[rock1_spec; abs(fftshift(rock1_vg))];
        rock2_vg=fft(g.*rock2_sample);
        rock2_spec=[rock2_spec; abs(fftshift(rock2_vg))];
        rock3_vg=fft(g.*rock3_sample);
        rock3_spec=[rock3_spec; abs(fftshift(rock3_vg))];
    end
    col_size=size(rock1_spec,1)*size(rock1_spec,2);
    rock1_spec=reshape(rock1_spec,[col_size,1]);
    rock1_total=[rock1_total rock1_spec];
    col_size=size(rock2_spec,1)*size(rock2_spec,2);
    rock2_spec=reshape(rock2_spec,[col_size,1]);
    rock2_total=[rock2_total rock2_spec];
    col_size=size(rock3_spec,1)*size(rock3_spec,2);
    rock3_spec=reshape(rock3_spec,[col_size,1]);
    rock3_total=[rock3_total rock3_spec];
end
jazz_total=[jazz1_total jazz2_total jazz3_total];
classical_total=[classical1_total classical2_total classical3_total];
rock_total=[rock1_total rock2_total rock3_total];

% SVD
X=[jazz_total classical_total rock_total];
[U,S,V] = svd(X,'econ');
Y=U'*X; 

%% LDA
q1=randperm(60);
q2=randperm(60);
q3=randperm(60);
xjazz=V(1:60,2:4);
xclassical=V(61:120,2:4);
xrock=V(121:180,2:4);
xtrain=[xjazz(q1(1:40),:); xclassical(q2(1:40),:); xrock(q3(1:40),:)];
xtest=[xjazz(q1(41:end),:); xclassical(q2(41:end),:); xrock(q3(41:end),:)];
ctrain=[ones(40,1);2*ones(40,1);3*ones(40,1)];
pre=classify(xtest,xtrain,ctrain);
bar(pre);
