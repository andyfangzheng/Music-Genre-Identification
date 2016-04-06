% test 1
clear all; close all; clc;

[black, Fs]=audioread('Soundgarden - Black Hole Sun.mp3');
[path, Fs]=audioread('pathetic Beethoven.mp3');
[beat, Fs]=audioread('Michael Jackson - Beat It.mp3');

%decimate the songs because the original version is too large and 
%my computer can't run the whole code. And I select the first channel
%of the audio
black=decimate(black(1:Fs*260,1),100);
path=decimate(path(1:Fs*260,1),100);
beat=decimate(beat(1:Fs*260,1),100);

% spectrogram
n=5*Fs/100;
t2=linspace(0,5,n+1); t=t2(1:n);
k=(2*pi/5)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);
tslide=0:0.5:5;

black_total=[]; %build empty spectrogram matrix
path_total=[];
beat_total=[];
for j=3:52 %5-second/sample, 50 samples for each song, start from 11s
    black_spec=[];%spectrogram for each sample
    path_spec=[];
    beat_spec=[];
    black_sample=black(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    path_sample=path(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    beat_sample=beat(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    for j2=1:length(tslide)
        g=exp(-10*(t-tslide(j2)).^2); %gabor filter
        black_vg=fft(g.*black_sample);
        black_spec=[black_spec; abs(fftshift(black_vg))];
        path_vg=fft(g.*path_sample);
        path_spec=[path_spec; abs(fftshift(path_vg))];
        beat_vg=fft(g.*beat_sample);
        beat_spec=[beat_spec; abs(fftshift(beat_vg))];
    end
    %put all spactrograms into one matrix for each song
    col_size=size(black_spec,1)*size(black_spec,2);
    black_spec=reshape(black_spec,[col_size,1]);
    black_total=[black_total black_spec];
    col_size=size(path_spec,1)*size(path_spec,2);
    path_spec=reshape(path_spec,[col_size,1]);
    path_total=[path_total path_spec];
    col_size=size(beat_spec,1)*size(beat_spec,2);
    beat_spec=reshape(beat_spec,[col_size,1]);
    beat_total=[beat_total beat_spec];
end

% SVD
X=[black_total path_total beat_total];
[U,S,V] = svd(X,'econ');
Y=U'*X; 

%% LDA
q1=randperm(50);
q2=randperm(50);
q3=randperm(50);
xblack=V(1:50,2:4);
xpath=V(51:100,2:4);
xbeat=V(101:150,2:4);
xtrain=[xblack(q1(1:30),:); xpath(q2(1:30),:); xbeat(q3(1:30),:)];
xtest=[xblack(q1(31:end),:); xpath(q2(31:end),:); xbeat(q3(31:end),:)];
ctrain=[ones(30,1);2*ones(30,1);3*ones(30,1)];
pre=classify(xtest,xtrain,ctrain);
bar(pre);




