% test 2
clear all; close all; clc;

[black, Fs]=audioread('Soundgarden - Black Hole Sun.mp3');
[alice, Fs]=audioread('Alice In Chains - Man In The Box.mp3');
[pearl, Fs]=audioread('Pearl Jam - Nothingman.mp3');

%decimate the songs because the original version is too large and 
%my computer can't run the whole code. And I select the first channel
%of the audio
black=decimate(black(1:Fs*260,1),100);
alice=decimate(alice(1:Fs*260,1),100);
pearl=decimate(pearl(1:Fs*260,1),100);

% spectrogram
n=5*Fs/100;
t2=linspace(0,5,n+1); t=t2(1:n);
k=(2*pi/5)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);
tslide=0:0.5:5;

black_total=[]; %build empty spectrogram matrix
alice_total=[];
pearl_total=[];
for j=3:52 %5-second/sample, 50 samples for each song, start from 11s
    black_spec=[];%spectrogram for each sample
    alice_spec=[];
    pearl_spec=[];
    black_sample=black(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    alice_sample=alice(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    pearl_sample=pearl(5*(j-1)*Fs/100+1:5*j*Fs/100)';
    for j2=1:length(tslide)
        g=exp(-10*(t-tslide(j2)).^2); %gabor filter
        black_vg=fft(g.*black_sample);
        black_spec=[black_spec; abs(fftshift(black_vg))];
        alice_vg=fft(g.*alice_sample);
        alice_spec=[alice_spec; abs(fftshift(alice_vg))];
        pearl_vg=fft(g.*pearl_sample);
        pearl_spec=[pearl_spec; abs(fftshift(pearl_vg))];
    end
    %put all spactrograms into one matrix for each song
    col_size=size(black_spec,1)*size(black_spec,2);
    black_spec=reshape(black_spec,[col_size,1]);
    black_total=[black_total black_spec];
    col_size=size(alice_spec,1)*size(alice_spec,2);
    alice_spec=reshape(alice_spec,[col_size,1]);
    alice_total=[alice_total alice_spec];
    col_size=size(pearl_spec,1)*size(pearl_spec,2);
    pearl_spec=reshape(pearl_spec,[col_size,1]);
    pearl_total=[pearl_total pearl_spec];
end

% SVD
X=[black_total alice_total pearl_total];
[U,S,V] = svd(X,'econ');
Y=U'*X; 

%% LDA
q1=randperm(50);
q2=randperm(50);
q3=randperm(50);
xblack=V(1:50,2:4);
xalice=V(51:100,2:4);
xpearl=V(101:150,2:4);
xtrain=[xblack(q1(1:30),:); xalice(q2(1:30),:); xpearl(q3(1:30),:)];
xtest=[xblack(q1(31:end),:); xalice(q2(31:end),:); xpearl(q3(31:end),:)];
ctrain=[ones(30,1);2*ones(30,1);3*ones(30,1)];
pre=classify(xtest,xtrain,ctrain);
bar(pre);
