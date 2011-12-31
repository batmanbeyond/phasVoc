
%[x,fs]=wavread('Roland-GR-1-Pick-Bass-2-C4.wav');

function [output] = PhaseVoc(x,fs, tCe ,wL, h)
% x is the signal
% fs is the sampling frequency 
% tCe is time compression expansion, x is intially tCe=1
% wL is the size of each window
% h is the hop size

%check wL to be power of two, if not find and adjust
if(wL==2^nextpow2(wL))
    wL=wL;
else 
    wL=2^nextpow2(wL);
end;

%check tCe, must be greater then 0
if(tCe<=0)
    error('tCe needs to be > 0');
end;

%Check input of hopsize, change if entered incorrectly 
if(h<=1)
    h= h;
else
    h=h/100;
end;

%=====================|
%Analysis Parameters  |
%=====================|

sigLen=length(x);             % length of input signal RSR
hW= hanning(wL);              % hanning window the length of parameter wL
r=1; %intialize counter
i=1;
hIa=floor((1-h)*wL); %hop in analysis
hOs=floor(tCe*hIa);  %hop out synthesis 

NumOfWin = ceil(sigLen/(wL*(1-h))); %number of windows, so space 2 b allcoated


PL= sigLen;%*tCe; %potential length of signal
z = zeros(1,(floor(PL))); %allocation of space given potential length of signal
                          %row of zeros of length PL

%for comparison in first loop

win = x(1:1+wL-1).*hW'; %first window for anaylysis
AmpS = fft(win,wL); %amplitude spectrum of first window
PhaS = angle(AmpS); %phase spectrum of first window
tP=2*pi; % make it easy to deal with radians 
k=fs/wL; %bins (still not sure i understand this one)


 
    while(i+wL<=sigLen)
        
%===========|
%Windowing  |
%===========|
      
      
      
       
           win = x(i:i+wL-1).*hW'; %window signal from r to sigLen
          
           subplot(2,1,1)
           plot(x) 
         
              AmpS=(fft(win, wL));  %amplitude spectrum 
              
              PhaS1=angle(AmpS);     %phase Spectrum
              
              AmpS=abs(fft(win,wL));%absolute value of aplitude spectrum
               
           
              
              
          
%===========|
%Encoding   |
%===========|       
    
    
 

 
 m=(((tP*k)/wL)*hIa+PhaS1); %%
 m= round(m/tP);            %%Calculating instantaneous freq of [s+1] 
 unP1= tP*m+PhaS1*wL;       %%
 
 DeltaP= (unP1-(PhaS*hIa))*[o:wL-1]; %Delta Phase ((((PhaS1 compared to PhaS))))

 iF= (tP^-1)*DeltaP*(hIa^-1)*fs; %Instantaneous Frequency 
    % update the phase values for calculating the next difference 
    
  % PhaS= PhaS1;
       
%===========|
%Decoding   |
%===========|
    
    DeltaP= (hOs/hIa)*DeltaP;%%%%%+(PhaS1*hIa);
    
    DeltaP= DeltaP +(PhaS*hIa);
    
  
    tSt = AmpS.*exp(1i*DeltaP);      %%%%iF1);
    xNew  = real(ifft(tSt)).*hW';
    
      PhaS= PhaS1;
   
    subplot(2,1,2) %CHECKING OUPUT (TEST)
    plot(xNew)     % 
  

  % z((i+1):(i+wL)) = z((i+1):(i+wL))+xNew;
z(i:wL+1)=xNew(i:wL+1);
i = i+(wL*(1-h)); %update loop
    
    end
    
sound(z,fs);
%plot(z(1:end));

end

  

