
%[x,fs]=wavread('Roland-GR-1-Pick-Bass-2-C4.wav');

function [z] = PhasVoc(x,fs, tCe ,wL, h)
% x is the signal
% fs is the sampling frequency 
% tCe is time compression expansion, x is intially tCe=1
% wL is the size of each window
% h is the hop size



%make the signal mono
x=x(:,1);

%check wL to be power of two, if not find and adjust
if(wL==2^nextpow2(wL))
    wL=wL;
else 
    wL=2^nextpow2(wL);
end

%check tCe, must be greater then 0
if(tCe<=0)
    error('tCe needs to be > 0');
end

%Check input of hopsize, change if entered incorrectly 
if(h<=1)
    h=h;
else
    h=h/100;
end

%=====================|
%Analysis Parameters  |
%=====================|

sigLen=length(x);             % length of input signal (Column)

hW= hanning(wL);              % hanning window the length of parameter wL


l=0;
r=0;


hIa=floor((1-h)*wL); %hop in analysis
hOs=floor(tCe*hIa);  %hop out synthesis 

PL= (sigLen*tCe)+wL; %potential length of signal
z = zeros(wL+ceil(PL),1); %allocation of space given potential length of signal
                               %column of zeros of length PL


 
    while l<length(x)-wL
        
%===========|
%Windowing  |
%===========|
    
              win = x(l+1:l+wL).*hW; %window, each element multiplied 
                                         %by hanning window equivalent
              ft=fft(win);  %amplitude  
              
              ft=abs(ft); %magnitude
              
              Phas=tCe*angle(ft);     %phase 
              
            
          
%===========|
%Encoding   |
%===========|       
    
    
 ft = ft.*exp(1i*Phas); 
 
 
%===========|
%Decoding   |
%===========|


    xNew  = real(ifft(ft)).*hW;
    
    
   
     
  

  z(r+1:r+wL) = z(r+1:r+wL)+xNew; %resythesize with hOs paraemter
 
  l = l+hIa; %update hIa parameter
  r = r+hOs; %update hOs parameter
    
     
    
    end
    
sound(z,fs);
plot(z(1:end));

end

  

