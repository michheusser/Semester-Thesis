close all
clear all
clc

Number_of_datasets = 18;
Fs = 44.1e3;
%% Feedforward (Chin-Ear) - Data Acquisition
%Ricardo
disp('Loading audio (.wav) FF files...')
Ricardo_FF_Sweep1_Data = wavread('.\Files\WAVs\RICARDO\FF\RICARDO_SWEEP1_Intensidad2.wav');
Ricardo_FF_Sweep1 = struct('Data',Ricardo_FF_Sweep1_Data,'Fs',Fs);
Ricardo_FF_Sweep2_Data = wavread('.\Files\WAVs\RICARDO\FF\RICARDO_SWEEP2_Intensidad2.wav');
Ricardo_FF_Sweep2 = struct('Data',Ricardo_FF_Sweep2_Data,'Fs',Fs);
Ricardo_FF_Noise_Data = wavread('.\Files\WAVs\RICARDO\FF\RICARDO_RUIDOBLANCOFILTRADO_Intensidad2.wav');
Ricardo_FF_Noise = struct('Data',Ricardo_FF_Noise_Data,'Fs',Fs);

%Patricio
Patricio_FF_Sweep1_Data = wavread('.\Files\WAVs\PATRICIO\FF\PATRICIO_SWEEP1_Intensidad2.wav');
Patricio_FF_Sweep1 = struct('Data',Patricio_FF_Sweep1_Data,'Fs',Fs);
Patricio_FF_Sweep2_Data = wavread('.\Files\WAVs\PATRICIO\FF\PATRICIO_SWEEP2_Intensidad2.wav');
Patricio_FF_Sweep2 = struct('Data',Patricio_FF_Sweep2_Data,'Fs',Fs);
Patricio_FF_Noise_Data = wavread('.\Files\WAVs\PATRICIO\FF\PATRICIO_RUIDOBLANCOFILTRADO_Intensidad2.wav');
Patricio_FF_Noise = struct('Data',Patricio_FF_Noise_Data,'Fs',Fs);

%Demian
Demian_FF_Sweep1_Data = wavread('.\Files\WAVs\DEMIAN\FF\DEMIAN_SWEEP1_Intensidad2.wav');
Demian_FF_Sweep1 = struct('Data',Demian_FF_Sweep1_Data,'Fs',Fs);
Demian_FF_Sweep2_Data = wavread('.\Files\WAVs\DEMIAN\FF\DEMIAN_SWEEP2_Intensidad2.wav');
Demian_FF_Sweep2 = struct('Data',Demian_FF_Sweep2_Data,'Fs',Fs);
Demian_FF_Noise_Data = wavread('.\Files\WAVs\DEMIAN\FF\DEMIAN_RUIDOBLANCOFILTRADO_Intensidad2.wav');
Demian_FF_Noise = struct('Data',Demian_FF_Noise_Data,'Fs',Fs);

disp('Audio (.wav) FF files loaded')
%% FF Struct
disp('Creating FF struct...')
ricardo_FF = struct('Sweep1',Ricardo_FF_Sweep1,'Sweep2',Ricardo_FF_Sweep2,'Noise',Ricardo_FF_Noise);
patricio_FF = struct('Sweep1',Patricio_FF_Sweep1,'Sweep2',Patricio_FF_Sweep2,'Noise',Patricio_FF_Noise);
demian_FF = struct('Sweep1',Demian_FF_Sweep1,'Sweep2',Demian_FF_Sweep2,'Noise',Demian_FF_Noise);

Data_FF = struct('ricardo',ricardo_FF,'patricio',patricio_FF,'demian',demian_FF);

disp('FF struct created.')
%% FB (Earphone-Ear) - Data Acquisition

disp('Loading audio (.wav) FB files...')
%Ricardo
Ricardo_FB_Sweep1_Data = wavread('.\Files\WAVs\RICARDO\EA\RICARDOELECTRO_SWEEP1_TOMA1.wav');
Ricardo_FB_Sweep1 = struct('Data',Ricardo_FB_Sweep1_Data,'Fs',Fs);
Ricardo_FB_Sweep2_Data = wavread('.\Files\WAVs\RICARDO\EA\RICARDOELECTRO_SWEEP2_TOMA1.wav');
Ricardo_FB_Sweep2 = struct('Data',Ricardo_FB_Sweep2_Data,'Fs',Fs);
Ricardo_FB_Noise_Data = wavread('.\Files\WAVs\RICARDO\EA\RICARDOELECTRO_RUIDOBLANCOFILTRADO_TOMA1.wav');
Ricardo_FB_Noise = struct('Data',Ricardo_FB_Noise_Data,'Fs',Fs);

%Patricio
Patricio_FB_Sweep1_Data = wavread('.\Files\WAVs\PATRICIO\EA\PATRICIOELECTRO_SWEEP1_TOMA2.wav');
Patricio_FB_Sweep1 = struct('Data',Patricio_FB_Sweep1_Data,'Fs',Fs);
Patricio_FB_Sweep2_Data = wavread('.\Files\WAVs\PATRICIO\EA\PATRICIOELECTRO_SWEEP2_TOMA2.wav');
Patricio_FB_Sweep2 = struct('Data',Patricio_FB_Sweep2_Data,'Fs',Fs);
Patricio_FB_Noise_Data = wavread('.\Files\WAVs\PATRICIO\EA\PATRICIOELECTRO_RUIDOBLANCOFILTRADO_TOMA2.wav');
Patricio_FB_Noise = struct('Data',Patricio_FB_Noise_Data,'Fs',Fs);

%Demian
Demian_FB_Sweep1_Data = wavread('.\Files\WAVs\DEMIAN\EA\RICARDOELECTRO_SWEEP1_TOMA1.wav');
Demian_FB_Sweep1 = struct('Data',Demian_FB_Sweep1_Data,'Fs',Fs);
Demian_FB_Sweep2_Data = wavread('.\Files\WAVs\DEMIAN\EA\RICARDOELECTRO_SWEEP2_TOMA1.wav');
Demian_FB_Sweep2 = struct('Data',Demian_FB_Sweep2_Data,'Fs',Fs);
Demian_FB_Noise_Data = wavread('.\Files\WAVs\DEMIAN\EA\RICARDOELECTRO_RUIDOBLANCOFILTRADO_TOMA1.wav');
Demian_FB_Noise = struct('Data',Demian_FB_Noise_Data,'Fs',Fs);

disp('Audio (.wav) FB files loaded')
%% FB Struct
disp('Creating FB struct...')
ricardo_FB = struct('Sweep1',Ricardo_FB_Sweep1,'Sweep2',Ricardo_FB_Sweep2,'Noise',Ricardo_FB_Noise);
patricio_FB = struct('Sweep1',Patricio_FB_Sweep1,'Sweep2',Patricio_FB_Sweep2,'Noise',Patricio_FB_Noise);
demian_FB = struct('Sweep1',Demian_FB_Sweep1,'Sweep2',Demian_FB_Sweep2,'Noise',Demian_FB_Noise);

Data_FB = struct('ricardo',ricardo_FB,'patricio',patricio_FB,'demian',demian_FB);
disp('FB struct created.')
%% Main Struct
disp('Creating main data struct')
Data = struct('FF',Data_FF,'FB',Data_FB);
save('Data_Struct.mat','Data');
disp('data main struct created.')
clear Data_FF Data_FB
%% Data Matrix
disp('Creating data matrix...')
%Ricardo
RFFS1 = Data.FF.ricardo.Sweep1;
RFFS2 = Data.FF.ricardo.Sweep2;
RFFN = Data.FF.ricardo.Noise;

RFBS1 = Data.FB.ricardo.Sweep1;
RFBS2 = Data.FB.ricardo.Sweep2;
RFBN = Data.FB.ricardo.Noise;

%Patricio
PFFS1 = Data.FF.patricio.Sweep1;
PFFS2 = Data.FF.patricio.Sweep2;
PFFN = Data.FF.patricio.Noise;

PFBS1 = Data.FB.patricio.Sweep1;
PFBS2 = Data.FB.patricio.Sweep2;
PFBN = Data.FB.patricio.Noise;

%Demian
DFFS1 = Data.FF.demian.Sweep1;
DFFS2 = Data.FF.demian.Sweep2;
DFFN = Data.FF.demian.Noise;

DFBS1 = Data.FB.demian.Sweep1;
DFBS2 = Data.FB.demian.Sweep2;
DFBN = Data.FB.demian.Noise;

clear Data


Max_Val = max([     length(RFFS1) length(RFFS2) length(RFFN)...
                    length(RFBS1) length(RFBS2) length(RFBN)...
                    length(PFFS1) length(PFFS2) length(PFFN)...
                    length(PFBS1) length(PFBS2) length(PFBN)...
                    length(DFFS1) length(DFFS2) length(DFFN)...
                    length(DFBS1) length(DFBS2) length(DFBN)]);

Data_Matrix = zeros(Max_Val,Number_of_datasets*2);
Data_Matrix_Length = zeros(Number_of_datasets,1);

%Ricardo
Data_Matrix(1:length(RFFS1),1)     = RFFS1(:,1);
Data_Matrix(1:length(RFFS1),2)     = RFFS1(:,2);
Data_Matrix_Length(1) = length(RFFS1);
Data_Matrix(1:length(RFFS2),3)     = RFFS2(:,1);
Data_Matrix(1:length(RFFS2),4)     = RFFS2(:,2);
Data_Matrix_Length(2) = length(RFFS2);
Data_Matrix(1:length(RFFN) ,5)     = RFFN(:,1);
Data_Matrix(1:length(RFFN) ,6)     = RFFN(:,2);
Data_Matrix_Length(3) = length(RFFN);

Data_Matrix(1:length(RFBS1),7)     = RFBS1(:,1);
Data_Matrix(1:length(RFBS1),8)     = RFBS1(:,2);
Data_Matrix_Length(4) = length(RFBS1);
Data_Matrix(1:length(RFBS2),9)     = RFBS2(:,1);
Data_Matrix(1:length(RFBS2),10)    = RFBS2(:,2);
Data_Matrix_Length(5) = length(RFBS2);
Data_Matrix(1:length(RFBN) ,11)    = RFBN(:,1);
Data_Matrix(1:length(RFBN) ,12)    = RFBN(:,2);
Data_Matrix_Length(6) = length(RFBN);

%Patricio
Data_Matrix(1:length(PFFS1),13)    = PFFS1(:,1);
Data_Matrix(1:length(PFFS1),14)    = PFFS1(:,2);
Data_Matrix_Length(7) = length(PFFS1);
Data_Matrix(1:length(PFFS2),15)    = PFFS2(:,1);
Data_Matrix(1:length(PFFS2),16)    = PFFS2(:,2);
Data_Matrix_Length(8) = length(PFFS2);
Data_Matrix(1:length(PFFN) ,17)    = PFFN(:,1);
Data_Matrix(1:length(PFFN) ,18)    = PFFN(:,2);
Data_Matrix_Length(9) = length(PFFN);

Data_Matrix(1:length(PFBS1),19)    = PFBS1(:,1);
Data_Matrix(1:length(PFBS1),20)    = PFBS1(:,2);
Data_Matrix_Length(10) = length(PFBS1);
Data_Matrix(1:length(PFBS2),21)    = PFBS2(:,1);
Data_Matrix(1:length(PFBS2),22)    = PFBS2(:,2);
Data_Matrix_Length(11) = length(PFBS2);
Data_Matrix(1:length(PFBN) ,23)    = PFBN(:,1);
Data_Matrix(1:length(PFBN) ,24)    = PFBN(:,2);
Data_Matrix_Length(12) = length(PFBN);

%Demian
Data_Matrix(1:length(DFFS1),25)    = DFFS1(:,1);
Data_Matrix(1:length(DFFS1),26)    = DFFS1(:,2);
Data_Matrix_Length(13) = length(DFFS1);
Data_Matrix(1:length(DFFS2),27)    = DFFS2(:,1);
Data_Matrix(1:length(DFFS2),28)    = DFFS2(:,2);
Data_Matrix_Length(14) = length(DFFS2);
Data_Matrix(1:length(DFFN) ,29)    = DFFN(:,1);
Data_Matrix(1:length(DFFN) ,30)    = DFFN(:,2);
Data_Matrix_Length(15) = length(DFFN);

Data_Matrix(1:length(DFBS1),31)    = DFBS1(:,1);
Data_Matrix(1:length(DFBS1),32)    = DFBS1(:,2);
Data_Matrix_Length(16) = length(DFBS1);
Data_Matrix(1:length(DFBS2),33)    = DFBS2(:,1);
Data_Matrix(1:length(DFBS2),34)    = DFBS2(:,2);
Data_Matrix_Length(17) = length(DFBS2);
Data_Matrix(1:length(DFBN) ,35)    = DFBN(:,1);
Data_Matrix(1:length(DFBN) ,36)    = DFBN(:,2);
Data_Matrix_Length(18) = length(DFBN);

save('Data_Matrix.mat','Data_Matrix');
save('Data_Matrix_Length.mat','Data_Matrix_Length');
disp('Data matrix created')

%% Data Struct Cell
Data_Cell = cell(6,1);

RFFS1 = Data.FF.ricardo.Sweep1;
RFFS2 = Data.FF.ricardo.Sweep2;
RFFN = Data.FF.ricardo.Noise;

RFBS1 = Data.FB.ricardo.Sweep1;
RFBS2 = Data.FB.ricardo.Sweep2;
RFBN = Data.FB.ricardo.Noise;

%Patricio
PFFS1 = Data.FF.patricio.Sweep1;
PFFS2 = Data.FF.patricio.Sweep2;
PFFN = Data.FF.patricio.Noise;

PFBS1 = Data.FB.patricio.Sweep1;
PFBS2 = Data.FB.patricio.Sweep2;
PFBN = Data.FB.patricio.Noise;

%Demian
DFFS1 = Data.FF.demian.Sweep1;
DFFS2 = Data.FF.demian.Sweep2;
DFFN = Data.FF.demian.Noise;

DFBS1 = Data.FB.demian.Sweep1;
DFBS2 = Data.FB.demian.Sweep2;
DFBN = Data.FB.demian.Noise;


