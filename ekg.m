% ---------------Part 1: ECG Fundamentals---------------

% Reading the healthy data into MATLAB
clear
healthydata = readtable("lab 3 healthy data.xlsx");
time = healthydata.Time;
L1 = healthydata.LeadI/1000; % Divided by 1000 since data is in microVolts
L2 = healthydata.LeadII/1000;
L3 = healthydata.LeadIII/1000;
avr = healthydata.aVR/1000;
avl = healthydata.aVL/1000;
avf = healthydata.aVF/1000;

% Plot 1: Plotting all 6 ECG leads from healthy data file
figure(1)
plotSixLeads(time,L1,L2,L3,avr,avl,avf)
subplot(6,1,1)
title('Healthy ECG (Original Data)')

% Deriving the four other leads using leads 1 and 2
calculated_L3 = L1 + L2;
calculated_avr = -(L1 + L2) / sqrt(3);
calculated_avl = ((2*L1) - L2) / sqrt(3);
calculated_avf = ((2*L2) - L1) / sqrt(3);

% Plot 2: Plotting leads 1, 2, and the four derived leads
figure(2)
plotSixLeads(time,L1,L2,calculated_L3,calculated_avr,calculated_avl,calculated_avf)
subplot(6,1,1)
title('Healthy ECG (Calculated Values)')


% --------------Part 2: Analysis of Healthy and Diseased Data--------------
% Calculating the healthy mean ECG then detrending it
meanekg = (L1 + L2 + L3)/3;
opol = 6;
[p,s,mu] = polyfit(time,meanekg,opol);
f_y = polyval(p,time,[],mu);
detrend = meanekg - f_y;

% Plot 3: Plotting the mean ECG signal (healthy)
figure(3)
plot(time,detrend);
xlabel('Time (seconds)');
ylabel('Voltage (mV)');
title('Mean ECG (Healthy)');

% Identification of QRS Complex, P and T waves on graph of mean ECG signal (healthy)
% P
Plocs = [0.6 1.425 2.225 3.04 3.875 4.75 5.635 6.535];
Ppks = detrend(Plocs*200);
text(Plocs,Ppks,'P');
% Q
Q = islocalmin(detrend,'MinProminence',0.15);
Qlocs = find(Q==1)./200;
Qpks = detrend(Qlocs*200);
text(Qlocs,Qpks,'Q')
% R
[Rpks, Rlocs] = findpeaks(detrend,time,'MinPeakHeight',0.25,'MinPeakDistance',0.5);
text(Rlocs,Rpks,'R');
% S
S=[164 330 481 650 813 987 1167 1344];
Slocs = S./200;
Spks = detrend(S);
text(Slocs,Spks,'S');
% T
T=[29 198 362 522 688 852 1024 1202 1385];
Tlocs = T./200;
Tpks = detrend(T);
text(Tlocs,Tpks,'T');

% Measured heart rate (healthy)
heart_rate_healthy = 60 / ((max(Rlocs) - min(Rlocs)) / (length(Rlocs) - 1));
% Maximum and minimum voltages (healthy)
max_voltage_healthy = max(detrend);
min_voltage_healthy = min(detrend);
% Average P-Q interval (healthy)
pq_sum = 0;
for i = 1:8
    pq_interval = Qlocs(i) - Plocs(i);
    pq_sum = pq_sum + pq_interval;
end
pq_average_healthy = pq_sum/8;
% Average P-R interval (healthy)
pr_sum = 0;
for j = 1:8
    pr_interval = Rlocs(j) - Plocs(j);
    pr_sum = pr_sum + pr_interval;
end
pr_average_healthy = pr_sum/8;
% Average Q-T interval (healthy)
qt_sum = 0;
for k = 1:8
    qt_interval = Tlocs(k+1) - Qlocs(k);
    qt_sum = qt_sum + qt_interval;
end
qt_average_healthy = qt_sum/8;

%Plot 4: Mean axis of depolarization (healthy), derived from two of the lead channels
max_L1 = max(L1);
max_L2 = max(L2);
L1_x = max_L1;
L1_y = 0;
L2_x = max(L2) * cosd(60);
L2_y = max(L2) * sind(60);
slopeL2= (L2_y/L2_x);
pL2= -slopeL2;
bL2= L2_y -(pL2*L2_x);
Inter= pL2*L1_x+bL2;
 
mea_angle= atand(Inter/L1_x);
figure(4)
grid on
set(gca,'visible','off')
title('Mean Electrical Axis (Healthy)')
pax = polaraxes;
theta = [0,mea_angle*(pi/180)];
rho = [0, max(Rpks)]; %find from height of R complex;
polarplot(theta, rho)
pax.ThetaDir = 'clockwise';
pax.FontSize = 12;


% Diseased Data
diseasedata=readtable("lab 3 disease data.xlsx");
time_3=diseasedata.Time;
dis_L1=diseasedata.LeadI; 
dis_L2=diseasedata.LeadII; 
dis_L3 = diseasedata.LeadIII;
dis_avr = diseasedata.aVR;
dis_avl = diseasedata.aVL;
dis_avf = diseasedata.aVF;

% Plot 5: Plotting all 6 ECG leads from diseased data file
figure(5)
plotSixLeads(time_3,dis_L1,dis_L2,dis_L3,dis_avr,dis_avl,dis_avf)
subplot(6,1,1)
title('Diseased ECG (Original Data)')

% Deriving the four other leads using leads 1 and 2
calculated_dis_L3 = dis_L1 + dis_L2;
calculated_dis_avr = -(dis_L1 + dis_L2) / sqrt(3);
calculated_dis_avl = ((2*dis_L1) - dis_L2) / sqrt(3);
calculated_dis_avf = ((2*dis_L2) - dis_L1) / sqrt(3);

% Plot 6: Plotting leads 1, 2, and the four derived leads
figure(6)
plotSixLeads(time_3,dis_L1,dis_L2,calculated_dis_L3,calculated_dis_avr,calculated_dis_avl,calculated_dis_avf)
subplot(6,1,1)
title('Diseased ECG (Calculated Values)')

%Plot 7: Plotting the mean ECG signal (diseased)
dis_meanekg = (dis_L1 + dis_L2 + dis_L3)/3;
dis_opol = 6;
[dis_p,dis_s,dis_mu] = polyfit(time_3,dis_meanekg,dis_opol);
dis_f_y = polyval(dis_p,time_3,[],dis_mu);
dis_detrend = dis_meanekg - dis_f_y;

figure(7)
plot(time_3,dis_detrend);
xlabel('Time (seconds)');
ylabel('Voltage (V)');
title('Mean ECG (Diseased)');

% Given the erratic nature of the data, only Q and R were able to be identified
% Q
Q_dis_a = islocalmin(dis_detrend,'MinProminence',0.15);
Q_dis = find(Q_dis_a==1);
Qlocs_dis = Q_dis./200;
Qpks_dis = dis_detrend(Q_dis);
text(Qlocs_dis,Qpks_dis,'Q')
% R
[Rpks_dis, Rlocs_dis] = findpeaks(dis_detrend,time_3,'MinPeakHeight',0.15,'MinPeakDistance',0.2);
text(Rlocs_dis,Rpks_dis,'R');
% Measured heart rate
heart_rate_dis = 60 / ((max(Rlocs_dis) - min(Rlocs_dis)) / (length(Rlocs_dis) - 1));
% Maximum and minimum voltages
max_voltage_dis = max(dis_detrend);
min_voltage_dis = min(dis_detrend);

% Plot 8: Mean Electrical Axis, Include lead channel vectors you used in your plot
max_dis_L1 = max(dis_L1);
max_dis_L2 = max(dis_L2);
dis_L1_x = max_dis_L1;
dis_L1_y = 0;
dis_L2_x = max(dis_L2) * cosd(60);
dis_L2_y = max(dis_L2) * sind(60);
dis_slopeL2= (dis_L2_y/dis_L2_x);
dis_pL2= -dis_slopeL2;
dis_bL2= dis_L2_y -(dis_pL2*dis_L2_x);
dis_Inter= dis_pL2*dis_L1_x+dis_bL2;
dis_mea_angle=atand(dis_Inter/dis_L1_x);
 
figure(8)
grid on
set(gca,'visible','off')
title('Mean Electrical Axis (Diseased)')
pax = polaraxes;
theta = [0,dis_mea_angle*(pi/180)];
rho = [0, max(Rpks_dis)]; %find from height of R complex;
polarplot(theta, rho)
pax.ThetaDir = 'clockwise';
pax.FontSize = 12;


% Table
Measurement = ["Measured Heart Rate (bpm)";"Maximum Voltage (V)";"Minimum Voltage (V)";...
    "Average P-Q Interval (s)";"Average P-R Interval (s)";"Average Q-T Interval (s)";"Mean Electrical Axis (Θ)"];
HealthyState = [heart_rate_healthy;max_voltage_healthy;min_voltage_healthy;...
    pq_average_healthy;pr_average_healthy;qt_average_healthy;mea_angle];
DiseasedState = [heart_rate_dis;max_voltage_dis;min_voltage_dis;...
    "N/A";"N/A";"N/A";dis_mea_angle];
table = table(Measurement,HealthyState,DiseasedState)

% ---------------Part 3: Testing the Code---------------
% Define thresholds for classification
normal_threshold = [60 100]; 
tachycardia_threshold = 100;
bradycardia_threshold = 60;
heart_rate2 = 60 / mean(dis_meanekg);

% Classify heart rate
if heart_rate2 < normal_threshold(1) 
    disp('The patient has bradycardia.');
elseif heart_rate2 > normal_threshold(2) && dis_mea_angle < -30 && dis_mea_angle > -90
    disp('The patient has a left axis deviation, potentially ventricular tachycardia.')
elseif heart_rate2 > normal_threshold(2) && dis_mea_angle > 90 && dis_mea_angle < 180
    disp('The patient has a right axis deviation, potentially ventricular tachycardia.')
elseif heart_rate2 > normal_threshold(2) && dis_mea_angle > -180 && dis_mea_angle < -90
    disp('The patient has tachycardia and extreme axis deviation, further testing required!')
elseif heart_rate2 > normal_threshold(2)
    disp('The patient has tachycardia.');
else
    disp('The heart rate is within the normal range.');
end



%% Graphing ECG function
function plotSixLeads(t,L1,L2,L3,aVR,aVL,aVF)
%This function generates 6 subplots in a column, intended to contain the 6
%ECG leads' data over time. Inputs are time array, Leads 1-3 arrays, AVR,
%AVL, AVF array
subplot(6,1,1)
plot(t,L1)
xlabel('Time (seconds)')
ylabel(['Lead 1', newline, 'Voltage (mV)'],'FontSize',10)
subplot(6,1,2)
plot(t,L2)
xlabel('Time (seconds)')
ylabel(['Lead 2', newline, 'Voltage (mV)'],'FontSize',10)
subplot(6,1,3)
plot(t,L3)
xlabel('Time (seconds)')
ylabel(['Lead 3', newline, 'Voltage (mV)'],'FontSize',10)
subplot(6,1,4)
plot(t,aVR)
xlabel('Time (seconds)')
ylabel(['aVR', newline, 'Voltage (mV)'],'FontSize',10)
subplot(6,1,5)
plot(t,aVL)
xlabel('Time (seconds)')
ylabel(['aVL', newline, 'Voltage (mV)'],'FontSize',10)
subplot(6,1,6)
plot(t,aVF)
xlabel('Time (seconds)')
ylabel(['aVF', newline, 'Voltage (mV)'],'FontSize',10)
end
