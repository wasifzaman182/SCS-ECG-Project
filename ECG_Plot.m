clc; clear; close all;

load('101m.mat');
ECGsignal = val(1,:) - mean(val(1,:));
Fs = 360; t = (0:length(ECGsignal)-1) / Fs;
threshold = max(ECGsignal)/3;

% QRS Detection
[pks,locs] = findpeaks(ECGsignal, 'MinPeakHeight', threshold, 'MinPeakDistance', Fs*0.2);
qrsWidth = round(0.1 * Fs);
q_locs = max(locs - qrsWidth / 2, 1); 
s_locs = min(locs + qrsWidth / 2, length(ECGsignal));

% Arrhythmia Analysis
RR_intervals = diff(locs)/Fs;
heart_rates = 60 ./ RR_intervals;
disp(['Average Heart Rate: ', num2str(mean(heart_rates)), ' BPM']);
analyzeArrhythmia(t, locs, RR_intervals);

% ST Segment and PR Interval Analysis
st_duration = s_locs - q_locs;
pr_duration = diff(locs) / Fs;
analyzeSTPR(t, locs, st_duration, pr_duration, Fs);
plotECG(t, ECGsignal, locs, q_locs, s_locs, pks, 'ECG Signal with QRS Complex Annotated');
analyzeHeartRate(locs, Fs);

t_wave_abnormalities = analyzeTWave(ECGsignal, t, locs, q_locs, Fs);

p_wave_abnormalities = analyzePWave(ECGsignal, t, locs, q_locs, s_locs, Fs);

displayResults(t, locs, t_wave_abnormalities, p_wave_abnormalities);

% Analyze and plot functions
function plotECG(t, ECGsignal, locs, q_locs, s_locs, pks, titleText)
    figure; plot(t, ECGsignal, t(locs), pks, 'ro', 'MarkerFaceColor', 'r');
    hold on; plot(t(q_locs), ECGsignal(q_locs), 'ks', 'MarkerFaceColor', 'b');
    plot(t(s_locs), ECGsignal(s_locs), 'ks', 'MarkerFaceColor', 'g');

    xlabel('Time (s)'); ylabel('Amplitude'); title(titleText);
    legend('ECG Signal', 'R Peaks', 'Estimated Q Points', 'Estimated S Points'); hold off;
end

function analyzeHeartRate(locs, Fs)
    heart_rates = 60 ./ diff(locs) * Fs;
    heart_rates = heart_rates(heart_rates > 40 & heart_rates < 180);
    disp(['Average Heart Rate: ', num2str(mean(heart_rates)), ' BPM']);
    % Add more analysis if needed
end

function displayResults(t, locs, t_wave_abnormalities, p_wave_abnormalities)
    disp('T-wave Abnormalities:');
    for i = 1:length(t_wave_abnormalities)
        disp(['T-wave after QRS complex at ' num2str(t(locs(i + 1))) 's - ' t_wave_abnormalities{i}]);
    end

    disp('P-wave Abnormalities:');
    for i = 1:length(p_wave_abnormalities)
        disp(['P-wave before QRS complex at ' num2str(t(locs(i))) 's - ' p_wave_abnormalities{i}]);
    end
end

function analyzeArrhythmia(t, locs, RR_intervals)
    tachycardia_threshold = 100; % BPM
    bradycardia_threshold = 60; % BPM
    heart_rates = 60 ./ RR_intervals;

    tachycardia_indices = find(heart_rates > tachycardia_threshold);
    bradycardia_indices = find(heart_rates < bradycardia_threshold);
    rr_std = std(RR_intervals);
    rr_irregularity_threshold = 0.1; % Customize as needed

    if ~isempty(tachycardia_indices)
        disp(['Tachycardia detected at times: ', num2str(t(locs(tachycardia_indices)))]);
    else
        disp('No tachycardia detected.');
    end

    if ~isempty(bradycardia_indices)
        disp(['Bradycardia detected at times: ', num2str(t(locs(bradycardia_indices)))]);
    else
        disp('No bradycardia detected.');
    end

    if rr_std > rr_irregularity_threshold
        disp('Possible arrhythmia (atrial fibrillation) detected.');
    else
        disp('No significant arrhythmia detected.');
    end
end

function analyzeSTPR(t, locs, st_duration, pr_duration, Fs)
    pr_abnormal_threshold = 0.2; % seconds
    st_abnormal_threshold = 0.1; % seconds

    for i = 1:length(pr_duration)
        if pr_duration(i) > pr_abnormal_threshold
            disp(['At ' num2str(t(locs(i))) 's: PR interval prolonged (possible heart block)']);
        elseif pr_duration(i) < 0.12
            disp(['At ' num2str(t(locs(i))) 's: PR interval too short (possible pre-excitation syndrome)']);
        end
    end

    for i = 1:length(st_duration)
        if st_duration(i) > st_abnormal_threshold * Fs
            disp(['At ' num2str(t(locs(i))) 's: ST segment abnormal (possible myocardial ischemia)']);
        end
    end
end

function t_wave_abnormalities = analyzeTWave(ECGsignal, t, locs, q_locs, Fs)
    % Initialize cell array to store T-wave analysis results
    t_wave_abnormalities = cell(length(locs) - 1, 1);
    
    % Define parameters for T-wave morphology analysis (you can adjust these)
    t_wave_duration_threshold = 0.3; % seconds
    t_wave_amplitude_threshold = 0.1; % amplitude units (adjust as needed)
    
    for i = 2:(length(locs)-1) % Iterate up to the second-to-last QRS complex
        % Extract the T-wave segment between current QRS complex and next QRS complex
        t_wave_start_index = q_locs(i);
        t_wave_end_index = q_locs(i + 1);
        t_wave_segment = ECGsignal(t_wave_start_index:t_wave_end_index);
        
        % Measure T-wave duration and amplitude
        t_wave_duration = (t_wave_end_index - t_wave_start_index) / Fs;
        t_wave_amplitude = max(t_wave_segment) - min(t_wave_segment);
        
        % Check for T-wave inversion (example criteria)
        if t_wave_duration > t_wave_duration_threshold && t_wave_amplitude < -t_wave_amplitude_threshold
            t_wave_abnormalities{i-1} = 'T-wave inversion detected';
        else
            t_wave_abnormalities{i-1} = 'T-wave within normal limits';
        end
    end
end

function p_wave_abnormalities = analyzePWave(ECGsignal, t, locs, q_locs, s_locs, Fs)
    % Initialize cell array to store P-wave analysis results
    p_wave_abnormalities = cell(length(locs) - 1, 1);
    
    % Define parameters for P-wave morphology analysis (you can adjust these)
    p_wave_duration_threshold = 0.12; % seconds
    p_wave_amplitude_threshold = 0.1; % amplitude units (adjust as needed)
    
    for i = 2:(length(locs)-1) % Iterate up to the second-to-last QRS complex
        % Extract the P-wave segment between previous S-wave and current QRS complex
        p_wave_start_index = s_locs(i - 1);
        p_wave_end_index = q_locs(i);
        p_wave_segment = ECGsignal(p_wave_start_index:p_wave_end_index);
        
        % Measure P-wave duration and amplitude
        p_wave_duration = (p_wave_end_index - p_wave_start_index) / Fs;
        p_wave_amplitude = max(p_wave_segment) - min(p_wave_segment);
        
        % Implement your P-wave morphology analysis logic here
        % Example: Check for P-wave shape, amplitude, and duration
        
        % Interpret P-wave abnormalities based on your analysis criteria
        if p_wave_duration > p_wave_duration_threshold && p_wave_amplitude < -p_wave_amplitude_threshold
            p_wave_abnormalities{i-1} = 'P-wave inversion detected';
        else
            p_wave_abnormalities{i-1} = 'P-wave within normal limits';
        end
    end
end