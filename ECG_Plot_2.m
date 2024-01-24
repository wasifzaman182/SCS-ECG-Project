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

% ST Segment and PR Interval Analysis
st_duration = s_locs - q_locs;
pr_duration = diff(locs) / Fs;
plotECG(t, ECGsignal, locs, q_locs, s_locs, pks, 'ECG Signal with QRS Complex Annotated',Fs);
disp(['Average Heart Rate: ', num2str(mean(heart_rates)), ' BPM']);
analyzeHeartRate(locs, Fs);
analyzeArrhythmia(t, locs, RR_intervals);
analyzeSTPR(t, locs, st_duration, pr_duration, Fs);

t_wave_abnormalities = analyzeTWave(ECGsignal, t, locs, q_locs, Fs);

p_wave_abnormalities = analyzePWave(ECGsignal, t, locs, q_locs, s_locs, Fs);

displayResults(t, locs, t_wave_abnormalities, p_wave_abnormalities);

% Analyze and plot functions
function plotECG(t, ECGsignal, locs, q_locs, s_locs, pks, titleText, Fs)
    figure;
    % Graph 1: ECG Waveform with QRS Complex Annotation
    subplot(2,1,1);
    xlabel('Time (s)'); ylabel('Amplitude'); title(titleText);
    hold on;

    % Graph 2: Cardiac Analysis Data
    subplot(2,1,2);
    xlabel('Time (s)'); ylabel('Metrics'); title('Cardiac Analysis Data');
    hold on;

    samplePeriod = 1 / Fs; % Time period of one sample
    bufferSize = 10; % Number of samples to buffer
    bufferTime = samplePeriod * bufferSize; % Total time of buffered samples

    accumulatedLocs = []; % Initialize outside the loop to accumulate R peaks

    for i = 1:bufferSize:length(ECGsignal)
        bufferEnd = min(i+bufferSize-1, length(ECGsignal));

        % Update Graph 1: ECG waveform
        subplot(2,1,1);
        plot(t(i:bufferEnd), ECGsignal(i:bufferEnd), 'b'); % Plot buffered ECG signal

        % Plotting R, Q, and S points within the buffer
        currentLocs = locs(locs >= i & locs <= bufferEnd);
        accumulatedLocs = [accumulatedLocs, currentLocs]; % Append new R peaks
        currentQ_locs = q_locs(q_locs >= i & q_locs <= bufferEnd);
        currentS_locs = s_locs(s_locs >= i & s_locs <= bufferEnd);

        if ~isempty(currentLocs)
            plot(t(currentLocs), pks(ismember(locs, currentLocs)), 'ro', 'MarkerFaceColor', 'r');
        end

        if ~isempty(currentQ_locs)
            plot(t(currentQ_locs), ECGsignal(currentQ_locs), 'ks', 'MarkerFaceColor', 'b');
        end

        if ~isempty(currentS_locs)
            plot(t(currentS_locs), ECGsignal(currentS_locs), 'ks', 'MarkerFaceColor', 'g');
        end

        % Update Graph 2: Analysis Data (e.g., Heart Rate)
        subplot(2,1,2);
        if ~isempty(accumulatedLocs)
           [current_heart_rate] =analyzeHeartRate(accumulatedLocs, Fs);
            %disp(['print heart beat ', num2str(current_heart_rate)]);
            analyzeArrhythmia(t, accumulatedLocs, diff(accumulatedLocs)/Fs);
            analyzeSTPR(t, accumulatedLocs, s_locs - q_locs, diff(accumulatedLocs)/Fs, Fs);
        end

        pause(bufferTime); % Pause for the duration of the buffered time
    end

    % Finalize Graph 1
    subplot(2,1,1);
    legend('ECG Signal', 'R Peaks', 'Estimated Q Points', 'Estimated S Points');

    % Finalize Graph 2
    subplot(2,1,2);
    legend('Analysis Data'); % Update with appropriate legend items
    hold off;
end



function current_heart_rate=analyzeHeartRate(locs, Fs)
    current_heart_rate = NaN;
    if length(locs) >= 2
        RR_intervals = diff(locs) / Fs;
        heart_rates = 60 ./ RR_intervals;
        current_heart_rate = heart_rates(end); % Latest heart rate based on the most recent RR interval
        times = locs(2:end) / Fs; % Times of heart rates (excluding the first R peak)
        subplot(2,1,2); % Select the second subplot for plotting
        plot(times, heart_rates, 'm-o'); % Plot heart rate over time
        ylabel('Heart Rate (BPM)');
        disp(['Current Heart Rate: ', num2str(current_heart_rate), ' BPM']);
    elseif length(locs) == 1
        % Only one R peak detected, insufficient for heart rate calculation
        disp('Awaiting additional data for heart rate calculation...');
    else
        % No R peaks detected yet
        disp('No heart rate data available yet.');
    end
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

function [tachycardia_indices, bradycardia_indices, arrhythmia_flag] = analyzeArrhythmia(t, locs, RR_intervals)
    tachycardia_threshold = 100; % BPM
    bradycardia_threshold = 60; % BPM
    heart_rates = 60 ./ RR_intervals;
    arrhythmia_flag = false;
    
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
        arrhythmia_flag=true;
    else
        disp('No significant arrhythmia detected.');
    end
end

function [prolonged_pr_indices, short_pr_indices, abnormal_st_indices] = analyzeSTPR(t, locs, st_duration, pr_duration, Fs)
    pr_abnormal_threshold = 0.2; % seconds
    st_abnormal_threshold = 0.1; % seconds

    prolonged_pr_indices = []; % Indices where PR interval is prolonged
    short_pr_indices = []; % Indices where PR interval is short
    abnormal_st_indices = []; % Indices where ST segment is abnormal

    for i = 1:length(pr_duration)
        if pr_duration(i) > pr_abnormal_threshold
            disp(['At ' num2str(t(locs(i))) 's: PR interval prolonged (possible heart block)']);
            prolonged_pr_indices = [prolonged_pr_indices, i];
        elseif pr_duration(i) < 0.12
            disp(['At ' num2str(t(locs(i))) 's: PR interval too short (possible pre-excitation syndrome)']);
            short_pr_indices = [short_pr_indices, i];
        end
    end

    for i = 1:length(st_duration)
        if st_duration(i) > st_abnormal_threshold * Fs
            disp(['At ' num2str(t(locs(i))) 's: ST segment abnormal (possible myocardial ischemia)']);
            abnormal_st_indices = [abnormal_st_indices, i];
        end
    end
end

function [t_wave_abnormalities, t_wave_abnormal_indices] = analyzeTWave(ECGsignal, t, locs, q_locs, Fs)
    % Initialize cell array to store T-wave analysis results and indices
    t_wave_abnormalities = cell(length(locs) - 1, 1);
    t_wave_abnormal_indices = [];

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
        
        % Check for T-wave inversion
        if t_wave_duration > t_wave_duration_threshold && t_wave_amplitude < -t_wave_amplitude_threshold
            t_wave_abnormalities{i-1} = 'T-wave inversion detected';
            t_wave_abnormal_indices = [t_wave_abnormal_indices, i-1];
        else
            t_wave_abnormalities{i-1} = 'T-wave within normal limits';
        end
    end

    % The function returns both the analysis results and indices of abnormalities
end


function [p_wave_abnormalities, p_wave_abnormal_indices] = analyzePWave(ECGsignal, t, locs, q_locs, s_locs, Fs)
    % Initialize cell array to store P-wave analysis results
    p_wave_abnormalities = cell(length(locs) - 1, 1);
    p_wave_abnormal_indices = [];

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
            p_wave_abnormal_indices = [p_wave_abnormal_indices, i-1];
        else
            p_wave_abnormalities{i-1} = 'P-wave within normal limits';
        end
    end
end