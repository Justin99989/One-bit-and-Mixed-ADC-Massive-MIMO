clear; close all; clc;
rng(1); % 設定隨機種子

% --- 模擬參數設定 ---
K = 10;           % 使用者數量
T = 200;          % 相干時間長度
% eta = K;        % eta 不再固定，將被優化
snr_dB = -20:1:20; % SNR 範圍
sigma2 = 1;       % 雜訊方差
trials = 1000;  % 暫時專注於理論優化，如果需要模擬，再取消註解

% --- 論文 Figure 1 中的不同配置 ---
mu_list = [0, 18/95, 25/55, 29/31];
M_list = [200, 95, 55, 31];
mu_frac_str = {'0', '18/95', '25/55', '29/31'};

% --- 結果儲存變數初始化 ---
sumSE_theory_optimized_eta = zeros(length(mu_list), length(snr_dB));
optimal_eta_values = zeros(length(mu_list), length(snr_dB)); % 儲存每個點的最佳eta

fprintf('開始進行模擬 (包含 eta 優化)...\n');

% --- eta 優化參數 ---
eta_min_factor = 1.0; % eta_min = K * eta_min_factor (至少為K)
eta_max_factor = 0.5; % eta_max = T * eta_max_factor (例如，不超過 T 的一半)

% --- 迴圈遍歷不同的 mu/M 配置 ---
for m_idx = 1:length(mu_list)
    mu = mu_list(m_idx);
    M = M_list(m_idx);
    M0 = round(mu * M);
    M1 = M - M0;

    fprintf('  正在執行配置: mu = %s, M = %d (M0=%d, M1=%d)\n', mu_frac_str{m_idx}, M, M0, M1);

    % --- 迴圈遍歷不同的 SNR 值 ---
    for s_idx = 1:length(snr_dB)
        SNR = snr_dB(s_idx);
        p = 10^(SNR/10);

        % --- 內部迴圈：優化 eta ---
        eta_search_min = ceil(K * eta_min_factor);
        eta_search_max = floor(T * eta_max_factor);
        if eta_search_max <= eta_search_min
            % 如果 T 太小，或者 K 太大，可能導致搜索範圍無效
            % 這裡可以選擇一個固定的 eta，例如 K，或者報錯
            warning('eta 搜索範圍無效 (T=%d, K=%d)。對此配置和SNR，將使用 eta = K。', T, K);
            current_optimal_eta = K;
            eta_range_current_optim = K;
        else
            eta_range_current_optim = eta_search_min:eta_search_max;
        end

        sumSE_for_eta_optimization = zeros(1, length(eta_range_current_optim));

        for eta_idx = 1:length(eta_range_current_optim)
            eta_current = eta_range_current_optim(eta_idx);

            if eta_current <= 0 || eta_current >= T % 確保 eta 在有效範圍內
                sumSE_for_eta_optimization(eta_idx) = -inf; % 無效的 eta 給予極低 SE
                continue;
            end

            % 理論總頻譜效率計算 (使用 eta_current)
            % sigma_0_sq 和 sigma_1_sq 依賴於 eta_current
            % 假設 p*K 中的 K 代表導頻資源量，現在是 eta_current
            effective_pilot_resource = p * eta_current; % 每個用戶的有效導頻能量/資源
            if effective_pilot_resource <= 0 % 避免 log(0) 或除以零
                 sigma_0_sq_eta = 0;
                 sigma_1_sq_eta = 0;
            else
                sigma_0_sq_eta = effective_pilot_resource / (effective_pilot_resource + sigma2);
                sigma_1_sq_eta = (2/pi) * effective_pilot_resource / (effective_pilot_resource + sigma2);
            end


            numerator_th = M * (mu * sigma_0_sq_eta + (1-mu) * sigma_1_sq_eta)^2;
            denominator_part1 = K + sigma2/p; % K 是用戶數
            denominator_part2 = mu * sigma_0_sq_eta + (1-mu) * (pi/2) * sigma_1_sq_eta;
            denominator_th = denominator_part1 * denominator_part2;

            if denominator_th > eps
                SINR_th_eta = numerator_th / denominator_th;
            else
                SINR_th_eta = 0;
            end

            % SE 計算，考慮訓練開銷 (1 - eta_current/T)
            SE_each_th_eta = (1 - eta_current/T) * log2(1 + max(SINR_th_eta, 0));
            sumSE_for_eta_optimization(eta_idx) = K * SE_each_th_eta;
        end

        % 找到當前 mu, M, SNR 配置下的最佳 eta
        [max_SE_current_config, best_eta_idx] = max(sumSE_for_eta_optimization);
        if isempty(max_SE_current_config) || max_SE_current_config == -inf
            % 如果所有 eta 都無效 (例如 T 太小)
            sumSE_theory_optimized_eta(m_idx, s_idx) = 0;
            optimal_eta_values(m_idx, s_idx) = NaN; % 標記為無有效 eta
        else
            sumSE_theory_optimized_eta(m_idx, s_idx) = max_SE_current_config;
            optimal_eta_values(m_idx, s_idx) = eta_range_current_optim(best_eta_idx);
        end

        if mod(s_idx, round(length(snr_dB)/4)) == 0 || s_idx == length(snr_dB)
            fprintf('    SNR = %d dB 完成. 優化得到的 eta = %d, Max SE = %.3f\n', ...
                SNR, optimal_eta_values(m_idx, s_idx), sumSE_theory_optimized_eta(m_idx, s_idx));
        end
    end % 結束 SNR 迴圈
end % 結束 mu/M 配置迴圈

fprintf('模擬完成。\n正在繪製結果...\n');

% --- 繪製結果圖 (優化 eta 後的理論頻譜效率) ---
figure;
hold on;
colors = {'k', 'b', 'r', 'g'};
lineStyles = {'-', '--', '-.', ':'};
% markers = {'o', 's', '^', 'd'}; % 理論曲線通常不加過多marker

legend_entries = {};

for m_idx = 1:length(mu_list)
    M_val = M_list(m_idx);
    frac_str = mu_frac_str{m_idx};
    plot(snr_dB, sumSE_theory_optimized_eta(m_idx,:), ...
        'Color', colors{m_idx}, ...
        'LineStyle', lineStyles{m_idx}, ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\mu = %s, M=%d (Optimized \\eta)', frac_str, M_val));
end

xlabel('SNR (dB)');
ylabel('Sum Spectral Efficiency (Optimized \eta) [bits/s/Hz]');
title('Mixed ADC Massive MIMO Performance with Optimized Pilot Length');
legend('show', 'Location', 'southeast');
grid on;
axis([min(snr_dB) max(snr_dB) 0 inf]);
hold off;

% --- (可選) 繪製最佳 eta 值隨 SNR 的變化 ---
figure;
hold on;
for m_idx = 1:length(mu_list)
    frac_str = mu_frac_str{m_idx};
    M_val = M_list(m_idx);
    plot(snr_dB, optimal_eta_values(m_idx,:), ...
        'Color', colors{m_idx}, ...
        'LineStyle', lineStyles{m_idx}, ...
        'Marker', '.', ... % 加小點標記
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\mu = %s, M=%d', frac_str, M_val));
end
xlabel('SNR (dB)');
ylabel('Optimal Pilot Length \eta^*');
title('Optimal Pilot Length vs. SNR for Different Configurations');
legend('show', 'Location', 'northwest');
grid on;
hold off;