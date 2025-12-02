clear; close all; clc;
rng(1); % 設定隨機種子，確保每次執行模擬結果一樣，方便除錯和比較

% --- 模擬參數設定 ---
K = 10;           % 使用者數量
T = 200;          % 相干時間長度 (一個通道保持不變的時間區塊大小，以符號數為單位)
eta = K;          % 導頻序列長度 (這裡設定導頻長度等於使用者數量，論文 Figure 1 的設定)
snr_dB = -20:1:20; % 訊號雜訊比 (SNR) 範圍，單位 dB
sigma2 = 1;       % 雜訊方差 (這裡假設雜訊功率為 1)
trials = 1000;    % 蒙地卡羅模擬的次數，越多結果越平滑，但計算時間越長

% --- 論文 Figure 1 中的不同配置 (mu 和 M) ---
% mu_list: 高解析度 ADC 天線佔總天線數的比例
% M_list: 對應的總天線數量 (這些 M 值是根據論文中的功耗預算和 mu 計算出來的)
mu_list = [0, 18/95, 25/55, 29/31];
M_list = [200, 95, 55, 31];
% 用於圖例顯示的 mu 分數格式字串
mu_frac_str = {'0', '18/95', '25/55', '29/31'};

% --- 結果儲存變數初始化 ---
% sumSE_theory: 儲存理論計算的總頻譜效率
% sumSE_sim: 儲存蒙地卡羅模擬得到的總頻譜效率
sumSE_theory = zeros(length(mu_list), length(snr_dB));
sumSE_sim = zeros(length(mu_list), length(snr_dB));

fprintf('開始進行模擬...\n');

% --- 迴圈遍歷不同的 mu/M 配置 ---
for m_idx = 1:length(mu_list)
    mu = mu_list(m_idx);
    M = M_list(m_idx);
    M0 = round(mu * M);   % 高解析度 ADC 的天線數量 (取整數)
    M1 = M - M0;          % 一比特 ADC 的天線數量

    fprintf('  正在執行配置: mu = %.3f, M = %d (M0=%d, M1=%d)\n', mu, M, M0, M1);

    % --- 迴圈遍歷不同的 SNR 值 ---
    for s_idx = 1:length(snr_dB)
        % 計算當前 SNR 對應的線性值和發射功率 p
        SNR = snr_dB(s_idx);
        p = 10^(SNR/10);  % 發射功率 (這裡假設每個使用者發射功率相同)

        % --- 理論總頻譜效率計算 (對應論文公式 16) ---
        % 論文公式 9: 通道估計的方差
        % sigma_0_sq: 高解析度 ADC 部分通道估計的方差
        sigma_0_sq = p * K / (p * K + sigma2);
        % sigma_1_sq: 一比特 ADC 部分通道估計的方差 (注意這裡有個 2/pi 的因子)
        sigma_1_sq = (2/pi) * p * K / (p * K + sigma2);

        % 論文公式 16 中的 SINR 項 (括號裡面的分數部分)
        % 分子部分 (對應論文公式 16 分子，並乘以 M，因為論文的 gamma 是針對單個使用者，總 SE 是 K*SE_k)
        numerator_th = M * (mu * sigma_0_sq + (1-mu) * sigma_1_sq)^2;
        % 分母部分 (對應論文公式 16 分母)
        denominator_part1 = K + sigma2/p;
        denominator_part2 = mu * sigma_0_sq + (1-mu) * (pi/2) * sigma_1_sq; % 注意這裡一比特部分有個 pi/2 的因子
        denominator_th = denominator_part1 * denominator_part2;

        % 計算理論 SINR (論文公式 16 括號裡面的部分)
        if denominator_th > eps % 避免除以零或非常小的數
            SINR_th = numerator_th / denominator_th;
        else
            SINR_th = 0;
        end

        % 計算每個使用者的理論頻譜效率 (對應論文公式 14 的 R(theta))
        % (1 - eta/T) 是前置因子，考慮了訓練時間的開銷
        SE_each_th = (1 - eta/T) * log2(1 + max(SINR_th, 0)); % log2(1+SINR)
        % 計算總頻譜效率 (K 個使用者總和)
        sumSE_theory(m_idx, s_idx) = K * SE_each_th;

        % --- 蒙地卡羅模擬計算總頻譜效率 ---
        % 這裡我們模擬系統的實際收發過程，然後統計計算 SINR
        cross_sum = zeros(K,1); % 用於累積 E[x_k^* * x_hat_k] 的結果
        pow_sum = zeros(K,1);   % 用於累積 E[|x_hat_k|^2] 的結果

        % 執行多次模擬試驗
        for t = 1:trials
            % 生成通道矩陣 H (假設 Rayleigh 衰落，獨立同分佈)
            H = (randn(M,K) + 1j*randn(M,K))/sqrt(2);
            H0 = H(1:M0, :);      % 高解析度 ADC 對應的通道部分
            H1 = H(M0+1:end, :);  % 一比特 ADC 對應的通道部分

            % --- 模擬訓練階段 ---
            % 生成訓練階段的雜訊
            n0_p = (randn(M0,K) + 1j*randn(M0,K))/sqrt(2) * sqrt(sigma2);
            n1_p = (randn(M1,K) + 1j*randn(M1,K))/sqrt(2) * sqrt(sigma2);
            % 生成訓練階段的接收訊號 (高解析度部分)
            Y0 = sqrt(p) * H0 + n0_p;
            % 生成訓練階段的接收訊號 (一比特部分，模擬類比訊號)
            Y1_analog = sqrt(p) * H1 + n1_p;
            % 對一比特部分的接收訊號進行量化 (取實部和虛部的符號)
            Y1_sign = (sign(real(Y1_analog)) + 1j*sign(imag(Y1_analog))) / sqrt(2); % 量化到 1/sqrt(2) * (±1 ±j)

            % --- 通道估計 ---
            % 這裡使用簡化的通道估計方法 (類似 LS 估計，並對一比特部分應用 Bussgang 增益)
            H0_est = Y0 / sqrt(p); % 高解析度部分直接除以 sqrt(p) (簡化處理)
            % 計算 Bussgang 增益 Ad (這裡使用資料階段的 Ad 形式，用導頻 SNR)
            effective_pilot_snr = p / sigma2;
            Ad = sqrt(2/pi) / sqrt(1 + 1/effective_pilot_snr); % 對應論文公式 10 下方的 Ad 定義
            H1_est = Ad * Y1_sign; % 一比特部分應用 Bussgang 增益 (簡化處理)
            % 注意：論文公式 8 使用的是更精確的 LMMSE 估計，這裡的模擬是簡化版本。

            % --- 構建 MRC 組合矩陣 W (對應論文公式 11) ---
            % 使用估計的通道構建 MRC 矩陣
            % 注意：論文公式 11 的 W 包含 A_d^{-1}，這裡的模擬 W 沒有包含，
            % 但後面的 SINR 計算方法會自動考慮這個影響，因為我們是直接計算 E[x_hat x^*] 和 E[|x_hat|^2]。
            W = [H0_est; H1_est]'; % 將高解析度估計和一比特估計合併並轉置

            % --- 模擬資料傳輸階段 ---
            % 生成待傳輸的資料符號 x (假設是 CN(0, I_K))
            x = (randn(K,1) + 1j*randn(K,1))/sqrt(2);
            % 生成資料階段的雜訊
            n0_d = (randn(M0,1) + 1j*randn(M0,1))/sqrt(2) * sqrt(sigma2);
            n1_d = (randn(M1,1) + 1j*randn(M1,1))/sqrt(2) * sqrt(sigma2);
            % 生成資料階段的接收訊號 (高解析度部分)
            y0 = sqrt(p) * H0 * x + n0_d;
            % 生成資料階段的接收訊號 (一比特部分，模擬類比訊號)
            y1_analog = sqrt(p) * H1 * x + n1_d;
            % 對一比特部分的接收訊號進行量化
            y1 = (sign(real(y1_analog)) + 1j*sign(imag(y1_analog))) / sqrt(2);
            % 合併接收訊號
            y = [y0; y1];
            % 應用 MRC 組合矩陣得到估計的資料符號
            x_hat = W * y; % 對應論文公式 12 的計算

            % --- 累積統計量以估計 SQINR ---
            % 累積 E[x_k^* * x_hat_k] 的分子部分 (用於計算信號功率)
            cross_sum = cross_sum + conj(x) .* x_hat;
            % 累積 E[|x_hat_k|^2] (用於計算總功率)
            pow_sum = pow_sum + abs(x_hat).^2;
        end % 結束蒙地卡羅試驗迴圈

        % --- 根據累積的統計量計算模擬的 SQINR 和總頻譜效率 ---
        cross_mean = cross_sum / trials; % 估計 E[x_k^* * x_hat_k]
        pow_mean = pow_sum / trials;   % 估計 E[|x_hat_k|^2]

        % 計算信號功率 (對應論文公式 15 分子中的 E[|h_hat_k^H h_k|^2] 項，但這裡直接從 x_hat 和 x 估計)
        signal_power = abs(cross_mean).^2;
        % 計算干擾加雜訊功率 (總功率減去信號功率)
        interf_noise_power = pow_mean - signal_power;
        interf_noise_power(interf_noise_power <= 0) = eps; % 防止出現非正數或零

        % 計算每個使用者的模擬 SINR (對應論文公式 15)
        SINR_sim_users = signal_power ./ interf_noise_power;
        % 計算每個使用者的模擬頻譜效率
        SE_sim_users = (1 - eta/T) * log2(1 + max(SINR_sim_users, 0));
        % 計算總頻譜效率 (所有使用者求和)
        sumSE_sim(m_idx, s_idx) = sum(SE_sim_users);

        % 顯示進度
        if mod(s_idx, 5) == 0 || s_idx == length(snr_dB)
            fprintf('    SNR = %d dB 完成。\n', SNR);
        end
    end % 結束 SNR 迴圈
end % 結束 mu/M 配置迴圈

fprintf('模擬完成。\n正在繪製結果...\n');

% --- 繪製結果圖 ---
figure;
hold on; % 允許在同一張圖上繪製多條曲線
colors = {'k', 'b', 'r', 'g'}; % 不同配置的顏色
lineStyles = {'-', '--', '-.', ':'}; % 不同配置的線條樣式
markers = {'o', 's', '^', 'd'}; % 不同配置的模擬點標記
legend_entries = {}; % 用於儲存圖例條目 (這裡使用了 DisplayName，這個變數其實可以不用)

% 迴圈繪製每種配置的理論曲線和模擬點
for m_idx = 1:length(mu_list)
    M_val = M_list(m_idx);
    mu_val = mu_list(m_idx);
    frac_str = mu_frac_str{m_idx};
    % Plot the theoretical curve with markers at intervals
    marker_interval = 5; % Show markers every 5 points
    plot(snr_dB, sumSE_theory(m_idx,:), ...
        'Color', colors{m_idx}, ...
        'LineStyle', lineStyles{m_idx}, ...
        'Marker', markers{m_idx}, ... % Add markers directly in the plot command
        'MarkerIndices', 1:marker_interval:length(snr_dB), ... % Specify marker positions
        'MarkerFaceColor', colors{m_idx}, ...
        'MarkerSize', 7, ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\mu = %s, M=%d', frac_str, M_val)); % Use DisplayName for legend
end

% 設定圖形標籤和標題
xlabel('SNR (dB)');
ylabel('sum spectral efficiency [bits/s/Hz]');
title('Mixed ADC Massive MIMO Performance');
legend('show', 'Location', 'southeast'); % 顯示圖例，位置在右下角
grid on; % 顯示網格
axis([min(snr_dB) max(snr_dB) 0 inf]); % 設定座標軸範圍
hold off; % 結束繪圖

