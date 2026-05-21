%% FD_V.m — 线性区自动检测：找出信号进入和离开线性区的边界
clear; clc;

%% ===== 用户配置 =====
fileV = 'your_data.csv';   % 输入数据 CSV 文件（请改为实际文件名）
rowStart = 2;               % 数据起始行（跳过表头）

% ---- 线性区检测参数 ----
R2Thresh = 0.9999;            % 线性拟合 R² 阈值，越高越严格（1=完美直线）
minLinearLen = 10;           % 最小线性区长度（点数），过短不视为线性区
minVoltageRange = 0.1;       % 最小电压变化量（V），排除平坦/近零伪线性段
smoothWin = 5;               % 差值/斜率平滑窗口（奇数），仅用于展示

% ---- 绘图选项 ----
plotOverview = true;         % 全波形 + 线性区标注
plotDiff     = true;         % 差值图 + 斜率图 + 线性区标注
plotR2Diag   = true;         % R² 诊断图（滑动窗口 R² + 逐点展开 R²）

%% ===== 读取数据 =====
optsV = detectImportOptions(fileV); optsV.DataLines = [rowStart Inf];
MV = readmatrix(fileV, optsV);

t = MV(:,1);
v = MV(:,2);

good = isfinite(t) & isfinite(v);
t = t(good);
v = v(good);

[t, idx] = sort(t);
v = v(idx);
[t, iu] = unique(t, 'last');
v = v(iu);

if numel(t) < 10
  error('有效数据点太少（至少10个）');
end

%% ===== 计算前后差值 与 斜率 =====
% dv(i) = v(i+1) - v(i)，斜率 = dv/dt (V/s) — 不受时间步长变化影响
dt = diff(t);
dv = diff(v);

bad = ~isfinite(dt) | (dt <= 0) | ~isfinite(dv);
dt(bad) = NaN; dv(bad) = NaN;

t_mid = t(1:end-1) + dt/2;

% 斜率：用于线性度判断（恒定斜率 = 线性行为）
slope = dv ./ dt;

% 平滑（展示用）
dv_smooth = movmedian(dv, smoothWin, 'omitnan');
slope_smooth = movmedian(slope, smoothWin, 'omitnan');

%% ===== 检测线性区：R² 扫描 =====
% 对所有连续子区间 [i:j] 做线性拟合，找 R² ≥ R2Thresh 的最长段
% 线性信号 v = a*t + b → R² ≈ 1

n = numel(t);
fprintf('\n===== R² 线性区扫描 =====\n');
fprintf('R² 阈值: %.6f, 最小长度: %d 点\n', R2Thresh, minLinearLen);
fprintf('正在扫描 %d 个点的所有子区间...\n', n);

bestLen = 0;
bestStart = 1;
bestEnd = 0;
totalChecked = 0;
flatSkipped = 0;

for i = 1:(n - minLinearLen)
  if (n - i + 1) < bestLen
    break;  % 剩余数据长度已经无法超越最优解
  end

  for j = (i + minLinearLen - 1):n
    winLen = j - i + 1;
    if winLen <= bestLen
      continue;  % 当前窗口长度无法超越最优解
    end

    tv = t(i:j);
    vv = v(i:j);

    % 跳过电压变化过小的平坦段（近零区伪线性）
    if max(vv) - min(vv) < minVoltageRange
      flatSkipped = flatSkipped + 1;
      continue;
    end

    totalChecked = totalChecked + 1;

    coeff = polyfit(tv, vv, 1);
    v_fit = polyval(coeff, tv);
    SS_res = sum((vv - v_fit).^2);
    SS_tot = sum((vv - mean(vv)).^2);

    if SS_tot > 0
      R2 = 1 - SS_res / SS_tot;
    else
      R2 = 0;
    end

    if R2 < R2Thresh
      break;  % 继续扩展只会让 R² 更低
    end

    if winLen > bestLen
      bestLen = winLen;
      bestStart = i;
      bestEnd = j;
    end
  end
end

fprintf('扫描完成，共检查 %d 个子区间（跳过 %d 个平坦段）。\n', totalChecked, flatSkipped);

if bestLen == 0
  error('未找到 R² ≥ %.4f 的线性区。请降低 R2Thresh 或 minLinearLen。', R2Thresh);
end

linearStart = bestStart;
linearEnd   = bestEnd;
linearLen   = bestEnd - bestStart + 1;

% 对最佳段做最终拟合
tv_best = t(linearStart:linearEnd);
vv_best = v(linearStart:linearEnd);
coeff_best = polyfit(tv_best, vv_best, 1);
v_fit_best = polyval(coeff_best, tv_best);
SS_res_best = sum((vv_best - v_fit_best).^2);
SS_tot_best = sum((vv_best - mean(vv_best)).^2);
R2_best = 1 - SS_res_best / SS_tot_best;

t_linear_start = t(linearStart);
t_linear_end   = t(linearEnd);

% 线性区统计
% linearStart/linearEnd 是原始 t/v 的索引；dv/slope 在第 linearStart:linearEnd-1 个点
dvIdxStart = linearStart;
dvIdxEnd   = min(linearEnd - 1, numel(dv));
if dvIdxEnd >= dvIdxStart
  dv_linear = dv(dvIdxStart:dvIdxEnd);
  dv_linear = dv_linear(isfinite(dv_linear));
  dv_mean = mean(dv_linear);
  dv_std  = std(dv_linear);

  slope_linear = slope(dvIdxStart:dvIdxEnd);
  slope_linear = slope_linear(isfinite(slope_linear));
  slope_mean = mean(slope_linear);
  slope_std  = std(slope_linear);
else
  dv_mean = NaN; dv_std = NaN;
  slope_mean = NaN; slope_std = NaN;
end

% 拟合直线参数
a = coeff_best(1);  % 斜率 (V/s)
b = coeff_best(2);  % 截距 (V)

%% ===== 命令行输出 =====
fprintf('\n========== 线性区检测结果 (R² 法) ==========\n');
fprintf('数据文件          : %s\n', fileV);
fprintf('总数据点数        : %d\n', n);
fprintf('R² 阈值           : %.6f\n', R2Thresh);
fprintf('线性区 R²         : %.6f\n', R2_best);
fprintf('拟合直线          : v = %.6e * t + %.6e\n', a, b);
fprintf('--------------------------------------\n');
fprintf('线性区起始        : idx=%d, t=%.6e s (%.3f ns), v=%.6f V\n', ...
  linearStart, t_linear_start, t_linear_start*1e9, v(linearStart));
fprintf('线性区结束        : idx=%d, t=%.6e s (%.3f ns), v=%.6f V\n', ...
  linearEnd, t_linear_end, t_linear_end*1e9, v(linearEnd));
fprintf('线性区长度        : %d 点\n', linearLen);
fprintf('线性区时间跨度    : %.6e s (%.3f ns)\n', ...
  t_linear_end - t_linear_start, (t_linear_end - t_linear_start)*1e9);
fprintf('电压变化          : %.6f V\n', v(linearEnd) - v(linearStart));
fprintf('--------------------------------------\n');
if ~isnan(slope_mean)
  fprintf('线性区内斜率均值  : %.6e V/s\n', slope_mean);
  fprintf('线性区内斜率标准差: %.6e V/s\n', slope_std);
  if slope_mean ~= 0
    fprintf('线性区内斜率波动  : ±%.2f%%\n', (slope_std / abs(slope_mean)) * 100);
  end
  fprintf('线性区内差值均值  : %.6e V/sample\n', dv_mean);
  fprintf('线性区内差值标准差: %.6e V/sample\n', dv_std);
end
fprintf('======================================\n\n');

%% ===== Figure 1：全波形 + 线性区标注 =====
if plotOverview
  figure('Color', 'w', 'Name', '全波形与线性区');
  plot(t*1e9, v*1e3, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8); hold on;
  plot(t(linearStart:linearEnd)*1e9, v(linearStart:linearEnd)*1e3, 'b-', 'LineWidth', 1.5);
  plot(t(linearStart:linearEnd)*1e9, v_fit_best*1e3, 'r--', 'LineWidth', 1.2);
  xline(t_linear_start*1e9, '--g', 'LineWidth', 1.2);
  xline(t_linear_end*1e9, '--r', 'LineWidth', 1.2);
  text(t_linear_start*1e9, max(v)*1e3*0.95, '进入线性区', 'Color', 'g', 'FontSize', 10);
  text(t_linear_end*1e9,   max(v)*1e3*0.85, '离开线性区', 'Color', 'r', 'FontSize', 10);

  legend({'原始数据', '线性区', '拟合直线', '起点', '终点'}, 'Location', 'best');
  xlabel('Time (ns)'); ylabel('Voltage (mV)');
  title(sprintf('全波形 — R²=%.6f, 线性区 idx=%d~%d', R2_best, linearStart, linearEnd));
  grid on; hold off;
end

%% ===== Figure 2：差值图 & 斜率图 + 线性区标注 =====
if plotDiff && dvIdxEnd >= dvIdxStart
  figure('Color', 'w', 'Name', '差值与斜率分析');

  subplot(2,2,1);
  plot(t_mid*1e9, dv, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on;
  plot(t_mid(dvIdxStart:dvIdxEnd)*1e9, dv(dvIdxStart:dvIdxEnd), 'b-', 'LineWidth', 1.2);
  xline(t_linear_start*1e9, '--g', 'LineWidth', 1.2);
  xline(t_linear_end*1e9, '--r', 'LineWidth', 1.2);
  xlabel('Time (ns)'); ylabel('dV (V/sample)');
  title(sprintf('前后差值 dv(i)=v(i+1)-v(i)'));
  legend({'dv原始', '线性区'}, 'Location', 'best');
  grid on; hold off;

  subplot(2,2,2);
  plot(t_mid*1e9, slope, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on;
  plot(t_mid(dvIdxStart:dvIdxEnd)*1e9, slope(dvIdxStart:dvIdxEnd), 'b-', 'LineWidth', 1.2);
  xline(t_linear_start*1e9, '--g', 'LineWidth', 1.2);
  xline(t_linear_end*1e9, '--r', 'LineWidth', 1.2);
  if ~isnan(slope_mean)
    yline(slope_mean, '--m', 'LineWidth', 1);
    yline(slope_mean + slope_std, ':m');
    yline(slope_mean - slope_std, ':m');
  end
  xlabel('Time (ns)'); ylabel('Slope (V/s)');
  title(sprintf('斜率 dv/dt'));
  legend({'slope原始', '线性区'}, 'Location', 'best');
  grid on; hold off;

  subplot(2,2,3);
  plot(t_mid*1e9, dv_smooth, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on;
  plot(t_mid(dvIdxStart:dvIdxEnd)*1e9, dv_smooth(dvIdxStart:dvIdxEnd), 'b-', 'LineWidth', 1.2);
  xline(t_linear_start*1e9, '--g', 'LineWidth', 1.2);
  xline(t_linear_end*1e9, '--r', 'LineWidth', 1.2);
  xlabel('Time (ns)'); ylabel('dV smooth (V/sample)');
  title(sprintf('平滑差值 (movmedian win=%d)', smoothWin));
  legend({'dv平滑', '线性区'}, 'Location', 'best');
  grid on; hold off;

  subplot(2,2,4);
  plot(t_mid*1e9, slope_smooth, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); hold on;
  plot(t_mid(dvIdxStart:dvIdxEnd)*1e9, slope_smooth(dvIdxStart:dvIdxEnd), 'b-', 'LineWidth', 1.2);
  xline(t_linear_start*1e9, '--g', 'LineWidth', 1.2);
  xline(t_linear_end*1e9, '--r', 'LineWidth', 1.2);
  xlabel('Time (ns)'); ylabel('Slope smooth (V/s)');
  title(sprintf('平滑斜率 (movmedian win=%d)', smoothWin));
  legend({'slope平滑', '线性区'}, 'Location', 'best');
  grid on; hold off;
end

%% ===== Figure 3：R² 诊断图 =====
if plotR2Diag
  figure('Color', 'w', 'Name', 'R² 线性度诊断');

  % 子图1：滑动窗口 R²（窗口大小 = minLinearLen*2）
  subplot(2,1,1);
  winR2 = minLinearLen * 2;
  if winR2 < 5, winR2 = 5; end
  if winR2 > n, winR2 = n; end

  r2_sliding = nan(n, 1);
  t_sliding  = nan(n, 1);
  for k = 1:(n - winR2 + 1)
    tv = t(k:k+winR2-1);
    vv = v(k:k+winR2-1);
    coeff = polyfit(tv, vv, 1);
    v_fit = polyval(coeff, tv);
    SS_res = sum((vv - v_fit).^2);
    SS_tot = sum((vv - mean(vv)).^2);
    if SS_tot > 0
      r2_sliding(k) = 1 - SS_res / SS_tot;
    end
    t_sliding(k) = t(k) + (t(k+winR2-1) - t(k))/2;
  end

  plot(t_sliding*1e9, r2_sliding, 'b-', 'LineWidth', 1); hold on;
  yline(R2Thresh, '--r', 'LineWidth', 1.2);
  xline(t_linear_start*1e9, '--g', 'LineWidth', 1.2);
  xline(t_linear_end*1e9, '--r', 'LineWidth', 1.2);
  xlabel('Time (ns)'); ylabel('Local R²');
  title(sprintf('滑动窗口 R² (窗口=%d点) — 阈值=%.4f', winR2, R2Thresh));
  legend({'局部R²', '阈值', '线性区起点', '线性区终点'}, 'Location', 'best');
  grid on; hold off;

  % 子图2：从线性区起点开始的逐点扩展 R²
  subplot(2,1,2);
  r2_expand = nan(n - linearStart + 1, 1);
  t_expand  = nan(n - linearStart + 1, 1);
  for k = linearStart:n
    tv = t(linearStart:k);
    vv = v(linearStart:k);
    if numel(tv) < 3, continue; end
    coeff = polyfit(tv, vv, 1);
    v_fit = polyval(coeff, tv);
    SS_res = sum((vv - v_fit).^2);
    SS_tot = sum((vv - mean(vv)).^2);
    if SS_tot > 0
      r2_expand(k - linearStart + 1) = 1 - SS_res / SS_tot;
    end
    t_expand(k - linearStart + 1) = t(k);
  end

  plot(t_expand*1e9, r2_expand, 'b-', 'LineWidth', 1); hold on;
  yline(R2Thresh, '--r', 'LineWidth', 1.2);
  xline(t_linear_end*1e9, '--r', 'LineWidth', 1.2);
  xlabel('Time (ns)'); ylabel('Cumulative R²');
  title(sprintf('从起点(idx=%d)逐点展开的 R² — 终点 idx=%d', linearStart, linearEnd));
  legend({'累积R²', '阈值', '检测终点'}, 'Location', 'best');
  grid on; hold off;
end
