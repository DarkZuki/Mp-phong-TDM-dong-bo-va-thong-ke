%% MÔ PHỎNG TDM ĐỒNG BỘ (STDM) VÀ TDM THỐNG KÊ (StatTDM) - 4 KÊNH
% Bài tập lớn môn Cơ sở thông tin số. Lưu file UTF-8; MATLAB R2016b trở lên, không cần toolbox.
% Chạy sẽ in toàn bộ số liệu ra Command Window và mở 3 cửa sổ hình:
%   Cửa sổ 1: so sánh STDM và StatTDM | Cửa sổ 2: riêng StatTDM | Cửa sổ 3: riêng STDM
% Nguyên tắc so sánh công bằng: cùng dữ liệu nguồn, mặt nạ hoạt động, khung 40 bit và cùng
% một mặt nạ lỗi BSC cho hai luồng. Giả định: DEMUX StatTDM biết trước active_mask và ranh giới
% khung; DEMUX STDM tự tìm FAS. Dữ liệu là bit ngẫu nhiên, không phải thoại G.711.

clear; clc; close all;

%% 1. THAM SỐ HỆ THỐNG (Bảng 3.1 báo cáo)
num_channels = 4;  bits_per_sample = 8;  frame_rate = 8000;      % khung/s
frame_period_us = 1e6/frame_rate;                                % 125 us
num_frames = 1000;  bits_per_frame = 40;                         % 1.000 khung = 40.000 bit
line_rate = bits_per_frame*frame_rate;                           % 320 kbps
fas = [1 0 0 1 1 0 1 1];  fas_len = numel(fas);                  % FAS 10011011
id_bits = 2;  rec_len = id_bits + bits_per_sample;               % Channel ID 2 bit, bản ghi 10 bit
ber_levels = [0 1e-4 1e-3 1e-2];  n_bers = numel(ber_levels);
n_trials = 50;  n_miss = 3;  seed_main = 42;                     % flywheel: 3 FAS sai liên tiếp mới mất ĐB
total_bits = num_frames*bits_per_frame;

scenario_names = {'A - 100% (4 kênh)','B - 75% (3 kênh)','C - 50% (2 kênh)','D - 25% (1 kênh)'};
load_labels    = {'100%','75%','50%','25%'};
load_names     = {'Tải 100%','Tải 75%','Tải 50%','Tải 25%'};
active_masks   = logical([1 1 1 1; 1 1 1 0; 1 1 0 0; 1 0 0 0]);
n_sc           = size(active_masks,1);
n_active_vec   = sum(active_masks,2);
method_names   = {'STDM - FSM cơ bản','STDM - FSM flywheel','StatTDM'};
P = struct('fas',fas,'fl',fas_len,'nf',num_frames,'bpf',bits_per_frame, ...
           'nch',num_channels,'bps',bits_per_sample,'nm',n_miss);     % gói tham số cho các hàm DEMUX

assert(bits_per_frame == fas_len + num_channels*bits_per_sample, 'Khung STDM không đủ 40 bit.');
assert(rec_len == 10 && 2^id_bits >= num_channels, 'Bản ghi/Channel ID không hợp lệ.');
assert(ber_levels(1) == 0, 'ber_levels phải bắt đầu bằng 0.');

%% 2. DỮ LIỆU NGUỒN, MUX VÀ KIỂM TRA ĐẾM BIT THEO LOẠI
rng(seed_main,'twister');
source_bits = randi([0 1], num_channels, num_frames*bits_per_sample);

[stdm_tx,stat_tx,lay_stdm,lay_stat,cnt_stdm,cnt_stat] = deal(cell(n_sc,1));
for s = 1:n_sc
    mk = active_masks(s,:);
    stdm_tx{s} = mux_stdm(source_bits,mk,num_frames,bits_per_sample,fas);
    stat_tx{s} = mux_statistical(source_bits,mk,num_frames,bits_per_sample,bits_per_frame);
    lay_stdm{s} = layout_stdm(mk,fas_len,bits_per_sample);
    lay_stat{s} = layout_stat(mk,bits_per_frame,id_bits,bits_per_sample);
    cnt_stdm{s} = check_and_count(stdm_tx{s},lay_stdm{s},source_bits,mk,num_frames,bits_per_sample,fas,true);
    cnt_stat{s} = check_and_count(stat_tx{s},lay_stat{s},source_bits,mk,num_frames,bits_per_sample,fas,false);
end

% Chỉ tiêu dung lượng (từ số bit đếm được, không nhập tay)
payload_bits = zeros(n_sc,1); ov_bits = zeros(n_sc,2); un_bits = zeros(n_sc,2);
for s = 1:n_sc
    payload_bits(s) = cnt_stdm{s}.n3;
    assert(cnt_stat{s}.n3 == payload_bits(s) && payload_bits(s) == n_active_vec(s)*bits_per_sample*num_frames, ...
        'Payload không đồng nhất giữa hai phương pháp.');
    ov_bits(s,:) = [cnt_stdm{s}.n1, cnt_stat{s}.n2];    % FAS ; Channel ID
    un_bits(s,:) = [cnt_stdm{s}.n4, cnt_stat{s}.n5];    % khe rỗng ; padding
end
eff_pct      = repmat(100*payload_bits/total_bits,1,2);   % hiệu suất payload (hai phương pháp)
overhead_pct = 100*ov_bits/total_bits;
unused_pct   = 100*un_bits/total_bits;
kappa_pct    = 100*ov_bits./repmat(payload_bits,1,2);     % overhead / payload

%% 3. MỘT LẦN CHẠY MINH HỌA (seed 42) - CÙNG MẶT NẠ LỖI CHO HAI LUỒNG
% m = 1: STDM FSM cơ bản | m = 2: STDM FSM flywheel | m = 3: StatTDM
rng(seed_main,'twister');
bit_err_1 = zeros(n_sc,3,n_bers);  ber_1 = bit_err_1;  lost_1 = bit_err_1;
sync_1 = zeros(n_sc,2,n_bers);  fl_1 = sync_1;  mis_1 = zeros(n_sc,n_bers);
rec0 = cell(n_sc,3);
for s = 1:n_sc
    mk = active_masks(s,:);
    for b = 1:n_bers
        noise = rand(1,total_bits) < ber_levels(b);        % một mặt nạ lỗi dùng chung
        [recs,vals,sy,fl,mi] = receive_all(xor(stdm_tx{s},noise),xor(stat_tx{s},noise),mk,P);
        sync_1(s,:,b) = sy;  fl_1(s,:,b) = fl;  mis_1(s,b) = mi;
        for m = 1:3
            [bit_err_1(s,m,b),ber_1(s,m,b)] = measure_payload_errors(source_bits,recs{m},vals{m},mk,num_frames,bits_per_sample);
            lost_1(s,m,b) = sum(~vals{m});
            if b == 1, rec0{s,m} = recs{m}; end
        end
    end
end

%% 4. CHẠY LẶP 50 LẦN (seed 1..50) - THỐNG KÊ BER, KHUNG MẤT, MẤT ĐỒNG BỘ, KHÓA NHẦM
ber_tr = zeros(n_sc,3,n_bers,n_trials);  lost_tr = ber_tr;
sync_tr = zeros(n_sc,2,n_bers,n_trials);  fl_tr = sync_tr;  mis_tr = zeros(n_sc,n_bers,n_trials);
fprintf('Đang chạy lặp %d lần (có thể mất vài chục giây)...\n',n_trials);
for t = 1:n_trials
    rng(t,'twister');   % seed t chỉ sinh mặt nạ lỗi; dữ liệu nguồn giữ nguyên (seed_main)
    for s = 1:n_sc
        mk = active_masks(s,:);
        for b = 1:n_bers
            noise = rand(1,total_bits) < ber_levels(b);
            if b == 1, continue; end      % kênh lý tưởng: kết quả tất định, đã kiểm tra ở mục 5
            [recs,vals,sy,fl,mi] = receive_all(xor(stdm_tx{s},noise),xor(stat_tx{s},noise),mk,P);
            for m = 1:3
                [~,ber_tr(s,m,b,t)] = measure_payload_errors(source_bits,recs{m},vals{m},mk,num_frames,bits_per_sample);
                lost_tr(s,m,b,t) = sum(~vals{m});
            end
            sync_tr(s,:,b,t) = sy;  fl_tr(s,:,b,t) = fl;  mis_tr(s,b,t) = mi;
        end
    end
    if mod(t,10) == 0, fprintf('  ... xong %d/%d lần\n',t,n_trials); end
end
ber_mean  = mean(ber_tr,4);  ber_std = std(ber_tr,0,4);
ber_min   = min(ber_tr,[],4);  ber_max = max(ber_tr,[],4);
lost_mean = mean(lost_tr,4);  sync_mean = mean(sync_tr,4);
fl_mean   = mean(fl_tr,4);  mis_mean = mean(mis_tr,3);

%% 5. KIỂM TRA KÊNH LÝ TƯỞNG (BER = 0) - PASS/FAIL TÍNH TỪ KẾT QUẢ
chk = false(n_sc,3);
for s = 1:n_sc
    mk = active_masks(s,:);
    for m = 1:3
        r = rec0{s,m};
        chk(s,m) = bit_err_1(s,m,1) == 0 && lost_1(s,m,1) == 0 && isequal(r(mk,:),source_bits(mk,:)) && ~any(any(r(~mk,:)));
        if m < 3, chk(s,m) = chk(s,m) && sync_1(s,m,1) == 0; end
    end
end

%% 6. IN KẾT QUẢ RA COMMAND WINDOW
p_plot = ber_levels(2:end);
fprintf('\n=================================================================\n');
fprintf(' THÔNG SỐ KHỞI TẠO\n');
fprintf('=================================================================\n');
fprintf('Số kênh                : %d (CH1..CH4)\n',num_channels);
fprintf('Mẫu mỗi kênh / khung   : %d bit\n',bits_per_sample);
fprintf('Tần số khung / chu kỳ  : %d khung/s ; %.0f us\n',frame_rate,frame_period_us);
fprintf('Độ dài khung           : %d bit (STDM: FAS %d + 4 x %d ; StatTDM: bản ghi %d bit + padding)\n',bits_per_frame,fas_len,bits_per_sample,rec_len);
fprintf('Tốc độ đường truyền    : %.0f kbps ; tốc độ danh định mỗi kênh %d kbps\n',line_rate/1000,bits_per_sample*frame_rate/1000);
fprintf('Từ mã FAS              : %s\n',sprintf('%d',fas));
fprintf('Số khung mô phỏng      : %d (%d bit, %.3f s)\n',num_frames,total_bits,num_frames/frame_rate);
fprintf('Mức BER kênh BSC       : %s\n',strtrim(sprintf('%g ',ber_levels)));
fprintf('Bộ sinh ngẫu nhiên     : rng(%d,''twister''); lặp %d lần với seed 1..%d\n',seed_main,n_trials,n_trials);
fprintf('\n=== BẢNG CHỈ TIÊU DUNG LƯỢNG (BER = 0, %d khung, %d bit) - dữ liệu Hình 5.1 ===\n',num_frames,total_bits);
fprintf('%-18s %-8s %9s %9s %9s %9s %7s %7s %7s\n','Kịch bản','PP','Payload','Overhead','Rỗng/Pad','Tổng bit','eta','rho_OH','kappa');
for s = 1:n_sc
    for m = 1:2
        nm = {'STDM','StatTDM'};
        fprintf('%-18s %-8s %9d %9d %9d %9d %6.1f%% %6.1f%% %6.1f%%\n',scenario_names{s},nm{m}, ...
            payload_bits(s),ov_bits(s,m),un_bits(s,m),total_bits,eff_pct(s,m),overhead_pct(s,m),kappa_pct(s,m));
    end
end
fprintf('(PP=phương pháp; Overhead: STDM = FAS, StatTDM = Channel ID; Rỗng/Pad: khe rỗng của STDM, padding của StatTDM)\n');
fprintf('Kiểm tra: payload + overhead + rỗng/pad = %d bit ở mọi kịch bản và cả hai phương pháp: PASS\n',total_bits);
fprintf('\n=== TÍNH LÝ THUYẾT VỚI N_f = 100 KHUNG (4.000 bit) - Bảng 5.1 và 5.2 ===\n');
nf_th = 100;
fprintf('%-18s %9s %9s %9s %9s %9s\n','Kịch bản','Payload','FAS','Rỗng STDM','ID Stat','Pad Stat');
for s = 1:n_sc
    n = n_active_vec(s);
    fprintf('%-18s %9d %9d %9d %9d %9d\n',scenario_names{s},n*8*nf_th,fas_len*nf_th,(4-n)*8*nf_th,2*n*nf_th,(40-10*n)*nf_th);
end
fprintf('Số lỗi kỳ vọng E[Ne] = M*p với M = %d bit payload:  ',4*8*nf_th);
fprintf('p=%g -> %.2f ;  ',[p_plot; 4*8*nf_th*p_plot]); fprintf('\n');
fprintf('\n=== MINH HỌA KHUNG ĐẦU TIÊN (KỊCH BẢN A, BER = 0) ===\n');
fprintf('STDM    TX (40 bit): %s\n',sprintf('%d',stdm_tx{1}(1:bits_per_frame)));
fprintf('StatTDM TX (40 bit): %s\n',sprintf('%d',stat_tx{1}(1:bits_per_frame)));
fprintf('\n16 bit đầu mỗi kênh: nguồn và sau DEMUX (kịch bản A, BER = 0)\n');
for ch = 1:num_channels
    idx = 1:16;
    fprintf('CH%d nguồn   : %s\n',ch,sprintf('%d',source_bits(ch,idx)));
    fprintf('CH%d STDM RX : %s\n',ch,sprintf('%d',rec0{1,1}(ch,idx)));
    fprintf('CH%d Stat RX : %s\n',ch,sprintf('%d',rec0{1,3}(ch,idx)));
end
fprintf('\n=== KIỂM TRA KHÔI PHỤC DỮ LIỆU TẠI BER = 0 ===\n');
fprintf('%-18s %-10s %-10s %-10s\n','Kịch bản','STDM cơ bản','STDM flywheel','StatTDM');
for s = 1:n_sc
    fprintf('%-18s %-10s %-13s %-10s\n',scenario_names{s},pass_fail(chk(s,1)),pass_fail(chk(s,2)),pass_fail(chk(s,3)));
end
fprintf('Kết luận STDM    : %s\n',pass_fail(all(chk(:,1)) && all(chk(:,2))));
fprintf('Kết luận StatTDM : %s\n',pass_fail(all(chk(:,3))));
fprintf('\n=== BER ĐẶT VÀO VÀ BER ĐO SAU DEMUX - MỘT LẦN CHẠY (seed %d) ===\n',seed_main);
fprintf('%-18s %-20s %9s %9s %11s %9s %9s %9s\n','Kịch bản','Phương pháp','p đặt','Bit lỗi','BER đo','Khung mất','Mất ĐB','Sai ID');
for b = 1:n_bers
    for s = 1:n_sc
        for m = 1:3
            if m < 3, sl = sync_1(s,m,b); else, sl = 0; end
            if m == 3, mi = mis_1(s,b); else, mi = 0; end      % Sai ID: số bản ghi StatTDM định tuyến sai
            fprintf('%-18s %-20s %9.1e %9d %11.4g %9d %9d %9d\n',scenario_names{s},method_names{m}, ...
                ber_levels(b),bit_err_1(s,m,b),ber_1(s,m,b),lost_1(s,m,b),sl,mi);
        end
    end
end
fprintf('\n=== THỐNG KÊ %d LẦN CHẠY - STDM, FSM cơ bản (Bảng 4.2) ===\n',n_trials);
fprintf('%-4s %8s %11s %20s %7s %10s %9s %9s\n','KB','p','BER TB','BER min - max','BER/p','Khung mất','Mất ĐB','Khóa nhầm');
for b = 2:n_bers
    for s = 1:n_sc
        fprintf('%-4s %8.0e %11.3e %9.1e - %-8.1e %7.1f %10.1f %9.1f %9.1f\n',scenario_names{s}(1),ber_levels(b), ...
            ber_mean(s,1,b),ber_min(s,1,b),ber_max(s,1,b),ber_mean(s,1,b)/ber_levels(b), ...
            lost_mean(s,1,b),sync_mean(s,1,b),fl_mean(s,1,b));
    end
end
fprintf('\n=== THỐNG KÊ %d LẦN CHẠY - STDM, FSM flywheel (n_miss = %d) ===\n',n_trials,n_miss);
fprintf('%-4s %8s %11s %7s %10s %9s %9s\n','KB','p','BER TB','BER/p','Khung mất','Mất ĐB','Khóa nhầm');
for b = 2:n_bers
    for s = 1:n_sc
        fprintf('%-4s %8.0e %11.3e %7.2f %10.1f %9.2f %9.2f\n',scenario_names{s}(1),ber_levels(b), ...
            ber_mean(s,2,b),ber_mean(s,2,b)/ber_levels(b),lost_mean(s,2,b),sync_mean(s,2,b),fl_mean(s,2,b));
    end
end
fprintf('\n=== THỐNG KÊ %d LẦN CHẠY - StatTDM (Bảng 4.5) ===\n',n_trials);
fprintf('%-4s %8s %11s %20s %7s %26s\n','KB','p','BER TB','BER min - max','BER/p','Bản ghi ID sai / tổng bản ghi');
for b = 2:n_bers
    for s = 1:n_sc
        nrec = n_active_vec(s)*num_frames;
        fprintf('%-4s %8.0e %11.3e %9.1e - %-8.1e %7.2f %14.1f / %d (%.2f%%)\n',scenario_names{s}(1),ber_levels(b), ...
            ber_mean(s,3,b),ber_min(s,3,b),ber_max(s,3,b),ber_mean(s,3,b)/ber_levels(b),mis_mean(s,b),nrec,100*mis_mean(s,b)/nrec);
    end
end
fprintf('\n=== STDM: BER LÝ THUYẾT (mục 4.2) VÀ BER ĐO - KỊCH BẢN A ===\n');
fprintf('%8s %12s %12s %20s %14s\n','p','q=1-(1-p)^8','BER lý thuyết','BER đo TB ± sigma (A)','BER đo TB (D)');
for b = 2:n_bers
    p = ber_levels(b); q = 1-(1-p)^fas_len;
    fprintf('%8.0e %12.3e %12.3e %12.3e ± %-7.1e %12.3e\n',p,q,q+(1-q)*p,ber_mean(1,1,b),ber_std(1,1,b),ber_mean(4,1,b));
end
fprintf('\n=== HIỆU SUẤT HIỆU DỤNG STDM = eta*(1 - BER) - KỊCH BẢN A (Bảng 4.4) ===\n');
fprintf('%8s %8s %12s %12s %12s %12s\n','p','eta','BER (cơ bản)','eta_hd','BER (flywheel)','eta_hd');
for b = 2:n_bers
    fprintf('%8.0e %7.1f%% %12.3e %11.2f%% %12.3e %11.2f%%\n',ber_levels(b),eff_pct(1,1), ...
        ber_mean(1,1,b),eff_pct(1,1)*(1-ber_mean(1,1,b)),ber_mean(1,2,b),eff_pct(1,1)*(1-ber_mean(1,2,b)));
end
fprintf('\nLưu ý: BER là đại lượng thống kê; kết quả 1 lần chạy (seed %d) chỉ để minh họa, dùng trung bình %d lần để kết luận.\n',seed_main,n_trials);

%% 7. DỮ LIỆU CHO ĐỒ THỊ
ber_basic_mean = reshape(ber_mean(:,1,2:end),n_sc,[]);
ber_fly_mean   = reshape(ber_mean(:,2,2:end),n_sc,[]);
ber_stat_mean  = reshape(ber_mean(:,3,2:end),n_sc,[]);
ber_basic_std  = reshape(ber_std(:,1,2:end),n_sc,[]);
ber_fly_std    = reshape(ber_std(:,2,2:end),n_sc,[]);
lost_basic     = reshape(lost_mean(:,1,2:end),n_sc,[]);
sync_basic     = reshape(sync_mean(:,1,2:end),n_sc,[]);
mis_plot       = reshape(mis_mean(:,2:end),n_sc,[]);
ber_basic_s42  = reshape(ber_1(1,1,2:end),1,[]);
ber_stat_s42   = reshape(ber_1(1,3,2:end),1,[]);

ber_colors = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
ber_legend = {'p = 10^{-4} (0,01%)','p = 10^{-3} (0,1%)','p = 10^{-2} (1%)'};
load_colors = lines(4);
tick_lab = load_labels;
xl_load = 'Mức tải (4, 3, 2, 1 kênh hoạt động)';
xl_ber = 'BER đặt vào kênh BSC (p)';
ss = get(0,'ScreenSize');

%% ==================== CỬA SỔ 1 - SO SÁNH STDM vs StatTDM ====================
f1 = figure('Name','Cửa sổ 1 - So sánh STDM và StatTDM','NumberTitle','off','Color','w', ...
    'Position',[ss(1)+ss(3)*0.01 ss(2)+ss(4)*0.03 ss(3)*0.88 ss(4)*0.90]);
sc_ex = 3;  mk = active_masks(sc_ex,:);  nb_show = 24;     % kịch bản minh họa C (2 kênh): có cả khe rỗng và padding
draw_frames(make_ax(f1,4,3,1,1,1,3),{stdm_tx{sc_ex}(1:bits_per_frame),stat_tx{sc_ex}(1:bits_per_frame)}, ...
    {lay_stdm{sc_ex},lay_stat{sc_ex}},{'STDM','StatTDM'},true, ...
    sprintf('Cấu trúc khung đầu tiên (bit thực tế do MUX tạo ra) - kịch bản %s',scenario_names{sc_ex}));
grids = {source_bits(:,1:nb_show),rec0{sc_ex,1}(:,1:nb_show),rec0{sc_ex,3}(:,1:nb_show)};
gtit  = {'Dữ liệu nguồn (24 bit đầu, BER = 0)','Sau DEMUX - STDM','Sau DEMUX - StatTDM'};
for k = 1:3, draw_bitgrid(make_ax(f1,4,3,2,k,1,1),grids{k},mk,gtit{k}); end

bar_panel(make_ax(f1,4,3,3,1,1,1),eff_pct,'Hình 5.2. Hiệu suất payload (BER = 0)','Hiệu suất payload (%)', ...
    {'STDM','StatTDM'},'northeast',[0 100],tick_lab,xl_load);
bar_panel(make_ax(f1,4,3,3,2,1,1),overhead_pct,'Hình 5.3. Tỷ lệ overhead (BER = 0)','Overhead / tổng dung lượng (%)', ...
    {'STDM - FAS','StatTDM - Channel ID'},'northeast',[0 35],tick_lab,xl_load);
bar_panel(make_ax(f1,4,3,3,3,1,1),unused_pct,'Hình 5.4. Khe rỗng (STDM) và padding (StatTDM)','Không mang payload (%)', ...
    {'STDM - khe rỗng','StatTDM - padding'},'northwest',[0 100],tick_lab,xl_load);

ax = make_ax(f1,4,3,4,1,1,1);
loglog(ax,p_plot,ber_basic_mean(1,:),'-o','LineWidth',1.5,'MarkerSize',6); hold(ax,'on');
loglog(ax,p_plot,ber_fly_mean(1,:),'-s','LineWidth',1.5,'MarkerSize',6);
loglog(ax,p_plot,ber_stat_mean(1,:),'-d','LineWidth',1.5,'MarkerSize',6);
loglog(ax,p_plot,p_plot,'--k','LineWidth',1.1);
grid(ax,'on'); box(ax,'on'); xlabel(ax,xl_ber); ylabel(ax,'BER payload sau DEMUX');
title(ax,sprintf('Hình 4.4. BER sau DEMUX - trung bình %d lần (4 kênh)',n_trials));
legend(ax,'STDM - FSM cơ bản','STDM - FSM flywheel','StatTDM','Y = X','Location','northwest');

ax = make_ax(f1,4,3,4,2,1,1);
semilogx(ax,p_plot,ber_basic_s42,'-o','LineWidth',1.5,'MarkerSize',6); hold(ax,'on');
semilogx(ax,p_plot,ber_stat_s42,'-s','LineWidth',1.5,'MarkerSize',6);
semilogx(ax,p_plot,p_plot,'--k','LineWidth',1.1);
xlim(ax,[min(p_plot) max(p_plot)]); ylim(ax,[0 max(0.1,1.15*max([ber_basic_s42 ber_stat_s42]))]);
grid(ax,'on'); box(ax,'on'); xlabel(ax,xl_ber); ylabel(ax,'BER payload sau DEMUX');
title(ax,sprintf('Hình 5.5. BER một lần chạy (seed %d, 4 kênh)',seed_main));
legend(ax,'STDM - FSM cơ bản','StatTDM','Y = X','Location','northwest');

bar_panel(make_ax(f1,4,3,4,3,1,1),kappa_pct,'Overhead trên mỗi bit payload \kappa_{OH} (BER = 0)','Overhead / payload (%)', ...
    {'STDM','StatTDM'},'northwest',[0 120],tick_lab,xl_load);

%% ==================== CỬA SỔ 2 - RIÊNG StatTDM ====================
f2 = figure('Name','Cửa sổ 2 - TDM thống kê (StatTDM)','NumberTitle','off','Color','w', ...
    'Position',[ss(1)+ss(3)*0.05 ss(2)+ss(4)*0.02 ss(3)*0.88 ss(4)*0.85]);
draw_frames(make_ax(f2,2,3,1,1,1,3),cellfun(@(x)x(1:bits_per_frame),stat_tx','UniformOutput',false),lay_stat',load_names,false, ...
    'StatTDM - khung đầu tiên theo mức tải: bản ghi 10 bit (Channel ID 2 bit + payload 8 bit), phần còn lại là padding');

ax = make_ax(f2,2,3,2,1,1,1);
for s = 1:n_sc
    loglog(ax,p_plot,ber_stat_mean(s,:),'-o','LineWidth',1.5,'MarkerSize',5,'Color',load_colors(s,:)); hold(ax,'on');
end
loglog(ax,p_plot,p_plot,'--k','LineWidth',1.0); grid(ax,'on'); box(ax,'on');
xlabel(ax,xl_ber); ylabel(ax,'BER payload StatTDM sau DEMUX');
title(ax,sprintf('BER StatTDM theo mức tải (TB %d lần)',n_trials));
legend(ax,[load_names {'Y = X'}],'Location','northwest');

grp_panel(make_ax(f2,2,3,2,2,1,1),mis_plot,'Bản ghi bị định tuyến sai do lỗi Channel ID', ...
    sprintf('Bản ghi định tuyến sai / %d khung',num_frames),tick_lab,xl_load,ber_colors,ber_legend);

ax = make_ax(f2,2,3,2,3,1,1);
eta_stat_eff = repmat(eff_pct(:,2),1,numel(p_plot)).*(1-ber_stat_mean);
for s = 1:n_sc
    semilogx(ax,p_plot,eta_stat_eff(s,:),'-o','LineWidth',1.5,'MarkerSize',5,'Color',load_colors(s,:)); hold(ax,'on');
    semilogx(ax,p_plot,eff_pct(s,2)*ones(size(p_plot)),':','Color',load_colors(s,:),'HandleVisibility','off');
end
grid(ax,'on'); box(ax,'on'); xlim(ax,[min(p_plot) max(p_plot)]); ylim(ax,[0 100]); set(ax,'XTick',p_plot);
xlabel(ax,xl_ber); ylabel(ax,'Hiệu suất hiệu dụng (%)');
title(ax,'Hiệu suất hiệu dụng StatTDM = \eta(1 - BER)'); legend(ax,load_names,'Location','southwest');

%% ==================== CỬA SỔ 3 - RIÊNG STDM ====================
f3 = figure('Name','Cửa sổ 3 - TDM đồng bộ (STDM)','NumberTitle','off','Color','w', ...
    'Position',[ss(1)+ss(3)*0.09 ss(2)+ss(4)*0.01 ss(3)*0.88 ss(4)*0.88]);
draw_frames(make_ax(f3,3,2,1,1,1,2),cellfun(@(x)x(1:bits_per_frame),stdm_tx','UniformOutput',false),lay_stdm',load_names,true, ...
    'STDM - khung đầu tiên theo mức tải: FAS 8 bit + 4 khe x 8 bit (kênh không hoạt động giữ khe rỗng)');

ax = make_ax(f3,3,2,2,1,1,1);
q_fas = 1-(1-p_plot).^fas_len;
ber_theory = q_fas + (1-q_fas).*p_plot;
loglog(ax,p_plot,ber_theory,'-k','LineWidth',1.5); hold(ax,'on');
loglog(ax,p_plot,p_plot,'--k','LineWidth',1.0);
loglog(ax,p_plot,ber_basic_s42,'rx','MarkerSize',9,'LineWidth',1.7);
errorbar(ax,p_plot,ber_basic_mean(1,:),min(ber_basic_std(1,:),0.99*ber_basic_mean(1,:)),ber_basic_std(1,:),'-o','LineWidth',1.3);
errorbar(ax,p_plot,ber_fly_mean(1,:),min(ber_fly_std(1,:),0.99*ber_fly_mean(1,:)),ber_fly_std(1,:),'-s','LineWidth',1.3);
grid(ax,'on'); box(ax,'on'); xlabel(ax,xl_ber); ylabel(ax,'BER payload sau DEMUX');
title(ax,'Hình 4.1. BER sau DEMUX của STDM (kịch bản A)');
legend(ax,'Lý thuyết FSM cơ bản (mục 4.2)','Tham chiếu Y = X',sprintf('FSM cơ bản - seed %d',seed_main), ...
    sprintf('FSM cơ bản - TB %d lần',n_trials),sprintf('FSM flywheel - TB %d lần',n_trials),'Location','northwest','FontSize',8);

grp_panel(make_ax(f3,3,2,2,2,1,1),lost_basic,sprintf('Hình 4.2. Số khung STDM bị mất (FSM cơ bản, TB %d lần)',n_trials), ...
    sprintf('Số khung bị mất / %d khung',num_frames),tick_lab,xl_load,ber_colors,ber_legend);

ax = make_ax(f3,3,2,3,1,1,1);
eta_A = eff_pct(1,1);
semilogx(ax,p_plot,eta_A*(1-ber_basic_mean(1,:)),'-o','LineWidth',1.5,'MarkerSize',6); hold(ax,'on');
semilogx(ax,p_plot,eta_A*(1-ber_fly_mean(1,:)),'-s','LineWidth',1.5,'MarkerSize',6);
semilogx(ax,p_plot,eta_A*ones(size(p_plot)),'--k','LineWidth',1.0);
grid(ax,'on'); box(ax,'on'); xlim(ax,[min(p_plot) max(p_plot)]); ylim(ax,[0 85]); set(ax,'XTick',p_plot);
xlabel(ax,xl_ber); ylabel(ax,'Hiệu suất hiệu dụng (%)');
title(ax,'Hiệu suất hiệu dụng STDM (tải 100%, Bảng 4.4)');
legend(ax,'FSM cơ bản','FSM flywheel','\eta khi p = 0','Location','southwest');

grp_panel(make_ax(f3,3,2,3,2,1,1),sync_basic,'Số lần FSM cơ bản của STDM mất đồng bộ', ...
    sprintf('Số lần mất đồng bộ / %d khung',num_frames),tick_lab,xl_load,ber_colors,ber_legend);

%% ============================ HÀM CỤC BỘ ============================
function ax = make_ax(fig,nr,nc,r,c,rs,cs)
% Tạo axes theo lưới nr x nc, bắt đầu ở ô (r,c), trải rs hàng và cs cột.
    L = 0.055; R = 0.015; B = 0.075; T = 0.05; gx = 0.05; gy = 0.095;
    w = (1-L-R-(nc-1)*gx)/nc;
    h = (1-B-T-(nr-1)*gy)/nr;
    W = cs*w + (cs-1)*gx;
    H = rs*h + (rs-1)*gy;
    x = L + (c-1)*(w+gx);
    y = 1 - T - (r-1)*(h+gy) - H;
    ax = axes('Parent',fig,'Position',[x y W H],'FontSize',9);
end

function col = field_colors()
% 1 FAS | 2 Channel ID | 3 payload | 4 khe rỗng | 5 padding
    col = [0.60 0.78 0.95; 0.86 0.72 0.95; 0.62 0.86 0.62; 0.90 0.90 0.90; 0.75 0.75 0.75];
end

function lay = layout_stdm(mask,faslen,bps)
% Loại của từng bit trong khung STDM (1 = FAS, 3 = payload, 4 = khe rỗng).
    lay = [ones(1,faslen) repelem(3+~mask(:)',bps)];
end

function lay = layout_stat(mask,bpf,idb,bps)
% Loại của từng bit trong khung StatTDM (2 = Channel ID, 3 = payload, 5 = padding).
    n = nnz(mask);
    lay = [repmat([2*ones(1,idb) 3*ones(1,bps)],1,n) 5*ones(1,bpf-n*(idb+bps))];
end

function c = check_and_count(tx,lay,source,mask,nf,bps,fas,is_stdm)
% Kiểm tra luồng TX thật đúng cấu trúc khung rồi đếm bit theo loại.
    types = repmat(lay,1,nf);
    assert(numel(tx) == numel(types),'Độ dài luồng TX không khớp cấu trúc khung.');
    act = find(mask); na = numel(act);
    exp_pay = reshape(permute(reshape(source(act,:),na,bps,nf),[2 1 3]),1,[]);
    assert(isequal(tx(types==3),exp_pay),'Payload trong luồng TX không khớp dữ liệu nguồn.');
    if is_stdm
        assert(isequal(tx(types==1),repmat(fas,1,nf)),'FAS trong luồng TX sai.');
    else
        exp_id = repmat(reshape([floor((act-1)/2); mod(act-1,2)],1,[]),1,nf);
        assert(isequal(tx(types==2),exp_id),'Channel ID trong luồng TX sai.');
    end
    assert(~any(tx(types==4 | types==5)),'Khe rỗng/padding phải bằng 0.');
    c.n1 = nnz(types==1); c.n2 = nnz(types==2); c.n3 = nnz(types==3);
    c.n4 = nnz(types==4); c.n5 = nnz(types==5);
    assert(c.n1+c.n2+c.n3+c.n4+c.n5 == numel(tx),'Tổng bit theo loại không bằng tổng bit truyền.');
end

function txt = pass_fail(tf)
    names = {'FAIL','PASS'}; txt = names{tf+1};
end

%% ---------- MUX ----------
function tx = mux_stdm(source_bits,mask,nf,bps,fas)
% Mỗi khung = FAS + payload các kênh (kênh im lặng giữ 8 bit 0).
    nch = size(source_bits,1);
    S = source_bits .* mask(:);                                           % bit của kênh không hoạt động = 0
    pay = reshape(permute(reshape(S,nch,bps,nf),[2 1 3]),nch*bps,nf);     % mỗi cột: payload một khung
    tx = reshape([repmat(fas(:),1,nf); pay],1,[]);
end

function tx = mux_statistical(source_bits,mask,nf,bps,bpf)
% Mỗi bản ghi = Channel ID (CH1=00, CH2=01, CH3=10, CH4=11) + payload; phần còn lại là padding bằng 0.
    act = find(mask); na = numel(act);
    rec = zeros(2+bps,na,nf);
    rec(1:2,:,:) = repmat([floor((act-1)/2); mod(act-1,2)],1,1,nf);
    rec(3:end,:,:) = permute(reshape(source_bits(act,:),na,bps,nf),[2 1 3]);
    tx = reshape([reshape(rec,na*(2+bps),nf); zeros(bpf-na*(2+bps),nf)],1,[]);
end

%% ---------- DEMUX ----------
function [recs,vals,sync,fl,mis] = receive_all(rx_a,rx_b,mk,P)
% Ba bộ thu: STDM FSM cơ bản, STDM flywheel (cùng luồng rx_a) và StatTDM (luồng rx_b).
    [r1,v1,l1,f1] = demux_stdm_fsm(rx_a,P);
    [r2,v2,l2,f2] = demux_stdm_flywheel(rx_a,P);
    [r3,v3,mis]   = demux_statistical(rx_b,mk,P);
    recs = {r1,r2,r3};  vals = {v1,v2,v3};  sync = [l1 l2];  fl = [f1 f2];
end

function [rec,valid,sl,fl] = demux_stdm_fsm(rx,P)
% FSM cơ bản (2 trạng thái HUNT/SYNC): khóa ngay khi thấy FAS, một FAS sai là về HUNT (n_miss = 1).
    rec = zeros(P.nch,P.nf*P.bps); valid = false(1,P.nf); sl = 0; fl = 0;
    state = 0; pos = 1; L = numel(rx); n = P.fl;          % 0 = HUNT, 1 = SYNC
    while pos+n-1 <= L
        hit = all(rx(pos:pos+n-1) == P.fas);
        if state == 0
            if hit
                state = 1;
                if mod(pos-1,P.bpf) ~= 0, fl = fl + 1; end            % khóa lệch đầu khung
            else
                pos = pos + 1;
            end
        elseif hit
            fi = round((pos-1)/P.bpf) + 1; ps = pos + n;
            if fi >= 1 && fi <= P.nf && ps+P.nch*P.bps-1 <= L
                rec(:,(fi-1)*P.bps+(1:P.bps)) = reshape(rx(ps:ps+P.nch*P.bps-1),P.bps,P.nch)';
                valid(fi) = true;
            end
            pos = pos + P.bpf;
        else
            sl = sl + 1; state = 0; pos = pos + 1;
        end
    end
end

function [rec,valid,sl,fl] = demux_stdm_flywheel(rx,P)
% FSM 4 trạng thái: HUNT -> CONFIRM -> SYNC <-> FLYWHEEL; mất đồng bộ khi P.nm FAS sai liên tiếp (G.706).
    rec = zeros(P.nch,P.nf*P.bps); valid = false(1,P.nf); sl = 0; fl = 0;
    HUNT = 1; CONFIRM = 2; SYNC = 3; FLYWHEEL = 4;
    state = HUNT; pos = 1; cand = 1; miss = 0; L = numel(rx); n = P.fl;
    while pos+n-1 <= L
        if state == HUNT
            if all(rx(pos:pos+n-1) == P.fas), cand = pos; state = CONFIRM; else, pos = pos + 1; end
        elseif state == CONFIRM
            cp = cand + P.bpf;
            if cp+n-1 <= L && all(rx(cp:cp+n-1) == P.fas)
                pos = cand; miss = 0; state = SYNC;
                if mod(cand-1,P.bpf) ~= 0, fl = fl + 1; end
            else
                pos = cand + 1; state = HUNT;
            end
        else    % SYNC hoặc FLYWHEEL
            if all(rx(pos:pos+n-1) == P.fas)
                miss = 0; state = SYNC;
            else
                if state == SYNC, miss = 1; else, miss = miss + 1; end
                state = FLYWHEEL;
                if miss >= P.nm
                    sl = sl + 1; state = HUNT; pos = pos + 1;
                    continue;
                end
            end
            fi = round((pos-1)/P.bpf) + 1; ps = pos + n;
            if fi >= 1 && fi <= P.nf && ps+P.nch*P.bps-1 <= L
                rec(:,(fi-1)*P.bps+(1:P.bps)) = reshape(rx(ps:ps+P.nch*P.bps-1),P.bps,P.nch)';
                valid(fi) = true;
            end
            pos = pos + P.bpf;
        end
    end
end

function [rec,valid,n_mis] = demux_statistical(rx,mask,P)
% DEMUX StatTDM: biết trước active_mask và ranh giới khung; mỗi bản ghi = ID 2 bit + payload; bỏ padding.
    nch = numel(mask); act = find(mask); rl = 2 + P.bps; n_mis = 0;
    rec = zeros(nch,P.nf*P.bps); valid = false(1,P.nf);
    for f = 1:P.nf
        cursor = (f-1)*P.bpf + 1;
        oi = (f-1)*P.bps + (1:P.bps);
        for r = 1:numel(act)
            if cursor+rl-1 > numel(rx), break; end
            ch = 2*rx(cursor) + rx(cursor+1) + 1;
            if ch ~= act(r), n_mis = n_mis + 1; end              % ID thu khác kênh phát
            if ch >= 1 && ch <= nch && mask(ch), rec(ch,oi) = rx(cursor+2:cursor+rl-1); end
            cursor = cursor + rl;
        end
        valid(f) = true;
    end
end

function [n_err,ber] = measure_payload_errors(source,recovered,valid,mask,nframes,bps)
% BER trên payload của các kênh hoạt động; khung bị mất tính sai toàn bộ payload.
    act = find(mask); na = numel(act);
    mism = source(act,:) ~= recovered(act,:);
    fe = sum(reshape(sum(mism,1),bps,nframes),1);       % số bit sai theo từng khung
    n_err = sum(fe(valid)) + na*bps*sum(~valid);
    ber = n_err/(na*nframes*bps);
end

%% ---------- VẼ ----------
function bar_panel(ax,Y,ttl,ylab,leg,loc,yl,tick_lab,xl)
% Biểu đồ cột theo mức tải (STDM vs StatTDM) kèm nhãn giá trị.
    bh = bar(ax,Y); grid(ax,'on'); box(ax,'on'); ylim(ax,yl);
    set(ax,'XTick',1:size(Y,1),'XTickLabel',tick_lab); xlabel(ax,xl); ylabel(ax,ylab);
    title(ax,ttl); legend(ax,leg{:},'Location',loc); label_bars(ax,bh,'%.0f',false);
end

function grp_panel(ax,Y,ttl,ylab,tick_lab,xl,cols,leg)
% Biểu đồ cột nhóm theo mức BER (p = 10^-4, 10^-3, 10^-2).
    bh = bar(ax,Y,'grouped'); grid(ax,'on'); box(ax,'on');
    for k = 1:numel(bh), bh(k).FaceColor = cols(k,:); end
    set(ax,'XTick',1:size(Y,1),'XTickLabel',tick_lab); xlabel(ax,xl); ylabel(ax,ylab);
    title(ax,ttl); legend(ax,bh,leg,'Location','northwest','FontSize',8); label_bars(ax,bh,'%.1f',true);
end

function draw_frames(ax,bits_rows,types_rows,names,slot8,ttl)
% Vẽ từng bit của khung đầu tiên, tô màu theo trường (FAS/ID/payload/khe rỗng/padding).
    cla(ax); hold(ax,'on');
    col = field_colors(); n = numel(names); nb = numel(bits_rows{1});
    for i = 1:n
        y = n - i + 1;
        bits = bits_rows{i}; ty = types_rows{i};
        for k = 1:nb
            rectangle('Parent',ax,'Position',[k-1 y-0.9 1 0.8],'FaceColor',col(ty(k),:),'EdgeColor',[0.5 0.5 0.5]);
            text(k-0.5,y-0.5,sprintf('%d',bits(k)),'Parent',ax,'HorizontalAlignment','center', ...
                'FontName','Courier New','FontSize',9,'FontWeight','bold','Interpreter','none');
        end
        for k = 1:nb-1
            if ty(k) ~= ty(k+1) || (slot8 && mod(k,8) == 0)
                plot(ax,[k k],[y-0.95 y-0.05],'k-','LineWidth',1.8);
            end
        end
        text(-0.6,y-0.5,names{i},'Parent',ax,'HorizontalAlignment','right','FontSize',9, ...
            'FontWeight','bold','Interpreter','none');
    end
    % chú giải màu vẽ ngay trong axes
    items = {1,'FAS';2,'Channel ID';3,'Payload';4,'Khe rỗng (STDM)';5,'Padding (StatTDM)'};
    x0 = 0;
    for j = 1:size(items,1)
        rectangle('Parent',ax,'Position',[x0 -1.0 1.2 0.7],'FaceColor',col(items{j,1},:),'EdgeColor','k');
        text(x0+1.6,-0.65,items{j,2},'Parent',ax,'FontSize',8,'Interpreter','none');
        x0 = x0 + 8.4;
    end
    xlim(ax,[-9 nb+0.5]); ylim(ax,[-1.3 n+0.1]);
    set(ax,'XTick',[],'YTick',[],'XColor','none','YColor','none','Box','off');
    title(ax,ttl,'FontSize',10,'Interpreter','none');
end

function draw_bitgrid(ax,M,mask,ttl)
% Lưới bit của 4 kênh; kênh không hoạt động tô xám và đánh dấu riêng.
    nb = size(M,2);
    C = double(M); C(~mask,:) = 2;
    imagesc(ax,C); hold(ax,'on');
    colormap(ax,[1 1 1; 0.62 0.86 0.62; 0.82 0.82 0.82]); caxis(ax,[0 2]);
    for ch = 1:size(M,1)
        for k = 1:nb
            if mask(ch)
                text(k,ch,sprintf('%d',M(ch,k)),'Parent',ax,'HorizontalAlignment','center','FontSize',7,'Interpreter','none');
            else
                text(k,ch,'-','Parent',ax,'HorizontalAlignment','center','FontSize',7,'Color',[0.45 0.45 0.45],'Interpreter','none');
            end
        end
    end
    for k = 8:8:nb-1, plot(ax,[k+0.5 k+0.5],[0.5 size(M,1)+0.5],'k-','LineWidth',1.2); end
    set(ax,'XTick',4.5:8:nb,'XTickLabel',arrayfun(@(v)sprintf('mẫu %d',v),1:ceil(nb/8),'UniformOutput',false), ...
        'YTick',1:size(M,1),'YTickLabel',{'CH1','CH2','CH3','CH4'},'TickLength',[0 0]);
    xlabel(ax,'Ô xám: kênh không hoạt động','FontSize',8);
    title(ax,ttl,'FontSize',10,'Interpreter','none');
end

function label_bars(ax,bh,fmt,adjust)
% Ghi giá trị lên đầu cột (bỏ qua nếu MATLAB cũ không có XEndPoints).
    if ~isprop(bh(1),'XEndPoints') || ~isprop(bh(1),'YEndPoints'), return; end
    ymax = 0;
    for j = 1:numel(bh), ymax = max(ymax,max(bh(j).YEndPoints)); end
    if ymax <= 0, return; end
    dy = 0.012*ymax;
    for j = 1:numel(bh)
        x = bh(j).XEndPoints; y = bh(j).YEndPoints;
        lab = arrayfun(@(v)sprintf(fmt,v),y,'UniformOutput',false);
        text(x,y+dy,lab,'Parent',ax,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7);
    end
    if adjust, ylim(ax,[0 ymax*1.5]); end   % chừa chỗ cho chú giải ở góc trên trái
end