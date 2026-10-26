function cfg = state_space_parameters()
%STATE_SPACE_PARAMETERS 状态空间复现的全部参数与实验配置。
%
% 使用方法：
%   result = run_state_space_repro();
% 或者：
%   cfg = state_space_parameters();
%   cfg.experiment.tableSCR = 10;
%   result = run_state_space_repro(cfg);
%
% 本文件只负责设置参数，不加载、不修改原 Simulink 模型。
% 硬件与结构来源：用户附件 GFMUPS100k_CPP_estimate.slx。
% 控制器与运行点采用当前文件中的值；保留你后续手动修改的参数。
% 论文：刘永奎等，《虚拟同步发电机型储能逆变器次同步振荡抑制策略》，
% 电力系统及其自动化学报，2026，38(2)：126-134。

%% 1. 图7、表2、表3与首次实验的 SCR 设置

cfg.experiment.firstSCR = 10;           % 首先计算的工况，不是只在图上标记
cfg.experiment.figure7SCR = 2:0.05:40;   % 图7扫描范围，可直接修改
cfg.experiment.tableSCR = 30;           % 表2：计算全部特征值的 SCR
cfg.experiment.participationSCR = 30;   % 表3：计算参与因子的 SCR

% false：只输出与你附件结构更接近的主模型。
% true ：为另一结构也生成完整SCR扫描、表2、表3和图形，并追加比较图。
% 不再仅计算另一结构的单个工况。
cfg.experiment.compareStructures = false;

% paper15：论文正文的15阶结构，代入用户参数；保留论文PCC电压前馈，
%          无功反馈为 v_pd，忽略电容ESR与虚拟阻抗电流预滤波。
% model17：在15阶骨架上增加虚拟阻抗的两个500Hz滤波状态，采用附件的
%          电压幅值反馈、无PCC电压前馈、滤波电感电阻和电容ESR。
% 注意：model17仍是连续平均近似，不包含RMS滑窗、PWM、采样保持与预同步。
cfg.model.structure = 'model17';    % 默认主模型：附件控制路径的17阶平均近似
% paper15只作为结构对照，不用于替代原附件模型的稳定性判断。

%% 2. 三相额定值与电网阻抗换算

cfg.base.Ull = 380;                 % 电网线电压有效值，V
cfg.base.Sn = 100e3;                % 三相额定视在功率基准，VA
% 100kVA由附件SCR=1/1.5/3/10/30支路交叉反推，约为100000VA。
% 此处是额定基准，不是实时输出功率。Pref=0不会使SCR无穷大。

cfg.base.f0 = 50;                   % 额定频率，Hz
cfg.base.w0 = 2*pi*cfg.base.f0;      % 额定角频率，rad/s
cfg.base.Eg = sqrt(2/3)*cfg.base.Ull; % 电网相电压峰值，V
cfg.base.Zb = cfg.base.Ull^2/cfg.base.Sn;
cfg.base.Ipk = sqrt(2)*cfg.base.Sn/(sqrt(3)*cfg.base.Ull);

cfg.grid.XR = 10;                  % X/R，附件SCR=10支路约为10
% 每一个SCR重新计算：
%   |Zg| = Ull^2 / (Sn*SCR)
%   Rg   = |Zg| / sqrt(1 + (X/R)^2)
%   Xg   = Rg*(X/R)，Lg = Xg/w0
% 不使用原模型支路名称来判断其真实SCR。

%% 3. 主电路与平均化范围

cfg.plant.Lf = 182e-6;             % LC滤波电感，每相，H
cfg.plant.Cf = 28e-6;              % LC滤波电容，每相，F
cfg.plant.Rf = 0.01;               % 电感串联电阻，每相，ohm
cfg.plant.Rc = 0.01;               % 电容串联ESR，每相，ohm，model17使用
cfg.plant.Udc = 1500;              % 直流母线电压，V，作为理想恒定电压
cfg.plant.fs = 20e3;               % PWM频率，Hz，仅用于说明平均模型范围
cfg.plant.Ts = 1/cfg.plant.fs;      % 内环控制周期，s
% 附件LC子系统内孤立的5uH支路没有物理端口连线，不计入Lf或Lg。
% 当前研究工况为三相平衡、并网、PCC无本地负载、预同步已退出。

%% 4. 有功环：避免把12000直接误当成论文 Dp

cfg.control.Pref = 90e3;            % 当前有功参考，W；采用你目前的90kW设置
cfg.control.Qref = 0;               % 无功参考，var
cfg.control.UrefLL = 380;           % 无功环原始电压参考，线电压RMS，V
cfg.control.Eref = sqrt(2/3)*cfg.control.UrefLL;

cfg.control.J = 1;                % 有功惯量
cfg.control.Komega = 16805;         % 频率偏差到功率的增益，保留你的当前设置
cfg.control.Ddirect = 52;           % 转矩通道直接阻尼，保留你的当前设置
cfg.control.Dp = cfg.control.Ddirect ...
    + cfg.control.Komega/cfg.base.w0;
% 附件实际方程：
%   J*d(domega)/dt = (Pref-Pf-Komega*domega)/(w0+domega)
%                    - Ddirect*domega。
% 当稳态omega=w0且Pf=Pref时，局部线性化与论文完全等效，
%   Dp = Ddirect + Komega/w0；随上面的原始参数自动重新计算。

%% 5. 无功环：线电压RMS到相电压峰值的单位换算

cfg.control.KqIntegralLL = 0.0128;  % 无功积分增益；保留你修改后的值
cfg.control.DqLL = 3000;           % 无功电压反馈增益；var/V(line RMS)

cfg.control.bq = sqrt(2/3)*cfg.control.KqIntegralLL;
cfg.control.K = 1/cfg.control.bq;    % 论文记法，dE/dt=(...)/K
cfg.control.Dq = sqrt(3/2)*cfg.control.DqLL;
% 因此：
%   dE/dt = bq*[Qref-Qf+Dq*(Eref-Ufeedback_peak)]。
% paper15：Ufeedback_peak=vpd。
% model17：Ufeedback_peak=sqrt(vpd^2+vpq^2)，近似附件RMS幅值反馈的基波。
% 附件RMS滑窗的动态没有被虚构成论文的功率低通滤波器。

%% 6. Pe、Qe输出处的功率滤波器

cfg.measurement.wp = 2*pi*10;       % Pe对应的一阶低通截止角频率，rad/s
cfg.measurement.wq = 2*pi*10;       % Qe对应的一阶低通截止角频率，rad/s
% 对应附件 坐标变换1/Subsystem2 与 /Subsystem1：
%   dPf/dt = wp*(Pe-Pf)，dQf/dt = wq*(Qe-Qf)。
% 上面2*pi*10表示截止频率10Hz，方程中使用的是62.8319rad/s。
% 瞬时功率使用PCC电压和滤波电感电流io，和论文式(1)一致。

%% 7. 虚拟阻抗与双闭环

%cfg.virtual.Rv = 0.04398;           % 虚拟电阻，ohm
%cfg.virtual.Lv = 1.4e-3;            % 虚拟电感，H
%cfg.virtual.wFilter = 2*pi*500;     % 虚拟阻抗电流预滤波，rad/s，model17使用

cfg.virtual.Rv = 0;           % 虚拟电阻，ohm
cfg.virtual.Lv = 0;            % 虚拟电感，H
cfg.virtual.wFilter = 2*pi*500;     % 虚拟阻抗电流预滤波，rad/s，model17使用

cfg.inner.kvp = 0.16;              % 电压环比例增益
cfg.inner.kvi = 500;               % 电压环积分增益
cfg.inner.kip = 1.2;              % 电流环比例增益
cfg.inner.kii = 526*2;             % 电流环积分增益=1052

%% 8. 表3的模态选择与参与因子定义

% auto：在指定频率区间，选择实部最大的正虚部特征值。
% manual：按照表2式排序后的正虚部编号，手动选择。
% 不把“第8、9个特征值”硬编码为附件模型的关键模态。
cfg.participation.selection = 'auto';
cfg.participation.frequencyHz = [0.05 50];
cfg.participation.manualIndex = 8;
% 表3同时输出p、abs(p)和abs(p)/sum(abs(p))三列。
% 论文表3的数值合计大于1，不应当不加说明地当成百分比。

%% 9. 可选：有功频率环与delta参与度实验

% 运行实验：
%   cfg = state_space_parameters();
%   cfg.activeExperiment.enabled = true;
%   result = run_state_space_repro(cfg);
% 普通复现仍按原流程执行；本实验不会自动覆盖基准控制参数或修改SLX。
cfg.activeExperiment.enabled = false;
cfg.activeExperiment.SCR = cfg.experiment.participationSCR;
cfg.activeExperiment.JFactors = [0.25 0.5 0.75 1 1.5 2 3 4];
cfg.activeExperiment.KomegaFactors = [0.25 0.5 1 1.5 2];
cfg.activeExperiment.DdirectFactors = [0 0.5 1 1.5 2 3 5];
% 补充保持Komega不变、联合改变J与等效阻尼的实验。
cfg.activeExperiment.DeqFactors = [0.6 0.7 0.75 0.8 0.85 0.9 0.95 1 1.1 1.25 1.5 2 3];
cfg.activeExperiment.continuationSteps = 12;
cfg.activeExperiment.maxTrackingStep = 0.75;
cfg.activeExperiment.validationSCR = [2 3 5 10 20 30 40];
cfg.activeExperiment.pdfName = '有功频率环_delta参与度实验报告.pdf';
% 同时跟踪共轭根对应的二维子空间：跨越复根合并时，不换另一对模态。
% 候选需满足：delta原始模值与占比下降、复根保留、阻尼比不降、
% 目标根实部不右移、全部根稳定，并通过validationSCR的全部根检查。

%% 10. 数值精度、绘图范围与文件输出

cfg.numeric.newtonTolerance = 1e-10;
cfg.numeric.maxIterations = 80;
cfg.numeric.complexStep = 1e-20;
cfg.numeric.finiteDifferenceStep = 1e-5;
cfg.numeric.stabilityTolerance = 1e-7;
cfg.numeric.conjugateTolerance = 1e-6;
cfg.numeric.jacobianTolerance = 1e-6;
cfg.numeric.localCheckAmplitude = 1e-5;

cfg.plot.visible = 'on';
cfg.plot.fontName = 'Microsoft YaHei';
cfg.plot.markerSize = 12;
cfg.plot.zoomReal = [-180 80];       % 图7低频局部图的实轴范围
cfg.plot.zoomImag = [-400 400];     % 图7低频局部图的虚轴范围，rad/s
cfg.plot.participationTopN = 8;

% 默认按主模型生成一份PDF；图和表留在MATLAB窗口。
% 打开compareStructures时，对照结构的完整图表也放进同一PDF。
% 不另存PNG/FIG/CSV/MAT；本程序不会修改SLX模型。
cfg.output.folder = fileparts(mfilename('fullpath'));
% 旧版PDF正在被外部阅读器占用，本次修订使用新文件名。
cfg.output.pdfName = '状态空间模型_图7与表2表3复现报告_修订版.pdf';
cfg.output.writePDF = true;
cfg.output.keepFigures = true;
end
