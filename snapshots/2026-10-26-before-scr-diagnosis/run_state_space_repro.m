function result = run_state_space_repro(cfg)
%RUN_STATE_SPACE_REPRO 按用户参数复现论文状态空间方法、图7、表2及表3。
%
% 直接运行：result = run_state_space_repro();
% 自定义配置：cfg = state_space_parameters(); 修改cfg后再调用本函数。
%
% 本文件包含所有运算、绘图、表格与PDF输出函数；无需Simulink或控制工具箱。
% 不读取/修改原模型，不运行开关仿真，不另存任何过程文件。
% 默认只绘制model17；paper15仅在参数文件中显式选择或开启对照时输出。
% result保存主模型的A、B、C、D、工作点、扫描、表格、参与因子和核验结果。
% 开启完整对照后，另一结构的同类结果保存在result.comparison中。
%
% 建模依据：论文式(4)-(12)。从非线性方程求Jacobian，而不是照抄附录A：
% 附录A的A(3,11)漏列了电流环比例通道kip*kvp/Lf，本文恢复该项。
% 使用相电压/相电流峰值dq，P=1.5*(vd*id+vq*iq)，Q=1.5*(vq*id-vd*iq)。

%% 1. 读取配置，并由原始设置重新计算派生参数

if nargin < 1
    cfg = state_space_parameters();
end
cfg = refresh_derived_parameters(cfg);
validate_configuration(cfg);

fprintf('\n=== 主模型：%s ===\n', cfg.model.structure);
fprintf('Pref=%.6g W；Qref=%.6g var；电压PI积分增益=%.6g\n', ...
    cfg.control.Pref, cfg.control.Qref, cfg.inner.kvi);
fprintf('结构名将同时显示在图形窗口、图题与PDF中。\n');

if isfield(cfg, 'activeExperiment') && cfg.activeExperiment.enabled
    % 实验模式只改变配置副本，跳过普通图7，以免重复打开两套图形。
    result = run_active_frequency_experiment(cfg);
    figures = make_active_experiment_figures(result, cfg);
    cfg.output.pdfName = cfg.activeExperiment.pdfName;
else
%% 2. 主模型：完整扫描、工作点、状态矩阵与两张表

result = calculate_experiment(cfg);
result.comparison = [];

%% 3. 可选对照：另一结构也执行同样的完整实验

if cfg.experiment.compareStructures
    cfgOther = cfg;
    cfgOther.experiment.compareStructures = false;
    if strcmp(cfg.model.structure, 'model17')
        cfgOther.model.structure = 'paper15';
    else
        cfgOther.model.structure = 'model17';
    end
    fprintf('\n=== 完整对照模型：%s（与主模型采用相同原始参数） ===\n', ...
        cfgOther.model.structure);
    result.comparison = calculate_experiment(cfgOther);
end

%% 4. 绘制主模型的图7、表2、表3和参与因子图

figures = make_experiment_figures(result, cfg);
if ~isempty(result.comparison)
    comparisonFigures = make_experiment_figures( ...
        result.comparison, result.comparison.cfg);
    figures = [figures; comparisonFigures; plot_structure_comparison(result, cfg)];
end
end

%% 5. 全部选中结构的图表放在同一个PDF中

if cfg.output.writePDF
    assert(isfolder(cfg.output.folder), '输出目录不存在。');
    pdfPath = fullfile(cfg.output.folder, cfg.output.pdfName);

    for k = 1:numel(figures)
        if k == 1
            exportgraphics(figures(k), pdfPath, ...
                'ContentType', 'vector', 'BackgroundColor', 'white');
        else
            exportgraphics(figures(k), pdfPath, 'Append', true, ...
                'ContentType', 'vector', 'BackgroundColor', 'white');
        end
    end
    fprintf('\nPDF报告：%s\n', pdfPath);
end

result.figures = figures;
if ~cfg.output.keepFigures
    close(figures);
end
fprintf('完成：主模型为%s；原SLX与基准控制参数未被修改。\n', cfg.model.structure);
fprintf('未生成MAT、CSV、PNG、FIG等过程性文件。\n');
end

function c = refresh_derived_parameters(c)
% 允许cfg=state_space_parameters()后再修改原始参数。
% 例如修改Komega后不能继续使用旧的Dp；所有派生量在此统一刷新。
c.base.w0 = 2*pi*c.base.f0;
c.base.Eg = sqrt(2/3)*c.base.Ull;
c.base.Zb = c.base.Ull^2/c.base.Sn;
c.base.Ipk = sqrt(2)*c.base.Sn/(sqrt(3)*c.base.Ull);
c.plant.Ts = 1/c.plant.fs;
c.control.Eref = sqrt(2/3)*c.control.UrefLL;
c.control.Dp = c.control.Ddirect+c.control.Komega/c.base.w0;
c.control.bq = sqrt(2/3)*c.control.KqIntegralLL;
c.control.K = 1/c.control.bq;
c.control.Dq = sqrt(3/2)*c.control.DqLL;
end

function r = calculate_experiment(c)
% 主模型与对照模型共用这一个流程，避免“只算了、没有画”的不一致。
fprintf('\n首次实验：SCR=%.6g，结构=%s\n', ...
    c.experiment.firstSCR, c.model.structure);
r = struct();
r.cfg = c;
r.first = analyse_case(c, c.experiment.firstSCR, []);
print_case(r.first);

scr = unique(c.experiment.figure7SCR(:).');
nCase = numel(scr);
r.scan = cell(1, nCase);
seed = [];

for k = 1:nCase
    r.scan{k} = analyse_case(c, scr(k), seed);
    seed = r.scan{k}.equilibrium.z;
end

r.SCR = scr;
r.spectralAbscissa = cellfun(@(a) max(real(a.lambda)), r.scan);
r.criticalSCR = refine_boundaries(c, r.scan);
r.tableCase = analyse_case(c, c.experiment.tableSCR, []);
r.participationCase = analyse_case(c, c.experiment.participationSCR, []);
r.table2 = make_eigenvalue_table(r.tableCase, c);
[r.table3, r.selectedMode] = make_participation_table(r.participationCase, c);

% 指定的次同步频段不一定包含整个系统最危险的模态。
% 另算全局最大实部模态，防止“低频模态稳定”被误读为“系统稳定”。
cfgGlobal = c;
cfgGlobal.participation.selection = 'manual';
[~, cfgGlobal.participation.manualIndex] = ...
    max(real(r.participationCase.lambda));
[r.globalTable3, r.globalMode] = make_participation_table( ...
    r.participationCase, cfgGlobal);

fprintf('\n表2：SCR = %.6g，%s\n', ...
    c.experiment.tableSCR, c.model.structure);
disp(r.table2);

fprintf('\n表3：SCR = %.6g，目标特征值 = %.8g %+.8gj\n', ...
    c.experiment.participationSCR, real(r.selectedMode.lambda), ...
    imag(r.selectedMode.lambda));
disp(r.table3);
fprintf('全局最大实部模态：%.8g %+.8gj；其参与因子见result.globalTable3。\n', ...
    real(r.globalMode.lambda), imag(r.globalMode.lambda));
fprintf('扫描点数=%d；最大实部范围=[%.8g, %.8g] 1/s。\n', ...
    nCase, min(r.spectralAbscissa), max(r.spectralAbscissa));
end

function figures = make_experiment_figures(r, c)
figures = gobjects(6, 1);
figures(1) = make_summary_page(r, c);
figures(2) = plot_figure7(r, c);
figures(3) = plot_eigenvalue_table(r, c);
figures(4) = plot_participation_table(r, c);
figures(5) = plot_participation_bars(r, c);
figures(6) = make_validation_page(r, c);
end

%% 参数检查：不将无可行工作点误判成动态失稳
function validate_configuration(c)
validatestring(c.model.structure, {'paper15', 'model17'});
assert(c.base.Sn > 0 && c.base.Ull > 0, '额定基准必须为正。');
assert(c.grid.XR > 0, 'X/R必须为正。');
assert(all(c.experiment.figure7SCR > 0), 'SCR扫描值必须为正。');
assert(c.experiment.firstSCR > 0 && c.experiment.tableSCR > 0 ...
    && c.experiment.participationSCR > 0, '工况SCR必须为正。');
assert(c.plant.Lf > 0 && c.plant.Cf > 0, 'LC参数必须为正。');
assert(c.control.J > 0 && c.control.bq > 0, '惯量与无功积分增益必须为正。');
assert(c.inner.kvi > 0 && c.inner.kii > 0, 'PI积分增益必须为正。');
assert(~isempty(c.experiment.figure7SCR), '图7扫描范围不能为空。');
end

%% 单工况分析：稳态 -> A/B/C/D -> 特征值 -> 左右特征向量
function a = analyse_case(c, scr, seed)
p = make_case_parameters(c, scr);
eq = solve_equilibrium(p, c, seed);
x0 = eq.x;
u0 = p.u0;

% 复步长法没有中心差分的相减消去误差。
% nonlinear_rhs内部必须保持解析运算，不可把abs、real、conj加入状态计算。
A = complex_jacobian(@(x) nonlinear_rhs(x, u0, p), ...
    x0, c.numeric.complexStep);
B = complex_jacobian(@(u) nonlinear_rhs(x0, u, p), ...
    u0, c.numeric.complexStep);
C = complex_jacobian(@(x) measured_outputs(x, u0, p), ...
    x0, c.numeric.complexStep);
D = complex_jacobian(@(u) measured_outputs(x0, u, p), ...
    u0, c.numeric.complexStep);

% 状态缩放只改善数值条件，不改变特征值与参与因子。
s = state_scales(p);
As = (A .* s.') ./ s;
[Vs, L] = eig(As);
lambda = diag(L);
order = order_modes(lambda, c.numeric.conjugateTolerance);
lambda = lambda(order);
Vs = Vs(:, order);

% W为左特征向量按行组成的矩阵：W*V=I。
% 这里不使用MATLAB的共轭转置'，参与因子需要代数转置.'。
Ws = Vs \ eye(numel(lambda));
V = s .* Vs;
W = Ws ./ s.';
P = Vs .* Ws.';

check = validate_linearization(x0, u0, p, A, B, Vs, Ws, lambda, c);
assert(check.jacobianRelativeError < c.numeric.jacobianTolerance, ...
    '状态Jacobian交叉核验未通过。');
assert(check.inputJacobianRelativeError < c.numeric.jacobianTolerance, ...
    '输入Jacobian交叉核验未通过。');

a = struct('SCR', scr, 'parameters', p, 'equilibrium', eq, ...
    'A', A, 'B', B, 'C', C, 'D', D, 'lambda', lambda, ...
    'rightEigenvectors', V, 'leftEigenvectors', W, ...
    'participation', P, 'validation', check);
a.outputs = measured_outputs(x0, u0, p);
a.stateNames = state_names(p);
end

function p = make_case_parameters(c, scr)
p = struct();
p.c = c;
p.structure = c.model.structure;
p.w0 = c.base.w0;
p.Lf = c.plant.Lf;
p.Cf = c.plant.Cf;
p.Rf = c.plant.Rf;
Zg = c.base.Ull^2/(c.base.Sn*scr);
p.Rg = Zg/sqrt(1+c.grid.XR^2);
p.Lg = p.Rg*c.grid.XR/p.w0;
p.Rc = 0;
p.voltageFeedforward = 1;
p.useMagnitude = false;
p.useVirtualFilter = false;

if strcmp(p.structure, 'model17')
    p.Rc = c.plant.Rc;
    p.voltageFeedforward = 0;
    p.useMagnitude = true;
    p.useVirtualFilter = true;
end

% 输入顺序与论文一致，但Eg与Eref均统一为相峰值。
p.u0 = [c.base.Eg; c.control.Pref; p.w0; ...
    c.control.Qref; c.control.Eref];
end

%% 工作点的六个独立代数条件
function r = equilibrium_residual(z, p)
% z=[vpd,vpq,igd,igq,delta,E]，全部为物理量。
v = z(1:2);
ig = z(3:4);
delta = z(5);
E = z(6);
J2 = [0 -1; 1 0];
vc = (eye(2)+p.Rc*p.w0*p.Cf*J2) \ v;
ic = p.w0*p.Cf*J2*vc;
io = ig+ic;
eg = p.u0(1)*[cos(delta); -sin(delta)];

ev = [E; 0] - p.c.virtual.Rv*io ...
    - p.w0*p.c.virtual.Lv*J2*io;
Pe = 1.5*(v(1)*io(1)+v(2)*io(2));
Qe = 1.5*(v(2)*io(1)-v(1)*io(2));

if p.useMagnitude
    Ufeedback = sqrt(v(1)^2+v(2)^2);
else
    Ufeedback = v(1);
end

r = [v-eg-p.Rg*ig-p.w0*p.Lg*J2*ig; ...
    Pe-p.u0(2); ...
    p.u0(4)-Qe+p.c.control.Dq*(p.u0(5)-Ufeedback); ...
    ev-v];
end

function eq = solve_equilibrium(p, c, seed)
% 只解六个独立条件，避免将15/17个混合单位状态一起交给黑箱求解器。
if isempty(seed)
    seed = [c.base.Eg; 0; 0; 0; 0; c.control.Eref];
end
zs = [c.base.Eg; c.base.Eg; c.base.Ipk; c.base.Ipk; 1; c.base.Eg];
rs = [c.base.Eg; c.base.Eg; c.base.Sn; c.base.Sn; ...
    c.base.Eg; c.base.Eg];
y = seed./zs;
fun = @(v) equilibrium_residual(v.*zs, p)./rs;

for k = 1:c.numeric.maxIterations
    residual = fun(y);
    if norm(residual, inf) < c.numeric.newtonTolerance
        break;
    end

    H = complex_jacobian(fun, y, c.numeric.complexStep);
    assert(rcond(H) > 1e-14, '工作点Jacobian近似奇异；请检查稳态可行性。');
    step = -H \ residual;
    alpha = 1;
    accepted = false;

    for j = 1:24
        trial = y+alpha*step;
        if trial(1) > 0 && trial(6) > 0 ...
                && norm(fun(trial), 2) < norm(residual, 2)
            y = trial;
            accepted = true;
            break;
        end
        alpha = alpha/2;
    end
    assert(accepted, '稳态求解线搜索失败，不输出虚假的稳定性结论。');
end

assert(norm(fun(y), inf) < c.numeric.newtonTolerance, ...
    '当前SCR下稳态求解未收敛。');
z = y.*zs;
v = z(1:2);
ig = z(3:4);
J2 = [0 -1; 1 0];
vc = (eye(2)+p.Rc*p.w0*p.Cf*J2) \ v;
io = ig+p.w0*p.Cf*J2*vc;
P0 = 1.5*(v(1)*io(1)+v(2)*io(2));
Q0 = 1.5*(v(2)*io(1)-v(1)*io(2));

% 电压环误差为0；积分器仍可有非零稳态值。
phi = (io-p.w0*p.Cf*J2*v)/c.inner.kvi;
gamma = ((1-p.voltageFeedforward)*v+p.Rf*io)/c.inner.kii;
x = [ig; io; vc; P0; Q0; 0; z(5); z(6); phi; gamma];
if p.useVirtualFilter
    x = [x; io];
end

dx = nonlinear_rhs(x, p.u0, p);
scaledResidual = norm(dx./state_scales(p), inf);
assert(scaledResidual < 1e-6, '完整状态方程的稳态残差过大。');
eq = struct('z', z, 'x', x, 'u', p.u0, ...
    'algebraicResidual', norm(fun(y), inf), ...
    'dynamicResidual', scaledResidual, 'iterations', k);
end

%% 非线性连续平均状态方程：论文式(4)-(12)的直接实现
function dx = nonlinear_rhs(x, u, p)
% 15阶状态顺序：
% [igd igq iod ioq vcd vcq Pf Qf domega delta E phi_d phi_q gamma_d gamma_q]'
% paper15中Rc=0，因此vc=vp，与论文完全相同。
% model17中再附加[ivd ivq]'两个虚拟阻抗滤波状态。
ig = x(1:2);
io = x(3:4);
vc = x(5:6);
Pf = x(7);
Qf = x(8);
domega = x(9);
delta = x(10);
E = x(11);
phi = x(12:13);
gamma = x(14:15);
w = p.w0+domega;
J2 = [0 -1; 1 0];

% 电容ESR使PCC电压与理想电容电压不同，不增加储能状态。
ic = io-ig;
vp = vc+p.Rc*ic;
eg = u(1)*[cos(delta); -sin(delta)];

if p.useVirtualFilter
    iv = x(16:17);
else
    iv = io;
end

% 虚拟阻抗：频率交叉项使用模型控制器中的额定w0。
ev = [E; 0] - p.c.virtual.Rv*iv ...
    - p.w0*p.c.virtual.Lv*J2*iv;
voltageError = ev-vp;

% 电压PI、电容电流交叉补偿、电流PI与电感交叉解耦。
iref = p.c.inner.kvp*voltageError+p.c.inner.kvi*phi ...
    + p.w0*p.Cf*J2*vp;
currentError = iref-io;
eo = p.c.inner.kip*currentError+p.c.inner.kii*gamma ...
    + p.w0*p.Lf*J2*io+p.voltageFeedforward*vp;

% 功率定义：PCC电压与变换器侧io；不能悄悄替换成网侧ig。
Pe = 1.5*(vp(1)*io(1)+vp(2)*io(2));
Qe = 1.5*(vp(2)*io(1)-vp(1)*io(2));

if p.useMagnitude
    Ufeedback = sqrt(vp(1)^2+vp(2)^2);
else
    Ufeedback = vp(1);
end

dx = zeros(size(x), 'like', x+u(1));
dx(1:2) = (vp-eg-p.Rg*ig)/p.Lg-w*J2*ig;
dx(3:4) = (eo-vp-p.Rf*io)/p.Lf-w*J2*io;
dx(5:6) = ic/p.Cf-w*J2*vc;
dx(7) = p.c.measurement.wp*(Pe-Pf);
dx(8) = p.c.measurement.wq*(Qe-Qf);

if p.useVirtualFilter
    % 附件的实际转矩方程；与论文在同步稳态附近一阶等效。
    dx(9) = ((u(2)-Pf-p.c.control.Komega*domega)/w ...
        - p.c.control.Ddirect*domega)/p.c.control.J;
else
    dx(9) = (u(2)-Pf)/(p.c.control.J*p.w0) ...
        - p.c.control.Dp*domega/p.c.control.J;
end

dx(10) = w-u(3);
dx(11) = p.c.control.bq*(u(4)-Qf ...
    + p.c.control.Dq*(u(5)-Ufeedback));
dx(12:13) = voltageError;
dx(14:15) = currentError;

if p.useVirtualFilter
    dx(16:17) = p.c.virtual.wFilter*(io-iv);
end
end

function y = measured_outputs(x, ~, p)
v = x(5:6)+p.Rc*(x(3:4)-x(1:2));
io = x(3:4);
y = [1.5*(v(1)*io(1)+v(2)*io(2)); ...
    1.5*(v(2)*io(1)-v(1)*io(2)); ...
    sqrt(3/2)*sqrt(v(1)^2+v(2)^2); ...
    p.w0+x(9); x(10)];
% 输出：[Pe(W), Qe(var), U_LL(RMS,V), omega(rad/s), delta(rad)]。
end

function J = complex_jacobian(fun, x, h)
y = fun(x);
J = zeros(numel(y), numel(x));
for k = 1:numel(x)
    trial = complex(x);
    trial(k) = trial(k)+1i*h;
    J(:, k) = imag(fun(trial))/h;
end
end

function s = state_scales(p)
c = p.c;
V = c.base.Eg;
I = c.base.Ipk;
s = [I; I; I; I; V; V; c.base.Sn; c.base.Sn; p.w0; 1; V; ...
    I/c.inner.kvi; I/c.inner.kvi; V/c.inner.kii; V/c.inner.kii];
if p.useVirtualFilter
    s = [s; I; I];
end
end

function names = state_names(p)
names = {'i_gd'; 'i_gq'; 'i_od'; 'i_oq'; 'v_pd'; 'v_pq'; ...
    'P_f'; 'Q_f'; 'omega'; 'delta'; 'E'; ...
    'phi_d'; 'phi_q'; 'gamma_d'; 'gamma_q'};
if p.useVirtualFilter
    names(5:6) = {'v_cd'; 'v_cq'};
    names = [names; {'i_vd'; 'i_vq'}];
end
end

%% 数值核验：独立差分、左右向量、共轭配对和局部线性误差
function check = validate_linearization(x, u, p, A, B, V, W, lambda, c)
s = state_scales(p);
As = (A.*s.')./s;
ys = x./s;
h = c.numeric.finiteDifferenceStep;
fd = zeros(size(A));
for k = 1:numel(x)
    e = zeros(size(x));
    e(k) = h;
    fp = nonlinear_rhs((ys+e).*s, u, p)./s;
    fm = nonlinear_rhs((ys-e).*s, u, p)./s;
    fd(:, k) = (fp-fm)/(2*h);
end
check.jacobianRelativeError = norm(As-fd, 'fro')/max(1, norm(As, 'fro'));

us = [c.base.Eg; c.base.Sn; p.w0; c.base.Sn; c.base.Eg];
Bscaled = (B.*us.')./s;
fdB = zeros(size(B));
for k = 1:numel(u)
    e = zeros(size(u));
    e(k) = h*us(k);
    fdB(:, k) = (nonlinear_rhs(x, u+e, p) ...
        - nonlinear_rhs(x, u-e, p))./(2*h*s);
end
check.inputJacobianRelativeError = ...
    norm(Bscaled-fdB, 'fro')/max(1, norm(Bscaled, 'fro'));
check.eigenResidual = norm(As*V-V*diag(lambda), 'fro') ...
    /max(1, norm(As, 'fro')*norm(V, 'fro'));
check.dualityError = norm(W*V-eye(size(A)), 'fro');
check.participationSumError = max(abs(sum(V.*W.', 1)-1));
check.eigenvectorCondition = cond(V);
check.conjugateError = max(arrayfun(@(z) min(abs(lambda-conj(z))), lambda));

% 相同方向上幅值减半，非线性余项应当近似缩小到1/4。
direction = sin((1:numel(x)).');
direction = direction/norm(direction);
amp = c.numeric.localCheckAmplitude;
error = zeros(2, 1);
f0 = nonlinear_rhs(x, u, p);
for k = 1:2
    d = s.*direction*amp/(2^(k-1));
    error(k) = norm((nonlinear_rhs(x+d, u, p)-f0-A*d)./s);
end
check.localRemainder = error;
check.localRemainderRatio = error(1)/max(realmin, error(2));
check.spectralAbscissa = max(real(lambda));
check.unstableModeCount = sum(real(lambda) > c.numeric.stabilityTolerance);
end

function order = order_modes(z, tol)
% 仅按频率与共轭关系排序；编号没有跨模型的物理不变性。
positive = find(imag(z) > tol);
[~, s] = sort(imag(z(positive)), 'descend');
positive = positive(s);
order = zeros(0, 1);
used = false(size(z));
for k = positive(:).'
    if used(k)
        continue;
    end
    order(end+1, 1) = k; %#ok<AGROW>
    used(k) = true;
    distance = abs(z-conj(z(k)));
    distance(used) = inf;
    [~, j] = min(distance);
    if distance(j) < tol*max(1, abs(z(k)))
        order(end+1, 1) = j; %#ok<AGROW>
        used(j) = true;
    end
end
remaining = find(~used);
[~, s] = sort(real(z(remaining)), 'ascend');
order = [order; remaining(s)];
end

function T = make_eigenvalue_table(a, c)
z = a.lambda;
Mode = (1:numel(z)).';
RealPart = real(z);
ImagPart = imag(z);
FrequencyHz = abs(imag(z))/(2*pi);
DampingRatio = -real(z)./max(abs(z), realmin);
Stable = real(z) < -c.numeric.stabilityTolerance;
T = table(Mode, RealPart, ImagPart, FrequencyHz, DampingRatio, Stable);
end

function [T, selected] = make_participation_table(a, c)
z = a.lambda;
if strcmp(c.participation.selection, 'manual')
    index = c.participation.manualIndex;
    assert(index >= 1 && index <= numel(z), '手动模态编号超出范围。');
else
    f = imag(z)/(2*pi);
    band = c.participation.frequencyHz;
    candidates = find(f >= band(1) & f <= band(2));
    assert(~isempty(candidates), ...
        '指定频段无正虚部振荡模态，请修改参数文件的模态频率区间。');
    [~, k] = max(real(z(candidates)));
    index = candidates(k);
end
raw = a.participation(:, index);
mag = abs(raw);
normalized = mag/sum(mag);
State = a.stateNames;
RealP = real(raw);
ImagP = imag(raw);
AbsP = mag;
Normalized = normalized;
T = table(State, RealP, ImagP, AbsP, Normalized);
selected = struct('index', index, 'lambda', z(index), ...
    'frequencyHz', abs(imag(z(index)))/(2*pi), ...
    'sumRaw', sum(raw), 'sumAbs', sum(mag));
end

function boundary = refine_boundaries(c, cases)
alpha = cellfun(@(a) max(real(a.lambda)), cases);
boundary = [];
for k = 1:numel(cases)-1
    if alpha(k)*alpha(k+1) >= 0
        continue;
    end
    lo = cases{k}.SCR;
    hi = cases{k+1}.SCR;
    aLo = alpha(k);
    for j = 1:30
        mid = (lo+hi)/2;
        current = analyse_case(c, mid, cases{k}.equilibrium.z);
        aMid = max(real(current.lambda));
        if aLo*aMid <= 0
            hi = mid;
        else
            lo = mid;
            aLo = aMid;
        end
    end
    boundary(end+1) = (lo+hi)/2; %#ok<AGROW>
end
end

function print_case(a)
fprintf('阶数=%d，Rg=%.9g ohm，Lg=%.9g H\n', ...
    size(a.A, 1), a.parameters.Rg, a.parameters.Lg);
fprintf('Pe=%.9g W，Qe=%.9g var，U_LL=%.9g V，delta=%.9g rad\n', ...
    a.outputs(1), a.outputs(2), a.outputs(3), a.outputs(5));
fprintf('最大特征根实部=%.9g 1/s，不稳定根数=%d\n', ...
    a.validation.spectralAbscissa, a.validation.unstableModeCount);
fprintf('工作点归一化动态残差=%.3e，Jacobian交叉误差=%.3e\n', ...
    a.equilibrium.dynamicResidual, a.validation.jacobianRelativeError);
end

%% 图表与报告页面：采用MATLAB原生图形，不依赖第三方报告工具
function f = new_page(c, name)
name = sprintf('[%s] %s', c.model.structure, name);
f = figure('Name', name, 'Color', 'white', 'NumberTitle', 'off', ...
    'Visible', c.plot.visible, 'Position', [80 60 1050 760]);
set(f, 'DefaultAxesFontName', c.plot.fontName, ...
    'DefaultTextFontName', c.plot.fontName);
end

function f = make_summary_page(r, c)
f = new_page(c, '模型与参数说明');
ax = axes(f, 'Position', [0.07 0.06 0.86 0.86]);
axis(ax, 'off');
title(ax, '状态空间模型与图7、表2、表3复现', 'FontSize', 20);
lines = {
    sprintf('当前结构：%s；首次SCR=%.4g；Pref=%.4g W；Qref=%.4g var', ...
        c.model.structure, c.experiment.firstSCR, c.control.Pref, c.control.Qref)
    '依据：刘永奎等，电力系统及其自动化学报，2026，38(2)：126-134。'
    '从正文式(4)-(12)建立非线性方程并求导；不硬填论文的特征值或临界SCR。'
    sprintf('U_LL=%.4g V，S_N=%.4g kVA，f0=%.4g Hz，X/R=%.4g', ...
        c.base.Ull, c.base.Sn/1e3, c.base.f0, c.grid.XR)
    sprintf('Lf=%.4g uH，Cf=%.4g uF，Rf=%.4g ohm，Rc=%.4g ohm(model17)', ...
        c.plant.Lf*1e6, c.plant.Cf*1e6, c.plant.Rf, c.plant.Rc)
    sprintf('J=%.5g；Dp=%.6g+%.6g/w0=%.8g；K=%.8g；Dq=%.8g', ...
        c.control.J, c.control.Ddirect, c.control.Komega, ...
        c.control.Dp, c.control.K, c.control.Dq)
    sprintf('Pe/Qe低通：wp=%.6g，wq=%.6g rad/s；截止频率%.6g / %.6g Hz。', ...
        c.measurement.wp, c.measurement.wq, ...
        c.measurement.wp/(2*pi), c.measurement.wq/(2*pi))
    sprintf('虚拟阻抗：Rv=%.6g ohm，Lv=%.6g mH；虚拟电流滤波%.6g Hz。', ...
        c.virtual.Rv, c.virtual.Lv*1e3, c.virtual.wFilter/(2*pi))
    sprintf('电压PI=(%.4g, %.4g)；电流PI=(%.4g, %.4g)', ...
        c.inner.kvp, c.inner.kvi, c.inner.kip, c.inner.kii)
    'paper15：保留论文电压前馈、v_pd无功反馈；忽略Rc和虚拟电流滤波。'
    'model17：采用附件无电压前馈、幅值反馈、Rc与两个虚拟电流滤波状态。'
    '两者均忽略PWM/采样延时、RMS滑窗、直流动态、限幅、预同步与保护切换。'
    '结果是局部小信号结论，不能直接当作完整开关模型的暂态稳定判定。'
    sprintf('本次工作点：Pe=%.7g W，Qe=%.7g var，U_LL=%.7g V。', ...
        r.first.outputs(1), r.first.outputs(2), r.first.outputs(3))
    '实部不为零是正常现象；任一特征根实部为正，才表示该线性模型失稳。'
    '参数只修改state_space_parameters.m；计算、图表与PDF全部位于主程序。'
    };
text(ax, 0, 0.98, lines, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 11, ...
    'Interpreter', 'none');
end

function f = plot_figure7(r, c)
f = new_page(c, '图7：SCR变化时的全部特征根');
layout = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('图7：%s，SCR = %.4g ~ %.4g', ...
    c.model.structure, min(r.SCR), max(r.SCR)), 'Interpreter', 'none');
n = numel(r.scan{1}.lambda);
Z = cell2mat(cellfun(@(a) a.lambda, r.scan, 'UniformOutput', false));
color = repmat(r.SCR, n, 1);
ax1 = nexttile(layout, [1 2]);
scatter(ax1, real(Z(:)), imag(Z(:)), c.plot.markerSize, color(:), 'filled');
hold(ax1, 'on');
xline(ax1, 0, '--k');
plot(ax1, real(r.first.lambda), imag(r.first.lambda), 'ro', ...
    'MarkerSize', 6, 'LineWidth', 1);
xlabel(ax1, '实部 / s^{-1}'); ylabel(ax1, '虚部 / rad s^{-1}');
grid(ax1, 'on'); colorbar(ax1); colormap(f, parula);
title(ax1, sprintf('全局分布；红圈为首次SCR=%.4g', c.experiment.firstSCR));

ax2 = nexttile(layout);
scatter(ax2, real(Z(:)), imag(Z(:)), c.plot.markerSize, color(:), 'filled');
xline(ax2, 0, '--k');
xlim(ax2, c.plot.zoomReal); ylim(ax2, c.plot.zoomImag);
xlabel(ax2, '实部 / s^{-1}'); ylabel(ax2, '虚部 / rad s^{-1}');
title(ax2, '低频模态局部放大'); grid(ax2, 'on');

ax3 = nexttile(layout);
plot(ax3, r.SCR, r.spectralAbscissa, 'LineWidth', 1.5);
yline(ax3, 0, '--k'); grid(ax3, 'on');
xlabel(ax3, 'SCR'); ylabel(ax3, 'max Re(lambda) / s^{-1}');
title(ax3, '全部模态的稳定性边界');
for value = r.criticalSCR
    xline(ax3, value, ':r', sprintf('%.4g', value));
end
end

function f = plot_eigenvalue_table(r, c)
f = new_page(c, '表2：状态矩阵的全部特征值');
ax = axes(f, 'Position', [0.06 0.05 0.88 0.84]);
axis(ax, 'off');
title(ax, sprintf('表2：%s，SCR=%.5g', ...
    c.model.structure, c.experiment.tableSCR), 'Interpreter', 'none');
T = r.table2;
header = {'编号', '特征值 / s^{-1}', '频率 / Hz', '阻尼比', '实部判定'};
rows = cell(height(T), 5);
for k = 1:height(T)
    rows{k, 1} = sprintf('%d', k);
    rows{k, 2} = sprintf('%.7g %+.7gj', T.RealPart(k), T.ImagPart(k));
    rows{k, 3} = sprintf('%.7g', T.FrequencyHz(k));
    rows{k, 4} = sprintf('%.7g', T.DampingRatio(k));
    if T.RealPart(k) > c.numeric.stabilityTolerance
        rows{k, 5} = '不稳定';
    elseif T.RealPart(k) < -c.numeric.stabilityTolerance
        rows{k, 5} = '稳定';
    else
        rows{k, 5} = '近虚轴，需进一步判断';
    end
end
draw_table(ax, header, rows, [0 0.13 0.56 0.71 0.85], c);
end

function f = plot_participation_table(r, c)
f = new_page(c, '表3：所选模态的参与因子');
ax = axes(f, 'Position', [0.06 0.06 0.88 0.82]);
axis(ax, 'off');
z = r.selectedMode.lambda;
title(ax, {sprintf('表3：%s，SCR=%.5g，指定频段目标模态', ...
    c.model.structure, c.experiment.participationSCR), ...
    sprintf('lambda=%.6g %+.6gj，f=%.6g Hz', ...
    real(z), imag(z), r.selectedMode.frequencyHz)}, 'Interpreter', 'none');
T = r.table3;
rows = cell(height(T), 4);
for k = 1:height(T)
    rows{k, 1} = T.State{k};
    rows{k, 2} = sprintf('%.6g %+.6gj', T.RealP(k), T.ImagP(k));
    rows{k, 3} = sprintf('%.7g', T.AbsP(k));
    rows{k, 4} = sprintf('%.5f%%', 100*T.Normalized(k));
end
draw_table(ax, {'状态', '原始p', '|p|', '|p|归一化占比'}, ...
    rows, [0 0.18 0.59 0.80], c);
text(ax, 0, -0.04, sprintf('sum(p)=%.8g %+.3gj；sum(|p|)=%.8g。', ...
    real(r.selectedMode.sumRaw), imag(r.selectedMode.sumRaw), ...
    r.selectedMode.sumAbs), 'Interpreter', 'none', 'FontSize', 10);
end

function draw_table(ax, header, rows, x, c)
n = size(rows, 1);
dy = 0.88/(n+1);
for j = 1:numel(header)
    text(ax, x(j), 0.98, header{j}, 'FontWeight', 'bold', ...
        'FontName', c.plot.fontName, 'FontSize', 11, 'Interpreter', 'none');
end
for k = 1:n
    for j = 1:size(rows, 2)
        text(ax, x(j), 0.98-k*dy, rows{k, j}, ...
            'FontName', c.plot.fontName, 'FontSize', 10, 'Interpreter', 'none');
    end
end
line(ax, [0 1], [1.02 1.02], 'Color', [0.1 0.2 0.35]);
line(ax, [0 1], [0.98-dy/2 0.98-dy/2], 'Color', [0.1 0.2 0.35]);
line(ax, [0 1], [0.98-(n+0.6)*dy 0.98-(n+0.6)*dy], 'Color', [0.1 0.2 0.35]);
xlim(ax, [0 1]); ylim(ax, [0 1.06]);
end

function f = plot_participation_bars(r, c)
f = new_page(c, '所选模态状态参与度');
[values, order] = sort(r.table3.Normalized, 'descend');
n = min(c.plot.participationTopN, numel(order));
ax = axes(f, 'Position', [0.15 0.16 0.76 0.70]);
barh(ax, 100*values(n:-1:1));
yticks(ax, 1:n);
yticklabels(ax, r.table3.State(order(n:-1:1)));
set(ax, 'TickLabelInterpreter', 'none');
xlabel(ax, '|p|归一化参与度 / %'); grid(ax, 'on');
title(ax, sprintf('%s：SCR=%.4g，振荡频率=%.6g Hz', ...
    c.model.structure, c.experiment.participationSCR, ...
    r.selectedMode.frequencyHz), 'Interpreter', 'none');
end

function f = make_validation_page(r, c)
f = new_page(c, '方程、核验与适用范围');
ax = axes(f, 'Position', [0.07 0.05 0.86 0.86]);
axis(ax, 'off');
title(ax, sprintf('%s：建模核验与结论边界', c.model.structure), ...
    'FontSize', 19, 'Interpreter', 'none');
v = r.first.validation;
lines = {
    '状态方程：dx/dt=f(x,u)；在独立求得的同步工作点线性化。'
    'A=df/dx，B=df/du，y=[Pe,Qe,U_LL,omega,delta]，C=dy/dx，D=dy/du。'
    'J2=[0,-1;1,0]，eg=Eg*[cos(delta),-sin(delta)]。'
    'dig/dt=(vp-eg-Rg*ig)/Lg-omega*J2*ig。'
    'dio/dt=(eo-vp-Rf*io)/Lf-omega*J2*io。'
    'dvc/dt=(io-ig)/Cf-omega*J2*vc；vp=vc+Rc*(io-ig)。'
    'ev=[E,0]-Rv*iv-w0*Lv*J2*iv；dphi/dt=ev-vp。'
    'iref=kvp*(ev-vp)+kvi*phi+w0*Cf*J2*vp。'
    'eo=kip*(iref-io)+kii*gamma+w0*Lf*J2*io+kff*vp。'
    'dgamma/dt=iref-io；dPf/dt=wp*(Pe-Pf)；dQf/dt=wq*(Qe-Qf)。'
    'paper15：domega/dt=(Pref-Pf)/(J*w0)-Dp*domega/J。'
    'model17：J*domega/dt=(Pref-Pf-Komega*domega)/(w0+domega)-Ddirect*domega。'
    'ddelta/dt=omega-wg；两种有功环在同步稳态附近一阶等效。'
    'dE/dt=bq*(Qref-Qf+Dq*(Eref-Ufeedback))。'
    sprintf('工作点归一化残差：%.3e；状态/输入Jacobian误差：%.3e / %.3e。', ...
        r.first.equilibrium.dynamicResidual, v.jacobianRelativeError, ...
        v.inputJacobianRelativeError)
    sprintf('特征方程残差：%.3e；W*V-I：%.3e；参与因子和误差：%.3e。', ...
        v.eigenResidual, v.dualityError, v.participationSumError)
    sprintf('局部扰动幅值减半后，二阶余项比值：%.7g，理论应趋近4。', ...
        v.localRemainderRatio)
    sprintf('首次工况max(Re(lambda))=%.8g /s；不稳定根数=%d。', ...
        v.spectralAbscissa, v.unstableModeCount)
    '论文表2的+0.7651实部已表示失稳，不能与“恰好临界”混用。'
    '论文表3不应直接解读为百分比；本代码明确给出原始值、模及归一化值。'
    '特征根虚部/(2*pi)是dq包络频率；abc电气频谱还需检查50Hz两侧边带。'
    '表3默认分析0.05-50Hz频段；系统稳定性以全部特征值的最大实部为准。'
    sprintf('全局最危险模态：%.6g %+.6gj，频率%.6g Hz；参与因子见globalTable3。', ...
        real(r.globalMode.lambda), imag(r.globalMode.lambda), r.globalMode.frequencyHz)
    };
text(ax, 0, 0.98, lines, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 10, 'Interpreter', 'none');
end

function f = plot_structure_comparison(r, c)
f = new_page(c, '两种结构的完整SCR扫描比较');
ax = axes(f, 'Position', [0.12 0.16 0.78 0.70]);
plot(ax, r.SCR, r.spectralAbscissa, 'LineWidth', 1.7, ...
    'DisplayName', r.cfg.model.structure);
hold(ax, 'on');
plot(ax, r.comparison.SCR, r.comparison.spectralAbscissa, ...
    'LineWidth', 1.7, 'DisplayName', r.comparison.cfg.model.structure);
yline(ax, 0, '--k', 'HandleVisibility', 'off');
grid(ax, 'on');
xlabel(ax, 'SCR'); ylabel(ax, 'max Re(lambda) / s^{-1}');
legend(ax, 'Location', 'best', 'Interpreter', 'none');
title(ax, '相同参数下的结构对照：稳定性以各自全部特征根为准');
end

%% 有功频率环实验：同一模态子空间跟踪，不改变基准参数
function r = run_active_frequency_experiment(c)
assert(strcmp(c.model.structure, 'model17'), '有功频率环实验请使用model17。');
e = c.activeExperiment;
a = analyse_case(c, e.SCR, []);
[~, mode] = make_participation_table(a, c);
assert(abs(imag(mode.lambda)) > c.numeric.conjugateTolerance, ...
    '实验基准必须选择一个振荡模态。');
s = state_scales(a.parameters);
ctx = struct();
ctx.base = a;
ctx.scale = s;
ctx.g0 = [c.control.J, c.control.Komega, c.control.Ddirect];
ctx.v0 = a.rightEigenvectors./s;
ctx.w0 = a.leftEigenvectors.*s.';
[~, conjugateIndex] = min(abs(a.lambda-conj(mode.lambda)));
ctx.pair0 = [mode.index, conjugateIndex];
ctx.projector0 = ctx.v0(:, ctx.pair0)*ctx.w0(ctx.pair0, :);
ctx.lambda0 = a.lambda(ctx.pair0);
baseline = track_active_case(ctx, ctx.g0, c);
baseline.Group = 'baseline';
rows = baseline;

% 三个物理参数的联合网格；所有其余参数和工作点保持相同。
for jf = e.JFactors
    for kf = e.KomegaFactors
        for df = e.DdirectFactors
            gain = ctx.g0.*[jf kf df];
            item = track_active_case(ctx, gain, c);
            item.Group = 'joint';
            rows(end+1) = item; %#ok<AGROW>
        end
    end
end

% 固定原有频差功率增益，只联合调节J与等效阻尼，便于直接试验。
for jf = e.JFactors
    for factor = e.DeqFactors
        directD = baseline.Deq*factor-c.control.Komega/c.base.w0;
        if directD < 0
            continue;
        end
        item = track_active_case(ctx, ...
            [c.control.J*jf, c.control.Komega, directD], c);
        item.Group = 'fixed_Komega';
        rows(end+1) = item; %#ok<AGROW>
    end
end

% 单参数扫描单独保存，避免把联合扫描的结果误当作单参数作用。
oneD = cell(3, 1);
factorSets = {e.JFactors, e.KomegaFactors, e.DdirectFactors};
for parameterIndex = 1:3
    values = factorSets{parameterIndex};
    axisRows = repmat(baseline, numel(values), 1);
    for k = 1:numel(values)
        gain = ctx.g0;
        gain(parameterIndex) = gain(parameterIndex)*values(k);
        axisRows(k) = track_active_case(ctx, gain, c);
    end
    oneD{parameterIndex} = axisRows;
end

% 不能只按参与比例选最小值：额外约束原始模值、阻尼、衰减速度与模态跟踪。
eligible = [rows.Oscillatory] & [rows.Alpha] < -c.numeric.stabilityTolerance ...
    & [rows.DampingRatio] >= baseline.DampingRatio-1e-6 ...
    & [rows.TargetReal] <= baseline.TargetReal+1e-6 ...
    & [rows.DeltaPercent] < baseline.DeltaPercent-1e-6 ...
    & [rows.DeltaAbs] < baseline.DeltaAbs-1e-6 ...
    & [rows.TrackingScore] < e.maxTrackingStep;
indices = find(eligible);
[~, order] = sort([rows(indices).DeltaPercent]);
indices = indices(order);

% 在额外SCR点检查所有特征根。该检查不宣称网格之间也必然稳定。
scr = unique(e.validationSCR(:).');
gridBase = cell(size(scr));
baseAlpha = zeros(size(scr));
for k = 1:numel(scr)
    gridBase{k} = analyse_case(c, scr(k), []);
    baseAlpha(k) = max(real(gridBase{k}.lambda));
end
chosen = [];
chosenAlpha = [];
validated = false(size(rows));
for index = indices
    gain = [rows(index).J, rows(index).Komega, rows(index).Ddirect];
    alpha = zeros(size(scr));
    for k = 1:numel(scr)
        A = active_matrix(gridBase{k}, gain, c.base.w0);
        alpha(k) = max(real(eig(A)));
    end
    if all(alpha < -c.numeric.stabilityTolerance)
        validated(index) = true;
        if isempty(chosen)
            chosen = rows(index);
            chosenAlpha = alpha;
        end
    end
end

% 验证把Komega与Ddirect重新分配、但保持Deq不变的等价性。
Kvalues = c.control.Komega*[0 0.5 1 1.5];
fixedDeq = zeros(numel(Kvalues), 3);
baseAs = (a.A.*s.')./s;
for k = 1:numel(Kvalues)
    D = baseline.Deq-Kvalues(k)/c.base.w0;
    other = active_matrix(a, [c.control.J Kvalues(k) D], c.base.w0);
    fixedDeq(k, :) = [Kvalues(k), D, norm(other-baseAs, 'fro')];
end

r = struct('cfg', c, 'first', a, 'baseline', baseline, 'cases', rows, ...
    'oneD', {oneD}, 'eligible', eligible, 'validated', validated, ...
    'candidate', chosen, 'validationSCR', scr, ...
    'validationBaselineAlpha', baseAlpha, 'validationCandidateAlpha', chosenAlpha, ...
    'fixedDeq', fixedDeq);
r.caseTable = active_case_table(rows);
[~, minIndex] = min([rows.DeltaPercent]);
r.minimumOnly = rows(minIndex);
if ~isempty(chosen)
    % 用完整非线性方程重新求导，再核验快速更新行的方法。
    candidateCfg = c;
    candidateCfg.control.J = chosen.J;
    candidateCfg.control.Komega = chosen.Komega;
    candidateCfg.control.Ddirect = chosen.Ddirect;
    candidateCfg = refresh_derived_parameters(candidateCfg);
    fullCase = analyse_case(candidateCfg, e.SCR, a.equilibrium.z);
    fullScaled = (fullCase.A.*s.')./s;
    fastScaled = active_matrix(a, [chosen.J chosen.Komega chosen.Ddirect], c.base.w0);
    error = norm(fullScaled-fastScaled, 'fro')/max(1, norm(fullScaled, 'fro'));
    assert(error < 1e-8, '快速参数更新与完整非线性Jacobian不一致。');
    r.candidateFullCase = fullCase;
    r.candidateJacobianError = error;
end
fprintf('\n=== 有功频率环delta参与度实验：SCR=%.6g ===\n', e.SCR);
fprintf('基准delta占比=%.6g%%，原始模值=%.6g，频率=%.6g Hz。\n', ...
    baseline.DeltaPercent, baseline.DeltaAbs, baseline.FrequencyHz);
fprintf('联合候选数=%d，满足局部约束=%d，通过额外SCR检查=%d。\n', ...
    numel(rows)-1, nnz(eligible), nnz(validated));
if isempty(chosen)
    fprintf('本次范围未找到同时满足参与度下降与稳定性约束的候选。\n');
else
    fprintf('候选：J=%.9g，Komega=%.9g，Ddirect=%.9g。\n', ...
        chosen.J, chosen.Komega, chosen.Ddirect);
    fprintf('delta占比=%.6g%%，原始模值=%.6g，目标根实部=%.6g。\n', ...
        chosen.DeltaPercent, chosen.DeltaAbs, chosen.TargetReal);
end
fprintf('基准配置未被覆盖，完整数值见result.caseTable。\n');
end

function As = active_matrix(a, gain, w0)
% 同步稳态Pf=Pref，三个有功参数不改变工作点。
% 因而只有频率状态这一行变化；以下是完整非线性方程的精确一阶导数。
J = gain(1);
Deq = gain(3)+gain(2)/w0;
assert(J > 0 && Deq > 0, 'J和等效阻尼必须为正。');
A = a.A;
A(9, :) = 0;
A(9, 7) = -1/(J*w0);
A(9, 9) = -Deq/J;
s = state_scales(a.parameters);
As = (A.*s.')./s;
end

function item = track_active_case(ctx, gain, c)
previousP = ctx.projector0;
previousZ = ctx.lambda0;
maxScore = 0;
steps = c.activeExperiment.continuationSteps;
for k = 1:steps
    currentGain = ctx.g0+(gain-ctx.g0)*(k/steps);
    As = active_matrix(ctx.base, currentGain, c.base.w0);
    [v, L] = eig(As);
    z = diag(L);
    w = v\eye(size(v));
    groups = active_mode_pairs(z, c.numeric.conjugateTolerance);
    scores = zeros(size(groups, 1), 1);
    for j = 1:size(groups, 1)
        ix = groups(j, :);
        P = v(:, ix)*w(ix, :);
        zz = z(ix);
        dP = norm(P-previousP, 'fro')/max(1, norm(previousP, 'fro'));
        dZ = min(max(abs(zz-previousZ)), max(abs(flipud(zz)-previousZ))) ...
            /max(1, max(abs(previousZ)));
        scores(j) = dP+0.2*dZ;
    end
    [score, selected] = min(scores);
    pair = groups(selected, :);
    previousP = v(:, pair)*w(pair, :);
    previousZ = z(pair);
    maxScore = max(maxScore, score);
end

% 对共轭对，此模值平均与取正虚部单根完全一致。
% 跨入两个实根时仍计算原二维子空间，防止换模态后虚假降低参与度。
raw = v.*w.';
pairAbs = mean(abs(raw(:, pair)), 2);
normalized = pairAbs/sum(pairAbs);
zz = z(pair);
item = struct('J', gain(1), 'Komega', gain(2), 'Ddirect', gain(3), ...
    'Deq', gain(3)+gain(2)/c.base.w0, 'DeltaPercent', 100*normalized(10), ...
    'DeltaAbs', pairAbs(10), 'Alpha', max(real(z)), ...
    'TargetReal', max(real(zz)), 'FrequencyHz', max(abs(imag(zz)))/(2*pi), ...
    'DampingRatio', min(-real(zz)./max(abs(zz), realmin)), ...
    'TrackingScore', maxScore, ...
    'Oscillatory', any(abs(imag(zz)) > c.numeric.conjugateTolerance), ...
    'Pair', zz, 'Normalized', normalized, 'Group', '');
end

function pairs = active_mode_pairs(z, tol)
pairs = zeros(0, 2);
for k = find(imag(z) > tol).'
    [~, j] = min(abs(z-conj(z(k))));
    pairs(end+1, :) = [k j]; %#ok<AGROW>
end
realIndices = find(abs(imag(z)) <= tol);
if numel(realIndices) >= 2
    pairs = [pairs; nchoosek(realIndices, 2)];
end
assert(~isempty(pairs), '没有可跟踪的模态对。');
end

function T = active_case_table(rows)
T = table([rows.J].', [rows.Komega].', [rows.Ddirect].', [rows.Deq].', ...
    [rows.DeltaPercent].', [rows.DeltaAbs].', [rows.Alpha].', ...
    [rows.TargetReal].', [rows.FrequencyHz].', [rows.DampingRatio].', ...
    [rows.Oscillatory].', [rows.TrackingScore].', ...
    'VariableNames', {'J','Komega','Ddirect','Deq','DeltaPercent','DeltaAbs', ...
    'MaxRealPart','TargetRealPart','FrequencyHz','DampingRatio', ...
    'Oscillatory','TrackingScore'});
end

function figures = make_active_experiment_figures(r, c)
figures = gobjects(4, 1);
figures(1) = new_page(c, '有功频率环实验摘要');
ax = axes(figures(1), 'Position', [0.07 0.07 0.86 0.84]); axis(ax, 'off');
title(ax, '有功频率环与delta参与度：保留基准参数的对照实验');
b = r.baseline;
lines = {
    sprintf('工况：SCR=%.6g，Pref=%.6g W，Pe/Qe低通=%.6g/%.6g Hz。', ...
        c.activeExperiment.SCR, c.control.Pref, c.measurement.wp/(2*pi), c.measurement.wq/(2*pi))
    sprintf('基准：J=%.8g，Komega=%.8g，Ddirect=%.8g，Deq=%.8g。', ...
        b.J, b.Komega, b.Ddirect, b.Deq)
    sprintf('基准目标根：%.8g %+.8gj，频率%.8g Hz，阻尼比%.8g。', ...
        real(b.Pair(1)), imag(b.Pair(1)), b.FrequencyHz, b.DampingRatio)
    sprintf('基准delta占比=%.8g%%；原始参与因子模=%.8g。', b.DeltaPercent, b.DeltaAbs)
    '沿参数连续变化跟踪同一二维模态子空间；不使用每次重新自动选择目标模态的方法。'
    '局部筛选要求：原始模与占比同时下降、共轭对保留、阻尼比不降、实部不右移、全部根稳定。'
    '还需通过参数文件validationSCR中各工况的全部特征根稳定性检查。'
    '参与百分比不是能量占比；降低delta参与度本身不等于提高系统稳定性。'
    'J和等效阻尼不改变当前同步工作点，所以只需更新频率状态行后再计算全部特征根。'
    '同J且Deq不变时，重新分配Komega与Ddirect不改变当前小信号状态矩阵。'
    '结果属于17阶连续平均模型；尚需原SLX的小扰动验证，未自动修改SLX。'
    };
if ~isempty(r.candidate)
    a = r.candidate;
    lines{end+1} = sprintf('候选：J=%.8g，Komega=%.8g，Ddirect=%.8g。', a.J, a.Komega, a.Ddirect);
    lines{end+1} = sprintf('候选delta=%.6g%%，原始模=%.6g，频率=%.6gHz，阻尼比=%.6g。', ...
        a.DeltaPercent, a.DeltaAbs, a.FrequencyHz, a.DampingRatio);
else
    lines{end+1} = '在本次范围和约束下没有找到合适候选，基准设置保持不变。';
end
text(ax, 0, 0.98, lines, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'Interpreter', 'none', 'FontSize', 11);

figures(2) = new_page(c, '有功环三个参数的单独扫描');
layout = tiledlayout(figures(2), 3, 3, 'TileSpacing', 'compact');
parameters = {'J', 'Komega', 'Ddirect'};
metrics = {'DeltaPercent', 'TargetReal', 'DampingRatio'};
labels = {'delta归一化占比 / %', '目标根最大实部 / s^{-1}', '目标阻尼比'};
for row = 1:3
    for column = 1:3
        ax = nexttile(layout);
        items = r.oneD{column};
        plot(ax, [items.(parameters{column})], [items.(metrics{row})], '-o', 'LineWidth', 1.2);
        yline(ax, b.(metrics{row}), '--k'); grid(ax, 'on');
        xlabel(ax, parameters{column}, 'Interpreter', 'none'); ylabel(ax, labels{row});
    end
end
title(layout, '单参数实验：跨入两个实根的结果仅作诊断，不用于共轭模态推荐');

figures(3) = new_page(c, '联合扫描与候选比较');
ax = axes(figures(3), 'Position', [0.08 0.09 0.86 0.80]); axis(ax, 'off');
items = [b, r.minimumOnly];
names = {'基准', '仅占比最小'};
if ~isempty(r.candidate)
    items = [b, r.candidate, r.minimumOnly];
    names = {'基准', '通过约束候选', '仅占比最小'};
end
rows = cell(numel(items), 7);
for k = 1:numel(items)
    a = items(k);
    rows(k, :) = {names{k}, sprintf('%.5g', a.J), sprintf('%.6g', a.Komega), ...
        sprintf('%.6g', a.Ddirect), sprintf('%.5f', a.DeltaPercent), ...
        sprintf('%.6g', a.TargetReal), sprintf('%.5g', a.DampingRatio)};
end
draw_table(ax, {'方案','J','Komega','Ddirect','delta/%','目标实部','阻尼比'}, ...
    rows, [0 0.20 0.29 0.43 0.57 0.73 0.88], c);
title(ax, sprintf('扫描%d组；局部通过%d组；额外SCR通过%d组', ...
    numel(r.cases)-1, nnz(r.eligible), nnz(r.validated)));

figures(4) = new_page(c, '候选的跨SCR稳定性检查');
ax = axes(figures(4), 'Position', [0.12 0.16 0.78 0.70]);
plot(ax, r.validationSCR, r.validationBaselineAlpha, '-o', 'LineWidth', 1.5, 'DisplayName', '基准');
hold(ax, 'on');
if ~isempty(r.candidate)
    plot(ax, r.validationSCR, r.validationCandidateAlpha, '-s', 'LineWidth', 1.5, 'DisplayName', '候选');
end
yline(ax, 0, '--k', 'HandleVisibility', 'off'); grid(ax, 'on');
xlabel(ax, 'SCR'); ylabel(ax, '全部特征根最大实部 / s^{-1}');
legend(ax, 'Location', 'best');
title(ax, '仅对列出的SCR采样点进行检查，不表示网格之间已获稳定性保证');
end
