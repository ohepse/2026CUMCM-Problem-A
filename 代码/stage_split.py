"""问题 2：预热平衡 / 恒温干燥分界点判定（纯温度判据）+ 可视化。

判据（只用温度，不用湿度）:
  烘房侧（附件 1 直接判定）:
    A1 分段变点：前段二次升温 + 后段常数平台，最小残差平方和
    A2 一阶指数饱和 T = Tf - (Tf-T0)e^{-t/tau}，取达 95 % 终值时刻
    A3 首次进入平台 +-2sigma 且此后不再离开
  药材侧（问题 1 模型延长积到 4 h）:
    B1 药材中心与烘房温差 < 1 K
    B2 药材表面-中心径向温差 < 1 K / < 0.5 K
    B3 药材表面与烘房温差 < 0.5 K
输出：判据表 + 两张图
  问题2_烘房温度变化曲线.png
  问题2_药材内部温度变化曲线.png
"""
import math
import openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import solve_q1 as S

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['SimSun', 'Times New Roman'],
    'font.size': 11, 'axes.unicode_minus': False,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 4, 'ytick.major.size': 4,
})
BLUE, ORANGE, GREEN = '#0072BD', '#D95319', '#77AC30'
RED = '#C00000'
TSTAR = 2.0                                   # 最终选定的分界点 / h


def load_air():
    rows = list(openpyxl.load_workbook(S.AIR, data_only=True)
                .active.iter_rows(min_row=2, values_only=True))
    return [r[0] for r in rows], [r[1] for r in rows]


def solve3(M, v):
    for c in range(3):
        p = max(range(c, 3), key=lambda r: abs(M[r][c]))
        M[c], M[p], v[c], v[p] = M[p], M[c], v[p], v[c]
        for r in range(c + 1, 3):
            f = M[r][c] / M[c][c]
            for cc in range(c, 3):
                M[r][cc] -= f * M[c][cc]
            v[r] -= f * v[c]
    x = [0.0] * 3
    for r in (2, 1, 0):
        x[r] = (v[r] - sum(M[r][cc] * x[cc] for cc in range(r + 1, 3))) / M[r][r]
    return x


def change_point(t, y):
    """前段二次多项式 + 后段常数，返回 (t*, 后段均值, R2)。"""
    n = len(t)
    best = None
    for i in range(20, n - 40):
        A = [[1.0, tt, tt * tt] for tt in t[:i]]
        M = [[sum(A[k][a] * A[k][b] for k in range(i)) for b in range(3)]
             for a in range(3)]
        v = [sum(A[k][a] * y[k] for k in range(i)) for a in range(3)]
        b = solve3(M, v)
        sse = sum((y[k] - sum(b[j] * A[k][j] for j in range(3))) ** 2
                  for k in range(i))
        seg = y[i:]
        m = sum(seg) / len(seg)
        sse += sum((z - m) ** 2 for z in seg)
        if best is None or sse < best[0]:
            best = (sse, t[i], m)
    tss = sum((z - sum(y) / n) ** 2 for z in y)
    return best[1], best[2], 1 - best[0] / tss


def exp_settle(t, y):
    """一阶指数饱和拟合，返回 (Tf, tau, 达 95 % 终值时刻)。"""
    best = None
    for i in range(960, 1041):                    # Tf = 48 .. 52 degC
        Tf = i * 0.05
        for j in range(20, 151):                  # tau = 400 .. 3000 s
            tau = j * 20.0
            e = sum((y[k] - (Tf - (Tf - y[0]) * math.exp(-t[k] / tau))) ** 2
                    for k in range(0, len(t), 2))
            if best is None or e < best[0]:
                best = (e, Tf, tau)
    return best[1], best[2], -best[2] * math.log(0.05)


def entry(t, y, i0, k=2.0):
    """平台段均值/标准差；返回首次进入 mean+-k*sigma 且不再离开的时刻。"""
    seg = y[i0:]
    m = sum(seg) / len(seg)
    sd = math.sqrt(sum((z - m) ** 2 for z in seg) / (len(seg) - 1))
    last = max((j for j in range(i0, len(t)) if abs(y[j] - m) > k * sd),
               default=i0 - 1)
    return t[min(last + 1, len(t) - 1)], m, sd


def material():
    """用问题 1 模型积到 4 h，返回 (t, T_center, T_surf, T_air)。"""
    S.AIRDATA = [(r[0], r[1], r[2]) for r in openpyxl.load_workbook(
        S.AIR, data_only=True).active.iter_rows(min_row=2, values_only=True)]
    S.setup(80)
    snaps, _, _ = S.solve(tmax=14400.0)
    n = S.N
    return ([s[0] for s in snaps], [s[1][0] for s in snaps],
            [s[1][n] for s in snaps], [S.air_at(s[0], 1) for s in snaps])


def first_stable(t, cond):
    """cond 此后恒成立的首个时刻。"""
    bad = [j for j in range(len(t)) if not cond(j)]
    return t[min(bad[-1] + 1, len(t) - 1)] if bad else t[0]


def draw_air(t, y, ts, m, sd):
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.axhspan(m - 2 * sd, m + 2 * sd, color='0.6', alpha=0.25, lw=0,
               label='恒温平台 $\\pm2\\sigma$（%.2f$^\\circ$C）' % (2 * sd))
    ax.plot([x / 3600 for x in t], y, color=BLUE, lw=1.6, label='附件 1 烘房温度')
    ax.axhline(m, color='0.35', lw=1.0, ls=':')
    ax.axvline(ts / 3600, color=ORANGE, lw=1.4, ls='--',
               label='温度变点 %.2f h' % (ts / 3600))
    ax.axvline(TSTAR, color=RED, lw=1.8, label='分界点 $t^*$ = %.1f h' % TSTAR)
    ax.set_xlim(0, 4)
    ax.set_ylim(27, 52)
    ax.set_xlabel('时间/h')
    ax.set_ylabel('烘房温度/$^\\circ$C')
    ax.set_title('烘房温度变化曲线', fontsize=12, pad=10)
    ax.legend(loc='lower right', fontsize=9, frameon=True, edgecolor='black',
              fancybox=False)
    fig.tight_layout()
    fig.savefig('问题2_烘房温度变化曲线.png', dpi=300)
    print('已保存 问题2_烘房温度变化曲线.png')


def draw_material(tm, tc, tsurf, tair, tb1):
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    th = [x / 3600 for x in tm]
    ax.fill_between(th, [v - 1 for v in tair], [v + 1 for v in tair],
                    color='0.6', alpha=0.25, lw=0, label='烘房温度 $\\pm$1 K')
    ax.plot(th, tair, color='0.35', lw=1.2, ls=':', label='烘房温度（参考）')
    ax.plot(th, tsurf, color=ORANGE, lw=1.6, label='药材表面 $r=2$ cm')
    ax.plot(th, tc, color=BLUE, lw=1.6, label='药材中心 $r=0$')
    ax.axvline(tb1 / 3600, color=GREEN, lw=1.4, ls='--',
               label='中心进入 $\\pm$1 K：%.2f h' % (tb1 / 3600))
    ax.axvline(TSTAR, color=RED, lw=1.8, label='分界点 $t^*$ = %.1f h' % TSTAR)
    ax.set_xlim(0, 4)
    ax.set_xlabel('时间/h')
    ax.set_ylabel('药材温度/$^\\circ$C')
    ax.set_title('药材内部温度变化曲线', fontsize=12, pad=10)
    ax.legend(loc='lower right', fontsize=9, frameon=True, edgecolor='black',
              fancybox=False)
    fig.tight_layout()
    fig.savefig('问题2_药材内部温度变化曲线.png', dpi=300)
    print('已保存 问题2_药材内部温度变化曲线.png')


if __name__ == '__main__':
    t, y = load_air()
    ts, m, r2 = change_point(t, y)
    Tf, tau, t95 = exp_settle(t, y)
    te, mm, sd = entry(t, y, t.index(ts))
    tm, tc, tsurf, tair = material()
    tb1 = first_stable(tm, lambda j: abs(tair[j] - tc[j]) < 1.0)
    tb2 = first_stable(tm, lambda j: abs(tsurf[j] - tc[j]) < 1.0)
    tb3 = first_stable(tm, lambda j: abs(tsurf[j] - tc[j]) < 0.5)
    tb4 = first_stable(tm, lambda j: abs(tair[j] - tsurf[j]) < 0.5)

    print('============ 烘房侧（附件 1 温度，纯数据） ============')
    print('A1 分段变点              t* = %5.0f s = %.2f h  (R2=%.4f, 平台 %.2f C)'
          % (ts, ts / 3600, r2, m))
    print('A2 指数饱和 Tf=%.2f C, tau=%4.0f s -> 达 95%% 终值 %.0f s = %.2f h'
          % (Tf, tau, t95, t95 / 3600))
    print('A3 首次稳定在平台 +-2sigma(%.2f C)      %.0f s = %.2f h'
          % (2 * sd, te, te / 3600))
    print('============ 药材侧（问题 1 模型积到 4 h） ============')
    print('B1 中心与烘房温差 < 1 K         %.0f s = %.2f h' % (tb1, tb1 / 3600))
    print('B2 表面-中心温差 < 1 K          %.0f s = %.2f h' % (tb2, tb2 / 3600))
    print('B3 表面-中心温差 < 0.5 K        %.0f s = %.2f h' % (tb3, tb3 / 3600))
    print('B4 表面与烘房温差 < 0.5 K       %.0f s = %.2f h' % (tb4, tb4 / 3600))
    print('============ 结论 ============')
    print('烘房 %.2f h 已恒温；药材 %.2f h 才处处与烘房差 < 1 K，故取 t* = %.1f h'
          % (ts / 3600, tb1 / 3600, TSTAR))

    assert 1.8 < tb1 / 3600 < 2.2 and 1.3 < ts / 3600 < 1.8, (ts / 3600, tb1 / 3600)
    draw_air(t, y, ts, m, sd)
    draw_material(tm, tc, tsurf, tair, tb1)
