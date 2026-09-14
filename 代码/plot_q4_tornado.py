"""问题 4 灵敏度图：① 龙卷风图 ② 空间网格收敛曲线。

龙卷风图数值来源：`sensitivity_q4.py` 实算结果（约 10 min 机时），
基准 t* = 51.23 h（N=160、实测插值 R(t)、Δt=10 s、处处判据）。
网格收敛曲线数值来源：论文表 17（q4.py 实算），基准取 N=320 的 51.12 h。
"""
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['SimSun', 'Times New Roman'],
    'font.size': 11, 'axes.unicode_minus': False,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 4, 'ytick.major.size': 4,
})
BASE = 51.23                                   # 龙卷风图基准 / h
BLUE, ORANGE, RED = '#0072BD', '#D95319', '#C00000'

# 龙卷风图：(因素, 最小 t_end, 最大 t_end, 是否高亮, 条端标注)
ROWS = [
    ('判据口径\n(表面/平均/处处)', 12.78, 51.23, False, '12.78 ~ 51.23 h'),
    ('空间网格 $N$\n(40 ~ 640)', 51.10, 54.48, True, '51.10 ~ 54.48 h'),
    ('收缩律\n(实测/单指数/Weibull)', 51.13, 51.23, False, '51.13 ~ 51.23 h'),
    ('时间步 $\\Delta t$\n(10 / 2.5 s)', 51.23, 51.23, False, '51.23 h（无变化）'),
]

# 表 17：问题四空间网格敏感性检验（基准 N=320）
GRID_N = [40, 80, 160, 320, 640]
GRID_T = [54.48, 51.73, 51.23, 51.12, 51.10]
T_CONV = 51.12


def draw_tornado():
    fig, ax = plt.subplots(figsize=(7.6, 3.8))
    for i, (name, lo, hi, hl, tag) in enumerate(ROWS):
        d0, d1 = 100 * (lo / BASE - 1), 100 * (hi / BASE - 1)
        w = max(d1 - d0, 1.0)                      # 最小可见宽度，精确值由标签给出
        c = ORANGE if hl else BLUE
        ax.barh(i, w, left=d0, height=0.52, color=c, alpha=0.9,
                edgecolor='black', lw=0.6)
        ax.text(max(d1, d0 + w) + 1.6, i, tag, va='center', ha='left', fontsize=9.5)
        ax.text(d0 - 1.6, i, '%.2f %%' % d0, va='center', ha='right',
                fontsize=9.5, color='0.25')
    ax.axvline(0, color='0.3', lw=1.3, ls='--')
    ax.set_yticks(range(len(ROWS)))
    ax.set_yticklabels([r[0] for r in ROWS], fontsize=10)
    ax.set_ylim(-0.6, len(ROWS) - 0.4)
    ax.set_xlim(-95, 30)
    ax.invert_yaxis()
    ax.set_xlabel('$t_{end}$ 相对基准的变化/%')
    ax.set_title('问题四 $t_{end}$ 的灵敏度龙卷风图（基准 $t^*$=%.2f h）' % BASE,
                 fontsize=12, pad=10)
    ax.annotate('基准取值：$N$=160、实测插值 $R(t)$、$\\Delta t$=10 s、处处判据',
                xy=(0.99, -0.30), xycoords='axes fraction', ha='right',
                va='top', fontsize=9, color='0.35')
    fig.tight_layout()
    fig.savefig('问题4_灵敏度龙卷风图.png', dpi=300)
    print('已保存 问题4_灵敏度龙卷风图.png')


def draw_grid():
    """表 17 数据：(a) 烘干时间 vs N；(b) 与基准差值 vs N（双对数，看收敛阶）。"""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    a1.axhline(T_CONV, color=RED, lw=1.3, ls='--',
               label='收敛值 %.2f h（$N$=320）' % T_CONV)
    a1.plot(GRID_N, GRID_T, 'o-', color=BLUE, lw=1.8, ms=6,
            mfc='white', mew=1.4, label='烘干时间')
    for x, y in zip(GRID_N, GRID_T):
        a1.annotate('%+.2f h' % (y - T_CONV), (x, y), textcoords='offset points',
                    xytext=(0, 10), ha='center', fontsize=9, color='0.25')
    a1.set_xscale('log', base=2)
    a1.set_xticks(GRID_N)
    a1.set_xticklabels([str(n) for n in GRID_N])
    a1.set_xlim(30, 800)
    a1.set_ylim(50.9, 55.2)
    a1.set_xlabel('空间区间数 $N$')
    a1.set_ylabel('烘干时间/h')
    a1.set_title('(a) 烘干时间随网格的收敛', fontsize=11)
    a1.legend(fontsize=9, frameon=True, edgecolor='black', fancybox=False)

    # (b) 相对 N=320 的偏差，双对数看收敛阶（取 N<=160 的三点）
    n = [GRID_N[i] for i in range(3)]
    d = [abs(GRID_T[i] - T_CONV) for i in range(3)]
    a2.loglog(n, d, 's-', color=ORANGE, lw=1.8, ms=6, mfc='white', mew=1.4,
              label='|$\\Delta t_{end}$|（表 17）')
    p = math.log(d[0] / d[-1]) / math.log(n[-1] / n[0])       # 实测收敛阶
    xr = [n[0], n[-1]]
    a2.loglog(xr, [d[0], d[0] * (xr[1] / xr[0]) ** (-2)], 'k--', lw=1.2,
              label='参考斜率 -2（二阶）')
    a2.set_xticks(n)
    a2.set_xticklabels([str(v) for v in n])
    a2.set_yticks([0.1, 0.2, 0.5, 1, 2, 5])
    a2.set_yticklabels(['0.1', '0.2', '0.5', '1', '2', '5'])   # 避开 SimSun 缺 U+2212
    a2.minorticks_off()
    a2.set_xlabel('空间区间数 $N$')
    a2.set_ylabel('与基准差值 $|\\Delta t_{end}|$/h')
    a2.set_title('(b) 偏差的收敛阶（实测 %.2f 阶）' % p, fontsize=11)
    a2.legend(fontsize=9, frameon=True, edgecolor='black', fancybox=False)
    fig.tight_layout()
    fig.savefig('问题4_网格收敛曲线.png', dpi=300)
    print('已保存 问题4_网格收敛曲线.png  (实测收敛阶 %.2f)' % p)
    return p


if __name__ == '__main__':
    draw_tornado()
    p = draw_grid()
    assert abs(100 * (54.48 / BASE - 1) - 6.34) < 0.05
    assert 2.0 < p < 3.0
