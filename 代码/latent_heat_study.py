"""水汽化潜热项能否忽略的定量论证（2026 A 题 问题 1 参数）。

模型（问题 1，一维轴对称，有限体积 + 显式时间推进）:
    dX/dt = r^-1 d/dr ( r D(X) dX/dr )
    rho*cp dT/dt = r^-1 d/dr ( k r dT/dr ) + rho_d*Lv dX/dt
    r=R:  k dT/dr = h (T_inf - T),   D dX/dr = hm (X_inf - X)
对照实验：同一套参数下令 Lv=0（即删掉潜热项），比较温度场。

自检：能量守恒残差 |对流供热 - (显热累积 + 潜热消耗)| / 对流供热 < 2%。
"""

import math

# ---------- 附录 2 参数（问题 1） ----------
RHO, CP, K = 820.0, 2600.0, 0.36      # kg/m^3, J/(kg K), W/(m K)
H, HM = 25.0, 8.0e-7                  # W/(m^2 K), m/s
X0, T0 = 2.55, 28.0                   # kg/kg, degC
R = 0.02                              # m
RHO_D = RHO / (1 + X0)                # 干基骨架密度 231.0 kg/m^3
LV = 2.43e6                           # J/kg, 30 degC 附近水的汽化潜热


def D_of(X):
    return 7e-9 * math.exp(-0.89 / X) if X > 1e-9 else 0.0


def solve(T_inf, X_inf, Lv, N=20, dt=0.5, tmax=1800.0, hm=HM):
    """返回 (T, X, 能量账)。T/X 为节点值；能量账单位 J/m（单位长度）。"""
    dr = R / N
    r = [i * dr for i in range(N + 1)]
    V, A = [], []                                     # 控制体体积、外侧面面积
    for i in range(N + 1):
        a = max(r[i] - dr / 2, 0.0)
        b = R if i == N else r[i] + dr / 2
        V.append(math.pi * (b * b - a * a))
        A.append(2 * math.pi * min(r[i] + dr / 2, R))

    T, X = [T0] * (N + 1), [X0] * (N + 1)
    heat_in = sensible = latent = 0.0
    for step in range(int(tmax / dt)):
        D = [D_of(x) for x in X]
        dX = [0.0] * (N + 1)
        dT = [0.0] * (N + 1)

        def Dh(a, b):                                  # 界面扩散系数：调和平均
            return 2 * a * b / (a + b) if a + b > 0 else 0.0

        # 水分
        dX[0] = A[0] * Dh(D[0], D[1]) * (X[1] - X[0]) / dr / V[0]
        for i in range(1, N):
            g_in = A[i - 1] * Dh(D[i - 1], D[i]) * (X[i] - X[i - 1]) / dr
            g_out = A[i] * Dh(D[i], D[i + 1]) * (X[i + 1] - X[i]) / dr
            dX[i] = (g_out - g_in) / V[i]
        dX[N] = (A[N] * hm * (X_inf - X[N])
                 - A[N - 1] * Dh(D[N - 1], D[N]) * (X[N] - X[N - 1]) / dr) / V[N]

        # 温度
        dT[0] = A[0] * K * (T[1] - T[0]) / dr / V[0] / (RHO * CP)
        for i in range(1, N):
            g_in = A[i - 1] * K * (T[i] - T[i - 1]) / dr
            g_out = A[i] * K * (T[i + 1] - T[i]) / dr
            dT[i] = (g_out - g_in) / V[i] / (RHO * CP)
        dT[N] = (A[N] * H * (T_inf - T[N])
                 - A[N - 1] * K * (T[N] - T[N - 1]) / dr) / V[N] / (RHO * CP)

        for i in range(N + 1):
            X[i] += dt * dX[i]
            T[i] += dt * (dT[i] + RHO_D * Lv * dX[i] / (RHO * CP))

        heat_in += dt * A[N] * H * (T_inf - T[N])

    sensible = sum(RHO * CP * (T[i] - T0) * V[i] for i in range(N + 1))
    latent = sum(RHO_D * Lv * (X0 - X[i]) * V[i] for i in range(N + 1))
    return T, X, (heat_in, sensible, latent)


def criterion(T_inf=60.0, X_inf=0.01, hm=HM):
    """判据 Pi = 最大潜热需求 / 最大对流供热（Pi<<1 才可忽略潜热）。"""
    return RHO_D * hm * LV * (X0 - X_inf) / (H * (T_inf - T0))


def sweep(T_inf=60.0):
    print(f"\n=== 灵敏度扫描 (T_inf={T_inf} degC) ===")
    print(f"{'hm/HM':>7} {'X_inf':>6} {'Pi':>7} {'潜热占供热':>10} {'最大温差K':>9}")
    for f in (0.05, 0.1, 0.5, 1.0, 2.0, 5.0):
        for x_inf in (0.01, 0.05):
            Ta, _, ea = solve(T_inf, x_inf, LV, hm=f * HM)
            Tb, _, _ = solve(T_inf, x_inf, 0.0, hm=f * HM)
            dT = max(Tb[i] - Ta[i] for i in range(21))
            share = ea[2] / ea[0]
            print(f"{f:7.2f} {x_inf:6.2f} {criterion(T_inf, x_inf, f*HM):7.3f}"
                  f" {share:10.1%} {dT:9.2f}")


def hm_threshold(T_inf=60.0, X_inf=0.01):
    """潜热项与显热项同量级的临界传质系数 (Pi=1)。"""
    return H * (T_inf - T0) / (RHO_D * LV * (X0 - X_inf))


def report(T_inf=60.0, X_inf=0.01):
    print(f"\n=== T_inf={T_inf} degC, X_inf={X_inf} kg/kg, t=1800 s ===")
    Ta, Xa, ea = solve(T_inf, X_inf, LV)
    Tb, Xb, eb = solve(T_inf, X_inf, 0.0)

    # 能量守恒自检
    q_in, sens, lat = ea
    err = abs(q_in - sens - lat) / q_in
    assert err < 0.02, f"能量不守恒，残差 {err:.2%}"
    print(f"供热 {q_in:.4g} J/m = 显热 {sens:.4g} + 潜热 {lat:.4g}"
          f"  (潜热占 {lat/q_in:.1%}, 残差 {err:.2%})")

    print("  r/cm     T(有潜热,degC)  T(无潜热,degC)   偏差K    X(有潜热)")
    for i in (0, 5, 10, 15, 20):
        print(f"  {i/10:4.1f}     {Ta[i]:8.3f}      {Tb[i]:8.3f}    {Tb[i]-Ta[i]:7.3f}"
              f"   {Xa[i]:7.4f}")
    dT = max(Tb[i] - Ta[i] for i in range(21))
    print(f"  最大温度偏差 {dT:.2f} K "
          f"(相当于空气-初始温差 {T_inf - T0:.0f} K 的 {dT/(T_inf-T0):.0%})")


if __name__ == "__main__":
    for t_inf in (50.0, 60.0, 70.0):
        report(t_inf)
    sweep()
    for t_inf in (50.0, 60.0, 70.0):
        print(f"T_inf={t_inf:.0f} degC: 潜热可忽略所需 hm < "
              f"{hm_threshold(t_inf):.2e} m/s (给定 {HM:.2e}, 为其 "
              f"{HM/hm_threshold(t_inf):.2f} 倍)")
    # 步长无关性
    a, _, _ = solve(60.0, 0.01, LV, dt=0.5)
    b, _, _ = solve(60.0, 0.01, LV, dt=0.25)
    print(f"\n步长无关性: dt=0.5 vs 0.25 s, 最大温差 "
          f"{max(abs(x-y) for x, y in zip(a, b)):.4f} K")
