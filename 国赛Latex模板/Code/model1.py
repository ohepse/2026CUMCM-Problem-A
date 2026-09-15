import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad
import matplotlib.pyplot as plt

# 解决Matplotlib中文显示问题
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

def module1_statics():
    # 系统基础参数初始化
    W, H, L, R, g, EA = 600.0, 15.5, 603.5, 1.2, 9.8, 1.5e7
    w0, xp = 80.0, 300.0
    
    # 1. 空置状态求解 [a, c1, c2]
    def eq_empty(vars):
        a, c1, c2 = vars
        eq1 = a * np.cosh(-c1 / a) + c2
        eq2 = a * np.cosh((W - c1) / a) + c2 - H
        # 被积函数：计算受自重拉伸后的形变积分
        int_val, _ = quad(lambda x: np.cosh((x - c1) / a) / (1 + a * R * g * np.cosh((x - c1) / a) / EA), 0, W)
        eq3 = int_val - L
        return [eq1, eq2, eq3]

    sol_empty = fsolve(eq_empty, [1000.0, W / 2, 0.0])
    
    # 2. 载人状态求解 [a_prime, c1L, c2L, c1R, c2R]
    def eq_loaded(vars):
        a_p, c1L, c2L, c1R, c2R = vars
        eq1 = a_p * np.cosh(-c1L / a_p) + c2L
        eq2 = a_p * np.cosh((W - c1R) / a_p) + c2R - H
        eq3 = (a_p * np.cosh((xp - c1L) / a_p) + c2L) - (a_p * np.cosh((xp - c1R) / a_p) + c2R)
        eq4 = a_p * R * g * (np.sinh((c1R - xp) / a_p) + np.sinh((xp - c1L) / a_p)) - w0 * g
        intL, _ = quad(lambda x: np.cosh((x - c1L) / a_p) / (1 + a_p * R * g * np.cosh((x - c1L) / a_p) / EA), 0, xp)
        intR, _ = quad(lambda x: np.cosh((x - c1R) / a_p) / (1 + a_p * R * g * np.cosh((x - c1R) / a_p) / EA), xp, W)
        eq5 = intL + intR - L
        return [eq1, eq2, eq3, eq4, eq5]

    sol_loaded = fsolve(eq_loaded, [1000.0, W / 2, 0.0, W / 2, 0.0])
    
    # --- 绘图逻辑 ---
    x_plot = np.linspace(0, W, 500)
    y_empty = sol_empty[0] * np.cosh((x_plot - sol_empty[1]) / sol_empty[0]) + sol_empty[2]
    
    x_L = np.linspace(0, xp, 250)
    x_R = np.linspace(xp, W, 250)
    y_L = sol_loaded[0] * np.cosh((x_L - sol_loaded[1]) / sol_loaded[0]) + sol_loaded[2]
    y_R = sol_loaded[0] * np.cosh((x_R - sol_loaded[3]) / sol_loaded[0]) + sol_loaded[4]
    
    plt.figure(figsize=(10, 5))
    plt.plot(x_plot, y_empty, 'b--', linewidth=1.5, label='空置状态理想悬链线')
    plt.plot(np.concatenate([x_L, x_R]), np.concatenate([y_L, y_R]), 'r-', linewidth=2, label='载人状态非对称耦合下陷')
    y_xp = sol_loaded[0] * np.cosh((xp - sol_loaded[1]) / sol_loaded[0]) + sol_loaded[2]
    plt.plot(xp, y_xp, 'ko', markersize=8, markerfacecolor='y', label='游客位置')
    
    plt.gca().invert_yaxis() # y轴向下为正
    plt.legend(loc='best')
    plt.title('溜索系统静力学空间形态')
    plt.xlabel('水平跨度 X (m)')
    plt.ylabel('下垂高度 Y (m)')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    module1_statics()