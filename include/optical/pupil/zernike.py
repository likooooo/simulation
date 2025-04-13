# import numpy as np
# import matplotlib.pyplot as plt
# import math
# def zernike_radial(n, m, rho):
#     """计算Zernike径向多项式 R_n^m(rho)"""
#     R = 0
#     for k in range(0, (n - abs(m)) // 2 + 1):
#         num = (-1)**k * math.factorial(n - k)
#         denom = math.factorial(k) * math.factorial((n + abs(m))//2 - k) * math.factorial((n - abs(m))//2 - k)
#         R += num / denom * rho**(n - 2*k)
#     return R

# def zernike(n, m, rho, theta):
#     """计算Zernike多项式 Z_n^m(rho, theta)"""
#     if m >= 0:
#         return zernike_radial(n, m, rho) * np.cos(m * theta)
#     else:
#         return zernike_radial(n, -m, rho) * np.sin(-m * theta)

# # 定义光瞳网格
# N = 256  # 分辨率
# x = np.linspace(-1, 1, N)
# y = np.linspace(-1, 1, N)
# X, Y = np.meshgrid(x, y)
# rho = np.sqrt(X**2 + Y**2)  # 归一化半径
# theta = np.arctan2(Y, X)    # 方位角

# # 定义光瞳（圆形孔径）
# pupil = (rho <= 1).astype(float)  # 1 表示透光，0 表示遮挡

# # 选择Zernike系数（示例：离焦 + 像散）
# c = {
#     (2, 0): 0.5,   # 离焦 (Defocus)
#     # (2, -2): 0.3,  # 像散 (Astigmatism)
# }

# # 计算波前像差
# W = np.zeros_like(rho)
# for (n, m), coeff in c.items():
#     W += coeff * zernike(n, m, rho, theta)

# # 构建光瞳函数（含像差）
# k = 2 * np.pi / 0.5  # 波长 0.5 单位
# P = pupil * np.exp(1j * k * W)

# # 可视化
# plt.figure(figsize=(10, 4))
# plt.subplot(1, 2, 1)
# plt.imshow(np.abs(P), cmap='gray')
# plt.title("Amplitude (Pupil)")
# plt.colorbar()

# plt.subplot(1, 2, 2)
# plt.imshow(np.angle(P), cmap='hsv')
# plt.title("Phase (Wavefront Aberration)")
# plt.colorbar()
# plt.show()

# import numpy as np
# from scipy.special import legendre

# # Zernike径向多项式 R_n^m(rho)
# def R_nm(n, m, rho):
#     R = 0
#     for k in range((n - m) // 2 + 1):
#         coeff = (-1)**k * np.math.factorial(n - k) / (np.math.factorial(k) * np.math.factorial((n + m) // 2 - k) * np.math.factorial((n - m) // 2 - k))
#         R += coeff * rho**(n - 2*k)
#     return R

# # 勒让德多项式 P_n(x)
# def P_n(n, x):
#     return legendre(n)(x)

# # 对比 R_2^0(rho) 和 P_2(2rho-1)
# rho = np.linspace(0, 1, 100)
# R20 = R_nm(2, 0, rho)       # Zernike R_2^0 = 2rho^2 - 1
# P2 = P_n(2, 2*(rho + 0.5) - 1)      # 勒让德 P_2(x) = (3x^2 - 1)/2

# import matplotlib.pyplot as plt
# plt.plot(rho, R20, label='Zernike $R_2^0(\\rho)$')
# plt.plot(rho, P2, '--', label='Legendre $P_2(2\\rho-1)$')
# plt.xlabel('$\\rho$'), plt.legend(), plt.title("Comparison of Zernike and Legendre")
# plt.show()
import numpy as np
from scipy.linalg import lstsq

def zernike_radial(n, m, rho):
    """计算Zernike径向多项式 R_n^m(rho)"""
    R = 0
    for k in range(0, (n - abs(m)) // 2 + 1):
        num = (-1)**k * np.math.factorial(n - k)
        denom = (np.math.factorial(k) * 
                np.math.factorial((n + abs(m)) // 2 - k) * 
                np.math.factorial((n - abs(m)) // 2 - k))
        R += num / denom * rho**(n - 2 * k)
    return R

def generate_zernike_polynomials(n_max, rho, theta):
    """生成Zernike多项式基（确保返回2D数组）"""
    Z = []
    indices = []
    for n in range(0, n_max + 1):
        for m in range(-n, n + 1, 2):
            if m >= 0:
                Z_nm = zernike_radial(n, m, rho) * np.cos(m * theta)
            else:
                Z_nm = zernike_radial(n, -m, rho) * np.sin(-m * theta)
            Z.append(Z_nm.flatten())
            indices.append((n, m))
    Z = np.column_stack(Z)  # 转换为2D矩阵
    return Z, indices

def fit_zernike(phi_unwrapped, rho, theta, n_max=6):
    """用Zernike多项式拟合相位（修复维度问题）"""
    Z, indices = generate_zernike_polynomials(n_max, rho, theta)
    phi_flat = phi_unwrapped.flatten()
    
    # 检查维度
    if Z.ndim != 2:
        Z = Z.reshape(-1, Z.shape[-1])  # 强制转为2D
    if phi_flat.ndim != 1:
        phi_flat = phi_flat.ravel()
    
    coeffs, _, _, _ = lstsq(Z, phi_flat)
    phi_fitted = Z @ coeffs
    return phi_fitted.reshape(phi_unwrapped.shape), coeffs, indices

# 测试用例
N = 100
x = np.linspace(-1, 1, N)
y = np.linspace(-1, 1, N)
X, Y = np.meshgrid(x, y)
rho = np.sqrt(X**2 + Y**2)
theta = np.arctan2(Y, X)
pupil = (rho <= 1)

# 模拟相位（离焦 + 像散）
phi_true = 0.5 * zernike_radial(2, 0, rho) + 0.3 * zernike_radial(2, 2, rho) * np.cos(2 * theta)
phi_true *= pupil

# 拟合Zernike
phi_fitted, coeffs, indices = fit_zernike(phi_true, rho, theta, n_max=4)
print("Fitted Coefficients:")
for (n, m), c in zip(indices, coeffs):
    if abs(c) > 0.01:
        print(f"Z_{n}^{m}: {c:.3f}")