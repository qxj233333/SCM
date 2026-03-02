# SI-BDF1 + SI-BDF2 MATLAB 实现（简洁版）

这个仓库提供了你给出的方程组的一个**清晰、易改**的 MATLAB 实现：

- `si_bdf1_bootstrap.m`：先用 SI-BDF1 跑一个大步长里的若干子步，得到 `psi^1`。
- `si_bdf2_solver.m`：主时间推进器（SI-BDF2）。
- `run_demo.m`：一键示例运行。

## 代码结构

- `run_demo.m`：设置参数、构造初值、调用求解器。
- `si_bdf2_solver.m`：完整 SI-BDF2 流程。
- `si_bdf1_bootstrap.m`：BDF2 的启动步骤。
- `build_fft_ops.m`：FFT 波数网格。
- `spectral_utils.m`：梯度 / 拉普拉斯（周期边界，谱方法）。
- `compute_velocity.m`：按公式计算速度场 `u`。
- `compute_pressure.m`：解压力泊松方程（零均值约束）。
- `solve_poisson_mean_zero.m`：周期 Poisson 求解器。
- `solve_helmholtz_fft.m`：周期 Helmholtz 求解器。
- `fixed_point_implicit_potential.m`：势能项隐式的简单不动点迭代。
- `normalize_pair.m`：归一化投影（`|psi1|^2+|psi2|^2=1`）。
- `dot_u_grad_psi.m` / `divergence_field.m`：常用算子。

## 使用方式

在 MATLAB 命令行执行：

```matlab
run_demo
```

运行后会输出：

- 最终时刻
- 归一化约束的均值检查
- 保存的时间步数量

## 说明

- 空间离散：3D FFT（天然支持周期边界）。
- 时间离散：
  - BDF2 处理时间导数；
  - 对流项 AB2 显式；
  - 色散项隐式；
  - 势能项采用“系数显式 + 波函数隐式”的思路，并用不动点迭代实现。
- 每一步后做两次投影：
  1. 归一化投影；
  2. 相位投影（通过解 `phi` 的 Poisson 方程）实现不可压约束。

如果你后续希望，我可以再帮你加：

1. 参数扫描脚本（批量跑不同 `epsilon`），
2. 涡结构可视化（isosurface），
3. 结果自动保存为 `.mat` 文件。
