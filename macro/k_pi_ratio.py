import matplotlib.pyplot as plt

# 데이터
mom = [645, 665, 685, 715, 735, 755, 790, 814, 842, 870, 933]
kbeam = [8, 19, 27, 30, 33, 36, 33, 36, 37, 36, 35]
pibeam = [320, 284, 360, 250, 200, 180, 120, 112, 89, 88, 115]

# ratio 계산
ratio = [k/p for k, p in zip(kbeam, pibeam)]

# figure 생성
fig, ax1 = plt.subplots()

# 왼쪽 y축 (kbeam)
color = 'tab:blue'
ax1.set_xlabel(r'Beam Momentum (MeV/$c$)')
ax1.set_ylabel(r'$K^-$ [k/spill]',color=color)
ax1.plot(mom, kbeam, marker='o', color=color, linewidth=2, label='K beam')
ax1.tick_params(axis='y', labelcolor=color)

# 오른쪽 y축 (ratio)
ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel(r'$K^-$ / $\pi^-$ ratio',color=color)
ax2.plot(mom, ratio, marker='s', linestyle='--', color=color, linewidth=2, label='Ratio')
ax2.tick_params(axis='y', labelcolor=color)

# 그래프 제목 & 그리드
#fig.suptitle('K Beam Count & K/Pi Ratio vs Momentum')
ax1.grid(True, linestyle='--', alpha=0.5)

# 출력
plt.tight_layout()
plt.savefig("k_pi_ratio.pdf")
plt.show()
plt.close()
