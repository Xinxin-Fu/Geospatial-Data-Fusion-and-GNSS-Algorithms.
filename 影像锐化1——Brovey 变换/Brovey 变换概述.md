# Brovey 变换（Brovey Transform）

## 1. 基本概念

Brovey 变换是一种常用的影像融合方法，其核心思想是将多光谱影像每个波段在总亮度中的比例作为权重，按比例将多光谱波段映射到高分辨率全色影像（PAN）上。这样可以在增强空间细节的同时，尽量保留原多光谱影像的光谱特征。  

**适用场景**：  
- 卫星影像融合（如 Landsat、SPOT、Sentinel）  
- 遥感图像可视化与解译  
- 精细地物分析和变化检测  

---

## 2. 基本步骤（文字描述）

1. **计算总亮度**  
   对每个像元，将多光谱影像所有波段的值加总得到总亮度。这一步是为了了解每个像元在整体光强中的贡献情况。  

2. **计算各波段比例**  
   将每个波段的值除以该像元的总亮度，得到该波段在总亮度中的比例。这个比例可以反映该波段在像元光谱信息中的相对重要性。  

3. **将比例映射到全色影像**  
   使用 Python 矩阵举例说明：假设有一个 3×3 的多光谱矩阵和对应的全色影像矩阵，我们直接按照比例将多光谱波段映射到全色影像上，得到融合后的高分辨率波段。  

---

## 3. Python 示例（3×3 矩阵）

```python
import numpy as np

# 定义多光谱影像矩阵（RGB 波段，3x3）
R = np.array([[50, 80, 60],
              [70, 60, 90],
              [30, 50, 40]], dtype=float)

G = np.array([[40, 60, 50],
              [60, 50, 80],
              [20, 40, 30]], dtype=float)

B = np.array([[30, 50, 40],
              [50, 40, 70],
              [10, 30, 20]], dtype=float)

# 定义全色影像矩阵
PAN = np.array([[200, 220, 210],
                [180, 200, 190],
                [150, 170, 160]], dtype=float)

# 计算总亮度
total_intensity = R + G + B

# 计算各波段比例
R_ratio = R / total_intensity
G_ratio = G / total_intensity
B_ratio = B / total_intensity

# Brovey 变换：将比例映射到全色影像
R_new = R_ratio * PAN
G_new = G_ratio * PAN
B_new = B_ratio * PAN

print("融合后的 R 波段：\n", R_new)
print("融合后的 G 波段：\n", G_new)
print("融合后的 B 波段：\n", B_new)
