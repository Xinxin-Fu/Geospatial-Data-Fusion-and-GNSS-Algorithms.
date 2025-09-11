# RTKLIB 常用坐标转换

在 RTKLIB 中，坐标转换主要涉及 **ECEF ↔ 地理经纬度（WGS84）** 和 **ECEF ↔ 本地 ENU/NED**。RTKLIB 已经提供了相应的函数，通常位于 `coord.c` 文件中。

---

## 1. ECEF ↔ 经纬度（WGS84）

### 函数

```c
// xyz2blh - ECEF (X,Y,Z) -> 地理坐标 (lat, lon, height)
extern void xyz2blh(const double *r, double *pos);

// blh2xyz - 地理坐标 (lat, lon, height) -> ECEF (X,Y,Z)
extern void blh2xyz(const double *pos, double *r);
```
## 示例

```c
#include "rtklib.h"
#include <stdio.h>

int main() {
    double xyz[3] = { 3875000.0, 332500.0, 5029000.0 }; // ECEF 坐标
    double blh[3]; // lat, lon, height

    // ECEF -> 地理坐标
    xyz2blh(xyz, blh);
    printf("Latitude: %f deg, Longitude: %f deg, Height: %f m\n",
           blh[0]*R2D, blh[1]*R2D, blh[2]);

    // 地理坐标 -> ECEF
    double xyz2[3];
    blh2xyz(blh, xyz2);
    printf("X: %f, Y: %f, Z: %f\n", xyz2[0], xyz2[1], xyz2[2]);

    return 0;
}
```

## 2. ECEF ↔ ENU（本地坐标系）

##  函数

```c
// ecef2pos - ECEF -> 地理坐标
void ecef2pos(const double *r, double *pos);

// ecef2enu - ECEF -> ENU
void ecef2enu(const double *r, const double *pos, double *enu);

// enu2ecef - ENU -> ECEF
void enu2ecef(const double *enu, const double *pos, double *r);
```
## 示例

```c
#include "rtklib.h"
#include <stdio.h>

int main() {
    double ref_xyz[3] = {3875000.0, 332500.0, 5029000.0}; // 接收机参考点
    double r[3] = {3875005.0, 332510.0, 5029002.0};       // 卫星或其他点
    double enu[3];

    // ECEF -> ENU
    ecef2enu(r, ref_xyz, enu);
    printf("E: %f, N: %f, U: %f\n", enu[0], enu[1], enu[2]);

    // ENU -> ECEF
    double r2[3];
    enu2ecef(enu, ref_xyz, r2);
    printf("X: %f, Y: %f, Z: %f\n", r2[0], r2[1], r2[2]);

    return 0;
}
```






