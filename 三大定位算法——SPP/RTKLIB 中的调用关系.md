# SPP 定位在 RTKLIB 中的调用关系
# # 主要 RTKLIB 函数调用关系图如下图
readobs() --> satpos() --> build_obs_equation() --> lsq() --> outsol()
