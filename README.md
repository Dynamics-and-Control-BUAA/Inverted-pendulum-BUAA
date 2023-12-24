# Inverted-pendulum-BUAA

针对一阶倒立摆和二阶倒立摆系统的稳定分析和控制

一、倒立摆控制主程序：
①一级倒立摆：main_order1.m

以一级倒立摆控制程序为例：
包含：
①主程序文件（main_order1.m）：针对一阶倒立摆的线性化模型的稳定性分析（可控性分析）、极点配置实验
②非线性动力学模型程序（IP_order1_dynamic.m）
③可视化绘图程序（IP_order1_draw.m）

二、Simulink仿真
①IP_order1_simscape.slx
②IP_order1_simscape_pre.slx（matlab低版本）
使用simulink搭建实体模型，再观察分析平衡点的稳定性，留有控制接口便于进一步的控制实验。

