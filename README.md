# SAR, Synthetic Aperture Radar 合成孔径雷达

这是我在研究生期间学习 SAR/InSAR/PolSAR 等相关领域时编写的Matlab代码（仿真、实验）以及一些总结。毕业后已经离开这个领域了，所以也不会再有新东西更新了（更新到此结束），哈哈哈。

记得之前读研时，网上这个领域的开源代码和材料等是很少的（论文倒是不少。不过SAR这个领域本身就挺小众，而且个人感觉入门门槛还是比较高的），什么东西都要自己一点点摸索。所以我把读研期间的一些东西整理并分享出来，如果能有一点参考价值，也就是这个项目的意义了。

1. 关于 SAR 成像算法的学习和仿真，主要是两个经典算法，距离多普勒算法即 RD 算法，chirp scaling算法即 CS 算法；仿真包括点目标仿真，以及对Radarsat-1数据的成像处理。简略版本可见仓库：https://github.com/denkywu/Simulation-of-SAR-Imaging-Algorithm （仅给出了RD和CS算法的代码）。本仓库中进行了分类，给出了一些实验报告和总结；
2. 关于 InSAR 即干涉SAR的学习和仿真。主要是自己仿真了（生成回波原始数据并进行成像和干涉处理）两个场景，一个是平地场景，一个是圆锥场景。简略版本可见仓库：https://github.com/denkywu/InSAR-Simulation-and-Studies 。本仓库中进行了分类，给出了一些实验报告和总结；
3. 关于 PolSAR 即极化SAR，我做的工作主要是极化定标。详见仓库：https://github.com/denkywu/PolSAR-Calibration 。由于极化处理中涉及到过多的实际数据处理和分析，不便于分享，所以也没有更多内容了。

由于本人水平有限，文中不免有错误，还望各位海涵。

谢谢！
