1. 关于 RD算法 的一些内容，还可以再参考目录“
SAR-Synthetic-Aperture-Radar/2-InSAR干涉SAR-人造场景仿真/RD算法，并采用二维频域相位相乘实现SRC/
”下的内容，地址 https://github.com/denkywu/SAR-Synthetic-Aperture-Radar/tree/master/2-InSAR%E5%B9%B2%E6%B6%89SAR-%E4%BA%BA%E9%80%A0%E5%9C%BA%E6%99%AF%E4%BB%BF%E7%9C%9F/RD%E7%AE%97%E6%B3%95%EF%BC%8C%E5%B9%B6%E9%87%87%E7%94%A8%E4%BA%8C%E7%BB%B4%E9%A2%91%E5%9F%9F%E7%9B%B8%E4%BD%8D%E7%9B%B8%E4%B9%98%E5%AE%9E%E7%8E%B0SRC 。

2. 原因如下：在做InSAR仿真时，固然需要先对场景进行成像（而后才是干涉处理），上述目录是我选择RD算法进行成像再处理的内容。由于InSAR仿真对成像的细节要求很高（峰值位置、保相性等），所以当我在做InSAR中的成像并基于RD算法时，就发现了以往RD算法（我的代码实现）中一些没注意到的点。由此进行了较深入的仿真和论证。上述目录下的成像文件是 RDA_imaging2_v3.m 和 RDA_imaging2_v4.m，且还有几个pdf总结供大家参考，《2015.01.12.一维LFM脉冲压缩峰值位置问题及解决.pdf》，《2015.01.12.一维LFM脉冲压缩结果的保相性问题及解决.pdf》，《2015.01.14. 在二维频域相位相乘实现SRC，RDA公式推导.pdf》，《2015.01.15.对于RDA算法的改进.pdf》（此改进并不是什么算法层面的改进，是对我以往代码实现中的一些细节进行了改进，或者说修正）。

水平有限，仅供参考。
谢谢！
