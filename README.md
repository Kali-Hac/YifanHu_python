# YifanHu_python
(Barnes_Hut算法)Quadtree和YifanHu layout的python实现 

**Implementation of the YifanHu Layout for representation of huge social network in Python3.5**

It is based on the article of [Efficient, High-Quality
Force-Directed Graph
Drawing](http://www.mathematica-journal.com/issue/v10i1/contents/graph_draw/graph_draw.pdf)

### The visualization and description are available in this project.You can also take a look at [my blog](https://yikunhaocong.com/)

## USAGE
```
from YifanHu import _layout
lo = _layout.YifanHu(graph=g, pos_array=pos)
while not lo.Converged:
	lo.run_layout()
```

# 第一个版本实验：
## s1=1, s2=1
- 运行前随机化坐标
![](result/origin_1.png)
- 生成布局
![](result/result_1.png)

## s1=10, s2=1
- 运行前随机化坐标
![](result/101_origin.png)
- 生成布局
![](result/101_result.png)


## s1=1, s2=10
- 运行前随机化坐标
![](result/110_origin.png)
- 生成布局
![](result/110_result.png)

## Barnes-Hut Visualization
- Insert Nodes 
![](http://op72m4y17.bkt.clouddn.com/Figure_1.png)
![](http://op72m4y17.bkt.clouddn.com/Figure_1-1.png)
![](http://op72m4y17.bkt.clouddn.com/Figure_1-2.png)
![](http://op72m4y17.bkt.clouddn.com/Figure_1-3.png)
![](http://op72m4y17.bkt.clouddn.com/Figure_1-4.png)
![](http://op72m4y17.bkt.clouddn.com/Figure_1-5.png)
![](http://op72m4y17.bkt.clouddn.com/Figure_1-6.png)
![](http://op72m4y17.bkt.clouddn.com/Figure_1-7.png)
![](http://op72m4y17.bkt.clouddn.com/Figure_1-8.png)
![](http://op72m4y17.bkt.clouddn.com/Figure_1-9.png)

- Get Mass Centers
![](http://op72m4y17.bkt.clouddn.com/node_100.png)
![](http://op72m4y17.bkt.clouddn.com/node_100_mass.png)

### Effectiveness of Model
![](http://op72m4y17.bkt.clouddn.com/eff_1.png)

### Need to optimize
See in my blog: [https://yikunhaocong.com/2018/04/29/barnes-hut/](https://yikunhaocong.com/2018/04/29/barnes-hut/)
![](http://op72m4y17.bkt.clouddn.com/eff_2.png)