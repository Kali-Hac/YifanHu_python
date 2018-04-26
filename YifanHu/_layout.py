#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
☆*°☆*°(∩^o^)~━━  2018/4/25 10:41        
      (ˉ▽￣～) ~~ 一捆好葱 (*˙︶˙*)☆*°
      Fuction：Build YifanHu layout √ ━━━━━☆*°☆*°
"""
import numpy as np
import networkx as nx
import random
import math
import datetime
from Barnes_Hut import _quadtree
import matplotlib.pyplot as plt
from gexf_transformer import *

class YifanHu:
	def __init__(self, graph, pos_array):
		assert graph is not None
		assert pos_array.shape[0] == graph.number_of_nodes() and pos_array.shape[1] == 2
		self.pos = pos_array  # np.array(shape(n, 2)) 所有节点的坐标(2D)
		self.G = graph  # 网络图 networkx graph
		self.node_num = pos_array.shape[0]  # 节点数量
		# self.fixed = fixed  # list 想要固定的节点的序号id数组
		self.Converged = False  # 是否完成迭代，跟convergenceThreshold相关
		self.delta = None  # np.array(shape=(n, n, 2)) 所有点对于所有点的位置偏移，原斥力公式将由负变正
		self.A = None  # np.array(shape=(n, n)) 邻接矩阵(对角线为该点的总度数)，矩阵计算加快速度
		self.distance = None  # np.array(shape=(n, n))  所有点与所有点的距离(对角线为0)
		self.AverageEdgeLength = None  # 图的平均边长
		self.relativeStrength = None  # electric斥力相对于spring引力的比例大小，即公式中的C default 0.2
		self.stepRatio = None  # 步长更新率的倒数 default 0.95
		self.quadTreeMaxLevel = None  # 四叉树的最大四分次数 default 20
		self.barnesHutTheta = None  # BH算法搜寻点的阈值(size / d < barnesHutTheta) default 1.2
		self.convergenceThreshold = None  # 能量收敛的比例阈值 default 1e-4
		self.adaptiveCooling = None  # 自适应退火 default True
		self.optimalDistance = None  # 最优目标距离 default relativeStrength ** (1.0 / 3) * getAverageEdgeLength()
		self.step = None  # 步长大小
		self.initialStep = None  # 初始步长 default optimalDistance / 5.0
		self.progress = None  # 正反馈系数，用于控制步长大小的更新 default 0 and progress < 5
		self.energy0 = None  # 初始能量
		self.energy = None  # 当前能量
		self.initPropertiesValues()

	# 需要调节的参数是s1(父节点对子节点的吸引力相对于spring大小),s2(子节点对父节点的排斥力相对于spring力大小)
	def getAverageEdgeLength(self):
		s1 = 3.2  # 已弃用，改成手动调节,Line 116
		s2 = -0.8  # 已弃用，改成手动调节,Line 117
		# 转换为邻接矩阵，对角线即为该点的度大小
		G = self.G.to_undirected()
		t = nx.to_numpy_matrix(self.G, weight='weight')
		self.A = nx.to_numpy_matrix(G, weight='weight')  # 已弃用，改成手动调节
		# nx.to_numpy_matrix确认邻接矩阵
		assert t.all() == self.A.all()
		self.A = np.asarray(self.A)
		# Gephi 0.9.1 JAVA源码考虑有向图，父节点加上引力，子节点减去引力
		# https://github.com/gephi/gephi/blob/master/modules/LayoutPlugin/src/main/java/org/gephi/layout/plugin/force/yifanHu/YifanHuLayout.java
		# Line 289-298
		# if (!e.getSource().equals(e.getTarget())) {
		# Node n1 = e.getSource();
		# Node n2 = e.getTarget();
		# ForceVector f1 = n1.getLayoutData();
		# ForceVector f2 = n2.getLayoutData();
		#
		# ForceVector f = getEdgeForce().calculateForce(n1, n2);
		# f1.add(f);  f即spring力大小，这里可以看成s1=1
		# f2.subtract(f); 这里可以看成s2=-1
		# }
		self.delta = np.zeros((self.node_num, self.node_num, self.pos.shape[1]), dtype=self.A.dtype)
		for i in range(self.pos.shape[1]):
			self.delta[:, :, i] = self.pos[:, i, None] - self.pos[:, i]
		self.distance = np.sqrt((self.delta ** 2).sum(axis=-1))
		# 乘以邻接矩阵，得到所有边长度的一个 n * n 矩阵， 这个矩阵的总和为所有边长和的两倍
		assert type(self.A) == type(self.distance)
		# 确认都是数组类型才能用 * 相乘，否则只能用
		# distance_x_2 = np.multiply(distance, A)
		distance_x_2 = self.distance * self.A
		edge_num = len(np.nonzero(distance_x_2)[0]) / 2
		self.AverageEdgeLength = np.sum(distance_x_2) / 2 / edge_num
		for u, v in self.G.edges():  # 已弃用，改成手动调节
			self.A[u, v] = s1  # 已弃用，改成手动调节 Line 116
			self.A[v, u] = s2  # 已弃用，改成手动调节 Line 117
		# 后面再转换
		# self.G = self.G.to_undirected()
		return self.AverageEdgeLength

	# 初始化所有参数，按照Gephi源码设置,其中比较可能需要改动的是relativeStrength，与optimalDistance相关
	def initPropertiesValues(self):
		self.progress = 0
		self.energy = 2 ** 31 - 1  # 初始能量默认为无穷大
		self.stepRatio = 0.95
		self.relativeStrength = 0.2
		self.quadTreeMaxLevel = 20
		self.barnesHutTheta = 1.2
		self.adaptiveCooling = True
		self.convergenceThreshold = 1e-4
		# Vital: 这里的最优距离要乘以初始坐标的放大倍数
		self.optimalDistance = self.relativeStrength ** (1.0 / 3) * self.getAverageEdgeLength()
		self.initialStep = self.optimalDistance / 5.0
		self.step = self.initialStep

	# 手动设置参数
	def resetPropertiesValues(self, stepRatio=0.95, relativeStrength=0.2, quadTreeMaxLevel=20,
	barnesHutTheta=1.2, convergenceThreshold=1e-4):
		assert 0.0 < stepRatio <= 1.0
		assert 0.0 < relativeStrength
		assert quadTreeMaxLevel > 10
		self.stepRatio = stepRatio
		self.relativeStrength = relativeStrength
		self.quadTreeMaxLevel = quadTreeMaxLevel
		self.barnesHutTheta = barnesHutTheta
		self.convergenceThreshold = convergenceThreshold
		self.optimalDistance = self.relativeStrength ** (1.0 / 3) * self.getAverageEdgeLength()
		self.initialStep = self.optimalDistance / 5.0
		self.step = self.initialStep

	def updateStep(self):
		if self.adaptiveCooling:
			if self.energy < self.energy0:
				self.progress += 1
				if self.progress >= 5:
					self.progress = 0
					self.step /= self.stepRatio
					# self.relativeStrength *= self.stepRatio
			else:
				self.progress = 0
				# self.step /= self.stepRatio
		else:
			self.step /= self.stepRatio

	def control_run(self):
		self.updateStep()
		print("Change ratio: %f" % (abs(self.energy - self.energy0) / self.energy))
		if abs(self.energy - self.energy0) / self.energy < self.convergenceThreshold:
			self.Converged = True

	def run_layout(self):
		min_x, min_y = min(self.pos.T[0]), min(self.pos.T[1])
		size = max(max(self.pos.T[0]) - min(self.pos.T[0]), max(self.pos.T[1]) - min(self.pos.T[1]))
		BH_tree = _quadtree.Quadtree(min_x=min_x, min_y=min_y, size=size, max_dep=self.quadTreeMaxLevel)
		for i in self.pos:
			BH_tree.insert_node(i.tolist())
		e_forces_move = []
		# 计算每一个点所受到的斥力偏移量，delta(受力点-施力点的位置偏移量)已包含在BH算法里面计算
		for i in self.pos:
			e_force_vector = [0.0, 0.0]
			BH_tree.cal_electrical_forces(node=i.tolist(), threshold=self.barnesHutTheta,
			c=self.relativeStrength, k=self.optimalDistance, e_force_vector=e_force_vector)
			e_forces_move.append(e_force_vector)
		#
		#  _fruchterman_reingold算法设置
		# self.optimalDistance = np.sqrt(100.0 / self.node_num)
		# self.distance = np.where(self.distance < 0.01, 0.01, self.distance)
		# spring_forces_move = np.transpose(np.transpose(self.delta) *
		# (self.optimalDistance * self.optimalDistance/self.distance ** 2 - self.A * self.distance / self.optimalDistance)).sum(axis=1)
		#

		electric_forces_move = np.array(e_forces_move)

		# 已弃用，改成手动调节
		# spring_forces_move = np.transpose(np.transpose(self.delta) *
		# (self.A * self.distance / self.optimalDistance)).sum(axis=1)

		# 手动调节，使用for循环来代替上面计算
		# s1代表父节点受子节点的引力相对于spring引力的大小
		# s2代表子节点受父节点的引力相对于spring引力的大小
		s1 = 1
		s2 = 1
		spring_forces_move = np.zeros(shape=(self.node_num, 2))
		for u, v in self.G.edges():
			spring_forces_move[u, :] += s1 * (self.pos[u, :] - self.pos[v, :]) * self.distance[u, v] / self.optimalDistance
			spring_forces_move[v, :] -= s2 * (self.pos[u, :] - self.pos[v, :]) * self.distance[u, v] / self.optimalDistance

		assert electric_forces_move.shape == spring_forces_move.shape

		# 此处取electric_forces_move - spring_forces_move，对应delta为受力点-施力点的位置偏移量
		# print(get_force_norm(electric_forces_move))
		# print(get_force_norm(spring_forces_move))
		# 实验一 都是引力，最大度被排斥，其余点聚集坍塌
		# displacement = -electric_forces_move + spring_forces_move
		# 实验二 都是斥力，最大度也被排斥，其余点保持一定距离
		# displacement = electric_forces_move + spring_forces_move
		# 实验三 斥力-引力
		displacement = electric_forces_move - spring_forces_move
		self.energy0 = self.energy
		max_force, self.energy = get_force_norm(displacement)

		# 分开绘制
		# self.pos += electric_forces_move * (self.step / max_force)
		# plt.clf()
		# nx.draw(self.G, self.pos, with_labels=False, node_size=1)
		# plt.pause(0.1)
		# #
		# self.pos -= spring_forces_move * (self.step / max_force)
		# plt.clf()
		# nx.draw(self.G, self.pos, with_labels=False, node_size=1)
		# plt.pause(0.1)

		assert max_force > 0
		self.pos += displacement * (self.step / max_force)
		self.control_run()
		print("Energy: %d" % self.energy)
		plt.clf()
		self.G = self.G.to_undirected()
		nx.draw(self.G, self.pos, with_labels=False, node_size=0.5, width=0.1)
		if self.Converged:
			plt.pause(100)
		plt.pause(0.1)

# 计算力矢量的大小总和及其最大值
def get_force_norm(forces):
	forces = np.sqrt((forces ** 2).sum(axis=1))
	# _fruchterman_reingold
	forces = np.where(forces < 0.01, 0.01, forces)
	return max(1.0, np.max(forces)), np.sum(forces)

# 测试计算平均边长的函数
def test_getAverageEdgeLength(node_num=10, pos_range=100):
	now = datetime.datetime.now()
	G = nx.Graph()
	for i in range(node_num):
		G.add_node(i)
	# 随机生成node的位置numpy二维数组(pos_range范围内)
	pos = np.random.randint(0, pos_range, (node_num, 2))
	# 随机加入无向边，边数设置为节点的两倍
	for i in range(node_num * 2):
		G.add_edge(random.randint(0, node_num-1), random.randint(0, node_num-1))
	distance_sums = 0.0
	for i in range(node_num):
		adj_nodes = list(G.adj[i].keys())
		for node in adj_nodes:
			distance_sums += math.sqrt((pos[i, 0] - pos[node, 0]) ** 2 + (pos[i, 1] - pos[node, 1]) ** 2)
	distance_sums /= 2 * node_num
	distance_sums = round(distance_sums, 6)
	test = YifanHu(graph=G, pos_array=pos)
	t = test.getAverageEdgeLength()
	t = round(t, 6)
	assert t == distance_sums
	end = datetime.datetime.now()
	print("## Function [getAverageEdgeLength] + %d nodes PASS" % node_num)
	print("### using time：%s" % str(end - now))
	# print(t)
	# print(distance_sums)

# 测试计算总力大小的函数
def test_get_force_norm(node_num=10000):
	t = np.random.random((100000, 2))
	forces_sum_1 = 0.0
	now = datetime.datetime.now()
	for i in range(t.shape[0]):
		forces_sum_1 += math.sqrt(t[i, 0] ** 2 + t[i, 1] ** 2)
	end = datetime.datetime.now()
	print("## ['for' iterations-calculation] using time: %s" % str(end - now))
	forces_sum_1 = round(forces_sum_1, 6)
	# print(forces_sum_1)
	now = datetime.datetime.now()
	forces_sum_2 = get_force_norm(t)
	forces_sum_2 = round(forces_sum_2, 6)
	# print(forces_sum_2)
	end = datetime.datetime.now()
	print("## [numpy-calculation] using time: %s" % str(end - now))
	print("## Function [getAverageEdgeLength] + %d nodes PASS" % node_num)
	assert forces_sum_1 == forces_sum_2

# 测试一个有多个簇的图
def test_run_layout():
	g = nx.Graph()
	for i in range(30):
		g.add_node(i)
	for i in range(1, 5):
		g.add_edge(0, i)
	for i in range(6, 10):
		g.add_edge(1, i)
	for i in range(10, 15):
		g.add_edge(2, i)
	for i in range(15, 20):
		g.add_edge(3, i)
	for i in range(20, 30):
		g.add_edge(4, i)
	# pos = np.array([[0, 0], [-1, 1], [1, 1], [-1, -1], [1, -1]], dtype=float)
	pos = np.asarray(np.random.random((g.number_of_nodes(), 2)) * 100, dtype=float)
	test = YifanHu(graph=g, pos_array=pos)
	while not test.Converged:
		test.run_layout()
	print(test.Converged)

if __name__ == "__main__":
	# test_getAverageEdgeLength(node_num=100, pos_range=100)
	# test_get_force_norm(node_num=10000)
	test_run_layout()
	#
	# g = gexf_to_networkx(r"../81214.gexf")
	# pos = np.asarray(np.random.random((g.number_of_nodes(), 2)) * 1000, dtype=float)
	# test = YifanHu(graph=g, pos_array=pos)
	# while not test.Converged:
	# 	test.run_layout()
	# print(test.Converged)