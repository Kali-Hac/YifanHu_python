#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
☆*°☆*°(∩^o^)~━━  2018/4/19 10:43        
      (ˉ▽￣～) ~~ 一捆好葱 (*˙︶˙*)☆*°
      Fuction： Build Quadtree to save nodes √ ━━━━━☆*°☆*°
生成四叉树可以有可视化输出过程，description = True 将生成描述语句1.2.3.4.5.6
"""
import math
import matplotlib.pyplot as plt
import numpy as np
import datetime
description = visualization = False
node_num = 1

class Quadtree:
	def __init__(self, min_x, min_y, size, max_dep):
		self.min_x = min_x  # 区域最小x坐标
		self.min_y = min_y  # 区域最小y坐标
		self.size = size  # 区域大小
		self.isleaf = True  # 是否为叶子节点
		self.max_dep = max_dep  # 该树最大四分深度
		self.mass = 0  # 区域重量
		self.children = []  # 保存子树
		self.mass_center = [0.0, 0.0]  # 区域重心

	# 将区域四等分，生成子树
	def quad_divide(self):
		global axes
		# 如果已经到最大分割层次，则不再切割，放置在原区域
		if self.max_dep == 0:
			if len(self.children) == 0:
				# change
				self.children = [self]
		else:
			# 将一个区域分割为四个子区域a
			children_size = self.size / 2.0
			if visualization:
				visualize_divide([self.min_x, self.min_x + self.size], [self.min_y + children_size, self.min_y + children_size])
				visualize_divide([self.min_x + children_size,self.min_x + children_size], [self.min_y, self.min_y + self.size])
			self.isleaf = False
			self.children.append(Quadtree(self.min_x, self.min_y + children_size, children_size, self.max_dep - 1))
			self.children.append(Quadtree(self.min_x + children_size, self.min_y + children_size, children_size, self.max_dep - 1))
			self.children.append(Quadtree(self.min_x, self.min_y, children_size, self.max_dep - 1))
			self.children.append(Quadtree(self.min_x + children_size, self.min_y, children_size, self.max_dep - 1))

	# 加入新节点
	def insert_node(self, node):
		if self.min_x <= node[0] <= self.min_x + self.size and self.min_y <= node[1] <= self.min_y + self.size:

			# 1.输出吸收节点的过程描述(前) 包括吸收的层级、吸收的区域坐标、吸收区域的重心
			if description:
				print('The quad-part that will assimilate (%.2f, %.2f):' % (node[0], node[1]))
				print("Max-depth:%d  [sw(%.2f, %.2f), nw(%.2f, %.2f),ne(%.2f, %.2f),se(%.2f, %.2f)] Mess_center: (%.2f, %.2f)\n"
				      % (self.max_dep, self.min_x, self.min_y, self.min_x, self.min_y + self.size, self.min_x + self.size, self.min_y + self.size,
				          self.min_x + self.size, self.min_y, self.mass_center[0], self.mass_center[1]))

			# 先把原重心保存，如果是叶子节点的话需要做进一步单一性判定
			ori_node = [i for i in self.mass_center]

			# 递归容纳进这个节点,并改变每个区域总质量和重心位置
			self.assimilateNode(node)

			# 2.输出吸收节点的过程描述(后)
			if description:
				print('After assimilating (%.2f, %.2f):' % (node[0], node[1]))
				print("Max-depth:%d  [sw(%.2f, %.2f), nw(%.2f, %.2f),ne(%.2f, %.2f),se(%.2f, %.2f)] Mess_center: (%.2f, %.2f)\n"
				      % (self.max_dep, self.min_x, self.min_y, self.min_x, self.min_y + self.size, self.min_x + self.size, self.min_y + self.size,
				          self.min_x + self.size, self.min_y, self.mass_center[0], self.mass_center[1]))

			if self.isleaf:
				self.quad_divide()
				# 3.新插入区域的重心，也就是新节点的坐标
				if description:
					print('New Children Mess-Center: (%.2f, %.2f)\n' % (self.mass_center[0], self.mass_center[1]))

				# 需要将新节点插入到当前节点不存在的区域
				ori_update = new_update = False
				for i in self.children:
					if ori_node[0] * ori_node[1] != 0.0 and i.min_x <= ori_node[0] <= i.min_x + i.size and i.min_y <= ori_node[1] <= i.min_y + i.size:
						# 4.该区域内包含之前的节点则递归切分再插入
						if description:
							print("## Located at the same part -> Start recurisve divisions again to insert")
						i.quad_insert(ori_node)
						if description:
							print("## Initialize this part by the existed node\n")
						ori_update = True
						if i.min_x <= node[0] <= i.min_x + i.size and i.min_y <= node[1] <= i.min_y + i.size:
							if i.insert_node(node):
								return True
					if i.quad_insert(node):
						new_update = True
					if ori_update and new_update:
						return True
			else:
				for i in self.children:
					if i.insert_node(node):
						return True
		else:
			return False

	# 节点插入子树
	def quad_insert(self, node):
		global node_num
		if self.min_x <= node[0] <= self.min_x + self.size and self.min_y <= node[1] <= self.min_y + self.size:

			# 该区域容纳进这个节点,并改变总质量和重心位置
			self.assimilateNode(node)
			if visualization:
				plt.scatter(node[0], node[1], marker='o', color='b', s=1)
			node_num += 1
			# 6.成功插入的提示，包括3递归过程中先插入原来已有的节点来更新该区域的重心
			if description:
				print("## (%.2f, %.2f) Successfully insert 1/4 of above part [sw(%.2f, %.2f), nw(%.2f, %.2f),ne(%.2f, %.2f),se(%.2f, %.2f)]\n" %
				      (node[0], node[1], self.min_x, self.min_y, self.min_x, self.min_y + self.size, self.min_x + self.size, self.min_y + self.size,
				          self.min_x + self.size, self.min_y))
			return True
		else:
			return False

	# 更新区域重心
	def assimilateNode(self, node):
		self.mass_center[0] = (self.mass * self.mass_center[0] + node[0]) / (self.mass + 1)
		self.mass_center[1] = (self.mass * self.mass_center[1] + node[1]) / (self.mass + 1)
		self.mass += 1

	# 获取与该点(距离/区域)大小不超过阈值的所有点或重心
	def get_farthest_mass(self, node, threshold, pos_arr):
		# 去除本身
		if node[0] == self.mass_center[0] and node[1] == self.mass_center[1]:
			return
		# 如果该区域大小/距离 小于阈值，则没有搜索的必要，直接计算重心给予的力
		if self.size / math.sqrt((node[0] - self.mass_center[0]) ** 2 + (node[1] - self.mass_center[1]) ** 2) < threshold\
				and not self.isleaf:
			# 返回重心所在区域的节点数，重心坐标
			pos_arr.append([self.mass, self.mass_center])
			return
		# 如果是叶子节点且存在点的话，说明已经搜索到该区域的尽头，直接计算该区域重心(即唯一的点)给予的力
		elif self.isleaf and self.mass != 0:
			pos_arr.append([self.mass, self.mass_center])
		for i in self.children:
			i.get_farthest_mass(node, threshold, pos_arr)

	# 获取所有不为空的区域的重心
	def get_all_mass_center(self, centers):
		if len(self.children) == 0:
			return
		for i in self.children:
			if i.mass_center[0] * i.mass_center[1] != 0.0:
				if visualization:
					plt.scatter(i.mass_center[0], i.mass_center[1], marker="D", color="r", s=2)
				# 返回重心所在区域的大小，重点所在区域的节点数，重心坐标
				centers.append([i.size, i.mass, i.mass_center])
				i.get_all_mass_center(centers)

	def cal_electrical_forces(self, node, threshold, c, k, e_force_vector):
		# e_force_vector即为displacement => node距离重心控制在threshold以内的所有点对其施力的偏移量
		# 去除本身
		if node[0] == self.mass_center[0] and node[1] == self.mass_center[1]:
			return
		# 如果该区域大小/距离 小于阈值，则没有搜索的必要，直接计算重心给予的力
		if self.size / math.sqrt(
								(node[0] - self.mass_center[0]) ** 2 + (node[1] - self.mass_center[1]) ** 2) < threshold \
				and not self.isleaf:
			# distance = math.sqrt((self.mass_center[0] - node[0]) ** 2 + (self.mass_center[1] - node[1]) ** 2)
			distance_sqr = (self.mass_center[0] - node[0]) ** 2 + (self.mass_center[1] - node[1]) ** 2
			# fr(i,j) = -C*K**2 / ||xi - xj||
			# 若按照公式来，fr(i,j)取负，那么方向为 施力点node(j) - 受力点node(i)
			# 但我们使用受力点向量与施力点向量的差 node(i) - node(j)，所以此处fr(i,j)取正
			# 力的大小与重心的大小self.mass成正比
			e_force_vector[0] += (node[0] - self.mass_center[0]) * (c * k ** 2 / distance_sqr) * self.mass
			e_force_vector[1] += (node[1] - self.mass_center[1]) * (c * k ** 2 / distance_sqr) * self.mass
			return
		# 如果是叶子节点且存在点的话，说明已经搜索到该区域的尽头，直接计算该区域重心(即唯一的点)给予的力
		elif self.isleaf and self.mass != 0:
			distance_sqr = (self.mass_center[0] - node[0]) ** 2 + (self.mass_center[1] - node[1]) ** 2
			e_force_vector[0] += (node[0] - self.mass_center[0]) * (c * k ** 2 / distance_sqr)
			e_force_vector[1] += (node[1] - self.mass_center[1]) * (c * k ** 2 / distance_sqr)
		# 达到最大分割次数的区域，不再遍历子树
		if self.children == [self]:
			return
		for i in self.children:
			i.cal_electrical_forces(node, threshold, c, k, e_force_vector)


def visualize_divide(node_1, node_2):
	# print(node_1, node_2)
	plt.plot(node_1, node_2, linewidth=0.5)
	# plt.pause(0.1)

class BHError(Exception):
	def __init__(self, value):
		self.value = value

	def __str__(self):
		return repr(self.value)

	def cal_error(self):
		pass

	def insert_error(self):
		pass

# 测试计算斥力得到displacement的函数
def test_electrical_forces(descri=False, visual=False):
	global description, visualization
	description = descri
	visualization = visual
	now = datetime.datetime.now()
	test = Quadtree(min_x=0, min_y=0, size=10, max_dep=100)
	node_list = [[1, 2], [8, 6], [8, 9], [8, 10], [6, 4], [1, 9]]
	for node in node_list:
		test.insert_node(node)
	t = []
	for node in node_list:
		force_vector = [0.0, 0.0]
		t = []
		test.cal_electrical_forces(node=node, threshold=0.5, c=1, k=1, e_force_vector=force_vector)
		test.get_farthest_mass(node=node, threshold=0.5, pos_arr=t)
		# print(t)
		x_force_sum = y_force_sum = 0.0
		# print(node)
		for i in t:
			distance_sqr = (node[0] - i[1][0]) ** 2 + (node[1] - i[1][1]) ** 2
			x_force_sum += (node[0] - i[1][0]) * 1 / distance_sqr * i[0]
			y_force_sum += (node[1] - i[1][1]) * 1 / distance_sqr * i[0]
		test_values = [round(x_force_sum, 6), round(y_force_sum, 6)]
		cal_values = [round(force_vector[0], 6), round(force_vector[1], 6)]
		assert test_values == cal_values
		# print([x_force_sum, y_force_sum])
		# print(force_vector)
	end = datetime.datetime.now()
	print("## Function [cal_electrical_forces] PASS")
	print("### using time：%s" % str(end - now))
	if visual:
		plt.pause(10)

# 测试插入节点的函数
def test_insert_node(node_num=1000, max_dep=20, visual=False, descri=False):
	global description, visualization
	visualization = visual
	description = descri
	test = Quadtree(min_x=0, min_y=0, size=node_num, max_dep=max_dep)
	arr = np.random.randint(0, node_num, (node_num, 2))
	now = datetime.datetime.now()
	for i in range(arr.shape[0]):
		test.insert_node([arr[i, 0], arr[i, 1]])
		# if i % 100 == 0:
		# 	print(i)
	end = datetime.datetime.now()
	print("## Function [insert_node] + %d PASS" % node_num)
	print("### using time：%s" % str(end - now))
	if visual:
		plt.pause(10)


if __name__ == "__main__":
	# 开启可视化会大幅增加运行时间
	test_insert_node(node_num=10000, max_dep=20, visual=False, descri=False)
	test_electrical_forces(visual=True, descri=False)


