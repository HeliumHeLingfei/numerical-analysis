# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt

x0 = []
y0 = []
x1 = []
y1 = []
x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
y = [-4.467, -0.452, 0.551, 0.048, -0.447, 0.549, 4.552]
count = 0
with open("data2.txt") as f0:
	with open("data.txt")as f1:
		for i in range(1, 4):
			count = 0
			while count < 250:
				s = f0.readline().split()
				x0.append(float(s[0]))
				y0.append(float(s[1]))
				t = f1.readline().split()
				x1.append(float(t[0]))
				y1.append(float(t[1]))
				count += 1
			plt.figure(i)
			plt.plot(x0, y0, 'r')
			plt.plot(x1, y1, 'g')
			plt.plot(x, y, '+', 'b')
			savename = str(i) + ".png"
			plt.savefig(savename)
			x0 = []
			y0 = []
			x1 = []
			y1 = []
		plt.show()
