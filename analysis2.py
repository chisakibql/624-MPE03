import matplotlib.pyplot as plt
import numpy as np
import math


def LeastSquareMethod(xlst, ylst):
    n = len(xlst)
    avgx = sum(xlst)/n
    avgy = sum(ylst)/n
    sumxy = sum([xlst[i]*ylst[i] for i in range(n)])
    sumx2 = sum([xi*xi for xi in xlst])
    sumy2 = sum([yi*yi for yi in ylst])
    
    k = (sumxy/n - avgx * avgy) / (sumx2/n - avgx * avgx)
    b = avgy - k * avgx
    r = (sumxy - n*avgx*avgy)/math.sqrt((sumx2-n*avgx*avgx)*(sumy2-n*avgy*avgy))
    return k,b,r


felst = []
nitrolst = []
fepeaklst = []
nitropeaklst = []

fp = open("Exp/SCA SWEEP.asc","r")
for i in range(512):
    s = fp.readline().split()
    felst.append(int(s[1]))
fp.close()

plt.plot([i for i in range(512)], felst)
plt.xlabel("x")
plt.ylabel("Counts")
plt.savefig("scaspec.png")
plt.savefig("scaspec.eps")
#plt.show()
plt.clf()


felst = []
nitrolst = []
fepeaklst = []
nitropeaklst = []
fp = open("Exp/FE ALPHA.asc","r")
for i in range(512):
    s = fp.readline().split()
    felst.append(int(s[1]))
fp.close()

plt.plot([i for i in range(512)], felst)
plt.xlabel("x")
plt.ylabel("Counts")
plt.savefig("feaspec.png")
plt.savefig("feaspec.eps")
#plt.show()
plt.clf()