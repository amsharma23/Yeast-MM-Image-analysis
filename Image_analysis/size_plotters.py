#size plotters
from numpy import genfromtxt
import matplotlib.pyplot as plt
import math

path = "/Users/amansharma/Documents/Data/test/20220531_ylb128_gal2p_fluor/Best_examples/1";
data = genfromtxt(path+"/manual.csv",delimiter=',');
l = (data[1:,-1 ]);

m_ax1 = l[0::4];
m_ax2 = l[1::4];

d_ax1 = l[2::4];
d_ax2 = l[3::4];


m_mj =[];
m_mn=[];

d_mj =[];
d_mn=[];


for i in range(len(m_ax2)):

	m_mj.append(max(m_ax1[i],m_ax2[i]));
	m_mn.append(min(m_ax1[i],m_ax2[i]));
	d_mj.append(max(d_ax1[i],d_ax2[i]));
	d_mn.append(min(d_ax1[i],d_ax2[i]));

ar_m = [math.pi * m_mj[i]*m_mn[i] for i  in range(18)];
ar_d = 	 [math.pi * d_mj[i]*d_mn[i] for i  in range(18)];

vol_m = [math.pi*m_mj[i]*m_mn[i]**2 for i in range(18)];
vol_d = [math.pi*d_mj[i]*d_mn[i]**2 for i in range(18)];


fig,ax = plt.subplots(3);
x = list(range(18));




ax[0].plot(x, d_mj,label='Daughter');
ax[0].scatter(x, d_mj,color='red');
ax[0].plot(x, m_mj,label='Mother');
ax[0].scatter(x, m_mj,color='red');

ax[0].legend();
ax[0].set_title('Major axis');

ax[1].plot(x, ar_d);
ax[1].scatter(x, ar_d,color='red');
ax[1].plot(x, ar_m);
ax[1].scatter(x, ar_m,color='red');

ax[1].set_title('Area');

ax[2].plot(x, vol_d);
ax[2].scatter(x, vol_d,color='red');
ax[2].plot(x, vol_m);
ax[2].scatter(x, vol_m,color='red');

ax[2].set_title('Volume');


plt.savefig(path+"/manual_mj.svg",format='svg');
fig.clf();
plt.close();
