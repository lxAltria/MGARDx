import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

SMALL_SIZE = 12
MEDIUM_SIZE = 20
BIGGER_SIZE = 32

plt.rc('font', size=24)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=28)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

coeff = [0.2, 0.8, 0.6, 1.2, 0.4]
phi = [0, 1, 0]
idx1 = [0, 1, 2, 3, 4]
idx2 = [0, 2, 4]



fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10,6))


p1, = ax1.plot(idx1, coeff, color = 'gray', marker = 'o', linestyle='-', linewidth=2, label='Function of interpolation difference')
p2, = ax1.plot(idx2, phi, color = 'black', linestyle='--', label='Nodal basis function')

ax1.annotate('', xy=(0, 0), xytext=(0, coeff[0]), arrowprops=dict(arrowstyle="|-|", color='blue', linewidth = 1))
ax1.annotate('$c_{2i-2}$', xy=(0, 0.5*(coeff[0])), xycoords='data', xytext=(5, 0), textcoords='offset points')
ax1.annotate('', xy=(2 - 0.1, 0), xytext=(2 - 0.1, coeff[2]), arrowprops=dict(arrowstyle="|-|", color='blue', linewidth = 1))
ax1.annotate('$c_{2i}$', xy=(2-0.1, 0.5*(coeff[2])), xycoords='data', xytext=(-15, 0), textcoords='offset points')
ax1.annotate('', xy=(4, 0), xytext=(4, coeff[4]), arrowprops=dict(arrowstyle="|-|", color='blue', linewidth = 1))
ax1.annotate('$c_{2i+2}$', xy=(4, 0.5*(coeff[4])), xycoords='data', xytext=(5, 0), textcoords='offset points')

ax1.annotate('', xy=(1, 0), xytext=(1, coeff[1]), arrowprops=dict(arrowstyle="|-|", color='green', linewidth = 1))
ax1.annotate('$c_{2i-1}$', xy=(1, 0.5*(coeff[1])), xycoords='data', xytext=(5, 0), textcoords='offset points')
ax1.annotate('', xy=(3, 0), xytext=(3, coeff[3]), arrowprops=dict(arrowstyle="|-|", color='green', linewidth = 1))
ax1.annotate('$c_{2i+1}$', xy=(3, 0.5*(coeff[3])), xycoords='data', xytext=(5, 0), textcoords='offset points')
ax1.annotate('', xy=(2, 0), xytext=(2, 1), arrowprops=dict(arrowstyle="|-|", color='green', linewidth = 1, ls='--'))
ax1.annotate('$1$', xy=(2, 0.5), xycoords='data', xytext=(5, 0), textcoords='offset points')

ax1.annotate('', xy=(0, 0), xytext=(1, 0), arrowprops=dict(arrowstyle="|-|", color='blue', linewidth = 1, ls='--'))
ax1.annotate('$h_l$', xy=(0.4, 0.02), xycoords='data', xytext=(5, 0), textcoords='offset points')
ax1.annotate('', xy=(1, 0), xytext=(2, 0), arrowprops=dict(arrowstyle="|-|", color='blue', linewidth = 1, ls='--'))
ax1.annotate('$h_l$', xy=(1.4, 0.02), xycoords='data', xytext=(5, 0), textcoords='offset points')
ax1.annotate('', xy=(2, 0), xytext=(3, 0), arrowprops=dict(arrowstyle="|-|", color='blue', linewidth = 1, ls='--'))
ax1.annotate('$h_l$', xy=(2.4, 0.02), xycoords='data', xytext=(5, 0), textcoords='offset points')
ax1.annotate('', xy=(3, 0), xytext=(4, 0), arrowprops=dict(arrowstyle="|-|", color='blue', linewidth = 1, ls='--'))
ax1.annotate('$h_l$', xy=(3.4, 0.02), xycoords='data', xytext=(5, 0), textcoords='offset points')

# ax1.annotate(r'Coefficient', xy=(1, 0.2+0.5*(coeff[1])), xytext=(60,25), 
#          textcoords='offset points', ha='center', va='bottom',color='blue',
#          bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
#          arrowprops=dict(arrowstyle='->', linestyle = '-', connectionstyle='arc3,rad=-0.3', 
#                             color='g'))

# ax1.annotate(r'Coefficient', xy=(3, 0.5), xytext=(60,25), 
#          textcoords='offset points', ha='center', va='bottom',color='blue',
#          bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
#          arrowprops=dict(arrowstyle='->', linestyle = '-', connectionstyle='arc3,rad=0.3', 
#                             color='g'))
ax1.get_yaxis().set_visible(False)


# ax1.annotate('', xy=(0, qu[0]), xytext=(0, after_correction[0]), arrowprops=dict(arrowstyle="<|-", color='blue', mutation_scale=15, linewidth = 2))
# ax1.annotate('', xy=(2, qu[2]), xytext=(2, after_correction[1]), arrowprops=dict(arrowstyle="<|-", color='blue', mutation_scale=15))
# ax1.annotate('', xy=(4, qu[4]), xytext=(4, after_correction[2]), arrowprops=dict(arrowstyle="<|-", color='blue', mutation_scale=15))

# ax1.annotate(r'Add Correction', xy=(0, 0.5*(qu[0] + after_correction[0])), xytext=(60,10), 
#          textcoords='offset points', ha='center', va='bottom',color='blue',
#          bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
#          arrowprops=dict(arrowstyle='->', linestyle = '-', connectionstyle='arc3,rad=-0.3', 
#                             color='b'))
# ax1.annotate(r'Add Correction', xy=(2, 0.5*(qu[2] + after_correction[1])), xytext=(-80,15), 
#          textcoords='offset points', ha='center', va='bottom',color='blue',
#          bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
#          arrowprops=dict(arrowstyle='->', linestyle = '-', connectionstyle='arc3,rad=-0.3', 
#                             color='b'))

# ax1.annotate(r'Add Correction', xy=(4, 0.5*(qu[4] + after_correction[2])), xytext=(-80,15), 
#          textcoords='offset points', ha='center', va='bottom',color='blue',
#          bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
#          arrowprops=dict(arrowstyle='->', linestyle = '-', connectionstyle='arc3,rad=-0.3', 
#                             color='b'))
plt.ylim(top=1.85)
ax1.set_xticks(idx1)
ax1.set_xticklabels(['$x^{2i-2}$', '$x^{2i-1}$','$x^{2i}$','$x^{2i+1}$','$x^{2i+2}$',])
ax1.tick_params(axis=u'both', which=u'both',length=0)
#ax1.set_xlabel("x")
#ax1.set_ylabel("y")
#ax1.legend(tuple([p1, p2, p3]), ['$Q_2u$ (Original Data)', 'Linear Interpolation', '$Q_1u$'])
ax1.legend(loc='upper left', bbox_to_anchor= (0, 1.02))
plt.tight_layout()
plt.box(False)
plt.savefig('load_vector_md.pdf', bbox_inches='tight')

