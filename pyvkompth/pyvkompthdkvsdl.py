import os, sys
import pyvkompthdl, pyvkompthdk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons  # import the Slider widget
from numpy import zeros, arange, where, sort, array, abs, argmin, pi

# Function smoother :
# input : phi is a numpy ndarray that contains the
# phase lags for the particular QPO at different energies
# (could be either the full grid or just the centers of the arbitrary energy bands)
def smoother(phi):
    phi = np.mod(phi+pi, 2.*pi)-pi-phi[0]
    ps = phi.size
    nphi = zeros(ps) + phi
    iis = arange(0, ps - 1)
    pivots = where(phi[iis] * phi[iis + 1] < 0)[0]
    pivots = sort(pivots)
    si = 0
    for pivot in pivots:
        patch1 = phi[si:pivot + 1]
        corrlist = array([6.*pi, 4.*pi, 2. * pi, -2. * pi, -4. * pi, -6.*pi, 0.])
        compcases = patch1[pivot] + corrlist
        # case1, case2, case3 = patch1 + 2.0 * pi, patch1 - 2.0 * pi, patch1 + 0.0 * pi
        compref = phi[pivot + 1]
        dist = abs(compcases - compref)  # np.array([np.abs(case1[-1] - comp), np.abs(case2[-1] - comp), np.abs(case3[-1] - comp)])
        mindist = argmin(dist)
        nphi[si:pivot + 1] = nphi[si:pivot + 1] + corrlist[mindist]

    # zeropos = np.where(nphi == 0)[0]
    # nphi[zeropos] = phi[zeropos]
    return nphi



kTes = [3.,50.]
taus = [0.5,15.0]
kTss = [0.01,2.0]
etas = [0.,0.999]
Lfs = [10.0,25000.0]
qpos = [0.1, 1000.]
afs = [10.0, 250.]
dls = [1e-3, 0.1]

af = 250.
kTe = 20.
tau = 1.3
kTs = 0.2
eta = 0.5
size = 7100.
qpo_freq = 4.5
dl = 1e-1

ne = 1000
ear = np.logspace(np.log10(0.1), np.log10(50.0), ne+1)
energies = 0.5*(ear[1:]+ear[:-1])

def make_plot_rms(plot):
    global kTe, size, eta, tau, kTs, qpo_freq, af, dHext, dTe_mod_dk, dTs_mod_dk, eta_int_dk, dTe_mod_dl, dTs_mod_dl, eta_int_dl, Hexo0_out, dl
    photar_dl, photer_dl, dTe_mod_dl, dTs_mod_dl, Hexo0_out, eta_int_dl = pyvkompthdl.pyvkompthdl(ear, ne, [kTs, kTe, tau, size, eta, qpo_freq, af, dl], 1)
    dHext = (kTe * 1.6021766208e-16) * qpo_freq * 6 * np.pi * (af**2 +2 * af * size +size**2) * dl / (Hexo0_out * (3 * af**2 +3 *af*size + size**2))
    photar_dk, photer_dk, dTe_mod_dk, dTs_mod_dk, Hexo0_out, eta_int_dk = pyvkompthdk.pyvkompthdk(ear, ne, [kTs, kTe, tau, size, eta, qpo_freq, af, dHext], 1)
    plt.plot(energies, photar_dk/(ear[1:]-ear[:-1]), lw=3)
    plt.plot(energies, photar_dl/(ear[1:]-ear[:-1]), lw=3)    
    plot.set_xscale('log'); plot.set_yscale('log'); #plot.set_ylim(0.001, 0.55)
    plot.set_xlabel('Energy [keV]'); plot.set_ylabel('fractional rms')
#    plot.set_yticks([0.001,0.01,0.1,0.4])
#    plot.set_yticklabels(['0.001','0.01','0.1',''])
    plot.set_xticks([0.1,1.0,10.0])
    plot.set_xticklabels(['0.1','1','10'])
    plot.grid()
    
    return plot


def make_plot_lag(plot):
    global kTe, size, eta, tau, kTs, qpo_freq, af, dl, dHext
    photar_dl, _, _, _, _, _ = pyvkompthdl.pyvkompthdl(ear, ne, [kTs, kTe, tau, size, eta, qpo_freq, af, dl], 2)
    photar_dk, _, _, _, _, _ = pyvkompthdk.pyvkompthdk(ear, ne, [kTs, kTe, tau, size, eta, qpo_freq, af, dHext], 2)
    plt.plot(energies, photar_dk/(ear[1:]-ear[:-1]), lw=3, label='dHext')
    plt.plot(energies, photar_dl/(ear[1:]-ear[:-1]), lw=3, label='dL')    
    plot.set_xscale('log'); #plot.set_ylim(-3.5,3.5)
    plot.set_ylabel('Phase Lag [rad]')
    plot.set_xticklabels('')
    plot.set_yticks([-3.4,-np.pi,-np.pi/2,0,np.pi/2,np.pi])
    plot.set_yticklabels(['','$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'])
    plot.legend()
    plot.grid()
    return plot


def update():
    global kTe, size, eta, tau, kTs, qpo_freq, af, dl, dTe_mod_dk, dTs_mod_dk, eta_int_dk, dTe_mod_dl, dTs_mod_dl, eta_int_dl, Hexo0_out, dHext

    plt.sca(plt_rms)
    plt_rms.cla()
    plt_plot_rms = make_plot_rms(plt_rms)

    plt.sca(plt_lag)
    plt_lag.cla()
    plt_plot_lag = make_plot_lag(plt_lag)

    plt.suptitle( r'QPO={:.1f}Hz  L={:.2f}km  $kT_e$={:.2f}keV  $kT_s$={:.2f}keV  $\tau$={:.2f}  $\eta$={:.2f}  $a_f$={:.1f}km $\delta L$={:.5f} $\delta Hext$={:.5f}'.format(qpo_freq,size, kTe, kTs, tau, eta, af, dl, dHext) + '\n' + r'$\eta_{{int,dH}}$={:.1f}%  $\delta T_{{e,dH}}$={:.1f}%  $\delta T_{{s,dH}}$={:.1f}%'.format(eta_int_dk*100, dTe_mod_dk*100, dTs_mod_dk*100) + r' $\eta_{{int,dL}}$={:.1f}%  $\delta T_{{e,dL}}$={:.1f}%  $\delta T_{{s,dL}}$={:.1f}%'.format(eta_int_dl*100, dTe_mod_dl*100, dTs_mod_dl*100), y=0.99) 

    plt.draw()


def update_af(a):
    global af
    af = 10**a
    af_slider.valtext.set_text('{:.1f} km'.format(af))
    update()

def update_qpo(a):
    global qpo_freq
    qpo_freq = 10**a
    qpo_slider.valtext.set_text('{:.1f} Hz'.format(qpo_freq))
    update()

def update_kTe(a):
    global kTe
    kTe = 10**a
    kTe_slider.valtext.set_text('{:.2f} keV'.format(kTe))
    update()

def update_L(a):
    global size
    size = 10**a
    L_slider.valtext.set_text('{:.1f} km'.format(size))
    update()

def update_eta(a):
    global eta
    eta = a
    update()

def update_tau(a):
    global tau
    tau = a
    update()

def update_kTs(a):
    global kTs
    kTs = a
    update()
    
def update_dl(a):
    global dl
    dl = 10**a
    dl_slider.valtext.set_text('{:.5f}'.format(dl))
    update()


plt.figure(figsize=(10,8))

plt_rms = plt.axes([0.15, 0.20, 0.8, 0.35])
plt_lag = plt.axes([0.15, 0.55, 0.8, 0.35])

slider_L = plt.axes([0.10, 0.07, 0.3, 0.02])
slider_kTe = plt.axes([0.10, 0.01, 0.3, 0.02])
slider_kTs = plt.axes([0.10, 0.04, 0.3, 0.02])

slider_tau = plt.axes([0.60, 0.07, 0.3, 0.02])
slider_eta = plt.axes([0.60, 0.04, 0.3, 0.02])
slider_qpo = plt.axes([0.60, 0.01, 0.1, 0.02])
slider_af = plt.axes([0.8, 0.01, 0.1, 0.02])
slider_dl = plt.axes([0.65, 0.1, 0.25, 0.02])

plt.sca(plt_rms)
plt.sca(plt_lag)

af_slider = Slider(slider_af, r'$a_f$', np.log10(afs[0]), np.log10(afs[-1]), valinit=np.log10(af), valfmt='%.1f km')
af_slider.valtext.set_text('{:.1f} km'.format(10**af_slider.val))
af_slider.on_changed(update_af)

qpo_slider = Slider(slider_qpo, r'QPO', np.log10(qpos[0]), np.log10(qpos[-1]), valinit=np.log10(qpo_freq), valfmt='%.1f Hz')
qpo_slider.valtext.set_text('{:.1f} Hz'.format(10**qpo_slider.val))
qpo_slider.on_changed(update_qpo)

kTe_slider = Slider(slider_kTe, r'$kT_e$', np.log10(kTes[0]), np.log10(kTes[-1]), valinit=np.log10(kTe), valfmt='%.2f keV')
kTe_slider.valtext.set_text('{:.2f} keV'.format(10**kTe_slider.val))
kTe_slider.on_changed(update_kTe)

L_slider = Slider(slider_L, r'$L$', np.log10(Lfs[0]), np.log10(Lfs[-1]), valinit=np.log10(size), valfmt='%.1f km')
L_slider.valtext.set_text('{:.1f} km'.format(10**L_slider.val))
L_slider.on_changed(update_L)

eta_slider = Slider(slider_eta, r'$\eta$', etas[0], etas[-1], valinit=eta, valfmt='%.2f')
eta_slider.on_changed(update_eta)

tau_slider = Slider(slider_tau, r'$\tau$', taus[0], taus[-1], valinit=tau, valfmt='%.2f')
tau_slider.on_changed(update_tau)

kTs_slider = Slider(slider_kTs, r'$kT_s$', kTss[0], kTss[-1], valinit=kTs, valfmt='%.2f keV')
kTs_slider.on_changed(update_kTs)

dl_slider = Slider(slider_dl, r'$\delta L$', np.log10(dls[0]), np.log10(dls[-1]), valinit=np.log10(dl), valfmt='%.5f')
dl_slider.valtext.set_text('{:.5f}'.format(10**dl_slider.val))
dl_slider.on_changed(update_dl)

update()
plt.show()
