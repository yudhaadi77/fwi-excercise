"""Render fwi_frames.npz into the 1080x1350 time-lapse video FWI_frequency_timelapse.mp4 (needs ffmpeg)."""
import numpy as np, matplotlib, os, subprocess
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import font_manager as fm
for f in ['Poppins-Regular.ttf','Poppins-Bold.ttf','Poppins-SemiBold.ttf']:
    p='/usr/share/fonts/truetype/google-fonts/'+f
    if os.path.exists(p): fm.fontManager.addfont(p)
plt.rcParams.update({'font.family':['Poppins','DejaVu Sans'],'font.size':13})
D = np.load('fwi_frames.npz'); S=D['snaps']; L=D['labels']; E=D['errs']; vt=D['vt']; v0=D['v0']; img=D['img']; water=D['water']
nx,nz = vt.shape; dx=25; ext=[0,(nx-1)*dx/1000,(nz-1)*dx/1000,0]
BLUE='#0038BF'; YEL='#FFEC3D'; GREY='#5B6B78'
cm = plt.get_cmap('RdBu_r')
g = img.copy(); g[water]=0; c=np.percentile(np.abs(g[~water]),98); g=np.clip(g/c,-1,1)
gray = plt.get_cmap('gray_r')((g.T+1)/2)[...,:3]*0.85+0.15*1  # light gray image
gray = 0.15+0.85*gray
def overlay(dv):
    a = np.clip(np.abs(dv)/12.0,0,1)**0.7*0.8    # opacity grows with |update|
    col = cm((np.clip(dv,-20,20)+20)/40)[...,:3]
    out = gray*(1-a[...,None]) + col*a[...,None]
    out[water.T] = [1,1,1]
    return out
dtrue = ((vt-v0)/v0*100).T
os.makedirs('frames',exist_ok=True)
# timeline: hold start, then between snapshots interpolate
seq=[]
seq += [(0,0.0)]*24
for i in range(len(S)-1):
    nsteps = 12
    if L[i+1][0]!=L[i][0] and L[i][0]!=0: seq += [(i,0.0)]*18   # pause at band change
    for k in range(nsteps): seq.append((i,k/nsteps))
seq += [(len(S)-1,0.0)]*60
bands=[3,5,8]
def frame(n,i,t):
    v = S[i] if t==0 else (1-t)*S[i]+t*S[i+1]
    j = i if t<0.5 else i+1
    f0,it = L[j]; err = (1-t)*E[i]+t*E[min(i+1,len(E)-1)] if t>0 else E[i]
    dv = ((v-v0)/v0*100).T
    fig = plt.figure(figsize=(10.8,13.5),dpi=100); fig.patch.set_facecolor('white')
    fig.text(0.06,0.955,'How FWI builds the model,',fontsize=30,fontweight='bold',color=BLUE)
    fig.text(0.06,0.92,'one frequency band at a time',fontsize=30,fontweight='bold',color=BLUE)
    fig.text(0.06,0.893,'Marmousi section · Devito acoustic FWI · 16 shots · Ricker 3, 5, 8 Hz, 12 updates per band',fontsize=13,color=GREY)
    # band chips
    for b,(bf) in enumerate(bands):
        x=0.06+b*0.16; act = (f0==bf); done = f0>bf
        fc = BLUE if act else ('#DCE7FB' if done else '#EEF1F6'); tc = 'white' if act else (BLUE if done else '#9AA3B5')
        fig.patches.append(matplotlib.patches.FancyBboxPatch((x,0.835),0.14,0.038,boxstyle='round,pad=0.004,rounding_size=0.01',transform=fig.transFigure,fc=fc,ec='none'))
        fig.text(x+0.07,0.854,f'{bf} Hz',ha='center',va='center',fontsize=17,fontweight='bold',color=tc)
    fig.text(0.56,0.854,'start model' if f0==0 else f'update {it} / 12',va='center',fontsize=17,color=BLUE,fontweight='bold')
    fig.text(0.80,0.854,f'error {err:.1f} %',va='center',fontsize=17,color='#D62828',fontweight='bold')
    ax1 = fig.add_axes([0.08,0.50,0.80,0.30]); ax2 = fig.add_axes([0.08,0.13,0.80,0.30])
    ax1.imshow(overlay(dv),extent=ext,aspect='auto',interpolation='bilinear'); ax2.imshow(overlay(dtrue),extent=ext,aspect='auto',interpolation='bilinear')
    for ax,tt in [(ax1,'FWI velocity update (what the inversion found)'),(ax2,'True update needed (known only in a synthetic test)')]:
        ax.set_title(tt,loc='left',fontsize=16,fontweight='bold',color=BLUE,pad=10)
        ax.set_ylabel('depth (km)'); ax.set_ylim(ext[2],0.2); ax.tick_params(colors=GREY)
        for s in ax.spines.values(): s.set_color('#9FB2D6')
    ax2.set_xlabel('distance (km)')
    cax = fig.add_axes([0.90,0.13,0.018,0.67])
    sm = matplotlib.cm.ScalarMappable(cmap=cm,norm=matplotlib.colors.Normalize(-20,20)); cb=fig.colorbar(sm,cax=cax)
    cb.set_label('velocity update (%)',color=GREY); cb.outline.set_visible(False)
    fig.text(0.08,0.075,'Colour: change from the starting model, over a grey 8 Hz RTM image (final model).',fontsize=12.5,color=GREY)
    fig.text(0.08,0.055,'Low frequencies fix the broad velocity; higher bands add layers and faults.',fontsize=12.5,color=GREY)
    fig.text(0.08,0.025,'Yudha Adi Putra · KFUPM Geophysics',fontsize=12.5,color=BLUE,fontweight='bold')
    fig.savefig(f'frames/f{n:04d}.png'); plt.close(fig)
import sys
only = int(sys.argv[1]) if len(sys.argv)>1 else None
if only is not None:
    i,t = seq[only]; frame(only,i,t)
else:
    for n,(i,t) in enumerate(seq): frame(n,i,t)
    subprocess.run(['ffmpeg','-y','-framerate','24','-i','frames/f%04d.png','-c:v','libx264','-pix_fmt','yuv420p','-crf','18','FWI_frequency_timelapse.mp4'],check=True)
print(len(seq))
