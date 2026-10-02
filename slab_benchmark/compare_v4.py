"""Compares v4 runs: 1M vs 2M (stability, 2 sigma) and splitting vs analog (bias, 0-80 cm),
and recomputes TVL1/TVLe vs NCRP 151. Run from the build folder: python3 compare_v4.py"""
import glob, re, numpy as np
E=np.array([0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1,1.5,2,3,4,5,6,8,10])
Hc=np.array([0.061,0.83,1.05,0.81,0.64,0.55,0.51,0.53,0.61,0.89,1.20,1.80,2.38,2.93,3.44,4.38,5.20,6.90,8.60,11.1,13.4,15.5,17.6,21.6,25.6])
h=lambda e: np.exp(np.interp(np.log(np.clip(e,0.01,10)),np.log(E),np.log(Hc)))
c=(np.arange(70)+0.5)*0.1
spec=np.loadtxt("spectrum_6MV.csv",delimiter=",",skiprows=1); H0=np.sum(spec[:,1]*h(spec[:,0]))/1e4
def load(folder,N):
    out={}
    for p in glob.glob(f"{folder}/bunker_h1_plane*.csv"):
        i=int(re.search(r"plane(\d+)",p).group(1))
        rows=[l for l in open(p) if not l.startswith("#")]
        a=np.array([r.split(",") for r in rows[1:]],dtype=float)[1:-1]
        H=np.sum(a[:,1]*h(c))/N; s=np.sqrt(np.sum(a[:,2]*h(c)**2))/N
        out[i*10]=(H/H0,s/H0)
    return out
def tvl(r):
    d=np.array(sorted(r)); lt=np.log10([r[x][0] for x in d])
    m=(d>=50)&(d<=190); return np.interp(1,-lt,d), -1/np.polyfit(d[m],lt[m],1)[0]
A=load("v4_1M",1e6); B=load("v4_2M",2e6); G=load("analog",5e6)
for name,r in (("v4_1M",A),("v4_2M",B)):
    t1,te=tvl(r); print(f"{name}: TVL1 = {t1:.1f} cm ({(t1/37-1)*100:+.1f} %)  TVLe = {te:.1f} cm ({(te/33-1)*100:+.1f} %)")
print("\nStability 1M vs 2M (|diff|/sigma, must be < 2):")
for x in sorted(A):
    dz=abs(A[x][0]-B[x][0])/np.hypot(A[x][1],B[x][1]); print(f"  {x:3d} cm  {dz:4.2f}  {'OK' if dz<2 else 'CHECK'}")
print("\nBias check splitting (2M) vs analog (5M), 0-80 cm:")
for x in sorted(G):
    if x>80: continue
    dz=abs(B[x][0]-G[x][0])/np.hypot(B[x][1],G[x][1]); print(f"  {x:3d} cm  ratio {B[x][0]/G[x][0]:.4f}  {dz:4.2f} sigma  {'OK' if dz<2 else 'CHECK'}")
