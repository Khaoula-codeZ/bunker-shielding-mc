"""v5 slab comparison (run from the build folder after run_v5.sh).
1. Consistency: per-history mean (hstat.csv) vs histogram mean (must agree to ~1e-10).
2. Variance: per-history sigma vs histogram sum-of-w2 sigma (ratio > 1 = correlation of split copies).
3. TVL1 / TVLe vs NCRP 151 (37 / 33 cm) with MC uncertainty (Gaussian resampling, per-history sigma).
4. Stability 1M vs 2M and bias splitting (2M) vs analog (5M, 0-80 cm), using per-history sigma.
Writes summary_v5.csv."""
import glob, re, numpy as np
E=np.array([0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1,1.5,2,3,4,5,6,8,10])
Hc=np.array([0.061,0.83,1.05,0.81,0.64,0.55,0.51,0.53,0.61,0.89,1.20,1.80,2.38,2.93,3.44,4.38,5.20,6.90,8.60,11.1,13.4,15.5,17.6,21.6,25.6])
h=lambda e: np.exp(np.interp(np.log(np.clip(e,0.01,10)),np.log(E),np.log(Hc)))
c=(np.arange(70)+0.5)*0.1
spec=np.loadtxt("spectrum_6MV.csv",delimiter=",",skiprows=1); H0=np.sum(spec[:,1]*h(spec[:,0]))/1e4

def load_hist(folder,N):
    out={}
    for p in glob.glob(f"{folder}/bunker_h1_plane*.csv"):
        i=int(re.search(r"plane(\d+)",p).group(1))
        rows=[l for l in open(p) if not l.startswith("#")]
        a=np.array([r.split(",") for r in rows[1:]],dtype=float)[1:-1]
        out[i*10]=(np.sum(a[:,1]*h(c))/N/H0, np.sqrt(np.sum(a[:,2]*h(c)**2))/N/H0)
    return out

def load_hstat(folder):
    a=np.loadtxt(f"{folder}/hstat.csv",delimiter=",",skiprows=1,ndmin=2)
    out={}; N=a[0,2]
    for _,d,n,s,s2 in a:
        m=s/n; var=max(s2/n-m*m,0.)/(n-1)
        out[int(round(d))]=(m/H0, np.sqrt(var)/H0)
    return out,N

def tvl(r):
    d=np.array(sorted(r)); lt=np.log10([r[x][0] for x in d])
    m=(d>=50)&(d<=190); return np.interp(1,-lt,d), -1/np.polyfit(d[m],lt[m],1)[0]

def tvl_unc(r,n=2000,seed=1):
    rng=np.random.default_rng(seed); d=sorted(r); t1=[];te=[]
    for _ in range(n):
        rr={x:(max(r[x][0]+rng.normal()*r[x][1],1e-300),0) for x in d}
        a,b=tvl(rr); t1.append(a); te.append(b)
    return np.std(t1), np.std(te)

R={}; rows=[]
for f in ("v5_1M","v5_2M","v5_analog"):
    P,N=load_hstat(f); Hh=load_hist(f,N); R[f]=P
    dev=max(abs(Hh[x][0]/P[x][0]-1) for x in P if P[x][0]>0)
    ratio=[P[x][1]/Hh[x][1] for x in P if Hh[x][1]>0]
    t1,te=tvl(P); u1,ue=tvl_unc(P)
    print(f"\n== {f} (N = {N:.0f})")
    print(f"  mean consistency hstat vs histogram: max rel. diff {dev:.1e}  {'OK' if dev<1e-8 else 'CHECK'}")
    print(f"  sigma ratio per-history / histogram: {min(ratio):.2f} - {max(ratio):.2f} (median {np.median(ratio):.2f})")
    print(f"  TVL1 = {t1:.1f} ± {u1:.1f} cm ({(t1/37-1)*100:+.1f} %)   TVLe = {te:.1f} ± {ue:.1f} cm ({(te/33-1)*100:+.1f} %)")
    for x in sorted(P): rows.append((f,x,P[x][0],P[x][1],Hh[x][1]))
    print("  depth  T(H*10)      rel.sigma(per-history)")
    for x in sorted(P): print(f"  {x:4d}   {P[x][0]:.4e}   {P[x][1]/P[x][0]*100:5.2f} %")

A,B,G=R["v5_1M"],R["v5_2M"],R["v5_analog"]
print("\nStability 1M vs 2M (|diff|/sigma, per-history, must be < 2):")
for x in sorted(A):
    dz=abs(A[x][0]-B[x][0])/np.hypot(A[x][1],B[x][1]); print(f"  {x:3d} cm  {dz:4.2f}  {'OK' if dz<2 else 'CHECK'}")
print("\nBias check splitting (2M) vs analog (5M), 0-80 cm, per-history sigma:")
for x in sorted(G):
    if x>80: continue
    dz=abs(B[x][0]-G[x][0])/np.hypot(B[x][1],G[x][1]); print(f"  {x:3d} cm  ratio {B[x][0]/G[x][0]:.4f}  {dz:4.2f} sigma  {'OK' if dz<2 else 'CHECK'}")
with open("summary_v5.csv","w") as fo:
    fo.write("run,depth_cm,T_Hstar10,sigma_perhistory,sigma_histogram\n")
    for r in rows: fo.write("%s,%d,%.6e,%.6e,%.6e\n" % r)
print("\nWrote summary_v5.csv")
