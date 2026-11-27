"""Synthetic clipped-start stress set from MITDB (no 118/119) and SVDB, as in PREREGISTRATION.md."""
import os,sys,json;W=os.getcwd();sys.path.insert(0,W+'/repos/tenhz-qrs/scripts');sys.path.insert(0,W+'/repos/tenhz-qrs/src')
import numpy as np,benchmark_three as b
rng=np.random.default_rng(42);OUT=W+'/data/synth-mat';os.makedirs(OUT,exist_ok=True)
def rail_segment(n,R,fs,rng):
    y=np.empty(n);i=0;s=rng.choice([-1,1]);tr=max(1,int(.02*fs))
    while i<n:
        L=max(1,int(rng.exponential(.5)*fs));y[i:i+L]=s*R;i+=L;s=-s
    # 20 ms linear transitions
    k=np.flatnonzero(np.diff(y)!=0)
    for j in k:
        a,c=max(0,j-tr//2),min(n,j+tr//2+1);y[a:c]=np.linspace(y[a],y[c-1],c-a)
    return y
def inject(x,a,e,R,tau,fs,rng):
    seg=rail_segment(e-a,R,fs,rng);last=seg[-1];x[a:e]=seg
    t=np.arange(len(x)-e)/fs;x[e:]+=last*np.exp(-t/tau)
man={};meta=[]
for db,rec,fs in b.records(W+'/data',['mitdb','svdb']):
    ecg,ann=b.load(W+'/data',db,rec);N=int(30*fs)
    for k in range(4):
        s0=int(rng.integers(0,len(ecg)-N));x=ecg[s0:s0+N].copy();a_=ann[(ann>=s0)&(ann<s0+N)]-s0
        mask=np.zeros(N,bool);art=rng.random()>=.25;info=dict(name=f'{db}_{rec}_{k}',art=bool(art))
        if art:
            A=np.percentile(np.abs(x-np.median(x)),99);m=rng.uniform(3,10);R=m*A
            D=rng.uniform(1,15);e=int(D*fs);tau=rng.uniform(.2,2)
            inject(x,0,e,R,tau,fs,rng);mask[:e]=True;info.update(clip_s=D,rail_mult=m,tau=tau)
            if rng.random()<.5:
                d=rng.uniform(.5,3);p=int(rng.integers(e,N-int(d*fs)));q=p+int(d*fs)
                inject(x,p,q,R,rng.uniform(.2,2),fs,rng);mask[p:q]=True;info.update(drop_s=d,drop_at=p/fs)
        a_=a_[~mask[a_]]
        nm=info['name'];np.savetxt(f'{OUT}/{nm}.txt',x[:,None],fmt='%.5f');np.savetxt(f'{OUT}/{nm}-ann.txt',(a_+1)[:,None],fmt='%d')
        d_=np.diff(np.r_[0,mask.astype(np.int8),0]);iv=np.c_[np.flatnonzero(d_==1)+1,np.flatnonzero(d_==-1)]
        np.savetxt(f'{OUT}/{nm}-mask.txt',iv.reshape(-1,2),fmt='%d')
        man[nm]=dict(fs=fs,status='OK',n_samples=N,n_beats=int(len(a_)));meta.append(info)
json.dump(man,open(f'{OUT}/_manifest.json','w'),indent=1);json.dump(meta,open(W+'/startup/synth_meta.json','w'))
print(len(man),'excerpts;',sum(m['art'] for m in meta),'with artifact;',sum('drop_s' in m for m in meta),'with dropout')
