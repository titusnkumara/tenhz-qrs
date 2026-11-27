import os,sys,warnings;warnings.filterwarnings('ignore')
W=os.getcwd();sys.path.insert(0,W+'/startup/pkg');sys.path.insert(0,W+'/repos/tenhz-qrs/scripts')
import numpy as np,pandas as pd
import tenhz;from tenhz import detectors as D
import benchmark_three as b
assert 'startup/pkg' in tenhz.__file__
VAR={'V0':dict(init='first',wd=None,rej=0.0,reinit=False),'V1':dict(init='robust',wd=None,rej=0.0,reinit=False),'V2':dict(init='robust',wd=5.0,rej=0.0,reinit=False),'V3':dict(init='robust',wd=5.0,rej=20.0,reinit=False),'V4':dict(init='robust',wd=5.0,rej=20.0,reinit=True)}
_cache={}
def run(job):
    db,rec,fs,methods=job
    ecg,ann=b.load(W+'/data',db,rec);mk=b.load_mask(W+'/data',db,rec,len(ecg));out=[]
    for v,m in methods:
        D.VARIANT.update(VAR[v]);d=np.asarray(tenhz.qrs_locations(ecg,fs,m),np.int64)
        if mk is not None: d=d[~mk[np.clip(d,0,len(ecg)-1)]]
        tp,fp,fn,_=b.match(d,ann,b.mround(.15*fs));out.append((db,rec,v,m,tp,fp,fn))
    return out
if __name__=='__main__':
    import multiprocessing as mp;from concurrent.futures import ProcessPoolExecutor
    dbs=sys.argv[1].split(',');tag=sys.argv[2]
    recs=b.records(W+'/data',dbs)
    meth=[(v,m) for v in (sys.argv[3].split(',') if len(sys.argv)>3 else VAR) for m in ['titen','fusion']]+([('V0','minel')] if len(sys.argv)<=3 else [])
    jobs=[(db,r,fs,meth) for db,r,fs in recs]
    jobs.sort(key=lambda j:j[0]!='ltafdb')
    with ProcessPoolExecutor(9,mp_context=mp.get_context('spawn')) as ex:
        res=[x for o in ex.map(run,jobs,chunksize=1) for x in o]
    X=pd.DataFrame(res,columns=['db','record','variant','method','tp','fp','fn']);X.to_csv(W+f'/startup/res_{tag}.csv',index=False)
    print(len(X),'rows')
