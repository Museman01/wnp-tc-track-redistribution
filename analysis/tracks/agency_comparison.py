'Cross-agency track, landfall and sensitivity analyses.'
from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import sys,json

from itertools import combinations
from pathlib import Path

import numpy as np,pandas as pd

from core.paths import PROJECT,output_path



import analysis.shared.record_composition as RC

from analysis.shared.inference_stats import trend_row, add_family_fdr, bh_fdr  # noqa: E402

from core.utils import load_config as project_config                     # noqa: E402

from analysis.shared.agency_catalogues import all_catalogues                             # noqa: E402



Z=_LAYOUT_ROOT / 'data/reference'

ZA=_LAYOUT_ROOT / 'data/reference/agency_comparison';ZV=_LAYOUT_ROOT / 'data/reference/environment_validation/supplemental_validation/results'

OUT=_LAYOUT_ROOT / 'results/agency_comparison'

AG4=['USA','TOKYO','CMA','HKO'];YEARS=np.arange(1966,2026)

def record(stage,info):
    print(stage,info,flush=True)



def with_hko_tables(cats):

    h=cats['HKO'];ts=h['ts_points']

    w=ts.loc[ts.wind.notna()];lmi=w.loc[w.groupby('sid')['wind'].idxmax(),['sid','season','iso_time','lat','lon','wind']].copy()

    lmi.columns=['sid','season','lmi_time','lmi_lat','lmi_lon','lmi_wind'];lmi['agency']='HKO'

    fr=ts.groupby('season')['sid'].nunique().rename('n_tc').reset_index();fr['agency']='HKO'

    h['lmi']=lmi;h['frequency']=fr;return cats



def common_four(cats):

    return set.intersection(*(set(cats[a]['ts_points'].loc[cats[a]['ts_points'].season.between(1966,2025),'sid']) for a in AG4))



def stage_core():

    core_out=OUT/'core'
    index_out=ZV if OUT.resolve()==(PROJECT/'results/agency_comparison').resolve() else OUT/'reference_index'
    targets=['trends_frequency_lmi_common.csv', 'lmi_common_four_agency_annual.csv', 'index_agreement_with_hko.csv']
    if (index_out/'index_hko_summary.csv').exists():
        raise FileExistsError(f'Refusing to overwrite existing index summary: {index_out}')
    if any((core_out/name).exists() for name in targets):
        raise FileExistsError(f'Refusing to overwrite existing core outputs: {core_out}')
    core_out.mkdir(parents=True,exist_ok=True)
    index_out.mkdir(parents=True,exist_ok=True)

    cfg=project_config();cats=with_hko_tables(all_catalogues());info={}

    common3=set.intersection(*(set(cats[a]['lmi']['sid']) for a in ['USA','TOKYO','CMA']))

    common4=common_four(cats);info['four_agency_common_storms']=len(common4)

    rows=[];lmi_rows=[]

    for fam_common,label in [(common3,'three'),(common4,'four')]:

        agencies=['USA','TOKYO','CMA'] if label=='three' else AG4

        tr=[]

        for a in agencies:

            f=cats[a]['frequency'].set_index('season')['n_tc'].reindex(YEARS,fill_value=0)

            l=cats[a]['lmi'];l=l.loc[l.season.between(1982,2025)]

            c=l.loc[l.sid.isin(fam_common)].groupby('season')['lmi_lat'].mean().reindex(range(1982,2026))

            tr.append(trend_row(f.index,f.values,label=a,family=f'frequency_{label}',cfg=cfg,extra={'analysis':'frequency','agency':a}))

            tr.append(trend_row(c.index,c.values,label=a,family=f'lmi_common_{label}',cfg=cfg,extra={'analysis':'lmi_common','agency':a}))

            if label=='four':

                for y,v in c.items():lmi_rows.append(dict(agency=a,season=y,mean_lmi_lat_common=v,n=int(l.loc[l.sid.isin(fam_common)&l.season.eq(y)].shape[0])))

        t=add_family_fdr(tr);t['common_set']=label;rows.append(t)

    t=pd.concat(rows,ignore_index=True)

    t.to_csv(core_out/'trends_frequency_lmi_common.csv',index=False);pd.DataFrame(lmi_rows).to_csv(core_out/'lmi_common_four_agency_annual.csv',index=False)

    ann=pd.read_csv(ZV/'wnp_tc_redistribution_index_annual.csv')

    h=cats['HKO']['ts_points'];h=h[h.season.between(1966,2025)&h.lon.between(100,180)&h.lat.between(0,40)]

    hf,*_=RC.annual_fields(h,YEARS,2.5);_,hoos,hpat=RC.projection_scores(hf,YEARS)

    hobs,hp=RC.block_permutation_projection(hf,YEARS,block=3,nperm=9999);htr=RC.sen_mk(YEARS,hoos)

    pd.DataFrame([dict(agency='HKO',weighting='track_point',late_minus_early_oos=hobs,block_permutation_p=hp,**htr)]).to_csv(index_out/'index_hko_summary.csv',index=False)

    pat=np.load(ZV/'wnp_tc_redistribution_pattern.npz')

    piv=ann[ann.weighting=='track_point'].pivot(index='year',columns='agency',values='index_oos');piv['HKO']=hoos

    agr=[]

    for a,b in combinations(['PRIMARY','USA','TOKYO','CMA','HKO'],2):

        pa=hpat if a=='HKO' else pat[f'{a}_track_point_pattern'];pb=hpat if b=='HKO' else pat[f'{b}_track_point_pattern']

        agr.append(dict(catalog_a=a,catalog_b=b,annual_index_correlation=float(piv[a].corr(piv[b])),change_pattern_correlation=float(np.corrcoef(pa,pb)[0,1])))

    agr=pd.DataFrame(agr);agr.to_csv(core_out/'index_agreement_with_hko.csv',index=False)

    info.update(hko_index=dict(late_minus_early=hobs,p=hp,sen=htr['slope_per_decade'],mk_p=htr['mk_p']),

        index_corr_range=[float(agr.annual_index_correlation.min()),float(agr.annual_index_correlation.max())],

        pattern_corr_range=[float(agr.change_pattern_correlation.min()),float(agr.change_pattern_correlation.max())])

    record('core',info)



# ---------------------------------------------------------------- landfall

import analysis.landfall.native_landfall_validation as LL                 # noqa: E402

import analysis.tracks.common_storm_validation as CS               # noqa: E402

import analysis.tracks.cutpoint_validation as CP              # noqa: E402

import analysis.landfall.coast_grouping_validation as CG        # noqa: E402

from analysis.shared.block_resampling import block_bootstrap_many, block_permutation_many, block_permutation_difference  # noqa: E402

CFG={'n_permutations':9999,'random_seed':202406,'n_bootstrap':4999}



def events4():

    f=pd.read_csv(_LAYOUT_ROOT / 'results/landfall/landfall_events_exact.csv');f=f.loc[f.pre_crossing_native_ts.astype(bool)]

    h=pd.read_csv(_LAYOUT_ROOT / 'results/agency_comparison/figure_series/hko_event_coast_classification.csv')

    h['coast']=np.where(h.nearest_distance_km.astype(float)<=50,h.nearest_segment,'Other');h['wind']=np.nan

    cols=['agency','sid','season','time','lat','lon','wind','coast','nearest_segment','nearest_distance_km']

    e=pd.concat([f[cols],h[cols]],ignore_index=True);e['time']=pd.to_datetime(e.time);e['coast']=e.coast.fillna('Other')

    return e.loc[e.season.between(1966,2025)].copy()



def stage_landfall():

    landfall_out=OUT/'landfall'
    if landfall_out.exists() and any(landfall_out.iterdir()):
        raise FileExistsError(f'Refusing to overwrite existing landfall outputs: {landfall_out}')
    landfall_out.mkdir(parents=True,exist_ok=True)

    info={};ev=events4()

    rows=[]

    for ag,g in ev.groupby('agency'):

        for dfn,fr in LL.definitions(g).items():

            for r in LL.annual_metrics(fr,YEARS):rows.append({'agency':ag,'definition':dfn,**r})

    ann=pd.DataFrame(rows);ann.to_csv(landfall_out/'landfall_annual_four_agency.csv',index=False)

    summ=LL.summarize(ann,CFG);summ.to_csv(landfall_out/'landfall_summary_four_agency.csv',index=False)

    # Leave-one-coast-out (first landfall, all events)

    lo=LL.leave_one_coast_out(ev,CFG);lo.to_csv(landfall_out/'leave_one_coast_out_four_agency.csv',index=False)

    # 15-yr segments of annual mean first-landfall latitude

    seg=[]

    for ag,g in ann[ann.definition.eq('first_landfall')].groupby('agency'):

        s=g.set_index('year').mean_lat

        v=[s.loc[a:a+14].mean() for a in (1966,1981,1996,2011)]

        seg.append(dict(agency=ag,s1966=v[0],s1981=v[1],s1996=v[2],s2011=v[3],early_internal=v[1]-v[0],late_retreat=v[3]-v[2]))

    pd.DataFrame(seg).to_csv(landfall_out/'landfall_segments_four_agency.csv',index=False)

    info['segments']=seg

    record('landfall',info)



def stage_common():

    result_out=OUT/'common'
    target=result_out/'common_storm_landfall_four_agency.csv'
    if target.exists():raise FileExistsError(f'Refusing to overwrite existing output: {target}')
    result_out.mkdir(parents=True,exist_ok=True)

    info={};ev=events4();cats=all_catalogues()

    early=np.flatnonzero(YEARS<=1995);late=np.flatnonzero(YEARS>=1996)

    def run(common,agencies):

        f=ev.loc[ev.sid.isin(common)&ev.agency.isin(agencies)].copy()

        first=f.sort_values('time').drop_duplicates(['agency','sid'],keep='first');rows=[]

        for ag in agencies:

            for dfn,fr in [('all_events',f),('first_landfall',first)]:

                part=fr.loc[fr.agency.eq(ag)]

                lat=part.groupby('season')['lat'].mean().reindex(YEARS)

                named=part.loc[part.coast.isin(CS.NORTH|CS.SOUTH)]

                n=named.coast.isin(CS.NORTH).groupby(named.season).sum().reindex(YEARS,fill_value=0)

                s_=(~named.coast.isin(CS.NORTH)).groupby(named.season).sum().reindex(YEARS,fill_value=0)

                share=(n/(n+s_)).where((n+s_)>0)

                matrix=np.column_stack([lat.values,share.values]).astype(float)

                ch,p=block_permutation_many(matrix,early,late,block=3,nperm=CFG['n_permutations'],seed=CFG['random_seed'])

                lo,hi=block_bootstrap_many(matrix[early],matrix[late],block=3,nboot=CFG['n_bootstrap'],seed=CFG['random_seed'])

                rows.append(dict(agency=ag,definition=dfn,late_minus_early=ch[0],ci_low=lo[0],ci_high=hi[0],block_permutation_p=p[0],n_storms=int(part.sid.nunique())))

        out=pd.DataFrame(rows);out['q_bh']=np.nan

        for d,idx in out.groupby('definition').groups.items():out.loc[idx,'q_bh']=bh_fdr(out.loc[idx,'block_permutation_p'])

        return out

    c4=common_four(cats);r4=run(c4,AG4);r4.to_csv(target,index=False)

    info['n_common4']=len(c4);record('common',info)



def stage_coast(thresholds, *, overwrite=False):

    sensitivity_out=OUT/'sensitivity'
    target=sensitivity_out/'coast_grouping_four_agency.csv'
    if target.exists() and not overwrite:
        raise FileExistsError(f'Refusing to overwrite sensitivity output: {target}')
    sensitivity_out.mkdir(parents=True,exist_ok=True)

    ev=events4();early,late=np.arange(30),np.arange(30,60);rows=[]

    taiwan={'taiwan_north':(CG.BASE_NORTH,CG.BASE_SOUTH),'taiwan_south':(CG.BASE_NORTH-{'Taiwan'},CG.BASE_SOUTH|{'Taiwan'}),'taiwan_excluded':(CG.BASE_NORTH-{'Taiwan'},CG.BASE_SOUTH)}

    for thr in thresholds:

        c=ev.copy();c['coast']=np.where(c.nearest_distance_km.astype(float)<=thr,c.nearest_segment,'Other')

        for rule in ['all_events','first_any','first_named']:

            fr=CG.select_events(c,rule)

            for ag in AG4:

                for opt,(no,so) in taiwan.items():

                    for den in ['named_only','include_other']:

                        sh=CG.annual_share(fr,ag,no,so,den)

                        eff,p=block_permutation_difference(sh,early,late,block=3,nperm=CFG['n_permutations'],seed=CFG['random_seed'])

                        rows.append(dict(agency=ag,event_rule=rule,taiwan_rule=opt,denominator=den,threshold_km=float(thr),change_percentage_points=eff*100,block_permutation_p=p))

    out=pd.DataFrame(rows)

    path=target

    old=pd.read_csv(path) if path.exists() else pd.DataFrame()

    if len(old):old=old[~old.threshold_km.isin([float(t) for t in thresholds])]

    out=pd.concat([old,out],ignore_index=True)

    out['fdr_family']=out.event_rule+'_'+out.taiwan_rule+'_'+out.denominator+'_'+out.threshold_km.map(lambda v:f'{v:g}km')

    out['q_bh']=np.nan

    for _,idx in out.groupby('fdr_family').groups.items():out.loc[idx,'q_bh']=bh_fdr(out.loc[idx,'block_permutation_p'])

    out.to_csv(path,index=False)

    record('coast_'+'_'.join(map(str,thresholds)),{'rows':int(len(out))})



def stage_cutpoint(schemes=('full_record','adjacent_25yr'),cuts=range(1991,2002), *, overwrite=False):

    sensitivity_out=OUT/'sensitivity'
    target=sensitivity_out/'cutpoint_four_agency.csv'
    if target.exists() and not overwrite:
        raise FileExistsError(f'Refusing to overwrite sensitivity output: {target}')
    sensitivity_out.mkdir(parents=True,exist_ok=True)

    ev=events4();cats=with_hko_tables(all_catalogues())

    from analysis.tracks.agency_composition_validation import annual_path_composition

    xe=np.arange(100,182.5,2.5);ye=np.arange(0,42.5,2.5)

    comp={a:annual_path_composition(cats[a]['ts_points'],YEARS,xe,ye) for a in AG4}

    land=CP.landfall_annual(ev.copy(),YEARS);core=pd.read_csv(ZA/'03_common_storms/core_crossagency_annual.csv')

    rows=[]

    for scheme in schemes:

        for cut in cuts:

            ey=np.arange(1966,cut) if scheme=='full_record' else np.arange(cut-25,cut);ly=np.arange(cut,2026) if scheme=='full_record' else np.arange(cut,cut+25)

            idx=np.array([np.flatnonzero(YEARS==y)[0] for y in np.r_[ey,ly]])

            for ag in AG4:

                tv,p=CP.field_test(comp[ag][idx],len(ey),block=3,nperm=CFG['n_permutations'],seed=CFG['random_seed'])

                rows.append(dict(scheme=scheme,cutpoint=cut,agency=ag,metric='path_total_variation',effect=tv,p_value=p))

                lf=land.loc[land.agency.eq(ag)]

                for col,met in [('first_landfall_mean_latitude','first_landfall_latitude_difference'),('north_share_all_events','north_share_difference')]:

                    e_,p_=CP.scalar_test(lf.season,lf[col],ey,ly,CFG);rows.append(dict(scheme=scheme,cutpoint=cut,agency=ag,metric=met,effect=e_,p_value=p_))

                if ag=='HKO':l=cats[ag]['lmi'];l=l.loc[l.season.between(1982,2025)].groupby('season')['lmi_lat'].mean()

                else:l=core.loc[core.agency.eq(ag)].set_index('season').mean_lmi_lat_full.dropna()

                e_,p_=CP.scalar_test(l.index,l.values,ey[ey>=1982],ly[ly>=1982],CFG);rows.append(dict(scheme=scheme,cutpoint=cut,agency=ag,metric='lmi_latitude_difference',effect=e_,p_value=p_))

    out=pd.DataFrame(rows);path=target

    if path.exists():

        old=pd.read_csv(path);old=old[~(old.scheme.isin(schemes)&old.cutpoint.isin(list(cuts)))];out=pd.concat([old,out],ignore_index=True)

    out.to_csv(path,index=False)

    record('cutpoint_'+'_'.join(schemes)+f'_{min(cuts)}_{max(cuts)}',{'rows':int(len(out))})



# ---------------------------------------------------------------- S06 / S07

def hko_composition(cats):

    from analysis.tracks.agency_composition_validation import annual_path_composition

    return annual_path_composition(cats['HKO']['ts_points'],YEARS,np.arange(100,182.5,2.5),np.arange(0,42.5,2.5))



def stage_pathdef(definitions, *, overwrite=False):

    pathdef_out=OUT/'path_definition'
    if definitions!=['check']:
        if not overwrite and pathdef_out.exists() and any(pathdef_out.iterdir()):
            raise FileExistsError(f'Refusing to overwrite existing path-definition outputs: {pathdef_out}')
        pathdef_out.mkdir(parents=True,exist_ok=True)

    import analysis.tracks.path_definition_validation as PD

    cats=all_catalogues();ev=events4()

    first=ev.sort_values('time').drop_duplicates(['agency','sid'])[['agency','sid','time']].rename(columns={'time':'first_landfall_time'})

    rows=[];changes={}

    defs=definitions

    for ag in ['HKO']:

        ts=cats[ag]['ts_points'].loc[cats[ag]['ts_points'].season.between(1966,2025)].copy();ts['agency']=ag

        full=cats[ag]['tracks'].loc[cats[ag]['tracks'].season.between(1966,2025)].copy();full['agency']=ag

        for d in defs:

            if d=='line_length':fields=PD.line_length_fields(ts,2.5)

            elif d=='pre_landfall_track_point':fields=PD.point_fields(PD.truncate_before_first_landfall(ts,first),2.5,'track_point')

            elif d=='full_life_track_point':fields=PD.point_fields(full,2.5,'track_point')

            else:fields=PD.point_fields(ts,2.5,d)

            PD.test_fields(rows,changes,ag,d,2.5,fields,[3],CFG,'native')

    out=pd.DataFrame(rows)

    path=pathdef_out/'path_definition_hko.csv'

    if path.exists():

        old=pd.read_csv(path);out=pd.concat([old[~old.definition.isin(defs)],out],ignore_index=True)

    out.to_csv(path,index=False)

    for (ag,d,w),ch in changes.items():np.save(pathdef_out/f'pathdef_change_{ag}_{d}.npy',ch)

    record('pathdef_'+'_'.join(defs),out.round(4).to_dict('records'))



def stage_climate():

    result_out=OUT/'climate'
    target=result_out/'climate_mode_four_agency.csv'
    if target.exists():raise FileExistsError(f'Refusing to overwrite existing output: {target}')
    result_out.mkdir(parents=True,exist_ok=True)

    import analysis.tracks.climate_adjustment_validation as CM

    cats=all_catalogues()

    z=np.load(ZA/'03_common_storms/annual_path_composition_2p5deg.npz');years=z['years'].astype(int)

    assert np.array_equal(years,YEARS)

    cl=pd.DataFrame({'season':years}).merge(pd.read_csv(PROJECT/'data/raw/indices/oni.csv'),on='season',how='left').merge(pd.read_csv(PROJECT/'data/raw/indices/pdo.csv'),on='season',how='left')

    cov=cl[['jas_oni','pdo']].to_numpy(float);cov=(cov-cov.mean(0))/cov.std(0,ddof=1)

    red=np.column_stack([np.ones(len(years)),cov]);per=(cl.season>=1996).astype(float).to_numpy()[:,None]

    comp={a:np.asarray(z[a],float) for a in ['USA','TOKYO','CMA']};comp['HKO']=hko_composition(cats)

    rows=[]

    for a in AG4:

        f,r2,p=CM.freedman_lane(np.sqrt(comp[a]),red,per,block=3,nperm=CFG['n_permutations'],seed=CFG['random_seed'])

        rows.append(dict(agency=a,pseudo_f=f,partial_r2=r2,block_permutation_p=p))

    out=pd.DataFrame(rows);out.to_csv(target,index=False)

    record('climate',out.round(4).to_dict('records'))



# ---------------------------------------------------------------- S09 / Table S3

def stage_speed():

    speed_out=OUT/'speed'
    if speed_out.exists() and any(speed_out.glob('speed_decomposition_*_four_agency.csv')):
        raise FileExistsError(f'Refusing to overwrite existing speed outputs: {speed_out}')
    speed_out.mkdir(parents=True,exist_ok=True)

    import analysis.tracks.speed_decomposition as SD

    from analysis.shared.inference_stats import add_family_fdr as fdr2

    from analysis.shared.agency_data import read_ibtracs_agencies, build_agency_catalog

    cfg=project_config();cfg['paths']={**cfg['paths'],'processed':str(PROJECT/'data/processed'),'raw':str(PROJECT/'data/raw')}

    raw=read_ibtracs_agencies(PROJECT/'data/raw/ibtracs_wp.csv',start=1980,end=2025)

    cats=all_catalogues()

    yrs=np.arange(1980,2024);annual_out=[];trend_rows=[];band_out=[]

    for ag in AG4:

        tr=build_agency_catalog(raw,ag)['tracks'] if ag!='HKO' else cats['HKO']['tracks'].loc[cats['HKO']['tracks'].season.between(1980,2025)]

        seg=SD._segments(tr,exact_six_hour=True)

        for st,col in [('total','speed_total'),('meridional','speed_meridional'),('zonal','speed_zonal')]:

            a,b=SD.decompose_segments(seg,yrs,lat_max=60,speed_column=col)

            a['agency']=ag;a['speed_type']=st;b['agency']=ag;b['speed_type']=st;annual_out.append(a);band_out.append(b)

            if st=='total':trend_rows.extend(SD._component_trends(a,cfg,f'{ag}_feng','total',1980,2023))

    pt=SD._project_tracks(cfg);a,b=SD.decompose_segments(SD._segments(pt,exact_six_hour=False),np.arange(1966,2026),lat_max=40,speed_column='speed_total')

    a['agency']='USA_project';a['speed_type']='total';annual_out.append(a)

    for st in (1966,1982):trend_rows.extend(SD._component_trends(a,cfg,f'PROJECT_ts_prelandfall_{st}','total',st,2025))

    tr=fdr2(trend_rows);tr.to_csv(speed_out/'speed_decomposition_trends_four_agency.csv',index=False)

    pd.concat(annual_out).to_csv(speed_out/'speed_decomposition_annual_four_agency.csv',index=False)

    pd.concat(band_out).to_csv(speed_out/'speed_decomposition_bands_four_agency.csv',index=False)

    qcol=[c for c in tr.columns if c.startswith('q')or c.endswith('q_bh') or c=='mk_q_bh']

    record('speed',{'qcols':qcol,'HKO':tr[tr.dataset.eq('HKO_feng')].round(5).to_dict('records')})



# ---------------------------------------------------------------- Table S7 / S10 sampling

ARC=PROJECT

def stage_sampling(agencies):

    sampling_out=OUT/'storm_sampling'
    target=sampling_out/'steering_sampling_four_agency.csv'
    if target.exists():raise FileExistsError(f'Refusing to overwrite existing sampling output: {target}')
    sampling_out.mkdir(parents=True,exist_ok=True)

    import importlib.util

    spec=importlib.util.spec_from_file_location('ss',_LAYOUT_ROOT / 'analysis/circulation/steering_sampling_decomposition.py');SS=importlib.util.module_from_spec(spec);spec.loader.exec_module(SS)

    cats=all_catalogues();rows=[]

    rng=np.random.default_rng(202406);blocks=np.arange(60).reshape(20,3)

    first=np.zeros((10000,60));second=np.zeros_like(first);first[0,:30]=1/30;second[0,30:]=1/30

    for k in range(1,10000):

        order=rng.permutation(20);first[k,blocks[order[:10]].ravel()]=1/30;second[k,blocks[order[10:]].ravel()]=1/30

    for layer in ['850-300','850-500']:

        U=dict(np.load(ARC/f'data/derived/sampling_{layer}.npz'))

        for ag in agencies:

            n=SS.distributions(cats[ag]['tracks']);N=n.reshape(60,-1)

            pool=n.mean(0);sp=pool/pool.sum((1,2),keepdims=True);nf=n.sum((2,3))[:,:,None,None]*sp[None]

            for c in ['u','v']:

                C=N@U[c].T;F=nf.reshape(60,-1)@U[c].T

                a=SS.values(C,first,second);f=SS.values(F,first,second)

                for name,v in {**a,'fixed_field':f['field'],'fixed_seasonal':f['sampling'],'fixed_cosampling':f['cosampling']}.items():

                    rows.append(dict(layer=layer,agency=ag,wind_component=c,term=name,change_ms=float(v[0]),two_sided_block_permutation_p=(int((abs(v[1:])>=abs(v[0])).sum())+1)/10000))

    out=pd.DataFrame(rows)

    ref=pd.read_csv(_LAYOUT_ROOT / 'results/circulation/steering_sampling_layers.csv')

    path=target

    allr=pd.concat([ref[out.columns],out[out.agency.eq('HKO')]],ignore_index=True).drop_duplicates(['layer','agency','wind_component','term'],keep='last')

    allr.to_csv(path,index=False)

    record('sampling_'+'_'.join(agencies),out[out.agency.eq('HKO')&out.term.isin(['total_change','field','sampling','cosampling','fixed_seasonal','fixed_field'])].round(4).to_dict('records') if 'HKO' in agencies else {})



SMS=_LAYOUT_ROOT / 'data/reference/storm_sampling'

def load_sts():

    import importlib.util

    spec=importlib.util.spec_from_file_location('sts',_LAYOUT_ROOT / 'analysis/circulation/storm_time_sampling.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m



def stage_s7perm(nperm=9999):

    """Block-permutation estimates for each agency."""
    target=OUT/'table_s7_permutation_four_agency.csv'
    if target.exists():raise FileExistsError(f'Refusing to overwrite existing output: {target}')
    OUT.mkdir(parents=True,exist_ok=True)
    sts=load_sts();cfg=sts.load_config();months=cfg['typhoon_season'];years=list(range(1966,2026))

    grid=dict(np.load(SMS/'steering_2p5deg_cache.npz'));U=sts.field_by_year_month(grid,months,years)

    assert all(np.isfinite(U[c]).all() for c in 'uv')

    cats=all_catalogues();rng=np.random.default_rng(cfg['statistics']['random_seed'])

    blocks=np.arange(60).reshape(20,3);rows=[]

    for ag in AG4:

        d,kept,_=sts.annual_point_distributions(cats[ag]['tracks'],cfg,months,years,grid)

        assert kept==years

        n=np.array([d[y] for y in years]);N=n.reshape(60,-1)

        pool=n.mean(0);sp=np.nan_to_num(pool/pool.sum((1,2),keepdims=True));nf=n.sum((2,3))[:,:,None,None]*sp[None]

        for comp in ('u','v'):

            W1=np.zeros((nperm+1,60));W2=np.zeros_like(W1);W1[0,:30]=1/30;W2[0,30:]=1/30

            for k in range(1,nperm+1):

                o=rng.permutation(20);W1[k,blocks[o[:10]].ravel()]=1/30;W2[k,blocks[o[10:]].ravel()]=1/30

            Uf=U[comp].reshape(60,-1);C=N@Uf.T;F=nf.reshape(60,-1)@Uf.T

            def vals(C):

                a=np.einsum('bi,ij,bj->b',W1,C,W1);b=np.einsum('bi,ij,bj->b',W1,C,W2);c=np.einsum('bi,ij,bj->b',W2,C,W1);dd=np.einsum('bi,ij,bj->b',W2,C,W2)

                delta=(W2-W1)@np.diag(C);return {'total_change':delta,'field':.5*(dd-c+b-a),'sampling':.5*(dd+c-b-a),'cosampling':delta-dd+a}

            a=vals(C);f=vals(F)

            terms={**a,'fixed_field':f['field'],'fixed_seasonal':f['sampling'],'fixed_cosampling':f['cosampling']}

            for name,v in terms.items():

                p=(int((abs(v[1:])>=abs(v[0])).sum())+1)/(nperm+1)

                rows.append(dict(agency=ag,wind_component=comp,term=name,change_ms=float(v[0]),p=p))

    out=pd.DataFrame(rows);out.to_csv(OUT/'table_s7_permutation_four_agency.csv',index=False)

    record('s7perm',{'HKO':out[out.agency.eq('HKO')].round(4).to_dict('records')})



def stage_s7sens():
    target=OUT/'table_s7_sensitivity_four_agency.csv'
    if target.exists():raise FileExistsError(f'Refusing to overwrite existing output: {target}')
    OUT.mkdir(parents=True,exist_ok=True)
    sts=load_sts();cfg=sts.load_config();months=cfg['typhoon_season']

    grid25=dict(np.load(SMS/'steering_2p5deg_cache.npz'));cats=all_catalogues();rows=[]

    import netCDF4

    with netCDF4.Dataset(PROJECT/'data/interim/steering.nc') as nc:

        lat=np.asarray(nc['latitude'][:],float);lon=np.asarray(nc['longitude'][:],float)

    lat=np.sort(lat[(lat>=0)&(lat<=40)]);lon=lon[(lon>=100)&(lon<=180)]

    e25=np.arange(0,42.5,2.5);l25=np.arange(100,182.5,2.5)

    wl=np.bincount(np.clip(np.digitize(lat,e25)-1,0,15),weights=np.cos(np.deg2rad(lat)),minlength=16)

    nl=np.bincount(np.clip(np.digitize(lon,l25)-1,0,31),minlength=32)

    den=wl[:,None]*nl[None,:]

    # 5-degree fields as exact weighted aggregates of the 2.5-degree cells

    def to5(a):

        num=(a*den).reshape(a.shape[0],8,2,16,2).sum((2,4));d=den.reshape(8,2,16,2).sum((1,3));return num/d

    grid5={'u':to5(grid25['u']),'v':to5(grid25['v']),'year':grid25['year'],'month':grid25['month'],'lon_edges':np.arange(100,185,5.0),'lat_edges':np.arange(0,45,5.0)}

    for case,grid,start in [('bin_5deg_1966',grid5,1966),('bin_2p5deg_1982',grid25,1982)]:

        years=list(range(start,2026));split=start+(2026-start)//2-1

        U=sts.field_by_year_month(grid,months,years);valid=np.isfinite(U['u']).all(0)&np.isfinite(U['v']).all(0)

        Uf={c:np.nan_to_num(U[c]) for c in 'uv'}

        for ag in AG4:

            tr=cats[ag]['tracks'];d,kept,_=sts.annual_point_distributions(tr,cfg,months,years,grid)

            a1=[y for y in years if y<=split and y in d];a2=[y for y in years if y>split and y in d]

            for comp in 'uv':

                r=sts.decompose(Uf[comp],d,years,a1,a2,valid)

                rows.append(dict(case=case,agency=ag,wind_component=comp,field=r['field'],fixed_seasonal=r['fixed_space']['seasonal'],fixed_cosampling=r['fixed_space']['cosampling'],total=r['total_change']))

    out=pd.DataFrame(rows);out.to_csv(OUT/'table_s7_sensitivity_four_agency.csv',index=False)

    record('s7sens',{'HKO':out[out.agency.eq('HKO')].round(4).to_dict('records')})



# ---------------------------------------------------------------- Table S6 OHC (HKO added)

def ohc_modules():
    from analysis.ocean import ohc_matching as matching
    from analysis.ocean import ohc_decomposition as decomposition
    from analysis.ocean import land_mask
    return matching, decomposition, land_mask


def hko_ohc_tracks():

    cats=all_catalogues();out=[]

    for d,key in [('ts_only','ts_points'),('eligible_full_track','tracks')]:

        x=cats['HKO'][key].loc[cats['HKO'][key].season.between(1966,2025),['sid','season','iso_time','lat','lon','wind']].copy()

        x.insert(0,'definition',d);x.insert(0,'agency','HKO');out.append(x)

    return pd.concat(out,ignore_index=True)



def stage_ohc_match(y0,y1):

    ocean_out=OUT/'ocean'
    target=ocean_out/f'ohc_hko_matched_{y0}_{y1}.csv.gz'
    if target.exists():raise FileExistsError(f'Refusing to overwrite existing OHC matches: {target}')
    ocean_out.mkdir(parents=True,exist_ok=True)

    M,D,LM=ohc_modules();cfg=json.loads((_LAYOUT_ROOT / 'config/ohc.json').read_text());dom=cfg['ohc']['domain']

    mpath=output_path('data/derived/ohc_land_mask.npz')

    if mpath.exists():mi=dict(np.load(mpath));mi={k:(v if k=='mask' else float(v)) for k,v in mi.items()}

    else:

        mi=LM.build_land_mask(PROJECT/cfg['matching']['gshhg_land_shapefile'],dom,float(cfg['matching']['land_mask_resolution_deg']));np.savez_compressed(mpath,**mi)

    prepared=_LAYOUT_ROOT / 'data/derived/ohc300'

    def match(points):

        points=points.copy();points['iso_time']=pd.to_datetime(points.iso_time)

        points['calendar_year']=points.iso_time.dt.year;points['match_status']='outside_domain'

        inside=points.lon.between(dom['west'],dom['east'])&points.lat.between(dom['south'],dom['north'])

        code=LM.classify_points(points.loc[inside,'lat'].to_numpy(float),points.loc[inside,'lon'].to_numpy(float),mi)

        points['point_is_land']=pd.NA;points.loc[inside,'point_is_land']=code.astype(bool)

        el=points.loc[inside].copy();el['ohc_year']=el.iso_time.dt.year;el['ohc_month']=el.iso_time.dt.month

        el=el.loc[el.ohc_year.between(1966,2025)]

        res=[M.match_group(g,prepared,'ohc300',int(y),int(m),6371.0088,75.0) for (y,m),g in el.groupby(['ohc_year','ohc_month'],sort=True)]

        return pd.concat(res) if res else el,points

    h=hko_ohc_tracks();h=h.loc[h.iso_time.dt.year.between(y0,y1)]

    mm,pts=match(h)

    cols=['ohc_year','ohc_month','grid_y','grid_x','grid_lat','grid_lon','match_km','ohc300','match_status']

    for c in cols:

        if c not in pts:pts[c]=np.nan

    pts.loc[mm.index,cols]=mm[cols].to_numpy()

    pts.to_csv(target,index=False)

    record(f'ohc_match_{y0}_{y1}',{'n':int(len(pts)),'matched':int(pts.match_status.eq('matched').sum())})



def stage_ohc_decomp(agency):

    ocean_out=OUT/'ocean'
    if agency=='HKO' and any((ocean_out/name).exists() for name in ('ohc_hko_summary.csv','ohc_hko_bootstrap.csv')):
        raise FileExistsError(f'Refusing to overwrite existing OHC results: {ocean_out}')
    ocean_out.mkdir(parents=True,exist_ok=True)

    M,D,LM=ohc_modules();cfg=json.loads((_LAYOUT_ROOT / 'config/ohc.json').read_text());dc=cfg['decomposition']

    if agency=='HKO':

        matched=pd.concat([pd.read_csv(f,parse_dates=['iso_time']) for f in sorted(ocean_out.glob('ohc_hko_matched_*.csv.gz'))],ignore_index=True)

    else:

        ref=pd.read_csv(_LAYOUT_ROOT / 'results/ocean/ohc300/matched_points.csv.gz',parse_dates=['iso_time']);matched=ref.loc[ref.agency.eq(agency)].copy()

    matched=matched.loc[matched.match_status.eq('matched')&matched.ohc300.notna()].copy()

    matched['season']=pd.to_numeric(matched.season).astype(int)

    for c in ['ohc_year','ohc_month','grid_y','grid_x']:matched[c]=pd.to_numeric(matched[c]).astype(int)

    matched['point_is_land']=matched.point_is_land.astype(str).str.strip().str.lower().map({'true':True,'false':False}).astype(bool)

    keys=D.state_table(matched);matched=D.attach_key_ids(matched,keys);matched['key_id']=matched.key_id.astype(int)

    years=list(range(1966,2026));yr={y:i for i,y in enumerate(years)}

    F=D.load_field_matrix(_LAYOUT_ROOT / 'data/derived/ohc300','ohc300',years,keys,float(cfg['ohc']['coordinate_equality_atol_degrees']))

    assert np.allclose(matched.ohc300.to_numpy(float),F[matched.ohc_year.map(yr).to_numpy(int),matched.key_id.to_numpy(int)],rtol=1e-7,atol=1e-3)

    srows=[];brows=[]

    for ya in dc['year_assignment_modes']:

        for lt in dc['land_treatments']:

            view=D.analysis_view(matched,ya,lt,1966,2025)

            for d in ['ts_only','eligible_full_track']:

                sub=view.loc[view.definition.eq(d)]

                for w in dc['weightings']:

                    W,_,_,_=D.annual_weight_matrix(sub,years,len(keys),w)

                    for comp in D.comparison_list(cfg):

                        r1=np.array([yr[y] for y in range(comp['period1'][0],comp['period1'][1]+1)]);r2=np.array([yr[y] for y in range(comp['period2'][0],comp['period2'][1]+1)])

                        summ,_=D.symmetric_decomposition(W[r1],F[r1],W[r2],F[r2],float(dc['closure_atol']),float(dc['closure_rtol']))

                        ident=dict(comparison=comp['label'],agency=agency,definition=d,year_assignment=ya,land_treatment=lt,weighting=w)

                        srows.append({**ident,**summ})

                        b=D.bootstrap_summary(W[r1],F[r1],W[r2],F[r2],int(dc['block_length_years']),int(dc['bootstrap_replicates']),int(dc['random_seed']),float(dc['closure_atol']),float(dc['closure_rtol']))

                        for k,v in ident.items():b[k]=v

                        brows.append(b)

    S=pd.DataFrame(srows);B=pd.concat(brows,ignore_index=True)

    if agency=='HKO':

        S.to_csv(ocean_out/'ohc_hko_summary.csv',index=False);B.to_csv(ocean_out/'ohc_hko_bootstrap.csv',index=False)

    main=S.loc[S.definition.eq('ts_only')&S.year_assignment.eq('calendar_year')&S.land_treatment.eq('ocean_points_only')&S.weighting.eq('storm_normalized_equal_year')]

    record(f'ohc_decomp_{agency}',main[['comparison','delta_exposure','redistribution_component','ocean_component','covariability_component']].round(0).to_dict('records'))



# ---------------------------------------------------------------- S01-S03 (records before 1966)

def stage_coverage():

    coverage_out=OUT/'coverage'
    if coverage_out.exists() and any(coverage_out.iterdir()):
        raise FileExistsError(f'Refusing to overwrite existing coverage outputs: {coverage_out}')
    coverage_out.mkdir(parents=True,exist_ok=True)

    from analysis.tracks.agency_composition_validation import annual_path_composition

    cats=all_catalogues(1945);rows=[]

    xe=np.arange(100,182.5,2.5);ye=np.arange(0,42.5,2.5)

    for ag in AG4:

        c=cats[ag];ts=c['ts_points'];tr=c['tracks']

        first=int(ts.season.min())

        for y in range(first,2026):

            t=tr.loc[tr.season.eq(y)];rows.append(dict(agency=ag,season=y,n_storms=int(ts.loc[ts.season.eq(y),'sid'].nunique()),

                n_ts_points=int(ts.season.eq(y).sum()),n_positions=int(len(t)),wind_available_pct=100*float(t.wind.notna().mean()) if len(t) else np.nan))

    cov=pd.DataFrame(rows);cov.to_csv(coverage_out/'coverage_four_agency.csv',index=False)

    proj=[];sy=[]

    for ag in AG4:

        ts=cats[ag]['ts_points'];yrs=np.arange(int(ts.season.min()),2026)

        F=annual_path_composition(ts,yrs,xe,ye)

        e=np.flatnonzero((yrs>=1966)&(yrs<=1995));l=np.flatnonzero(yrs>=1996)

        em,lm=F[e].mean(0),F[l].mean(0);pat=lm-em;mid=.5*(lm+em)

        sc=2*((F-mid)@pat)/float(pat@pat)

        for y,v in zip(yrs,sc):proj.append(dict(agency=ag,year=int(y),projection=v))

        for start in (1945,1966,1982):

            if start==1945 and ag=='HKO':continue  # 1945-start comparison kept to the three records available before 1961

            ee=np.flatnonzero((yrs>=start)&(yrs<=1995));d=lm-F[ee].mean(0)

            sy.append(dict(agency=ag,start=start,correlation_with_1966=float(np.corrcoef(d,pat)[0,1]),total_variation=float(.5*np.abs(d).sum()),n_early=len(ee)))

    pd.DataFrame(proj).to_csv(coverage_out/'projection_1945_2025_four_agency.csv',index=False)

    S=pd.DataFrame(sy);S.to_csv(coverage_out/'start_year_four_agency.csv',index=False)

    w=S.pivot(index='agency',columns='start',values='total_variation');info={'tv_drop_1945_pct':(100*(1-w[1945]/w[1966])).round(2).dropna().to_dict(),

        'corr_1945':S[S.start.eq(1945)].set_index('agency').correlation_with_1966.round(3).to_dict(),'first_year':cov.groupby('agency').season.min().to_dict(),

        'hko_wind_first':int(cov[cov.agency.eq('HKO')&cov.wind_available_pct.gt(0)].season.min())}

    record('coverage',info)



if __name__=='__main__':

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')
    if '--output-root' in sys.argv:
        index=sys.argv.index('--output-root')
        OUT=(PROJECT/Path(sys.argv[index+1])).resolve()
        if not OUT.is_relative_to(PROJECT.resolve()):raise ValueError('Output root must remain within the project')
        del sys.argv[index:index+2]
    st=sys.argv[1]

    if st=='coast':
        args=sys.argv[2:];stage_coast([int(v) for v in args if v!='--overwrite'], overwrite='--overwrite' in args)

    elif st=='pathdef':
        args=sys.argv[2:];stage_pathdef([arg for arg in args if arg!='--overwrite'], overwrite='--overwrite' in args)

    elif st=='climate':stage_climate()

    elif st=='speed':stage_speed()

    elif st=='ohc_match':stage_ohc_match(int(sys.argv[2]),int(sys.argv[3]))

    elif st=='ohc_decomp':stage_ohc_decomp(sys.argv[2])

    elif st=='coverage':stage_coverage()

    elif st=='s7perm':stage_s7perm(int(sys.argv[2]) if len(sys.argv)>2 else 9999)

    elif st=='s7sens':stage_s7sens()

    elif st=='sampling':stage_sampling(sys.argv[2:])

    elif st=='cutpoint':
        args=[v for v in sys.argv[2:] if v!='--overwrite'];overwrite='--overwrite' in sys.argv[2:]
        if args:
            if len(args)!=3:raise ValueError('cutpoint requires scheme, start year and end year')
            stage_cutpoint((args[0],),range(int(args[1]),int(args[2])+1),overwrite=overwrite)
        else:stage_cutpoint(overwrite=overwrite)

    else:{'core':stage_core,'landfall':stage_landfall,'common':stage_common,'cutpoint':stage_cutpoint}[st]()

