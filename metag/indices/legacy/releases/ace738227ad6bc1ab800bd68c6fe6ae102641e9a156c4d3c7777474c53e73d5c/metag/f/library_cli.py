"""Build the optimized group index inside one shared, standalone library."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from common import CODE_ROOT, load_config
from index_selection import INDEX_NAMES, select_index
from library_paths import checked_catalog


def parser():
    p=argparse.ArgumentParser(description=__doc__,allow_abbrev=False)
    p.add_argument('action',nargs='?',default='build',choices=('build','group','check','status','kraken','human','sylph'))
    p.add_argument('--index-use',choices=INDEX_NAMES,default=None)
    p.add_argument('--config',default=str(CODE_ROOT/'config.json'))
    p.add_argument('--root')
    p.add_argument('--max-gib',type=float,default=None,help='Optimized hash cap; config default 384 GiB')
    p.add_argument('--estimate',choices=('auto','reuse','native','saved'),default=None)
    p.add_argument('--submit',choices=('hpc','local','auto'),default='auto')
    p.add_argument('--plan-only',action='store_true')
    p.add_argument('--full',action='store_true')
    p.add_argument('--queue',default='1t81c')
    p.add_argument('--cores',type=int,default=32)
    p.add_argument('--threads',type=int,default=32)
    p.add_argument('--mem-gib',type=int,default=640)
    p.add_argument('--hours',type=int,default=120)
    p.add_argument('--replace',choices=('TRUE','FALSE'),default='TRUE',help='Replace only an owned optimized index when its inputs/settings changed')
    return p


def main(argv=None):
    # Preserve established expert/reference-audit routes and all their arguments.
    argv = list(sys.argv[1:] if argv is None else argv)
    legacy = {'kraken': ['bash',str(CODE_ROOT/'f/entry.sh'),'library'],
              'human': ['bash',str(CODE_ROOT/'f/build_human.sh')],
              'sylph': ['bash',str(CODE_ROOT/'f/build_sylph.sh')]}
    if argv and argv[0] in legacy:
        if any(x == '--index-use' or x.startswith('--index-use=') for x in argv[1:]):
            raise SystemExit('Expert builders use --config, not --index-use.')
        # A bundle must never be converted in place into a full-library build.
        if argv[0] == 'kraken' and not any(x in ('-h','--help','--help-full') for x in argv):
            mini=argparse.ArgumentParser(add_help=False);mini.add_argument('--config');mini.add_argument('--root')
            options,_=mini.parse_known_args(argv[1:]);base=load_config(options.config,options.root)['library_dir']
            if (Path(base)/'library.json').exists() or Path(base).name=='lib0':
                raise SystemExit('Choose a separate new library_dir in --config for a full rebuild.')
        command=legacy[argv[0]]+argv[1:];os.execvp(command[0],command)
    if argv and argv[0]=='check' and '--references' in argv:
        command=['bash',str(CODE_ROOT/'f/analysis_entry.sh'),*argv];os.execvp(command[0],command)
    p=parser();a=p.parse_args(argv)
    cfg=select_index(load_config(a.config,a.root),a.index_use)
    name=cfg['index_use']
    if a.full and a.action!='check':p.error('--full is only valid with check')
    if a.action in ('build','group'):
        from partitioned_build import database_snapshot
        from group_library import source_catalog,compact_snapshot
        database_snapshot(cfg);source_catalog(cfg)
        if name=='legacy':
            compact_snapshot(cfg)
            print('Reusing the preserved legacy index: '+cfg['group_library_dir']);return
        if name=='full':
            print('Reusing the complete category indexes: '+cfg['db_dir']);return
        cap=a.max_gib if a.max_gib is not None else float(cfg.get('group_max_gib',384))
        estimate=a.estimate or cfg.get('group_estimate','auto')
        if min(a.cores,a.threads,a.mem_gib,a.hours)<1 or a.threads>a.cores or not 0<cap<a.mem_gib-8:
            p.error('Invalid resources: threads <= cores; hash cap plus reserve must fit memory')
        print('Library: '+cfg['library_dir']+'\nBuild: '+cfg['group_library_dir'],flush=True)
        cmd=[sys.executable,str(CODE_ROOT/'f/qc_submit.py'),'group_library',
             '--config',a.config,'--index-use',name,'--max-gib',str(cap),
             '--estimate',estimate,'--threads',str(a.threads),'--replace',a.replace,
             '--submit',a.submit,'--queue',a.queue,'--cores',str(a.cores),
             '--mem-gib',str(a.mem_gib),'--hours',str(a.hours),'--job-name','metag.index.'+name]
        if a.root:cmd+=['--root',a.root]
        if a.plan_only:cmd.append('--plan-only')
        subprocess.run(cmd,check=True);return
    if a.action=='check' and a.full:
        cmd=[sys.executable,str(CODE_ROOT/'f/qc_submit.py'),'library_check',
             '--config',a.config,'--index-use',name,'--full','--submit',a.submit,
             '--queue',a.queue,'--cores','2','--mem-gib','16','--hours',str(a.hours)]
        if a.root:cmd+=['--root',a.root]
        if a.plan_only:cmd.append('--plan-only')
        subprocess.run(cmd,check=True);return
    if a.plan_only:
        print('Metadata-only '+a.action+': '+cfg['library_dir']);return
    if a.action=='status':
        base=Path(cfg['library_dir'])
        result=dict(library=str(base),default_index=cfg['index_use'],indexes={})
        for index in INDEX_NAMES:
            path=base/('index' if index=='full' else 'group/'+index)/'database.json'
            if path.exists():
                data=checked_catalog(path);entry=dict(status=data['status'],identity=data['identity'],directory=str(path.parent))
                if index!='full':entry['selected_hash_fraction']=data['downsampling']['native_options']['selected_hash_fraction']
            else:entry=dict(status='not_built',directory=str(path.parent))
            result['indexes'][index]=entry
        print(json.dumps(result,indent=2));return
    from library_check import check
    check(cfg)


if __name__=='__main__':main()
