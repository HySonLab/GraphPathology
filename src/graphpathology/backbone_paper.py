"""Build the maintained matched-backbone manuscript without research execution."""
import json
import re
import shutil
import subprocess
from pathlib import Path

from .assets import sha
from .publication import compile_documents

PACKAGE='manuscript/journal_revision/transformer_egnn'


def build(root,output):
    source=root/PACKAGE;out=Path(output).resolve()
    if out.exists():raise FileExistsError(out)
    shutil.copytree(source,out,ignore=shutil.ignore_patterns('*.aux','*.log','*.out','*.fls','*.blg','*.pass*.txt','__pycache__'))
    reports=compile_documents(out,bibliography=True)
    layout=json.loads((out/'SOURCE_LAYOUT_MANIFEST.json').read_text())
    names={'main':'GraphPathology_Combined_Manuscript_Transformer_EGNN.pdf','article_only':'GraphPathology_Article_Only_Transformer_EGNN.pdf','supplementary':'GraphPathology_Supplementary_Transformer_EGNN.pdf'}
    for stem,name in names.items():shutil.copy2(out/(stem+'.pdf'),out/name)
    text=subprocess.check_output(['pdftotext','-layout',str(out/'main.pdf'),'-'],text=True)
    figures=set(re.findall(r'\bFigure\s+(S?\d+)\s*:',text));tables=set(re.findall(r'\bTable\s+(S?\d+)\s*:',text))
    expected_figures={str(i) for i in range(1,layout['main_figures']+1)}|{'S'+str(i) for i in range(1,layout['supplementary_figures']+1)}
    expected_tables={str(i) for i in range(1,layout['main_tables']+1)}|{'S'+str(i) for i in range(1,layout['supplementary_tables']+1)}
    assert figures==expected_figures,(figures,expected_figures)
    assert tables==expected_tables,(tables,expected_tables)
    refs=len(re.findall(r'\\bibitem',(out/'main.bbl').read_text()));assert refs==layout['references'],refs
    contents=(out/'article_content.tex').read_text()
    assert 'fig:extension_architecture' not in contents and 'Graph predictors and information boundaries' not in text
    assert contents.index('\\section{Introduction}')<contents.index('\\section{Related Work}')<contents.index('\\section{Results}')
    related=contents.split('\\section{Related Work}',1)[1].split('\\section{Results}',1)[0]
    assert related.count('\\noindent\\textbf')==9
    intro=contents.split('\\section{Introduction}',1)[1].split('\\section{Related Work}',1)[0]
    assert intro.count('\\item ')==4
    for i,page in enumerate(text.split('\f')):
        if page.strip():assert page.strip().splitlines()[-1].strip()==str(i+1)
    report={'status':'COMPLETE_DOCUMENT_BUILD','research_status':layout['research_status'],'documents':reports,'figures':len(figures),'tables':len(tables),'references':refs,'continuous_pages':True,'no_final_tex_diagnostics':True,'removed_schematic_absent':True,'intro_contributions':4,'related_work_themes':9,'outputs':{name:sha(out/name) for name in names.values()}}
    (out/'EXTENSION_BUILD_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def package(root,output,check_directory):
    from .editorial import package as create_package
    return create_package(root,output,check_directory,report_path=root/'experiments/journal_revision/backbone_graph_comparison/CLEAN_SOURCE_REPORT.json')


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);sub=p.add_subparsers(dest='action',required=True)
    b=sub.add_parser('build');b.add_argument('--output',type=Path,required=True)
    z=sub.add_parser('package');z.add_argument('--output',type=Path,required=True);z.add_argument('--check-directory',type=Path,required=True)
    a=p.parse_args();r=build(a.root.resolve(),a.output) if a.action=='build' else package(a.root.resolve(),a.output,a.check_directory);print(json.dumps(r,indent=2))
