"""Track deliberate display changes without replacing historical target ledgers."""
import csv,json,re,subprocess
from pathlib import Path

PACKAGE='manuscript/journal_revision/transformer_egnn'
START='ab0811134bb1d62719206ab11a43bd7c4dd5d813'


def inventory(read):
    def expand(text):
        return re.sub(r'\\input\{([^}]+)\}',lambda m:expand(read(m.group(1) if m.group(1).endswith('.tex') else m.group(1)+'.tex')),text)
    result={}
    for source in ['article_content.tex','supplementary_content.tex']:
        text=expand(read(source))
        for match in re.finditer(r'\\begin\{(figure|table|longtable)\}(.*?)\\end\{\1\}',text,re.S):
            labels=re.findall(r'\\label\{([^}]+)\}',match.group(2))
            if not labels:
                title=re.search(r'\\caption\{\\textbf\{([^}]+)',match.group(2)).group(1)
                labels=[{'Primary gene panel.':'tab:primary_gene_panel','Patient-resolved reconstruction metrics.':'tab:patient_resolved_metrics'}[title]]
            assert len(labels)==1,(source,labels)
            label=labels[0];assert label not in result
            result[label]={'kind':'figure' if match.group(1)=='figure' else 'table','source':source,'supplementary':source.startswith('supplementary')}
    return result


def baseline(root):
    path=root/PACKAGE/'DISPLAY_BASELINE.json'
    if path.exists():return json.loads(path.read_text())
    def read(name):return subprocess.check_output(['git','show',f'{START}:{PACKAGE}/{name}'],cwd=root,text=True)
    displays=inventory(read);assert len(displays)==62
    result={'starting_commit':START,'displays':displays,'implicit_label_mapping':{'tab:primary_gene_panel':'Prior caption Primary gene panel.; explicit label added in current source without changing numbering.','tab:patient_resolved_metrics':'Prior caption Patient-resolved reconstruction metrics.; explicit label added in current source without changing numbering.'}};path.write_text(json.dumps(result,indent=2)+'\n');return result


def check(root,build):
    source=root/PACKAGE;old=baseline(root)['displays'];current=inventory(lambda name:(source/name).read_text())
    removed=set(old)-set(current);assert removed=={'fig:extension_architecture'},removed
    allowed={'tab:backbone_resources','tab:matched_backbone','fig:matched_backbone_effects'}
    added=set(current)-set(old);assert added<=allowed,added
    for label in set(old)&set(current):assert old[label]==current[label],label
    aux=(build/'main.aux').read_text();numbers={label:(number,page) for label,number,page in re.findall(r'\\newlabel\{([^}]+)\}\{\{([^}]+)\}\{([^}]+)\}',aux)}
    rows=[]
    for label,entry in current.items():
        assert label in numbers,label
        number,page=numbers[label];rows.append({'label':label,**entry,'number':number,'combined_page':page,'change':'ADDED' if label in added else 'RETAINED'})
    with (build/'CURRENT_DISPLAY_MAP.csv').open('w') as f:
        w=csv.DictWriter(f,lineterminator='\n',fieldnames=['label','kind','source','supplementary','number','combined_page','change']);w.writeheader();w.writerows(rows)
    return {'prior_displays':len(old),'current_displays':len(current),'removed_labels':sorted(removed),'added_labels':sorted(added),'empirical_displays_removed':0,'all_other_labels_and_roles_retained':True}

if __name__=='__main__':
    print(json.dumps(baseline(Path(__file__).resolve().parents[2]),indent=2))
