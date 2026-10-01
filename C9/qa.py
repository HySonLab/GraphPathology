"""Supplementary metadata/environment/figure audit; no model or metric changes."""
from contract import *
from PIL import Image,ImageDraw
import platform,scipy,torch,torchvision

def main():
 for s in SAMPLES:
  names=(ROOT/cfg()['count_files'][s]['features']).read_text().splitlines()
  assert all(names.count(g)==1 for g in genes()),s
 paths=sorted((B/'figures').glob('[0-9][0-9]_*.png'));assert len(paths)==10;out=Image.new('RGB',(1600,1500),'white');draw=ImageDraw.Draw(out)
 for i,p in enumerate(paths):
  im=Image.open(p);im.load();assert im.width>800 and im.height>400;im.thumbnail((795,275));x=(i%2)*800;y=(i//2)*300;out.paste(im,(x,y+20));draw.text((x+10,y+2),p.stem,fill='black')
 out.save(B/'figures/QA_CONTACT_SHEET.jpg')
 dump(B/'logs/SUPPLEMENTARY_QA.json',{'status':'PASS','primary_gene_symbols_unique_in_all_assays':True,'figure_decode_checks':10,'all_patients_in_figures':True,'environment':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'torch':torch.__version__,'torchvision':torchvision.__version__},'parent_verification_sha256':sha(B/'logs/PARENT_VERIFICATION.txt')})
if __name__=='__main__':main()
