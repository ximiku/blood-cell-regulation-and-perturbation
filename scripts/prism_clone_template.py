"""Clone native Prism graph formatting while replacing only documented XML table data.

Usage: python prism_clone_template.py TEMPLATE_PZFX DATA_PZFX OUTPUT_DIR STEM
The compressed Template payload remains byte-for-byte unchanged. Graphs still
require Prism to load, render, save and export; generated projects are inputs.
"""
from pathlib import Path
import argparse,copy,hashlib,json
import xml.etree.ElementTree as E

def local(tag):return tag.split('}')[-1]
def unnamespace(e):
    e.tag=local(e.tag)
    for child in e:unnamespace(child)

def main():
    a=argparse.ArgumentParser();a.add_argument('template');a.add_argument('data');a.add_argument('output');a.add_argument('stem');args=a.parse_args()
    base=E.parse(args.template).getroot();data=E.parse(args.data).getroot();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    templates=[e for e in base if local(e.tag)=='Template'];assert len(templates)==1
    native_payload=templates[0].text;original_tables=[e for e in base if local(e.tag)=='Table'];assert len(original_tables)==1,'Use a native project with exactly one table.'
    table_id=original_tables[0].get('ID');manifest=[]
    for i,t in enumerate([e for e in data if local(e.tag)=='Table'],1):
        r=copy.deepcopy(base)
        for e in list(r):
            if local(e.tag) in ['Table','Info','InfoSequence']:r.remove(e)
        nt=copy.deepcopy(t);unnamespace(nt);nt.set('ID',table_id)
        r.insert(list(r).index(next(e for e in r if local(e.tag)=='Template')),nt)
        created=next((e for e in r if local(e.tag)=='Created'),None)
        if created is not None:
            for e in list(created):
                if local(e.tag)=='MostRecentVersion':created.remove(e)
            E.SubElement(created,'MostRecentVersion',CreatedByProgram='HW2 native-template data substitution',CreatedByVersion='1.0')
        assert next(e for e in r if local(e.tag)=='Template').text==native_payload
        p=out/f'{args.stem}_{i}.pzfx';E.indent(r);E.ElementTree(r).write(p,encoding='utf-8',xml_declaration=True)
        manifest.append({'path':str(p),'source_table':t.get('ID'),'native_template':args.template,'native_payload_sha256':hashlib.sha256(native_payload.encode()).hexdigest(),'rendered_in_Prism':False})
    (out/f'{args.stem}_preparation.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');print(json.dumps(manifest,ensure_ascii=False))
if __name__=='__main__':main()
