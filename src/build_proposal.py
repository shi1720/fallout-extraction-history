"""Build the AIP proposal using the official template's fields and character caps."""
import json
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'data/proposal.json').read_text())
assert len(d['abstract'])<=1300 and len(d['biography'])<=1300
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Top',fontName='Helvetica-Bold',fontSize=13,leading=17,alignment=TA_CENTER,spaceAfter=20))
styles.add(ParagraphStyle(name='Field',fontName='Helvetica',fontSize=11,leading=15,spaceAfter=10))
styles.add(ParagraphStyle(name='Label',fontName='Helvetica-Bold',fontSize=11,leading=15,spaceBefore=10,spaceAfter=5))
styles.add(ParagraphStyle(name='TextBody',fontName='Helvetica',fontSize=11,leading=15,spaceAfter=12))
out=ROOT/'paper/AIP_ECC_2027_Proposal_Shivam_Gupta.pdf'
doc=SimpleDocTemplate(str(out),pagesize=A4,leftMargin=55,rightMargin=55,topMargin=48,bottomMargin=48,title=d['title'],author=d['name'])
story=[Paragraph('7th AIP Early-Career Conference for<br/>Historians of the Physical Sciences',styles['Top'])]
for label,key in [('Name','name'),('E-mail address','email'),('Institutional affiliation','affiliation'),('Presentation title','title')]:
 story.append(Paragraph('<b>'+label+':</b> '+escape(d[key]),styles['Field']))
for label,key in [('Abstract','abstract'),('Short biography','biography')]:
 story += [Paragraph(label,styles['Label']),Paragraph(escape(d[key]),styles['TextBody'])]
doc.build(story)
print(out)
print('Character counts, including spaces:',len(d['abstract']),len(d['biography']))
