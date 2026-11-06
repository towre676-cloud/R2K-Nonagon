"""Render the engineering report and plots from retained numerical results."""
from pathlib import Path
import html,json,re
from fractions import Fraction
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,Image,Spacer,KeepTogether
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

HERE=Path(__file__).resolve().parent
OUT=HERE/'report';OUT.mkdir(exist_ok=True)
INK='#183747';TEAL='#087f79';GOLD='#b87c2e';GRAY='#61727b'


def figures():
    d=json.loads((HERE/'fault_calibration_results.json').read_text())
    keys=[('greedy_fixed_intervention',6,'Selected local\n6 readings'),
          ('canonical_intervention',21,'Canonical local\n21 readings'),
          ('greedy_fixed_reciprocal_scan',6,'Reciprocal scan\n6 readings')]
    fig,ax=plt.subplots(figsize=(8,3.1))
    x=np.arange(3)
    for shift,eta,label,col in [(-.18,'1/100000','10 microvolt response bound',TEAL),(.18,'1/10000','100 microvolt response bound',GOLD)]:
        y=[next(r['all_error_box_identifiable'] for r in d['summaries'] if r['eta']==eta and r['policy']==p and r['budget']==n) for p,n,_ in keys]
        bars=ax.bar(x+shift,y,.32,label=label,color=col)
        for b,v in zip(bars,y):ax.text(b.get_x()+b.get_width()/2,v+1,str(v),ha='center',fontsize=10,color=INK)
    ax.set_xticks(x,[k[2] for k in keys]);ax.set_ylim(0,118);ax.set_ylabel('Robustly distinguished faults / 100')
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    ax.legend(frameon=False,fontsize=9,loc='lower right');fig.tight_layout();fig.savefig(OUT/'fault_comparison.png',dpi=240,bbox_inches='tight');plt.close(fig)

    d=json.loads((HERE/'spice_results.json').read_text())
    case=next(x for x in d['transient']['cases'] if x['name']=='random_03_n6')
    data=np.loadtxt(HERE/case['raw_trace'],skiprows=1);t=data[:,0];voltage=data[:,1]
    k=np.array([[float(Fraction(x)) for x in r] for r in case['exact_K_normalized']])
    scale=float(Fraction(case['gscale_S']));current=float(Fraction(case['current_A']));cap=case['capacitance_each_F']
    eigen,basis=np.linalg.eigh(k*scale/cap);steady=np.linalg.solve(k*scale,np.eye(len(k))[0]*current)
    ref=steady[0]-(basis[0,:]*(basis.T@steady))@np.exp(-eigen[:,None]*t[None,:])
    fig,axes=plt.subplots(1,2,figsize=(8,2.7))
    axes[0].plot(t*1000,voltage*1000,color=TEAL,lw=2,label='ngspice')
    axes[0].plot(t[::110]*1000,ref[::110]*1000,'o',mfc='white',mec=GOLD,ms=3,label='Independent reference')
    axes[0].set_ylabel('Port voltage (mV)');axes[0].legend(frameon=False,fontsize=7,loc='lower right')
    axes[1].plot(t*1000,(voltage-ref)*1e9,color=TEAL,lw=1);axes[1].set_ylabel('Simulator - reference (nV)')
    for ax in axes:ax.set_xlabel('Time (ms)');ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
    fig.tight_layout();fig.savefig(OUT/'rc_transient.png',dpi=240,bbox_inches='tight');plt.close(fig)


def main():
    figures()
    for name,file in [('Body','DejaVuSerif.ttf'),('Bold','DejaVuSerif-Bold.ttf'),('Sans','DejaVuSans.ttf'),('Mono','DejaVuSansMono.ttf')]:
        pdfmetrics.registerFont(TTFont(name,'/usr/share/fonts/truetype/dejavu/'+file))
    pdfmetrics.registerFontFamily('Body',normal='Body',bold='Bold',italic='Body',boldItalic='Bold')
    base=dict(fontName='Body',textColor=colors.HexColor(INK))
    body=ParagraphStyle('Body',**base,fontSize=10,leading=15,spaceAfter=9)
    title=ParagraphStyle('Title',fontName='Bold',textColor=colors.HexColor(INK),fontSize=31,leading=36,spaceAfter=14)
    heading=ParagraphStyle('Heading',fontName='Bold',textColor=colors.HexColor(INK),fontSize=15,leading=19,spaceBefore=17,spaceAfter=10,keepWithNext=True)
    caption=ParagraphStyle('Caption',fontName='Sans',textColor=colors.HexColor(GRAY),fontSize=8,leading=11,spaceAfter=10)
    cell=ParagraphStyle('Cell',fontName='Sans',textColor=colors.HexColor(INK),fontSize=8.3,leading=11)
    header=ParagraphStyle('Header',parent=cell,textColor=colors.white)
    text=(HERE/'ENGINEERING_REPORT.md').read_text().replace('\u2014',' - ').replace('\u2013','-').replace('\u2011','-')
    def rich(s):
        s=html.escape(s);s=re.sub(r'\*\*(.*?)\*\*',r'<b>\1</b>',s)
        return re.sub(r'`(.*?)`',r'<font name="Mono">\1</font>',s)
    story=[Paragraph('NONAGON / ENGINEERING QUALIFICATION',caption)]
    blocks=re.split(r'\n\s*\n',text.strip())
    for b in blocks:
        if b.startswith('# '):
            story.append(Paragraph(rich(b[2:]),title));story.append(Paragraph('Cody · 30 October 2026 · Part VIII engineering extension',caption))
        elif b.startswith('## '):
            if b[3:]=='Calibration determines the useful inference task':
                story.extend([Image(str(OUT/'rc_transient.png'),width=480,height=162),Paragraph('Figure 1. External ngspice transient and an independent matrix-exponential calculation for the declared six-node RC fixture, with 1 microfarad at every site.',caption)])
            story.append(Paragraph(rich(b[3:]),heading))
        elif b.startswith('|'):
            lines=[line for line in b.splitlines() if not re.match(r'^\|[ :\-|]+\|$',line)]
            rows=[[x.strip() for x in line.strip('|').split('|')] for line in lines]
            widths={2:[315,165],4:[240,80,80,80]}.get(len(rows[0]),[215,45,110,110])
            cells=[[Paragraph(rich(v),header if i==0 else cell) for v in row] for i,row in enumerate(rows)]
            tab=Table(cells,colWidths=widths,repeatRows=1,hAlign='LEFT')
            tab.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor(TEAL)),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,1),(-1,-1),.4,colors.HexColor('#d9e3e5')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#f3f7f7'),colors.white])]))
            story.extend([KeepTogether([tab]),Spacer(1,11)])
            if len(rows[0])==4 and 'Scalar budget' in rows[0]:
                story.extend([Image(str(OUT/'fault_comparison.png'),width=480,height=186),Paragraph('Figure 2. Worst-case prediction-box discrimination on the 100-case calibrated confirmation catalog. Source current and conductance scaling are both 1 milli-unit, so the response bound is numerically in volts.',caption)])
        else:story.append(Paragraph(rich(b.replace('\n',' ')),body))
    def footer(c,d):
        c.setStrokeColor(colors.HexColor('#d9e3e5'));c.line(56,43,556,43);c.setFont('Sans',7.5);c.setFillColor(colors.HexColor(GRAY))
        c.drawString(56,30,'NONAGON · CIRCUIT VALIDATION AND FAULT DIAGNOSIS');c.drawRightString(556,30,str(d.page))
    path=OUT/'Nonagon_Engineering_Tests.pdf'
    SimpleDocTemplate(str(path),pagesize=(612,792),leftMargin=64,rightMargin=64,topMargin=49,bottomMargin=57,title='Nonagon Engineering Tests',author='Cody').build(story,onFirstPage=footer,onLaterPages=footer)
    print(json.dumps({'pdf':str(path),'words':len(text.split())}))


if __name__=='__main__':main()
