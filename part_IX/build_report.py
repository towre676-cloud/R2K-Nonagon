"""Render the Part VIII chapter, equations and retained results to PDF and HTML."""
from pathlib import Path
import html
import json
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.utils import ImageReader
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, KeepTogether, PageBreak
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

HERE = Path(__file__).resolve().parent
OUT = HERE / 'report'
OUT.mkdir(exist_ok=True)
FIG = OUT / 'figures'
FIG.mkdir(exist_ok=True)
INK = '#152b3c'; TEAL = '#087e78'; GOLD = '#b97928'


def figures():
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.35), gridspec_kw={'width_ratios':[1,1.15]})
    ax = axes[0]; ax.axis('off')
    pos = np.array([[.12,.5],[.47,.8],[.83,.5]])
    for shift, label, edges in [(0,'Path',[(0,1),(1,2)]),(-.44,'Triangle',[(0,1),(1,2),(0,2)])]:
        p = pos.copy(); p[:,1] += shift
        for i,j in edges:ax.plot(p[[i,j],0],p[[i,j],1],color=TEAL,lw=2)
        ax.scatter(p[:,0],p[:,1],s=185,c='white',edgecolors=INK,zorder=3)
        for j in range(3):ax.text(*p[j],str(j+1),ha='center',va='center',color=INK,fontsize=9)
        ax.text(.04,p[:,1].max()+.075,label,color=INK,fontsize=11,fontweight='bold')
    ax.set_xlim(-.02,1);ax.set_ylim(-.08,1)
    ax=axes[1]
    baseline=5/18; path=7/26; twin=1088/3989
    ax.plot([0,1],[baseline,path],'-o',color=TEAL,lw=2,label='Path')
    ax.plot([0,1],[baseline,twin],'-s',color=GOLD,lw=2,label='Triangle')
    ax.annotate('5/18', (0,baseline), xytext=(12,3),textcoords='offset points',fontsize=10)
    ax.annotate('7/26', (1,path), xytext=(-48,-5),textcoords='offset points',fontsize=10,color=TEAL)
    ax.annotate('1088/3989', (1,twin), xytext=(-83,7),textcoords='offset points',fontsize=10,color=GOLD)
    ax.set_xticks([0,1],['Original','Unit load at node 2']);ax.set_xlim(-.12,1.14)
    ax.set_ylabel('Static compliance (normalized)');ax.legend(frameon=False,loc='lower left',fontsize=9)
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.18)
    fig.tight_layout();fig.savefig(FIG/'twins.png',dpi=240,bbox_inches='tight');plt.close(fig)

    data=json.loads((HERE/'offset_audit.json').read_text())['corrected_archived_BK_rows']
    fig,axs=plt.subplots(1,2,figsize=(9,3.25))
    cc=[TEAL,GOLD,'#5569a5']
    for row,c in zip(data,cc):
        L=[r['L'] for r in row['rows']]
        axs[0].plot(L,[r['old_BK_minus_observed'] for r in row['rows']],'-o',color=c,ms=3,label=f"Index {row['n0']:,}")
        axs[1].plot(L,[r['corrected_residual_over_archived_block_se'] for r in row['rows']],'-o',color=c,ms=3)
    x=np.linspace(.5,16,100);axs[0].plot(x,2e-4*x-1e-8,'--',color=INK,lw=1.5,label='Omitted endpoint')
    axs[0].set_ylabel('Archived prediction minus observation')
    axs[1].axhspan(-2,2,color=TEAL,alpha=.06);axs[1].axhline(0,color=INK,lw=.6)
    axs[1].set_ylabel('Corrected residual / archived block SE');axs[1].set_ylim(-2.1,2.1)
    for ax in axs:
        ax.set_xlabel('Window length L');ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
    axs[0].legend(frameon=False,fontsize=7);fig.tight_layout()
    fig.savefig(FIG/'offset.png',dpi=240,bbox_inches='tight');plt.close(fig)


def equation(tex,index):
    p=FIG/f'equation_{index}.png'
    fig=plt.figure(figsize=(9,.75))
    fig.text(.5,.5,'$'+tex+'$',ha='center',va='center',fontsize=16,color=INK)
    fig.savefig(p,dpi=260,bbox_inches='tight',pad_inches=.1,transparent=True)
    plt.close(fig);return p


def rich(s):
    s=html.escape(s)
    s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s)
    # In the PDF, descriptive source labels keep long URLs from spilling margins.
    if s[:2].strip().rstrip('.').isdigit() and 'https://' in s:
        label,url=s.split('https://',1)
        url='https://'+url.split(' ; ')[0]
        return label+f'<link href="{url}"><font color="{TEAL}">manufacturer source</font></link>'
    return s


def main():
    figures()
    for face,file in [('Body','DejaVuSerif.ttf'),('BodyBold','DejaVuSerif-Bold.ttf'),('Sans','DejaVuSans.ttf')]:
        pdfmetrics.registerFont(TTFont(face,'/usr/share/fonts/truetype/dejavu/'+file))
    pdfmetrics.registerFontFamily('Body',normal='Body',bold='BodyBold',italic='Body',boldItalic='BodyBold')
    styles={
      'body':ParagraphStyle('body',fontName='Body',fontSize=9.5,leading=14.3,textColor=colors.HexColor(INK),spaceAfter=8.5),
      'h1':ParagraphStyle('h1',fontName='BodyBold',fontSize=29,leading=34,textColor=colors.HexColor(INK),spaceAfter=13),
      'h2':ParagraphStyle('h2',fontName='BodyBold',fontSize=14,leading=19,textColor=colors.HexColor(INK),spaceBefore=13,spaceAfter=8,keepWithNext=True),
      'caption':ParagraphStyle('caption',fontName='Sans',fontSize=8,leading=11,textColor=colors.HexColor('#526571'),spaceAfter=12),
      'eyebrow':ParagraphStyle('eyebrow',fontName='Sans',fontSize=9,leading=14,textColor=colors.HexColor(TEAL),spaceAfter=12)}
    story=[Paragraph('NONAGON / PART VIII',styles['eyebrow'])]
    text=(HERE/'CHAPTER.md').read_text(); blocks=re.split(r'\n\s*\n',text.strip())
    eq=0
    for block in blocks:
        block=block.strip()
        if block.startswith('# '):
            story.append(Paragraph(rich(block[2:]),styles['h1']))
            story.append(Paragraph('Cody · 30 October 2026',styles['eyebrow']))
            story.append(Paragraph('Controlled interventions reveal a network, check its reconstructed connections, and select changes that meet a response target.',styles['body']))
        elif block.startswith('## '):
            title=block[3:]
            if title=='The number-variance offset: isolate the common numerical operation':story.append(PageBreak())
            story.append(Paragraph(rich(title),styles['h2']))
        elif block.startswith('$$'):
            eq+=1;p=equation(block.strip('$').strip(),eq);w,h=ImageReader(str(p)).getSize();width=min(460,w/2.6)
            image=Image(str(p),width=width,height=width*h/w)
            if eq==4:
                previous=story.pop();story.append(KeepTogether([previous,image,Spacer(1,7)]))
            else:
                story.append(image);story.append(Spacer(1,7))
        else:
            for para in block.split('\n'):
                story.append(Paragraph(rich(para),styles['body']))
            if block.startswith('All six canonical'):
                story.append(Image(str(FIG/'twins.png'),width=465,height=173))
                story.append(Paragraph('Figure 1. The same labeled intervention separates the v7 twins at one fixed readout.',styles['caption']))
            if block.startswith('The supplied audit separates'):
                story.append(Image(str(FIG/'offset.png'),width=465,height=169))
                story.append(Paragraph('Figure 2. The omitted endpoint follows the large-L discrepancy. Right: endpoint-corrected residuals use the archived block standard errors; full BK re-evaluation remains a separate run.',styles['caption']))
    def page(c,d):
        w,h=d.pagesize;c.setStrokeColor(colors.HexColor('#d4dfe4'));c.line(48,43,w-48,43)
        c.setFont('Sans',7.3);c.setFillColor(colors.HexColor('#526571'))
        c.drawString(48,30,'NONAGON · ACTIVE SPECTRAL TOMOGRAPHY');c.drawRightString(w-48,30,str(d.page))
    target=OUT/'Nonagon_Part_VIII_Active_Spectral_Tomography.pdf'
    doc=SimpleDocTemplate(str(target),pagesize=(612,792),leftMargin=62,rightMargin=62,topMargin=48,bottomMargin=58,title='Active Spectral Tomography and Inverse Design',author='Cody')
    doc.build(story,onFirstPage=page,onLaterPages=page)
    # Self-contained report page uses equation images, allowing offline reading.
    pageparts=[];eq=0
    for block in blocks:
        if block.startswith('# '):pageparts.append('<h1>'+html.escape(block[2:])+'</h1><p class="byline">Cody · Nonagon Part VIII · 30 October 2026</p>')
        elif block.startswith('## '):pageparts.append('<h2>'+html.escape(block[3:])+'</h2>')
        elif block.startswith('$$'):
            import base64
            eq+=1;src=base64.b64encode((FIG/f'equation_{eq}.png').read_bytes()).decode();pageparts.append('<div class="equation"><img alt="'+html.escape(block.strip('$'))+'" src="data:image/png;base64,'+src+'"></div>')
        else:pageparts.append('<p>'+re.sub(r'\*\*(.+?)\*\*',r'<strong>\1</strong>',html.escape(block))+'</p>')
    css='body{margin:0;background:#f5f3ed;color:#193443;font:17px/1.78 Georgia,serif}main{max-width:850px;margin:0 auto;padding:55px 26px}h1{font-size:clamp(34px,6vw,54px);line-height:1.12;letter-spacing:-.04em}h2{font-size:27px;line-height:1.3;margin-top:42px}p{margin:0 0 20px;overflow-wrap:anywhere}.byline{font:13px/1.7 system-ui;color:#087e78}.equation{margin:26px 0;background:white;padding:16px}.equation img{display:block;max-width:100%;height:auto;margin:auto}a{color:#087e78}nav{font:14px system-ui;margin-bottom:32px}@media print{body{background:white}main{padding:0;max-width:none}h2{break-after:avoid}}'
    content='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Nonagon Part VIII · Research Chapter</title><style>'+css+'</style><main><nav><a href="./">← Open the interactive laboratory</a></nav>'+''.join(pageparts)+'</main></html>'
    (OUT/'chapter.html').write_text(content)
    print(json.dumps({'pdf':str(target),'html':str(OUT/'chapter.html'),'words':len(text.split())}))


if __name__=='__main__':main()
