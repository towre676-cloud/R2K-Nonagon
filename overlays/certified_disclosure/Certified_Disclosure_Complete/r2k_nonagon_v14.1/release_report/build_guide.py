#!/usr/bin/env python3
"""Build the integrated release guide from reviewed facts and retained receipts."""
from pathlib import Path
import json
import shutil
import html
import re
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak, KeepTogether
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / 'Nonagon_Integrated_Release.pdf'
INK = colors.HexColor('#17303c')
TEAL = colors.HexColor('#087d76')
MUTED = colors.HexColor('#526773')
PALE = colors.HexColor('#edf4f3')
LINE = colors.HexColor('#cbdadf')
S = {}


def setup():
    for name, file in [('Body','DejaVuSerif.ttf'),('BodyBold','DejaVuSerif-Bold.ttf'),
                       ('Sans','DejaVuSans.ttf'),('SansBold','DejaVuSans-Bold.ttf'),('Mono','DejaVuSansMono.ttf')]:
        pdfmetrics.registerFont(TTFont(name, '/usr/share/fonts/truetype/dejavu/' + file))
    pdfmetrics.registerFontFamily('Body', normal='Body', bold='BodyBold', italic='Body', boldItalic='BodyBold')
    pdfmetrics.registerFontFamily('Sans', normal='Sans', bold='SansBold', italic='Sans', boldItalic='SansBold')
    S.update({
        'kicker':ParagraphStyle('kicker',fontName='SansBold',fontSize=8.4,leading=12,textColor=TEAL,spaceAfter=11),
        'title':ParagraphStyle('title',fontName='BodyBold',fontSize=29,leading=33,textColor=INK,spaceAfter=15),
        'subtitle':ParagraphStyle('subtitle',fontName='Sans',fontSize=11.2,leading=16,textColor=MUTED,spaceAfter=16),
        'h1':ParagraphStyle('h1',fontName='BodyBold',fontSize=22,leading=27,textColor=INK,spaceAfter=15),
        'h2':ParagraphStyle('h2',fontName='SansBold',fontSize=10.2,leading=14,textColor=TEAL,spaceBefore=11,spaceAfter=6,keepWithNext=True),
        'body':ParagraphStyle('body',fontName='Body',fontSize=9.8,leading=14.7,textColor=INK,spaceAfter=9),
        'small':ParagraphStyle('small',fontName='Sans',fontSize=8.1,leading=11.5,textColor=MUTED,spaceAfter=7),
        'cell':ParagraphStyle('cell',fontName='Sans',fontSize=8.1,leading=11.4,textColor=INK),
        'head':ParagraphStyle('head',fontName='SansBold',fontSize=8.1,leading=11.4,textColor=colors.white),
        'code':ParagraphStyle('code',fontName='Mono',fontSize=8.0,leading=12,textColor=INK,backColor=PALE,borderPadding=9,spaceAfter=12),
    })


def rich(s):
    s=html.escape(s)
    s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s)
    for tag in ['sub','/sub','sup','/sup','br/']:
        s=s.replace('&lt;'+tag+'&gt;','<'+tag+'>')
    return s.replace('\n','<br/>')


def p(s, style='body'):
    return Paragraph(rich(s),S[style])


def eq(s, label):
    d=HERE/'figures';d.mkdir(exist_ok=True)
    fig=plt.figure(figsize=(8,0.62));fig.patch.set_alpha(0)
    fig.text(.5,.5,'$'+s+'$',ha='center',va='center',fontsize=14,color='#17303c')
    path=d/(label+'.png');fig.savefig(path,dpi=210,transparent=True,bbox_inches='tight',pad_inches=.08);plt.close(fig)
    w,h=ImageReader(str(path)).getSize();width=min(480,w/2.8)
    return Image(str(path),width=width,height=width*h/w)


def table(headers, rows, widths):
    data=[[p(x,'head') for x in headers]]+[[p(str(x),'cell') for x in row] for row in rows]
    t=Table(data,colWidths=widths,hAlign='LEFT')
    t.setStyle(TableStyle([
        ('BACKGROUND',(0,0),(-1,0),INK),('VALIGN',(0,0),(-1,-1),'TOP'),
        ('LEFTPADDING',(0,0),(-1,-1),9),('RIGHTPADDING',(0,0),(-1,-1),9),
        ('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,PALE]),
        ('LINEBELOW',(0,-1),(-1,-1),.5,LINE),
    ]))
    return t


def section(kicker,title):
    return [p(kicker.upper(),'kicker'),p(title,'h1')]


def main():
    global OUT
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=OUT)
    args=ap.parse_args()
    OUT=args.output.resolve()
    setup()
    facts=json.loads((HERE/'release_facts.json').read_text())
    if facts.get('status')!='reviewed':
        raise RuntimeError('Release facts must be reviewed before producing the final guide')
    story=[]
    story += [p('NONAGON / INTEGRATED RELEASE 9.0','kicker'),
              p('From a fixed port to a<br/>verified physical network','title'),
              p('The mathematical inverse, its resolution limits, and the calibrated acquisition program in one reproducible release.','subtitle')]
    story += [p('**Controlled grounding is the expansion that makes the Nonagon an experimental instrument.** A fixed source and readout remain at one port. Known positive shunts are routed to labeled sites. Their scalar responses determine a connected grounded network, predict unused settings, and support targeted fault diagnosis and component design.'),
              p('Part VIII now supplies the canonical theory. Part IX supplies the strict observation format, exact reconstruction, bounded-error propagation, engineering qualification, and acquisition extensions. The two independent inverse implementations agree on **all 103 retained networks**, using the same **2,430 acquisition responses and 309 held-out responses**. No synthetic replacement readings are used in that comparison.')]
    story += [table(['Layer','What is established'],[
        ['Finite theory','Exact n(n+1)/2 acquisition; qualified minimum count; Kron boundary; statistical and bounded-error resolution obstructions.'],
        ['Production acceptance','Strict rational input contract and certificate checks before the independent reference inverse runs.'],
        ['Engineering evidence','External SPICE, calibrated synthetic faults, measured-CSV import, and fresh finite-load and transient studies.'],
    ],[119,369]),Spacer(1,12)]
    story += [p('**The first physical platform is a passive resistor fixture.** Its nodal conductance matrix is exactly the required operator. Existing source-measure units, two programmable resistor channels and relay routing implement the operations. The release contains zero physical hardware measurements; the engineering results are circuit simulation and declared synthetic calibration experiments.'),
              p('Cody · 30 October 2026\nRelease guide. Detailed proofs, executable contracts, raw records and independent replay receipts remain in the corresponding chapter folders.','small')]

    story += [PageBreak()]+section('01 / Mathematical contract','What one fixed port can reveal')
    story += [p('Let K be a symmetric positive definite grounded network matrix, with nonpositive off-diagonal entries, nonnegative row sums and connected underlying graph. Its inverse G is entrywise positive. At fixed port 0, an added nonnegative diagonal load Θ gives the response below.'),
              eq(r'h(\Theta)=e_0^{\mathsf{T}}(K+\Theta)^{-1}e_0,\qquad G=K^{-1}','response'),
              p('Two single-site strengths determine G<sub>0j</sub><sup>2</sup> and G<sub>jj</sub>; positivity fixes the sign of G<sub>0j</sub>. A paired intervention gives a quadratic for G<sub>jk</sub>. The inverse-M-matrix path-product inequality selects its physical smaller root. Inverting the recovered G gives K. Neither simple modes, cyclicity nor known masses are needed for this static result.'),
              eq(r'1+2(n-1)+\binom{n-1}{2}=\frac{n(n+1)}{2}','count'),
              p('**The count is sharp under the stated design regularity.** Part VIII proves the lower bound on a full-dimensional open grounded cone for fixed protocols and protocols whose adaptive queries depend real-analytically on earlier readings. The release does not extend that proof to unrestricted discontinuous adaptive rules.'),
              p('**Access determines the identifiable object.** If only a subset A of sites can receive interventions, all responses depend on the Kron reduction K<sub>AA</sub> - K<sub>AH</sub>K<sub>HH</sub><sup>-1</sup>K<sub>HA</sub>. A hidden star and its equivalent triangle can therefore remain indistinguishable. Structural priors may change the inverse problem; unchanged accessible data alone cannot recover arbitrary hidden topology.'),
              p('**Resolution depends on the error model.** Independent Gaussian errors yield a two-model tail-probability lower bound at scale σ/(√m G<sub>0j</sub><sup>2</sup>). It is a worst-case statistical obstruction, not a pointwise standard-deviation floor. Arbitrary per-reading bounded errors give a separate worst-case grounding-error obstruction η/G<sub>0j</sub><sup>2</sup>, without a √m averaging gain. An unresolved sufficient interval calculation is not automatically an impossibility theorem.'),
              p('The number-variance repair is consolidated as well: the omitted interval [0, 10<sup>-4</sup>] contributes the leading correction 2La - a², equal to 0.00319999 at L = 16. Endpoint-corrected archived comparisons lie within 1.454 archived block standard errors. High-precision samples corroborate the local remainder; a uniform production-formula remainder certificate remains open.','small')]

    story += [PageBreak()]+section('02 / Finite acquisition design','Choose strengths for the readings actually taken')
    story += [p('The clamp-limit optimum and the implemented finite protocol solve different acquisition problems. A finite measurement paired with an ideal clamp has the first-order optimum θG<sub>jj</sub> = 1/√2 under its response-drop convention. The production protocol measures the raw baseline h<sub>0</sub>, h<sub>j</sub>(t), and h<sub>j</sub>(2t); both drops share the same baseline error.'),
              p('With q = G<sub>0j</sub><sup>2</sup>, b = G<sub>jj</sub> and u = tb, independently bounded raw errors of magnitude η give the first-order relative-q bound'),
              eq(r'\frac{|\delta q|}{q}\;\lesssim\;\eta\frac{b}{q}\left(\frac{4}{u}+8+4u\right)','bounded_design'),
              p('This coefficient is minimized at **u = 1**, with factor 16. For independent raw Gaussian noise, the corresponding squared coefficient is 6.5/u² + 24/u + 34 + 24u + 8u², minimized at u ≈ 0.9720604. These are local extraction results with known loads; whole-matrix inversion, resistor calibration and pilot cost still determine practical performance.'),
              p(facts['planner_summary']),
              table(facts['planner_table_headers'],facts['planner_table_rows'],[244,122,122]),Spacer(1,11),
              p(facts['planner_interpretation']),
              p('The planner keeps the existing observation contract: a single common t, two single strengths t and 2t, and pair strengths t. It does not relabel unequal site-specific settings into a fixed-strength certificate. Pilot readings and any reused configurations are counted explicitly. The independent checks exercise the acquisition Jacobian and uncertainty arithmetic. Exact edge absence uses the declared minimum-positive-edge gap 1/32.','small')]

    story += [PageBreak()]+section('03 / Dynamic acquisition','Recover the pencil from finite transients')
    story += [p('Part VIII applies the same inverse to P(s) = K + sC + s²M at positive real Laplace arguments. Grounded passive damping C and positive diagonal masses M preserve the M-matrix structure. Recovering P at three arguments separates the coefficients:'),
              eq(r'M=\frac{(P(s_2)-K)/s_2-(P(s_1)-K)/s_1}{s_2-s_1},\qquad C=\frac{P(s_1)-K}{s_1}-s_1M','dynamic'),
              p('The sufficient count 3n(n+1)/2 is not proved minimal: the parameter count n(n+2) leaves a gap of n(n-1)/2. Finite physical transients also require an explicit acquisition-error contract.'),
              p(facts['dynamic_summary']),
              table(facts['dynamic_table_headers'],facts['dynamic_table_rows'],[275,213]),Spacer(1,11),
              p(facts['dynamic_interpretation']),
              p('The earlier SPICE study checks four RC settling transients. This separate second-order demonstration carries finite-record displacement uncertainty into K, C and M, under its declared model and calibration premises.','small')]

    story += [PageBreak()]+section('04 / Engineering qualification','The useful first task is calibrated diagnosis')
    story += [p('**A known fault catalog is presently easier to qualify than an arbitrary unknown graph.** Six catalog-selected fixed-port readings distinguish all 100 fresh calibrated fault specimens at a 10 µV response envelope, including the declared resistor-setting and relay uncertainties. At 100 µV, 97 remain uniquely diagnosed with positive robust separation. A strong six-reading conventional reciprocal scan achieves 100/100 at both envelopes; the fixed-port method offers a different access arrangement, not a demonstrated universal accuracy or speed advantage.'),
              table(['Retained engineering check','Result'],[
                  ['External ngspice resistor fixtures / DC settings','14 / 286'],
                  ['Largest recovered normalized matrix error','2.46 × 10⁻¹²'],
                  ['Independent engineering replay','682,370 arithmetic predicates'],
                  ['Fresh calibrated catalog faults, six readings','100/100 at 10 µV; 97/100 at 100 µV'],
                  ['Full-graph combined-error qualification','72/72 unresolved'],
                  ['Physical hardware acquisitions','0'],
              ],[275,213]),Spacer(1,12),
              p('The 72 unresolved full-graph cases use 120 ppm of measured voltage plus 5 µV absolute error, 200 ppm current error, 0.1% shunt uncertainty and up to 0.05 Ω relay resistance. The same tested envelope still leaves a **3.326 mV gap** between the v7 twin predictions at 1 mA. Candidate separation, matrix enclosure and complete graph disclosure are therefore reported as distinct conclusions.'),
              p('**Existing equipment already matches the first interface.** A source-measure unit injects calibrated current and reads port voltage. Two isolated programmable resistors supply the maximum simultaneous shunt count. A relay matrix routes them to known pads. Examples documented in the engineering chapter are Keithley 2450, Pickering 40-261 and Keysight 34980A/34932A. These are capability mappings, not combined accuracy specifications for an assembled rack.'),
              p('Fault trials fix source gain; bench gain and drift must fit the response envelope. The importer normalizes V/I, encloses current and shunt calibration, records relay resistance and settling, and supports current reversal to cancel constant thermoelectric offsets. A six-node canonical plan has 21 logical scalar settings and 42 raw voltage samples with reversal, before repeats. Its completed examples are synthetic; the first measured coupon should use the same CSV contract.','small')]

    story += [PageBreak()]+section('05 / Reproduce and extend','One release, separate acceptance boundaries')
    story += [p('**Run from the release root.** The default replay checks preserved originals, production certificates, engineering arithmetic, the independent inverse bridge, corrected theory and both fresh studies. Numerical acquisition extensions need NumPy and SciPy; the rational production core uses the standard library.'),
              p('py -3.11 -m pip install numpy scipy\npy -3.11 replay_all.py','code'),
              p(facts['replay_summary']),
              table(['Path','Role'],[
                  ['part_VIII/','Canonical mathematical chapter and separately implemented reference inverse.'],
                  ['part_IX/','Production exact and bounded inverses, hardware plans and original engineering qualification.'],
                  ['part_IX/load_planning/','Finite common-load planner and fresh prospective comparison.'],
                  ['part_IX/dynamic_acquisition/','Finite transient-to-Laplace enclosures and joint coefficient study.'],
                  ['cross_checks/','Strict-gated bridge, independent controls and integrated replay receipts.'],
                  ['research/claude_v8/','Preserved source materials, with their original namespace and historical outputs.'],
              ],[206,282]),Spacer(1,12),
              p('The two verify_intervention.py files retain separate roles. **Part IX is the production acceptance gate** for serialized observations, physical class, calibration, extraction witnesses and independent forward equations. The Part VIII algebraic oracle runs after that gate; it independently replays every acquisition and held-out response.'),
              p('The first measured confirmation should use a coded resistor coupon, current-reversed calibration, a frozen fault catalog and a conventional-scan comparison. Report contacts, raw samples, relay settings, settling, elapsed time and the smallest reliable defect.'),
              p('Evidence: part_VIII/body8.html; part_IX/README.md; engineering_tests/ENGINEERING_REPORT.md; cross_checks/integration_receipt.json; acquisition-extension READMEs. Archived records retain their declared scope.','small')]

    def page(c,d):
        c.setStrokeColor(LINE);c.line(58,45,554,45)
        c.setFont('Sans',7.3);c.setFillColor(MUTED)
        c.drawString(58,31,'NONAGON 9.0  /  THEORY, ACQUISITION AND VERIFICATION')
        c.drawRightString(554,31,str(d.page))
    OUT.parent.mkdir(parents=True,exist_ok=True)
    doc=SimpleDocTemplate(str(OUT),pagesize=(612,792),leftMargin=62,rightMargin=62,topMargin=48,bottomMargin=60,
                          title='Nonagon 9.0: Integrated Theory and Engineering',author='Cody')
    doc.build(story,onFirstPage=page,onLaterPages=page)
    if OUT != HERE / OUT.name:
        shutil.copy2(OUT,HERE/OUT.name)
    print(json.dumps({'pdf':str(OUT),'bytes':OUT.stat().st_size}))


if __name__=='__main__':
    main()
