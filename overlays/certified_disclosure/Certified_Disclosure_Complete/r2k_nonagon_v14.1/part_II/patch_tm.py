s=open('body2.html').read()
def rep(a,b):
    global s
    assert a in s, a[:60]; s=s.replace(a,b)
rep('<tr><td>N1</td><td>cost of rigorous moments at height</td><td>open, measured</td><td>2.2 s at $t=14$; over 600 s at zero 3000</td></tr>',
    '<tr><td>N1</td><td>cost of rigorous moments at height</td><td>closed for the wrapping barrier; point-evaluation cost open</td><td>Lehmer tile: over 100 min → 19 s</td></tr>')
rep('<tr><td>P2</td><td>discriminant meter at Lehmer pairs</td><td>closed</td><td>collision $t_c=-7.112\\times10^{-4}$</td></tr>',
    '<tr><td>P2</td><td>discriminant meter at Lehmer pairs</td><td>closed</td><td>rigorous gap $0.0376984977261402\\pm3.1\\times10^{-17}$; $t_c=-7.112\\times10^{-4}$</td></tr>')
rep('<strong>P2</strong>: the discriminant meter gives a rigorous enclosure of the gap of the classic Lehmer pair,',
    '<strong>P2</strong>: the discriminant meter gives a rigorous enclosure of the gap of the classic Lehmer pair, $0.0376984977261402\\pm3.1\\times10^{-17}$, without locating either zero,')
rep('The Davenport–Heilbronn off-line pair is proved by a rigorously negative $2\\times2$ minor.',
    'The Davenport–Heilbronn off-line pair is proved by a rigorously negative $2\\times2$ minor. A new <strong>Taylor-model contour integrator</strong> removes the wrapping effect that stalls ball-arithmetic quadrature at height: Arb\'s adaptive integrator ran for more than 100 minutes on a tile at $t=7005$ without finishing, and the Taylor-model integrator certifies the same tile in 19 seconds with a height-independent segment count.')
rep('<h3>The margin law</h3>', open('tm_section.html').read() + '\n<h3>The margin law</h3>')
rep('<!--LEHMER-->','<!--LEHMERTM-->')
a=s.index('<p><strong>N1, the cost of rigorous moments at height</strong>'); b=s.index('<p><strong>N2,')
s=s[:a]+open('n1.html').read()+'\n\n'+s[b:]
rep('<code>rig_p1.py</code> (Arb certificates; the 19 tiles in about 60 s, boxes in 207 s and 179 s, Davenport–Heilbronn in 63 s, the Lehmer tile at height 7005 the longest)',
    '<code>rig_p1.py</code> (Arb certificates; the 19 tiles in about 60 s, boxes in 207 s and 179 s, Davenport–Heilbronn in 63 s), <code>tmint.py</code> (the Taylor-model integrator of Theorem II.1.7), <code>lehmer_tm.py</code> (Lehmer certificate, 19 s), <code>tmcost.py</code> (height sweep)')
open('body2.html','w').write(s); print('patched')
