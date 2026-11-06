(()=>{
const root=document.documentElement,theme=document.getElementById('theme');
theme.onclick=()=>{root.dataset.theme=root.dataset.theme==='dark'?'light':'dark';theme.textContent=root.dataset.theme==='dark'?'Light paper':'Dark paper';draw()};
let size=18;document.getElementById('smaller').onclick=()=>{size=Math.max(15,size-1);root.style.setProperty('--read',size+'px')};
document.getElementById('larger').onclick=()=>{size=Math.min(24,size+1);root.style.setProperty('--read',size+'px')};
document.getElementById('print').onclick=()=>window.print();
function progress(){const max=document.documentElement.scrollHeight-innerHeight;document.getElementById('reading-progress').style.width=(max>0?scrollY/max*100:0)+'%'}
addEventListener('scroll',progress,{passive:true});progress();
const NS='http://www.w3.org/2000/svg',svg=document.getElementById('edge-svg');
function elt(type,attrs,text){const e=document.createElementNS(NS,type);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;svg.appendChild(e);return e}
function draw(){
const w=+document.getElementById('edge-w').value,t=+document.getElementById('edge-t').value;
document.getElementById('edge-w-label').textContent=w.toFixed(2);document.getElementById('edge-t-label').textContent=t.toFixed(1);
const c=getComputedStyle(root),ox=c.getPropertyValue('--ox').trim(),teal=c.getPropertyValue('--teal').trim(),rule=c.getPropertyValue('--rule').trim();
svg.replaceChildren();const W=Math.max(280,svg.clientWidth),H=290,L=54,R=18,T=24,B=53;svg.setAttribute('viewBox','0 0 '+W+' '+H);const x=v=>L+v/20*(W-L-R),y=v=>H-B-(v-.30)/.22*(H-T-B);
for(const yy of [.30,.35,.40,.45,.50]){elt('line',{x1:L,x2:W-R,y1:y(yy),y2:y(yy),stroke:rule,'stroke-width':.6});elt('text',{x:L-10,y:y(yy)+4,'text-anchor':'end'},yy.toFixed(2))}
for(const xx of (W<440?[0,10,20]:[0,5,10,15,20]))elt('text',{x:x(xx),y:H-B+23,'text-anchor':'middle'},String(xx));
elt('text',{x:(W+L-R)/2,y:H-8,'text-anchor':'middle'},'Endpoint ground strength θ');
elt('text',{x:16,y:140,transform:'rotate(-90 16 140)','text-anchor':'middle'},'Port response h');
const F=z=>(4+4*w+z*(2+w))/(8+8*w+z*(5+3*w));
elt('line',{x1:L,x2:W-R,y1:y(.5),y2:y(.5),stroke:teal,'stroke-width':1.5,'stroke-dasharray':'5 4'});
const pts=Array.from({length:161},(_,i)=>{let q=i/8;return(i?'L':'M')+x(q)+','+y(F(q))}).join(' ');
elt('path',{d:pts,fill:'none',stroke:ox,'stroke-width':2.2});
elt('circle',{cx:x(t),cy:y(F(t)),r:5,fill:ox});elt('text',{x:W-R-4,y:y(.5)-7,'text-anchor':'end'},'Unloaded: h = 1/2');
let gamma=1/(5+3*w),der=-t*t/(8+8*w+t*(5+3*w))**2;
document.getElementById('edge-values').textContent='Loaded response '+F(t).toPrecision(10)+' · edge sensitivity ∂h/∂w = '+der.toExponential(5)+' · optimal contrast Γ = '+gamma.toPrecision(8)+' · clamp response '+((2+w)/(5+3*w)).toPrecision(8);
}
new ResizeObserver(draw).observe(svg);
document.getElementById('edge-w').addEventListener('input',draw);document.getElementById('edge-t').addEventListener('input',draw);draw();
document.querySelectorAll('[data-download]').forEach(btn=>btn.onclick=()=>{
const source=document.getElementById(btn.dataset.download),blob=new Blob([source.textContent],{type:btn.dataset.mime||'text/plain'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=btn.dataset.filename;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
});
})();
