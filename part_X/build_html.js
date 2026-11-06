/* Build a single offline scholarly HTML, with prerendered mathematics and sources. */
const fs=require('fs'),path=require('path');
const here=__dirname,root=path.dirname(here),katex=require(process.env.KATEX_PATH||'katex');
const esc=s=>s.replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
let input=fs.readFileSync(path.join(here,'paper_source.html'),'utf8'),eq=0,inline=0;
input=input.replace(/\$\$([\s\S]*?)\$\$|\$([^$\n]+?)\$/g,(m,display,small)=>{
 const tex=(display||small).trim();let math=katex.renderToString(tex,{displayMode:!!display,throwOnError:true,strict:'warn',output:'htmlAndMathml',trust:false});
 if(display){eq++;return '<div class="equation" id="eq-'+eq+'"><div class="math">'+math+'</div><span class="eqno">(X.'+eq+')</span></div>'}
 inline++;return math;
});
input=input.replace(/<figure data-figure="([^"]+)">/g,(m,name)=>{
 let svg=fs.readFileSync(path.join(here,'figures',name+'.svg'),'utf8');svg=svg.slice(svg.indexOf('<svg'));
 svg=svg.replace('<svg ','<svg role="img" aria-label="'+name.replaceAll('-',' ')+'" ');
 return m+'<div class="plot-scroll" data-plot="'+name+'">'+svg+'</div>';
});
const programs=['replay.py','frame.py','frame_formalism.py','occupation.py','ground_submodularity.py','forests.py','dynamic.py','certification.py','nested_boxes.py','origin_cauchy.py','design.py','frequency_for_pairs.py','claude_uniform_rank.py','lossless_passivity.py','basis_policy.py','compensated_edge.py','requirements.txt'];
let files=[...programs,...fs.readdirSync(here).filter(x=>x.endsWith('_results.json')), 'replay_receipt.json','viii_stability_reference.json','union_7005_contour.json','same_input_comparison.json','ten_zero_fixture_provenance.json'];
let appendix=files.filter(f=>fs.existsSync(path.join(here,f))).map((name,i)=>{
 let text=fs.readFileSync(path.join(here,name),'utf8'),mime=name.endsWith('.json')?'application/json':'text/plain';
 return '<details class="source-details"><summary>'+esc(name)+' · '+(Buffer.byteLength(text)/1024).toFixed(1)+' KiB</summary><div class="source-actions"><button data-download="source-'+i+'" data-filename="'+esc(name)+'" data-mime="'+mime+'">Download '+esc(name)+'</button></div><pre><code id="source-'+i+'">'+esc(text)+'</code></pre></details>';
}).join('\n');
input=input.replace('<div id="source-entries"></div>','<div id="source-entries">'+appendix+'</div>');
const dist=path.dirname(require.resolve(process.env.KATEX_PATH||'katex'));
let mathCSS=fs.readFileSync(path.join(dist,'katex.min.css'),'utf8');
mathCSS=mathCSS.replace(/src:[^}]*?url\((fonts\/[^)]+\.woff2)\)[^}]*?(?=;font-weight)/g,(m,font)=>{
const data=fs.readFileSync(path.join(dist,font)).toString('base64');return 'src:url(data:font/woff2;base64,'+data+') format("woff2")';
});
/* Match distribution CSS robustly if source ordering changes. */
mathCSS=mathCSS.replace(/url\((fonts\/[^)]+)\)/g,(m,font)=>{
const ext=path.extname(font).slice(1),type=ext==='woff2'?'font/woff2':ext==='woff'?'font/woff':'font/ttf';return 'url(data:'+type+';base64,'+fs.readFileSync(path.join(dist,font)).toString('base64')+')';
});
const css=fs.readFileSync(path.join(here,'paper.css'),'utf8'),js=fs.readFileSync(path.join(here,'paper.js'),'utf8');
const html='<!doctype html>\n<html lang="en" data-theme="light"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="Part X of the Rank-Two Nonagon: elimination calculus, damped tomography, exact edge resolution, and composable contour certificates."><title>The Elimination Calculus · Nonagon X</title><style>'+mathCSS+'\n'+css+'</style></head><body><div id="reading-progress"></div>'+input+'<script>'+js+'</script></body></html>';
fs.writeFileSync(path.join(root,'nonagon_X.html'),html);
const prose=input.split('<section id="source-appendix">')[0].replace(/<annotation[\s\S]*?<\/annotation>/g,'').replace(/<math[\s\S]*?<\/math>/g,'').replace(/<[^>]*>/g,' ').split(/\s+/).filter(Boolean).length;
console.log(JSON.stringify({display_equations:eq,inline_equations:inline,embedded_sources:files.length,bytes:Buffer.byteLength(html),rendered_text_words_approx:prose,output:path.join(root,'nonagon_X.html')},null,2));
