/* Copyright (C) 2026 Andrew D. Smith
 *
 * FalcoAmadeus fork modifications Copyright (C) 2026 Joachim Wolff
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *   
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *   
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef SRC_HTML_TEMPLATE_HPP_
#define SRC_HTML_TEMPLATE_HPP_

#include "FalcoConfig.hpp"
#include <string>

constexpr auto falco_config_html_template = R"(<html>
<head>
    <meta charset="utf-8">
    <script>(function(){var d=document.documentElement,t;try{t=localStorage.getItem('falco-theme');}catch(e){}if(!t)t='dark';d.setAttribute('data-theme',t);})();</script>
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">

	<title>
     {{filename}} - report
	</title>
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
<link href="https://stackpath.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css" rel="stylesheet" integrity="sha384-wvfXpqpZZVQGK6TAh5PVlGOfQNHSoD2xbE+QkPxCAFlNEevoEH3Sl0sibVcOQVnN" crossorigin="anonymous">
<style type="text/css">
:root {
  --bg:#f5f5f7; --surface:#ffffff; --surface-2:#f5f5f7; --border:#e3e3e8;
  --text:#1d1d1f; --muted:#6e6e73; --brand:#0071e3; --brand-2:#42a5ff;
  --pass:#34c759; --warn:#ff9f0a; --fail:#ff3b30;
  --header-bg:rgba(255,255,255,.72);
  --shadow:0 8px 30px rgba(0,0,0,.06),0 1px 2px rgba(0,0,0,.05);
  --shadow-lg:0 18px 50px rgba(0,0,0,.10);
  --radius:18px; --sidebar-w:19em;
}
* { box-sizing:border-box; }
html { scroll-behavior:smooth; }
body {
  color:var(--text); margin:0; padding:0; background:var(--bg);
  font-family:-apple-system,BlinkMacSystemFont,'SF Pro Display','SF Pro Text','Segoe UI',Roboto,Helvetica,Arial,sans-serif;
  font-size:15px; line-height:1.55; -webkit-font-smoothing:antialiased;
  letter-spacing:-.006em;
  -webkit-print-color-adjust:exact; print-color-adjust:exact;
}
a { color:var(--brand); text-decoration:none; }
a:hover { text-decoration:underline; }
div.header {
  position:sticky; top:0; z-index:20; display:flex; align-items:center;
  justify-content:space-between; gap:1em; flex-wrap:wrap; width:100%; margin:0;
  background:var(--header-bg); border-bottom:1px solid var(--border);
  -webkit-backdrop-filter:saturate(180%) blur(20px); backdrop-filter:saturate(180%) blur(20px);
  padding:.9em 1.6em; font-size:100%;
}
#header_title {
  float:none; display:inline-flex; align-items:center;
  font-size:1.6rem; font-weight:700; letter-spacing:-.022em; color:var(--text);
}
#header_filename {
  float:none; text-align:right; font-size:.82rem; color:var(--muted);
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
}
div.header h3 { font-size:50%; margin:0; }
div.summary {
  position:fixed; top:4.6em; left:0; width:var(--sidebar-w);
  height:calc(100vh - 6.6em); overflow-y:auto; margin:0; padding:1em .9em;
}
div.summary h2 {
  font-size:.72rem; text-transform:uppercase; letter-spacing:.08em;
  color:var(--muted); font-weight:700; margin:.2em 0 .7em .5em;
}
div.summary ul { padding:0; margin:0; list-style:none; }
div.summary ul li img { display:none; }
div.summary ul li a {
  display:flex; align-items:center; gap:.6em; padding:.5em .7em;
  margin-bottom:2px; border-radius:8px; color:var(--text);
  font-size:.9rem; font-weight:500; transition:background .12s;
}
div.summary ul li a:hover { background:var(--surface-2); }
div.summary ul li a::before {
  content:""; flex:0 0 auto; width:9px; height:9px; border-radius:50%; background:var(--muted);
}
div.summary ul li a.pass::before { background:var(--pass); }
div.summary ul li a.warn::before { background:var(--warn); }
div.summary ul li a.fail::before { background:var(--fail); }
div.main { margin-left:var(--sidebar-w); padding:1.5em 2em 4.5em; max-width:1180px; }
div.module {
  background:var(--surface); border:1px solid var(--border); border-radius:var(--radius);
  box-shadow:var(--shadow); padding:1.35em 1.7em 1.6em; margin-bottom:1.5em;
  transition:box-shadow .2s ease, transform .2s ease;
}
div.module:hover { box-shadow:var(--shadow-lg); transform:translateY(-2px); }
div.module h2 {
  display:flex; align-items:center; gap:.55em; flex-wrap:wrap; clear:none;
  font-size:1.12rem; font-weight:650; color:var(--text);
  margin:0 0 1em; padding:0 0 .7em; border-bottom:1px solid var(--border);
}
div.module h2::before {
  content:""; flex:0 0 auto; width:11px; height:11px; border-radius:50%; background:var(--muted);
}
div.module h2.pass::before { background:var(--pass); }
div.module h2.warn::before { background:var(--warn); }
div.module h2.fail::before { background:var(--fail); }
div.module h2.pass, div.module h2.warn, div.module h2.fail { color:var(--text); }
div.module h2.collapsible { cursor:pointer; user-select:none; }
div.module h2.collapsible::after {
  content:""; margin-left:auto; flex:0 0 auto; width:.5em; height:.5em; margin-top:-.15em;
  border-right:2px solid var(--muted); border-bottom:2px solid var(--muted);
  transform:rotate(45deg); transition:transform .15s;
}
div.module.collapsed h2.collapsible { margin-bottom:0; padding-bottom:0; border-bottom:0; }
div.module.collapsed h2.collapsible::after { transform:rotate(-45deg); margin-top:.15em; }
div.module.collapsed { padding-bottom:1em; }
div.module.collapsed > :not(h2) { display:none; }
a.toggle-all {
  display:inline-block; margin:.1em 0 .7em .5em; font-size:.72rem;
  color:var(--brand); cursor:pointer; font-weight:600;
}
p.module-desc {
  margin:-.3em 0 1.15em; color:var(--muted); font-size:.86rem;
  line-height:1.5; max-width:78ch;
}
p.module-desc a { font-weight:600; white-space:nowrap; }
div.module > div[id] { width:100%; }
.pass { color:var(--pass); } .warn { color:var(--warn); } .fail { color:var(--fail); }
table {
  margin:.4em 0; border-collapse:collapse; min-width:55%; text-align:left;
  font-size:.87rem; border:1px solid var(--border); border-radius:10px; overflow:hidden;
}
th {
  background:var(--surface-2); color:var(--muted); text-align:left;
  text-transform:uppercase; font-size:.71rem; letter-spacing:.04em; font-weight:700;
  padding:.6em .9em; border-bottom:1px solid var(--border);
}
td {
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace; text-align:left;
  background:var(--surface); color:var(--text); padding:.52em .9em;
  border-bottom:1px solid var(--border);
}
tr:last-child td { border-bottom:0; }
tr:hover td { background:var(--surface-2); }
div.footer {
  position:fixed; bottom:0; left:0; width:100%; z-index:20; margin:0;
  background:var(--header-bg); border-top:1px solid var(--border);
  -webkit-backdrop-filter:saturate(180%) blur(20px); backdrop-filter:saturate(180%) blur(20px);
  color:var(--muted); font-size:.8rem; padding:.55em 1.6em;
}
div.summary ul li a.pass::before { box-shadow:0 0 0 3px rgba(52,199,89,.16); }
div.summary ul li a.warn::before { box-shadow:0 0 0 3px rgba(255,159,10,.16); }
div.summary ul li a.fail::before { box-shadow:0 0 0 3px rgba(255,59,48,.16); }
div.module h2.pass::before { box-shadow:0 0 0 4px rgba(52,199,89,.15); }
div.module h2.warn::before { box-shadow:0 0 0 4px rgba(255,159,10,.15); }
div.module h2.fail::before { box-shadow:0 0 0 4px rgba(255,59,48,.15); }
img.indented { margin-left:3em; } p { margin-top:0; }
@media (max-width:820px) {
  div.summary { position:static; width:auto; height:auto; margin:.8em;
    background:var(--surface); border:1px solid var(--border); border-radius:var(--radius); }
  div.summary ul { display:flex; flex-wrap:wrap; gap:.25em; }
  div.main { margin-left:0; padding:1em; }
}
@media print {
  div.summary, div.footer { position:static; width:auto; height:auto; }
  div.main { margin-left:0; } div.module { break-inside:avoid; box-shadow:none; }
  .modebar { display:none !important; }
}
:root[data-theme="dark"] {
  --bg:#000000; --surface:#1c1c1e; --surface-2:#2c2c2e; --border:#38383a;
  --text:#f5f5f7; --muted:#98989d; --brand:#0a84ff; --brand-2:#64d2ff;
  --pass:#30d158; --warn:#ff9f0a; --fail:#ff453a;
  --header-bg:rgba(28,28,30,.72);
  --shadow:0 8px 30px rgba(0,0,0,.5),0 1px 2px rgba(0,0,0,.6);
  --shadow-lg:0 18px 50px rgba(0,0,0,.6);
}
button.theme-toggle {
  border:1px solid var(--border); background:var(--surface-2); color:var(--text);
  border-radius:999px; width:2.1em; height:2.1em; font-size:1rem; line-height:1;
  cursor:pointer; display:inline-flex; align-items:center; justify-content:center;
  padding:0; margin-left:.15em; flex:0 0 auto; transition:background .15s, transform .15s;
}
button.theme-toggle:hover { background:var(--brand); color:#fff; transform:scale(1.08); }
select.palette-select {
  border:1px solid var(--border); background:var(--surface-2); color:var(--text);
  border-radius:999px; height:2.1em; font-size:.82rem; font-weight:600; line-height:1;
  cursor:pointer; padding:0 1.9em 0 .9em; margin-left:.35em; flex:0 0 auto;
  font-family:inherit; -webkit-appearance:none; appearance:none;
  background-image:url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='10' height='6'><path d='M1 1l4 4 4-4' fill='none' stroke='%236e6e73' stroke-width='1.6'/></svg>");
  background-repeat:no-repeat; background-position:right .8em center;
}
select.palette-select:hover { border-color:var(--brand); }
</style>
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
<div class="header">
	<div id="header_title">FalcoAmadeus</div>
  <div id="header_filename">{{date}}<br/> {{filename}}
	</div>
</div>
<div class="summary"><h2>Summary</h2>
<ul>
  {{basicstatisticscs}}
  <li><a class="{{passbasicstatistics}}" href="#basicstatistics">
    {{basicstatisticsname}}
  </a></li>
  {{basicstatisticsce}}

  {{perbasesequencequalitycs}}
	<li><a class="{{passperbasesequencequality}}" href="#perbasesequencequality">
    {{perbasesequencequalityname}}</a></li>
  {{perbasesequencequalityce}}

  {{pertilesequencequalitycs}}
	<li><a class="{{passpertilesequencequality}}" href="#pertilesequencequality">{{pertilesequencequalityname}}</a></li>
  {{pertilesequencequalityce}}

  {{persequencequalityscorescs}}
	<li><a class="{{passpersequencequalityscores}}" href="#persequencequalityscores">{{persequencequalityscoresname}}</a></li>
  {{persequencequalityscoresce}}

  {{perbasesequencecontentcs}}
	<li><a class="{{passperbasesequencecontent}}" href="#perbasesequencecontent">{{perbasesequencecontentname}}</a></li>
  {{perbasesequencecontentce}}

  {{persequencegccontentcs}}
	<li><a class="{{passpersequencegccontent}}" href="#persequencegccontent">{{persequencegccontentname}}</a></li>
  {{persequencegccontentce}}

  {{perbasencontentcs}}
	<li><a class="{{passperbasencontent}}" href="#perbasencontent">{{perbasencontentname}}</a></li>
  {{perbasencontentce}}

  {{sequencelengthdistributioncs}}
	<li><a class="{{passsequencelengthdistribution}}" href="#sequencelengthdistribution">{{sequencelengthdistributionname}}</a></li>
  {{sequencelengthdistributionce}}

  {{sequenceduplicationlevelscs}}
	<li><a class="{{passsequenceduplicationlevels}}" href="#sequenceduplicationlevels">{{sequenceduplicationlevelsname}}</a></li>
  {{sequenceduplicationlevelsce}}

  {{overrepresentedsequencescs}}
	<li><a class="{{passoverrepresentedsequences}}" href="#overrepresentedsequences">{{overrepresentedsequencesname}}</a></li>
  {{overrepresentedsequencesce}}

  {{adaptercontentcs}}
	<li><a class="{{passadaptercontent}}" href="#adaptercontent">{{adaptercontentname}}</a></li>
  {{adaptercontentce}}

  {{kmercontentcs}}
	<li><a class="{{passkmercontent}}" href="#kmercontent">{{kmercontentname}}</a></li>
  {{kmercontentce}}


</ul>
</div>
<div class="main">
{{combinednote}}
<div class="module">
  <h2 class="{{passbasicstatistics}}" id="basicstatistics">
    {{basicstatisticsname}}: {{passbasicstatistics}}
  </h2>
  {{basicstatisticsdata}}
</div>

{{perbasesequencequalitycs}}
<div class="module">
	<h2 class="{{passperbasesequencequality}}" id="perbasesequencequality">
    {{perbasesequencequalityname}}: {{passperbasesequencequality}}</h2>
 	<div id="seqbasequalityboxplot"></div>
</div>
{{perbasesequencequalityce}}


{{pertilesequencequalitycs}}
<div class="module">
	<h2 class="{{passpertilesequencequality}}" id="pertilesequencequality">
    {{pertilesequencequalityname}} : {{passpertilesequencequality}}
  </h2>
 	<div id="tilequalityheatmap"></div>
</div>
{{pertilesequencequalityce}}

{{persequencequalityscorescs}}
<div class="module">
	<h2 class="{{passpersequencequalityscores}}" id="persequencequalityscores">
    {{persequencequalityscoresname}} : {{passpersequencequalityscores}}
  </h2>
 	<div id="seqqualitylineplot"></div>
</div>
{{persequencequalityscoresce}}

{{perbasesequencecontentcs}}
<div class="module">
	<h2 class="{{passperbasesequencecontent}}" id="perbasesequencecontent">
    {{perbasesequencecontentname}} : {{passperbasesequencecontent}}
  </h2>
 	<div id="basesequencecontentlineplot"></div>
</div>
{{perbasesequencecontentce}}

{{persequencegccontentcs}}
<div class="module">
	<h2 class="{{passpersequencegccontent}}" id="persequencegccontent">
    {{persequencegccontentname}}: {{passpersequencegccontent}}
  </h2>
 	<div id="sequencegccontentlineplot"></div>
</div>
{{persequencegccontentce}}


{{perbasencontentcs}}
<div class="module">
	<h2 class="{{passperbasencontent}}" id="perbasencontent">
    {{perbasencontentname}} : {{passperbasencontent}}
  </h2>
 	<div id="basencontentlineplot"></div>
</div>
{{perbasencontentce}}

{{sequencelengthdistributioncs}}
<div class="module">
	<h2 class="{{passsequencelengthdistribution}}" id="sequencelengthdistribution">
    {{sequencelengthdistributionname}} : {{passsequencelengthdistribution}}
  </h2>
 	<div id="sequencelengthdistributionlineplot"></div>
</div>
{{sequencelengthdistributionce}}

{{sequenceduplicationlevelscs}}
<div class="module">
	<h2 class="{{passsequenceduplicationlevels}}" id="sequenceduplicationlevels">
    {{sequenceduplicationlevelsname}} : {{passsequenceduplicationlevels}}
  </h2>
  {{pairednote}}
 	<div id="seqduplevelslineplot"></div>
</div>
{{sequenceduplicationlevelsce}}

{{overrepresentedsequencescs}}
<div class="module">
	<h2 class="{{passoverrepresentedsequences}}" id="overrepresentedsequences">
    {{overrepresentedsequencesname}} : {{passoverrepresentedsequences}}</h2>
  {{overrepresentedsequencesdata}}
</div>
{{overrepresentedsequencesce}}

{{adaptercontentcs}}
<div class="module">
  <h2 class="{{passadaptercontent}}" id="adaptercontent">
    {{adaptercontentname}} : {{passadaptercontent}}
  </h2>
 	<div id="adapterlineplot"></div>
</div>
{{adaptercontentce}}

{{kmercontentcs}}
<div class="module">
  <h2 class="{{passkmercontent}}" id="kmercontent">
    {{kmercontentname}} : {{passkmercontent}}
  </h2>
 	<div id="kmerlineplot"></div>
</div>
{{kmercontentce}}


</div>
<div class="footer">FalcoAmadeus {{version}}
</div>
</body>
<script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" 
integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" 
crossorigin="anonymous"></script>

<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js" 
integrity="sha384-UO2eT0CpHqdSJQ6hJty5KVphtPhzWj9WO1clHTMGa3JDZwrnQq4sF86dIHNDz0W1"
crossorigin="anonymous"></script>

<script src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js"
integrity="sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM"
crossorigin="anonymous"></script>

<script>
// [FalcoAmadeus] apply one clean, consistent theme to every Plotly plot by
// wrapping newPlot -- keeps transparent backgrounds, light grids, a tidy font,
// slimmer boxes, an on-hover toolbar and no logo, without touching each call.
(function(){
  if (!window.Plotly) return;
  // Selectable colorblind-safe palettes, with a variant per background so no
  // colour ever disappears into it (the near-black in light mode becomes a light
  // neutral in dark mode, dark hues are brightened). Each palette still never
  // mixes red WITH green: 'ocean' has no red, 'sunset'/'ink' have no green.
  // Each palette: categorical light/dark arrays, plus a colourblind-safe
  // diverging pair (warm,cool) used for heatmaps so every plot follows the
  // chosen scheme (lines, box and heatmap alike).
  window.FALCO_PALETTES={
    ocean:{light:['#0071e3','#00c7be','#5856d6','#ff9500','#1d1d1f','#8e8e93'],
           dark: ['#4aa3ff','#40e0d6','#a29bff','#ffb340','#f2f2f7','#c7c7cc'],
           divLight:['#ff9500','#0071e3'],divDark:['#ffb340','#4aa3ff']},
    sunset:{light:['#ff9500','#ff2d55','#af52de','#a2845e','#1d1d1f','#0071e3'],
            dark: ['#ffb340','#ff6482','#d68bff','#d6b48c','#f2f2f7','#64d2ff'],
            divLight:['#ff9500','#af52de'],divDark:['#ffb340','#d68bff']},
    ink:{light:['#1d1d1f','#0071e3','#ff9500','#af52de','#8a6d1f','#6e6e73'],
         dark: ['#f2f2f7','#4aa3ff','#ffb340','#d68bff','#d4b06a','#aeaeb2'],
         divLight:['#ff9500','#0071e3'],divDark:['#ffb340','#4aa3ff']}
  };
  function palColors(name){
    var p=window.FALCO_PALETTES[name]||window.FALCO_PALETTES.ocean;
    return (document.documentElement.getAttribute('data-theme')==='dark')?p.dark:p.light;
  }
  var PAL='ocean'; try{var _p=localStorage.getItem('falco-palette'); if(_p&&window.FALCO_PALETTES[_p])PAL=_p;}catch(e){}
  var FONT={family:'-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif',size:12,color:'#3a4553'};
  // Grid is visible by default (a clear reference grid on every plot).
  var AX={gridcolor:'#d6dbe3',zerolinecolor:'#c4cbd4',linecolor:'#b7bec9',ticks:'outside',tickcolor:'#b7bec9',tickfont:{size:11,color:'#697586'},gridwidth:1,automargin:true};
  // Legend as a horizontal row BELOW the figure so it never collides with the
  // modebar (camera/zoom toolbar), which Plotly pins to the top-right corner.
  var LEGEND={orientation:'h',xanchor:'center',x:0.5,yanchor:'top',y:-0.22,borderwidth:0,font:{size:11}};
  // Force a solid, readable hover label instead of Plotly's default (which tints
  // it with the trace colour and can render same-colour text on same-colour bg).
  var HOVER={bgcolor:'#ffffff',bordercolor:'#d6dbe3',font:{family:FONT.family,size:12,color:'#1d1d1f'}};
  // 'x unified' collapses all series at a position into ONE box that honours the
  // hoverlabel bg/font, instead of per-point boxes whose colours fight the theme.
  var THEME={paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',colorway:palColors(PAL).slice(),boxgap:0.28,boxgroupgap:0.2,legend:LEGEND,hoverlabel:HOVER,hovermode:'x unified'};
  var CFG={displaylogo:false,responsive:true,displayModeBar:'hover',modeBarButtonsToRemove:['select2d','lasso2d','autoScale2d','toggleSpikelines','hoverClosestCartesian','hoverCompareCartesian']};
  var orig=Plotly.newPlot;
  Plotly.newPlot=function(id,data,layout,config){
    layout=layout||{};
    var L=Object.assign({},THEME,layout);
    L.font=Object.assign({},FONT,layout.font);
    L.xaxis=Object.assign({},AX,layout.xaxis);
    L.yaxis=Object.assign({},AX,layout.yaxis);
    L.legend=Object.assign({},LEGEND,layout.legend);
    L.margin=Object.assign({t:34,r:18,b:70,l:60},layout.margin);
    return orig.call(Plotly,id,data,L,Object.assign({},CFG,config));
  };
})();

if (document.getElementById('seqbasequalityboxplot') !== null) {
  Plotly.newPlot('seqbasequalityboxplot', [
   {{perbasesequencequalitydata}}
  ], {
    showlegend: false,
    boxgap: 0.45,
    xaxis : {title : 'Base position'},
    yaxis : {title : 'Phred quality', range : [0, 42]},
    shapes : [
      {type:'rect',xref:'paper',yref:'y',x0:0,x1:1,y0:28,y1:42,fillcolor:'rgba(52,199,89,0.20)',line:{width:0},layer:'below'},
      {type:'rect',xref:'paper',yref:'y',x0:0,x1:1,y0:20,y1:28,fillcolor:'rgba(255,159,10,0.22)',line:{width:0},layer:'below'},
      {type:'rect',xref:'paper',yref:'y',x0:0,x1:1,y0:0,y1:20,fillcolor:'rgba(255,59,48,0.20)',line:{width:0},layer:'below'}
    ]
  });
}
if (document.getElementById('tilequalityheatmap') !== null) {
  Plotly.newPlot('tilequalityheatmap', [
   {{pertilesequencequalitydata}}
  ], {
    margin: { t: 0 }, 
    showlegend: false,
    xaxis : {title : 'Base position'},
    yaxis : {title : 'tile', type: 'category'} 
  });
}

if (document.getElementById('seqqualitylineplot') !== null) {
  Plotly.newPlot('seqqualitylineplot', [
   {{persequencequalityscoresdata}}
  ], {
    margin: { t: 0 },
    showlegend: true,
    xaxis : {title : 'Phread quality'},
    yaxis : {title : 'Density'}
  } );
}
if (document.getElementById('basesequencecontentlineplot') !== null) {
  Plotly.newPlot('basesequencecontentlineplot', [
   {{perbasesequencecontentdata}}
  ], {
    margin: { t: 0 },
    showlegend: true,
    xaxis : {title : 'Base position'},
    yaxis : {title : '% sequence content'}
  } );
}
if (document.getElementById('sequencegccontentlineplot') !== null) {
  Plotly.newPlot('sequencegccontentlineplot', [
   {{persequencegccontentdata}}
  ], {
    margin: { t: 0 }, 
    showlegend: true,
    xaxis : {title : '% GC'},
    yaxis : {title : 'Density'}
  } );
}
if (document.getElementById('basencontentlineplot') !== null) {
  Plotly.newPlot('basencontentlineplot', [
   {{perbasencontentdata}}
  ], {
    margin: { t: 0 }, 
    showlegend: true,
    xaxis : {title : 'Base position'},
    yaxis : {title : '% N'}
  } );
}

if (document.getElementById('sequencelengthdistributionlineplot') !== null) {
  Plotly.newPlot('sequencelengthdistributionlineplot', [
   {{sequencelengthdistributiondata}}
  ], {
    margin: { t: 0 }, 
    showlegend: true,
    xaxis : {title : 'Sequence length'},
    yaxis : {title : 'Number of sequences'}
  } );
}
if (document.getElementById('seqduplevelslineplot') !== null) {
  Plotly.newPlot('seqduplevelslineplot', [
   {{sequenceduplicationlevelsdata}}
  ], {
    margin: { t: 0 }, 
    showlegend: true,
    xaxis : {title : 'Duplication rate',
             tickvals : [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],
             ticktext : ['1','2','3','4','5','6','7','8','9','10+','50+','100+','500+','1k+','5k+','10k+']},
    yaxis : {title : '% of sequences'}
  } );
}

if (document.getElementById('adapterlineplot') !== null) {
  Plotly.newPlot('adapterlineplot', [
   {{adaptercontentdata}}
  ], {
    margin: { t: 0 }, 
    showlegend: true,
    xaxis : {title : 'Base position'},
    yaxis : {title : '% sequences with adapter before position'}
  } );
}
if (document.getElementById('kmerlineplot') !== null) {
  Plotly.newPlot('kmerlineplot', [
   {{kmercontentdata}}
  ], {
    margin: { t: 0 }, 
    showlegend: true,
    xaxis : {title : 'Base position'},
    yaxis : {title : 'log2(obs/ exp max)'}
  } );
}

// [FalcoAmadeus] make each module collapsible (click its header to toggle).
// Runs after the plots are drawn, so every plot renders at full width first.
(function(){
  var BASE="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/";
  var INFO={
   basicstatistics:["Summary counts for the file - total reads, read-length range, %GC and the inferred quality encoding. Check these match what you expect for the library.","1%20Basic%20Statistics.html"],
   perbasesequencequality:["Box plots of the quality score at each base position. Quality usually declines toward the end of the read; medians in the low range mean you should consider quality trimming.","2%20Per%20Base%20Sequence%20Quality.html"],
   persequencequalityscores:["Distribution of the mean quality per read. One peak at high quality is good; a low-quality shoulder means a subset of reads is poor.","3%20Per%20Sequence%20Quality%20Scores.html"],
   perbasesequencecontent:["Proportion of A/C/G/T at each position - expected flat and parallel for random libraries. Strong bias (often at the read start) is normal for RNA-seq, amplicon and bisulfite data; see the --assay option.","4%20Per%20Base%20Sequence%20Content.html"],
   persequencegccontent:["Observed GC distribution across reads versus a modeled normal. A shifted or bimodal curve can mean contamination, or is expected for bisulfite and targeted libraries.","5%20Per%20Sequence%20GC%20Content.html"],
   perbasencontent:["Percentage of ambiguous N base calls per position. More than a fraction of a percent, especially if localized, indicates a sequencing issue.","6%20Per%20Base%20N%20Content.html"],
   sequencelengthdistribution:["Distribution of read lengths. A single length is typical of untrimmed Illumina reads; a spread is expected after trimming or for small-RNA / long-read libraries.","7%20Sequence%20Length%20Distribution.html"],
   sequenceduplicationlevels:["Fraction of reads seen more than once - a proxy for PCR/optical duplication and library complexity. High duplication is expected for amplicon, small-RNA and highly-expressed RNA-seq. For paired-end data use the paired report; single-end overcounts.","8%20Duplicate%20Sequences.html"],
   overrepresentedsequences:["Individual sequences that make up an unusually large share of the library - often adapters or contaminants, or (expected) the amplicons/miRNAs themselves.","9%20Overrepresented%20Sequences.html"],
   adaptercontent:["Cumulative fraction of reads containing known adapters by position. A rise near the read end means adapter read-through - trim before downstream analysis.","10%20Adapter%20Content.html"],
   kmercontent:["Positional enrichment of short k-mers. Spikes flag adapters, primers or other position-specific sequence bias.","11%20Kmer%20Content.html"],
   pertilesequencequality:["Quality deviation per flow-cell tile relative to the run average. Persistent blue/red bands flag a localized sequencing problem on specific tiles.","12%20Per%20Tile%20Sequence%20Quality.html"]
  };
  var mods=document.querySelectorAll("div.module");
  mods.forEach(function(m){
    var h=m.querySelector("h2"); if(!h) return;
    var info=INFO[h.id];
    if(info){
      var p=document.createElement("p"); p.className="module-desc";
      p.innerHTML=info[0]+" <a href='"+BASE+info[1]+"' target='_blank' rel='noopener'>FastQC help ↗</a>";
      h.parentNode.insertBefore(p,h.nextSibling);
    }
    h.classList.add("collapsible");
    h.addEventListener("click",function(){ m.classList.toggle("collapsed"); });
  });
  var sum=document.querySelector("div.summary h2");
  if(sum){
    var btn=document.createElement("a"); btn.className="toggle-all"; btn.textContent="collapse all";
    btn.addEventListener("click",function(){
      var collapse=btn.textContent==="collapse all";
      mods.forEach(function(m){ m.classList.toggle("collapsed",collapse); });
      btn.textContent=collapse?"expand all":"collapse all";
    });
    sum.parentNode.insertBefore(btn,sum.nextSibling);
  }
})();

// [FalcoAmadeus] dark-mode toggle (persisted) that also re-themes the plots
(function(){
  var root=document.documentElement;
  // Theme the plot chrome: a clear reference grid on every plot by default,
  // recoloured for dark mode. (No separate contrast mode -- the grid is default.)
  function applyPlotStyle(){
    if(!window.Plotly) return;
    var dark=root.getAttribute('data-theme')==='dark';
    var fg=dark?'#e6ebf2':'#2a3441', grid=dark?'#3a4451':'#d6dbe3',
        line=dark?'#4a5462':'#b7bec9', tick=dark?'#9aa6b2':'#697586';
    // readable hover label: solid surface bg + contrasting text for this theme
    var hbg=dark?'#2c2c2e':'#ffffff', hbd=dark?'#48484a':'#d6dbe3',
        htx=dark?'#f5f5f7':'#1d1d1f';
    document.querySelectorAll('.js-plotly-plot').forEach(function(d){
      Plotly.relayout(d,{'font.color':fg,'xaxis.gridcolor':grid,'yaxis.gridcolor':grid,
        'xaxis.gridwidth':1,'yaxis.gridwidth':1,
        'xaxis.linecolor':line,'yaxis.linecolor':line,
        'xaxis.tickcolor':line,'yaxis.tickcolor':line,
        'xaxis.tickfont.color':tick,'yaxis.tickfont.color':tick,'xaxis.zerolinecolor':grid,
        'yaxis.zerolinecolor':grid,'xaxis.title.font.color':fg,'yaxis.title.font.color':fg,
        'legend.font.color':fg,
        'hoverlabel.bgcolor':hbg,'hoverlabel.bordercolor':hbd,'hoverlabel.font.color':htx});
    });
  }
  // Recolour every categorical plot from the chosen palette, cycling colours
  // across each plot's solid line/bar traces (dashed reference lines, the box
  // plot and heatmaps keep their own colours).
  function hexToRgba(h,a){var m=(h||'').replace('#','');
    if(m.length===3)m=m[0]+m[0]+m[1]+m[1]+m[2]+m[2];
    var r=parseInt(m.substr(0,2),16),g=parseInt(m.substr(2,2),16),b=parseInt(m.substr(4,2),16);
    return 'rgba('+r+','+g+','+b+','+a+')';}
  function applyPalette(name){
    if(!window.Plotly||!window.FALCO_PALETTES) return;
    // theme-aware: pick the light/dark variant for the active background
    var dark=root.getAttribute('data-theme')==='dark';
    var p=window.FALCO_PALETTES[name]||window.FALCO_PALETTES.ocean;
    var pal=dark?p.dark:p.light, dv=dark?p.divDark:p.divLight,
        neutral=dark?'#3a3f47':'#eef0f3';
    document.querySelectorAll('.js-plotly-plot').forEach(function(gd){
      if(!gd.data) return;
      var k=0;
      gd.data.forEach(function(tr,i){
        var t=tr.type||'scatter';
        if(t==='heatmap'||t==='histogram2d'||t==='histogram2dcontour'||t==='contour'){
          // diverging scale from the palette (warm -> neutral -> cool)
          try{Plotly.restyle(gd,{colorscale:[[[0,dv[0]],[0.5,neutral],[1,dv[1]]]]},[i]);}catch(e){}
        } else if(t==='box'){
          var bc=pal[0];
          try{Plotly.restyle(gd,{fillcolor:[hexToRgba(bc,0.55)],'line.color':[bc],'marker.color':[bc]},[i]);}catch(e){}
        } else {
          var dashed=tr.line&&tr.line.dash&&tr.line.dash!=='solid';
          if(dashed) return;
          var c=pal[k%pal.length]; k++;
          try{Plotly.restyle(gd,{'line.color':[c],'marker.color':[c]},[i]);}catch(e){}
        }
      });
    });
  }
  var btn=document.createElement('button'); btn.className='theme-toggle';
  btn.title='Toggle light / dark';
  function set(t){ root.setAttribute('data-theme',t);
    try{localStorage.setItem('falco-theme',t);}catch(e){}
    btn.textContent=(t==='dark')?'☀':'☾'; applyPlotStyle();
    // re-map the palette to the new background so no colour is lost in it
    if(typeof sel!=='undefined'&&sel) applyPalette(sel.value); }
  btn.onclick=function(){ set(root.getAttribute('data-theme')==='dark'?'light':'dark'); };
  // palette picker
  var sel=document.createElement('select'); sel.className='palette-select';
  sel.title='Plot colour scheme';
  var NAMES={ocean:'Ocean',sunset:'Sunset',ink:'Ink'};
  ['ocean','sunset','ink'].forEach(function(k){
    var o=document.createElement('option'); o.value=k; o.textContent=NAMES[k]; sel.appendChild(o); });
  var curPal='ocean'; try{var sp=localStorage.getItem('falco-palette'); if(sp&&NAMES[sp])curPal=sp;}catch(e){}
  sel.value=curPal;
  sel.onchange=function(){ try{localStorage.setItem('falco-palette',sel.value);}catch(e){} applyPalette(sel.value); };
  var title=document.querySelector('#header_title');
  if(title){ title.parentNode.insertBefore(btn,title.nextSibling); btn.parentNode.insertBefore(sel,btn.nextSibling); }
  set(root.getAttribute('data-theme')||'dark');
  window.addEventListener('load',function(){ applyPlotStyle(); applyPalette(sel.value); });
  // Expose theme/palette setters so an embedding page (the paired combined
  // report) can drive this report from one shared header.
  window.falcoSetTheme=function(t){ set(t); };
  window.falcoSetPalette=function(p){ if(!NAMES[p])return; sel.value=p;
    try{localStorage.setItem('falco-palette',p);}catch(e){} applyPalette(p); };
  window.falcoTheme=function(){ return root.getAttribute('data-theme')||'dark'; };
  window.falcoPalette=function(){ return sel.value; };
})();

</script>
</html>;)";

#endif
