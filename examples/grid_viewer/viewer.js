"use strict";
(() => {
  const entries = JSON.parse(document.getElementById("scene-data").textContent);
  const $ = id => document.getElementById(id);
  const ns = "http://www.w3.org/2000/svg";
  let selected = Math.max(0, entries.findIndex(e => e.scene.technique === "ThreeHouseMatching"));
  let frameIndex = 0;
  const completeLibrary = new Set(entries.filter(e=>e.proof_bundle).map(e=>e.scene.technique)).size === 48;
  if(!completeLibrary)$("coverage-note").textContent="This file shows selected steps from one unique puzzle. Each step has downloadable independent evidence.";
  const statusLabels = {checked_certificate:"✓ Certificate verified", checked_deduction:"✓ Deduction verified", checked_pattern_and_deduction:"✓ Pattern and deduction verified", metadata_only:"? Metadata only", rejected:"! Evidence rejected", rules_satisfied:"✓ Placed values valid", rule_violation:"! Rule violation"};
  const palette = {
    pattern:["#e6edff","#244dae"],base:["#d9e7ff","#224b9b"],cover:["#e4e9fc","#4d3b94"],
    fin:["#ffebbf","#825100"],floor:["#d9e7ff","#224b9b"],roof:["#ffebbf","#825100"],
    on:["#d9e7ff","#224b9b"],off:["#ffe0e4","#a1273e"],assumption:["#ffebbf","#825100"],
    elimination:["#ffe0e4","#a1273e"],placement:["#c9f0df","#076246"],conflict:["#ffcbd4","#9b1730"]
  };
  const priority = ["pattern","base","cover","floor","roof","on","off","fin","placement","elimination","assumption","conflict"];
  function svg(tag, attrs, text) {
    const el = document.createElementNS(ns, tag);
    for (const [k,v] of Object.entries(attrs)) el.setAttribute(k, v);
    if (text !== undefined) el.textContent = text;
    return el;
  }
  function topRole(marks) { return marks.reduce((r,h) => priority.indexOf(h.role) > priority.indexOf(r) ? h.role : r, "pattern"); }
  function point(candidate) {
    return [24+(candidate.cell%9)*64+((candidate.digit-1)%3+.5)*64/3,
      24+Math.floor(candidate.cell/9)*64+(Math.floor((candidate.digit-1)/3)+.5)*64/3];
  }
  function inspect(cell, state, marks) {
    $("cell-title").textContent = `r${Math.floor(cell/9)+1}c${cell%9+1}`;
    const roles = marks.map(h => `${h.digit ? "Digit "+h.digit+": " : "Cell: "}${h.role.replaceAll("_"," ")}`);
    $("cell-description").textContent = [state.value ? `Value ${state.value}${state.given ? " (given)" : ""}.` : `Candidates: ${state.candidates.join(", ") || "none"}.`, ...roles].join(" ");
  }
  function draw(frame, original) {
    const board = $("grid"); board.replaceChildren();
    const cells = original ? Array.from(entries[selected].original_puzzle, c => ({value:/[1-9]/.test(c)?Number(c):null,given:/[1-9]/.test(c),candidates:[]})) : frame.cells;
    const marks = original ? [] : frame.highlights;
    board.setAttribute("aria-label", `${original ? "Original givens" : frame.title}. Sudoku grid. Select a cell for its values and roles.`);
    for (let n=0;n<9;n++) {
      board.append(svg("text",{x:56+n*64,y:14,"text-anchor":"middle",fill:"#5a687b","font-size":11},`c${n+1}`));
      board.append(svg("text",{x:12,y:60+n*64,"text-anchor":"middle",fill:"#5a687b","font-size":11},`r${n+1}`));
    }
    cells.forEach((state,cell) => {
      const x=24+(cell%9)*64, y=24+Math.floor(cell/9)*64;
      const annotations=marks.filter(h=>h.cell===cell);
      const cellMarks=annotations.filter(h=>h.digit===null);
      const color=cellMarks.length ? palette[topRole(cellMarks)][0] : "#fff";
      const g=svg("g",{tabindex:0,role:"button","aria-label":`r${Math.floor(cell/9)+1}c${cell%9+1}, ${state.value ? "value "+state.value : "candidates "+(state.candidates.join(" ")||"none")}`});
      g.append(svg("rect",{class:"cell-bg",x,y,width:64,height:64,fill:color,stroke:"#d4dbe7","stroke-width":.7}));
      if (state.value) {
        const valueMarks=annotations.filter(h=>h.digit===state.value);
        const role=valueMarks.length ? topRole(valueMarks) : null;
        if(role) g.append(svg("rect",{x:x+8,y:y+8,width:48,height:48,rx:6,fill:palette[role][0]}));
        g.append(svg("text",{x:x+32,y:y+44,"text-anchor":"middle","font-size":34,"font-weight":state.given?650:400,fill:role?palette[role][1]:state.given?"#1d2b3d":"#244dae"},state.value));
      } else {
        for(let d=1;d<=9;d++) {
          const candidateMarks=annotations.filter(h=>h.digit===d);
          const present=state.candidates.includes(d);
          if(!present&&!candidateMarks.length) continue;
          const [cx,cy]=point({cell,digit:d});
          const role=candidateMarks.length?topRole(candidateMarks):null;
          if(role) g.append(svg("rect",{x:cx-9,y:cy-9,width:18,height:18,rx:4,fill:palette[role][0],stroke:role==="assumption"?palette[role][1]:"none","stroke-dasharray":role==="assumption"?"2 2":""}));
          g.append(svg("text",{x:cx,y:cy+4.5,"text-anchor":"middle","font-size":13,"font-weight":role?700:400,fill:role?palette[role][1]:"#58677b",opacity:present?1:.75},d));
          if(["elimination","off"].includes(role)) g.append(svg("path",{d:`M${cx-7} ${cy+7}L${cx+7} ${cy-7}`,stroke:palette[role][1],"stroke-width":1.6}));
          if(role==="placement") g.append(svg("path",{d:`M${cx-7} ${cy+6}L${cx-4} ${cy+8}L${cx+6} ${cy-7}`,fill:"none",stroke:palette[role][1],"stroke-width":1.3}));
        }
      }
      if(cellMarks.some(h=>h.role==="conflict")) {
        g.append(svg("rect",{x:x+2,y:y+2,width:60,height:60,fill:"none",stroke:"#9b1730","stroke-width":2}));
        g.append(svg("text",{x:x+58,y:y+13,"text-anchor":"end",fill:"#9b1730","font-size":12,"font-weight":800},"!"));
      }
      g.addEventListener("click",()=>inspect(cell,state,annotations));
      g.addEventListener("keydown",e=>{if(["Enter"," "].includes(e.key)){e.preventDefault();inspect(cell,state,annotations);}});
      board.append(g);
    });
    for(let n=0;n<=9;n++) {
      const width=n%3===0?2.3:.6;
      board.append(svg("line",{x1:24+n*64,y1:24,x2:24+n*64,y2:600,stroke:"#23354d","stroke-width":width,"pointer-events":"none"}));
      board.append(svg("line",{x1:24,y1:24+n*64,x2:600,y2:24+n*64,stroke:"#23354d","stroke-width":width,"pointer-events":"none"}));
    }
    if(!original) for(const link of frame.links) {
      const [x1,y1]=point(link.from),[x2,y2]=point(link.to);
      board.append(svg("line",{x1,y1,x2,y2,stroke:link.kind==="strong"?"#244dae":"#a1273e","stroke-width":2,opacity:.7,"stroke-dasharray":link.kind==="weak"?"5 4":"","pointer-events":"none"}));
    }
    $("grid-text").textContent=cells.map((s,i)=>`r${Math.floor(i/9)+1}c${i%9+1}: ${s.value || "["+s.candidates.join(",")+"]"}`).join("\n");
  }
  function renderFrame() {
    const entry=entries[selected],scene=entry.scene,frame=scene.frames[frameIndex];
    $("frame-kind").textContent=frame.kind.toUpperCase(); $("frame-title").textContent=frame.title;
    $("frame-description").textContent=frame.description; $("frame-count").textContent=`${frameIndex+1} / ${scene.frames.length}`;
    $("timeline").value=frameIndex; $("frame-select").value=frameIndex;
    $("previous").disabled=frameIndex===0; $("next").disabled=frameIndex===scene.frames.length-1;
    $("board-mode").textContent=$("original").checked?"Original givens · proof highlights hidden":"Current proof frame";
    $("cell-title").textContent="Inspect a cell";$("cell-description").textContent="Select a cell on the board.";
    draw(frame,$("original").checked);
  }
  function library() {
    const query=$("search").value.toLowerCase();$("examples").replaceChildren();let count=0;
    entries.forEach((entry,i)=>{
      if(!`${entry.title} ${entry.scene.title} ${entry.origin} ${entry.scene.evidence.status}`.toLowerCase().includes(query))return;
      count++;
      const button=document.createElement("button");button.type="button";button.setAttribute("aria-current",String(i===selected));
      const title=document.createElement("span");title.textContent=entry.title;button.append(title);
      const meta=document.createElement("small");meta.textContent=statusLabels[entry.scene.evidence.status];button.append(meta);
      button.addEventListener("click",()=>select(i));$("examples").append(button);
    });
    const techniques=new Set(entries.map(e=>e.scene.technique).filter(Boolean));
    $("coverage").textContent=`${count} examples shown · ${techniques.size} / 48 techniques illustrated. ${completeLibrary ? "Every technique has a real puzzle example." : "Selected puzzle steps; all technique IDs are supported."}`;
  }
  function select(i) {
    selected=i;frameIndex=0;const entry=entries[i],scene=entry.scene;
    $("technique").textContent=scene.title;$("example-title").textContent=entry.title;
    $("status").className=`status ${scene.evidence.status}`;$("status").textContent=statusLabels[scene.evidence.status];
    $("evidence").textContent=scene.evidence.message;$("origin").textContent=entry.origin;
    $("source").hidden=!entry.source_url;
    if(entry.source_url)$("source").href=entry.source_url;
    $("export-proof").disabled=!entry.proof_bundle;
    $("original").checked=false;$("original").disabled=entry.original_puzzle.length!==81;
    $("timeline").max=scene.frames.length-1;$("frame-select").replaceChildren();
    scene.frames.forEach((f,index)=>{const option=document.createElement("option");option.value=index;option.textContent=`${index+1}. ${f.title}`;$("frame-select").append(option);});
    $("truncated").hidden=!scene.trace_truncated;library();renderFrame();
  }
  function move(delta){frameIndex=Math.max(0,Math.min(entries[selected].scene.frames.length-1,frameIndex+delta));renderFrame();}
  $("previous").addEventListener("click",()=>move(-1));$("next").addEventListener("click",()=>move(1));
  for(const id of ["timeline","frame-select"])$(id).addEventListener("input",()=>{frameIndex=Number($(id).value);renderFrame();});
  $("original").addEventListener("change",renderFrame);$("search").addEventListener("input",library);
  document.addEventListener("keydown",e=>{if(["INPUT","SELECT","BUTTON","SUMMARY"].includes(document.activeElement.tagName))return;if(e.key==="ArrowRight"){e.preventDefault();move(1);}if(e.key==="ArrowLeft"){e.preventDefault();move(-1);}});
  $("export").addEventListener("click",()=>{const blob=new Blob([JSON.stringify(entries[selected],null,2)],{type:"application/json"});const url=URL.createObjectURL(blob);const a=document.createElement("a");a.href=url;a.download="sudoku-scene.json";a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
  $("export-proof").addEventListener("click",()=>{const proof=entries[selected].proof_bundle;if(!proof)return;const blob=new Blob([JSON.stringify(proof,null,2)],{type:"application/json"});const url=URL.createObjectURL(blob);const a=document.createElement("a");a.href=url;a.download="sudoku-proof.json";a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
  select(selected);
})();
