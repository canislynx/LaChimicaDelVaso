(function(){
  "use strict";

  /* ============================================================
     MODELLO CHIMICO
     Tutte le concentrazioni interne sono in mol/L (M).
     n_i(h)  = C_i * Ka_i / (Ka_i + [H+])           contributo di una coppia
     n_w(h)  = Kw/[H+] - [H+]                        contributo dell'acqua
     n_s(h)  = ∓C                                    acido/base forte (costante)
     beta    = d n / d pH   (larghezza del vaso, sempre ≥ 0)
     Un vaso poliprotico somma i contributi di più coppie (pKa può essere un array).
     Ogni vaso porta con sé una quantità assoluta e conservata di titolante già aggiunto
     (v.tauExtra, mol/L equivalenti dalla sua "nascita" a pH 7). Il livello di un gruppo
     connesso risolve:
       sum_i n_i(h)  =  sum_i [ n_i(7) + tauExtra_i ]
     così unire due gruppi già titolati separatamente SOMMA le loro quantità — chimicamente
     corretto — invece di ripartire da un nuovo zero per il gruppo appena unito.
     ============================================================ */
  const KW = 1e-14;
  const H = pH => Math.pow(10, -pH);
  const RT_LN10_KJ = 5.708; // kJ/mol per unità di pH, a 298 K

  function pairN(Ka, C, pH){ const h = H(pH); return C * Ka / (Ka + h); }
  function pairBeta(Ka, C, pH){ const h = H(pH); return Math.LN10 * C * Ka * h / Math.pow(Ka + h, 2); }
  function n_w(pH){ const h = H(pH); return KW/h - h; }
  function beta_w(pH){ const h = H(pH); return Math.LN10 * (h + KW/h); }

  // HCl e NaOH sono vasi normalissimi, nessun caso speciale: la stessa formula di una
  // qualunque coppia debole. HCl ha pKa −7, fuori dalla finestra 0-14, quindi è
  // "sempre dissociato" nell'intervallo visibile (beta≈0 lì) con un vero rigonfiamento
  // appena fuori scala. NaOH ha invece pKa 14 — il vero pKa dell'equilibrio H₂O/OH⁻,
  // proprio al bordo della finestra: non essendo spostato fuori scala come HCl, nel suo
  // caso la formula mostra fedelmente anche la parte debole dell'equilibrio (poco
  // dissociato a pH bassi), non un comportamento da base forte su tutto 0-14.
  // diluire un vaso (v.dil, fattore ≥ 1: quante volte il suo volume originale è
  // stato annacquato) non si limita a restringere la sua C: la soluzione diventa
  // sempre più simile ad acqua pura, finché a diluizione "infinita" il suo pH deve
  // tendere a 7. f è la frazione di soluzione originale rimasta (1 = non diluito,
  // →0 = quasi solo acqua aggiunta): il contributo del vaso è una media pesata fra
  // la sua chimica propria e quella dell'acqua, non solo la sua chimica rimpicciolita.
  function vesselN(v, pH){
    if(v.water) return n_w(pH);
    const pKas = Array.isArray(v.pKa) ? v.pKa : [v.pKa];
    let s = 0;
    for(const pk of pKas) s += pairN(Math.pow(10,-pk), v.C, pH);
    const f = 1/(v.dil||1);
    return f===1 ? s : f*s + (1-f)*n_w(pH);
  }
  function vesselBeta(v, pH){
    if(v.water) return beta_w(pH);
    const pKas = Array.isArray(v.pKa) ? v.pKa : [v.pKa];
    let s = 0;
    for(const pk of pKas) s += pairBeta(Math.pow(10,-pk), v.C, pH);
    const f = 1/(v.dil||1);
    return f===1 ? s : f*s + (1-f)*beta_w(pH);
  }
  function primaryPka(v){ return Array.isArray(v.pKa) ? v.pKa[0] : v.pKa; }

  function groupF(members, pH){
    let s = 0;
    for(const v of members) s += vesselN(v, pH);
    return s;
  }

  // il titolante NON è più un unico numero globale condiviso da tutto il banco: ogni
  // vaso porta con sé una quantità assoluta e conservata (v.tauExtra, mol/L equivalenti
  // aggiunti dalla sua "nascita" a pH 7). Il bersaglio di un gruppo è semplicemente la
  // somma dei bersagli dei suoi membri — così unire due gruppi già titolati separatamente
  // SOMMA correttamente le loro quantità (chimicamente corretto), invece di ripartire da
  // un nuovo zero a pH 7 per il gruppo appena unito.
  function groupTarget(members){
    let s = 0;
    for(const v of members) s += vesselN(v, 7) + (v.tauExtra||0);
    return s;
  }

  // dominio di bisezione di un gruppo: normalmente 0-14, esteso a EXTRA_LO/EXTRA_HI se
  // il gruppo contiene solo vasi con pKa fuori scala (HCl, NaOH isolati o fra loro) —
  // così, isolati con una valvola chiusa, possono davvero riempirsi fino al loro pKa reale
  function groupDomain(members){
    let lo = 0, hi = 14;
    members.forEach(v=>{
      if(v.water) return;
      const pKas = Array.isArray(v.pKa) ? v.pKa : [v.pKa];
      pKas.forEach(pk=>{ if(pk<0) lo = EXTRA_LO; if(pk>14) hi = EXTRA_HI; });
    });
    return {lo, hi};
  }

  // risolve groupF(h) = groupTarget(members), per bisezione sul dominio del gruppo
  function solveLevel(members){
    const target = groupTarget(members);
    const {lo:LO, hi:HI} = groupDomain(members);
    let lo = LO, hi = HI;
    let flo = groupF(members, lo) - target;
    let fhi = groupF(members, hi) - target;
    if(flo >= 0) return lo;
    if(fhi <= 0) return hi;
    for(let i=0;i<60;i++){
      const mid = (lo+hi)/2;
      const fm = groupF(members, mid) - target;
      if(fm > 0) hi = mid; else lo = mid;
    }
    return (lo+hi)/2;
  }
  // dato un livello pH desiderato per un gruppo isolato, il tauExtra TOTALE (sommato sui
  // membri) che lo produce, partendo da tauExtra=0 per tutti
  function tauForLevel(members, level){
    return groupF(members, level) - members.reduce((s,v)=>s+vesselN(v,7),0);
  }

  // ogni vaso (HCl/NaOH inclusi, ora coppie come le altre) ha sempre un pH definito
  function isGrounded(members){
    return true;
  }

  /* ============================================================
     STATO
     ============================================================ */
  const PRESETS = [
    {key:"acetico",   label:"Acetico / Acetato", short:"Acetico", pKa:4.76, C:1.0, color:"var(--c-1)", cat:"buffer", group:"Tamponi semplici",
      species:[{acid:"CH₃COOH", base:"CH₃COO⁻"}]},
    {key:"ammonio",    label:"Ammonio / Ammoniaca", short:"Ammonio", pKa:9.25, C:1.0, color:"var(--c-3)", cat:"buffer", group:"Tamponi semplici",
      species:[{acid:"NH₄⁺", base:"NH₃"}]},

    {key:"fosfato1",   label:"Fosfato: 1ª coppia (H₃PO₄/H₂PO₄⁻)", short:"Fosfato 1ª", pKa:2.15, C:1.0, color:"var(--c-2)", cat:"buffer", group:"Fosfato", menu:false,
      species:[{acid:"H₃PO₄", base:"H₂PO₄⁻"}]},
    {key:"fosfato2",   label:"Fosfato: 2ª coppia (H₂PO₄⁻/HPO₄²⁻)", short:"Fosfato 2ª", pKa:7.20, C:1.0, color:"var(--c-2)", cat:"buffer", group:"Fosfato", menu:false,
      species:[{acid:"H₂PO₄⁻", base:"HPO₄²⁻"}]},
    {key:"fosfato3",   label:"Fosfato: 3ª coppia (HPO₄²⁻/PO₄³⁻)", short:"Fosfato 3ª", pKa:12.35, C:1.0, color:"var(--c-2)", cat:"buffer", group:"Fosfato", menu:false,
      species:[{acid:"HPO₄²⁻", base:"PO₄³⁻"}]},
    {key:"fosfatoCompleto", label:"Fosfato: sistema completo (3 coppie)", short:"Fosfato (3)", pKa:[2.15,7.20,12.35], C:1.0, color:"var(--c-2)", cat:"buffer", group:"Fosfato",
      species:[{acid:"H₃PO₄", base:"H₂PO₄⁻"},{acid:"H₂PO₄⁻", base:"HPO₄²⁻"},{acid:"HPO₄²⁻", base:"PO₄³⁻"}]},

    {key:"carbonico1",  label:"Carbonico: 1ª coppia (H₂CO₃/HCO₃⁻)", short:"Carbonico 1ª", pKa:6.35, C:1.0, color:"var(--c-4)", cat:"buffer", group:"Carbonato", menu:false,
      species:[{acid:"H₂CO₃", base:"HCO₃⁻"}]},
    {key:"carbonico2",  label:"Carbonico: 2ª coppia (HCO₃⁻/CO₃²⁻)", short:"Carbonico 2ª", pKa:10.33, C:1.0, color:"var(--c-4)", cat:"buffer", group:"Carbonato", menu:false,
      species:[{acid:"HCO₃⁻", base:"CO₃²⁻"}]},
    {key:"carbonicoCompleto", label:"Carbonico: sistema completo (2 coppie)", short:"Carbonico (2)", pKa:[6.35,10.33], C:1.0, color:"var(--c-4)", cat:"buffer", group:"Carbonato",
      species:[{acid:"H₂CO₃", base:"HCO₃⁻"},{acid:"HCO₃⁻", base:"CO₃²⁻"}]},

    {key:"glicinaCOOH", label:"Glicina: gruppo −COOH (pKa 2.34)", short:"Glicina COOH", pKa:2.34, C:1.0, color:"var(--c-5)", cat:"buffer", group:"Amminoacidi", menu:false,
      species:[{acid:"H₃N⁺CH₂COOH", base:"H₃N⁺CH₂COO⁻"}]},
    {key:"glicinaNH3", label:"Glicina: gruppo −NH₃⁺ (pKa 9.60)", short:"Glicina NH₃⁺", pKa:9.60, C:1.0, color:"var(--c-5)", cat:"buffer", group:"Amminoacidi", menu:false,
      species:[{acid:"H₃N⁺CH₂COO⁻", base:"H₂NCH₂COO⁻"}]},
    {key:"glicinaCompleta", label:"Glicina: sistema completo (2 coppie)", short:"Glicina", pKa:[2.34,9.60], C:1.0, color:"var(--c-5)", cat:"buffer", group:"Amminoacidi",
      species:[{acid:"H₃N⁺CH₂COOH", base:"H₃N⁺CH₂COO⁻"},{acid:"H₃N⁺CH₂COO⁻", base:"H₂NCH₂COO⁻"}]},

    {key:"fenolftaleina", label:"Indicatore: fenolftaleina", short:"Fenolftaleina", pKa:9.10, C:0.00005, color:"var(--c-6)", cat:"indicator", group:"Indicatori",
      acidColor:[223,220,221], baseColor:[196,43,135],
      species:[{acid:"HIn (incolore)", base:"In²⁻ (fucsia)"}]},
    {key:"bromotimolo",   label:"Indicatore: blu di bromotimolo", short:"Bromotimolo", pKa:7.10, C:0.00005, color:"var(--c-7)", cat:"indicator", group:"Indicatori",
      acidColor:[223,185,43], baseColor:[37,101,173],
      species:[{acid:"HIn (giallo)", base:"In⁻ (blu)"}]},
    {key:"metilarancio",  label:"Indicatore: metilarancio", short:"Metilarancio", pKa:3.70, C:0.00005, color:"var(--c-1)", cat:"indicator", group:"Indicatori",
      acidColor:[200,40,40], baseColor:[235,146,31],
      species:[{acid:"HIn⁺ (rosso)", base:"In (giallo-arancio)"}]},
    {key:"rossometile",   label:"Indicatore: rosso di metile", short:"Rosso metile", pKa:5.10, C:0.00005, color:"var(--c-5)", cat:"indicator", group:"Indicatori",
      acidColor:[200,40,40], baseColor:[230,195,45],
      species:[{acid:"HIn (rosso)", base:"In⁻ (giallo)"}]},

    {key:"hcl",  label:"HCl (acido forte, pKa −7 stimato)", short:"HCl", pKa:-7, C:1.0, color:"var(--c-strong)", cat:"strong", group:"Acidi e basi forti",
      species:[{acid:"HCl", base:"Cl⁻"}]},
    {key:"naoh", label:"NaOH (base, pKa 14: equilibrio H₂O/OH⁻)", short:"NaOH", pKa:14, C:1.0, color:"var(--c-strong)", cat:"strong", group:"Acidi e basi forti",
      species:[{acid:"H₂O", base:"OH⁻"}]}
  ];

  let vessels = []; // sempre almeno un vaso (di norma l'acqua)
  let valves = [];  // length vessels.length-1 ; true = aperta (collegata)
  let stacked = []; // length vessels.length-1 ; true = i due vasi condividono la stessa colonna visiva
  let nextId = 1;

  // --- titolante: quantità assoluta per vaso (v.tauExtra), NON più un unico numero
  // globale — vedi groupTarget più sopra. v.tauBase è solo un offset di VISUALIZZAZIONE
  // (non entra in nessun calcolo chimico): permette allo slider di un gruppo di "ripartire
  // da 0" dopo una preparazione o una separazione, pur mantenendo lo stesso pH reale.
  function groupDisplayTau(members){
    return members.reduce((s,v)=> s + (v.tauExtra||0) - (v.tauBase||0), 0);
  }
  function setGroupDisplayTau(members, mmolValue){
    const delta = mmolValue/1000 - groupDisplayTau(members);
    members[0].tauExtra = (members[0].tauExtra||0) + delta;
  }
  // diluizione dell'INTERO gruppo (v.dil, ≥1 su ogni vaso non-acqua del gruppo): non
  // è il titolante a cambiare slider, è la stessa quantità già aggiunta che — dissolta
  // in sempre più acqua — pesa sempre meno, mentre la chimica di ciascun vaso (in
  // vesselN/vesselBeta) si mescola sempre più con l'acqua pura. Combinate, le due cose
  // fanno sì che il pH tenda a 7 quando dil→∞, qualunque titolante fosse stato lasciato.
  function groupDil(members){
    const v = members.find(v=>!v.water);
    return v ? (v.dil||1) : 1;
  }
  function setGroupDil(members, newDil){
    members.forEach(v=>{
      if(v.water) return;
      const oldDil = v.dil||1;
      const ratio = oldDil/newDil; // f_new/f_old
      v.tauExtra = (v.tauExtra||0) * ratio;
      v.tauBase = (v.tauBase||0) * ratio;
      v.dil = newDil;
    });
  }
  // estremo (mol/L, simmetrico) dello slider di UN gruppo: quanto titolante serve
  // davvero a portarlo ai bordi del proprio dominio (0/14, esteso se ha un pKa fuori
  // scala), con un margine — così quei bordi restano sempre raggiungibili.
  function groupBound(members){
    // un gruppo CON acqua resta sempre mirato a 0/14: a pH così estremi (21, 28...) il
    // contributo di [H+]/[OH-] dell'acqua da solo esploderebbe la barra a valori assurdi
    // (all'estremo NaOH, [OH-] varrebbe ~10¹⁴ mol/L) senza bisogno reale di arrivarci —
    // solo un gruppo isolato SENZA acqua può davvero mirare al proprio pKa fuori scala.
    const hasWater = members.some(v=>v.water);
    const {lo,hi} = hasWater ? {lo:0,hi:14} : groupDomain(members);
    const baseSum = members.reduce((s,v)=>s+(v.tauBase||0),0);
    const rawLo = tauForLevel(members, lo)*1.15 - baseSum;
    const rawHi = tauForLevel(members, hi)*1.15 - baseSum;
    return Math.max(1, Math.abs(rawLo), Math.abs(rawHi)); // mol/L, minimo 1 di sicurezza
  }
  // consolida un (sotto)gruppo appena isolato su un unico vaso rappresentante: gli altri
  // restano a contributo zero (conta solo la somma), il rappresentante assorbe la
  // quantità assoluta che riproduce ESATTAMENTE "level" da solo (continuità del pH), e il
  // suo tauBase è impostato uguale, così il display riparte da 0 proprio da qui.
  function consolidateGroup(members, level){
    const rawNeeded = tauForLevel(members, level);
    members.forEach(v=>{ v.tauExtra = 0; v.tauBase = 0; });
    members[0].tauExtra = rawNeeded;
    members[0].tauBase = rawNeeded;
  }
  // unica porta d'accesso per aprire/chiudere una valvola (di rotaia o interna a una
  // pila): chiudendola, il gruppo attuale si spacca in due e ciascuna metà viene
  // consolidata (vedi sopra) così il pH di ognuna resta esattamente quello di un attimo
  // prima — nessun salto, e da quel momento le due parti si titolano indipendentemente.
  function toggleValve(i){
    if(valves[i]){
      const groups = groupsFromValves();
      const group = groups.find(g => g.includes(vessels[i]) && g.includes(vessels[i+1]));
      const level = solveLevel(group);
      valves[i] = false;
      const newGroups = groupsFromValves();
      const subA = newGroups.find(g=>g.includes(vessels[i]));
      const subB = newGroups.find(g=>g.includes(vessels[i+1]));
      consolidateGroup(subA, level);
      if(subB!==subA) consolidateGroup(subB, level);
    } else {
      valves[i] = true;
    }
    renderAll();
  }

  function addVesselFromPreset(key){
    const p = PRESETS.find(p=>p.key===key) || PRESETS[0];
    const v = { id:"v"+(nextId++), label:p.label, short:p.short, color:p.color, tauExtra:0, tauBase:0, dil:1 };
    v.pKa = Array.isArray(p.pKa) ? p.pKa.slice() : p.pKa;
    v.C = p.C;
    v.species = p.species;
    if(p.acidColor && p.baseColor){ v.acidColor = p.acidColor; v.baseColor = p.baseColor; }
    vessels.push(v);
    valves.push(true);
    stacked.push(false);
  }
  function resetDefault(){
    vessels = [{ id:"water", water:true, label:"Acqua pura", short:"Acqua", color:"var(--c-water)", tauExtra:0, tauBase:0 }];
    valves = [];
    stacked = [];
    nextId = 1;
    mode = "libero";
  }
  let mode = "libero";
  resetDefault();

  function groupsFromValves(){
    const groups = [];
    let cur = [vessels[0]];
    for(let i=0;i<valves.length;i++){
      if(valves[i]){
        cur.push(vessels[i+1]);
      } else {
        groups.push(cur);
        cur = [vessels[i+1]];
      }
    }
    groups.push(cur);
    return groups;
  }

  // colonne: vasi consecutivi impilati (stacked[i]===true) condividono la stessa
  // colonna visiva. Ogni colonna è un array di INDICI in vessels[] (non oggetti),
  // usato sia dal banco (renderRig) sia dal pannello "Vasi sul banco".
  function computeColumns(){
    const columns = [];
    let cur = [0];
    for(let i=0;i<valves.length;i++){
      if(stacked[i]) cur.push(i+1);
      else { columns.push(cur); cur = [i+1]; }
    }
    columns.push(cur);
    return columns;
  }

  function computeLevels(){
    const groups = groupsFromValves();
    const levelByVesselId = {};
    const groupInfo = [];
    for(const g of groups){
      const grounded = isGrounded(g);
      const level = grounded ? solveLevel(g) : null;
      const betaTotal = grounded ? g.reduce((s,v)=>s+vesselBeta(v, level), 0) : 0;
      groupInfo.push({members:g, level, betaTotal, grounded});
      for(const v of g) levelByVesselId[v.id] = level;
    }
    return {levelByVesselId, groupInfo};
  }

  /* ============================================================
     COLORE DEL LIQUIDO — scala tipo indicatore universale
     (funzione del VALORE del pH, non della posizione sul banco)
     ============================================================ */
  const PH_STOPS = [
    [0,  [176,36,43]],
    [2,  [199,56,38]],
    [4,  [222,120,30]],
    [6,  [190,157,35]],
    [7,  [70,148,80]],
    [8,  [31,148,150]],
    [10, [45,105,183]],
    [12, [70,72,168]],
    [14, [118,58,150]]
  ];
  function pHColor(pH){
    pH = Math.max(0, Math.min(14, pH));
    for(let i=0;i<PH_STOPS.length-1;i++){
      const [p0,c0] = PH_STOPS[i], [p1,c1] = PH_STOPS[i+1];
      if(pH>=p0 && pH<=p1){
        const t = (pH-p0)/(p1-p0);
        const r = Math.round(c0[0]+(c1[0]-c0[0])*t);
        const g = Math.round(c0[1]+(c1[1]-c0[1])*t);
        const b = Math.round(c0[2]+(c1[2]-c0[2])*t);
        return `rgb(${r},${g},${b})`;
      }
    }
    return "rgb(120,120,120)";
  }
  // per gli indicatori acido-base, il colore "reale" delle due forme (quello che
  // si vede davvero in laboratorio) conta più della scala universale: si
  // interpola fra colore della forma acida e della base coniugata in base alla
  // stessa frazione di dissociazione già usata per la chimica (pairN/C).
  function indicatorColor(v, pH){
    const Ka = Math.pow(10, -v.pKa);
    const frac = pairN(Ka, v.C, pH) / v.C; // 0 = tutto HIn, 1 = tutto In-
    const [r0,g0,b0] = v.acidColor, [r1,g1,b1] = v.baseColor;
    const r = Math.round(r0+(r1-r0)*frac);
    const g = Math.round(g0+(g1-g0)*frac);
    const b = Math.round(b0+(b1-b0)*frac);
    return `rgb(${r},${g},${b})`;
  }

  /* ============================================================
     RENDER: BANCO DI LAVORO (vasi)
     Asse pH invertito: pH 0 in alto, pH 14 in basso — allineato al
     potenziale chimico del protone, che cresce verso l'alto.
     ============================================================ */
  const svgNS = "http://www.w3.org/2000/svg";
  const rigSvg = document.getElementById("rigSvg");

  // scala FISSA (non dipende dal numero o dalla composizione dei vasi correnti):
  // calibrata su un tampone 1 M esattamente al proprio pKa.
  const REF_MAX_BETA = Math.LN10 * 1.0 / 4;
  const FIXED_HALF_PX = 26;
  const MIN_HALF_PX = 2;
  // il taglio (clamp) vero e proprio dei bordi del vaso avviene molto più in là di
  // FIXED_HALF_PX: si allarga fin quasi a toccare il vaso vicino (CLIP_GAP_PX di
  // margine), entro un tetto assoluto ragionevole. La CALIBRAZIONE del disegno (quanti
  // pixel vale 1 M esattamente al proprio pKa) resta invece ancorata a FIXED_HALF_PX,
  // quindi l'aspetto dei vasi "normali" non cambia: solo le concentrazioni molto più
  // alte del riferimento continuano a crescere prima di venire davvero tagliate.
  const CLIP_GAP_PX = 6;
  const ABS_MAX_HALF_PX = 75;

  const RIG = { x0:112, x1:940, yTop:120, yBottom:470, railY:492, stubTop:58, stubBot:572, viewW:980, viewH:610 };
  const RIG_MAX_COLS = 5; // al massimo 5 vasi/colonne affiancati senza scorrimento orizzontale
  // HCl/NaOH sono vasi reali con pKa fuori da 0-14 (-7 e 21): quando sono sul banco si
  // guadagna spazio sopra/sotto per mostrarne il vero rigonfiamento, con la scala che
  // dallo 0 (o dal 14) salta direttamente a quel valore, invece di restare vuota nel
  // mezzo. EXTRA_LO/EXTRA_HI sono i domini virtuali di quello spazio (il vero pKa
  // cade esattamente a metà), EXTRA_VESSEL_H i pixel guadagnati.
  const EXTRA_VESSEL_H = 130;
  const EXTRA_LO = -14;
  const EXTRA_HI = 28;

  function sampleH(fn, steps, lo, hi){
    if(lo===undefined){ lo=0; hi=14; }
    const out = [];
    for(let i=0;i<=steps;i++){
      const pH = lo + (hi-lo)*i/steps;
      out.push([pH, fn(pH)]);
    }
    return out;
  }

  function el(tag, attrs, parent){
    const e = document.createElementNS(svgNS, tag);
    for(const k in attrs) e.setAttribute(k, attrs[k]);
    if(parent) parent.appendChild(e);
    return e;
  }

  function resolveVar(name){
    return getComputedStyle(document.documentElement).getPropertyValue(name.replace(/^var\(|\)$/g,"")).trim();
  }
  function colorOf(v){
    if(v.color && v.color.startsWith("var(")) return resolveVar(v.color);
    return v.color;
  }

  // pH0 -> yTop, pH14 -> yBottom; oltre 0 o 14 continua nello spazio guadagnato sopra/sotto
  // (vedi EXTRA_LO/EXTRA_HI) — un'unica funzione continua, così un vaso con pKa fuori
  // scala (HCl, NaOH) si disegna come un vaso qualunque, solo su un tratto d'asse più lungo.
  const yOf = pH => {
    if(pH >= 0 && pH <= 14) return RIG.yTop + (pH/14)*(RIG.yBottom-RIG.yTop);
    if(pH < 0) return -EXTRA_VESSEL_H + ((pH-EXTRA_LO)/(0-EXTRA_LO)) * EXTRA_VESSEL_H;
    return RIG.viewH + ((pH-14)/(EXTRA_HI-14)) * EXTRA_VESSEL_H;
  };

  function renderRig(){
    rigSvg.innerHTML = "";
    const defs = el("defs", {}, rigSvg);

    const columns = computeColumns();
    const nCols = columns.length;
    // scala per colonna: fissa, calibrata su al massimo 5 colonne affiancate. Con più di 5
    // colonne le successive NON si rimpiccioliscono: il banco si allunga e scorre in orizzontale.
    const slotW = (RIG.x1-RIG.x0)/Math.min(nCols, RIG_MAX_COLS);
    const minHalf = MIN_HALF_PX;
    // calibrazione: un tampone 1 M esattamente al proprio pKa arriva qui (FISSA, non
    // dipende da quante colonne ci sono davvero sul banco)
    const calHalf = Math.min(FIXED_HALF_PX, slotW*0.38);
    const pxPerBeta = (calHalf-minHalf)/REF_MAX_BETA;
    // tetto di disegno vero e proprio: molto più largo della calibrazione, si ferma
    // solo a un pelo (CLIP_GAP_PX) dalla colonna vicina, o al tetto assoluto se il
    // banco ha pochissime colonne
    const maxHalf = Math.min(ABS_MAX_HALF_PX, Math.max(calHalf, slotW/2 - CLIP_GAP_PX));
    const plotX1 = RIG.x0 + slotW*nCols;
    const colCenters = columns.map((_,ci)=> RIG.x0 + slotW*(ci+0.5));
    const vesselColIdx = {};
    columns.forEach((idxs,ci)=> idxs.forEach(i=> vesselColIdx[i]=ci));

    // se un vaso ha un pKa sotto 0 (HCl) e/o sopra 14 (NaOH), si guadagna spazio SOPRA
    // e/o SOTTO per mostrarne il vero rigonfiamento: per tutti gli altri casi il canvas
    // resta esattamente quello di sempre.
    const allPkas = vessels.flatMap(v=> v.water ? [] : (Array.isArray(v.pKa)?v.pKa:[v.pKa]));
    const hasLowPka = allPkas.some(pk=>pk<0);
    const hasHighPka = allPkas.some(pk=>pk>14);
    const topExtra = hasLowPka ? EXTRA_VESSEL_H : 0;
    const botExtra = hasHighPka ? EXTRA_VESSEL_H : 0;
    const viewBoxY = -topExtra;
    const viewBoxH = RIG.viewH + topExtra + botExtra;

    // larghezza dell'SVG: al massimo 5 colonne riempiono il contenitore come prima; oltre,
    // la larghezza in pixel resta quella "per colonna" già in uso (misurata al contenitore),
    // e il banco scorre in orizzontale invece di rimpicciolire i vasi già presenti.
    const rigWrap = rigSvg.parentElement;
    if(nCols <= RIG_MAX_COLS){
      rigSvg.setAttribute("viewBox", `0 ${viewBoxY} ${RIG.viewW} ${viewBoxH}`);
      rigSvg.style.width = "100%";
      rigWrap.style.overflowX = "visible";
    } else {
      rigSvg.style.width = "100%"; // misura la larghezza naturale del contenitore prima di fissarla
      const containerW = rigSvg.getBoundingClientRect().width || RIG.viewW;
      const pxPerUnit = containerW / RIG.viewW;
      const totalUnits = plotX1 + (RIG.viewW - RIG.x1); // stesso margine destro della vista a 5 colonne
      rigSvg.setAttribute("viewBox", `0 ${viewBoxY} ${totalUnits} ${viewBoxH}`);
      rigSvg.style.width = Math.round(pxPerUnit*totalUnits) + "px";
      rigWrap.style.overflowX = "auto";
    }

    const {levelByVesselId, groupInfo} = computeLevels();

    /* --- assi: pH (scende) e potenziale chimico (sale) --- */
    const axisG = el("g", {}, rigSvg);
    for(let p=0;p<=14;p+=2){
      const y = yOf(p);
      el("line", {x1:RIG.x0-10,x2:plotX1+8,y1:y,y2:y,stroke:"var(--line)","stroke-width":1,"shape-rendering":"crispEdges"}, axisG);
      const t = el("text", {x:RIG.x0-18,y:y+4,"text-anchor":"end","font-family":"IBM Plex Mono, monospace","font-size":13,fill:"var(--ink-faint)"}, axisG);
      t.textContent = p;
      const mu = (-RT_LN10_KJ*p).toFixed(1);
      const tm = el("text", {x:RIG.x0-58,y:y+4,"text-anchor":"end","font-family":"IBM Plex Mono, monospace","font-size":12,fill:"var(--ink-faint)",opacity:0.85}, axisG);
      tm.textContent = mu;
    }
    // scala interrotta: dallo 0 si salta direttamente a -7 (HCl) e/o dal 14 a 21 (NaOH),
    // con un simbolo a zig-zag sull'asse a segnare il salto, invece di riempire lo spazio
    // vuoto fra -1 e -6 (o 15 e 20) con una scala enorme e inutile
    [[hasLowPka,-7],[hasHighPka,21]].forEach(([has,p])=>{
      if(!has) return;
      const y = yOf(p);
      const yBreak = yOf(p<0 ? 0 : 14);
      el("line", {x1:RIG.x0-10,x2:plotX1+8,y1:y,y2:y,stroke:"var(--line)","stroke-width":1,"shape-rendering":"crispEdges"}, axisG);
      el("text", {x:RIG.x0-18,y:y+4,"text-anchor":"end","font-family":"IBM Plex Mono, monospace","font-size":13,fill:"var(--ink-faint)"}, axisG).textContent = p;
      const midY = (y+yBreak)/2;
      [-5,5].forEach(dx=>{
        el("line", {x1:RIG.x0-14+dx, x2:RIG.x0-4+dx, y1:midY-4, y2:midY+4, stroke:"var(--ink-faint)", "stroke-width":1.4}, axisG);
      });
    });
    // etichette con freccia direzionale
    el("text", {x:RIG.x0-18, y:RIG.yTop-38, "text-anchor":"end","font-family":"IBM Plex Sans, sans-serif","font-size":14,"font-weight":600,fill:"var(--ink-soft)"}, axisG).textContent = "pH ↓";
    el("text", {x:RIG.x0-58, y:RIG.yTop-38, "text-anchor":"end","font-family":"IBM Plex Sans, sans-serif","font-size":13,"font-weight":600,fill:"var(--ink-soft)"}, axisG).textContent = "μ(H⁺) ↑";
    el("text", {x:RIG.x0-58, y:RIG.yTop-22, "text-anchor":"end","font-family":"IBM Plex Sans, sans-serif","font-size":11,fill:"var(--ink-faint)"}, axisG).textContent = "−μ° in kJ/mol";

    /* --- rotaia sottile --- */
    el("line", {x1:RIG.x0, x2:plotX1, y1:RIG.railY, y2:RIG.railY, stroke:"var(--line-strong)","stroke-width":3,"stroke-linecap":"round"}, rigSvg);

    /* --- valvole / giunti impilati --- */
    for(let i=0;i<valves.length;i++){
      const cxA = colCenters[vesselColIdx[i]];
      const cxB = colCenters[vesselColIdx[i+1]];
      const open = valves[i];
      const g = el("g", {class:"valve-btn", "data-valve":i, style:"cursor:pointer;"}, rigSvg);
      if(stacked[i]){
        el("circle", {cx:cxA, cy:RIG.railY, r:9, fill:"var(--accent-soft)", stroke:"var(--accent)", "stroke-width":2}, g);
        el("text", {x:cxA, y:RIG.railY+4, "text-anchor":"middle","font-family":"IBM Plex Sans, sans-serif","font-size":12.5,"font-weight":700, fill:"var(--accent-strong)"}, g).textContent = "⋮";
        el("title", {}, g).textContent = "Vasi impilati nella stessa colonna — clic per separare";
      } else {
        const x = (cxA+cxB)/2;
        el("line", {x1:cxA, x2:cxB, y1:RIG.railY, y2:RIG.railY, stroke: open ? "var(--accent)" : "var(--warn)", "stroke-width":3, "stroke-linecap":"round", opacity: open?1:0.6}, rigSvg);
        el("circle", {cx:x, cy:RIG.railY, r:10, fill:"var(--surface)", class:"valve-ring", stroke: open ? "var(--accent)" : "var(--warn)", "stroke-width":2}, g);
        if(open){
          el("circle", {cx:x, cy:RIG.railY, r:3.2, fill:"var(--accent)"}, g);
        } else {
          el("line", {x1:x-4,x2:x+4,y1:RIG.railY-4,y2:RIG.railY+4, stroke:"var(--warn)", "stroke-width":2, "stroke-linecap":"round"}, g);
          el("line", {x1:x-4,x2:x+4,y1:RIG.railY+4,y2:RIG.railY-4, stroke:"var(--warn)", "stroke-width":2, "stroke-linecap":"round"}, g);
        }
        el("title", {}, g).textContent = open ? "Valvola aperta — clic per chiudere" : "Valvola chiusa — clic per aprire";
      }
    }

    /* --- linea di equilibrio tratteggiata rossa, per ciascun gruppo con pH definito ---
       disegnata PRIMA dei vasi e dei rubinetti interni, così questi ultimi restano sempre
       sopra (e quindi cliccabili) anche quando cadono quasi alla stessa altezza --- */
    groupInfo.forEach(g=>{
      if(!g.grounded) return;
      const gIdxs = g.members.map(v=>vessels.indexOf(v));
      const colIdxs = gIdxs.map(i=>vesselColIdx[i]);
      const xL = colCenters[Math.min(...colIdxs)] - maxHalf - 14;
      const xR = colCenters[Math.max(...colIdxs)] + maxHalf + 14;
      const y = yOf(g.level);
      el("line", {x1:xL, x2:xR, y1:y, y2:y, stroke:"var(--eq-line)", "stroke-width":1.6, "stroke-dasharray":"6,4"}, rigSvg);
      el("text", {x:xR+8, y:y+4.5, "font-family":"IBM Plex Mono, monospace","font-size":14.5, "font-weight":600, fill:"var(--eq-line)"}, rigSvg).textContent = `pH ${g.level.toFixed(2)}`;
    });

    columns.forEach((idxs, ci)=>{
      const cx = colCenters[ci];
      const members = idxs.map(i=>vessels[i]);
      const widthFn = pH => members.reduce((s,vv)=> s + vesselBeta(vv,pH), 0);
      // dominio di questa colonna: normalmente 0-14, esteso se uno dei vasi impilati ha
      // un pKa fuori scala (HCl/NaOH) — un'unica forma continua, non un pezzo a parte
      const colPkas = members.flatMap(v=> v.water?[]:(Array.isArray(v.pKa)?v.pKa:[v.pKa]));
      const domainLo = colPkas.some(pk=>pk<0) ? EXTRA_LO : 0;
      const domainHi = colPkas.some(pk=>pk>14) ? EXTRA_HI : 14;
      const steps = Math.round(64 * (domainHi-domainLo)/14);
      const pts = sampleH(widthFn, steps, domainLo, domainHi);
      const left = [], right = [];
      for(const [pH,b] of pts){
        const half = Math.min(maxHalf, Math.max(minHalf, b*pxPerBeta));
        const y = yOf(pH);
        left.push([cx-half, y]);
        right.push([cx+half, y]);
      }
      const pathD = "M " + left.map(p=>p[0]+","+p[1]).join(" L ") +
                    " L " + right.slice().reverse().map(p=>p[0]+","+p[1]).join(" L ") + " Z";

      const clipId = "clip-col"+ci;
      const clip = el("clipPath", {id:clipId}, defs);
      el("path", {d:pathD}, clip);

      const gradId = "grad-col"+ci;
      const grad = el("linearGradient", {id:gradId, x1:"0%", y1:"0%", x2:"100%", y2:"0%"}, defs);
      el("stop", {offset:"0%", "stop-color":"#ffffff", "stop-opacity":0.35}, grad);
      el("stop", {offset:"18%", "stop-color":"#ffffff", "stop-opacity":0.05}, grad);
      el("stop", {offset:"85%", "stop-color":"#000000", "stop-opacity":0.03}, grad);
      el("stop", {offset:"100%", "stop-color":"#000000", "stop-opacity":0.12}, grad);

      const vg = el("g", {}, rigSvg);

      // un vaso impilato può contenere PIÙ gruppi indipendenti (valvole interne chiuse):
      // ognuno riempie SOLO il proprio tratto di colonna, dal proprio livello fino al
      // confine col vicino — esattamente come vasi separati, solo disegnati insieme.
      // Il confine fra un membro e il successivo è il punto medio fra i loro pKa (l'acqua
      // vale 7, un punto neutro di riferimento) cosicché ognuno "possiede" la porzione di
      // vaso più vicina al proprio pKa.
      const anchors = members.map(v=> v.water ? 7 : primaryPka(v));
      const boundaries = [domainLo];
      for(let k=0;k<members.length-1;k++) boundaries.push((anchors[k]+anchors[k+1])/2);
      boundaries.push(domainHi);

      let subStart = 0;
      for(let k=0;k<=members.length-1;k++){
        const isLastMember = k===members.length-1;
        const boundaryOpen = isLastMember ? false : valves[idxs[k]];
        if(isLastMember || !boundaryOpen){
          const subEnd = k;
          const subV = members[subStart];
          const level = levelByVesselId[subV.id];
          if(level !== null){
            const yTopSub = yOf(boundaries[subStart]);
            const yBotSub = yOf(boundaries[subEnd+1]);
            const levelY = Math.max(yOf(level), yTopSub);
            const liquidColor = (subV.acidColor && subV.baseColor) ? indicatorColor(subV, level) : pHColor(level);
            const fillH = Math.max(0, yBotSub-levelY);
            if(fillH>0){
              el("rect", {x:cx-maxHalf-4, y:levelY, width:(maxHalf+4)*2, height:fillH, fill:liquidColor, "clip-path":`url(#${clipId})`, opacity:0.86}, vg);
              el("rect", {x:cx-maxHalf-4, y:levelY-1.5, width:(maxHalf+4)*2, height:3, fill:liquidColor, "clip-path":`url(#${clipId})`, opacity:1}, vg);
            }
          }
          subStart = k+1;
        }
      }

      const outlineColor = colorOf(members[0]);
      el("path", {d:pathD, fill:`url(#${gradId})`, stroke: outlineColor, "stroke-width":2, "clip-path":`url(#${clipId})`}, vg);
      el("path", {d:pathD, fill:"none", stroke: outlineColor, "stroke-width":2}, vg);

      // rubinetti interni: uno per ogni confine fra vasi impilati, indipendenti dal
      // rubinetto/giunto sulla rotaia in basso — permettono di isolare un tratto (col suo
      // proprio livello, come un vaso a sé) o ricollegarlo al resto della colonna
      for(let k=0;k<members.length-1;k++){
        const gi = idxs[k];
        const open = valves[gi];
        const y = yOf(boundaries[k+1]);
        const ig = el("g", {class:"stack-valve-btn", "data-valve":gi, style:"cursor:pointer;"}, rigSvg);
        el("circle", {cx, cy:y, r:12.5, fill:"var(--surface)", opacity:0.97, stroke: open ? "var(--accent)" : "var(--warn)", "stroke-width":2.6}, ig);
        if(open){ el("circle", {cx, cy:y, r:4, fill:"var(--accent)"}, ig); }
        else {
          el("line", {x1:cx-4.5,x2:cx+4.5,y1:y-4.5,y2:y+4.5, stroke:"var(--warn)", "stroke-width":2.2, "stroke-linecap":"round"}, ig);
          el("line", {x1:cx-4.5,x2:cx+4.5,y1:y+4.5,y2:y-4.5, stroke:"var(--warn)", "stroke-width":2.2, "stroke-linecap":"round"}, ig);
        }
        el("title", {}, ig).textContent = open ? "Rubinetto interno aperto — clic per chiudere" : "Rubinetto interno chiuso — clic per aprire";
      }

      members.forEach(v=>{
        // pKa tick: ben visibile, ancorato alla parete reale della colonna in quel punto
        // (funziona anche per pKa fuori scala come HCl/NaOH, che ora sono vasi come gli altri)
        if(!v.water){
          const pKas = Array.isArray(v.pKa) ? v.pKa : [v.pKa];
          pKas.forEach(pk=>{
            const ypk = yOf(pk);
            const halfAtPka = Math.min(maxHalf, Math.max(minHalf, widthFn(pk)*pxPerBeta));
            el("line", {x1:cx-halfAtPka-12, x2:cx-halfAtPka-2, y1:ypk, y2:ypk, stroke:colorOf(v), "stroke-width":3, "stroke-linecap":"round"}, vg);
          });
        }
        // modalità "Specie predominanti": a destra della colonna, la specie che
        // predomina in ciascuna zona di pH separata dai pKa (acida più in alto, via
        // via più basica scendendo, come i pKa stessi lungo l'asse)
        if(mode==="specie" && !v.water && v.species){
          const pKas = Array.isArray(v.pKa) ? v.pKa : [v.pKa];
          const bounds = [domainLo, ...pKas, domainHi];
          const zoneLabels = [v.species[0].acid];
          v.species.forEach(s=> zoneLabels.push(s.base));
          zoneLabels.forEach((label, zi)=>{
            const loB = bounds[zi], hiB = bounds[zi+1];
            if(hiB-loB < 0.15) return;
            const mid = (loB+hiB)/2;
            const ymid = yOf(mid);
            const halfAtMid = Math.min(maxHalf, Math.max(minHalf, widthFn(mid)*pxPerBeta));
            el("text", {x:cx+halfAtMid+12, y:ymid+5, "text-anchor":"start","font-family":"IBM Plex Sans, sans-serif","font-size":13.5,"font-weight":600,"font-style":"italic", fill:colorOf(v)}, rigSvg).textContent = label;
          });
        }
      });

      // etichetta sotto la colonna: nomi brevi (per non sovrapporsi ai vicini) + dettaglio
      const labelFont = Math.max(12, Math.min(16, slotW/7));
      const subFont = Math.max(10.5, Math.min(13.5, slotW/8.2));
      // staccato dal tubo con le valvole (railY) di più del semplice bordo del vaso
      // (yBottom), così il nome non resta "appiccicato" alla rotaia
      const labelY = RIG.railY + 26;
      const lbl = el("text", {x:cx, y:labelY, "text-anchor":"middle","font-family":"IBM Plex Sans, sans-serif","font-size":labelFont, "font-weight":600, fill:"var(--ink)"}, rigSvg);
      lbl.textContent = members.map(v=>v.short||v.label).join(" + ");
      let sub;
      if(members.length===1){
        const v = members[0];
        if(v.water) sub = "senza tampone";
        else { const pKas = Array.isArray(v.pKa)?v.pKa:[v.pKa]; sub = "pKa " + pKas.map(p=>p.toFixed(2)).join("/"); }
      } else {
        const withPka = members.filter(v=>!v.water);
        const pkaPart = withPka.map(v=> (Array.isArray(v.pKa)?v.pKa:[v.pKa]).map(p=>p.toFixed(2)).join("/")).join(" · ");
        sub = members.some(v=>v.water) ? ("pKa " + pkaPart + " · acqua") : ("pKa " + pkaPart);
      }
      el("text", {x:cx, y:labelY+15, "text-anchor":"middle","font-family":"IBM Plex Mono, monospace","font-size":subFont, fill:"var(--ink-faint)"}, rigSvg).textContent = sub;
    });

    // click sulle valvole / giunti impilati (ciclo: impilato -> collegato -> chiuso -> collegato)
    rigSvg.querySelectorAll(".valve-btn").forEach(g=>{
      g.addEventListener("click", ()=>{
        const i = +g.getAttribute("data-valve");
        // impilato: clic separa in due colonne (la valvola resta quella che era, si vede
        // e si comanda separatamente col rubinetto sulla rotaia); non impilato: clic
        // apre/chiude direttamente la valvola fra le due colonne affiancate
        if(stacked[i]){ stacked[i] = false; renderAll(); }
        else toggleValve(i);
      });
    });
    rigSvg.querySelectorAll(".stack-valve-btn").forEach(g=>{
      g.addEventListener("click", ()=>{
        const i = +g.getAttribute("data-valve");
        toggleValve(i);
      });
    });

    return {levelByVesselId, groupInfo};
  }

  /* ============================================================
     RENDER: STATISTICHE
     ============================================================ */
  function renderStats(groupInfo){
    const row = document.getElementById("statRow");
    row.innerHTML = "";
    groupInfo.forEach((g, gi)=>{
      const names = g.members.map(v=>v.short||v.label).join(" + ");
      const rep = g.members.find(v=>!v.water && !v.strong) || g.members[g.members.length-1];
      const swatchColor = colorOf(rep);
      const span = document.createElement("span");
      span.className = "gsum";
      span.innerHTML = g.grounded
        ? `<span class="swatch" style="background:${swatchColor}"></span><span class="muted">${names}:</span> pH <b>${g.level.toFixed(2)}</b> <span class="muted">·</span> &beta; <b>${(g.betaTotal*1000).toFixed(2)}</b><span class="muted"> mmol/(L&middot;pH)</span>`
        : `<span class="swatch" style="background:${swatchColor}"></span><span class="muted">${names}:</span> <span class="muted">non in soluzione</span>`;
      row.appendChild(span);
    });
  }

  /* ============================================================
     RENDER: PANNELLO VASI (controlli)
     ============================================================ */
  const vesselList = document.getElementById("vesselList");
  const GRAB = `<span class="vgrab" title="Trascinate SOLO da qui per riordinare o impilare">⠿</span>`;

  function buildVesselCard(v){
    const card = document.createElement("div");
    card.className = "vcard";
    card.dataset.vesselId = v.id;
    if(v.water){
      card.innerHTML = `
        <div class="vhead">
          ${GRAB}
          <span class="vdot" style="background:${colorOf(v)}"></span>
          <span class="vname" style="cursor:default;">${v.label}</span>
          <button class="vdel" draggable="false" data-id="${v.id}" title="Rimuovi vaso">&times;</button>
        </div>
        <div class="fixed-note">Riferimento: solo [H⁺]/[OH⁻], nessun tampone.</div>
      `;
      wireDragForCard(card, v);
      return card;
    }
    const pKas = Array.isArray(v.pKa) ? v.pKa : [v.pKa];
    card.innerHTML = `
      <div class="vhead">
        ${GRAB}
        <span class="vdot" style="background:${colorOf(v)}"></span>
        <input class="vname" draggable="false" value="${v.label}" data-id="${v.id}" data-field="label">
        <button class="vdel" draggable="false" data-id="${v.id}" title="Rimuovi vaso">&times;</button>
      </div>
      <div class="pka-readout">pK<sub>a</sub> <span class="mono">${pKas.map(p=>p.toFixed(2)).join(" / ")}</span></div>
      <div class="row-slider">
        <span>Conc.</span>
        <input type="range" draggable="false" min="0" max="100" step="1" value="${cToSlider(v.C)}" data-id="${v.id}" data-field="C">
        <output>${fmtConc(v.C)}</output>
      </div>
    `;
    wireDragForCard(card, v);
    return card;
  }

  /* Le schede sono raggruppate per colonna (vasi impilati = stessa colonna, vedi
     computeColumns): un contenitore "a mazzo di carte" per i gruppi impilati, e un
     connettore etichettato (collegati/separati) fra una colonna e la successiva.
     Sia il connettore di impilamento sia quello di collegamento sono cliccabili,
     esattamente come le valvole/i giunti sul banco. */
  function renderVesselControls(){
    vesselList.innerHTML = "";
    const columns = computeColumns();
    columns.forEach((idxs, ci)=>{
      const group = document.createElement("div");
      group.className = "vstack-group" + (idxs.length>1 ? " multi" : "");
      idxs.forEach((vi, k)=>{
        group.appendChild(buildVesselCard(vessels[vi]));
        if(k < idxs.length-1){
          const open = valves[vi];
          const tag = document.createElement("div");
          tag.className = "stack-tag";
          tag.innerHTML = `<span class="stack-valve-mini ${open?"open":""}" data-idx="${vi}" title="${open?"Rubinetto interno aperto — clic per chiudere":"Rubinetto interno chiuso — clic per aprire"}">${open?"●":"✕"}</span><span class="stack-sep" data-idx="${vi}" title="Vasi impilati nella stessa colonna — clic per separare"><span>⋮</span> stessa colonna</span>`;
          group.appendChild(tag);
        }
      });
      vesselList.appendChild(group);
      if(ci < columns.length-1){
        const boundaryIdx = idxs[idxs.length-1];
        const open = valves[boundaryIdx];
        const conn = document.createElement("div");
        conn.className = "col-conn" + (open ? " open" : " closed");
        conn.dataset.idx = boundaryIdx;
        conn.title = open ? "Valvola aperta, stesso pH — clic per chiudere" : "Valvola chiusa, pH indipendenti — clic per aprire";
        conn.innerHTML = `<span class="stack-valve-mini conn-dot ${open?"open":""}">${open?"●":"✕"}</span> ${open?"collegati":"separati"}`;
        vesselList.appendChild(conn);
      }
    });

    vesselList.querySelectorAll('input[data-field="label"]').forEach(inp=>{
      inp.addEventListener("input", e=>{
        const v = vessels.find(v=>v.id===e.target.dataset.id);
        v.label = e.target.value || v.label;
        renderRigAndStats();
      });
    });
    vesselList.querySelectorAll('input[data-field="C"]').forEach(inp=>{
      inp.addEventListener("input", e=>{
        const v = vessels.find(v=>v.id===e.target.dataset.id);
        if(!v) return;
        // diluire/concentrare non tocca solo la capacità tampone (v.C): anche il
        // titolante già aggiunto (v.tauExtra, v.tauBase) è disciolto nella stessa
        // soluzione e va perciò riscalato dallo stesso fattore — altrimenti
        // resterebbe una quantità assoluta fissa che, relativamente a un v.C
        // sempre più piccolo, farebbe divergere il pH invece di lasciarlo
        // costante (comportamento tipico di un tampone diluito).
        const newC = sliderToC(+e.target.value);
        if(v.C > 0){
          const ratio = newC / v.C;
          v.tauExtra = (v.tauExtra||0) * ratio;
          v.tauBase = (v.tauBase||0) * ratio;
        }
        v.C = newC;
        e.target.nextElementSibling.textContent = fmtConc(v.C);
        renderTitrantControls();
        renderRigAndStats();
      });
    });
    vesselList.querySelectorAll(".vdel").forEach(btn=>{
      btn.addEventListener("click", e=>{
        if(vessels.length<=1) return; // resta sempre almeno un vaso sul banco
        const id = e.target.dataset.id;
        const idx = vessels.findIndex(v=>v.id===id);
        if(idx<0) return;
        vessels.splice(idx,1);
        const bIdx = Math.max(0, idx-1);
        valves.splice(bIdx,1);
        stacked.splice(bIdx,1);
        renderAll();
      });
    });
    vesselList.querySelectorAll(".stack-sep").forEach(t=>{
      t.addEventListener("click", ()=>{ stacked[+t.dataset.idx] = false; renderAll(); });
    });
    vesselList.querySelectorAll(".stack-valve-mini").forEach(t=>{
      t.addEventListener("click", ()=>{ const i=+t.dataset.idx; if(!isNaN(i)) toggleValve(i); });
    });
    vesselList.querySelectorAll(".col-conn").forEach(c=>{
      c.addEventListener("click", ()=>{ const i=+c.dataset.idx; toggleValve(i); });
    });
  }

  /* ---- trascina e impila: riordina i vasi trascinando le loro schede ----
     Il trascinamento parte SOLO dalla maniglia ⠿ (non dall'intera scheda): così
     regolare lo slider di concentrazione, il nome o il pulsante elimina non
     sposta mai la scheda per errore. Trascinare vaso A su vaso B: A si sposta
     subito dopo B. Se entrambi sono tamponi/indicatori/acqua (mai un vaso forte),
     il nuovo confine diventa anche "impilato" (stessa colonna nel banco). */
  let dragSrcId = null;
  function wireDragForCard(card, v){
    const grab = card.querySelector(".vgrab");
    grab.setAttribute("draggable", "true");
    grab.addEventListener("dragstart", e=>{
      dragSrcId = v.id;
      card.classList.add("dragging");
      if(e.dataTransfer){ e.dataTransfer.effectAllowed = "move"; try{ e.dataTransfer.setData("text/plain", v.id); }catch(err){} }
    });
    grab.addEventListener("dragend", ()=>{ card.classList.remove("dragging"); dragSrcId = null; });
    card.addEventListener("dragover", e=>{
      if(!dragSrcId || dragSrcId===v.id) return;
      e.preventDefault();
      card.classList.add("drop-target");
    });
    card.addEventListener("dragleave", ()=> card.classList.remove("drop-target"));
    card.addEventListener("drop", e=>{
      e.preventDefault();
      card.classList.remove("drop-target");
      if(!dragSrcId || dragSrcId===v.id) return;
      moveAndMaybeStack(dragSrcId, v.id);
      dragSrcId = null;
    });
  }
  function moveAndMaybeStack(srcId, targetId){
    // ricorda lo stato di ogni confine per coppia di id (non per indice, che cambia col
    // riordino): così pile di 3+ vasi già costruite sopravvivono a un nuovo trascinamento
    const boundaryState = new Map();
    for(let i=0;i<vessels.length-1;i++){
      boundaryState.set(vessels[i].id+"|"+vessels[i+1].id, {valve:valves[i], stacked:stacked[i]});
    }
    // se src era già subito dopo target, quel confine va SEMPRE trattato come "nuovo"
    // (si impila per via del trascinamento appena fatto), non riletto dal vecchio stato
    boundaryState.delete(targetId+"|"+srcId);

    const srcIdx = vessels.findIndex(v=>v.id===srcId);
    if(srcIdx<0) return;
    const src = vessels[srcIdx];
    vessels.splice(srcIdx,1);
    const tgtIdx = vessels.findIndex(v=>v.id===targetId);
    if(tgtIdx<0){ vessels.splice(srcIdx,0,src); return; }
    vessels.splice(tgtIdx+1, 0, src);

    // qualunque coppia di vasi può impilarsi (HCl/NaOH inclusi: sono vasi normali, con
    // pKa fuori scala ma non più un caso speciale) — anche 4 o più, uno sopra l'altro
    const newValves = [], newStacked = [];
    for(let i=0;i<vessels.length-1;i++){
      const key = vessels[i].id+"|"+vessels[i+1].id;
      if(boundaryState.has(key)){
        const st = boundaryState.get(key);
        newValves.push(st.valve);
        newStacked.push(st.stacked);
      } else if(i===tgtIdx){
        newValves.push(true); newStacked.push(true);
      } else {
        newValves.push(true); newStacked.push(false);
      }
    }
    valves = newValves;
    stacked = newStacked;

    // i confini "di proprietà" usati per il riempimento indipendente (quando una valvola
    // interna è chiusa) seguono l'ordine reale lungo l'asse (per pKa), non l'ordine in
    // cui i vasi sono stati trascinati: se il nuovo vaso non è nel punto giusto della
    // pila a cui si è unito, si riordina l'intera pila per pKa (l'acqua vale 7).
    const anchorOf = v => v.water ? 7 : primaryPka(v);
    let runStart = tgtIdx, runEnd = tgtIdx+1;
    while(runStart>0 && stacked[runStart-1]) runStart--;
    while(runEnd<vessels.length-1 && stacked[runEnd]) runEnd++;
    if(runEnd>runStart){
      const runVessels = vessels.slice(runStart, runEnd+1);
      const sorted = runVessels.slice().sort((a,b)=>anchorOf(a)-anchorOf(b));
      const alreadySorted = runVessels.every((v,i)=>v===sorted[i]);
      if(!alreadySorted){
        for(let i=0;i<sorted.length;i++) vessels[runStart+i] = sorted[i];
        for(let i=runStart;i<runEnd;i++){ valves[i]=true; stacked[i]=true; }
      }
    }
    renderAll();
  }

  function cToSlider(C){ return ((Math.log10(C)+1)/2)*100; } // 0.1..10 M -> 0..100 (metà = 1 M)
  function sliderToC(s){ return Math.pow(10, -1 + (s/100)*2); }
  function fmtConc(C){
    const mM = C*1000;
    if(mM<1) return (mM*1000).toFixed(0)+" µM";
    return mM>=1000 ? (mM/1000).toFixed(2)+" M" : mM.toFixed(mM<10?2:0)+" mM";
  }

  // fattore di diluizione (v.dil, ≥1): scala log fino a un "quasi infinito" di 10⁷
  // volte il volume originale — abbastanza da avvicinare per davvero il pH a 7, senza
  // mai dichiarare una diluizione infinita vera e propria.
  const DIL_MAX = 1e7;
  function dilToSlider(d){ return (Math.log10(d)/Math.log10(DIL_MAX))*100; }
  function sliderToDil(s){ return Math.pow(10, (s/100)*Math.log10(DIL_MAX)); }
  function fmtDil(d){
    if(d<=1.02) return "non diluito";
    if(d>=1e6) return "× "+(d/1e6).toFixed(d>=1e7?0:1)+"M";
    if(d>=1000) return "× "+Math.round(d/1000)+"k";
    return "× "+Math.round(d);
  }

  /* ============================================================
     GRAFICO beta(pH)  — asse pH orizzontale, convenzione standard
     ============================================================ */
  const chartSvg = document.getElementById("chartSvg");
  const CH = {x0:56, x1:876, y0:210, y1:26};
  function renderChart(levelByVesselId){
    chartSvg.innerHTML = "";
    const buffers = vessels.filter(v=>!v.strong);

    const buffers_ = buffers.filter(v=>!v.water);
    const scaleSource = buffers_.length ? buffers_ : buffers;
    let maxBeta = 1e-9;
    scaleSource.forEach(v=>{
      for(let i=0;i<=112;i++){
        const pH = 14*i/112;
        maxBeta = Math.max(maxBeta, vesselBeta(v,pH));
      }
    });
    maxBeta *= 1.12;

    const X = pH => CH.x0 + (pH/14)*(CH.x1-CH.x0);
    const Y = b => Math.max(CH.y1, CH.y0 - (b/maxBeta)*(CH.y0-CH.y1));

    const g = el("g",{},chartSvg);
    for(let p=0;p<=14;p+=2){
      el("line",{x1:X(p),x2:X(p),y1:CH.y0,y2:CH.y1,stroke:"var(--line)","stroke-width":1}, g);
      const t = el("text",{x:X(p),y:CH.y0+18,"text-anchor":"middle","font-family":"IBM Plex Mono, monospace","font-size":11,fill:"var(--ink-faint)"}, g);
      t.textContent = p;
    }
    const xlab = el("text",{x:(CH.x0+CH.x1)/2, y:CH.y0+38, "text-anchor":"middle","font-family":"IBM Plex Sans, sans-serif","font-size":11,fill:"var(--ink-faint)"}, g);
    xlab.textContent = "pH";
    el("line",{x1:CH.x0,x2:CH.x0,y1:CH.y0,y2:CH.y1,stroke:"var(--line-strong)","stroke-width":1}, g);
    const ylab = el("text",{x:CH.x0-40, y:(CH.y0+CH.y1)/2, "text-anchor":"middle","font-family":"IBM Plex Sans, sans-serif","font-size":11,fill:"var(--ink-faint)","transform":`rotate(-90 ${CH.x0-40} ${(CH.y0+CH.y1)/2})`}, g);
    ylab.textContent = "capacità tampone β";

    buffers.forEach(v=>{
      const pts = [];
      for(let i=0;i<=112;i++){
        const pH = 14*i/112;
        pts.push(`${X(pH)},${Y(vesselBeta(v,pH))}`);
      }
      el("polyline", {points:pts.join(" "), fill:"none", stroke:colorOf(v), "stroke-width": v.water?1.6:2.2, "stroke-dasharray": v.water?"3,3":"none", opacity: v.water?0.75:0.95, "stroke-linejoin":"round"}, chartSvg);
      const level = levelByVesselId[v.id];
      const cx = X(level), cy = Y(vesselBeta(v, level));
      el("circle", {cx, cy, r:4.5, fill:colorOf(v), stroke:"var(--surface)", "stroke-width":1.5}, chartSvg).appendChild(
        (()=>{ const t=document.createElementNS(svgNS,"title"); t.textContent = `${v.label}: pH ${level.toFixed(2)}, β ${(vesselBeta(v,level)*1000).toFixed(2)} mmol/(L·pH)`; return t; })()
      );
    });

    const legend = document.getElementById("chartLegend");
    legend.innerHTML = "";
    buffers.forEach(v=>{
      const item = document.createElement("span");
      item.className = "item";
      const pKas = Array.isArray(v.pKa) ? v.pKa : (v.water?null:[v.pKa]);
      item.innerHTML = `<span class="sw" style="background:${colorOf(v)}${v.water?';opacity:.75':''}"></span>${v.label}${pKas?" (pKa "+pKas.map(p=>p.toFixed(2)).join("/")+")":""}`;
      legend.appendChild(item);
    });
  }

  /* ============================================================
     TITOLANTE: un cursore per ciascun gruppo di vasi collegati
     ============================================================ */
  const titrantGroupsEl = document.getElementById("titrantGroups");

  const STEPS = [0.001, 0.01, 0.1, 1, 10, 100, 1000];
  const stepSelect = document.getElementById("stepSelect");
  STEPS.forEach((s,i)=>{
    const opt = document.createElement("option");
    opt.value = s; opt.textContent = "passo " + s + " mmol/L";
    if(i===2) opt.selected = true;
    stepSelect.appendChild(opt);
  });

  // ricostruisce l'intero pannello (una card per gruppo): usato quando la topologia dei
  // gruppi può essere cambiata (vaso aggiunto/tolto, valvola aperta/chiusa, concentrazione
  // di un vaso modificata, ecc.) — MAI durante il trascinamento del cursore di un singolo
  // gruppo, altrimenti lo si ricrea sotto al puntatore e si interrompe il drag.
  function renderTitrantControls(){
    titrantGroupsEl.innerHTML = "";
    const groups = groupsFromValves();
    groups.forEach(members=>{
      const boundMmol = groupBound(members)*1000;
      const dispMmol = groupDisplayTau(members)*1000;
      const label = members.map(v=>v.short||v.label).join(" + ");
      const dilutable = members.some(v=>!v.water);

      const wrap = document.createElement("div");
      wrap.className = "titrant-group";
      wrap.innerHTML = `
        <div class="titrant-group-label">${label}</div>
        <div class="titrant-row">
          <span class="mono" style="font-size:.68rem;color:var(--ink-faint);">H⁺</span>
          <input type="range" class="tg-slider" min="${-boundMmol}" max="${boundMmol}" step="0.1" value="${Math.max(-boundMmol,Math.min(boundMmol,dispMmol)).toFixed(3)}">
          <span class="mono" style="font-size:.68rem;color:var(--ink-faint);">OH⁻</span>
        </div>
        <div class="step-row" style="margin-top:0;">
          <button type="button" class="tg-minus" title="Togli un passo">−</button>
          <input class="titrant-num mono tg-num" type="number" step="0.001" value="${dispMmol.toFixed(3)}">
          <button type="button" class="tg-plus" title="Aggiungi un passo">+</button>
          <span class="mono" style="font-size:.66rem;color:var(--ink-faint);">mmol/L eq.</span>
        </div>
        ${dilutable ? `
        <div class="row-slider" style="margin-top:8px;">
          <span>Diluizione</span>
          <input type="range" class="tg-dil" draggable="false" min="0" max="100" step="1" value="${dilToSlider(groupDil(members))}">
          <output>${fmtDil(groupDil(members))}</output>
        </div>` : ""}
      `;
      const slider = wrap.querySelector(".tg-slider");
      const num = wrap.querySelector(".tg-num");
      slider.addEventListener("input", ()=>{
        const v = +slider.value;
        setGroupDisplayTau(members, v);
        num.value = v.toFixed(3);
        renderRigAndStats();
      });
      num.addEventListener("input", ()=>{
        const v = parseFloat(num.value);
        if(!isFinite(v)) return;
        setGroupDisplayTau(members, v);
        slider.value = Math.max(+slider.min, Math.min(+slider.max, v));
        renderRigAndStats();
      });
      wrap.querySelector(".tg-minus").addEventListener("click", ()=>{
        setGroupDisplayTau(members, groupDisplayTau(members)*1000 - (+stepSelect.value));
        renderAll();
      });
      wrap.querySelector(".tg-plus").addEventListener("click", ()=>{
        setGroupDisplayTau(members, groupDisplayTau(members)*1000 + (+stepSelect.value));
        renderAll();
      });
      const dilSlider = wrap.querySelector(".tg-dil");
      if(dilSlider){
        dilSlider.addEventListener("input", ()=>{
          const newDil = sliderToDil(+dilSlider.value);
          setGroupDil(members, newDil);
          dilSlider.nextElementSibling.textContent = fmtDil(newDil);
          const boundMmolNow = groupBound(members)*1000;
          slider.min = -boundMmolNow; slider.max = boundMmolNow;
          const dispNow = groupDisplayTau(members)*1000;
          slider.value = Math.max(-boundMmolNow, Math.min(boundMmolNow, dispNow)).toFixed(3);
          num.value = dispNow.toFixed(3);
          renderRigAndStats();
        });
      }
      titrantGroupsEl.appendChild(wrap);
    });
  }

  /* ============================================================
     PRESET / AGGIUNTA VASI
     ============================================================ */
  // raggruppato per famiglia (Fosfato, Carbonato, Amminoacidi...) invece di un unico elenco
  // piatto: con tutte le singole coppie disponibili la lista sarebbe altrimenti troppo lunga
  // da scorrere per trovare quella giusta.
  function fillGroupedSelect(select, presets){
    select.innerHTML = "";
    const groups = [];
    presets.forEach(p=>{
      let g = groups.find(g=>g.name===p.group);
      if(!g){ g = {name:p.group, items:[]}; groups.push(g); }
      g.items.push(p);
    });
    groups.forEach(g=>{
      const og = document.createElement("optgroup");
      og.label = g.name;
      g.items.forEach(p=>{
        const opt = document.createElement("option");
        opt.value = p.key; opt.textContent = p.label;
        og.appendChild(opt);
      });
      select.appendChild(og);
    });
  }

  const presetSelect = document.getElementById("presetSelect");
  fillGroupedSelect(presetSelect, PRESETS.filter(p=>p.menu!==false));
  document.getElementById("addBtn").addEventListener("click", ()=>{
    addVesselFromPreset(presetSelect.value); // aggiunto in coda, collegato al vaso precedente
    renderAll();
  });

  document.getElementById("resetBtn").addEventListener("click", ()=>{
    resetDefault();
    setMode("libero");
    renderAll();
  });

  /* ============================================================
     MODALITÀ: libero / preparazione / titolazione
     ============================================================ */
  const modeTabs = document.getElementById("modeTabs");
  const titrPanel = document.getElementById("titrPanel");
  const benchRow = document.getElementById("benchRow");
  const vesselCard = document.getElementById("vesselCard");
  const prepSideCard = document.getElementById("prepSideCard");

  const prepSystem = document.getElementById("prepSystem");
  fillGroupedSelect(prepSystem, PRESETS.filter(p=>p.cat==="buffer" && !Array.isArray(p.pKa)));

  function setMode(m){
    mode = m;
    modeTabs.querySelectorAll(".mode-tab").forEach(b=>b.classList.toggle("active", b.dataset.mode===m));
    titrPanel.hidden = m!=="titolazione";
    vesselCard.hidden = m==="prep";
    prepSideCard.hidden = m!=="prep";
    benchRow.classList.toggle("titr-layout", m==="titolazione");
    renderRig();
    if(m==="titolazione") renderTitrationChart();
  }
  modeTabs.querySelectorAll(".mode-tab").forEach(b=>{
    b.addEventListener("click", ()=> setMode(b.dataset.mode));
  });

  document.getElementById("prepGoBtn").addEventListener("click", ()=>{
    const key = prepSystem.value;
    const method = document.querySelector('input[name="prepMethod"]:checked').value;
    const p = PRESETS.find(p=>p.key===key);
    if(!p) return;
    const water = { id:"water", water:true, label:"Acqua pura", short:"Acqua", color:"var(--c-water)", tauExtra:0, tauBase:0 };
    const buf = {id:"v"+(nextId++), label:p.label, short:p.short, pKa:p.pKa, C:p.C, color:p.color, species:p.species, tauExtra:0, tauBase:0};
    vessels = [water, buf];
    valves = [true];
    stacked = [false];
    const pk = p.pKa;
    // "dall'acido"/"dalla base" partono dalla forma pressoché pura (±4 unità di pH dal
    // pKa, oltre il 99.99% in quella forma): il titolante riparte da 0 proprio lì, così
    // aggiungendolo si vede subito quanto ne serve per arrivare al tampone desiderato
    let startLevel;
    if(method==="mix") startLevel = pk;
    else if(method==="fromAcid") startLevel = Math.max(0.05, pk-4);
    else startLevel = Math.min(13.95, pk+4);
    consolidateGroup([water, buf], startLevel);
    renderAll();
  });

  /* ============================================================
     TITOLAZIONE: curva pH vs titolante, trascinabile
     ============================================================ */
  const titrSvg = document.getElementById("titrSvg");
  const TT = {x0:112, x1:940, y0:560, y1:80, viewW:980};
  // valore ipotetico di un gruppo se il suo titolante mostrato fosse dispMmol: usato per
  // tracciare la curva (che esplora l'intero asse) senza toccare lo stato reale del gruppo
  function levelAtDisplay(members, dispMmol){
    const saved = members[0].tauExtra;
    members[0].tauExtra += dispMmol/1000 - groupDisplayTau(members);
    const lvl = solveLevel(members);
    members[0].tauExtra = saved;
    return lvl;
  }
  let titrGroupsCache = []; // aggiornata a ogni render, usata per capire quale curva si sta trascinando
  function renderTitrationChart(){
    titrSvg.innerHTML = "";
    const groups = groupsFromValves();
    titrGroupsCache = groups;
    // scala condivisa dell'asse: il massimo fra quanto serve a ciascun gruppo per
    // raggiungere i propri estremi, così ogni curva ci sta comodamente
    const B = Math.max(1, ...groups.map(g=>groupBound(g)*1000));
    const Xt = t => TT.x0 + ((t+B)/(2*B))*(TT.x1-TT.x0);
    const Yt = pH => TT.y0 - (pH/14)*(TT.y0-TT.y1);

    const g = el("g",{},titrSvg);
    for(let pH=0; pH<=14; pH+=2){
      const y = Yt(pH);
      el("line",{x1:TT.x0,x2:TT.x1,y1:y,y2:y,stroke:"var(--line)","stroke-width":1}, g);
      el("text",{x:TT.x0-10,y:y+4.5,"text-anchor":"end","font-family":"IBM Plex Mono, monospace","font-size":13,fill:"var(--ink-faint)"}, g).textContent = pH;
    }
    // quattro tacche eque più lo zero: si adattano all'estremo corrente della scala
    [-B, -B/2, 0, B/2, B].forEach(t=>{
      const x = Xt(t);
      el("text",{x, y:TT.y0+24,"text-anchor":"middle","font-family":"IBM Plex Mono, monospace","font-size":13,fill:"var(--ink-faint)"}, g).textContent = Math.round(t);
    });
    el("text",{x:(TT.x0+TT.x1)/2, y:TT.y0+46,"text-anchor":"middle","font-family":"IBM Plex Sans, sans-serif","font-size":13,fill:"var(--ink-faint)"}, g).textContent = "titolante τ di ciascun gruppo (mmol/L, H⁺ ← 0 → OH⁻)";
    el("text",{x:TT.x0-52, y:(TT.y0+TT.y1)/2,"text-anchor":"middle","font-family":"IBM Plex Sans, sans-serif","font-size":13,fill:"var(--ink-faint)","transform":`rotate(-90 ${TT.x0-52} ${(TT.y0+TT.y1)/2})`}, g).textContent = "pH";

    const legend = document.getElementById("titrLegend");
    legend.innerHTML = "";

    groups.forEach((members, gi)=>{
      const rep = members.find(v=>!v.water && !v.strong) || members[members.length-1];
      const color = colorOf(rep);
      const pts = [];
      for(let i=0;i<=140;i++){
        const t = -B + 2*B*i/140;
        const lvl = levelAtDisplay(members, t);
        pts.push(`${Xt(t)},${Yt(lvl)}`);
      }
      el("polyline", {points:pts.join(" "), fill:"none", stroke:color, "stroke-width":3, "stroke-linejoin":"round"}, titrSvg);

      const curDisp = groupDisplayTau(members)*1000;
      const curLevel = solveLevel(members);
      el("circle", {cx:Xt(curDisp), cy:Yt(curLevel), r:6.5, fill:color, stroke:"var(--surface)", "stroke-width":2, class:"titr-dot", "data-gi":gi}, titrSvg);

      const item = document.createElement("span");
      item.className = "item";
      item.innerHTML = `<span class="sw" style="background:${color}"></span>Gruppo ${gi+1}: ${members.map(v=>v.short||v.label).join(" + ")}`;
      legend.appendChild(item);
    });

    el("rect", {id:"titrHit", x:TT.x0, y:TT.y1, width:TT.x1-TT.x0, height:TT.y0-TT.y1, fill:"transparent"}, titrSvg);
  }

  // listener persistenti sull'elemento <svg> stesso (non ricreato a ogni render, a differenza del
  // rettangolo "hit" interno): così pointerdown/move/up e la pointer capture restano validi durante
  // tutto il trascinamento, anche se renderTitrationChart ridisegna il contenuto a ogni passo.
  let titrDragging = false;
  let titrDragMembers = null; // il gruppo (array di vasi) attualmente trascinato
  function titrPickGroupAt(clientX, clientY){
    const rect = titrSvg.getBoundingClientRect();
    const svgX = (clientX - rect.left) * (TT.viewW / rect.width);
    const svgY = (clientY - rect.top) * (TT.viewW / rect.width); // stessa scala, viewBox quadrata in unità
    let best = null, bestDist = Infinity;
    titrSvg.querySelectorAll(".titr-dot").forEach(dot=>{
      const dx = +dot.getAttribute("cx") - svgX, dy = +dot.getAttribute("cy") - svgY;
      const dist = dx*dx + dy*dy;
      if(dist < bestDist){ bestDist = dist; best = +dot.getAttribute("data-gi"); }
    });
    return best!==null ? titrGroupsCache[best] : (titrGroupsCache[0] || null);
  }
  function titrSetFromClientX(clientX){
    if(!titrDragMembers) return;
    const rect = titrSvg.getBoundingClientRect();
    const svgX = (clientX - rect.left) * (TT.viewW / rect.width);
    const B = Math.max(1, ...titrGroupsCache.map(g=>groupBound(g)*1000));
    let t = ((svgX - TT.x0)/(TT.x1-TT.x0))*(2*B) - B;
    t = Math.max(-B, Math.min(B, t));
    setGroupDisplayTau(titrDragMembers, t);
    renderRigAndStats();
    renderTitrantControls();
  }
  titrSvg.addEventListener("pointerdown", e=>{
    titrDragging = true;
    titrDragMembers = titrPickGroupAt(e.clientX, e.clientY);
    titrSvg.setPointerCapture(e.pointerId);
    titrSetFromClientX(e.clientX);
  });
  titrSvg.addEventListener("pointermove", e=>{ if(titrDragging) titrSetFromClientX(e.clientX); });
  titrSvg.addEventListener("pointerup", ()=>{ titrDragging = false; });
  titrSvg.addEventListener("pointercancel", ()=>{ titrDragging = false; });

  /* ============================================================
     ISTRUZIONI (dialog)
     ============================================================ */
  const instrDialog = document.getElementById("instrDialog");
  document.getElementById("instrBtn").addEventListener("click", ()=> instrDialog.showModal());
  document.getElementById("instrClose").addEventListener("click", ()=> instrDialog.close());
  instrDialog.addEventListener("click", e=>{ if(e.target===instrDialog) instrDialog.close(); });

  // le due illustrazioni nelle istruzioni sono vere istantanee SVG del banco/della curva
  // (non screenshot rasterizzati): si costruisce uno stato dimostrativo, si cattura il
  // markup, poi si ripristina lo stato di partenza dell'app.
  function buildInstrFigures(){
    const savedVessels = vessels, savedValves = valves, savedStacked = stacked, savedMode = mode;

    const water = { id:"water", water:true, label:"Acqua pura", short:"Acqua", color:"var(--c-water)", tauExtra:0, tauBase:0 };
    const acetico = PRESETS.find(p=>p.key==="acetico");
    const fosfato = PRESETS.find(p=>p.key==="fosfato2");
    const demo1 = {id:"demo1", label:acetico.label, short:acetico.short, pKa:acetico.pKa, C:acetico.C, color:acetico.color, species:acetico.species, tauExtra:0, tauBase:0};
    const demo2 = {id:"demo2", label:fosfato.label, short:fosfato.short, pKa:fosfato.pKa, C:fosfato.C, color:fosfato.color, species:fosfato.species, tauExtra:0, tauBase:0};

    vessels = [water, demo1];
    valves = [true];
    stacked = [false];
    mode = "libero";
    renderRig();
    const benchSvg = document.getElementById("rigSvg").outerHTML
      .replace('viewBox="0 0 980 610"', 'viewBox="0 0 980 560"')
      .replace(' id="rigSvg"', "")
      .replace(/clip-col/g, "instr-bench-clip-col")
      .replace(/grad-col/g, "instr-bench-grad-col");
    document.getElementById("figBench").insertAdjacentHTML("afterbegin", benchSvg);

    vessels = [water, demo1, demo2];
    valves = [true, true];
    stacked = [false, false];
    mode = "specie";
    renderRig();
    // "clip-col*"/"grad-col*" sono id generici riusati a ogni render: rinominati qui perché
    // non collidano con quelli del banco vero (un url(#id) duplicato punterebbe al primo
    // trovato nel documento, quindi al banco vero anziché a questa istantanea statica)
    const specieSvg = document.getElementById("rigSvg").outerHTML
      .replace('viewBox="0 0 980 610"', 'viewBox="0 0 980 560"')
      .replace(' id="rigSvg"', "")
      .replace(/clip-col/g, "instr-clip-col")
      .replace(/grad-col/g, "instr-grad-col");
    document.getElementById("figSpecie").insertAdjacentHTML("afterbegin", specieSvg);

    mode = "titolazione";
    setGroupDisplayTau([water, demo1, demo2], -groupBound([water, demo1, demo2])*1000*0.22);
    renderRig();
    renderTitrationChart();
    const titrSvgMarkup = document.getElementById("titrSvg").outerHTML
      .replace(' id="titrSvg"', "")
      .replace(' id="titrHit"', "");
    document.getElementById("figTitr").insertAdjacentHTML("afterbegin", titrSvgMarkup);

    vessels = savedVessels; valves = savedValves; stacked = savedStacked; mode = savedMode;
    setMode(mode);
    renderAll();
  }
  buildInstrFigures();
  // le istruzioni sono la prima cosa che si vede aprendo la pagina: si chiudono e si inizia
  instrDialog.showModal();

  /* ============================================================
     TEMA
     ============================================================ */
  const themeBtn = document.getElementById("themeToggle");
  const THEMES = ["system","light","dark"];
  const THEME_LABEL = {system:"Tema: sistema", light:"Tema: chiaro", dark:"Tema: scuro"};
  let themeIdx = 0;
  themeBtn.addEventListener("click", ()=>{
    themeIdx = (themeIdx+1)%THEMES.length;
    const t = THEMES[themeIdx];
    if(t==="system") document.documentElement.removeAttribute("data-theme");
    else document.documentElement.setAttribute("data-theme", t);
    themeBtn.textContent = THEME_LABEL[t];
    renderAll();
  });

  /* ============================================================
     LOOP DI RENDER
     ============================================================ */
  function renderRigAndStats(){
    const {levelByVesselId, groupInfo} = renderRig();
    renderStats(groupInfo);
    renderChart(levelByVesselId);
    if(mode==="titolazione") renderTitrationChart();
  }
  function renderAll(){
    renderVesselControls();
    renderTitrantControls();
    renderRigAndStats();
  }

  // se la finestra viene ridimensionata mentre il banco è oltre le 5 colonne (larghezza fissa
  // in pixel, misurata al contenitore), rimisura e ridisegna così la larghezza resta coerente.
  let resizeT = null;
  window.addEventListener("resize", ()=>{
    clearTimeout(resizeT);
    resizeT = setTimeout(renderRig, 120);
  });

  setMode("libero");
  renderAll();
})();
