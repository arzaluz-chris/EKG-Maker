"use strict";
/* =======================================================================
   EKG MAKER — simulador de electrocardiograma
   Módulos: Util · Morf · Derivaciones · Ritmo · Senal · Animacion · Papel · Presets · UI
   ======================================================================= */

/* ============================== UTIL ============================== */
function mulberry32(a){
  return function(){
    a|=0; a=a+0x6D2B79F5|0;
    let t=Math.imul(a^a>>>15,1|a);
    t=t+Math.imul(t^t>>>7,61|t)^t;
    return ((t^t>>>14)>>>0)/4294967296;
  };
}
const clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
const grados=g=>g*Math.PI/180;
const copia=o=>JSON.parse(JSON.stringify(o));
/* mezcla profunda: parte de la plantilla y sobrescribe con lo que traiga el archivo,
   de modo que un JSON incompleto no deje claves indefinidas */
function fusionarEstado(base,dato){
  if(dato===null||dato===undefined) return copia(base);
  const objeto=v=>v&&typeof v==="object"&&!Array.isArray(v);
  if(!objeto(base)||!objeto(dato)) return copia(dato);
  const salida=copia(base);
  for(const k of Object.keys(dato)) salida[k]=fusionarEstado(base[k],dato[k]);
  return salida;
}

/* ============================== MORFOLOGÍA ============================== */
const E0=Math.exp(-4.5);
/* campana asimétrica normalizada a 0 en los extremos; apex 0..1 */
function campana(dt,dur,apex){
  if(dt<=0||dt>=dur) return 0;
  const pico=dur*apex;
  const s=(dt<pico?pico:dur-pico)/3;
  if(s<=0) return 0;
  const e=Math.exp(-0.5*Math.pow((dt-pico)/s,2));
  return Math.max(0,(e-E0)/(1-E0));
}
/* interpolación lineal sobre una lista de vértices [[t,mV],...] */
function interpVert(v,dt){
  if(dt<=v[0][0]||dt>=v[v.length-1][0]) return 0;
  for(let i=1;i<v.length;i++){
    if(dt<=v[i][0]){
      const a=v[i-1],b=v[i],f=(dt-a[0])/(b[0]-a[0]||1e-9);
      return a[1]+(b[1]-a[1])*f;
    }
  }
  return 0;
}
/* muestreo suavizado: da el aspecto empastado de los complejos anchos */
function vertSuave(v,dt,suave){
  if(!suave) return interpVert(v,dt);
  let s=0;const n=5;
  for(let i=0;i<n;i++) s+=interpVert(v,dt+(i/(n-1)-0.5)*2*suave);
  return s/n;
}

/* ============================== DERIVACIONES ============================== */
const DERIVACIONES=["I","II","III","aVR","aVL","aVF","V1","V2","V3","V4","V5","V6"];
const ANGULO={I:0,II:60,III:120,aVR:-150,aVL:-30,aVF:90};
const ES_PRECORDIAL=d=>d.charAt(0)==="V";
const IDX_V=d=>+d.slice(1)-1;

/* factores empíricos por derivación precordial, relativos a la amplitud base */
const PERFILES={
  normal:{
    P:[0.30,0.50,0.50,0.50,0.45,0.35],
    Q:[0.00,0.00,0.10,0.25,0.50,0.60],
    R:[0.18,0.40,0.70,1.05,1.10,0.85],
    S:[1.40,1.50,0.95,0.45,0.22,0.12],
    T:[-0.15,0.45,0.60,0.60,0.50,0.38],
    U:[0.30,0.40,0.40,0.35,0.30,0.25]
  },
  mala:{
    P:[0.30,0.50,0.50,0.50,0.45,0.35],
    Q:[0,0,0,0.10,0.25,0.35],
    R:[0.10,0.13,0.20,0.35,0.60,0.80],
    S:[1.30,1.45,1.30,1.00,0.60,0.30],
    T:[-0.15,0.35,0.45,0.50,0.45,0.35],
    U:[0.30,0.40,0.40,0.35,0.30,0.25]
  },
  hvi:{
    P:[0.35,0.55,0.55,0.55,0.50,0.40],
    Q:[0,0,0.10,0.25,0.50,0.60],
    R:[0.15,0.30,0.60,1.25,1.95,1.70],
    S:[2.20,2.40,1.40,0.60,0.25,0.15],
    T:[-0.10,0.30,0.35,0.10,-0.45,-0.40],
    U:[0.30,0.40,0.40,0.35,0.30,0.25]
  },
  hvd:{
    P:[0.55,0.60,0.50,0.45,0.40,0.30],
    Q:[0,0,0,0.10,0.20,0.30],
    R:[1.20,0.85,0.55,0.40,0.35,0.30],
    S:[0.30,0.60,0.90,1.00,0.90,0.70],
    T:[-0.45,-0.40,-0.25,0.15,0.30,0.30],
    U:[0.25,0.30,0.30,0.25,0.20,0.20]
  }
};

/* territorios para desplazamiento del ST (valor = fracción de la elevación) */
const TERRITORIOS={
  ninguno:{},
  septal:{V1:1,V2:1,V3:.4,II:-.3,III:-.4,aVF:-.35},
  anteroseptal:{V1:.5,V2:1,V3:1,V4:.35,II:-.28,III:-.32,aVF:-.30},
  anterior:{V1:.8,V2:1,V3:1,V4:.9,aVL:.4,I:.35,II:-.4,III:-.5,aVF:-.45},
  anterior_extenso:{V1:.7,V2:1,V3:1,V4:1,V5:.8,V6:.6,I:.6,aVL:.6,II:-.4,III:-.6,aVF:-.5},
  lateral:{I:1,aVL:1,V5:.9,V6:.9,II:-.3,III:-.6,aVF:-.5,V1:-.3},
  inferior:{II:.9,III:1,aVF:1,I:-.5,aVL:-.7,V1:-.2,V2:-.2},
  posterior:{V1:-1,V2:-1,V3:-.7,II:.35,III:.4,aVF:.4,V6:.3},
  vd:{V1:1,III:.9,aVF:.6,II:.5,I:-.4,aVL:-.5},
  difusa:{I:.7,II:.9,III:.3,aVL:.5,aVF:.7,V2:.6,V3:.8,V4:.9,V5:.8,V6:.7,aVR:-.9}
};

/* ============================== ESTADO ============================== */
const ESTADO_BASE={
  semilla:12345,
  paciente:{nombre:"Paciente Demo",edad:"58",sexo:"Masculino",expediente:"CASO-001",
            fecha:"",hora:"",interpretacion:""},
  ondas:{
    P :{amp:0.15,dur:100,apex:0.50},
    PR:{dur:160},
    Q :{amp:-0.08,dur:20},
    R :{amp:1.10,dur:48},
    S :{amp:-0.25,dur:30},
    ST:{nivel:0.00,dur:90,curva:0.0},
    T :{amp:0.35,dur:170,apex:0.62},
    TP:{dur:0},
    U :{amp:0.04,dur:90}
  },
  ejes:{p:50,qrs:45,t:40},
  ritmo:{
    origenAuricular:"sinusal",   // sinusal | flutter | fibrilacion | ninguno
    fcA:70, fcV:70,
    conduccion:"normal",         // normal | wenckebach | mobitz2 | completo | independiente
    relacion:2,                  // para mobitz2 / flutter
    irregularidad:0,             // 0..1
    variabilidad:0.02,           // variabilidad sinusal
    escape:"union",              // union | ventricular
    ectopia:"ninguna",           // ninguna | ev | bigeminismo | trigeminismo | dupletas | ea
    cadaEctopico:6,
    modo:"normal"                // normal | fv | asistolia | torsades
  },
  morfologia:{
    qrs:"normal",                // normal | brd | bri | ventricular | wpw | marcapasos
    perfil:"normal",             // progresión precordial
    hbai:false, hbpi:false,
    osborn:false, delta:false
  },
  isquemia:{territorio:"ninguno",elevacion:0.30,fase:"aguda"}, // aguda|hiperaguda|q|cronica
  /* desnivel del ST propio de cada derivación, en mV; se suma al nivel global y al territorio.
     Sirve para la discordancia secundaria del BRI, el marcapasos o los ritmos ventriculares. */
  stDeriv:{I:0,II:0,III:0,aVR:0,aVL:0,aVF:0,V1:0,V2:0,V3:0,V4:0,V5:0,V6:0},
  /* ajuste opcional y disperso de la morfología de una derivación (solo las que lo necesitan):
     {V2:{q:mV que se suma a la Q (negativo = más profunda), qDur:ms, r:× amplitud de R, rDur:ms,
          s:× amplitud de S, t:mV que se suma a la T}}. Vacío = sin cambios; los casos antiguos no lo traen. */
  ajusteDeriv:{},
  artefactos:{basal:0,muscular:0,red:0},
  papel:{tamano:"carta-h",formato:"tira",duracion:10,velocidad:25,ganancia:10,tiraRitmo:true,
         mostrarMarca:true,mostrarAviso:true,mostrarTitulo:true},
  /* formato: "tira" = una derivación por fila | "4x3" = las 12 en una hoja, 4 columnas × 3 filas */
  monitor:{canales:1,deriv:"II",velocidad:25,ganancia:10}
};
let E=copia(ESTADO_BASE);

/* ============================== MORFOLOGÍA POR DERIVACIÓN ============================== */
/* factor de proyección para las derivaciones de miembros */
function proyMiembro(deriv,eje){ return Math.cos(grados(ANGULO[deriv]-eje)); }

/* para precordiales la tabla es multiplicador directo sobre la amplitud base */
function ampPrecordial(comp,deriv,e){
  return {P:e.ondas.P.amp,T:e.ondas.T.amp,U:e.ondas.U.amp}[comp]*PERFILES[e.morfologia.perfil][comp][IDX_V(deriv)]*2;
}
function ampOnda(comp,deriv,e){
  if(ES_PRECORDIAL(deriv)) return ampPrecordial(comp,deriv,e);
  const eje={P:e.ejes.p,T:e.ejes.t,U:e.ejes.t}[comp];
  return {P:e.ondas.P.amp,T:e.ondas.T.amp,U:e.ondas.U.amp}[comp]*proyMiembro(deriv,eje);
}

/* ---- construcción del QRS como lista de vértices (segundos, mV) ---- */
function verticesQRS(deriv,e,opciones){
  const o=e.ondas, m=e.morfologia, op=opciones||{};
  const tipo=op.tipo||m.qrs;
  const s=1/1000;
  let q,r,sw,dq,dr,ds,suave=0,vert=[],dur;

  /* factores por derivación */
  let fq,fr,fs;
  if(ES_PRECORDIAL(deriv)){
    const i=IDX_V(deriv),P=PERFILES[m.perfil];
    fq=P.Q[i]; fr=P.R[i]; fs=P.S[i];
  }else{
    let k=proyMiembro(deriv,e.ejes.qrs);
    /* hemibloqueos: rotan el eje pero también dejan qR / rS típicos */
    fq=k*0.9; fr=k; fs=k*0.9;
  }

  dq=o.Q.dur*s; dr=o.R.dur*s; ds=o.S.dur*s;
  q=o.Q.amp*Math.abs(fq)*(ES_PRECORDIAL(deriv)?1:Math.sign(fq)||1);
  r=o.R.amp*fr;
  sw=o.S.amp*Math.abs(fs)*(ES_PRECORDIAL(deriv)?1:Math.sign(fr)||1);
  const aj=ajusteDe(deriv,e);
  if(aj){
    if(typeof aj.q==="number") q+=aj.q;
    if(typeof aj.qDur==="number") dq=aj.qDur*s;
    if(typeof aj.r==="number") r*=aj.r;
    if(typeof aj.rDur==="number") dr=aj.rDur*s;
    if(typeof aj.s==="number") sw*=aj.s;
  }

  const armar=comps=>{
    let t=0;const v=[[0,0]];
    for(const c of comps){
      if(c.dur<=0) continue;
      if(c.mellado){
        v.push([t+c.dur*0.28,c.amp]);
        v.push([t+c.dur*0.50,c.amp*0.72]);
        v.push([t+c.dur*0.72,c.amp*1.0]);
      }else{
        v.push([t+c.dur*(c.apex||0.5),c.amp]);
      }
      t+=c.dur; v.push([t,0]);
    }
    return {v,dur:t};
  };

  if(tipo==="brd"){
    suave=0.006;
    if(deriv==="V1"||deriv==="V2"){
      ({v:vert,dur}=armar([{amp:0.25,dur:0.022},{amp:-0.30,dur:0.028},{amp:0.85,dur:0.075}]));
    }else if(deriv==="V3"){
      ({v:vert,dur}=armar([{amp:0.45,dur:0.03},{amp:-0.55,dur:0.035},{amp:0.35,dur:0.055}]));
    }else if(deriv==="I"||deriv==="aVL"||deriv==="V5"||deriv==="V6"){
      ({v:vert,dur}=armar([{amp:q,dur:dq},{amp:r*0.9,dur:dr},{amp:-0.32,dur:0.062}]));
    }else{
      ({v:vert,dur}=armar([{amp:q,dur:dq},{amp:r,dur:dr},{amp:sw*1.4,dur:0.055}]));
    }
  }else if(tipo==="bri"){
    suave=0.008;
    if(deriv==="V1"||deriv==="V2"||deriv==="V3"){
      ({v:vert,dur}=armar([{amp:0.06,dur:0.022},{amp:-1.25,dur:0.11,apex:0.55}]));
    }else if(deriv==="I"||deriv==="aVL"||deriv==="V5"||deriv==="V6"){
      ({v:vert,dur}=armar([{amp:1.35,dur:0.132,mellado:true}]));
    }else{
      const k=ES_PRECORDIAL(deriv)?0.5:proyMiembro(deriv,e.ejes.qrs);
      ({v:vert,dur}=armar([{amp:1.1*k,dur:0.132,mellado:true}]));
    }
  }else if(tipo==="ventricular"){
    suave=0.010;
    /* origen ventricular izquierdo: QRS ancho, negativo en precordiales derechas */
    let k;
    if(ES_PRECORDIAL(deriv)){ k=[-1,-1,-0.7,0.2,0.9,1.0][IDX_V(deriv)]; }
    else k=proyMiembro(deriv,e.ejes.qrs-70);
    ({v:vert,dur}=armar([{amp:1.6*k,dur:0.14,mellado:true}]));
  }else if(tipo==="marcapasos"){
    suave=0.009;
    let k;
    if(ES_PRECORDIAL(deriv)) k=[-1,-1,-0.8,0.1,0.8,0.9][IDX_V(deriv)];
    else k=proyMiembro(deriv,e.ejes.qrs-60);
    ({v:vert,dur}=armar([{amp:1.4*k,dur:0.13,mellado:true}]));
    /* espiga del marcapasos */
    vert=[[0,0],[0.002,2.2*(ES_PRECORDIAL(deriv)?0.6:1)],[0.005,0]].concat(vert.map(p=>[p[0]+0.008,p[1]]));
    dur+=0.008;
  }else if(tipo==="wpw"||m.delta){
    /* preexcitación: la onda delta es un empastamiento progresivo (cóncavo hacia arriba) del inicio
       del QRS que se funde con la rama ascendente de la R, sin volver a la línea de base; sustituye a
       la q septal y alarga el QRS a ~130 ms */
    suave=0.006;
    const dDelta=0.050, hDelta=0.30;
    vert=[[0,0]];
    for(let k=1;k<=6;k++){ const u=k/6; vert.push([dDelta*u,r*hDelta*Math.pow(u,1.4)]); }
    vert.push([dDelta+dr*0.55,r],[dDelta+dr,0]);
    if(ds>0) vert.push([dDelta+dr+ds*0.5,sw],[dDelta+dr+ds,0]);
    dur=dDelta+dr+Math.max(0,ds);
  }else{
    ({v:vert,dur}=armar([{amp:q,dur:dq},{amp:r,dur:dr},{amp:sw,dur:ds}]));
  }

  /* hemibloqueos: modifican el eje, ya cubierto por e.ejes.qrs; onda J de Osborn */
  if(m.osborn){
    const k=ES_PRECORDIAL(deriv)?[0.4,0.5,0.7,0.9,1,0.9][IDX_V(deriv)]:Math.max(0,proyMiembro(deriv,e.ejes.qrs));
    vert=vert.concat([[dur+0.012,0.28*k],[dur+0.030,0.15*k],[dur+0.045,0]]);
    dur+=0.045;
  }
  return {vert,dur,suave};
}

/* ajuste propio de la derivación (ver ESTADO_BASE.ajusteDeriv), o null */
function ajusteDe(deriv,e){
  const a=e.ajusteDeriv&&e.ajusteDeriv[deriv];
  return (a&&typeof a==="object")?a:null;
}
/* desplazamiento del punto J por isquemia, en mV, para una derivación */
function nivelST(deriv,e){
  const is=e.isquemia;
  let n=e.ondas.ST.nivel;
  if(e.stDeriv&&typeof e.stDeriv[deriv]==="number") n+=e.stDeriv[deriv]; // desnivel propio de la derivación
  if(is.territorio!=="ninguno"){
    const f=TERRITORIOS[is.territorio][deriv];
    if(f){
      let mult=1;
      if(is.fase==="hiperaguda") mult=0.45;
      else if(is.fase==="q") mult=0.55;
      else if(is.fase==="cronica") mult=0.1;
      n+=is.elevacion*f*mult;
    }
  }
  return n;
}
/* Discordancia secundaria del ST: con el QRS ancho (BRI, marcapasos, origen ventricular)
   el ST se desvía en sentido opuesto al predominio del QRS y de forma proporcional a él.
   Devuelve un juego de desniveles por derivación, en mV. */
const MORFOLOGIA_ANCHA=["bri","brd","ventricular","marcapasos"];
function discordanciaSecundaria(e,fraccion){
  const f=(typeof fraccion==="number")?fraccion:0.15;
  const salida={};
  for(const d of DERIVACIONES){
    const v=verticesQRS(d,e,{tipo:e.morfologia.qrs}).vert.map(p=>p[1]);
    const neto=Math.max(...v)+Math.min(...v);          // predominio positivo o negativo del QRS
    salida[d]=Math.max(-0.5,Math.min(0.5,+(-neto*f).toFixed(3)));
  }
  return salida;
}

/* onda T modificada por la fase evolutiva del infarto */
function ampT(deriv,e){
  let t=ampOnda("T",deriv,e);
  const is=e.isquemia;
  if(is.territorio!=="ninguno"){
    const f=TERRITORIOS[is.territorio][deriv];
    if(f){
      if(is.fase==="hiperaguda") t+=Math.abs(is.elevacion)*f*1.6;
      else if(is.fase==="q")     t-=is.elevacion*f*1.2;
      else if(is.fase==="cronica")t-=is.elevacion*f*2.2;
    }
  }
  const aj=ajusteDe(deriv,e);
  if(aj&&typeof aj.t==="number") t+=aj.t;
  return t;
}
/* onda Q patológica en fase Q / crónica */
function factorQPatologica(deriv,e){
  const is=e.isquemia;
  if(is.territorio==="ninguno"||is.fase==="aguda"||is.fase==="hiperaguda") return 0;
  const f=TERRITORIOS[is.territorio][deriv];
  return (f&&f>0.5)?(is.fase==="q"?1:1.2):0;
}

/* ============================== RITMO ============================== */
/* Devuelve lista de eventos ordenados: {t, tipo:'P'|'V', ...} */
function generarEventos(e,duracion){
  const r=e.ritmo, ev=[], rnd=mulberry32(e.semilla);
  const margen=1.5;                       // se genera algo antes y después del rango visible
  const t0=-margen, t1=duracion+margen;

  if(r.modo==="fv"||r.modo==="asistolia") return ev;

  const pAuric=[];
  if(r.origenAuricular==="sinusal"){
    let t=t0, i=0;
    const rrBase=60/Math.max(15,r.fcA);
    while(t<t1){
      const esEA=(r.ectopia==="ea"&&i>0&&i%Math.max(2,r.cadaEctopico)===0);
      if(esEA){
        const tEA=t-rrBase*0.32;
        pAuric.push({t:tEA,ectopico:true});
        t=tEA+rrBase;                       // pausa no compensatoria: se reinicia el nodo sinusal
      }else{
        pAuric.push({t,ectopico:false});
        t+=rrBase*(1+(rnd()*2-1)*r.variabilidad);
      }
      i++;
    }
  }

  const vent=[];
  if(r.conduccion==="normal"||r.conduccion==="wenckebach"||r.conduccion==="mobitz2"){
    const N=Math.max(2,Math.round(r.relacion));
    let i=0;
    for(const p of pAuric){
      let pr=e.ondas.PR.dur/1000, conduce=true;
      if(r.conduccion==="wenckebach"){
        const k=i%N;
        if(k===N-1) conduce=false; else pr=pr+k*0.055;
      }else if(r.conduccion==="mobitz2"){
        if(i%N===N-1) conduce=false;
      }
      if(p.ectopico) pr=pr*1.15;
      if(conduce) vent.push({t:p.t+pr,pr,ectopico:p.ectopico,morf:e.morfologia.qrs});
      i++;
    }
  }else{
    /* ventrículo independiente: FA, flutter con conducción variable, TV, escape, marcapasos */
    let t=t0;
    const rrBase=60/Math.max(15,r.fcV);
    let morf=e.morfologia.qrs;
    if(r.conduccion==="completo") morf=(r.escape==="ventricular")?"ventricular":e.morfologia.qrs;
    while(t<t1){
      vent.push({t,pr:null,ectopico:false,morf});
      let rr=rrBase;
      if(r.irregularidad>0) rr=rrBase*(1+(rnd()*2-1)*r.irregularidad);
      else rr=rrBase*(1+(rnd()*2-1)*r.variabilidad);
      t+=Math.max(0.16,rr);
    }
  }

  /* extrasístoles ventriculares */
  if(["ev","bigeminismo","trigeminismo","dupletas"].includes(r.ectopia)&&vent.length>2){
    const paso={bigeminismo:2,trigeminismo:3,dupletas:5,ev:Math.max(2,r.cadaEctopico)}[r.ectopia];
    for(let i=1;i<vent.length;i++){
      if(i%paso===0){
        const rr=vent[i].t-vent[i-1].t;
        vent[i]={t:vent[i].t-rr*0.36,pr:null,ectopico:true,morf:"ventricular"};
        if(r.ectopia==="dupletas"&&i+1<vent.length){
          vent[i+1]={t:vent[i].t+0.34,pr:null,ectopico:true,morf:"ventricular"};
        }
      }
    }
  }

  for(const p of pAuric) ev.push({t:p.t,tipo:"P",ectopico:p.ectopico});
  for(const v of vent)   ev.push(Object.assign({tipo:"V"},v));
  ev.sort((a,b)=>a.t-b.t);
  return ev;
}

/* frecuencia ventricular medida sobre los eventos generados */
function medirFC(eventos,duracion){
  const t=eventos.filter(x=>x.tipo==="V"&&x.t>=0&&x.t<=duracion).map(x=>x.t);
  if(t.length<2) return 0;
  const rr=(t[t.length-1]-t[0])/(t.length-1);
  return rr>0?Math.round(60/rr):0;
}

/* ============================== SEÑAL ============================== */
/* Parámetros precalculados por derivación (evita recalcular en cada muestra) */
function paramsDeriv(deriv,e){
  const o=e.ondas;
  const qrsNormal=verticesQRS(deriv,e,{tipo:e.morfologia.qrs});
  const qrsVent  =verticesQRS(deriv,e,{tipo:"ventricular"});
  const qPat=factorQPatologica(deriv,e);
  if(qPat>0){
    /* onda Q patológica: profunda (>25 % de la R) y ancha (≥40 ms) */
    const rMax=Math.max(...qrsNormal.vert.map(p=>p[1]));
    const prof=-Math.max(0.25*Math.abs(rMax),0.25)*qPat;
    const desp=0.042;
    qrsNormal.vert=[[0,0],[desp*0.55,prof],[desp,0]].concat(
      qrsNormal.vert.filter(p=>p[0]>0).map(p=>[p[0]+desp,p[1]*(qPat>1?0.45:0.7)]));
    qrsNormal.dur+=desp;
  }
  return {
    ampP:ampOnda("P",deriv,e), durP:o.P.dur/1000, apexP:o.P.apex,
    ampT:ampT(deriv,e),        durT:o.T.dur/1000, apexT:o.T.apex,
    ampU:ampOnda("U",deriv,e), durU:o.U.dur/1000,
    st:nivelST(deriv,e), durST:o.ST.dur/1000, curvaST:o.ST.curva,
    qrsNormal, qrsVent,
    esV1:deriv==="V1",
    fFlutter:(["II","III","aVF"].includes(deriv))?-1:(deriv==="V1"?0.7:(deriv==="aVR"?0.8:0.35)),
    fFib:(deriv==="V1")?1:(["II","III","aVF"].includes(deriv)?0.5:0.3)
  };
}

/* complejo ventricular completo: QRS + ST + T + U */
function valorComplejo(dt,p,ev){
  const ectop=ev.ectopico&&ev.morf==="ventricular";
  const q=(ev.morf==="ventricular")?p.qrsVent:p.qrsNormal;
  const durQ=q.dur;
  if(dt<0) return 0;

  let st=p.st, aT=p.ampT, dT=p.durT;
  if(ectop){                       /* T discordante y ST discordante en la extrasístole */
    const pico=q.vert.reduce((m,x)=>Math.abs(x[1])>Math.abs(m)?x[1]:m,0);
    const sg=pico>=0?-1:1;
    aT=sg*Math.max(0.35,Math.abs(p.ampT)*1.4);
    st=sg*0.12; dT=0.20;
  }
  /* la rama terminal del QRS aterriza en el punto J desplazado, sin escalón */
  if(dt<durQ){
    const r0=durQ*0.6;
    const rampa=dt>r0?st*(dt-r0)/(durQ-r0):0;
    return vertSuave(q.vert,dt,q.suave)+rampa;
  }
  const durST=Math.max(0.02,p.durST);
  /* ST elevado con T positiva: el segmento ST y la T forman una sola curva que sale del punto J,
     sube hasta el vértice de la T y baja con la rama descendente de la T. La curvatura del ST
     decide la forma de la subida: >0 convexa (lesión; «lápida» cuando el ST domina), 0 rectilínea,
     <0 cóncava (pericarditis). Así no aparece una joroba del ST separada de la T. */
  if(st>0.08&&aT>=0&&!ectop){
    /* el vértice cae algo antes que el de la T aislada y la rama descendente es más larga y suave,
       como en el registro real; si el ST domina, el vértice queda apenas por encima del punto J */
    const xp=durST+dT*p.apexT*0.75, W=xp+dT*(1-p.apexT)*1.5, x=dt-durQ;
    if(x<W){
      const pico=Math.max(st*1.08,st*0.6+aT);
      if(x<=xp){
        const u=x/xp, c=Math.max(-1,Math.min(1,p.curvaST));
        const h=(1-c)*u+c*Math.sin(u*Math.PI/2);
        return st+(pico-st)*h;
      }
      const v=(x-xp)/(W-xp);
      return pico*(1+Math.cos(Math.PI*v))/2;
    }
    const u0=durQ+W+0.06;
    if(dt<u0+p.durU) return p.ampU*campana(dt-u0,p.durU,0.5);
    return 0;
  }
  if(dt<durQ+durST){
    const f=(dt-durQ)/durST;
    /* curva: 0 recto, >0 convexo (lesión), <0 cóncavo (pericarditis) */
    const c=p.curvaST;
    return st*(1-f*0.15)+st*c*Math.sin(Math.PI*f)*0.5;
  }
  const t0=durQ+durST;
  if(dt<t0+dT){
    const f=(dt-t0)/dT;
    return st*(0.85)*(1-f)+aT*campana(dt-t0,dT,p.apexT);
  }
  const u0=t0+dT+0.06;
  if(dt<u0+p.durU) return p.ampU*campana(dt-u0,p.durU,0.5);
  return 0;
}
function duracionComplejo(p){
  return p.qrsNormal.dur+p.durST+p.durT+0.06+p.durU+0.02;
}

/* genera un arreglo de muestras para una derivación */
function generarSenal(e,eventos,deriv,duracion,fs){
  const n=Math.round(duracion*fs), out=new Float32Array(n);
  const p=paramsDeriv(deriv,e), r=e.ritmo;
  const rnd=mulberry32(e.semilla+DERIVACIONES.indexOf(deriv)*977+1);
  const durC=duracionComplejo(p);

  /* modos caóticos */
  if(r.modo==="fv"){
    const f1=4.5+rnd()*2, f2=6.5+rnd()*2, fase=rnd()*6;
    for(let i=0;i<n;i++){
      const t=i/fs;
      const amp=(0.45+0.35*Math.sin(2*Math.PI*0.35*t+fase))*(deriv==="aVR"?-0.8:1);
      out[i]=amp*(Math.sin(2*Math.PI*f1*t)+0.6*Math.sin(2*Math.PI*f2*t+1.3)+0.3*Math.sin(2*Math.PI*(f1+3.1)*t));
    }
  }else if(r.modo==="asistolia"){
    for(let i=0;i<n;i++) out[i]=0;
  }else{
    for(const ev of eventos){
      if(ev.tipo==="P"){
        if(r.origenAuricular!=="sinusal") continue;
        const amp=p.ampP*(ev.ectopico?0.75:1)*(ev.ectopico&&p.esV1?-1:1);
        const i0=Math.max(0,Math.ceil(ev.t*fs)), i1=Math.min(n,Math.ceil((ev.t+p.durP)*fs));
        for(let i=i0;i<i1;i++){
          const dt=i/fs-ev.t;
          let v=campana(dt,p.durP,p.apexP);
          if(p.esV1) v=campana(dt,p.durP*0.5,0.5)-0.55*campana(dt-p.durP*0.5,p.durP*0.5,0.5); // P bifásica
          out[i]+=amp*v;
        }
        /* descenso del PR en la pericarditis */
        if(e.isquemia.territorio==="difusa"){
          const f=TERRITORIOS.difusa[deriv]||0;
          const j0=Math.max(0,Math.ceil((ev.t+p.durP)*fs)), j1=Math.min(n,Math.ceil((ev.t+p.durP+0.09)*fs));
          for(let i=j0;i<j1;i++) out[i]-=0.08*f;
        }
      }else{
        const i0=Math.max(0,Math.ceil(ev.t*fs)), i1=Math.min(n,Math.ceil((ev.t+durC+0.2)*fs));
        for(let i=i0;i<i1;i++) out[i]+=valorComplejo(i/fs-ev.t,p,ev);
      }
    }
    /* torsades: modulación del eje que da el giro de puntas */
    if(r.modo==="torsades"){
      for(let i=0;i<n;i++) out[i]*=Math.sin(2*Math.PI*(i/fs)/2.4)*0.9+0.15;
    }
    /* actividad auricular continua */
    if(r.origenAuricular==="flutter"){
      const f=r.fcA/60;
      for(let i=0;i<n;i++){
        const x=((i/fs)*f)%1;
        out[i]+=0.20*p.fFlutter*(x<0.7?(x/0.7)*2-1:1-((x-0.7)/0.3)*2);
      }
    }else if(r.origenAuricular==="fibrilacion"){
      const fr=[];
      for(let k=0;k<6;k++) fr.push({f:5+rnd()*6,a:rnd(),ph:rnd()*6.28});
      for(let i=0;i<n;i++){
        let s=0; for(const o of fr) s+=o.a*Math.sin(2*Math.PI*o.f*(i/fs)+o.ph);
        out[i]+=0.045*p.fFib*s;
      }
    }
  }

  /* artefactos */
  const a=e.artefactos;
  if(a.basal>0){
    for(let i=0;i<n;i++) out[i]+=a.basal*0.25*Math.sin(2*Math.PI*0.22*(i/fs)+deriv.length);
  }
  if(a.muscular>0){
    let s=0;
    for(let i=0;i<n;i++){ s=s*0.55+(rnd()*2-1); out[i]+=a.muscular*0.06*s; }
  }
  if(a.red>0){
    for(let i=0;i<n;i++) out[i]+=a.red*0.035*Math.sin(2*Math.PI*60*(i/fs));
  }
  return out;
}

/* ============================== MEDICIONES ============================== */
function mediciones(e,eventos,duracion){
  const o=e.ondas;
  const fc=medirFC(eventos,duracion)||e.ritmo.fcV;
  const rr=fc>0?60/fc:1;
  const p2=paramsDeriv("II",e);
  const qrs=Math.round(p2.qrsNormal.dur*1000);
  const qt=Math.round((p2.qrsNormal.dur+p2.durST+p2.durT)*1000);
  const qtc=Math.round(qt/Math.sqrt(rr));
  let pr=o.PR.dur;
  if(e.ritmo.conduccion==="completo"||e.ritmo.conduccion==="independiente") pr=0;
  return {fc,pr,qrs,qt,qtc,ejeP:e.ejes.p,ejeQRS:e.ejes.qrs,ejeT:e.ejes.t,
          fcA:e.ritmo.origenAuricular==="sinusal"?e.ritmo.fcA:(e.ritmo.origenAuricular==="flutter"?e.ritmo.fcA:0)};
}

/* ============================== PRESETS ============================== */
/* discordancia: fracción del QRS que se aplica como desnivel opuesto del ST en cada derivación */
function pre(desc,cambios,discordancia){ return {desc,cambios,discordancia}; }
const PRESETS={
  "— Ritmo —":null,
  "Sinusal normal":pre("Ritmo sinusal normal, 70 lpm, ejes e intervalos dentro de límites.",{}),
  "Taquicardia sinusal":pre("Sinusal a 130 lpm; PR y QT se acortan con la frecuencia.",
    {ritmo:{fcA:130,fcV:130},ondas:{PR:{dur:130},T:{amp:0.3,dur:130,apex:0.6},ST:{nivel:0,dur:55,curva:0}}}),
  "Bradicardia sinusal":pre("Sinusal a 45 lpm; onda U más visible.",
    {ritmo:{fcA:45,fcV:45},ondas:{U:{amp:0.07,dur:110},ST:{nivel:0,dur:120,curva:0}}}),
  "Arritmia sinusal respiratoria":pre("Variación fásica del RR con la respiración; normal en jóvenes.",
    {ritmo:{fcA:72,fcV:72,variabilidad:0.16}}),
  "Fibrilación auricular con respuesta rápida":pre("Sin ondas P, ondulación de la línea de base y RR irregularmente irregular a 140 lpm.",
    {ritmo:{origenAuricular:"fibrilacion",conduccion:"independiente",fcV:140,irregularidad:0.30}}),
  "Fibrilación auricular controlada":pre("FA con respuesta ventricular de 80 lpm.",
    {ritmo:{origenAuricular:"fibrilacion",conduccion:"independiente",fcV:80,irregularidad:0.32}}),
  "Flutter auricular 2:1":pre("Ondas F en diente de sierra a 300/min, negativas en II, III y aVF; respuesta 150 lpm.",
    {ritmo:{origenAuricular:"flutter",conduccion:"independiente",fcA:300,fcV:150,irregularidad:0.02}}),
  "Flutter auricular 4:1":pre("Flutter con conducción 4:1; respuesta ventricular de 75 lpm.",
    {ritmo:{origenAuricular:"flutter",conduccion:"independiente",fcA:300,fcV:75,irregularidad:0.02}}),
  "Taquicardia supraventricular (TSV)":pre("QRS estrecho y regular a 180 lpm sin P visibles.",
    {ritmo:{origenAuricular:"ninguno",conduccion:"independiente",fcV:180,irregularidad:0},
     ondas:{ST:{nivel:-0.08,dur:45,curva:0},T:{amp:0.25,dur:110,apex:0.6}}}),
  "Ritmo de la unión":pre("Escape de la unión a 50 lpm, QRS estrecho sin P conducida.",
    {ritmo:{origenAuricular:"ninguno",conduccion:"independiente",fcV:50}}),
  "Taquicardia ventricular monomórfica":pre("QRS ancho, regular, a 180 lpm, con disociación AV.",
    {ritmo:{origenAuricular:"ninguno",conduccion:"independiente",fcV:180},
     morfologia:{qrs:"ventricular"},ondas:{ST:{nivel:0,dur:40,curva:0},T:{amp:-0.4,dur:120,apex:0.5}}}),
  "Torsades de pointes":pre("TV polimórfica con giro de puntas sobre la línea de base.",
    {ritmo:{origenAuricular:"ninguno",conduccion:"independiente",fcV:230,modo:"torsades"},
     morfologia:{qrs:"ventricular"}}),
  "Fibrilación ventricular":pre("Actividad caótica sin complejos identificables.",
    {ritmo:{modo:"fv",fcV:300,origenAuricular:"ninguno",conduccion:"independiente"}}),
  "Ritmo idioventricular acelerado":pre("QRS ancho a 70 lpm, frecuente tras la reperfusión.",
    {ritmo:{origenAuricular:"ninguno",conduccion:"independiente",fcV:70},morfologia:{qrs:"ventricular"}},0.15),
  "Asistolia":pre("Ausencia de actividad eléctrica.",{ritmo:{modo:"asistolia"}}),
  "Marcapasos ventricular":pre("Espiga seguida de QRS ancho con morfología de bloqueo de rama izquierda y repolarización discordante.",
    {ritmo:{origenAuricular:"ninguno",conduccion:"independiente",fcV:70},morfologia:{qrs:"marcapasos"}},0.15),
  "Extrasístoles ventriculares aisladas":pre("EV con QRS ancho, T discordante y pausa compensatoria.",
    {ritmo:{ectopia:"ev",cadaEctopico:5}}),
  "Bigeminismo ventricular":pre("Un latido sinusal alternando con una extrasístole ventricular.",
    {ritmo:{ectopia:"bigeminismo"}}),
  "Extrasístoles auriculares":pre("P prematura de distinta morfología con pausa no compensatoria.",
    {ritmo:{ectopia:"ea",cadaEctopico:5}}),

  "— Conducción —":null,
  "Bloqueo AV de primer grado":pre("PR de 280 ms, constante, todas las P conducen.",
    {ondas:{PR:{dur:280}},ritmo:{fcA:68,fcV:68}}),
  "Bloqueo AV 2.º Mobitz I (Wenckebach)":pre("Alargamiento progresivo del PR hasta una P bloqueada, con agrupamiento de latidos.",
    {ritmo:{conduccion:"wenckebach",relacion:4,fcA:75},ondas:{PR:{dur:170}}}),
  "Bloqueo AV 2.º Mobitz II":pre("PR constante con caída súbita de un QRS; alto riesgo de progresión.",
    {ritmo:{conduccion:"mobitz2",relacion:4,fcA:75},ondas:{PR:{dur:180}}}),
  "Bloqueo AV 2:1":pre("Una de cada dos P conduce.",
    {ritmo:{conduccion:"mobitz2",relacion:2,fcA:80},ondas:{PR:{dur:180}}}),
  "Bloqueo AV completo":pre("Disociación AV: aurículas a 90/min y escape ventricular a 35 lpm.",
    {ritmo:{conduccion:"completo",fcA:90,fcV:35,escape:"ventricular"},morfologia:{qrs:"ventricular"},
     ejes:{qrs:100}}),   /* eje del escape: con 45° el QRS sale casi isoeléctrico en DII y la tira de ritmo no lo enseña */
  "Bloqueo de rama derecha":pre("QRS ≥120 ms, rSR' en V1-V2, S ancha en I y V6, con ST y T discordantes en V1-V2.",
    {morfologia:{qrs:"brd"}},0.12),
  "Bloqueo de rama izquierda":pre("QRS ≥120 ms, QS en V1-V3 y R ancha mellada en V5-V6 sin q septal, con ST y T discordantes.",
    {morfologia:{qrs:"bri"},ondas:{T:{amp:-0.35,dur:180,apex:0.6},ST:{nivel:0,dur:80,curva:0}}},0.15),
  "Hemibloqueo anterior izquierdo":pre("Eje desviado a −60°, qR en I y aVL, rS en II, III y aVF.",
    {ejes:{qrs:-60},morfologia:{hbai:true}}),
  "Hemibloqueo posterior izquierdo":pre("Eje desviado a +120° con rS en I y qR en III.",
    {ejes:{qrs:120},morfologia:{hbpi:true}}),
  "Wolff-Parkinson-White":pre("PR corto (<120 ms), onda delta y QRS ancho por preexcitación.",
    {ondas:{PR:{dur:90},P:{dur:85},T:{amp:0.25}},morfologia:{qrs:"wpw"}},0.08),

  "— Isquemia e infarto —":null,
  "IAMCEST anteroseptal":pre("Elevación del ST en V2-V3 (con V1 y V4 menores) e imagen recíproca inferior.",
    {isquemia:{territorio:"anteroseptal",elevacion:0.45,fase:"aguda"},ondas:{ST:{nivel:0,dur:90,curva:0.6}}}),
  "IAMCEST anterior extenso":pre("Elevación de V1 a V6, I y aVL; oclusión proximal de la descendente anterior.",
    {isquemia:{territorio:"anterior_extenso",elevacion:0.50,fase:"aguda"},ondas:{ST:{nivel:0,dur:90,curva:0.6}}}),
  "IAMCEST inferior":pre("Elevación en II, III y aVF con descenso recíproco en I y aVL.",
    {isquemia:{territorio:"inferior",elevacion:0.40,fase:"aguda"},ondas:{ST:{nivel:0,dur:90,curva:0.6}},
     ejes:{qrs:70},ritmo:{fcA:58,fcV:58}}),
  "IAMCEST lateral":pre("Elevación en I, aVL, V5 y V6.",
    {isquemia:{territorio:"lateral",elevacion:0.40,fase:"aguda"},ondas:{ST:{nivel:0,dur:90,curva:0.6}}}),
  "IAMCEST sobre BRI (criterio de Sgarbossa)":pre(
    "BRI previo con elevación del ST concordante ≥1 mm en I, aVL, V5 y V6, por encima de la discordancia propia del bloqueo.",
    {morfologia:{qrs:"bri"},isquemia:{territorio:"lateral",elevacion:0.42,fase:"aguda"},
     ondas:{ST:{nivel:0,dur:90,curva:0.5},T:{amp:-0.35,dur:180,apex:0.6}},
     ritmo:{fcA:78,fcV:78},ejes:{qrs:-20}},0.15),
  /* posterior aislado (circunfleja): solo la imagen en espejo de V1-V3, sin elevación inferior.
     Se arma con desniveles y ajustes por derivación porque el territorio «posterior» del motor
     incluye la cara inferior (inferoposterior). */
  "IAM posterior":pre("IAM posterior aislado: imagen en espejo en V1-V3 con descenso horizontal del ST, R alta y ancha (R/S >1 en V2) y T positiva y alta; sin elevación inferior.",
    {ondas:{ST:{nivel:0,dur:90,curva:0}},
     stDeriv:{V1:-0.20,V2:-0.25,V3:-0.20},
     ajusteDeriv:{V1:{r:4.0,rDur:60,s:0.35,t:0.55},V2:{r:2.6,rDur:60,s:0.45,t:0.45},V3:{r:1.5,s:0.65,t:0.30}}}),
  "IAM de ventrículo derecho":pre("Elevación en V1 y III mayor que en II; acompaña al infarto inferior.",
    {isquemia:{territorio:"vd",elevacion:0.40,fase:"aguda"},ritmo:{fcA:55,fcV:55}}),
  "T hiperagudas (fase muy precoz)":pre("Ondas T altas y anchas antes de la elevación franca del ST.",
    {isquemia:{territorio:"anterior",elevacion:0.35,fase:"hiperaguda"}}),
  /* fase «cronica» del motor: Q ≥40 ms y >25 % de la R, ST prácticamente isoeléctrico (residuo 0.03 mV)
     y T negativa; el ajuste de T en II, III y aVF anula la T alta «recíproca» que no corresponde a un infarto antiguo */
  "Infarto con onda Q establecido":pre("Infarto anterior establecido: Q patológicas (≥40 ms, >25 % de la R) en V1-V4, ST isoeléctrico y T negativa simétrica en el territorio.",
    {isquemia:{territorio:"anterior",elevacion:0.30,fase:"cronica"},ondas:{T:{apex:0.5}},
     ajusteDeriv:{II:{t:-0.26},III:{t:-0.33},aVF:{t:-0.30}}}),
  "Isquemia subendocárdica":pre("Descenso horizontal del ST con T negativas, sin elevación.",
    {ondas:{ST:{nivel:-0.22,dur:100,curva:-0.2},T:{amp:-0.30,dur:180,apex:0.55}},ritmo:{fcA:95,fcV:95}}),
  "Pericarditis aguda":pre("Elevación cóncava difusa del ST con descenso del PR y elevación del ST en aVR invertida.",
    {isquemia:{territorio:"difusa",elevacion:0.22,fase:"aguda"},ondas:{ST:{nivel:0,dur:95,curva:-0.7}},
     ritmo:{fcA:100,fcV:100}}),

  "— Metabólico y otros —":null,
  "Hiperkalemia leve":pre("Ondas T picudas, altas y simétricas.",
    {ondas:{T:{amp:0.95,dur:120,apex:0.5},P:{amp:0.10,dur:100,apex:0.5}}}),
  "Hiperkalemia grave":pre("P plana o ausente, QRS ensanchado y T picudas: patrón preagónico.",
    {ondas:{T:{amp:1.20,dur:120,apex:0.5},P:{amp:0.02,dur:120,apex:0.5},R:{amp:0.9,dur:80},S:{amp:-0.45,dur:60}},
     morfologia:{qrs:"bri"},ritmo:{fcA:58,fcV:58}}),
  "Hipokalemia":pre("T aplanada, onda U prominente y descenso del ST.",
    {ondas:{T:{amp:0.09,dur:190,apex:0.6},U:{amp:0.28,dur:150},ST:{nivel:-0.12,dur:110,curva:-0.2}}}),
  "Hipocalcemia (QT largo)":pre("Segmento ST alargado con QT prolongado.",
    {ondas:{ST:{nivel:0,dur:210,curva:0},T:{amp:0.30,dur:180,apex:0.6}}}),
  "Hipercalcemia (QT corto)":pre("ST muy corto y QT acortado.",
    {ondas:{ST:{nivel:0,dur:25,curva:0},T:{amp:0.35,dur:150,apex:0.6}}}),
  "Síndrome de QT largo":pre("QT prolongado con T ancha y bífida; riesgo de torsades.",
    {ondas:{ST:{nivel:0,dur:200,curva:0},T:{amp:0.28,dur:260,apex:0.7}},ritmo:{fcA:60,fcV:60}}),
  "Hipertrofia ventricular izquierda":pre("Criterios de voltaje de Sokolow-Lyon con patrón de sobrecarga.",
    {morfologia:{perfil:"hvi"},ondas:{R:{amp:1.6,dur:55},S:{amp:-0.35,dur:40}},ejes:{qrs:20,t:210}}),
  "Hipertrofia ventricular derecha":pre("R dominante en V1 con eje derecho y sobrecarga derecha.",
    {morfologia:{perfil:"hvd"},ejes:{qrs:110,t:250}}),
  /* solo V1-V2: el resto del trazo queda normal. La curvatura del ST solo actúa donde hay desnivel. */
  "Brugada tipo 1":pre("Elevación del ST ≥2 mm en domo (coved) en V1-V2 seguida de T negativa solo en V1-V2, con r' de BRDHH incompleto; resto de derivaciones normales.",
    {ondas:{ST:{nivel:0,dur:80,curva:0.3}},
     stDeriv:{V1:0.30,V2:0.35},
     ajusteDeriv:{V1:{s:0.7,t:-0.30},V2:{s:0.7,t:-0.70}}}),
  "Miocardiopatía hipertrófica":pre("Q profundas y estrechas («en daga», <40 ms) en I, aVL, II, III, aVF, V5 y V6, voltaje de HVI y T negativas laterales.",
    {morfologia:{perfil:"hvi"},ondas:{R:{amp:1.6,dur:55},S:{amp:-0.35,dur:40}},ejes:{qrs:45,t:210},
     stDeriv:{I:-0.05,aVL:-0.05,V5:-0.08,V6:-0.08},
     ajusteDeriv:{I:{q:-0.40,qDur:28},aVL:{q:-0.30,qDur:28},II:{q:-0.45,qDur:28},III:{q:-0.35,qDur:28},
                  aVF:{q:-0.40,qDur:28},V5:{q:-0.55,qDur:28},V6:{q:-0.50,qDur:28}}}),
  "Tromboembolia pulmonar (S1Q3T3)":pre("S en I, Q y T negativa en III, taquicardia sinusal y sobrecarga derecha.",
    {ritmo:{fcA:115,fcV:115},ejes:{qrs:100,t:280},morfologia:{perfil:"hvd"},ondas:{S:{amp:-0.45,dur:45}}}),
  "Hipotermia (onda J de Osborn)":pre("Deflexión J al final del QRS con bradicardia.",
    {ritmo:{fcA:45,fcV:45},morfologia:{osborn:true}})
};

function aplicarPreset(nombre){
  const p=PRESETS[nombre];
  if(!p) return;
  E=copia(ESTADO_BASE);
  E.paciente=copia(estadoPacienteActual());
  fusionar(E,p.cambios);
  if(p.discordancia) E.stDeriv=discordanciaSecundaria(E,p.discordancia);
  E.presetActual=nombre;
  const sel=document.getElementById("sel-preset");
  if(sel.value!==nombre) sel.value=nombre;
  document.getElementById("preset-desc").textContent=p.desc;
  sincronizarUI();
  regenerar();
}
function fusionar(dest,src){
  for(const k in src){
    if(src[k]&&typeof src[k]==="object"&&!Array.isArray(src[k])){
      if(!dest[k]) dest[k]={};
      fusionar(dest[k],src[k]);
    }else dest[k]=src[k];
  }
}
function estadoPacienteActual(){ return E.paciente; }

/* ============================== ONDAS: definición de la UI ============================== */
const COLOR_ONDA={P:"#5ee6a8",PR:"#4da3ff",F:"#ffffff",Q:"#ff7b8a",R:"#ffd166",
                  S:"#a78bfa",ST:"#f59e0b",T:"#f472d0",TP:"#67e8f9",U:"#a3e635"};
const DEF_ONDAS=[
  {id:"P", campos:[["amp","Amplitud","mV",-0.5,0.6,0.01],["dur","Duración","ms",40,200,5],["apex","Posición del ápex","",0.1,0.9,0.01]]},
  {id:"PR",campos:[["dur","Intervalo PR","ms",80,400,5]]},
  {id:"Q", campos:[["amp","Amplitud","mV",-1.2,0,0.01],["dur","Duración","ms",0,80,1]]},
  {id:"R", campos:[["amp","Amplitud","mV",0,3,0.05],["dur","Duración","ms",20,140,2]]},
  {id:"S", campos:[["amp","Amplitud","mV",-2,0,0.05],["dur","Duración","ms",0,100,2]]},
  {id:"ST",campos:[["nivel","Desnivel del ST","mV",-0.5,0.8,0.01],["dur","Duración","ms",20,260,5],["curva","Convexidad","",-1,1,0.05]]},
  {id:"T", campos:[["amp","Amplitud","mV",-1,1.5,0.01],["dur","Duración","ms",80,320,5],["apex","Simetría (ápex)","",0.3,0.85,0.01]]},
  {id:"U", campos:[["amp","Amplitud","mV",-0.2,0.5,0.01],["dur","Duración","ms",40,200,5]]}
];
let ondaSel="P";

/* ============================== ANIMACIÓN ============================== */
const Anim={
  buf:{}, dur:16, fs:500, t:0, x:0, corriendo:true, ultimo:0, ultimoQRS:-1, ctxAudio:null
};
function bufferDeriv(d){
  if(!Anim.buf[d]) Anim.buf[d]=generarSenal(E,eventosAnim,d,Anim.dur,Anim.fs);
  return Anim.buf[d];
}
let eventosAnim=[];

function derivsMonitor(){
  const c=E.monitor.canales;
  if(c===1) return [E.monitor.deriv];
  if(c===3) return ["II","V1","V5"];
  return DERIVACIONES;
}
function dibujarRejilla(ctx,w,h,filas,mmPx){
  ctx.fillStyle="#05080d"; ctx.fillRect(0,0,w,h);
  ctx.lineWidth=1;
  ctx.strokeStyle="rgba(70,120,90,.20)";
  ctx.beginPath();
  for(let x=0;x<w;x+=mmPx){ ctx.moveTo(Math.round(x)+.5,0); ctx.lineTo(Math.round(x)+.5,h); }
  for(let y=0;y<h;y+=mmPx){ ctx.moveTo(0,Math.round(y)+.5); ctx.lineTo(w,Math.round(y)+.5); }
  ctx.stroke();
  ctx.strokeStyle="rgba(70,160,110,.38)";
  ctx.beginPath();
  for(let x=0;x<w;x+=mmPx*5){ ctx.moveTo(Math.round(x)+.5,0); ctx.lineTo(Math.round(x)+.5,h); }
  for(let y=0;y<h;y+=mmPx*5){ ctx.moveTo(0,Math.round(y)+.5); ctx.lineTo(w,Math.round(y)+.5); }
  ctx.stroke();
}
let ctxMon,ctxMonFondo,lienzoMon;
function prepararMonitor(){
  lienzoMon=document.getElementById("lienzo-monitor");
  const dpr=window.devicePixelRatio||1;
  const w=lienzoMon.clientWidth||600;
  const filas=derivsMonitor().length;
  const h=filas<=1?260:(filas===3?390:840);
  lienzoMon.style.height=h+"px";
  lienzoMon.width=w*dpr; lienzoMon.height=h*dpr;
  ctxMon=lienzoMon.getContext("2d");
  ctxMon.setTransform(dpr,0,0,dpr,0,0);
  Anim.w=w; Anim.h=h; Anim.filas=filas;
  Anim.mmPx=Math.max(3.2,w/(25*6.5));   // tamaño de cuadro fijo: a 50 mm/s cabe la mitad de tiempo
  Anim.pxPorSeg=Anim.mmPx*E.monitor.velocidad;
  Anim.altoFila=h/filas;
  dibujarRejilla(ctxMon,w,h,filas,Anim.mmPx);
  Anim.x=0; Anim.t=0; Anim.prev=null;
}
function bip(){
  if(!document.getElementById("chk-bip").checked) return;
  try{
    if(!Anim.ctxAudio) Anim.ctxAudio=new (window.AudioContext||window.webkitAudioContext)();
    const c=Anim.ctxAudio,o=c.createOscillator(),g=c.createGain();
    o.frequency.value=880; g.gain.value=0.05;
    o.connect(g); g.connect(c.destination); o.start();
    g.gain.exponentialRampToValueAtTime(0.0001,c.currentTime+0.09);
    o.stop(c.currentTime+0.1);
  }catch(err){}
}
function marco(ts){
  requestAnimationFrame(marco);
  if(!Anim.corriendo||!ctxMon){ Anim.ultimo=ts; return; }
  const dt=Math.min(0.05,(ts-Anim.ultimo)/1000||0); Anim.ultimo=ts;
  const derivs=derivsMonitor();
  const px=Anim.pxPorSeg*dt;
  const x0=Anim.x, x1=Anim.x+px;
  const gan=E.monitor.ganancia*Anim.mmPx;

  /* borrado por delante del barrido */
  const gap=Math.max(10,Anim.pxPorSeg*0.08);
  ctxMon.save();
  ctxMon.beginPath(); ctxMon.rect(x1,0,gap,Anim.h); ctxMon.clip();
  dibujarRejilla(ctxMon,Anim.w,Anim.h,Anim.filas,Anim.mmPx);
  ctxMon.restore();

  ctxMon.lineWidth=1.6; ctxMon.lineJoin="round"; ctxMon.lineCap="round";
  ctxMon.strokeStyle="#3ddc97";
  for(let f=0;f<derivs.length;f++){
    const buf=bufferDeriv(derivs[f]);
    const base=Anim.altoFila*(f+0.5);
    ctxMon.beginPath();
    const pasos=Math.max(1,Math.ceil(px));
    for(let i=0;i<=pasos;i++){
      const x=x0+px*i/pasos;
      const t=(Anim.t+dt*i/pasos)%Anim.dur;
      const v=buf[Math.floor(t*Anim.fs)]||0;
      const y=base-v*gan;
      if(i===0){ ctxMon.moveTo(x,y); } else ctxMon.lineTo(x,y);
    }
    ctxMon.stroke();
    if(f===0&&derivs.length>1){}
  }
  /* etiquetas de derivación */
  if(derivs.length>1){
    ctxMon.fillStyle="#9aa7b8"; ctxMon.font="600 11px sans-serif";
    for(let f=0;f<derivs.length;f++) ctxMon.fillText(derivs[f],6,Anim.altoFila*f+14);
  }

  /* bip en cada QRS */
  const tAnt=Anim.t, tNue=Anim.t+dt;
  for(const ev of eventosAnim){
    if(ev.tipo==="V"){
      const te=ev.t%Anim.dur;
      if(te>tAnt&&te<=tNue) bip();
    }
  }

  Anim.t=(Anim.t+dt)%Anim.dur;
  Anim.x=x1;
  if(Anim.x>=Anim.w) Anim.x=0;   /* el tiempo sigue corriendo: cada barrido muestra un tramo distinto */
}
requestAnimationFrame(marco);

/* ---- vista previa de la onda seleccionada ---- */
function dibujarPreviaOnda(){
  const c=document.getElementById("lienzo-onda");
  const dpr=window.devicePixelRatio||1, w=c.clientWidth||400, h=150;
  c.width=w*dpr; c.height=h*dpr;
  const g=c.getContext("2d"); g.setTransform(dpr,0,0,dpr,0,0);
  const mm=Math.max(3,w/60);
  dibujarRejilla(g,w,h,1,mm);
  const fs=500, dur=1.2;
  const evs=[{t:0.22,tipo:"P",ectopico:false},
             {t:0.22+E.ondas.PR.dur/1000,tipo:"V",morf:E.morfologia.qrs,ectopico:false}];
  const eLimpio=copia(E);
  eLimpio.artefactos={basal:0,muscular:0,red:0};
  eLimpio.ritmo.origenAuricular="sinusal";
  eLimpio.ritmo.modo="normal";
  const s=generarSenal(eLimpio,evs,E.monitor.deriv,dur,fs);
  g.strokeStyle=COLOR_ONDA[ondaSel]||"#3ddc97"; g.lineWidth=1.8; g.beginPath();
  for(let i=0;i<s.length;i++){
    const x=i/s.length*w, y=h*0.62-s[i]*mm*10;
    i?g.lineTo(x,y):g.moveTo(x,y);
  }
  g.stroke();
}

/* ============================== PAPEL (SVG en mm) ============================== */
const PAPELES={
  "carta-h":{w:279.4,h:215.9},"carta-v":{w:215.9,h:279.4},
  "a4-h":{w:297,h:210},"a4-v":{w:210,h:297}
};
function geometriaPapel(){
  const p=PAPELES[E.papel.tamano], m=6;
  const anchoUtil=p.w-m*2, altoUtil=p.h-m*2;
  const anchoEtiqueta=13;
  const anchoTrazo=anchoUtil-anchoEtiqueta;
  const durMax=anchoTrazo/E.papel.velocidad;
  const altoFila=Math.max(16,2.4*E.papel.ganancia);
  /* la interpretación se reparte en varias líneas; el encabezado crece con ellas */
  const lineasInterp=partirTexto(E.paciente.interpretacion||"",Math.max(20,Math.floor(anchoUtil/1.5)));
  const encabezado=lineasInterp.length?24+(lineasInterp.length-1)*4:20;
  const pie=6;
  /* filas que caben en la primera hoja (lleva encabezado) y en las siguientes */
  const filasPrimera=Math.max(1,Math.floor((altoUtil-encabezado-pie)/altoFila));
  const filasResto  =Math.max(1,Math.floor((altoUtil-4-pie)/altoFila));
  return {p,m,anchoUtil,altoUtil,anchoEtiqueta,anchoTrazo,durMax,altoFila,encabezado,pie,filasPrimera,filasResto,lineasInterp};
}
/* corta un texto en líneas de como máximo `maxCar` caracteres, sin partir palabras */
function partirTexto(txt,maxCar){
  const palabras=String(txt||"").trim().split(/\s+/).filter(Boolean);
  if(!palabras.length) return [];
  const lineas=[];
  let linea="";
  for(const pal of palabras){
    if(!linea) linea=pal;
    else if((linea+" "+pal).length<=maxCar) linea+=" "+pal;
    else{ lineas.push(linea); linea=pal; }
  }
  if(linea) lineas.push(linea);
  return lineas;
}
function rejillaSVG(x,y,w,h){
  let fino="",grueso="";
  for(let i=0;i<=Math.floor(w);i++){
    const gx=(x+i).toFixed(2);
    (i%5===0?(grueso+=`M${gx} ${y.toFixed(2)}V${(y+h).toFixed(2)}`):(fino+=`M${gx} ${y.toFixed(2)}V${(y+h).toFixed(2)}`));
  }
  for(let j=0;j<=Math.floor(h);j++){
    const gy=(y+j).toFixed(2);
    (j%5===0?(grueso+=`M${x.toFixed(2)} ${gy}H${(x+w).toFixed(2)}`):(fino+=`M${x.toFixed(2)} ${gy}H${(x+w).toFixed(2)}`));
  }
  return `<path d="${fino}" stroke="var(--papel-linea)" stroke-width="0.12" fill="none" opacity="0.75"/>`+
         `<path d="${grueso}" stroke="var(--papel-linea2)" stroke-width="0.28" fill="none" opacity="0.85"/>`;
}
function trazoSVG(senal,fs,x0,yBase,vel,gan,dur){
  const paso=Math.max(1,Math.round(fs/(vel*8)));   // ~8 puntos por milímetro
  let d="";
  for(let i=0;i<senal.length;i+=paso){
    const t=i/fs; if(t>dur) break;
    const x=(x0+t*vel).toFixed(2);
    const y=(yBase-clamp(senal[i],-6,6)*gan).toFixed(2);
    d+=(i===0?"M":"L")+x+" "+y;
  }
  return `<path d="${d}" fill="none" stroke="var(--tinta)" stroke-width="0.32" stroke-linejoin="round" stroke-linecap="round"/>`;
}
function pulsoCalibracion(x,yBase,gan){
  const a=gan*1.0;                     // 1 mV
  return `<path d="M${x} ${yBase}h2v${-a}h5v${a}h2" fill="none" stroke="var(--tinta)" stroke-width="0.32"/>`;
}
function esc(s){ return String(s==null?"":s).replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c])); }

/* ---- encabezado de identificación, común a los dos formatos ---- */
function encabezadoSVG(g,pa,med){
  const altoCab=g.m+g.encabezado-2;
  let s=`<rect x="0" y="0" width="${g.p.w}" height="${altoCab}" fill="#ffffff"/>`;
  s+=`<line x1="0" y1="${altoCab}" x2="${g.p.w}" y2="${altoCab}" stroke="var(--papel-linea2)" stroke-width="0.3"/>`;
  s+=`<g font-family="Helvetica,Arial,sans-serif" fill="var(--tinta)">`;
  s+=`<text x="${g.m}" y="${g.m+5}" font-size="4.2" font-weight="bold">${esc(pa.nombre||"—")}</text>`;
  s+=`<text x="${g.m}" y="${g.m+10.5}" font-size="3">Edad: ${esc(pa.edad||"—")}    Sexo: ${esc(pa.sexo||"—")}    Exp.: ${esc(pa.expediente||"—")}    ${esc(pa.fecha||"")} ${esc(pa.hora||"")}</text>`;
  s+=`<text x="${g.m}" y="${g.m+15.5}" font-size="3">FC ${med.fc} lpm   PR ${med.pr||"—"} ms   QRS ${med.qrs} ms   QT/QTc ${med.qt}/${med.qtc} ms   Ejes P ${med.ejeP}° QRS ${med.ejeQRS}° T ${med.ejeT}°</text>`;
  g.lineasInterp.forEach((ln,i)=>{
    s+=`<text x="${g.m}" y="${g.m+20.5+i*4}" font-size="3" font-style="italic">${esc(ln)}</text>`;
  });
  /* bloque derecho: cada línea se puede ocultar desde la barra del papel */
  const pp=E.papel, lineasDer=[];
  if(pp.mostrarMarca!==false)
    lineasDer.push({t:"EKG Maker",tam:3.6,neg:true,op:1});
  if(pp.mostrarAviso!==false)
    lineasDer.push({t:"SIMULACIÓN CON FINES EDUCATIVOS — NO ES UN ESTUDIO REAL",tam:2.6,neg:false,op:.75});
  if(pp.mostrarTitulo!==false)
    lineasDer.push({t:E.presetActual||"Trazo personalizado",tam:2.6,neg:false,op:.75});
  let yDer=g.m+5;
  for(const l of lineasDer){
    s+=`<text x="${g.p.w-g.m}" y="${yDer.toFixed(2)}" font-size="${l.tam}" text-anchor="end"`+
       (l.neg?` font-weight="bold"`:``)+(l.op<1?` opacity="${l.op}"`:``)+`>${esc(l.t)}</text>`;
    yDer+=l.tam>3?5:4.5;
  }
  s+=`</g>`;
  return s;
}

/* ---- escalado en pantalla y aviso de impresión, comunes a los dos formatos ---- */
function ajustarHojas(g,recorte){
  const cont=document.getElementById("hojas");
  const zoom=document.getElementById("chk-zoom").checked;
  const anchoDisp=cont.clientWidth-40;
  document.querySelectorAll("#hojas .hoja").forEach(h=>{
    if(zoom){
      const mmPx=96/25.4, k=Math.min(1,anchoDisp/(g.p.w*mmPx));
      h.style.transform=`scale(${k})`;
      h.style.marginBottom=(g.p.h*mmPx*(k-1))+"px";
    }else{ h.style.transform=""; h.style.marginBottom=""; }
  });
  const aviso=document.querySelector(".aviso");
  if(recorte){
    aviso.innerHTML=`<b>Aviso:</b> a ${E.papel.velocidad} mm/s este papel sólo admite ${g.durMax.toFixed(1)} s por fila; el trazo se recortó. Usa papel horizontal o baja la duración.`;
  }else{
    aviso.innerHTML=`En el diálogo de impresión elige <b>Escala 100 %</b> (no «Ajustar al área imprimible») y márgenes «Ninguno» para conservar la escala real de ${E.papel.velocidad} mm/s y ${E.papel.ganancia} mm/mV.`;
  }
}

function construirHojas(){
  if(E.papel.formato==="4x3") construirHoja4x3();
  else construirHojasTira();
}

/* deja la barra del papel en sintonía con el estado (p. ej. tras importar un caso) */
function sincronizarBarraPapel(){
  const v=(id,val)=>{const el=document.getElementById(id); if(el) el.value=val;};
  v("sel-papel",E.papel.tamano);
  v("sel-formato",E.papel.formato||"tira");
  v("sel-dur",String(E.papel.duracion));
  v("sel-vel-papel",String(E.papel.velocidad));
  v("sel-gan-papel",String(E.papel.ganancia));
  const marca=(id,val)=>{const el=document.getElementById(id); if(el) el.checked=val!==false;};
  marca("chk-ritmo-abajo",E.papel.tiraRitmo);
  marca("chk-marca",E.papel.mostrarMarca);
  marca("chk-aviso-educativo",E.papel.mostrarAviso);
  marca("chk-titulo-caso",E.papel.mostrarTitulo);
  /* en el formato 4×3 la duración la fija el ancho de la hoja */
  const gd=document.getElementById("grupo-duracion");
  if(gd) gd.style.display=(E.papel.formato==="4x3")?"none":"";
}

/* ---- formato clásico de 12 derivaciones en una hoja: 4 columnas × 3 filas + tira de ritmo ---- */
function construirHoja4x3(){
  const g=geometriaPapel(), pa=E.paciente;
  const fs=500;
  const COLUMNAS=[["I","II","III"],["aVR","aVL","aVF"],["V1","V2","V3"],["V4","V5","V6"]];
  const anchoCol=g.anchoTrazo/4;
  const segDur=anchoCol/E.papel.velocidad;      // segundos por columna
  const durTotal=segDur*4;
  const med=mediciones(E,eventosPapel,durTotal);
  const nFilas=3+(E.papel.tiraRitmo?1:0);

  const yTrazos=g.m+g.encabezado;
  const disponible=g.altoUtil-g.encabezado-g.pie;
  const altoFila=disponible/nFilas;
  const x0=g.m+g.anchoEtiqueta;

  const cont=document.getElementById("hojas");
  cont.innerHTML="";
  document.getElementById("estilo-pagina").textContent=
    `@page{size:${g.p.w}mm ${g.p.h}mm;margin:0}`;

  let svg=`<svg class="hoja" xmlns="http://www.w3.org/2000/svg" width="${g.p.w}mm" height="${g.p.h}mm" viewBox="0 0 ${g.p.w} ${g.p.h}">`;
  svg+=`<rect width="${g.p.w}" height="${g.p.h}" fill="var(--papel)"/>`;
  svg+=encabezadoSVG(g,pa,med);
  svg+=rejillaSVG(g.m,yTrazos,g.anchoUtil,Math.min(nFilas*altoFila,g.p.h-yTrazos-g.m));

  /* tres filas de cuatro derivaciones, cada una con su tramo del mismo registro continuo */
  for(let f=0;f<3;f++){
    const yBase=yTrazos+altoFila*(f+0.55);
    svg+=pulsoCalibracion(g.m+3,yBase,E.papel.ganancia);
    for(let c=0;c<4;c++){
      const d=COLUMNAS[c][f];
      const senal=senalPapel(d,durTotal,fs);
      const i0=Math.round(c*segDur*fs), i1=Math.round((c+1)*segDur*fs);
      const x=x0+c*anchoCol;
      if(c>0) svg+=`<path d="M${x.toFixed(2)} ${(yBase-7).toFixed(2)}v14" stroke="var(--tinta)" stroke-width="0.22" opacity=".6" fill="none"/>`;
      svg+=trazoSVG(senal.slice(i0,i1),fs,x,yBase,E.papel.velocidad,E.papel.ganancia,segDur);
      svg+=`<text x="${(x+1).toFixed(2)}" y="${(yBase+6).toFixed(2)}" font-size="3.4" font-weight="bold" font-family="Helvetica,Arial,sans-serif" fill="var(--tinta)">${d}</text>`;
    }
  }

  /* tira de ritmo: DII completa, del ancho de la hoja */
  if(E.papel.tiraRitmo){
    const yBase=yTrazos+altoFila*3.55;
    svg+=pulsoCalibracion(g.m+3,yBase,E.papel.ganancia);
    svg+=trazoSVG(senalPapel("II",durTotal,fs),fs,x0,yBase,E.papel.velocidad,E.papel.ganancia,durTotal);
    svg+=`<text x="${(x0+1).toFixed(2)}" y="${(yBase+6).toFixed(2)}" font-size="3.4" font-weight="bold" font-family="Helvetica,Arial,sans-serif" fill="var(--tinta)">II</text>`;
  }

  svg+=`<text x="${g.m}" y="${g.p.h-g.m+1}" font-size="2.8" font-family="Helvetica,Arial,sans-serif" fill="var(--tinta)">${E.papel.velocidad} mm/s    ${E.papel.ganancia} mm/mV    0.05–150 Hz    4 de ${segDur.toFixed(1)} s${E.papel.tiraRitmo?" + 1 deriv. de ritmo":""}    Hoja 1 de 1</text>`;
  svg+=`</svg>`;
  cont.insertAdjacentHTML("beforeend",svg);

  ajustarHojas(g,false);
}

function construirHojasTira(){
  const g=geometriaPapel(), pa=E.paciente, med=mediciones(E,eventosPapel,E.papel.duracion);
  const dur=Math.min(E.papel.duracion,g.durMax);
  const fs=500;
  const filas=DERIVACIONES.slice();
  if(E.papel.tiraRitmo) filas.push("II*");
  const paginas=[];
  for(let i=0;i<filas.length;){
    const n=paginas.length===0?g.filasPrimera:g.filasResto;
    paginas.push(filas.slice(i,i+n)); i+=n;
  }

  const cont=document.getElementById("hojas");
  cont.innerHTML="";
  document.getElementById("estilo-pagina").textContent=
    `@page{size:${g.p.w}mm ${g.p.h}mm;margin:0}`;

  paginas.forEach((leads,ip)=>{
    const primera=ip===0;
    const yTrazos=g.m+(primera?g.encabezado:4);
    /* reparte el espacio sobrante entre las filas de la hoja, sin separarlas en exceso */
    const disponible=g.altoUtil-(primera?g.encabezado:4)-g.pie;
    const altoFila=Math.min(g.altoFila*1.5,disponible/leads.length);
    const altoRejilla=leads.length*altoFila;
    let svg=`<svg class="hoja" xmlns="http://www.w3.org/2000/svg" width="${g.p.w}mm" height="${g.p.h}mm" viewBox="0 0 ${g.p.w} ${g.p.h}">`;
    svg+=`<rect width="${g.p.w}" height="${g.p.h}" fill="var(--papel)"/>`;

    if(primera) svg+=encabezadoSVG(g,pa,med);

    svg+=rejillaSVG(g.m,yTrazos,g.anchoUtil,Math.min(altoRejilla,g.p.h-yTrazos-g.m));

    leads.forEach((lead,i)=>{
      const d=lead==="II*"?"II":lead;
      const yBase=yTrazos+altoFila*(i+0.55);
      const senal=senalPapel(d,dur,fs);
      svg+=`<text x="${g.m+0.8}" y="${yBase-1}" font-size="3.6" font-weight="bold" font-family="Helvetica,Arial,sans-serif" fill="var(--tinta)">${lead==="II*"?"II":lead}</text>`;
      svg+=pulsoCalibracion(g.m+6.5,yBase,E.papel.ganancia);
      svg+=trazoSVG(senal,fs,g.m+g.anchoEtiqueta,yBase,E.papel.velocidad,E.papel.ganancia,dur);
    });

    svg+=`<text x="${g.m}" y="${g.p.h-g.m+1}" font-size="2.8" font-family="Helvetica,Arial,sans-serif" fill="var(--tinta)">${E.papel.velocidad} mm/s    ${E.papel.ganancia} mm/mV    0.05–150 Hz    Hoja ${ip+1} de ${paginas.length}</text>`;
    svg+=`</svg>`;
    cont.insertAdjacentHTML("beforeend",svg);
  });

  ajustarHojas(g,E.papel.duracion>g.durMax+0.01);
}
let eventosPapel=[], cachePapel={};
function senalPapel(d,dur,fs){
  const k=d+"|"+dur+"|"+fs;
  if(!cachePapel[k]) cachePapel[k]=generarSenal(E,eventosPapel,d,dur,fs);
  return cachePapel[k];
}

/* ---- exportación PNG a 300 ppp ---- */
function exportarPNG(){
  const svgs=document.querySelectorAll("#hojas .hoja");
  if(!svgs.length) return;
  const g=geometriaPapel();
  const ppp=300, kf=ppp/25.4;
  svgs.forEach((svg,i)=>{
    const clon=svg.cloneNode(true);
    /* las variables CSS no viajan dentro del SVG serializado: se fijan literales */
    let txt=new XMLSerializer().serializeToString(clon)
      .replace(/var\(--papel\)/g,"#f7cfc6").replace(/var\(--papel-linea2\)/g,"#c8624f")
      .replace(/var\(--papel-linea\)/g,"#e08a7d").replace(/var\(--tinta\)/g,"#12100f");
    const img=new Image();
    const url=URL.createObjectURL(new Blob([txt],{type:"image/svg+xml;charset=utf-8"}));
    img.onload=()=>{
      const c=document.createElement("canvas");
      c.width=Math.round(g.p.w*kf); c.height=Math.round(g.p.h*kf);
      const cx=c.getContext("2d");
      cx.fillStyle="#f7cfc6"; cx.fillRect(0,0,c.width,c.height);
      cx.drawImage(img,0,0,c.width,c.height);
      URL.revokeObjectURL(url);
      c.toBlob(b=>{
        const a=document.createElement("a");
        a.href=URL.createObjectURL(b);
        a.download=`ECG_${(E.paciente.expediente||"caso").replace(/\s+/g,"_")}_hoja${i+1}.png`;
        a.click();
        setTimeout(()=>URL.revokeObjectURL(a.href),4000);
      },"image/png");
    };
    img.src=url;
  });
}

/* ============================== REGENERACIÓN ============================== */
function regenerar(){
  eventosAnim=generarEventos(E,Anim.dur);
  eventosPapel=generarEventos(E,Math.max(12,E.papel.duracion));
  Anim.buf={}; cachePapel={};
  prepararMonitor();
  dibujarPreviaOnda();
  pintarMediciones();
  if(document.getElementById("vista-papel").classList.contains("abierta")) construirHojas();
}

function pintarMediciones(){
  const m=mediciones(E,eventosAnim,Anim.dur);
  document.getElementById("fc-monitor").textContent=m.fc||"—";
  const al=(v,a,b)=>(v<a||v>b)?" alerta":"";
  const cel=[
    ["FC",m.fc+"",""],
    ["PR",(m.pr?m.pr:"—")+"",m.pr?al(m.pr,120,200):""],
    ["QRS",m.qrs+"",al(m.qrs,60,110)],
    ["QT",m.qt+"",""],
    ["QTc",m.qtc+"",al(m.qtc,0,450)],
    ["Eje QRS",m.ejeQRS+"°",al(m.ejeQRS,-30,100)],
    ["Eje P",m.ejeP+"°",""],
    ["FC auric.",(m.fcA||"—")+"",""]
  ];
  document.getElementById("medidas").innerHTML=cel.map(c=>
    `<div class="medida${c[2]}"><b>${c[1]}</b><small>${c[0]}</small></div>`).join("");
}

/* ============================== UI ============================== */
/* constructor genérico de sliders con enlace directo */
function campoRango(obj,clave,etiqueta,unidad,min,max,paso,dec,cb){
  const id="r"+(campoRango.n=(campoRango.n||0)+1);
  const v=obj[clave];
  const div=document.createElement("div");
  div.className="campo";
  div.innerHTML=`<label>${etiqueta}<b id="${id}v">${dec!=null?(+v).toFixed(dec):v}${unidad?" "+unidad:""}</b></label>
    <input type="range" id="${id}" min="${min}" max="${max}" step="${paso}" value="${v}">`;
  setTimeout(()=>{
    const el=document.getElementById(id);
    if(!el) return;
    el.addEventListener("input",()=>{
      const nv=parseFloat(el.value); obj[clave]=nv;
      document.getElementById(id+"v").textContent=(dec!=null?nv.toFixed(dec):nv)+(unidad?" "+unidad:"");
      cb&&cb();
    });
  },0);
  return div;
}
function selectCampo(etiqueta,opciones,valor,cb){
  const div=document.createElement("div");
  div.className="campo";
  const id="s"+(selectCampo.n=(selectCampo.n||0)+1);
  div.innerHTML=`<label>${etiqueta}</label><select id="${id}">`+
    opciones.map(o=>`<option value="${o[0]}"${o[0]===valor?" selected":""}>${o[1]}</option>`).join("")+
    `</select>`;
  setTimeout(()=>{
    const el=document.getElementById(id);
    if(el) el.addEventListener("change",()=>cb(el.value));
  },0);
  return div;
}
function checkCampo(etiqueta,valor,cb){
  const lab=document.createElement("label");
  lab.className="check";
  lab.innerHTML=`<input type="checkbox"${valor?" checked":""}> ${etiqueta}`;
  lab.querySelector("input").addEventListener("change",e=>cb(e.target.checked));
  return lab;
}

/* ---- lista de ondas ---- */
function pintarOndas(){
  const chips=document.getElementById("chips");
  chips.innerHTML=DEF_ONDAS.map(d=>d.id).map(k=>
    `<div class="chip${k===ondaSel?" on":""}" data-onda="${k}" style="background:${COLOR_ONDA[k]}">${k}</div>`).join("");
  chips.querySelectorAll(".chip").forEach(c=>c.addEventListener("click",()=>{
    ondaSel=c.dataset.onda; pintarOndas(); pintarControlesOnda();
  }));
  const lista=document.getElementById("lista-ondas");
  lista.innerHTML="";
  for(const d of DEF_ONDAS){
    const o=E.ondas[d.id];
    const amp=("amp" in o)?o.amp.toFixed(2)+" mV":(("nivel" in o)?o.nivel.toFixed(2)+" mV":"—");
    const b=document.createElement("button");
    b.className="onda"+(d.id===ondaSel?" sel":"");
    b.innerHTML=`<span class="bola" style="background:${COLOR_ONDA[d.id]}">${d.id}</span>
      <span class="v1">${amp}</span><span class="v2">${o.dur} ms</span>`;
    b.addEventListener("click",()=>{ondaSel=d.id;pintarOndas();pintarControlesOnda();});
    lista.appendChild(b);
  }
}
function pintarControlesOnda(){
  const d=DEF_ONDAS.find(x=>x.id===ondaSel)||DEF_ONDAS[0];
  document.getElementById("titulo-onda").textContent="Controles de la onda "+d.id;
  const cont=document.getElementById("controles-onda");
  cont.innerHTML="";
  const o=E.ondas[d.id];
  for(const c of d.campos){
    const dec=(c[2]==="mV"||c[1].indexOf("ápex")>=0||c[1].indexOf("Simetría")>=0||c[1]==="Convexidad")?2:0;
    cont.appendChild(campoRango(o,c[0],c[1],c[2],c[3],c[4],c[5],dec,()=>{pintarOndas();regenerar();}));
  }
}

/* ---- controles de ritmo ---- */
function pintarRitmo(){
  const c=document.getElementById("controles-ritmo"); c.innerHTML="";
  const r=E.ritmo;
  c.appendChild(selectCampo("Actividad auricular",[
    ["sinusal","Ondas P sinusales"],["flutter","Flutter (ondas F)"],
    ["fibrilacion","Fibrilación (ondas f)"],["ninguno","Sin actividad auricular"]],
    r.origenAuricular,v=>{r.origenAuricular=v;pintarRitmo();regenerar();}));
  if(r.origenAuricular==="sinusal"||r.origenAuricular==="flutter"){
    c.appendChild(campoRango(r,"fcA",r.origenAuricular==="flutter"?"Frecuencia auricular (flutter)":"Frecuencia auricular","/min",
      r.origenAuricular==="flutter"?200:25, r.origenAuricular==="flutter"?400:220,1,0,()=>regenerar()));
  }
  c.appendChild(selectCampo("Conducción AV",[
    ["normal","Normal (1:1)"],["wenckebach","Mobitz I (Wenckebach)"],
    ["mobitz2","Mobitz II / relación fija"],["completo","Bloqueo AV completo"],
    ["independiente","Ventrículo independiente"]],
    r.conduccion,v=>{r.conduccion=v;pintarRitmo();regenerar();}));
  if(r.conduccion==="wenckebach"||r.conduccion==="mobitz2")
    c.appendChild(campoRango(r,"relacion","Relación de conducción (N:N−1)","",2,6,1,0,()=>regenerar()));
  if(r.conduccion!=="normal"&&r.conduccion!=="wenckebach"&&r.conduccion!=="mobitz2"){
    c.appendChild(campoRango(r,"fcV","Frecuencia ventricular","lpm",20,260,1,0,()=>regenerar()));
    c.appendChild(campoRango(r,"irregularidad","Irregularidad del RR","",0,0.5,0.01,2,()=>regenerar()));
  }
  if(r.conduccion==="completo")
    c.appendChild(selectCampo("Foco de escape",[["union","De la unión (QRS estrecho)"],["ventricular","Ventricular (QRS ancho)"]],
      r.escape,v=>{r.escape=v;regenerar();}));
  c.appendChild(campoRango(r,"variabilidad","Variabilidad sinusal","",0,0.3,0.01,2,()=>regenerar()));
  c.appendChild(selectCampo("Extrasístoles",[
    ["ninguna","Ninguna"],["ev","Ventriculares aisladas"],["bigeminismo","Bigeminismo ventricular"],
    ["trigeminismo","Trigeminismo ventricular"],["dupletas","Dupletas ventriculares"],["ea","Auriculares"]],
    r.ectopia,v=>{r.ectopia=v;pintarRitmo();regenerar();}));
  if(r.ectopia==="ev"||r.ectopia==="ea")
    c.appendChild(campoRango(r,"cadaEctopico","Una cada N latidos","",3,15,1,0,()=>regenerar()));
  c.appendChild(selectCampo("Modo especial",[["normal","Ninguno"],["torsades","Torsades de pointes"],
    ["fv","Fibrilación ventricular"],["asistolia","Asistolia"]],r.modo,v=>{r.modo=v;regenerar();}));
}

/* ---- ejes y morfología ---- */
function pintarEjes(){
  const c=document.getElementById("controles-ejes"); c.innerHTML="";
  c.appendChild(campoRango(E.ejes,"qrs","Eje del QRS","°",-120,180,5,0,()=>regenerar()));
  c.appendChild(campoRango(E.ejes,"p","Eje de la P","°",-60,120,5,0,()=>regenerar()));
  c.appendChild(campoRango(E.ejes,"t","Eje de la T","°",-90,300,5,0,()=>regenerar()));
  c.appendChild(selectCampo("Morfología del QRS",[
    ["normal","Normal"],["brd","Bloqueo de rama derecha"],["bri","Bloqueo de rama izquierda"],
    ["ventricular","Origen ventricular (ancho)"],["wpw","Preexcitación (WPW)"],["marcapasos","Marcapasos"]],
    E.morfologia.qrs,v=>{E.morfologia.qrs=v;regenerar();}));
  c.appendChild(selectCampo("Progresión de R en precordiales",[
    ["normal","Normal"],["mala","Mala progresión de R"],["hvi","Hipertrofia izquierda"],["hvd","Hipertrofia derecha"]],
    E.morfologia.perfil,v=>{E.morfologia.perfil=v;regenerar();}));
  c.appendChild(checkCampo("Onda J de Osborn (hipotermia)",E.morfologia.osborn,v=>{E.morfologia.osborn=v;regenerar();}));
}

/* ---- isquemia ---- */
function pintarIsquemia(){
  const c=document.getElementById("controles-isquemia"); c.innerHTML="";
  c.appendChild(selectCampo("Territorio",[
    ["ninguno","Ninguno"],["septal","Septal (V1-V2)"],["anteroseptal","Anteroseptal (V2-V3)"],["anterior","Anterior (V1-V4)"],
    ["anterior_extenso","Anterior extenso"],["lateral","Lateral (I, aVL, V5-V6)"],
    ["inferior","Inferior (II, III, aVF)"],["posterior","Posterior (espejo V1-V3)"],
    ["vd","Ventrículo derecho"],["difusa","Difusa (pericarditis)"]],
    E.isquemia.territorio,v=>{E.isquemia.territorio=v;regenerar();}));
  c.appendChild(campoRango(E.isquemia,"elevacion","Magnitud del desnivel","mV",0,1,0.05,2,()=>regenerar()));
  c.appendChild(selectCampo("Fase evolutiva",[
    ["hiperaguda","T hiperagudas"],["aguda","Elevación del ST"],
    ["q","Onda Q con T negativa"],["cronica","Crónica (Q y T invertida)"]],
    E.isquemia.fase,v=>{E.isquemia.fase=v;regenerar();}));
  const n=document.createElement("div");
  n.className="nota";
  n.textContent="Los cambios recíprocos en las derivaciones opuestas se calculan automáticamente.";
  c.appendChild(n);
  c.appendChild(bloqueStPorDerivacion());
}

/* ---- desnivel del ST propio de cada derivación ---- */
function bloqueStPorDerivacion(){
  const det=document.createElement("details");
  det.className="st-deriv";
  det.open=DERIVACIONES.some(d=>E.stDeriv[d]);
  det.innerHTML=`<summary>Desnivel del ST por derivación</summary>`;

  const nota=document.createElement("div");
  nota.className="nota";
  nota.textContent="Se suma al nivel global de la onda ST y al territorio isquémico. 0.10 mV = 1 mm. "+
                   "Úsalo para la discordancia secundaria del BRI, del marcapasos o del ritmo ventricular.";
  det.appendChild(nota);

  const rej=document.createElement("div");
  rej.className="rejilla-st";
  for(const d of DERIVACIONES){
    rej.appendChild(campoRango(E.stDeriv,d,"ST "+d,"mV",-0.5,0.5,0.025,3,()=>regenerar()));
  }
  det.appendChild(rej);

  const fila=document.createElement("div");
  fila.className="fila-botones";

  const bAuto=document.createElement("button");
  bAuto.className="btn";
  bAuto.textContent="Discordancia automática";
  bAuto.title="Calcula el desnivel opuesto al QRS de cada derivación (15 % de su amplitud)";
  bAuto.addEventListener("click",()=>{
    if(!MORFOLOGIA_ANCHA.includes(E.morfologia.qrs)){
      alert("La discordancia secundaria solo aplica con QRS ancho.\n\n"+
            "Elige antes una morfología de bloqueo de rama, marcapasos u origen ventricular.");
      return;
    }
    E.stDeriv=discordanciaSecundaria(E,0.15);
    pintarIsquemia(); regenerar();
  });
  fila.appendChild(bAuto);

  const bCero=document.createElement("button");
  bCero.className="btn";
  bCero.textContent="Poner todo a cero";
  bCero.addEventListener("click",()=>{
    for(const d of DERIVACIONES) E.stDeriv[d]=0;
    pintarIsquemia(); regenerar();
  });
  fila.appendChild(bCero);

  det.appendChild(fila);
  return det;
}

/* ---- artefactos ---- */
function pintarArtefactos(){
  const c=document.getElementById("controles-artefactos"); c.innerHTML="";
  c.appendChild(campoRango(E.artefactos,"basal","Vagabundeo de la línea base","",0,1,0.05,2,()=>regenerar()));
  c.appendChild(campoRango(E.artefactos,"muscular","Temblor muscular","",0,1,0.05,2,()=>regenerar()));
  c.appendChild(campoRango(E.artefactos,"red","Interferencia de 60 Hz","",0,1,0.05,2,()=>regenerar()));
}

/* ---- datos del paciente ---- */
function pintarPaciente(){
  const f=document.getElementById("form-paciente");
  const campos=[["nombre","Nombre",1],["expediente","Expediente",0],["edad","Edad",0],
                ["sexo","Sexo",0],["fecha","Fecha",0],["hora","Hora",0],["interpretacion","Interpretación / diagnóstico",1]];
  f.innerHTML=campos.map(c=>
    `<div class="${c[2]?"ancho":""}"><label>${c[1]}</label>
     <input type="text" data-p="${c[0]}" value="${esc(E.paciente[c[0]])}"></div>`).join("");
  f.querySelectorAll("input").forEach(i=>i.addEventListener("input",()=>{
    E.paciente[i.dataset.p]=i.value;
    if(document.getElementById("vista-papel").classList.contains("abierta")) construirHojas();
  }));
}

/* ---- presets ---- */
function pintarPresets(){
  const s=document.getElementById("sel-preset");
  let html="",abierto=false;
  for(const k of Object.keys(PRESETS)){
    if(PRESETS[k]===null){
      if(abierto) html+="</optgroup>";
      html+=`<optgroup label="${esc(k.replace(/—/g,"").trim())}">`; abierto=true;
    }else html+=`<option value="${esc(k)}">${esc(k)}</option>`;
  }
  if(abierto) html+="</optgroup>";
  s.innerHTML=html;
  s.value="Sinusal normal";
  s.addEventListener("change",()=>aplicarPreset(s.value));
  document.getElementById("preset-desc").textContent=PRESETS["Sinusal normal"].desc;
}

/* ---- derivación del monitor ---- */
function pintarSelDeriv(){
  const s=document.getElementById("sel-deriv");
  s.innerHTML=DERIVACIONES.map(d=>`<option${d===E.monitor.deriv?" selected":""}>${d}</option>`).join("");
  s.addEventListener("change",()=>{E.monitor.deriv=s.value;Anim.buf={};prepararMonitor();dibujarPreviaOnda();});
}

function sincronizarUI(){
  pintarOndas(); pintarControlesOnda(); pintarRitmo(); pintarEjes();
  pintarIsquemia(); pintarArtefactos(); pintarPaciente();
}

/* ---- casos guardados ---- */
const LLAVE="ekgmaker.casos";
function leerCasos(){ try{return JSON.parse(localStorage.getItem(LLAVE)||"{}");}catch(e){return {};} }
function escribirCasos(c){
  try{ localStorage.setItem(LLAVE,JSON.stringify(c)); return true; }
  catch(e){ alert("No se pudieron guardar los casos en este navegador. Usa «Exportar .json»."); return false; }
}
function pintarCasos(){
  const casos=leerCasos(), cont=document.getElementById("lista-casos");
  const nombres=Object.keys(casos);
  if(!nombres.length){ cont.innerHTML=`<div class="nota">Aún no hay casos guardados.</div>`; return; }
  cont.innerHTML="";
  for(const n of nombres){
    const d=document.createElement("div");
    d.className="caso";
    d.innerHTML=`<span title="${esc(n)}">${esc(n)}</span><button data-a="c">Cargar</button><button class="del" data-a="b">✕</button>`;
    d.querySelector('[data-a="c"]').addEventListener("click",()=>{
      E=copia(casos[n]); sincronizarUI(); regenerar();
      document.getElementById("preset-desc").textContent="Caso guardado: "+n;
    });
    d.querySelector('[data-a="b"]').addEventListener("click",()=>{
      const c=leerCasos(); delete c[n]; escribirCasos(c); pintarCasos();
    });
    cont.appendChild(d);
  }
}

/* ============================== EVENTOS DE LA INTERFAZ ============================== */
function seg(idSeg,cb){
  const c=document.getElementById(idSeg);
  c.addEventListener("click",e=>{
    const b=e.target.closest("button"); if(!b) return;
    c.querySelectorAll("button").forEach(x=>x.classList.remove("on"));
    b.classList.add("on"); cb(b.dataset);
  });
}
function iniciar(){
  const ahora=new Date();
  E.paciente.fecha=ahora.toLocaleDateString("es-MX");
  E.paciente.hora=ahora.toLocaleTimeString("es-MX",{hour:"2-digit",minute:"2-digit"});

  pintarPresets(); pintarSelDeriv(); sincronizarUI(); pintarCasos(); regenerar();

  seg("seg-canales",d=>{E.monitor.canales=+d.c;Anim.buf={};prepararMonitor();});
  seg("seg-vel",d=>{E.monitor.velocidad=+d.v;prepararMonitor();});
  seg("seg-gan",d=>{E.monitor.ganancia=+d.g;prepararMonitor();});
  document.getElementById("btn-pausa").addEventListener("click",e=>{
    Anim.corriendo=!Anim.corriendo; e.target.textContent=Anim.corriendo?"Pausar":"Reanudar";
  });

  /* papel */
  const vista=document.getElementById("vista-papel");
  document.getElementById("btn-papel").addEventListener("click",()=>{
    vista.classList.add("abierta"); sincronizarBarraPapel(); construirHojas();
  });
  document.getElementById("btn-cerrar-papel").addEventListener("click",()=>vista.classList.remove("abierta"));
  document.getElementById("sel-papel").addEventListener("change",e=>{E.papel.tamano=e.target.value;construirHojas();});
  document.getElementById("sel-formato").addEventListener("change",e=>{
    E.papel.formato=e.target.value; sincronizarBarraPapel(); construirHojas();
  });
  document.getElementById("sel-dur").addEventListener("change",e=>{
    E.papel.duracion=+e.target.value; eventosPapel=generarEventos(E,Math.max(12,E.papel.duracion));
    cachePapel={}; construirHojas();
  });
  document.getElementById("sel-vel-papel").addEventListener("change",e=>{E.papel.velocidad=+e.target.value;construirHojas();});
  document.getElementById("sel-gan-papel").addEventListener("change",e=>{E.papel.ganancia=+e.target.value;construirHojas();});
  document.getElementById("chk-ritmo-abajo").addEventListener("change",e=>{E.papel.tiraRitmo=e.target.checked;construirHojas();});
  document.getElementById("chk-marca").addEventListener("change",e=>{E.papel.mostrarMarca=e.target.checked;construirHojas();});
  document.getElementById("chk-aviso-educativo").addEventListener("change",e=>{E.papel.mostrarAviso=e.target.checked;construirHojas();});
  document.getElementById("chk-titulo-caso").addEventListener("change",e=>{E.papel.mostrarTitulo=e.target.checked;construirHojas();});
  document.getElementById("chk-zoom").addEventListener("change",construirHojas);
  document.getElementById("btn-imprimir").addEventListener("click",()=>window.print());
  document.getElementById("btn-png").addEventListener("click",exportarPNG);

  /* guardar / exportar */
  function guardarCaso(){
    const campo=document.getElementById("nombre-caso");
    const n=(campo.value||"").trim()||
      (E.paciente.nombre+" — "+(E.presetActual||"personalizado"));
    const c=leerCasos(); c[n]=copia(E);
    if(escribirCasos(c)) campo.value="";
    pintarCasos();
  }
  document.getElementById("btn-guardar").addEventListener("click",guardarCaso);
  document.getElementById("btn-guardar-caso").addEventListener("click",guardarCaso);
  document.getElementById("nombre-caso").addEventListener("keydown",e=>{if(e.key==="Enter")guardarCaso();});
  document.getElementById("btn-json").addEventListener("click",()=>{
    const b=new Blob([JSON.stringify(E,null,2)],{type:"application/json"});
    const a=document.createElement("a");
    a.href=URL.createObjectURL(b);
    a.download="ecg_"+(E.paciente.expediente||"caso").replace(/\s+/g,"_")+".json";
    a.click(); setTimeout(()=>URL.revokeObjectURL(a.href),4000);
  });
  document.getElementById("btn-importar").addEventListener("click",()=>document.getElementById("archivo-json").click());
  document.getElementById("archivo-json").addEventListener("change",e=>{
    const ent=e.target;
    const f=ent.files[0];
    ent.value="";                       // permite reimportar el mismo archivo
    if(!f) return;
    const r=new FileReader();
    r.onerror=()=>alert("No se pudo leer el archivo «"+f.name+"».");
    r.onload=()=>{
      let datos;
      try{
        datos=JSON.parse(r.result);
      }catch(err){
        alert("El archivo no es un JSON válido.\n\n"+err.message);
        return;
      }
      if(!datos||typeof datos!=="object"||Array.isArray(datos)){
        alert("El archivo no contiene un caso de ECG (se esperaba un objeto JSON).");
        return;
      }
      const previo=copia(E);
      try{
        E=fusionarEstado(ESTADO_BASE,datos); // completa las claves faltantes
        sincronizarUI(); regenerar();
        document.getElementById("preset-desc").textContent=
          datos.presetActual||("Caso importado: "+f.name);
      }catch(err){
        E=previo; sincronizarUI(); regenerar();
        alert("Error al aplicar el caso importado; se conservó el trazo anterior.\n\n"+err.message);
      }
    };
    r.readAsText(f);
  });

  window.addEventListener("resize",()=>{prepararMonitor();dibujarPreviaOnda();});
}
document.addEventListener("DOMContentLoaded",iniciar);