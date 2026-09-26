"use client";

import { useState } from "react";
import { ResultView } from "../../components/result-view";
import { SequenceBox } from "../../components/sequence-box";
import { Sidebar } from "../../components/sidebar";
import { ThemeToggle } from "../../components/theme-toggle";

type Locale = "en" | "pt" | "es";
type Tool = "home" | "dnaRna" | "dnaProtein" | "rnaProtein" | "rnaDna" | "protein" | "mutation";
type Options = { full: boolean; count: boolean; pi: boolean; charge: boolean; aromaticity: boolean; secondary: boolean; weight: boolean; hydrophobic: boolean; composition: boolean; extinction: boolean };

const COPY = {
  en: { brand:"BioSTAR", eyebrow:"Bioinformatics Software for Targeted Analysis and Research", start:"Start", gene:"Gene Section", proteinSection:"Protein Section", experimental:"Experimental Section", help:"Help Section", results:"Results", compare:"Compare", save:"Save JSON", upload:"Upload", optimize:"Optimize Codons", howTo:"How To Use", github:"Visit on GitHub", title:"Targeted biological sequence analysis", subtitle:"Analyze DNA, RNA and proteins in one focused workspace.", sequenceType:"Sequence type", dna:"DNA", rna:"RNA", protein:"Protein", peptideThreshold:"Peptide size threshold", sequence:"Sequence", reference:"Reference sequence", analyze:"Analyze", mutation:"Mutation compare", otherTools:"Other Tools", codonTable:"Codon Table", fasta:"FASTA Format Tool", started:"How To Get Started?", viewer:"3D Protein Viewer", sampleDna:"Load DNA Sample", sampleProtein:"Load Protein Sample", options:"Analysis options", full:"Full test suite", count:"Amino-acid count", pi:"Isoelectric point", charge:"Charge at pH", aromaticity:"Aromaticity", secondary:"Secondary structure", weight:"Molecular weight", hydrophobic:"Hydrophobic index", composition:"Composition ratio", extinction:"Extinction coefficient", pH:"pH", empty:"Enter a sequence to begin.", limit:"GET conversions: 1,000 nt · POST analysis: 10,000 nt", copy:"Copy", copied:"Copied", saved:"JSON saved", mutationHint:"Compare two coding DNA sequences of equal length.", footer:"Developed for BioSTAR", local:"Local-first · API powered", coming:"Coming soon" },
  pt: { brand:"BioSTAR", eyebrow:"Software de Bioinformática para Análise e Pesquisa Direcionada", start:"Início", gene:"Seção de Genes", proteinSection:"Seção de Proteínas", experimental:"Seção Experimental", help:"Ajuda", results:"Resultados", compare:"Comparar", save:"Salvar JSON", upload:"Enviar", optimize:"Otimizar Códons", howTo:"Como Usar", github:"Visitar no GitHub", title:"Análise direcionada de sequências biológicas", subtitle:"Analise DNA, RNA e proteínas em um único espaço de trabalho.", sequenceType:"Tipo de sequência", dna:"DNA", rna:"RNA", protein:"Proteína", peptideThreshold:"Limite de tamanho do peptídeo", sequence:"Sequência", reference:"Sequência de referência", analyze:"Analisar", mutation:"Comparar mutações", otherTools:"Outras Ferramentas", codonTable:"Tabela de Códons", fasta:"Formatador FASTA", started:"Como Começar?", viewer:"Visualizador 3D de Proteínas", options:"Opções de análise", full:"Suite completa", count:"Contagem de aminoácidos", pi:"Ponto isoelétrico", charge:"Carga em pH", aromaticity:"Aromaticidade", secondary:"Estrutura secundária", weight:"Peso molecular", hydrophobic:"Índice hidrofóbico", composition:"Razão de composição", extinction:"Coeficiente de extinção", pH:"pH", empty:"Insira uma sequência para começar.", limit:"Conversões GET: 1.000 nt · Análise POST: 10.000 nt", copy:"Copiar", copied:"Copiado", saved:"JSON salvo", mutationHint:"Compare duas sequências de DNA codificante com o mesmo tamanho.", footer:"Desenvolvido para o BioSTAR", local:"Local-first · API", coming:"Em breve" },
  es: { brand:"BioSTAR", eyebrow:"Software de Bioinformática para Análisis e Investigación Dirigida", start:"Inicio", gene:"Sección de Genes", proteinSection:"Sección de Proteínas", experimental:"Sección Experimental", help:"Ayuda", results:"Resultados", compare:"Comparar", save:"Guardar JSON", upload:"Subir", optimize:"Optimizar Codones", howTo:"Cómo Usar", github:"Visitar en GitHub", title:"Análisis dirigido de secuencias biológicas", subtitle:"Analiza ADN, ARN y proteínas en un único espacio de trabajo.", sequenceType:"Tipo de secuencia", dna:"ADN", rna:"ARN", protein:"Proteína", peptideThreshold:"Umbral de tamaño del péptido", sequence:"Secuencia", reference:"Secuencia de referencia", analyze:"Analizar", mutation:"Comparar mutaciones", otherTools:"Otras Herramientas", codonTable:"Tabla de Codones", fasta:"Formateador FASTA", started:"¿Cómo Empezar?", viewer:"Visor 3D de Proteínas", sampleDna:"Cargar Muestra de ADN", sampleProtein:"Cargar Muestra de Proteína", options:"Opciones de análisis", full:"Suite completa", count:"Conteo de aminoácidos", pi:"Punto isoeléctrico", charge:"Carga a pH", aromaticity:"Aromaticidad", secondary:"Estructura secundaria", weight:"Peso molecular", hydrophobic:"Índice hidrofóbico", composition:"Razón de composición", extinction:"Coeficiente de extinción", pH:"pH", empty:"Introduce una secuencia para comenzar.", limit:"Conversiones GET: 1.000 nt · Análisis POST: 10.000 nt", copy:"Copiar", copied:"Copiado", saved:"JSON guardado", mutationHint:"Compara dos secuencias de ADN codificante del mismo tamaño.", footer:"Desarrollado para BioSTAR", local:"Local-first · API", coming:"Próximamente" }
} as const;

type Copy = (typeof COPY)[Locale];

const defaults: Options = { full:false, count:true, pi:true, charge:false, aromaticity:true, secondary:true, weight:true, hydrophobic:true, composition:false, extinction:false };

export default function Page() {
  const [locale, setLocale] = useState<Locale>("en");
  const [tool, setTool] = useState<Tool>("home");
  const [sequence, setSequence] = useState("");
  const [reference, setReference] = useState("");
  const [result, setResult] = useState<unknown>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const [copied, setCopied] = useState(false);
  const [saved, setSaved] = useState(false);
  const [options, setOptions] = useState(defaults);
  const [pH, setPH] = useState(7);
  const t = COPY[locale];

  const selectTool = (next: Tool) => { setTool(next); setResult(null); setError(""); setSaved(false); };
  const toggle = (key: keyof Options) => setOptions(o => ({ ...o, [key]: !o[key] }));

  async function run() {
    setLoading(true); setError(""); setResult(null); setSaved(false);
    try {
      let response: Response;
      if (tool === "protein") response = await fetch("/api/protein", { method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify({ sequence, get_full_test_results:options.full, get_aminoacids_count:options.count, get_isoelectric_point:options.pi, get_charge_at_pH:options.charge ? Number(pH) : null, get_aromaticity:options.aromaticity, get_secondary_structure_propensity:options.secondary, get_molecular_weight:options.weight, get_hydrophobic_index:options.hydrophobic, get_composition_ratio:options.composition, get_extinction_coefficient:options.extinction }) });
      else if (tool === "mutation") response = await fetch("/api/mutation_compare", { method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify({ reference, sequence }) });
      else { const paths: Record<string,string> = { dnaRna:"dna-rna", dnaProtein:"dna-protein", rnaProtein:"rna-protein", rnaDna:"rna-dna" }; response = await fetch(`/api/${paths[tool]}?sequence=${encodeURIComponent(sequence)}`); }
      const data = await response.json(); if (!response.ok) throw new Error(data.detail || "The analysis could not be completed."); setResult(data);
    } catch (e) { setError(e instanceof Error ? e.message : "The analysis could not be completed."); } finally { setLoading(false); }
  }

  async function copyResult() { if (!result) return; await navigator.clipboard.writeText(JSON.stringify(result,null,2)); setCopied(true); setTimeout(()=>setCopied(false),1500); }
  function saveJson() { if (!result) return; const payload = { application:"BioSTAR", version:"3", saved_at:new Date().toISOString(), tool, sequence:tool === "mutation" ? {reference,sequence} : sequence, options:tool === "protein" ? options : undefined, result }; const url=URL.createObjectURL(new Blob([JSON.stringify(payload,null,2)],{type:"application/json"})); const a=document.createElement("a"); a.href=url; a.download=`biostar-${tool}-${new Date().toISOString().replace(/[:.]/g,"-")}.json`; a.click(); URL.revokeObjectURL(url); setSaved(true); setTimeout(()=>setSaved(false),1800); }
  function loadDnaSample() { setTool("dnaProtein"); setSequence(">DNA Sample 1\\nGATCTTTGAGAAAGGGGATTTTAATGGTCAGATGCATGAGACCACGGAAGACTGCCCTTCCATCATGGAGCAGTTCCACATGCGGGAGGTCCACTCCTGTAAGGTGCTGGAGGGCGCCTGGATCTTCTATGAGCTGCCCAACTACCGAGGCAGGCAGTACCTGCTGGACAAGAAGGAGTACCGGAAGCCCGTCGACTGGGGTGCAGCTTCCCCAGCTGTCCAGTCTTTCCGCCGCATTGTGGAGTGATGATACAGATGCGGCCAAAC"); }
  function loadProteinSample() { setTool("protein"); setSequence(">Protein Sample\\nMVLSPADKTNVKAAW"); }

  const label: Record<Tool,string> = { home:t.start, dnaRna:"DNA → RNA", dnaProtein:"DNA → Protein", rnaProtein:"RNA → Protein", rnaDna:"RNA → DNA", protein:t.protein, mutation:t.mutation };
  const canRun = tool === "mutation" ? Boolean(sequence && reference) : Boolean(sequence);

  return <main className="app-shell"><Sidebar t={t} tool={tool} onTool={selectTool}/><div className="content-wrap">
    <header className="topbar"><span>{t.brand} <b>|</b> {t.eyebrow}</span><div className="topbar-actions"><div className="lang">{(["en","pt","es"] as Locale[]).map(l=><button key={l} className={locale===l?"active":""} onClick={()=>setLocale(l)}>{l.toUpperCase()}</button>)}</div><ThemeToggle/></div></header>
    <div className="content">
      {tool === "home" ? <Home t={t} onTool={selectTool} loadDnaSample={loadDnaSample} loadProteinSample={loadProteinSample}/> : <>
        <div className="page-heading"><div><span className="eyebrow">BioSTAR / ANALYSIS</span><h1>{String(label[tool])}</h1><p>{t.limit}</p></div>{result && <div className="heading-actions"><button onClick={saveJson}>{saved?t.saved:t.save}</button><button onClick={copyResult}>{copied?t.copied:t.copy}</button></div>}</div>
        <section className="analysis-card">{tool === "mutation" ? <div className="two-inputs"><SequenceBox label={t.reference} value={reference} onChange={setReference}/><SequenceBox label={t.sequence} value={sequence} onChange={setSequence}/></div> : <SequenceBox label={t.sequence} value={sequence} onChange={setSequence} fasta={tool === "protein"}/>} {tool === "protein" && <ProteinOptions t={t} options={options} toggle={toggle} pH={pH} setPH={setPH}/>}<button className="analyze" disabled={loading||!canRun} onClick={run}>{loading?"…":tool === "mutation"?t.compare:t.analyze}<span>→</span></button></section>
        <section className="result-section"><div className="section-title"><span className="eyebrow">OUTPUT</span><h2>{t.results}</h2></div><div className="result-box">{error?<div className="error">{error}</div>:result?<ResultView result={result}/>:<div className="empty"><div>⌁</div><p>{t.empty}</p></div>}</div></section>
      </>}
    </div><footer><span>{t.footer}</span><span>{t.local}</span></footer>
  </div></main>;
}

function Home({ t, onTool, loadDnaSample, loadProteinSample }: { t: Copy; onTool:(tool:Tool)=>void; loadDnaSample:()=>void; loadProteinSample:()=>void }) {
  return <><div className="home-grid"><section className="welcome"><span className="eyebrow">BIOINFORMATICS WORKSPACE</span><h1>{t.title}</h1><p>{t.subtitle}</p><div className="home-actions"><button className="primary" onClick={()=>onTool("dnaProtein")}>{t.analyze} DNA → Protein</button><button onClick={()=>onTool("protein")}>{t.protein} Analysis</button></div></section><section className="sequence-launch"><div className="launch-title">Quick samples</div><p className="launch-copy">Start with a known sequence or jump directly into an analysis.</p><div className="sample-actions"><button onClick={loadDnaSample}>{t.sampleDna}</button><button onClick={loadProteinSample}>{t.sampleProtein}</button></div></section></div><div className="home-lower"><section><div className="section-title"><span className="eyebrow">TOOLS</span><h2>{t.gene} / {t.proteinSection}</h2></div><div className="tool-cards"><button onClick={()=>onTool("dnaRna")}>DNA → RNA</button><button onClick={()=>onTool("dnaProtein")}>DNA → Protein</button><button onClick={()=>onTool("rnaProtein")}>RNA → Protein</button><button onClick={()=>onTool("rnaDna")}>RNA → DNA</button><button onClick={()=>onTool("mutation")}>{t.mutation}</button><button onClick={()=>onTool("protein")}>{t.protein} Analysis</button></div></section><section className="other-links"><div className="section-title"><span className="eyebrow">UTILITY</span><h2>{t.otherTools}</h2></div><button>{t.codonTable}</button><button>{t.fasta}</button><button>{t.started}</button><button disabled>{t.viewer} · {t.coming}</button></section></div></>;
}

function ProteinOptions({ t, options, toggle, pH, setPH }: { t: Copy; options:Options; toggle:(key:keyof Options)=>void; pH:number; setPH:(v:number)=>void }) {
  const fields:[keyof Options,string][] = [["count",t.count],["pi",t.pi],["charge",t.charge],["aromaticity",t.aromaticity],["secondary",t.secondary],["weight",t.weight],["hydrophobic",t.hydrophobic],["composition",t.composition],["extinction",t.extinction]];
  return <div className="protein-options"><div className="options-head"><h3>{t.options}</h3><button onClick={()=>toggle("full")}>{options.full?"Clear":"Enable full suite"}</button></div><label className="full-option"><input type="checkbox" checked={options.full} onChange={()=>toggle("full")}/>{t.full}</label><div className="options-grid">{fields.map(([key,text])=><label key={key}><input type="checkbox" checked={options[key]} onChange={()=>toggle(key)}/><span>{text}</span>{key==="charge"&&options.charge&&<input className="ph" type="number" step="0.1" value={pH} onChange={e=>setPH(Number(e.target.value))}/>}</label>)}</div></div>;
}
