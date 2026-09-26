"use client";

import { useMemo, useState } from "react";

const COPY = {
  en: { brand: "BioSTAR", eyebrow: "Biological sequence analysis", title: "Turn sequences into insight.", subtitle: "A focused workspace for DNA, RNA and protein analysis — powered by the BioSTAR engine.", dnaRna: "DNA → RNA", dnaProtein: "DNA → Protein", rnaProtein: "RNA → Protein", rnaDna: "RNA → DNA", protein: "Protein analysis", mutation: "Mutation compare", sequence: "Sequence", reference: "Reference sequence", analyze: "Analyze", compare: "Compare", results: "Results", clear: "Clear", copy: "Copy", copied: "Copied", options: "Analysis options", full: "Full test suite", count: "Amino-acid count", pi: "Isoelectric point", charge: "Charge at pH", aromaticity: "Aromaticity", secondary: "Secondary structure", weight: "Molecular weight", hydrophobic: "Hydrophobic index", composition: "Composition ratio", extinction: "Extinction coefficient", pH: "pH", empty: "Enter a sequence to begin.", limit: "GET conversions are limited to 1,000 nt.", mutationHint: "Compare two coding DNA sequences of equal length.", footer: "BioSTAR · Open sequence analysis" },
  pt: { brand: "BioSTAR", eyebrow: "Análise de sequências biológicas", title: "Transforme sequências em informação.", subtitle: "Um espaço focado para análise de DNA, RNA e proteínas — powered pelo motor BioSTAR.", dnaRna: "DNA → RNA", dnaProtein: "DNA → Proteína", rnaProtein: "RNA → Proteína", rnaDna: "RNA → DNA", protein: "Análise de proteína", mutation: "Comparar mutações", sequence: "Sequência", reference: "Sequência de referência", analyze: "Analisar", compare: "Comparar", results: "Resultados", clear: "Limpar", copy: "Copiar", copied: "Copiado", options: "Opções de análise", full: "Suite completa", count: "Contagem de aminoácidos", pi: "Ponto isoelétrico", charge: "Carga em pH", aromaticity: "Aromaticidade", secondary: "Estrutura secundária", weight: "Peso molecular", hydrophobic: "Índice hidrofóbico", composition: "Razão de composição", extinction: "Coeficiente de extinção", pH: "pH", empty: "Insira uma sequência para começar.", limit: "Conversões GET são limitadas a 1.000 nt.", mutationHint: "Compare duas sequências de DNA codificante com o mesmo tamanho.", footer: "BioSTAR · Análise de sequências" },
  es: { brand: "BioSTAR", eyebrow: "Análisis de secuencias biológicas", title: "Convierte secuencias en información.", subtitle: "Un espacio enfocado para analizar ADN, ARN y proteínas — impulsado por el motor BioSTAR.", dnaRna: "ADN → ARN", dnaProtein: "ADN → Proteína", rnaProtein: "ARN → Proteína", rnaDna: "ARN → ADN", protein: "Análisis de proteína", mutation: "Comparar mutaciones", sequence: "Secuencia", reference: "Secuencia de referencia", analyze: "Analizar", compare: "Comparar", results: "Resultados", clear: "Limpiar", copy: "Copiar", copied: "Copiado", options: "Opciones de análisis", full: "Suite completa", count: "Conteo de aminoácidos", pi: "Punto isoeléctrico", charge: "Carga a pH", aromaticity: "Aromaticidad", secondary: "Estructura secundaria", weight: "Peso molecular", hydrophobic: "Índice hidrofóbico", composition: "Razón de composición", extinction: "Coeficiente de extinción", pH: "pH", empty: "Introduce una secuencia para comenzar.", limit: "Las conversiones GET están limitadas a 1.000 nt.", mutationHint: "Compara dos secuencias de ADN codificante del mismo tamaño.", footer: "BioSTAR · Análisis de secuencias" }
} as const;

type Locale = keyof typeof COPY;
type Tool = "dnaRna" | "dnaProtein" | "rnaProtein" | "rnaDna" | "protein" | "mutation";

export default function Page({ params }: { params: Promise<{ locale: string }> }) {
  const [locale, setLocale] = useState<Locale>("en");
  const [tool, setTool] = useState<Tool>("protein");
  const [sequence, setSequence] = useState("");
  const [reference, setReference] = useState("");
  const [result, setResult] = useState<unknown>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");
  const [copied, setCopied] = useState(false);
  const [options, setOptions] = useState({ full: false, count: true, pi: true, charge: false, aromaticity: true, secondary: true, weight: true, hydrophobic: true, composition: false, extinction: false });
  const [pH, setPH] = useState(7);
  const t = COPY[locale];

  const tools = useMemo(() => [
    ["dnaRna", t.dnaRna], ["dnaProtein", t.dnaProtein], ["rnaProtein", t.rnaProtein], ["rnaDna", t.rnaDna], ["protein", t.protein], ["mutation", t.mutation]
  ] as [Tool, string][], [t]);

  async function run() {
    setLoading(true); setError(""); setResult(null); setCopied(false);
    try {
      let response: Response;
      if (tool === "protein") {
        response = await fetch("/api/protein", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ sequence, get_full_test_results: options.full, get_aminoacids_count: options.count, get_isoelectric_point: options.pi, get_charge_at_pH: options.charge ? Number(pH) : null, get_aromaticity: options.aromaticity, get_secondary_structure_propensity: options.secondary, get_molecular_weight: options.weight, get_hydrophobic_index: options.hydrophobic, get_composition_ratio: options.composition, get_extinction_coefficient: options.extinction }) });
      } else if (tool === "mutation") {
        response = await fetch("/api/mutation_compare", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ reference, sequence }) });
      } else {
        const paths: Record<string, string> = { dnaRna: "dna-rna", dnaProtein: "dna-protein", rnaProtein: "rna-protein", rnaDna: "rna-dna" };
        response = await fetch(`/api/${paths[tool]}?sequence=${encodeURIComponent(sequence)}`);
      }
      const data = await response.json();
      if (!response.ok) throw new Error(data.detail || "The analysis could not be completed.");
      setResult(data);
    } catch (e) { setError(e instanceof Error ? e.message : "The analysis could not be completed."); }
    finally { setLoading(false); }
  }

  function toggle(key: keyof typeof options) { setOptions((o) => ({ ...o, [key]: !o[key] })); }
  async function copyResult() { await navigator.clipboard.writeText(JSON.stringify(result, null, 2)); setCopied(true); setTimeout(() => setCopied(false), 1500); }

  return <main className="shell">
    <header className="topbar">
      <div className="brand"><div className="brand-mark">B</div><div><strong>{t.brand}</strong><span>{t.eyebrow}</span></div></div>
      <div className="lang"><button className={locale === "en" ? "active" : ""} onClick={() => setLocale("en")}>EN</button><button className={locale === "pt" ? "active" : ""} onClick={() => setLocale("pt")}>PT</button><button className={locale === "es" ? "active" : ""} onClick={() => setLocale("es")}>ES</button></div>
    </header>

    <section className="hero"><div className="hero-copy"><div className="kicker">{t.eyebrow}</div><h1>{t.title}</h1><p>{t.subtitle}</p></div><div className="hero-orbit"><div className="orbit-dot dot-a">A</div><div className="orbit-dot dot-t">T</div><div className="orbit-dot dot-g">G</div><div className="orbit-dot dot-c">C</div><div className="core">DNA</div></div></section>

    <section className="workspace">
      <aside className="sidebar"><div className="side-label">Tools</div>{tools.map(([id, label]) => <button key={id} className={`tool ${tool === id ? "selected" : ""}`} onClick={() => { setTool(id); setResult(null); setError(""); }}><span className="tool-dot" />{label}</button>)}</aside>
      <div className="panel">
        <div className="panel-head"><div><span className="eyebrow">{t.brand}</span><h2>{tools.find(([id]) => id === tool)?.[1]}</h2></div><span className="limit">{tool === "protein" || tool === "mutation" ? "POST · 10,000" : "GET · 1,000 nt"}</span></div>
        <div className="editor-grid">
          <div className="input-card"><label>{tool === "mutation" ? t.reference : t.sequence}</label><textarea value={tool === "mutation" ? reference : sequence} onChange={(e) => tool === "mutation" ? setReference(e.target.value) : setSequence(e.target.value)} placeholder={tool === "protein" ? ">my-protein\nMKWVTFISLLFLFSSAYSR" : tool === "mutation" ? "ATGGCCGAA..." : "ATGGCCGAA..."} spellCheck={false} />{tool === "mutation" && <><label className="second-label">{t.sequence}</label><textarea className="small-textarea" value={sequence} onChange={(e) => setSequence(e.target.value)} placeholder="ATGGTCGAA..." spellCheck={false} /></> }<div className="input-foot"><span>{tool === "mutation" ? t.mutationHint : t.limit}</span><span>{(tool === "mutation" ? Math.max(sequence.length, reference.length) : sequence.length).toLocaleString()} nt</span></div></div>
          {tool === "protein" && <div className="options-card"><div className="card-title">{t.options}<button className="select-all" onClick={() => setOptions((o) => ({ ...o, full: !o.full }))}>{t.full}</button></div><div className="option-list">{([['count','count'],['pi','pi'],['charge','charge'],['aromaticity','aromaticity'],['secondary','secondary'],['weight','weight'],['hydrophobic','hydrophobic'],['composition','composition'],['extinction','extinction']] as [keyof typeof options, keyof typeof t][]).map(([key,label]) => <label className="check" key={key}><input type="checkbox" checked={options[key]} onChange={() => toggle(key)} /><span>{t[label]}</span>{key === "charge" && options.charge && <input className="ph" type="number" step="0.1" value={pH} onChange={(e) => setPH(Number(e.target.value))} aria-label={t.pH} />}</label>)}</div></div>}
        </div>
        <button className="run" disabled={loading || (tool === "mutation" ? !sequence || !reference : !sequence)} onClick={run}>{loading ? "…" : tool === "mutation" ? t.compare : t.analyze}<span>↗</span></button>
        <div className="results-head"><div><span className="eyebrow">Output</span><h3>{t.results}</h3></div>{result && <button className="copy" onClick={copyResult}>{copied ? t.copied : t.copy}</button>}</div>
        <div className="result-box">{error ? <div className="error">{error}</div> : result ? <ResultView result={result} /> : <div className="empty"><div className="empty-icon">⌁</div><p>{t.empty}</p></div>}</div>
      </div>
    </section>
    <footer>{t.footer}<span>Local-first · API powered</span></footer>
  </main>;
}

function ResultView({ result }: { result: unknown }) {
  if (typeof result !== "object" || result === null) return <pre>{String(result)}</pre>;
  const entries = Object.entries(result as Record<string, unknown>);
  return <div className="result-grid">{entries.map(([key, value]) => <div className="metric" key={key}><span>{key.replaceAll("_", " ")}</span>{typeof value === "object" ? <pre>{JSON.stringify(value, null, 2)}</pre> : <strong>{String(value)}</strong>}</div>)}</div>;
}
