"use client";

type Tool = "home" | "dnaRna" | "dnaProtein" | "rnaProtein" | "rnaDna" | "protein" | "mutation";

type Props = { t: Record<string, string>; tool: Tool; onTool: (tool: Tool) => void };

function Section({ title, icon, children }: { title: string; icon: string; children: React.ReactNode }) {
  return <div className="nav-section"><div className="nav-section-title"><span>{icon}</span>{title}</div>{children}</div>;
}

export function Sidebar({ t, tool, onTool }: Props) {
  const item = (next: Tool, label: string) => <button className={`nav-sub ${tool === next ? "active-sub" : ""}`} onClick={() => onTool(next)}>↳ {label}</button>;
  return <aside className="app-sidebar">
    <div className="brand"><div className="brand-mark">★</div><div><strong>{t.brand}</strong><span>{t.eyebrow}</span></div></div>
    <nav>
      <button className={`nav-item ${tool === "home" ? "active" : ""}`} onClick={() => onTool("home")}><span>★</span>{t.start}</button>
      <Section title={t.gene} icon="DNA">{item("dnaProtein", t.results)}{item("mutation", t.compare)}<button className="nav-sub disabled">↳ {t.save}</button><button className="nav-sub disabled">↳ {t.upload}</button></Section>
      <Section title={t.proteinSection} icon="◈">{item("protein", t.results)}<button className="nav-sub disabled">↳ {t.save}</button><button className="nav-sub disabled">↳ {t.upload}</button></Section>
      <Section title={t.experimental} icon="⚗"><button className="nav-sub disabled">↳ {t.optimize}</button><button className="nav-sub disabled">↳ {t.save}</button><button className="nav-sub disabled">↳ {t.upload}</button></Section>
      <Section title={t.help} icon="?"><button className="nav-sub disabled">↳ {t.howTo}</button></Section>
      <a className="nav-item github" href="https://github.com/NicolasMartins23/BioSTAR" target="_blank" rel="noopener noreferrer"><span>⌘</span>{t.github}</a>
    </nav>
  </aside>;
}
