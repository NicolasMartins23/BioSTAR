"use client";

export function SequenceBox({ label, value, onChange, fasta = false }: { label: string; value: string; onChange: (value: string) => void; fasta?: boolean }) {
  const length = value.replace(/^>[^\n]*\n?/, "").replace(/\s/g, "").length;
  return <div className="sequence-box">
    <div className="sequence-label-row"><label>{label}</label><span>{length.toLocaleString()} nt / aa</span></div>
    <textarea value={value} onChange={(e) => onChange(e.target.value)} placeholder={fasta ? ">Protein name\nSEQUENCE" : "Enter sequence (FASTA format supported)"} spellCheck={false} />
    <div className="sequence-meta"><span>{fasta ? "FASTA" : "DNA / RNA / FASTA"}</span><span>Paste or type your sequence</span></div>
  </div>;
}
