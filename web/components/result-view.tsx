export function ResultView({ result }: { result: unknown }) {
  if (!result || typeof result !== "object") return <pre>{String(result)}</pre>;
  const entries = Object.entries(result as Record<string, unknown>);
  return <div className="result-grid">{entries.map(([key, value]) => <div className="metric" key={key}>
    <span>{key.replaceAll("_", " ")}</span>
    {typeof value === "object" ? <pre>{JSON.stringify(value, null, 2)}</pre> : <strong>{String(value)}</strong>}
  </div>)}</div>;
}
