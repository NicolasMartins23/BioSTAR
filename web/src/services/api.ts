const API_BASE_URL: string = import.meta.env.VITE_API_BASE_URL ?? "/api";

export interface HealthResponse { status: string; }
export interface SequenceResponse { sequence: string; }
export interface ProteinAnalysisRequest {
  sequence: string;
  get_full_test_results: boolean;
  get_aminoacids_count: boolean;
  get_isoelectric_point: boolean;
  get_charge_at_pH: number | null;
  get_aromaticity: boolean;
  get_secondary_structure_propensity: boolean;
  get_molecular_weight: boolean;
  get_hydrophobic_index: boolean;
  get_composition_ratio: boolean;
  get_extinction_coefficient: boolean;
}
export interface ProteinAnalysisResponse {
  sequence: string;
  length: number;
  aminoacids_count?: Record<string, number>;
  isoelectric_point?: number;
  charge_at_pH?: { pH: number; charge: number };
  aromaticity?: number;
  secondary_structure_propensity?: Record<string, number> | number;
  molecular_weight?: number;
  hydrophobic_index?: number;
  composition_ratio?: Record<string, number>;
  extinction_coefficient?: number | Record<string, number>;
}
const request = async <T>(path: string, init?: RequestInit): Promise<T> => {
  const response: Response = await fetch(\`\${API_BASE_URL}\${path}\`, init);
  if (!response.ok) {
    let detail: string = \`API request failed with status \${response.status}\`;
    try {
      const body: unknown = await response.json();
      if (typeof body === "object" && body !== null && "detail" in body && typeof body.detail === "string") {
        detail = body.detail;
      }
    } catch {}
    throw new Error(detail);
  }
  return response.json() as Promise<T>;
};
export const getHealth = async (): Promise<HealthResponse> => request<HealthResponse>("/v1/health");
export const convertDnaToRna = async (sequence: string): Promise<SequenceResponse> => request<SequenceResponse>(\`/dna-rna?sequence=\${encodeURIComponent(sequence)}\`);
export const convertDnaToProtein = async (sequence: string): Promise<SequenceResponse> => request<SequenceResponse>(\`/dna-protein?sequence=\${encodeURIComponent(sequence)}\`);
export const convertRnaToProtein = async (sequence: string): Promise<SequenceResponse> => request<SequenceResponse>(\`/rna-protein?sequence=\${encodeURIComponent(sequence)}\`);
export const convertRnaToDna = async (sequence: string): Promise<SequenceResponse> => request<SequenceResponse>(\`/rna-dna?sequence=\${encodeURIComponent(sequence)}\`);
export const analyzeProtein = async (requestData: ProteinAnalysisRequest): Promise<ProteinAnalysisResponse> => request<ProteinAnalysisResponse>("/protein", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify(requestData),
});
