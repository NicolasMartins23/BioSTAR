const API_BASE_URL: string = import.meta.env.VITE_API_BASE_URL ?? "/api";

export interface MessageResource {
  code: string;
  message: string;
}

export class APIError extends Error {
  public readonly status: number;
  public readonly messageResource: MessageResource | null;

  constructor(status: number, messageResource: MessageResource | null) {
    super(messageResource?.message ?? `API request failed with status ${status}`);
    this.name = "APIError";
    this.status = status;
    this.messageResource = messageResource;
  }
}

export interface APIResponse<T> {
  data: T | null;
  message: MessageResource | null;
}

export interface HealthResponse {
  status: string;
}

export interface SequenceResponse {
  sequence: string;
}

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

const request = async <T>(path: string, init?: RequestInit): Promise<APIResponse<T>> => {
  const response: Response = await fetch(`${API_BASE_URL}${path}`, init);

  let body: APIResponse<T>;

  try {
    body = await response.json() as APIResponse<T>;
  } catch {
    throw new Error(`API request failed with status ${response.status}`);
  }

  if (!response.ok) {
    throw new APIError(response.status, body.message);
  }

  return body;
};

export const getHealth = async (): Promise<APIResponse<HealthResponse>> => {
  return request<HealthResponse>("/v1/health");
};

export const convertDnaToRna = async (
  sequence: string,
): Promise<APIResponse<SequenceResponse>> => {
  return request<SequenceResponse>(`/dna-rna?sequence=${encodeURIComponent(sequence)}`);
};

export const convertDnaToProtein = async (
  sequence: string,
): Promise<APIResponse<SequenceResponse>> => {
  return request<SequenceResponse>(`/dna-protein?sequence=${encodeURIComponent(sequence)}`);
};

export const convertRnaToProtein = async (
  sequence: string,
): Promise<APIResponse<SequenceResponse>> => {
  return request<SequenceResponse>(`/rna-protein?sequence=${encodeURIComponent(sequence)}`);
};

export const convertRnaToDna = async (
  sequence: string,
): Promise<APIResponse<SequenceResponse>> => {
  return request<SequenceResponse>(`/rna-dna?sequence=${encodeURIComponent(sequence)}`);
};

export const analyzeProtein = async (
  requestData: ProteinAnalysisRequest,
): Promise<APIResponse<ProteinAnalysisResponse>> => {
  return request<ProteinAnalysisResponse>("/protein", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(requestData),
  });
};
