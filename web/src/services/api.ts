const API_BASE_URL: string = import.meta.env.VITE_API_BASE_URL ?? "/api";

export interface HealthResponse {
  status: string;
}

export interface SequenceResponse {
  sequence: string;
}

const request = async <T>(path: string, init?: RequestInit): Promise<T> => {
  const response: Response = await fetch(`${API_BASE_URL}${path}`, init);

  if (!response.ok) {
    let detail: string = `API request failed with status ${response.status}`;

    try {
      const body: unknown = await response.json();
      if (
        typeof body === "object"
        && body !== null
        && "detail" in body
        && typeof body.detail === "string"
      ) {
        detail = body.detail;
      }
    } catch {
      // Keep the status-based error when the response is not JSON.
    }

    throw new Error(detail);
  }

  return response.json() as Promise<T>;
};

export const getHealth = async (): Promise<HealthResponse> => {
  return request<HealthResponse>("/v1/health");
};

export const convertDnaToRna = async (sequence: string): Promise<SequenceResponse> => {
  return request<SequenceResponse>(`/dna-rna?sequence=${encodeURIComponent(sequence)}`);
};

export const convertDnaToProtein = async (sequence: string): Promise<SequenceResponse> => {
  return request<SequenceResponse>(`/dna-protein?sequence=${encodeURIComponent(sequence)}`);
};

export const convertRnaToProtein = async (sequence: string): Promise<SequenceResponse> => {
  return request<SequenceResponse>(`/rna-protein?sequence=${encodeURIComponent(sequence)}`);
};

export const convertRnaToDna = async (sequence: string): Promise<SequenceResponse> => {
  return request<SequenceResponse>(`/rna-dna?sequence=${encodeURIComponent(sequence)}`);
};
