const API_BASE_URL: string = import.meta.env.VITE_API_BASE_URL ?? "/api";

export interface HealthResponse {
  status: string;
}

const request = async <T>(path: string, init?: RequestInit): Promise<T> => {
  const response: Response = await fetch(`${API_BASE_URL}${path}`, init);

  if (!response.ok) {
    throw new Error(`API request failed with status ${response.status}`);
  }

  return response.json() as Promise<T>;
};

export const getHealth = async (): Promise<HealthResponse> => {
  return request<HealthResponse>("/v1/health");
};
