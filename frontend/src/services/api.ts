/**
 * API Service for GenoScope
 * Handles all API communications with the backend
 */

const API_BASE_URL = process.env.REACT_APP_API_URL || 'http://localhost:8000/api/v1';
const WS_BASE_URL = process.env.REACT_APP_WS_URL || 'ws://localhost:8000/ws';

// Token management
let accessToken: string | null = localStorage.getItem('access_token');
let refreshToken: string | null = localStorage.getItem('refresh_token');

/**
 * Base API request handler
 */
async function apiRequest(
  endpoint: string,
  options: RequestInit = {}
): Promise<any> {
  const url = `${API_BASE_URL}${endpoint}`;
  
  const headers = {
    'Content-Type': 'application/json',
    ...(accessToken && { Authorization: `Bearer ${accessToken}` }),
    ...options.headers,
  };

  try {
    const response = await fetch(url, {
      ...options,
      headers,
    });

    if (response.status === 401) {
      // Try to refresh token
      const refreshed = await refreshAccessToken();
      if (refreshed) {
        // Retry request with new token
        return apiRequest(endpoint, options);
      } else {
        // Redirect to login
        window.location.href = '/login';
        throw new Error('Authentication failed');
      }
    }

    if (!response.ok) {
      const error = await response.json();
      throw new Error(error.detail || 'API request failed');
    }

    return response.json();
  } catch (error) {
    console.error('API request failed:', error);
    throw error;
  }
}

/**
 * Refresh access token using refresh token
 */
async function refreshAccessToken(): Promise<boolean> {
  if (!refreshToken) return false;

  try {
    const response = await fetch(`${API_BASE_URL}/auth/refresh`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ refresh_token: refreshToken }),
    });

    if (response.ok) {
      const data = await response.json();
      accessToken = data.access_token;
      if (accessToken) {
        localStorage.setItem('access_token', accessToken);
      }
      return true;
    }
  } catch (error) {
    console.error('Token refresh failed:', error);
  }

  return false;
}

/**
 * Authentication API
 */
export const authAPI = {
  async login(email: string, password: string) {
    const response = await fetch(`${API_BASE_URL}/auth/login`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ email, password }),
    });

    if (!response.ok) {
      throw new Error('Invalid credentials');
    }

    const data = await response.json();
    accessToken = data.access_token;
    refreshToken = data.refresh_token;
    if (accessToken) {
      localStorage.setItem('access_token', accessToken);
    }
    if (refreshToken) {
      localStorage.setItem('refresh_token', refreshToken);
    }
    
    return data;
  },

  async register(userData: {
    email: string;
    password: string;
    username: string;
    full_name: string;
  }) {
    return apiRequest('/auth/register', {
      method: 'POST',
      body: JSON.stringify(userData),
    });
  },

  logout() {
    accessToken = null;
    refreshToken = null;
    localStorage.removeItem('access_token');
    localStorage.removeItem('refresh_token');
    window.location.href = '/login';
  },

  async getCurrentUser() {
    return apiRequest('/auth/me');
  },

  async createApiKey() {
    return apiRequest('/auth/api-key', { method: 'POST' });
  },

  async revokeApiKey(apiKey: string) {
    return apiRequest(`/auth/api-key/${apiKey}`, { method: 'DELETE' });
  },
};

/**
 * Analysis API
 */
export const analysisAPI = {
  async createAnalysis(data: {
    name: string;
    pipeline: string;
    files: string[];
    parameters: any;
  }) {
    return apiRequest('/analyses', {
      method: 'POST',
      body: JSON.stringify(data),
    });
  },

  async getAnalyses(params?: {
    status?: string;
    pipeline?: string;
    limit?: number;
    offset?: number;
  }) {
    const queryParams = new URLSearchParams(params as any).toString();
    return apiRequest(`/analyses?${queryParams}`);
  },

  async getAnalysis(id: string) {
    return apiRequest(`/analyses/${id}`);
  },

  async getAnalysisStatus(id: string) {
    return apiRequest(`/analyses/${id}/status`);
  },

  async cancelAnalysis(id: string) {
    return apiRequest(`/analyses/${id}/cancel`, { method: 'POST' });
  },

  async deleteAnalysis(id: string) {
    return apiRequest(`/analyses/${id}`, { method: 'DELETE' });
  },

  async getAnalysisLogs(id: string) {
    return apiRequest(`/analyses/${id}/logs`);
  },

  async retryAnalysis(id: string) {
    return apiRequest(`/analyses/${id}/retry`, { method: 'POST' });
  },
};

/**
 * Results API
 */
export const resultsAPI = {
  async getResults(analysisId: string) {
    return apiRequest(`/results/${analysisId}`);
  },

  async getVariants(analysisId: string, params?: {
    chromosome?: string;
    gene?: string;
    impact?: string;
    pathogenic?: boolean;
    limit?: number;
    offset?: number;
  }) {
    const queryParams = new URLSearchParams(params as any).toString();
    return apiRequest(`/results/${analysisId}/variants?${queryParams}`);
  },

  async getVariantDetails(analysisId: string, variantId: string) {
    return apiRequest(`/results/${analysisId}/variants/${variantId}`);
  },

  async getCoverageStats(analysisId: string) {
    return apiRequest(`/results/${analysisId}/coverage`);
  },

  async getQCMetrics(analysisId: string) {
    return apiRequest(`/results/${analysisId}/qc`);
  },

  async exportResults(analysisId: string, format: 'vcf' | 'csv' | 'json') {
    return apiRequest(`/results/${analysisId}/export?format=${format}`);
  },

  async getIGVData(analysisId: string) {
    return apiRequest(`/results/${analysisId}/igv`);
  },
};

/**
 * Pipeline API
 */
export const pipelineAPI = {
  async getPipelines() {
    return apiRequest('/pipelines');
  },

  async getPipeline(id: string) {
    return apiRequest(`/pipelines/${id}`);
  },

  async createCustomPipeline(data: {
    name: string;
    description: string;
    steps: any[];
  }) {
    return apiRequest('/pipelines', {
      method: 'POST',
      body: JSON.stringify(data),
    });
  },

  async updatePipeline(id: string, data: any) {
    return apiRequest(`/pipelines/${id}`, {
      method: 'PUT',
      body: JSON.stringify(data),
    });
  },

  async deletePipeline(id: string) {
    return apiRequest(`/pipelines/${id}`, { method: 'DELETE' });
  },

  async validatePipeline(data: any) {
    return apiRequest('/pipelines/validate', {
      method: 'POST',
      body: JSON.stringify(data),
    });
  },
};

/**
 * Files API
 */
export const filesAPI = {
  async uploadFile(file: File, onProgress?: (progress: number) => void) {
    const formData = new FormData();
    formData.append('file', file);

    return new Promise((resolve, reject) => {
      const xhr = new XMLHttpRequest();

      xhr.upload.addEventListener('progress', (e) => {
        if (e.lengthComputable && onProgress) {
          const progress = (e.loaded / e.total) * 100;
          onProgress(progress);
        }
      });

      xhr.addEventListener('load', () => {
        if (xhr.status === 200) {
          resolve(JSON.parse(xhr.responseText));
        } else {
          reject(new Error('Upload failed'));
        }
      });

      xhr.addEventListener('error', () => {
        reject(new Error('Upload failed'));
      });

      xhr.open('POST', `${API_BASE_URL}/files/upload`);
      if (accessToken) {
        xhr.setRequestHeader('Authorization', `Bearer ${accessToken}`);
      }
      xhr.send(formData);
    });
  },

  async getFiles() {
    return apiRequest('/files');
  },

  async getFile(id: string) {
    return apiRequest(`/files/${id}`);
  },

  async deleteFile(id: string) {
    return apiRequest(`/files/${id}`, { method: 'DELETE' });
  },

  async getPresignedUrl(fileId: string, operation: 'get' | 'put' = 'get') {
    return apiRequest(`/files/${fileId}/presigned-url?operation=${operation}`);
  },
};

/**
 * Reports API
 */
export const reportsAPI = {
  async generateReport(data: {
    name: string;
    template: string;
    analysisIds: string[];
    sections: any;
    format: string;
  }) {
    return apiRequest('/reports/generate', {
      method: 'POST',
      body: JSON.stringify(data),
    });
  },

  async getReports() {
    return apiRequest('/reports');
  },

  async getReport(id: string) {
    return apiRequest(`/reports/${id}`);
  },

  async downloadReport(id: string) {
    const response = await fetch(`${API_BASE_URL}/reports/${id}/download`, {
      headers: {
        ...(accessToken && { Authorization: `Bearer ${accessToken}` }),
      },
    });

    if (!response.ok) {
      throw new Error('Download failed');
    }

    const blob = await response.blob();
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `report_${id}.pdf`;
    document.body.appendChild(a);
    a.click();
    window.URL.revokeObjectURL(url);
    document.body.removeChild(a);
  },

  async shareReport(id: string, emails: string[]) {
    return apiRequest(`/reports/${id}/share`, {
      method: 'POST',
      body: JSON.stringify({ emails }),
    });
  },

  async deleteReport(id: string) {
    return apiRequest(`/reports/${id}`, { method: 'DELETE' });
  },
};

/**
 * Projects API
 */
export const projectsAPI = {
  async getProjects() {
    return apiRequest('/projects');
  },

  async getProject(id: string) {
    return apiRequest(`/projects/${id}`);
  },

  async createProject(data: {
    name: string;
    description: string;
    collaborators?: string[];
  }) {
    return apiRequest('/projects', {
      method: 'POST',
      body: JSON.stringify(data),
    });
  },

  async updateProject(id: string, data: any) {
    return apiRequest(`/projects/${id}`, {
      method: 'PUT',
      body: JSON.stringify(data),
    });
  },

  async deleteProject(id: string) {
    return apiRequest(`/projects/${id}`, { method: 'DELETE' });
  },

  async addCollaborator(projectId: string, email: string, role: string) {
    return apiRequest(`/projects/${projectId}/collaborators`, {
      method: 'POST',
      body: JSON.stringify({ email, role }),
    });
  },

  async removeCollaborator(projectId: string, userId: string) {
    return apiRequest(`/projects/${projectId}/collaborators/${userId}`, {
      method: 'DELETE',
    });
  },
};

/**
 * Samples API
 */
export const samplesAPI = {
  async getSamples() {
    return apiRequest('/samples');
  },

  async getSample(id: string) {
    return apiRequest(`/samples/${id}`);
  },

  async createSample(data: {
    name: string;
    type: string;
    metadata: any;
  }) {
    return apiRequest('/samples', {
      method: 'POST',
      body: JSON.stringify(data),
    });
  },

  async updateSample(id: string, data: any) {
    return apiRequest(`/samples/${id}`, {
      method: 'PUT',
      body: JSON.stringify(data),
    });
  },

  async deleteSample(id: string) {
    return apiRequest(`/samples/${id}`, { method: 'DELETE' });
  },

  async linkFiles(sampleId: string, fileIds: string[]) {
    return apiRequest(`/samples/${sampleId}/files`, {
      method: 'POST',
      body: JSON.stringify({ file_ids: fileIds }),
    });
  },
};

/**
 * Settings API
 */
export const settingsAPI = {
  async getSettings() {
    return apiRequest('/settings');
  },

  async updateSettings(data: any) {
    return apiRequest('/settings', {
      method: 'PUT',
      body: JSON.stringify(data),
    });
  },

  async getNotificationSettings() {
    return apiRequest('/settings/notifications');
  },

  async updateNotificationSettings(data: any) {
    return apiRequest('/settings/notifications', {
      method: 'PUT',
      body: JSON.stringify(data),
    });
  },

  async getIntegrations() {
    return apiRequest('/settings/integrations');
  },

  async configureIntegration(integration: string, config: any) {
    return apiRequest(`/settings/integrations/${integration}`, {
      method: 'PUT',
      body: JSON.stringify(config),
    });
  },
};

/**
 * WebSocket connection for real-time updates
 */
export class WebSocketService {
  private ws: WebSocket | null = null;
  private reconnectAttempts = 0;
  private maxReconnectAttempts = 5;
  private reconnectDelay = 1000;
  private listeners: Map<string, Set<Function>> = new Map();

  connect() {
    const wsUrl = `${WS_BASE_URL}?token=${accessToken}`;
    this.ws = new WebSocket(wsUrl);

    this.ws.onopen = () => {
      console.log('WebSocket connected');
      this.reconnectAttempts = 0;
    };

    this.ws.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data);
        this.emit(data.type, data.payload);
      } catch (error) {
        console.error('Failed to parse WebSocket message:', error);
      }
    };

    this.ws.onerror = (error) => {
      console.error('WebSocket error:', error);
    };

    this.ws.onclose = () => {
      console.log('WebSocket disconnected');
      this.reconnect();
    };
  }

  private reconnect() {
    if (this.reconnectAttempts < this.maxReconnectAttempts) {
      this.reconnectAttempts++;
      console.log(`Reconnecting... (${this.reconnectAttempts}/${this.maxReconnectAttempts})`);
      setTimeout(() => this.connect(), this.reconnectDelay * this.reconnectAttempts);
    }
  }

  disconnect() {
    if (this.ws) {
      this.ws.close();
      this.ws = null;
    }
  }

  send(type: string, payload: any) {
    if (this.ws && this.ws.readyState === WebSocket.OPEN) {
      this.ws.send(JSON.stringify({ type, payload }));
    } else {
      console.error('WebSocket is not connected');
    }
  }

  on(event: string, callback: Function) {
    if (!this.listeners.has(event)) {
      this.listeners.set(event, new Set());
    }
    this.listeners.get(event)!.add(callback);
  }

  off(event: string, callback: Function) {
    const callbacks = this.listeners.get(event);
    if (callbacks) {
      callbacks.delete(callback);
    }
  }

  private emit(event: string, data: any) {
    const callbacks = this.listeners.get(event);
    if (callbacks) {
      callbacks.forEach(callback => callback(data));
    }
  }
}

// Export WebSocket instance
export const wsService = new WebSocketService();

// GraphQL client
export const graphqlQuery = async (query: string, variables?: any) => {
  const response = await fetch(`${API_BASE_URL.replace('/v1', '')}/graphql`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      ...(accessToken && { Authorization: `Bearer ${accessToken}` }),
    },
    body: JSON.stringify({ query, variables }),
  });

  if (!response.ok) {
    throw new Error('GraphQL query failed');
  }

  const data = await response.json();
  if (data.errors) {
    throw new Error(data.errors[0].message);
  }

  return data.data;
};

export default {
  auth: authAPI,
  analysis: analysisAPI,
  results: resultsAPI,
  pipeline: pipelineAPI,
  files: filesAPI,
  reports: reportsAPI,
  projects: projectsAPI,
  samples: samplesAPI,
  settings: settingsAPI,
  ws: wsService,
  graphql: graphqlQuery,
};
