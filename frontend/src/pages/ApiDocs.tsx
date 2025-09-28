import React from 'react';
import {
  Box,
  Grid,
  Paper,
  Typography,
  Card,
  CardContent,
  Button,
  Chip,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
} from '@mui/material';
import {
  Dna,
  ChevronDown,
  Copy,
  Check,
  ExternalLink,
} from 'lucide-react';

const ApiDocs: React.FC = () => {
  const [copiedEndpoint, setCopiedEndpoint] = React.useState<string | null>(null);

  const endpoints = [
    {
      method: 'GET',
      path: '/health',
      description: 'Health Check',
      color: '#5B8FFF',
    },
    {
      method: 'POST',
      path: '/demo',
      description: 'Demo Task',
      color: '#10B981',
    },
    {
      method: 'GET',
      path: '/report/{rid}',
      description: 'Get Report',
      color: '#5B8FFF',
    },
    {
      method: 'POST',
      path: '/files',
      description: 'Upload File',
      color: '#10B981',
    },
    {
      method: 'GET',
      path: '/jobs/{job_id}',
      description: 'Job Status',
      color: '#5B8FFF',
    },
    {
      method: 'GET',
      path: '/reports/{job_id}',
      description: 'Download Report',
      color: '#5B8FFF',
    },
    {
      method: 'POST',
      path: '/datasets/upload',
      description: 'Datasets Upload',
      color: '#10B981',
    },
    {
      method: 'POST',
      path: '/variants/filter',
      description: 'Variants Filter',
      color: '#10B981',
    },
  ];

  const schemas = [
    {
      name: 'Body_datasets_upload_datasets_upload_post',
      description: 'Upload dataset schema',
      expand: 'object',
    },
    {
      name: 'Body_upload_file_files_post',
      description: 'File upload schema',
      expand: 'object',
    },
    {
      name: 'EvidenceCard',
      description: 'Evidence card data structure',
      expand: 'object',
    },
    {
      name: 'FilterResponse',
      description: 'Filter response format',
      expand: 'object',
    },
    {
      name: 'HTTPValidationError',
      description: 'Validation error schema',
      expand: 'object',
    },
    {
      name: 'UploadResponse',
      description: 'Upload response format',
      expand: 'object',
    },
    {
      name: 'ValidationError',
      description: 'Validation error details',
      expand: 'object',
    },
  ];

  const handleCopy = (text: string) => {
    navigator.clipboard.writeText(text);
    setCopiedEndpoint(text);
    setTimeout(() => setCopiedEndpoint(null), 2000);
  };

  return (
    <Box sx={{ p: 3, bgcolor: '#0f1419', minHeight: '100vh' }}>
      {/* Header */}
      <Box sx={{ display: 'flex', alignItems: 'center', mb: 4 }}>
        <Dna size={32} style={{ marginRight: 12, color: '#10B981' }} />
        <Typography variant="h4" sx={{ fontWeight: 600, color: '#10B981' }}>
          GenoScope API
        </Typography>
        
        <Box sx={{ ml: 2, display: 'flex', gap: 1 }}>
          <Chip label="v1.0" size="small" sx={{ bgcolor: '#2D3748' }} />
          <Chip label="OAS 3.1" size="small" sx={{ bgcolor: '#10B981', color: 'white' }} />
        </Box>
        
        <Box sx={{ ml: 'auto' }}>
          <Typography variant="body2" color="text.secondary">
            /openapi.json
          </Typography>
        </Box>
      </Box>

      {/* Endpoints Section */}
      <Box sx={{ mb: 4 }}>
        <Accordion
          defaultExpanded
          sx={{
            bgcolor: '#1a1f2e',
            '&:before': { display: 'none' },
            boxShadow: 'none',
            border: '1px solid #2D3748',
          }}
        >
          <AccordionSummary
            expandIcon={<ChevronDown />}
            sx={{ '&:hover': { bgcolor: 'rgba(45, 55, 72, 0.2)' } }}
          >
            <Typography variant="h6">default</Typography>
          </AccordionSummary>
          <AccordionDetails>
            <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
              {endpoints.map((endpoint, index) => (
                <Paper
                  key={index}
                  sx={{
                    p: 2,
                    display: 'flex',
                    alignItems: 'center',
                    gap: 2,
                    cursor: 'pointer',
                    transition: 'all 0.2s',
                    '&:hover': {
                      bgcolor: 'rgba(91, 143, 255, 0.05)',
                      transform: 'translateX(4px)',
                    },
                  }}
                  onClick={() => handleCopy(endpoint.path)}
                >
                  <Chip
                    label={endpoint.method}
                    size="small"
                    sx={{
                      bgcolor: endpoint.color,
                      color: 'white',
                      fontWeight: 600,
                      minWidth: 60,
                    }}
                  />
                  <Typography variant="body2" sx={{ fontFamily: 'monospace', flex: 1 }}>
                    {endpoint.path}
                  </Typography>
                  <Typography variant="body2" color="text.secondary">
                    {endpoint.description}
                  </Typography>
                  {copiedEndpoint === endpoint.path ? (
                    <Check size={16} style={{ color: '#10B981' }} />
                  ) : (
                    <Copy size={16} style={{ color: '#9CA3AF' }} />
                  )}
                </Paper>
              ))}
            </Box>
          </AccordionDetails>
        </Accordion>
      </Box>

      {/* Schemas Section */}
      <Accordion
        sx={{
          bgcolor: '#1a1f2e',
          '&:before': { display: 'none' },
          boxShadow: 'none',
          border: '1px solid #2D3748',
        }}
      >
        <AccordionSummary
          expandIcon={<ChevronDown />}
          sx={{ '&:hover': { bgcolor: 'rgba(45, 55, 72, 0.2)' } }}
        >
          <Typography variant="h6">Schemas</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
            {schemas.map((schema, index) => (
              <Paper
                key={index}
                sx={{
                  p: 2,
                  display: 'flex',
                  alignItems: 'center',
                  gap: 2,
                  cursor: 'pointer',
                  transition: 'all 0.2s',
                  '&:hover': {
                    bgcolor: 'rgba(91, 143, 255, 0.05)',
                  },
                }}
              >
                <Typography variant="body2" sx={{ fontFamily: 'monospace', fontWeight: 600 }}>
                  {schema.name}
                </Typography>
                <ChevronDown size={16} />
                <Typography variant="caption" color="text.secondary">
                  Expand all
                </Typography>
                <Chip
                  label={schema.expand}
                  size="small"
                  sx={{ bgcolor: '#2D3748', ml: 'auto' }}
                />
              </Paper>
            ))}
          </Box>
        </AccordionDetails>
      </Accordion>

      {/* Example Request */}
      <Box sx={{ mt: 4 }}>
        <Typography variant="h6" gutterBottom>
          Example Request
        </Typography>
        <Paper sx={{ p: 3, bgcolor: '#0a0d13', border: '1px solid #2D3748' }}>
          <Typography
            variant="body2"
            sx={{
              fontFamily: 'monospace',
              color: '#10B981',
              whiteSpace: 'pre-wrap',
            }}
          >
{`curl -X POST "http://localhost:8000/api/demo" \\
  -H "Content-Type: application/json" \\
  -d '{
    "analysis_type": "variant_calling",
    "reference": "hg38",
    "quality_threshold": 30
  }'`}
          </Typography>
        </Paper>
      </Box>
    </Box>
  );
};

export default ApiDocs;
