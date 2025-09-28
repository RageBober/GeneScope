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
} from '@mui/material';
import {
  Dna,
  FileText,
  FolderOpen,
  Beaker,
  Zap,
  BarChart3,
  Filter,
  Link,
  HelpCircle,
  BookOpen,
  Code,
  Users,
  Server,
  Database,
  HardDrive,
} from 'lucide-react';
import { useNavigate } from 'react-router-dom';

const Help: React.FC = () => {
  const navigate = useNavigate();

  const features = [
    {
      icon: <FileText size={32} />,
      title: 'Supported File Formats',
      description: 'VCF, BAM, CSV, TSV, FASTA, GFF, HDFS, Excel files. Upload single files or process in batches.',
      color: '#5B8FFF',
    },
    {
      icon: <Beaker size={32} />,
      title: 'Analysis Presets',
      description: 'Pre-configured analysis pipelines for rare diseases, oncology, and pharmacogenomics research.',
      color: '#EF4444',
    },
    {
      icon: <BarChart3 size={32} />,
      title: 'Interactive Reports',
      description: 'Generate comprehensive HTML/PDF reports with interactive visualizations and detailed statistics.',
      color: '#10B981',
    },
    {
      icon: <Filter size={32} />,
      title: 'Advanced Filtering',
      description: 'Filter variants by allele frequency, quality scores, genes, chromosomes, and custom criteria.',
      color: '#F59E0B',
    },
    {
      icon: <Zap size={32} />,
      title: 'Batch Processing',
      description: 'Process multiple files simultaneously with consistent parameters and automated reporting.',
      color: '#8B5CF6',
    },
    {
      icon: <Link size={32} />,
      title: 'API Access',
      description: 'Programmatic access via REST API for integration with your bioinformatics workflows.',
      color: '#EC4899',
    },
  ];

  return (
    <Box sx={{ p: 3, bgcolor: '#0f1419', minHeight: '100vh' }}>
      {/* Header */}
      <Box sx={{ display: 'flex', alignItems: 'center', mb: 4 }}>
        <Dna size={32} style={{ marginRight: 12, color: '#10B981' }} />
        <Typography variant="h4" sx={{ fontWeight: 600, color: '#10B981' }}>
          GenoScope
        </Typography>
        
        <Box sx={{ ml: 'auto', display: 'flex', gap: 2 }}>
          <Button variant="text" sx={{ color: 'text.primary' }}>Analysis</Button>
          <Button variant="text" sx={{ color: 'text.primary' }}>Results</Button>
          <Button variant="text" sx={{ color: 'text.primary', fontWeight: 600 }}>Help</Button>
          <Button variant="text" sx={{ color: 'text.primary' }}>API Docs</Button>
        </Box>
      </Box>

      <Grid container spacing={3}>
        {/* Quick Actions */}
        <Grid item xs={12} md={3}>
          <Grid container spacing={2}>
            <Grid item xs={12}>
              <Card>
                <CardContent>
                  <Typography variant="h6" gutterBottom sx={{ mb: 3 }}>
                    Quick Actions
                  </Typography>
                  
                  <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                    <Button
                      fullWidth
                      variant="contained"
                      startIcon={<Zap size={18} />}
                      onClick={() => navigate('/analysis')}
                      sx={{
                        bgcolor: '#5B8FFF',
                        '&:hover': { bgcolor: '#4169E1' },
                        py: 1.5,
                      }}
                    >
                      New Analysis
                    </Button>
                    
                    <Button
                      fullWidth
                      variant="outlined"
                      startIcon={<FolderOpen size={18} />}
                      onClick={() => navigate('/results')}
                      sx={{
                        borderColor: '#2D3748',
                        color: 'text.primary',
                        py: 1.5,
                        '&:hover': {
                          borderColor: '#5B8FFF',
                          bgcolor: 'rgba(91, 143, 255, 0.1)',
                        },
                      }}
                    >
                      View Results
                    </Button>
                    
                    <Button
                      fullWidth
                      variant="outlined"
                      startIcon={<FileText size={18} />}
                      sx={{
                        borderColor: '#2D3748',
                        color: 'text.primary',
                        py: 1.5,
                        '&:hover': {
                          borderColor: '#EF4444',
                          bgcolor: 'rgba(239, 68, 68, 0.1)',
                        },
                      }}
                    >
                      Load Example
                    </Button>
                  </Box>
                </CardContent>
              </Card>
            </Grid>

            {/* System Status */}
            <Grid item xs={12}>
              <Card>
                <CardContent>
                  <Typography variant="h6" gutterBottom>
                    System Status
                  </Typography>
                  
                  <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
                    <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                      <Typography variant="body2">API Server</Typography>
                      <Chip
                        label="Ready"
                        size="small"
                        sx={{
                          bgcolor: 'rgba(16, 185, 129, 0.2)',
                          color: '#10B981',
                        }}
                      />
                    </Box>
                    
                    <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                      <Typography variant="body2">Database</Typography>
                      <Chip
                        label="Connected"
                        size="small"
                        sx={{
                          bgcolor: 'rgba(16, 185, 129, 0.2)',
                          color: '#10B981',
                        }}
                      />
                    </Box>
                    
                    <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                      <Typography variant="body2">Cache</Typography>
                      <Chip
                        label="Active"
                        size="small"
                        sx={{
                          bgcolor: 'rgba(16, 185, 129, 0.2)',
                          color: '#10B981',
                        }}
                      />
                    </Box>
                  </Box>
                </CardContent>
              </Card>
            </Grid>

            {/* Recent Files */}
            <Grid item xs={12}>
              <Card>
                <CardContent>
                  <Typography variant="h6" gutterBottom>
                    Recent Files
                  </Typography>
                  
                  <Typography variant="body2" color="text.secondary">
                    No recent files
                  </Typography>
                </CardContent>
              </Card>
            </Grid>
          </Grid>
        </Grid>

        {/* Help & Documentation */}
        <Grid item xs={12} md={9}>
          <Box sx={{ mb: 3 }}>
            <Typography variant="h5" gutterBottom sx={{ color: 'text.primary' }}>
              Help & Documentation
            </Typography>
            <Typography variant="body2" color="text.secondary">
              Learn how to use GenoScope effectively
            </Typography>
          </Box>

          <Grid container spacing={3}>
            {features.map((feature, index) => (
              <Grid item xs={12} md={6} key={index}>
                <Card
                  sx={{
                    height: '100%',
                    transition: 'transform 0.2s, box-shadow 0.2s',
                    '&:hover': {
                      transform: 'translateY(-4px)',
                      boxShadow: '0 8px 24px rgba(0,0,0,0.3)',
                    },
                  }}
                >
                  <CardContent>
                    <Box sx={{ display: 'flex', alignItems: 'flex-start', gap: 2 }}>
                      <Box
                        sx={{
                          p: 1.5,
                          borderRadius: 2,
                          bgcolor: `${feature.color}20`,
                          color: feature.color,
                        }}
                      >
                        {feature.icon}
                      </Box>
                      <Box sx={{ flex: 1 }}>
                        <Typography variant="h6" gutterBottom>
                          {feature.title}
                        </Typography>
                        <Typography variant="body2" color="text.secondary">
                          {feature.description}
                        </Typography>
                      </Box>
                    </Box>
                  </CardContent>
                </Card>
              </Grid>
            ))}
          </Grid>
        </Grid>
      </Grid>
    </Box>
  );
};

export default Help;
