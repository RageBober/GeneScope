import React, { useState, useCallback } from 'react';
import {
  Box,
  Paper,
  Typography,
  Stepper,
  Step,
  StepLabel,
  Button,
  Grid,
  Card,
  CardContent,
  TextField,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  Chip,
  Alert,
  LinearProgress,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  IconButton,
  Tooltip,
} from '@mui/material';
import { useNavigate } from 'react-router-dom';
import { useDropzone } from 'react-dropzone';
import {
  Upload,
  FileText,
  X,
  PlayCircle,
  Settings,
  Database,
  Info,
  CheckCircle,
} from 'lucide-react';

const Analysis: React.FC = () => {
  const navigate = useNavigate();
  const [activeStep, setActiveStep] = useState(0);
  const [analysisType, setAnalysisType] = useState('');
  const [uploadedFiles, setUploadedFiles] = useState<File[]>([]);
  const [parameters, setParameters] = useState<any>({});
  const [referenceGenome, setReferenceGenome] = useState('');
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [projectName, setProjectName] = useState('');
  const [description, setDescription] = useState('');

  const steps = ['Select Analysis Type', 'Upload Data', 'Configure Parameters', 'Review & Submit'];

  const analysisTypes = [
    {
      id: 'rna-seq',
      name: 'RNA-Seq Analysis',
      description: 'Quantify gene expression levels',
      icon: '🧬',
    },
    {
      id: 'wgs',
      name: 'Whole Genome Sequencing',
      description: 'Complete genome analysis',
      icon: '🧬',
    },
    {
      id: 'chip-seq',
      name: 'ChIP-Seq Analysis',
      description: 'Protein-DNA interaction analysis',
      icon: '🔬',
    },
    {
      id: 'variant-calling',
      name: 'Variant Calling',
      description: 'Identify genetic variants',
      icon: '🔍',
    },
  ];

  const onDrop = useCallback((acceptedFiles: File[]) => {
    setUploadedFiles([...uploadedFiles, ...acceptedFiles]);
  }, [uploadedFiles]);

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: {
      'application/x-gzip': ['.gz'],
      'text/plain': ['.fastq', '.fq', '.fa', '.fasta'],
    },
  });

  const handleNext = () => {
    setActiveStep((prevActiveStep) => prevActiveStep + 1);
  };

  const handleBack = () => {
    setActiveStep((prevActiveStep) => prevActiveStep - 1);
  };

  const handleSubmit = async () => {
    setIsSubmitting(true);
    // Simulate submission
    setTimeout(() => {
      navigate('/dashboard');
    }, 2000);
  };

  const removeFile = (index: number) => {
    const newFiles = [...uploadedFiles];
    newFiles.splice(index, 1);
    setUploadedFiles(newFiles);
  };

  const renderStepContent = (step: number) => {
    switch (step) {
      case 0:
        return (
          <Box>
            <Typography variant="h6" sx={{ mb: 3 }}>
              Choose your analysis type
            </Typography>
            <Grid container spacing={3}>
              {analysisTypes.map((type) => (
                <Grid item xs={12} sm={6} key={type.id}>
                  <Card
                    sx={{
                      cursor: 'pointer',
                      border: analysisType === type.id ? '2px solid' : '1px solid',
                      borderColor: analysisType === type.id ? 'primary.main' : 'grey.300',
                      '&:hover': {
                        borderColor: 'primary.main',
                        transform: 'translateY(-2px)',
                        transition: 'all 0.3s',
                      },
                    }}
                    onClick={() => setAnalysisType(type.id)}
                  >
                    <CardContent>
                      <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
                        <Typography variant="h4" sx={{ mr: 2 }}>
                          {type.icon}
                        </Typography>
                        <Typography variant="h6" sx={{ fontWeight: 600 }}>
                          {type.name}
                        </Typography>
                      </Box>
                      <Typography variant="body2" color="text.secondary">
                        {type.description}
                      </Typography>
                    </CardContent>
                  </Card>
                </Grid>
              ))}
            </Grid>
          </Box>
        );

      case 1:
        return (
          <Box>
            <Typography variant="h6" sx={{ mb: 3 }}>
              Upload your data files
            </Typography>
            <Box
              {...getRootProps()}
              sx={{
                border: '2px dashed',
                borderColor: isDragActive ? 'primary.main' : 'grey.400',
                borderRadius: 2,
                p: 4,
                textAlign: 'center',
                cursor: 'pointer',
                backgroundColor: isDragActive ? 'primary.light' : 'grey.50',
                mb: 3,
              }}
            >
              <input {...getInputProps()} />
              <Upload size={48} />
              <Typography variant="h6" sx={{ mt: 2 }}>
                {isDragActive
                  ? 'Drop the files here...'
                  : 'Drag & drop files here, or click to select'}
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
                Supported formats: FASTQ, FASTA, GZ
              </Typography>
            </Box>

            {uploadedFiles.length > 0 && (
              <Box>
                <Typography variant="subtitle1" sx={{ mb: 2 }}>
                  Uploaded Files ({uploadedFiles.length})
                </Typography>
                <List>
                  {uploadedFiles.map((file, index) => (
                    <ListItem
                      key={index}
                      sx={{
                        backgroundColor: 'grey.50',
                        borderRadius: 1,
                        mb: 1,
                      }}
                    >
                      <ListItemIcon>
                        <FileText />
                      </ListItemIcon>
                      <ListItemText
                        primary={file.name}
                        secondary={`${(file.size / 1024 / 1024).toFixed(2)} MB`}
                      />
                      <IconButton onClick={() => removeFile(index)}>
                        <X size={20} />
                      </IconButton>
                    </ListItem>
                  ))}
                </List>
              </Box>
            )}
          </Box>
        );

      case 2:
        return (
          <Box>
            <Typography variant="h6" sx={{ mb: 3 }}>
              Configure analysis parameters
            </Typography>
            <Grid container spacing={3}>
              <Grid item xs={12}>
                <TextField
                  fullWidth
                  label="Project Name"
                  value={projectName}
                  onChange={(e) => setProjectName(e.target.value)}
                />
              </Grid>
              <Grid item xs={12}>
                <TextField
                  fullWidth
                  multiline
                  rows={3}
                  label="Description"
                  value={description}
                  onChange={(e) => setDescription(e.target.value)}
                />
              </Grid>
              <Grid item xs={12} sm={6}>
                <FormControl fullWidth>
                  <InputLabel>Reference Genome</InputLabel>
                  <Select
                    value={referenceGenome}
                    onChange={(e) => setReferenceGenome(e.target.value)}
                    label="Reference Genome"
                  >
                    <MenuItem value="hg38">Human (GRCh38/hg38)</MenuItem>
                    <MenuItem value="hg19">Human (GRCh37/hg19)</MenuItem>
                    <MenuItem value="mm10">Mouse (GRCm38/mm10)</MenuItem>
                    <MenuItem value="dm6">Drosophila (dm6)</MenuItem>
                  </Select>
                </FormControl>
              </Grid>
              <Grid item xs={12} sm={6}>
                <TextField
                  fullWidth
                  type="number"
                  label="Quality Score Threshold"
                  defaultValue="30"
                />
              </Grid>
              <Grid item xs={12} sm={6}>
                <TextField
                  fullWidth
                  type="number"
                  label="Minimum Read Length"
                  defaultValue="50"
                />
              </Grid>
              <Grid item xs={12} sm={6}>
                <FormControl fullWidth>
                  <InputLabel>Adapter Trimming</InputLabel>
                  <Select defaultValue="auto" label="Adapter Trimming">
                    <MenuItem value="auto">Auto-detect</MenuItem>
                    <MenuItem value="illumina">Illumina Universal</MenuItem>
                    <MenuItem value="nextera">Nextera</MenuItem>
                    <MenuItem value="none">None</MenuItem>
                  </Select>
                </FormControl>
              </Grid>
            </Grid>
          </Box>
        );

      case 3:
        return (
          <Box>
            <Typography variant="h6" sx={{ mb: 3 }}>
              Review your analysis configuration
            </Typography>
            <Paper sx={{ p: 3, backgroundColor: 'grey.50' }}>
              <Grid container spacing={2}>
                <Grid item xs={12}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Analysis Type
                  </Typography>
                  <Typography variant="body1" sx={{ fontWeight: 600 }}>
                    {analysisTypes.find((t) => t.id === analysisType)?.name}
                  </Typography>
                </Grid>
                <Grid item xs={12}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Project Name
                  </Typography>
                  <Typography variant="body1" sx={{ fontWeight: 600 }}>
                    {projectName || 'Untitled Project'}
                  </Typography>
                </Grid>
                <Grid item xs={12}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Files
                  </Typography>
                  <Typography variant="body1" sx={{ fontWeight: 600 }}>
                    {uploadedFiles.length} files uploaded
                  </Typography>
                </Grid>
                <Grid item xs={12}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Reference Genome
                  </Typography>
                  <Typography variant="body1" sx={{ fontWeight: 600 }}>
                    {referenceGenome || 'Not selected'}
                  </Typography>
                </Grid>
              </Grid>
            </Paper>
            {isSubmitting && (
              <Box sx={{ mt: 3 }}>
                <Alert severity="info">
                  Submitting your analysis job...
                </Alert>
                <LinearProgress sx={{ mt: 2 }} />
              </Box>
            )}
          </Box>
        );

      default:
        return null;
    }
  };

  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        New Analysis
      </Typography>

      <Paper sx={{ p: 3 }}>
        <Stepper activeStep={activeStep} sx={{ mb: 4 }}>
          {steps.map((label) => (
            <Step key={label}>
              <StepLabel>{label}</StepLabel>
            </Step>
          ))}
        </Stepper>

        <Box sx={{ minHeight: 400 }}>{renderStepContent(activeStep)}</Box>

        <Box sx={{ display: 'flex', justifyContent: 'space-between', mt: 4 }}>
          <Button
            disabled={activeStep === 0}
            onClick={handleBack}
            sx={{ mr: 1 }}
          >
            Back
          </Button>
          <Box>
            {activeStep === steps.length - 1 ? (
              <Button
                variant="contained"
                onClick={handleSubmit}
                startIcon={<PlayCircle />}
                disabled={isSubmitting}
              >
                Submit Analysis
              </Button>
            ) : (
              <Button
                variant="contained"
                onClick={handleNext}
                disabled={
                  (activeStep === 0 && !analysisType) ||
                  (activeStep === 1 && uploadedFiles.length === 0)
                }
              >
                Next
              </Button>
            )}
          </Box>
        </Box>
      </Paper>
    </Box>
  );
};

export default Analysis;
