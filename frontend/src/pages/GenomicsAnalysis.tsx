import React, { useState, useEffect, useRef } from 'react';
import {
  Box,
  Card,
  CardContent,
  Typography,
  TextField,
  Button,
  Alert,
  CircularProgress,
  Chip,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
  Grid,
  IconButton,
  Tooltip,
  Dialog,
  DialogTitle,
  DialogContent,
  Tab,
  Tabs,
  List,
  ListItem,
  ListItemText,
  Divider,
} from '@mui/material';
import {
  Search,
  Dna,
  AlertCircle,
  CheckCircle,
  HelpCircle,
  Eye,
  Download,
  Filter,
  Info,
  Upload,
} from 'lucide-react';
import axios from 'axios';

// VCF Analysis Component - Improved
const VCFAnalysis: React.FC = () => {
  const [vcfFile, setVcfFile] = useState<File | null>(null);
  const [loading, setLoading] = useState(false);
  const [vcfStats, setVcfStats] = useState<any>(null);
  const [annotations, setAnnotations] = useState<any>(null);
  const [error, setError] = useState<string | null>(null);

  const handleVcfUpload = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    // Check file extension
    if (!file.name.endsWith('.vcf') && !file.name.endsWith('.vcf.gz')) {
      setError('Please upload a VCF file (.vcf or .vcf.gz)');
      return;
    }

    setVcfFile(file);
    setLoading(true);
    setError(null);

    const formData = new FormData();
    formData.append('file', file);

    try {
      // Upload VCF
      const uploadResponse = await axios.post('http://localhost:8000/api/genomics/upload/vcf', formData, {
        headers: { 'Content-Type': 'multipart/form-data' }
      });
      
      setVcfStats(uploadResponse.data.statistics);
      
      // Automatically annotate the uploaded file
      if (uploadResponse.data.file) {
        const annotResponse = await axios.post('http://localhost:8000/api/genomics/annotate', {
          vcf_path: uploadResponse.data.file,
        });
        setAnnotations(annotResponse.data);
      }
    } catch (error: any) {
      console.error('Upload failed:', error);
      setError(error.response?.data?.detail || 'Failed to process VCF file');
    } finally {
      setLoading(false);
    }
  };

  return (
    <Card>
      <CardContent>
        <Typography variant="h6" gutterBottom>
          VCF File Analysis
        </Typography>
        
        <Box sx={{ mb: 3 }}>
          <input
            type="file"
            accept=".vcf,.vcf.gz"
            onChange={handleVcfUpload}
            style={{ display: 'none' }}
            id="vcf-upload"
          />
          <label htmlFor="vcf-upload">
            <Button
              variant="contained"
              component="span"
              disabled={loading}
              startIcon={loading ? <CircularProgress size={20} /> : <Upload />}
            >
              {loading ? 'Processing...' : 'Upload VCF File'}
            </Button>
          </label>
          
          {vcfFile && (
            <Typography variant="body2" sx={{ mt: 1 }}>
              File: {vcfFile.name}
            </Typography>
          )}
        </Box>

        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        {vcfStats && (
          <Paper sx={{ p: 2, mb: 2 }}>
            <Typography variant="subtitle1" gutterBottom>
              📊 VCF Statistics
            </Typography>
            <Grid container spacing={2}>
              <Grid item xs={6} md={3}>
                <Typography variant="caption" color="text.secondary">
                  Total Variants
                </Typography>
                <Typography variant="h6">{vcfStats.total_variants || 0}</Typography>
              </Grid>
              <Grid item xs={6} md={3}>
                <Typography variant="caption" color="text.secondary">
                  Samples
                </Typography>
                <Typography variant="h6">{vcfStats.samples?.length || 0}</Typography>
              </Grid>
              <Grid item xs={12} md={6}>
                <Typography variant="caption" color="text.secondary">
                  Variant Types
                </Typography>
                <Box sx={{ mt: 1 }}>
                  {Object.entries(vcfStats.variant_types || {}).map(([type, count]) => (
                    <Chip
                      key={type}
                      label={`${type}: ${count}`}
                      size="small"
                      sx={{ mr: 1, mb: 1 }}
                      color={type === 'SNP' ? 'primary' : 'default'}
                    />
                  ))}
                </Box>
              </Grid>
            </Grid>
          </Paper>
        )}

        {annotations && (
          <Paper sx={{ p: 2 }}>
            <Typography variant="subtitle1" gutterBottom>
              🧬 Annotation Results
            </Typography>
            <Grid container spacing={2} sx={{ mb: 2 }}>
              <Grid item xs={4}>
                <Alert severity="info" icon={<Info />}>
                  <Typography variant="body2">
                    ClinVar Matches: <strong>{annotations.clinvar_matches}</strong>
                  </Typography>
                </Alert>
              </Grid>
              <Grid item xs={4}>
                <Alert severity="success" icon={<CheckCircle />}>
                  <Typography variant="body2">
                    dbSNP Matches: <strong>{annotations.dbsnp_matches}</strong>
                  </Typography>
                </Alert>
              </Grid>
              <Grid item xs={4}>
                <Alert severity="error" icon={<AlertCircle />}>
                  <Typography variant="body2">
                    Pathogenic: <strong>{annotations.pathogenic_count}</strong>
                  </Typography>
                </Alert>
              </Grid>
            </Grid>

            {annotations.annotated_variants?.length > 0 && (
              <>
                <Typography variant="subtitle2" gutterBottom>
                  Annotated Variants (First 10)
                </Typography>
                <TableContainer component={Paper} variant="outlined">
                  <Table size="small">
                    <TableHead>
                      <TableRow>
                        <TableCell>Position</TableCell>
                        <TableCell>Change</TableCell>
                        <TableCell>Type</TableCell>
                        <TableCell>Quality</TableCell>
                        <TableCell>Clinical Significance</TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {annotations.annotated_variants.slice(0, 10).map((item: any, idx: number) => (
                        <TableRow key={idx}>
                          <TableCell>
                            {item.variant.chromosome}:{item.variant.position}
                          </TableCell>
                          <TableCell>
                            <Chip 
                              label={`${item.variant.ref}→${item.variant.alt}`}
                              size="small"
                              variant="outlined"
                            />
                          </TableCell>
                          <TableCell>{item.variant.type || 'SNP'}</TableCell>
                          <TableCell>{item.variant.quality?.toFixed(1) || '-'}</TableCell>
                          <TableCell>
                            {item.clinvar?.clinical_significance ? (
                              <Chip
                                size="small"
                                label={item.clinvar.clinical_significance}
                                color={
                                  item.clinvar.clinical_significance.toLowerCase().includes('pathogenic')
                                    ? 'error'
                                    : item.clinvar.clinical_significance.toLowerCase().includes('benign')
                                    ? 'success'
                                    : 'warning'
                                }
                              />
                            ) : (
                              '-'
                            )}
                          </TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                </TableContainer>
              </>
            )}
          </Paper>
        )}
      </CardContent>
    </Card>
  );
};

// ClinVar Search Component - Improved
const ClinVarSearch: React.FC = () => {
  const [searchTerm, setSearchTerm] = useState('');
  const [searchType, setSearchType] = useState<'rsid' | 'position'>('rsid');
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<any>(null);
  const [error, setError] = useState<string | null>(null);

  // Example searches for user guidance
  const examples = {
    rsid: ['rs80357906', 'rs121913343', 'rs1799966'],
    position: ['17:43044295 G>A', '7:117559590 ATCT>A', '13:32911888 A>G']
  };

  const handleSearch = async () => {
    if (!searchTerm) {
      setError('Please enter a search term');
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null);

    try {
      let response;
      if (searchType === 'rsid') {
        // Clean up rsID
        const rsid = searchTerm.startsWith('rs') ? searchTerm : `rs${searchTerm}`;
        response = await axios.get(`http://localhost:8000/api/genomics/clinvar/${rsid}`);
      } else {
        // Parse position format: chr:position ref>alt
        const match = searchTerm.match(/(\w+):(\d+)\s+([ACGT]+)>([ACGT]+)/i);
        if (!match) {
          throw new Error('Invalid format. Use: chr:position ref>alt (e.g., 17:43044295 G>A)');
        }
        const [, chr, pos, ref, alt] = match;
        response = await axios.post('http://localhost:8000/api/genomics/variants/query', {
          chromosome: chr,
          position: parseInt(pos),
          ref,
          alt,
        });
      }
      setResult(response.data);
    } catch (err: any) {
      setError(err.response?.data?.detail || err.message || 'Search failed');
    } finally {
      setLoading(false);
    }
  };

  const handleExampleClick = (example: string) => {
    setSearchTerm(example);
  };

  const getSignificanceColor = (significance: string) => {
    const sig = significance?.toLowerCase() || '';
    if (sig.includes('pathogenic')) return 'error';
    if (sig.includes('benign')) return 'success';
    if (sig.includes('uncertain')) return 'warning';
    return 'default';
  };

  const getSignificanceIcon = (significance: string) => {
    const sig = significance?.toLowerCase() || '';
    if (sig.includes('pathogenic')) return <AlertCircle size={16} />;
    if (sig.includes('benign')) return <CheckCircle size={16} />;
    return <HelpCircle size={16} />;
  };

  return (
    <Card>
      <CardContent>
        <Typography variant="h6" gutterBottom>
          ClinVar Database Search
        </Typography>
        
        <Alert severity="info" sx={{ mb: 2 }}>
          <Typography variant="body2">
            Search for genetic variants in ClinVar to find clinical significance and associated conditions.
          </Typography>
        </Alert>

        <Tabs value={searchType} onChange={(e, v) => setSearchType(v)} sx={{ mb: 2 }}>
          <Tab label="Search by rsID" value="rsid" />
          <Tab label="Search by Position" value="position" />
        </Tabs>

        <Box sx={{ mb: 2 }}>
          <TextField
            fullWidth
            label={searchType === 'rsid' ? 'Enter rsID (e.g., rs80357906)' : 'Enter position (e.g., 17:43044295 G>A)'}
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
            onKeyPress={(e) => e.key === 'Enter' && handleSearch()}
            sx={{ mb: 1 }}
          />
          
          <Box sx={{ mb: 2 }}>
            <Typography variant="caption" color="text.secondary">
              Try these examples:
            </Typography>
            <Box sx={{ mt: 0.5 }}>
              {examples[searchType].map((example) => (
                <Chip
                  key={example}
                  label={example}
                  size="small"
                  onClick={() => handleExampleClick(example)}
                  sx={{ mr: 1, mb: 1, cursor: 'pointer' }}
                  variant="outlined"
                />
              ))}
            </Box>
          </Box>

          <Button
            variant="contained"
            onClick={handleSearch}
            disabled={loading || !searchTerm}
            startIcon={loading ? <CircularProgress size={20} /> : <Search />}
            fullWidth
          >
            {loading ? 'Searching...' : 'Search ClinVar'}
          </Button>
        </Box>

        {error && (
          <Alert severity="error" sx={{ mb: 2 }}>
            {error}
          </Alert>
        )}

        {result && (
          <Box>
            <Paper sx={{ p: 2, mb: 2 }}>
              <Grid container spacing={2}>
                <Grid item xs={12}>
                  <Typography variant="subtitle1" gutterBottom>
                    Search Results
                  </Typography>
                </Grid>
                <Grid item xs={12} md={6}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Clinical Significance
                  </Typography>
                  <Chip
                    icon={getSignificanceIcon(result.clinical_significance)}
                    label={result.clinical_significance || 'Not specified'}
                    color={getSignificanceColor(result.clinical_significance)}
                    sx={{ mt: 1 }}
                  />
                </Grid>
                <Grid item xs={12} md={6}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Gene
                  </Typography>
                  <Typography variant="h6">{result.gene || 'N/A'}</Typography>
                </Grid>
                <Grid item xs={12}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Associated Conditions
                  </Typography>
                  <Box sx={{ mt: 1 }}>
                    {result.conditions?.length > 0 ? (
                      result.conditions.map((condition: string, idx: number) => (
                        <Chip 
                          key={idx} 
                          label={condition} 
                          size="small" 
                          sx={{ mr: 1, mb: 1 }}
                          variant="outlined"
                        />
                      ))
                    ) : (
                      <Typography variant="body2">No conditions listed</Typography>
                    )}
                  </Box>
                </Grid>
                <Grid item xs={12}>
                  <Typography variant="subtitle2" color="text.secondary">
                    Review Status
                  </Typography>
                  <Typography variant="body2">{result.review_status || 'Not reviewed'}</Typography>
                </Grid>
                {result.url && (
                  <Grid item xs={12}>
                    <Button
                      variant="outlined"
                      href={result.url}
                      target="_blank"
                      startIcon={<Eye />}
                      size="small"
                    >
                      View in ClinVar Database
                    </Button>
                  </Grid>
                )}
              </Grid>
            </Paper>

            {result.interpretation && (
              <Paper sx={{ p: 2, backgroundColor: 'background.default' }}>
                <Typography variant="subtitle2" gutterBottom>
                  Clinical Interpretation
                </Typography>
                <List dense>
                  <ListItem>
                    <ListItemText 
                      primary="Category"
                      secondary={result.interpretation.category}
                    />
                  </ListItem>
                  <ListItem>
                    <ListItemText 
                      primary="Actionable"
                      secondary={result.interpretation.actionable ? 'Yes - Medical action may be needed' : 'No'}
                    />
                  </ListItem>
                  <ListItem>
                    <ListItemText 
                      primary="Follow-up Required"
                      secondary={result.interpretation.requires_follow_up ? 'Yes - Consult genetics specialist' : 'No'}
                    />
                  </ListItem>
                </List>
              </Paper>
            )}
          </Box>
        )}
      </CardContent>
    </Card>
  );
};

// Simple Genome Browser Component (without external dependencies)
const SimpleGenomeBrowser: React.FC = () => {
  const [chromosome, setChromosome] = useState('1');
  const [position, setPosition] = useState('1000000');
  
  const chromosomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'];
  
  return (
    <Card>
      <CardContent>
        <Typography variant="h6" gutterBottom>
          Genome Browser
        </Typography>
        
        <Alert severity="info" sx={{ mb: 2 }}>
          <Typography variant="body2">
            Navigate through the human genome. Enter a chromosome and position to explore.
          </Typography>
        </Alert>
        
        <Grid container spacing={2}>
          <Grid item xs={12} md={4}>
            <TextField
              select
              fullWidth
              label="Chromosome"
              value={chromosome}
              onChange={(e) => setChromosome(e.target.value)}
              SelectProps={{ native: true }}
            >
              {chromosomes.map(chr => (
                <option key={chr} value={chr}>chr{chr}</option>
              ))}
            </TextField>
          </Grid>
          <Grid item xs={12} md={4}>
            <TextField
              fullWidth
              label="Position"
              value={position}
              onChange={(e) => setPosition(e.target.value)}
              type="number"
            />
          </Grid>
          <Grid item xs={12} md={4}>
            <Button 
              variant="contained" 
              fullWidth 
              sx={{ height: '56px' }}
              startIcon={<Search />}
            >
              Navigate
            </Button>
          </Grid>
        </Grid>
        
        <Box sx={{ mt: 3, p: 2, bgcolor: 'background.default', borderRadius: 1 }}>
          <Typography variant="body2" color="text.secondary">
            Current Location: <strong>chr{chromosome}:{parseInt(position).toLocaleString()}</strong>
          </Typography>
          <Typography variant="caption" color="text.secondary">
            Note: For full genome visualization, external genome browser tools like UCSC Genome Browser or Ensembl are recommended.
          </Typography>
        </Box>
        
        <Box sx={{ mt: 2 }}>
          <Typography variant="subtitle2" gutterBottom>
            Quick Links to External Browsers:
          </Typography>
          <Button
            size="small"
            variant="outlined"
            sx={{ mr: 1 }}
            href={`https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr${chromosome}:${position}`}
            target="_blank"
          >
            UCSC Browser
          </Button>
          <Button
            size="small"
            variant="outlined"
            href={`https://www.ensembl.org/Homo_sapiens/Location/View?r=${chromosome}:${position}`}
            target="_blank"
          >
            Ensembl
          </Button>
        </Box>
      </CardContent>
    </Card>
  );
};

// Main Genomics Analysis Page
const GenomicsAnalysis: React.FC = () => {
  const [activeTab, setActiveTab] = useState(0);

  return (
    <Box>
      <Typography variant="h4" gutterBottom>
        Genomics Analysis Platform
      </Typography>

      <Box sx={{ borderBottom: 1, borderColor: 'divider', mb: 3 }}>
        <Tabs value={activeTab} onChange={(e, v) => setActiveTab(v)}>
          <Tab label="ClinVar Search" icon={<Search />} />
          <Tab label="VCF Analysis" icon={<Dna />} />
          <Tab label="Genome Browser" icon={<Eye />} />
        </Tabs>
      </Box>

      {activeTab === 0 && <ClinVarSearch />}
      {activeTab === 1 && <VCFAnalysis />}
      {activeTab === 2 && <SimpleGenomeBrowser />}
    </Box>
  );
};

export default GenomicsAnalysis;
