import React, { useState, useEffect, useContext } from 'react';
import {
  Grid,
  Card,
  CardContent,
  Typography,
  Box,
  Button,
  LinearProgress,
  Chip,
  Avatar,
  IconButton,
  Menu,
  MenuItem,
  Paper,
  List,
  ListItem,
  ListItemAvatar,
  ListItemText,
  ListItemSecondaryAction,
} from '@mui/material';
import { useNavigate } from 'react-router-dom';
import {
  Activity,
  FileText,
  Database,
  TrendingUp,
  MoreVertical,
  PlayCircle,
  Clock,
  CheckCircle,
  AlertCircle,
  ArrowUpRight,
  ArrowDownRight,
  Folder,
  TestTube2,
  Dna,
} from 'lucide-react';
import { UserContext } from '../App';

const Dashboard: React.FC = () => {
  const navigate = useNavigate();
  const { user } = useContext(UserContext);
  const [loading, setLoading] = useState(false);
  const [stats, setStats] = useState({
    totalAnalyses: 156,
    activeJobs: 3,
    completedToday: 12,
    storageUsed: 75,
  });

  const [recentAnalyses, setRecentAnalyses] = useState([
    {
      id: 1,
      name: 'RNA-Seq Analysis - Sample_001',
      status: 'completed',
      progress: 100,
      time: '2 hours ago',
      type: 'RNA-Seq',
    },
    {
      id: 2,
      name: 'Variant Calling - Patient_XY123',
      status: 'running',
      progress: 65,
      time: 'Started 30 min ago',
      type: 'WGS',
    },
    {
      id: 3,
      name: 'ChIP-Seq Peak Detection',
      status: 'queued',
      progress: 0,
      time: 'Queued',
      type: 'ChIP-Seq',
    },
  ]);

  const [recentProjects, setRecentProjects] = useState([
    {
      id: 1,
      name: 'Cancer Genomics Study',
      samples: 45,
      lastModified: '1 day ago',
      status: 'active',
    },
    {
      id: 2,
      name: 'COVID-19 Variant Analysis',
      samples: 128,
      lastModified: '3 days ago',
      status: 'active',
    },
    {
      id: 3,
      name: 'Microbiome Diversity Project',
      samples: 67,
      lastModified: '1 week ago',
      status: 'completed',
    },
  ]);

  useEffect(() => {
    // Fetch dashboard data
    fetchDashboardData();
  }, []);

  const fetchDashboardData = async () => {
    setLoading(true);
    // Simulate API call
    setTimeout(() => {
      setLoading(false);
    }, 1000);
  };

  const getStatusColor = (status: string): "success" | "primary" | "warning" | "error" | "default" | "secondary" | "info" => {
    switch (status) {
      case 'completed':
        return 'success';
      case 'running':
        return 'primary';
      case 'queued':
        return 'warning';
      case 'failed':
        return 'error';
      default:
        return 'default';
    }
  };

  const getStatusIcon = (status: string) => {
    switch (status) {
      case 'completed':
        return <CheckCircle size={16} />;
      case 'running':
        return <PlayCircle size={16} />;
      case 'queued':
        return <Clock size={16} />;
      case 'failed':
        return <AlertCircle size={16} />;
      default:
        return undefined;
    }
  };

  return (
    <Box>
      {/* Welcome Section */}
      <Box sx={{ mb: 4 }}>
        <Typography variant="h4" sx={{ fontWeight: 700, mb: 1 }}>
          Welcome back, {user?.name || 'User'}!
        </Typography>
        <Typography variant="body1" color="text.secondary">
          Here's an overview of your genomic analysis platform
        </Typography>
      </Box>

      {/* Stats Cards */}
      <Grid container spacing={3} sx={{ mb: 4 }}>
        <Grid item xs={12} sm={6} md={3}>
          <Card>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
                <Avatar sx={{ bgcolor: 'primary.light', mr: 2 }}>
                  <Activity size={24} />
                </Avatar>
                <Box>
                  <Typography color="text.secondary" variant="body2">
                    Total Analyses
                  </Typography>
                  <Typography variant="h5" sx={{ fontWeight: 700 }}>
                    {stats.totalAnalyses}
                  </Typography>
                </Box>
              </Box>
              <Box sx={{ display: 'flex', alignItems: 'center' }}>
                <ArrowUpRight size={16} color="#10B981" />
                <Typography variant="body2" color="success.main" sx={{ ml: 0.5 }}>
                  +12% from last month
                </Typography>
              </Box>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} sm={6} md={3}>
          <Card>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
                <Avatar sx={{ bgcolor: 'warning.light', mr: 2 }}>
                  <PlayCircle size={24} />
                </Avatar>
                <Box>
                  <Typography color="text.secondary" variant="body2">
                    Active Jobs
                  </Typography>
                  <Typography variant="h5" sx={{ fontWeight: 700 }}>
                    {stats.activeJobs}
                  </Typography>
                </Box>
              </Box>
              <Chip
                label="View all"
                size="small"
                onClick={() => navigate('/analysis')}
              />
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} sm={6} md={3}>
          <Card>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
                <Avatar sx={{ bgcolor: 'success.light', mr: 2 }}>
                  <CheckCircle size={24} />
                </Avatar>
                <Box>
                  <Typography color="text.secondary" variant="body2">
                    Completed Today
                  </Typography>
                  <Typography variant="h5" sx={{ fontWeight: 700 }}>
                    {stats.completedToday}
                  </Typography>
                </Box>
              </Box>
              <Box sx={{ display: 'flex', alignItems: 'center' }}>
                <ArrowUpRight size={16} color="#10B981" />
                <Typography variant="body2" color="success.main" sx={{ ml: 0.5 }}>
                  Above average
                </Typography>
              </Box>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} sm={6} md={3}>
          <Card>
            <CardContent>
              <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
                <Avatar sx={{ bgcolor: 'info.light', mr: 2 }}>
                  <Database size={24} />
                </Avatar>
                <Box>
                  <Typography color="text.secondary" variant="body2">
                    Storage Used
                  </Typography>
                  <Typography variant="h5" sx={{ fontWeight: 700 }}>
                    {stats.storageUsed}%
                  </Typography>
                </Box>
              </Box>
              <LinearProgress
                variant="determinate"
                value={stats.storageUsed}
                sx={{ height: 6, borderRadius: 3 }}
              />
            </CardContent>
          </Card>
        </Grid>
      </Grid>

      {/* Main Content Grid */}
      <Grid container spacing={3}>
        {/* Recent Analyses */}
        <Grid item xs={12} md={8}>
          <Paper sx={{ p: 3 }}>
            <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 3 }}>
              <Typography variant="h6" sx={{ fontWeight: 600 }}>
                Recent Analyses
              </Typography>
              <Button
                variant="text"
                size="small"
                onClick={() => navigate('/reports')}
              >
                View All
              </Button>
            </Box>
            <List>
              {recentAnalyses.map((analysis) => {
                const statusIcon = getStatusIcon(analysis.status);
                return (
                  <ListItem
                    key={analysis.id}
                    sx={{
                      mb: 2,
                      backgroundColor: 'grey.50',
                      borderRadius: 1,
                      '&:hover': {
                        backgroundColor: 'grey.100',
                        cursor: 'pointer',
                      },
                    }}
                    onClick={() => navigate(`/results/${analysis.id}`)}
                  >
                    <ListItemAvatar>
                      <Avatar sx={{ bgcolor: `${getStatusColor(analysis.status)}.light` }}>
                        <Dna size={24} />
                      </Avatar>
                    </ListItemAvatar>
                    <ListItemText
                      primary={
                        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                          <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
                            {analysis.name}
                          </Typography>
                          <Chip
                            label={analysis.type}
                            size="small"
                            variant="outlined"
                          />
                          {statusIcon ? (
                            <Chip
                              icon={statusIcon}
                              label={analysis.status}
                              size="small"
                              color={getStatusColor(analysis.status)}
                            />
                          ) : (
                            <Chip
                              label={analysis.status}
                              size="small"
                              color={getStatusColor(analysis.status)}
                            />
                          )}
                        </Box>
                      }
                      secondary={
                        <Box>
                          <Typography variant="body2" color="text.secondary">
                            {analysis.time}
                          </Typography>
                          {analysis.status === 'running' && (
                            <LinearProgress
                              variant="determinate"
                              value={analysis.progress}
                              sx={{ mt: 1, height: 4, borderRadius: 2 }}
                            />
                          )}
                        </Box>
                      }
                    />
                    <ListItemSecondaryAction>
                      <IconButton edge="end">
                        <MoreVertical size={20} />
                      </IconButton>
                    </ListItemSecondaryAction>
                  </ListItem>
                );
              })}
            </List>
            <Box sx={{ textAlign: 'center', mt: 2 }}>
              <Button
                variant="contained"
                startIcon={<PlayCircle />}
                onClick={() => navigate('/analysis')}
              >
                Start New Analysis
              </Button>
            </Box>
          </Paper>
        </Grid>

        {/* Recent Projects */}
        <Grid item xs={12} md={4}>
          <Paper sx={{ p: 3 }}>
            <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 3 }}>
              <Typography variant="h6" sx={{ fontWeight: 600 }}>
                Recent Projects
              </Typography>
              <Button
                variant="text"
                size="small"
                onClick={() => navigate('/projects')}
              >
                View All
              </Button>
            </Box>
            <List>
              {recentProjects.map((project) => (
                <ListItem
                  key={project.id}
                  sx={{
                    mb: 2,
                    backgroundColor: 'grey.50',
                    borderRadius: 1,
                    '&:hover': {
                      backgroundColor: 'grey.100',
                      cursor: 'pointer',
                    },
                  }}
                  onClick={() => navigate(`/projects/${project.id}`)}
                >
                  <ListItemAvatar>
                    <Avatar sx={{ bgcolor: 'primary.light' }}>
                      <Folder size={24} />
                    </Avatar>
                  </ListItemAvatar>
                  <ListItemText
                    primary={
                      <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>
                        {project.name}
                      </Typography>
                    }
                    secondary={
                      <Box>
                        <Typography variant="caption" color="text.secondary">
                          {project.samples} samples • {project.lastModified}
                        </Typography>
                      </Box>
                    }
                  />
                  <Chip
                    label={project.status}
                    size="small"
                    color={project.status === 'active' ? 'success' : 'default'}
                  />
                </ListItem>
              ))}
            </List>
          </Paper>
        </Grid>
      </Grid>
    </Box>
  );
};

export default Dashboard;
