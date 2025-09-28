import React, { useState, useEffect } from 'react';
import {
  Box,
  Card,
  CardContent,
  Typography,
  Grid,
  Button,
  Chip,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  Alert,
  CircularProgress,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  TextField,
  Paper,
  LinearProgress,
  Divider,
} from '@mui/material';
import {
  Check,
  X,
  CreditCard,
  TrendingUp,
  Package,
  Zap,
  Shield,
  Users,
  Star,
  AlertCircle,
  ChevronRight,
} from 'lucide-react';
import axios from 'axios';

interface Plan {
  id: string;
  name: string;
  price: number | string;
  price_display: string;
  recommended?: boolean;
  features: {
    analyses_per_month: number | string;
    storage: string;
    api_calls: string;
    max_file_size: string;
    features_list: string[];
  };
}

interface Usage {
  month: string;
  analyses_count: number;
  storage_mb: number;
  api_calls: number;
  estimated_cost: number;
}

const Billing: React.FC = () => {
  const [plans, setPlans] = useState<Plan[]>([]);
  const [currentPlan, setCurrentPlan] = useState<string>('free');
  const [usage, setUsage] = useState<Usage | null>(null);
  const [loading, setLoading] = useState(true);
  const [upgradeDialog, setUpgradeDialog] = useState(false);
  const [selectedPlan, setSelectedPlan] = useState<Plan | null>(null);
  const [paymentProcessing, setPaymentProcessing] = useState(false);

  useEffect(() => {
    fetchPlansAndUsage();
  }, []);

  const fetchPlansAndUsage = async () => {
    try {
      // Fetch available plans
      const plansResponse = await axios.get('http://localhost:8000/api/billing/plans');
      setPlans(plansResponse.data.plans);

      // Fetch usage (mock user ID for demo)
      const usageResponse = await axios.get('http://localhost:8000/api/billing/usage?user_id=demo_user');
      setUsage(usageResponse.data.usage);
      setCurrentPlan(usageResponse.data.plan?.toLowerCase() || 'free');
    } catch (error) {
      console.error('Failed to fetch billing data:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleUpgrade = (plan: Plan) => {
    setSelectedPlan(plan);
    setUpgradeDialog(true);
  };

  const processUpgrade = async () => {
    if (!selectedPlan) return;

    setPaymentProcessing(true);
    try {
      const response = await axios.post('http://localhost:8000/api/billing/subscribe', {
        plan_id: selectedPlan.id,
        user_email: 'user@example.com',
        user_name: 'John Doe',
      });

      if (response.data.status === 'success') {
        alert(`Successfully upgraded to ${selectedPlan.name} plan!`);
        setCurrentPlan(selectedPlan.id);
        setUpgradeDialog(false);
      } else if (response.data.status === 'contact_sales') {
        alert('Please contact our sales team for enterprise pricing.');
      }
    } catch (error) {
      console.error('Upgrade failed:', error);
      alert('Failed to process upgrade. Please try again.');
    } finally {
      setPaymentProcessing(false);
    }
  };

  const getPlanColor = (planId: string) => {
    switch (planId) {
      case 'free':
        return 'default';
      case 'professional':
        return 'primary';
      case 'team':
        return 'secondary';
      case 'enterprise':
        return 'warning';
      default:
        return 'default';
    }
  };

  const getFeatureIcon = (feature: string) => {
    if (feature.includes('analysis') || feature.includes('analyses')) return <TrendingUp size={16} />;
    if (feature.includes('storage')) return <Package size={16} />;
    if (feature.includes('api')) return <Zap size={16} />;
    if (feature.includes('team') || feature.includes('collaboration')) return <Users size={16} />;
    if (feature.includes('support') || feature.includes('priority')) return <Shield size={16} />;
    return <Check size={16} />;
  };

  if (loading) {
    return (
      <Box sx={{ display: 'flex', justifyContent: 'center', p: 4 }}>
        <CircularProgress />
      </Box>
    );
  }

  return (
    <Box>
      <Typography variant="h4" gutterBottom>
        Billing & Subscription
      </Typography>

      {/* Current Usage */}
      {usage && (
        <Card sx={{ mb: 4 }}>
          <CardContent>
            <Typography variant="h6" gutterBottom>
              Current Usage - {usage.month}
            </Typography>
            
            <Grid container spacing={3}>
              <Grid item xs={12} md={3}>
                <Box>
                  <Typography variant="caption" color="text.secondary">
                    Analyses This Month
                  </Typography>
                  <Typography variant="h4">
                    {usage.analyses_count}
                  </Typography>
                  <LinearProgress 
                    variant="determinate" 
                    value={Math.min((usage.analyses_count / 50) * 100, 100)}
                    sx={{ mt: 1 }}
                  />
                </Box>
              </Grid>
              
              <Grid item xs={12} md={3}>
                <Box>
                  <Typography variant="caption" color="text.secondary">
                    Storage Used
                  </Typography>
                  <Typography variant="h4">
                    {(usage.storage_mb / 1024).toFixed(1)} GB
                  </Typography>
                  <LinearProgress 
                    variant="determinate" 
                    value={Math.min((usage.storage_mb / 102400) * 100, 100)}
                    sx={{ mt: 1 }}
                  />
                </Box>
              </Grid>
              
              <Grid item xs={12} md={3}>
                <Box>
                  <Typography variant="caption" color="text.secondary">
                    API Calls
                  </Typography>
                  <Typography variant="h4">
                    {usage.api_calls.toLocaleString()}
                  </Typography>
                  <LinearProgress 
                    variant="determinate" 
                    value={Math.min((usage.api_calls / 10000) * 100, 100)}
                    sx={{ mt: 1 }}
                  />
                </Box>
              </Grid>
              
              <Grid item xs={12} md={3}>
                <Box>
                  <Typography variant="caption" color="text.secondary">
                    Estimated Additional Charges
                  </Typography>
                  <Typography variant="h4">
                    ${usage.estimated_cost}
                  </Typography>
                  <Typography variant="caption" color="text.secondary">
                    For usage over plan limits
                  </Typography>
                </Box>
              </Grid>
            </Grid>
          </CardContent>
        </Card>
      )}

      {/* Subscription Plans */}
      <Typography variant="h5" gutterBottom sx={{ mb: 3 }}>
        Available Plans
      </Typography>

      <Grid container spacing={3}>
        {plans.map((plan) => (
          <Grid item xs={12} md={6} lg={3} key={plan.id}>
            <Card 
              sx={{ 
                height: '100%',
                position: 'relative',
                border: currentPlan === plan.id ? 2 : 0,
                borderColor: 'primary.main',
              }}
            >
              {plan.recommended && (
                <Chip
                  label="RECOMMENDED"
                  color="primary"
                  size="small"
                  sx={{
                    position: 'absolute',
                    top: 10,
                    right: 10,
                  }}
                  icon={<Star size={14} />}
                />
              )}
              
              {currentPlan === plan.id && (
                <Chip
                  label="CURRENT PLAN"
                  color="success"
                  size="small"
                  sx={{
                    position: 'absolute',
                    top: 10,
                    left: 10,
                  }}
                />
              )}
              
              <CardContent>
                <Typography variant="h5" gutterBottom sx={{ mt: 2 }}>
                  {plan.name}
                </Typography>
                
                <Typography variant="h3" gutterBottom>
                  {plan.price_display}
                </Typography>
                
                <Divider sx={{ my: 2 }} />
                
                <Typography variant="subtitle2" gutterBottom>
                  Features:
                </Typography>
                
                <List dense>
                  <ListItem>
                    <ListItemIcon sx={{ minWidth: 32 }}>
                      <TrendingUp size={16} />
                    </ListItemIcon>
                    <ListItemText 
                      primary={`${plan.features.analyses_per_month} analyses/month`}
                    />
                  </ListItem>
                  
                  <ListItem>
                    <ListItemIcon sx={{ minWidth: 32 }}>
                      <Package size={16} />
                    </ListItemIcon>
                    <ListItemText 
                      primary={`${plan.features.storage} storage`}
                    />
                  </ListItem>
                  
                  <ListItem>
                    <ListItemIcon sx={{ minWidth: 32 }}>
                      <Zap size={16} />
                    </ListItemIcon>
                    <ListItemText 
                      primary={`${plan.features.api_calls} API calls`}
                    />
                  </ListItem>
                  
                  {plan.features.features_list.slice(0, 3).map((feature, idx) => (
                    <ListItem key={idx}>
                      <ListItemIcon sx={{ minWidth: 32 }}>
                        {getFeatureIcon(feature)}
                      </ListItemIcon>
                      <ListItemText 
                        primary={feature.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}
                      />
                    </ListItem>
                  ))}
                </List>
                
                {plan.features.features_list.length > 3 && (
                  <Typography variant="caption" color="text.secondary">
                    +{plan.features.features_list.length - 3} more features
                  </Typography>
                )}
                
                <Box sx={{ mt: 3 }}>
                  {currentPlan === plan.id ? (
                    <Button 
                      fullWidth 
                      variant="outlined" 
                      disabled
                    >
                      Current Plan
                    </Button>
                  ) : (
                    <Button 
                      fullWidth 
                      variant={plan.recommended ? "contained" : "outlined"}
                      onClick={() => handleUpgrade(plan)}
                      endIcon={<ChevronRight />}
                    >
                      {plan.id === 'free' ? 'Downgrade' : 
                       plan.id === 'enterprise' ? 'Contact Sales' : 'Upgrade'}
                    </Button>
                  )}
                </Box>
              </CardContent>
            </Card>
          </Grid>
        ))}
      </Grid>

      {/* Billing History */}
      <Card sx={{ mt: 4 }}>
        <CardContent>
          <Typography variant="h6" gutterBottom>
            Billing History
          </Typography>
          
          <Alert severity="info" sx={{ mb: 2 }}>
            Your next billing date is <strong>January 1, 2025</strong>
          </Alert>
          
          <List>
            <ListItem>
              <ListItemText 
                primary="December 1, 2024"
                secondary="Professional Plan - $99.00"
              />
              <Button size="small">Download Invoice</Button>
            </ListItem>
            <Divider />
            <ListItem>
              <ListItemText 
                primary="November 1, 2024"
                secondary="Professional Plan - $99.00"
              />
              <Button size="small">Download Invoice</Button>
            </ListItem>
          </List>
          
          <Box sx={{ mt: 2 }}>
            <Button 
              variant="outlined" 
              startIcon={<CreditCard />}
            >
              Manage Payment Methods
            </Button>
          </Box>
        </CardContent>
      </Card>

      {/* Upgrade Dialog */}
      <Dialog open={upgradeDialog} onClose={() => setUpgradeDialog(false)} maxWidth="sm" fullWidth>
        <DialogTitle>
          Upgrade to {selectedPlan?.name} Plan
        </DialogTitle>
        <DialogContent>
          <Alert severity="info" sx={{ mb: 2 }}>
            <Typography variant="body2">
              You're upgrading from <strong>{currentPlan}</strong> to <strong>{selectedPlan?.name}</strong> plan.
            </Typography>
          </Alert>
          
          {selectedPlan?.id === 'enterprise' ? (
            <Box>
              <Typography variant="body1" gutterBottom>
                Enterprise plan offers custom pricing based on your needs:
              </Typography>
              <List>
                <ListItem>
                  <ListItemIcon><Check /></ListItemIcon>
                  <ListItemText primary="Unlimited analyses and storage" />
                </ListItem>
                <ListItem>
                  <ListItemIcon><Check /></ListItemIcon>
                  <ListItemText primary="Dedicated infrastructure" />
                </ListItem>
                <ListItem>
                  <ListItemIcon><Check /></ListItemIcon>
                  <ListItemText primary="24/7 priority support" />
                </ListItem>
                <ListItem>
                  <ListItemIcon><Check /></ListItemIcon>
                  <ListItemText primary="Custom integrations" />
                </ListItem>
              </List>
              
              <TextField
                fullWidth
                label="Company Name"
                sx={{ mt: 2 }}
              />
              <TextField
                fullWidth
                label="Contact Email"
                sx={{ mt: 2 }}
              />
              <TextField
                fullWidth
                label="Expected Monthly Volume"
                sx={{ mt: 2 }}
                multiline
                rows={2}
              />
            </Box>
          ) : (
            <Box>
              <Typography variant="body1" gutterBottom>
                New plan features:
              </Typography>
              <List dense>
                <ListItem>
                  <ListItemText 
                    primary={`${selectedPlan?.features.analyses_per_month} analyses per month`}
                  />
                </ListItem>
                <ListItem>
                  <ListItemText 
                    primary={`${selectedPlan?.features.storage} storage`}
                  />
                </ListItem>
                <ListItem>
                  <ListItemText 
                    primary={`${selectedPlan?.features.api_calls} API calls`}
                  />
                </ListItem>
              </List>
              
              <Alert severity="success" sx={{ mt: 2 }}>
                <Typography variant="body2">
                  🎉 Start with a <strong>14-day free trial</strong>! No charge until trial ends.
                </Typography>
              </Alert>
              
              <Typography variant="h6" sx={{ mt: 2 }}>
                Total: {selectedPlan?.price_display}
              </Typography>
            </Box>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setUpgradeDialog(false)}>
            Cancel
          </Button>
          <Button 
            variant="contained" 
            onClick={processUpgrade}
            disabled={paymentProcessing}
            startIcon={paymentProcessing ? <CircularProgress size={16} /> : <CreditCard />}
          >
            {selectedPlan?.id === 'enterprise' ? 'Contact Sales' : 'Confirm Upgrade'}
          </Button>
        </DialogActions>
      </Dialog>
    </Box>
  );
};

export default Billing;
