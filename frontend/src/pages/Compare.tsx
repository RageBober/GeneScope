import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Compare: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Compare Analyses
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Compare analyses page</Typography>
      </Paper>
    </Box>
  );
};

export default Compare;
