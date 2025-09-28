import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Reports: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Reports
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Reports page content will be here</Typography>
      </Paper>
    </Box>
  );
};

export default Reports;
