import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Datasets: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Datasets
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Datasets management page</Typography>
      </Paper>
    </Box>
  );
};

export default Datasets;
