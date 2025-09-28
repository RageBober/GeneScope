import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const References: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Reference Genomes
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Reference genomes page</Typography>
      </Paper>
    </Box>
  );
};

export default References;
