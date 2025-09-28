import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Samples: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Samples
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Samples management page</Typography>
      </Paper>
    </Box>
  );
};

export default Samples;
