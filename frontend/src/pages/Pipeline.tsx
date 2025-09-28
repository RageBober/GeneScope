import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const Pipeline: React.FC = () => {
  return (
    <Box>
      <Typography variant="h4" sx={{ fontWeight: 700, mb: 3 }}>
        Pipeline Builder
      </Typography>
      <Paper sx={{ p: 3 }}>
        <Typography>Pipeline Builder page</Typography>
      </Paper>
    </Box>
  );
};

export default Pipeline;
